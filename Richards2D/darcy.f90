module darcy_mod

    use general
    use constants
    use geometry2d
    use inletflow_mod
    use unsaturated

    implicit none

    abstract interface
        subroutine iqin(t, q)
            import dp
            real(dp) :: t
            real(dp), dimension(:) :: q
        end subroutine
    end interface



    type darcy
        type(domain), pointer :: dom
        type(vangenuchten), dimension(m, n) :: umdl
        type(vangenuchten), dimension(m,n+1) :: umdlx
        type(vangenuchten), dimension(m+1,n) :: umdly
        real(dp), dimension(m, n + 1) :: jx
        real(dp), dimension(m + 1, n) :: jy
        real(dp), dimension(m, n+1) :: scondx
        real(dp), dimension(m+1,n) :: scondy
        real(dp), dimension(m, n) :: pres
        real(dp), dimension(m, n) :: sat
        real(dp), dimension(m,n ) :: dtheta

    end type


contains

    subroutine init_uniform_darcy(dom, umdl, dtheta, scond, drc)
        type(darcy), intent(inout) :: drc
        type(domain), intent(in), target :: dom
        type(vangenuchten), intent(in) :: umdl
        real(dp), intent(in) :: dtheta, scond

        drc%dom => dom

        drc%umdl = umdl
        drc%umdlx = umdl
        drc%umdly = umdl
        drc%dtheta = dtheta
        drc%scondx = scond
        drc%scondy = scond
    end subroutine init_uniform_darcy

    subroutine create_lensing(i, jbeg, jend, kmin, kmax, drc)
        integer, intent(in) :: i, jbeg, jend
        real(dp), intent(in) :: kmin, kmax
        type(darcy), intent(inout) :: drc

        drc%scondx(i, 1:jbeg) = kmin
        drc%scondx(i, jbeg+1:jend-1) = kmax
        drc%scondx(i, jend:n+1) = kmin

        drc%scondy(i:i+1, 1:jbeg-1) = kmin
        drc%scondy(i:i+1, jbeg:jend) = kmax
        drc%scondy(i:i+1, jend+1:n) = kmin

    end subroutine create_lensing

    subroutine create_channel(ibeg, iend, j, kmax, drc)
        integer, intent(in) :: ibeg, iend, j
        real(dp), intent(in) :: kmax
        type(darcy), intent(inout) :: drc

        drc%scondy(ibeg:iend, j) = kmax

    end subroutine create_channel


    subroutine run_darcy(drc, qin, sat0, tpoints, satout, qina, qouta)
        type(darcy), intent(inout) :: drc
        real(dp), dimension(:,:), intent(in) :: sat0
        real(dp), dimension(:), intent(in) :: tpoints
        real(dp), dimension(:,:,:), intent(out), target :: satout
        real(dp), dimension(:,:), intent(out) :: qina, qouta
        class(inletflow), intent(inout) :: qin
        external dlsoda

        !INTERFACE
        !    SUBROUTINE DLSODA(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
        !            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
        !        USE GENERAL
        !        EXTERNAL :: F, JAC
        !        INTEGER :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, JT, IWORK
        !        REAL(DP) :: T, TOUT, RTOL, Y, ATOL, RWORK
        !        DIMENSION Y(*), RWORK(LRW), IWORK(LIW)
        !    END SUBROUTINE DLSODA
        !END INTERFACE

        integer :: mn, i, j, ierr, itrys, it, itcur
        real(dp) :: tcur, tnext
        ! real(dp), dimension(:), pointer :: psat
        ! LSODA DATA
        real(dp), dimension(:), allocatable :: rwork
        integer, dimension(:), allocatable :: iwork
        integer :: neq, itol, itask, istate, iopt, lrw, liw, jt
        real(dp) :: rtol, atol
        ! EXTERNAL JDUM


        mn = m * n
        it = size(tpoints)
        ! PSAT(1:MN) => DRC%SAT

        ! Allocate storage for LSODA
        allocate (iwork(20+mn), rwork(100 + 22 + mn * max(16, mn + 9)), stat=ierr)
        ! IF (ERR > 0) GO TO 100
        ! Init LSODA data
        rwork(:) = -5D0                         ! Useful for debugging
        iwork(:) = 15                           ! Useful for debugging
        iwork(1) = m                            ! Upper and lower bands
        iwork(2) = m                            !   in Jacobian
        neq = mn                                ! Number of equations
        rtol = 1.0D-5                           ! Relative tolerance
        atol = 1.0D-5                           ! Absolute tolerance
        itol = 1                                ! ATOL is scalar
        itask = 1                               ! Initial task, always 1
        istate = 1                              ! Initial state
        iopt = 0                                ! No options
        lrw = 22 + 100 + mn * max(16, mn + 9)   ! Size of RWORK
        liw = 20 + mn                           ! Size of IWORK
        jt = 5                                  ! No Jacobian, use banded 5
        drc%sat(:,:) = sat0(:,:)                ! Set initial value
        satout(:,:,1) = sat0(:,:)               ! Store initial value
        call qin%q(0D0, qina(:,1))             ! The first value of flow in
        itrys = 0

        do itcur=1,it-1
            tcur = tpoints(itcur)
            tnext = tpoints(itcur+1)
            associate(sat => drc%sat)
10              call dlsoda(fdarcy, neq, sat, tcur, tnext, itol, rtol, atol, itask, &
                    istate, iopt, rwork, lrw, iwork, liw, jdum, jt)

                if (istate == -1 .and. itrys < 20) then
                    write (*,*) " "
                    write(*,*) "Too much work in DLSODA. Trying to restart the step"
                    write(*, 501) tcur, tnext
                    501 format ('TCUR = ', f10.4, '  TNEXT = ',f10.4)
                    istate = 2
                    itrys = itrys + 1
                    go to 10
                end if
                if (itrys == 20) then
                    write (*,*) "Ran out of restarts, giving up"
                end if
            end associate
            if (istate < 0) go to 110               ! Error has occurred
            satout(:, :, itcur+1) = drc%sat(:,:)    ! Record saturation
            call qin%q(tnext, qina(:, itcur+1))
            ! write (*,'(A)', advance='NO') "="     ! Make a mark
            !20      FORMAT(' At T =',D8.0,'   S(2) =',D14.6)
        end do
        forall (j=1:n, i=1:it)
            qouta(j, i) = drc%scondy(m+1, j) * drc%umdly(m+1, j)%rcond(satout(m, j, i))
        end forall

110     deallocate(rwork, iwork, stat=ierr)


    contains

        subroutine fdarcy(neq, t, sat, dsat)
            integer, intent(in) :: neq
            real(dp), intent(in) :: t
            real(dp), dimension(neq), intent(in), target :: sat
            real(dp), dimension(neq), intent(out), target :: dsat

            real(dp), dimension(:,:), pointer, contiguous :: s
            real(dp), dimension(:, :), pointer, contiguous :: ds
            integer :: i, j
            real(dp) :: x, h, sb

            s(1:m, 1:n) => sat(1:neq)
            ds(1:m,1:n) => dsat(1:neq)

            ! Setup all the quantities on centroids
            forall(i = 1 : m, j = 1 : n)
                drc%pres(i,j) = drc%umdl(i,j)%cpres(s(i,j))
                ! DRC%COND(I,J) = DRC%UMDL(I,J)%RCOND(S(I,J))
            end forall
            ! Calculate fluxes
            ! Boundary fluxes:
            drc%jx(:,1) = 0D0
            drc%jx(:, n+1) = 0D0
            call qin%q(t, drc%jy(1,:))
            forall (j=1:n) drc%jy(m+1,j) = drc%scondy(m+1, j) * drc%umdly(m+1, j)%rcond(s(m, j))
            do j=2,n
                do i=1,m
                    h = drc%dom%dcenx(i,j-1)
                    x = drc%dom%dcenfacx(i,j-1)
                    ! h = drc%dom%cen(i,j)%x - drc%dom%cen(i,j-1)%x
                    ! x = drc%dom%facx(i,j)%x - drc%dom%cen(i,j-1)%x
                    ! Saturation on the boundary
                    sb = (s(i,j) - s(i,j-1))*x/h + s(i,j-1)
                    drc%jx(i,j) = - drc%scondx(i,j) * drc%umdlx(i,j)%rcond(sb) * (drc%pres(i,j) - drc%pres(i,j-1)) / h
                end do
            end do

            do j=1,n
                do i=2,m
                    h = drc%dom%dceny(i-1,j)
                    x = drc%dom%dcenfacy(i-1,j)
!                    h = drc%dom%cen(i,j)%y - drc%dom%cen(i-1,j)%y
!                    x = drc%dom%facy(i,j)%y - drc%dom%cen(i-1,j)%y
                    sb = (s(i,j) - s(i-1,j))*x/h + s(i-1,j)
                    drc%jy(i,j) = - drc%scondy(i,j) * drc%umdly(i,j)%rcond(sb) * ((drc%pres(i,j) - drc%pres(i-1,j)) / h - 1D0)
                end do
            end do
            ! Mass balance for each cell
            do j=1,n
                do i=1,m
                    x = drc%dom%dfacx(i,j)
                    h = drc%dom%dfacy(i,j)
!                    x = drc%dom%facx(i,j+1)%x - drc%dom%facx(i,j)%x
!                    h = drc%dom%facy(i+1,j)%y - drc%dom%facy(i,j)%y
                    ds(i,j) = -1D0/drc%dtheta(i,j) * ( (drc%jx(i,j+1) - drc%jx(i,j))/x + &
                                                       (drc%jy(i+1,j) - drc%jy(i,j))/h )
                end do
            end do

        end subroutine fdarcy

    end subroutine run_darcy

    subroutine jdum
    end subroutine

end module darcy_mod
