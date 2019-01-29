module toem
    use iso_c_binding, only: c_double, c_int
    use general
    use unsaturated
    implicit none

contains

    subroutine run_toem(scond, vgl, vga, vgm, vgn, dtheta, height, qin, s0, &
                        nt, tout, sout) bind(c)

        !DEC$ ATTRIBUTES DLLEXPORT :: run_toem

        real(c_double), intent(in) :: scond, vgl, vga, vgm, vgn, dtheta, &
            height, qin, s0
        integer(c_int) :: nt
        real(c_double), dimension(nt), intent(in) :: tout
        real(c_double), dimension(nt), intent(out) :: sout

        type(vangenuchten) :: vg
        real(dp) :: denom
        integer :: itries, it
        ! LSODA parameters
        integer, dimension(64) ::  iwork
        integer :: neq, itol, itask, istate, iopt, lrw, liw, jt
        real(dp) :: tcur, tnext, rtol, atol
        real(dp), dimension(64) :: rwork
        external lsoda

        vg = vangenuchten(vgl, vga, vgn, vgm)
        denom = dtheta * height

        rwork(:) = -5D0                         ! Useful for debugging
        iwork(:) = 15                           ! Useful for debugging
        neq = 1                                 ! Number of equations
        rtol = 1.0D-7                           ! Relative tolerance
        atol = 1.0D-7                           ! Absolute tolerance
        itol = 1                                ! ATOL is scalar
        itask = 1                               ! Initial task, always 1
        istate = 1                              ! Initial state
        iopt = 0                                ! No options
        lrw = 64                                ! Size of RWORK
        liw = 64                                ! Size of IWORK
        jt = 2                                  ! No Jacobian
        sout(1) = s0
        sout(2) = s0
        itries = 0

        do it = 2, nt
            tcur = tout(it-1)
            tnext = tout(it)
101         call dlsoda(feq, neq, sout(it), tcur, tnext, itol, rtol, atol, itask, &
                istate, iopt, rwork, lrw, iwork, liw, jdum, jt)

            if (istate == -1 .and. itries < 20) then
                istate = 2
                itries = itries + 1
                go to 101
            end if
            if (istate < 0) go to 201

            sout(it+1) = sout(it)
        end do

        return

201     write (*,*) 'Error occurred'
        return


    contains

        subroutine feq(neq, t, s, ds)
            integer neq
            real(dp), intent(in) :: t
            real(dp), dimension(neq), intent(in) :: s
            real(dp), dimension(neq), intent(out) :: ds

            ds(1) = (qin - scond * vg%rcond(s(1))) / denom
        end subroutine feq

        subroutine jdum
        end subroutine jdum


    end subroutine run_toem

end module toem
