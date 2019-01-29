!>@brief            Unsaturated models
module unsaturated

    use general
    use constants
    implicit none

    !>@brief                Van Genuchten model representation
    type, public :: vangenuchten
        !>@brief    Mualem parameter, usually 0.5
        real(dp) :: l
        !>@brief    Inverse of bubbling pressure, alpha
        real(dp) :: a
        !>@brief    Van Genuchten fitting parameter
        real(dp) :: n
        !>@brief    Van Genuchten fitting parameter
        !!
        !! This parameter is usually set as `1-1/n`. Here, however,
        !! it is allowed to be set independently.
        real(dp) :: m
    contains
        procedure, pass(this), public :: average => vg_average
        procedure, pass(this), public :: cpres => vg_cpres
        procedure, pass(this), public :: rcond => vg_rcond
        procedure, pass(this), public :: sat => vg_sat
        procedure, pass(this), public :: copyto => vg_copyto
        procedure, pass(recip), public :: vg_assignment
        generic ::  assignment (=) => vg_assignment
    end type vangenuchten


    !>@brief    Generic van Genuchten model constructor (factory)
    !!
    !! Description
    !!   Convenience function to construct van Genuchten model.
    !!   Has three forms:
    !!   1.
    !!   \code{.f90}
    !!       MDL = MKVANGENUCHTEN(L, A, N[, M])
    !!   \endcode
    !!   with L, A, N, M being REAL(DP) - constructs 1-dimensional model
    !!   2.
    !!       MDL = MKVANGENUCHTEN(I, L, A, N[, M])
    !!   the same as 1, but additionally I is an INTEGER defining the
    !!   dimensionality of the model.
    !!   3.
    !!       MDL = MKVANGENUCHTEN(VL, VA, VN, VM)
    !!   with VL, VA, VN, VM being vectors of REAL(DP) - constructs
    !!   van Genuchten model with given vector parameters.
    !!
    !!   If the form 1 or 2 is used, DELUMODEL must be used at the end of
    !!   computation.
    interface mkvangenuchten
        module procedure mkvangenuchten1, mkvangenuchtenn, mkvangenuchtenv
    end interface


contains

    subroutine vg_assignment(recip, source)
        class(vangenuchten), intent(out) :: recip
        class(vangenuchten), intent(in) :: source
        select type (source)
        type is (vangenuchten)
            recip%l = source%l
            recip%a = source%a
            recip%n = source%n
            recip%m = source%m
        end select
    end subroutine



    elemental real(dp) function vg_rcond(this, sat)
        ! VG_RCOND
        !   Calculates relative conductivity for the van Genuchten model
        !
        ! Arguments
        !   SAT         Input saturation
        !   MDL         (Input) model
        !   KOND        Output: relative conductivities
        !
        class(vangenuchten), intent(in) :: this
        real(dp), intent(in) :: sat

        real(dp) :: s

        if (sat >= 1D0) then
            vg_rcond = 1D0
            return
        end if

        s = max(0.001D0, sat)

        vg_rcond = (s ** this%l) * &
            (1.0D0 - (1.0D0 - s ** (1.0D0 / this%m)) ** this%m) ** 2

    end function vg_rcond


    elemental real(dp) function vg_cpres(this, sat)
        ! VG_CPRES
        !   Calculates capillary pressure for van Genuchten model
        !
        ! Arguments
        !   SAT         Input saturation
        !   MDL         (Input) model
        !   PRES        Output: capillary pressure
        !
        real(dp), intent(in) :: sat
        class(vangenuchten), intent(in) :: this

        real(dp) :: a, n, l, m
        real(dp) :: s

        a = this%a
        n = this%n
        l = this%l
        m = this%m


        if (sat >= 1D0) then
            vg_cpres = 0D0
            return
        end if

        s = max(0.001D0, sat)

        vg_cpres = - 1.0D0 / this%a * &
            (s ** (-1.0D0 / this%m) - 1.0D0) ** (1.0D0 / this%n)

    end function vg_cpres


    elemental real(dp) function vg_sat(this, pres)
        ! V_SAT
        !   Calculates saturation from pressure of van Genuchten model
        !
        ! Arguments
        !   PRES        Input capillary pressure
        !   MDL         (Input) van Genuchten model
        !   PRES        Output: capillary pressure
        !
        real(dp), intent(in) :: pres
        class(vangenuchten), intent(in) :: this

        vg_sat = (1.0D0 / (1.0D0 + abs(this%a * pres) ** this%n)) ** this%m

    end function vg_sat


    function mkvangenuchten1(l, a, n, m) result(vg)
        ! MKVANGENUCHTEN1
        !   Constructs one-dimensional van Genuchten model.
        !
        ! Arguments
        !   L, A, N, M  Input REAL(DP) model parameters
        !                   M is optional, if not provided: M = 1-1/N
        !
        ! Result
        !   van Genuchten model (TYPE(VANGENUCHTEN))
        !
        real(dp), intent(in) :: l, a, n
        real(dp), intent(in), optional :: m

        type(vangenuchten) :: vg

        real(dp) :: m1

        if (.not. present(m)) then
            m1 = 1.0D0 - 1.0D0 / n
        else
            m1 = m
        end if

        vg%l = l; vg%a = a; vg%n = n; vg%m = m1

    end function mkvangenuchten1


    function mkvangenuchtenn(i, l, a, n, m) result(vg)
        ! MKVANGENUCHTENN
        !   Constructs I-dimensional van Genuchten model.
        !   It allocates space for the model in a vector - DELUMODEL
        !   must be used at the end of computation.
        !
        ! Arguments
        !   I           Model dimensionality
        !   L, A, N, M  Input REAL(DP) model parameters
        !                   M is optional, if not provided: M = 1-1/N
        !
        ! Result
        !   Vector of van Genuchten models (TYPE(VANGENUCHTEN), DIMENSION(I))
        !
        integer, intent(in) :: i
        real(dp), intent(in) :: l, a, n
        real(dp), intent(in), optional :: m

        type(vangenuchten), dimension(:), allocatable :: vg

        integer :: err
        real(dp) :: m1

        allocate (vg(i), stat=err)

        if (err > 0) then
            write (*,*) "MKVANGENUCHTENN: Cannot allocate space. Stopping"
            stop
        end if

        if (.not. present(m)) then
            m1 = 1.0D0 - 1.0D0 / n
        else
            m1 = m
        end if

        vg(:)%l = l; vg(:)%a = a; vg(:)%n = n; vg(:)%m = m1

    end function mkvangenuchtenn


    function mkvangenuchtenv(l, a, n, m) result(vg)
        ! MKVANGENUCHTENV
        !   Constructs general multi-dimensional van Genuchten model.
        !   Allocates the space for the model - use DELUMODEL at the end of
        !   computation.
        !
        ! Arguments
        !   L, A, N, M  Input REAL(DP) model parameters
        !               NOTE: M is compulsary in this form!
        !
        ! Result
        !   Vector of van Genuchten models (TYPE(VANGENUCHTEN))
        !
        real(dp), dimension(:), intent(in) :: l, a, n
        real(dp), dimension(:), optional, intent(in) :: m

        type(vangenuchten), dimension(:), allocatable :: vg

        integer :: i, err

        i = size(l)
        allocate (vg(i), stat=err)
        if (err > 0) then
            write (*,*) "VANGENUCHTENV: Cannot allocate space for the model. Stopping"
            stop
        end if

        vg(:)%l = l(:); vg(:)%a = a(:); vg(:)%n = n(:)
        if (present(m)) then
            vg(:)%m = m(:)
        else
            vg(:)%m = 1D0 - 1D0/n(:)
        end if

    end function mkvangenuchtenv


    subroutine vg_average(this, mdl2, mdlres)
        class(vangenuchten), intent(in) :: this
        class(vangenuchten), intent(in) :: mdl2
        class(vangenuchten), intent(inout) :: mdlres

        select type (mdl2)
        type is (vangenuchten)
            select type (mdlres)
            type is (vangenuchten)
                mdlres%l = 0.5d0*(this%l + mdl2%l)
                mdlres%a = 0.5d0*(this%a + mdl2%a)
                mdlres%n = 0.5d0*(this%n + mdl2%n)
                mdlres%m = 0.5d0*(this%m + mdl2%m)
            end select
        end select

    end subroutine vg_average

    elemental subroutine vg_copyto(this, dest)
        class(vangenuchten), intent(in) :: this
        class(vangenuchten), intent(inout) :: dest
        select type (dest)
        type is (vangenuchten)
            dest%l = this%l
            dest%a = this%a
            dest%n = this%n
            dest%m = this%m
        end select
    end subroutine


end module unsaturated
