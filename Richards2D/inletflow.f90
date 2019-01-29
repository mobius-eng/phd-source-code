!> \brief Inlet flow related functions
module inletflow_mod

    use general
    use constants
    use growing_array_mod
    implicit none


    !> \brief   Base type for all inlet flow types. Defaults to no flow.
    !!
    !! Inlet flow is represented by a polymorphic data type to provide
    !! common and extensible interface to different inlet floew regimes.
    type inletflow
    contains
        !> \brief Inelt flow function
        !!
        !! General form:
        !! ```
        !!      result = flow%q(t)
        !! ```
        !! Note, that `flow` may change as a result of calling `q`.
        procedure :: q => default_q
    end type


    !> \brief    Constant inlet flow
    type, extends(inletflow) :: const_inletflow
        !> \brief   The value of the flow rate (assumed in m/s)
        !INTEGER, LEN :: N
        real(dp) :: constq
        logical, dimension(n) :: mask
    contains
        !> \brief   Overriding flow function
        procedure :: q => const_inletflow_q
    end type


contains

    !> \brief   Default zero flow for base `inletflow` type
    subroutine default_q(qin, t, q)
        class(inletflow), intent(inout) :: qin
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(out) :: q
        q(:) = 0D0
    end subroutine default_q


    !> \brief   Constant inlet flow for `const_inletflow` type
    subroutine const_inletflow_q(qin, t, q)
        class(const_inletflow), intent(inout) :: qin
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(out) :: q

        ! logical, dimension(n) :: mask

        integer :: i

        do i=1,n
            if (qin%mask(i)) then
                q(i) = qin%constq
            else
                q(i) = 0.0_dp
            end if
        end do

    end subroutine const_inletflow_q


end module inletflow_mod
