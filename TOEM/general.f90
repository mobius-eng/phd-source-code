!> \brief General utility functions
!!
module general
    implicit none

    !> \brief Kind of double precision
    integer, parameter :: dp = kind(1.0D0)

    !> \brief Generic function producing a vector with equally spaced items
    !!
    !! Possible calls:
    !!
    !! ```fortran
    !!      result = linspace(x0, xend, [n], xout)
    !! ```
    !!
    !! If `n` is not provided, `xout` must be already allocated.
    !! If `n` is provided, a new space for `xout` will allocated.
    !! `result` is `.true.` if  successful and is `.false.` if
    !! allocation didn't succeed.
    interface linspace
        module procedure linspace_x, linspace_n
    end interface

contains

    !> \brief Produces vector with equally spaced items
    !!
    !! \param[in]  x0       Initial (first) value
    !! \param[in]  xend     End (last) value
    !! \param[in]  nx       Number of items to produce
    !! \param[in,out] xout  Output vector, will be allocated
    !! \result              `.true.` if successful and `.false.`
    !!                      if error has occurred.
    !!
    !! This function allocates the space for the result `xout`. If the
    !! allocation was unsuccessful, it will return `.false.`
    logical function linspace_x(x0, xend, nx, xout)

        real(dp), intent(in) :: x0, xend
        integer, intent(in) :: nx
        real(dp), dimension(:), allocatable, intent(out) :: xout
        integer :: aerr

        allocate(xout(nx), stat=aerr)

        if (aerr > 0) then
            linspace_x = .false.
            return
        end if

        linspace_x = linspace_n(x0,xend,xout)

    end function linspace_x

    !> \brief Produces vector with equally spaced items
    !!
    !! \param[in]  x0       Initial (first) value
    !! \param[in]  xend     End (last) value
    !! \param[in,out] xout  Output vector
    !! \result              Always returns `.true.`
    !!
    !! Fills the input vector `xout` with equally spaced items
    !! from `x0` to `xend` (both included).
    logical function linspace_n(x0,xend,xout)
        real(dp), intent(in) :: x0,xend
        real(dp), dimension(:), intent(out) :: xout

        integer ::  i, nx
        real(dp) :: dx

        nx = size(xout)
        dx = (xend - x0) / nx

        xout(1) = x0
        do i=1,nx-2
            xout(i+1) = xout(i) + dx
        end do
        xout(nx) = xend

        linspace_n = .true.
    end function

end module general
