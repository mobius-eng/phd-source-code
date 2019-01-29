module geometry2d

    use general
    use constants
    implicit none

    type point
        real(dp) :: x, y
    end type

    type domain
        type(point), dimension(m,n) :: cen
        type(point), dimension(m,n+1) :: facx
        type(point), dimension(m+1,n) :: facy
        real(dp), dimension(m,n) :: vol
        real(dp), dimension(m,n+1) :: ax
        real(dp), dimension(m+1,n) :: ay
        real(dp), dimension(m,n-1) :: dcenx
        real(dp), dimension(m-1,n) :: dceny
        real(dp), dimension(m,n) :: dcenfacx, dcenfacy
        real(dp), dimension(m,n) :: dfacx, dfacy
    end type


contains

    subroutine init_rect(dx, dy, dz, dom)
        real(dp), intent(in) :: dx, dy, dz
        type(domain), intent(out) :: dom

        real(dp) :: vol, ax, ay, cx, cy, dx2, dy2
        integer :: i, j
        real(dp) :: dfx, dfy
        real(dp) :: x1, x2, y1, y2

        vol = dx * dy * dz
        ax = dy * dz
        ay = dx * dz

        dx2 = dx / 2
        dy2 = dy / 2

        dom%vol(:,:) = vol
        dom%ax(:,:)  = ax
        dom%ay(:,:)  = ay

        cx = dx2
        cy = dy2
        do j=1,n
            cy = dy2
            do i=1,m
                dom%cen(i,j) = point(cx, cy)
                dom%facx(i,j) = point(cx - dx2, cy)
                dom%facy(i,j) = point(cx, cy - dy2)
                cy = cy + dy
            end do
            dom%facy(m+1,j) = point(cx, cy - dy2)
            cx = cx + dx
        end do
        forall (i=1:m) dom%facx(i,n+1) = point(cx - dx2, dy2+dy*(i-1))

        do j=1,n-1
            do i=1,m-1
                dom%dcenx(i,j) = dom%cen(i,j+1)%x - dom%cen(i,j)%x
                dom%dceny(i,j) = dom%cen(i+1,j)%y - dom%cen(i,j)%y

                dom%dcenfacx(i,j) = dom%facx(i, j+1)%x - dom%cen(i,j)%x
                dom%dcenfacy(i,j) = dom%facy(i+1, j)%y - dom%cen(i,j)%y

                dom%dfacx(i,j) = dom%facx(i,j+1)%x - dom%facx(i,j)%x
                dom%dfacy(i,j) = dom%facy(i+1,j)%y - dom%facy(i,j)%y
            end do
            dom%dcenx(m,j) = dom%cen(i,j+1)%x - dom%cen(i,j)%x
            dom%dcenfacx(m,j) = dom%facx(m, j+1)%x - dom%cen(m,j)%x
            dom%dfacx(m,j) = dom%facx(m,j+1)%x - dom%facx(m,j)%x
            dom%dfacy(m,j) = dom%facy(m+1,j)%y - dom%facy(m,j)%y
        end do
        do i=1,m-1
            dom%dceny(i,n) = dom%cen(i+1,n)%y - dom%cen(i,n)%y
            dom%dcenfacy(i,n) = dom%facy(i+1, n)%y - dom%cen(i, n)%y
            dom%dfacx(i,n) = dom%facx(i,n+1)%x - dom%facx(i,n)%x
            dom%dfacy(i,n) = dom%facy(i+1,n)%y - dom%facy(i,n)%y
        end do
        dom%dcenfacx(m,n) = dom%facx(m, n+1)%x - dom%cen(m,n)%x
        dom%dcenfacy(m,n) = dom%facy(m+1, n)%y - dom%cen(m, n)%y
        dom%dfacx(m,n) = dom%facx(m,n+1)%x - dom%facx(m,n)%x
        dom%dfacy(m,n) = dom%facy(m+1,n)%y - dom%facy(m,n)%y

        dfx = dom%dfacx(1,n)
        dfy = dom%dfacy(m,1)
        x1 = dom%facx(1,n)%x
        x2 = dom%facx(1,n+1)%x
        y1 = dom%facy(m,1)%y
        y2 = dom%facy(m+1,1)%y
        return

    end subroutine init_rect

    function cendist(dom, i, j, di, dj)
        type(domain), intent(in) :: dom
        integer, intent(in) :: i, j, di, dj
        real(dp) :: cendist

        if (di * dj == 0) then
            cendist = abs(dom%cen(i,j)%x - dom%cen(i+di,j+dj)%x) + &
                abs(dom%cen(i,j)%y - dom%cen(i+di,j+dj)%y)
        else
            cendist = sqrt((dom%cen(i,j)%x - dom%cen(i+di,j+dj)%x)**2 + &
                (dom%cen(i,j)%y - dom%cen(i+di,j+dj)%y)**2)
        end if
    end function cendist

    function facxdist(dom, i, j, di, dj)
        type(domain), intent(in) :: dom
        integer, intent(in) :: i, j, di, dj
        real(dp) :: facxdist

        if (di * dj == 0) then
            facxdist = abs(dom%facx(i,j)%x - dom%facx(i+di,j+dj)%x) + &
                abs(dom%facx(i,j)%y - dom%facx(i+di,j+dj)%y)
        else
            facxdist = sqrt((dom%facx(i,j)%x - dom%facx(i+di,j+dj)%x)**2 + &
                (dom%facx(i,j)%y - dom%facx(i+di,j+dj)%y)**2)
        end if
    end function facxdist

        function facydist(dom, i, j, di, dj)
        type(domain), intent(in) :: dom
        integer, intent(in) :: i, j, di, dj
        real(dp) :: facydist

        if (di * dj == 0) then
            facydist = abs(dom%facy(i,j)%x - dom%facy(i+di,j+dj)%x) + &
                abs(dom%facy(i,j)%y - dom%facy(i+di,j+dj)%y)
        else
            facydist = sqrt((dom%facy(i,j)%x - dom%facy(i+di,j+dj)%x)**2 + &
                (dom%facy(i,j)%y - dom%facy(i+di,j+dj)%y)**2)
        end if
    end function facydist

end module geometry2d
