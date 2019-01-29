module driver

    use iso_c_binding, only: c_double, c_int
    use general
    use constants
    use geometry2d
    use darcy_mod
    use unsaturated
    use inletflow_mod
    implicit none

contains

    subroutine init_and_run_model(height, width, l, a, vgm, vgn, scond, qin, s0, dtheta, tfinal, &
        isde, ilensing, ichannelling, &
        jdrip, &
        ilens, jlensbeg, jlensend, kminlens, kmaxlens, &
        ichanbeg, ichanend, jchan, kmaxchan, &
        tpoints, sat, qina, qouta) bind(c)

        !DEC$ ATTRIBUTES DLLEXPORT :: init_and_run_model

        real(c_double) :: height, width, l, a, vgm, vgn, scond, qin, s0, dtheta, tfinal, &
            kminlens, kmaxlens, kmaxchan
        integer(c_int) :: isde, ilensing, ichannelling, jdrip, ilens, jlensbeg, jlensend, &
            ichanbeg, ichanend, jchan
        real(c_double), dimension(nt), intent(inout) :: tpoints
        real(c_double), dimension(m,n,nt), intent(inout) :: sat
        real(c_double), dimension(n, nt), intent(inout) :: qina, qouta

        type(domain) :: dom
        real(dp) :: dx, dy, dz = 1.0_dp
        type(darcy) :: drc
        type(vangenuchten) :: vg
        type(const_inletflow) :: q
        logical, dimension(n) :: mask
        logical :: lresult
        real(dp), dimension(m,n) :: sat0

        vg = vangenuchten(l, a, vgn, vgm)

        dx = width / n
        dy = height / m
        call init_rect(dx, dy, dz, dom)

        tpoints(1) = 0.0_dp
        lresult = linspace(30.0_dp, tfinal, tpoints(2:nt))

        if (isde == 0) then
            mask(:) = .true.
        else
            mask(:) = .false.
            mask(jdrip) = .true.
        end if
        q%constq = qin
        q%mask(:) = mask(:)

        call init_uniform_darcy(dom, vg, dtheta, scond, drc)

        if (ilensing /= 0) then
            call create_lensing(ilens, jlensbeg, jlensend, kminlens, kmaxlens, drc)
        end if

        if (ichannelling /= 0) then
            call create_channel(ichanbeg, ichanend, jchan, kmaxchan, drc)
        end if

        sat0 = s0

        call run_darcy(drc, q, sat0, tpoints, sat, qina, qouta)

    end subroutine


end module driver
