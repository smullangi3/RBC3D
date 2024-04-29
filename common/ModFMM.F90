module ModFMM

  use ModDataTypes
  use ModData
  use ModDataStruct
  use MPI

  implicit none

  public :: FMM_AddVelocity

contains

  subroutine FMM_AddVelocity(slist, tlist, v, c1, c2)
    type(t_SourceList) :: slist
    type(t_TargetList) :: tlist
    real(WP) :: v(:, :) ! is dimension (3, tlist_length)
    real(WP) :: c1, c2

    !variables relevant for Stokeslet FMM3D
    integer :: nd, nsource, ntarg
    real(WP) :: eps
    integer :: ifstrslet, ifstoklet
    real(WP), allocatable :: source(:, :), targets(:, :)
    real(WP), allocatable :: stoklet(:, :, :), strslet(:, :, :), strsvec(:, :, :)
    integer ::ifppreg, ifppregtarg
    real(WP), allocatable :: pot(:, :, :), pre(:, :), grad(:, :, :, :)
    real(WP), allocatable :: pottarg(:, :, :), pretarg(:, :), gradtarg(:, :, :, :)
    integer :: ier

    !Temp vars for prepping data
    real(WP), allocatable :: source_T(:, :), stoklet_T(:, :)
    real(WP), allocatable :: strslet_T(:, :), strsvec_T(:, :)
    integer :: i, ii

    real(WP) :: FMM_clock

    eps = 1E-6
    ifstoklet = merge(1, 0, c1 .ne. 0)
    ifstrslet = merge(1, 0, c2 .ne. 0)

    nd = 1
    nsource = slist%nPoint
    ntarg = size(tlist%x)/3

    allocate (source_T(nsource, 3), source(3, nsource), targets(3, ntarg))
    allocate (stoklet_T(nsource, 3), stoklet(nd, 3, nsource))
    allocate (strslet_T(nsource, 3), strslet(nd, 3, nsource))
    allocate (strsvec_T(nsource, 3), strsvec(nd, 3, nsource))
    allocate (pot(nd, 3, nsource), pre(nd, nsource), grad(nd, 3, 3, nsource))
    allocate (pottarg(nd, 3, ntarg), pretarg(nd, ntarg), gradtarg(nd, 3, 3, ntarg))

    forall (i=1:nsource)
      source_T(i, :) = slist_rbc%x(i, :)
      stoklet_T(i, :) = slist_rbc%f(i, :)
      strslet_T(i, :) = slist_rbc%g(i, :)*slist_rbc%Bcoef(i)
      strsvec_T(i, :) = slist_rbc%a3(i, :)
    end forall

    source = transpose(source_T)
    stoklet(nd, :, :) = transpose(stoklet_T)
    strslet(nd, :, :) = transpose(strslet_T)
    strsvec(nd, :, :) = transpose(strsvec_T)
    targets = transpose(tlist%x)

    ! carry the coefficients inward
    stoklet = c1*stoklet
    strslet = c2*strslet

    ifppreg = 0
    ifppregtarg = 1

    pottarg = 0
    FMM_clock = MPI_WTime()
    call stfmm3d(nd, eps, nsource, source, ifstoklet, stoklet, &
                 ifstrslet, strslet, strsvec, ifppreg, pot, pre, grad, &
                 ntarg, targets, ifppregtarg, pottarg, pretarg, gradtarg, ier)

    ! if (rootWorld) write(*,*) "FMM Time Cost: ", (MPI_WTime() - FMM_clock), " Error Code: " , ier, &
    !     " Max Velocity: ", maxval(pottarg), &
    !     " Stoklet: ", maxval(stoklet), " Strslet", maxval(strslet), " Strsvec", maxval(strsvec), &
    !     " If Strslet: ", ifstrslet, " If Stoklet: ", ifstoklet

    do i = 1, tlist%nPoint
      if (tlist%active(i)) then
        v(i, :) = v(i, :) + (pottarg(nd, :, i)/tlist%Acoef(i))
      end if
    end do

    deallocate (source_T, source, targets)
    deallocate (stoklet_T, stoklet)
    deallocate (strslet_T, strslet)
    deallocate (strsvec_T, strsvec)
    deallocate (pot, pre, grad, pottarg, pretarg, gradtarg)
  end subroutine FMM_AddVelocity

end module ModFMM
