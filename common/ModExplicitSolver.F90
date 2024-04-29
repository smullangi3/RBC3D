module ModDirectSolver

  use ModDataStruct
  use ModDataTypes
  use ModData
  use ModBasicMath
  use MPI

  implicit none

  public :: DirectSolve
contains

  subroutine DirectSolve(slist, tlist, results)
    type(t_SourceList) :: slist
    type(t_TargetList) :: tlist
    real(WP) :: results(:, :) ! dimension of (tlist%nPoint, 3)
    real(WP), allocatable :: v_tmp(:, :)
    integer :: i, j, l, ierr
    integer :: s_it, t_it
    real(WP) :: S(3, 3), T(3, 3, 3)
    real(WP) :: sp(3), tp(3), f(3), h(3), v(3), u(3)
    real(WP) :: A, B
    real(WP) :: dist

    integer :: nodeNum, numNodes
    call MPI_Comm_Rank(MPI_COMM_WORLD, nodeNum, ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, numNodes, ierr)
    nodeNum = nodeNum + 1

    allocate (v_tmp(tlist%nPoint, 3))
    v_tmp = 0

    do t_it = 1, tlist%nPoint
      if (.not. tlist%active(t_it)) cycle
      tp = tlist%x(t_it, :)
      u = 0
      A = tlist%Acoef(t_it)
      do s_it = nodeNum, slist%nPoint, numNodes
        sp = slist%x(s_it, :)
        f = slist%f(s_it, :)
        h = slist%g(s_it, :)
        v = slist%a3(s_it, :)
        B = slist%Bcoef(s_it)

        dist = VecNorm(tp - sp)
        if (dist .eq. 0) cycle
        S = 0
        T = 0

        do i = 1, 3
        do j = 1, 3
          !stokeslet
          S(i, j) = (tp(i) - sp(i))*(tp(j) - sp(j))/(dist**3)
          if (i .eq. j) then
            S(i, j) = S(i, j) + 1/dist
          end if

          do l = 1, 3
            T(i, j, l) = (tp(i) - sp(i))*(tp(j) - sp(j))*(tp(l) - sp(l))/(dist**5)
            u(i) = u(i) + T(i, j, l)*h(j)*v(l)*B/(4*PI*A)
          end do !l
          u(i) = u(i) + S(i, j)*f(j)/(-4*PI*A)
        end do !j
        end do !i

      end do !s_it
      ! write(*,*) "Velocity: " , u
      v_tmp(t_it, :) = v_tmp(t_it, :) + u
    end do !t_it

    ! MPI allreduce to collect velocities
    call MPI_Allreduce(MPI_IN_PLACE, v_tmp, tlist%nPoint*3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
    results = results + v_tmp; 
    deallocate (v_tmp)
  end subroutine DirectSolve

end module ModDirectSolver
