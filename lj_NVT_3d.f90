  !   =================================================
  !   Standard 3D NVT simulation of Lennard-Jones fluid
  !                   October, 2020
  !   =================================================

  !   ===============================================================================
  !   Author: Anze Hubman
  !           FCCT Univ. of Ljubljana & National Institute of Chemistry/Theory dept.
  !   ===============================================================================

  !   =====================================================     
  !   A simple 3D Monte Carlo simulation of a Lennard-Jones
  !   fluid in NVT ensemble. All parameters are given in
  !   reduced units. A simple potential truncation at r=rc
  !   is used and then the total energy is adjusted using a
  !   tail correction. The only output is the total energy
  !   of a system per particle, namely U*/N. The user can
  !   modify the code freely to produce the desired amount
  !   of output. The code reproduces NIST data given at:
  !   https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
  !   =====================================================


module params
  integer, parameter :: N    = 500        ! number of particles
  real*8,  parameter :: rc   = 3.0d0      ! cutoff radius
  real*8,  parameter :: T    = 0.85d0     ! temperature
  real*8,  parameter :: den  = 0.776d0    ! density
  integer, parameter :: nstp = 5000000    ! number of steps
  real*8,  parameter :: dr   = 0.1d0      ! maximum displacement
end module params

program mc
  use params
  implicit none

  real*8, allocatable :: x(:), y(:), z(:)
  real*8 :: L, acr
  integer :: acc, m

  allocate(x(1:N), y(1:N), z(1:N))

  L = (dble(N)/den)**(1.0d0/3.0d0)
  print *, 'Cell edge:', L
  open(unit=31,file="energy")

  ! initialise particles
  call initial(x,y,z,L)

  ! main run
  do m = 1, nstp

     !===========================
     if (mod(m,10000) == 0) then
        print *, 'Step:', m
     end if
     !===========================
     
     call energy(x,y,z,L)
     call attempt_move(x,y,z,L,acc)
  end do

  ! acceptance ratio
  acr = (dble(acc)/dble(nstp))*100.0d0
  print *, 'Acc. ratio in %:', acr
end program mc

!    =========================================================
!    Initial random arrangement of particles; overlaps allowed
!    =========================================================
subroutine initial(x,y,z,L)
  use params
  implicit none

  real*8, intent(inout) :: x(1:N), y(1:N), z(1:N), L
  integer :: j
  real*8 :: xj, yj, zj

  do j = 1, N
     call random_number(xj)
     call random_number(yj)
     call random_number(zj)
     x(j) = L*xj
     y(j) = L*yj
     z(j) = L*zj
  end do

end subroutine initial

!    ==================================
!    Compute current energy of a system
!    ==================================
subroutine energy(x,y,z,L)
  use params
  implicit none

  real*8, intent(inout) :: x(1:N), y(1:N), z(1:N), L
  integer :: i, j
  real*8 :: dx, dy, dz, dU, U, xi, yi, zi, xj, yj, zj, r2, r2c, halfL, pi, utail

  dU = 0.0d0
  U  = 0.0d0
  r2c = rc**2
  halfL = 0.5d0*L
  pi = 4.0d0*atan(1.d0)
  utail = (8.0d0/3.0d0)*pi*den*(((1.0d0/3.0d0)*((1.0d0/rc)**9)) - ((1.0d0/rc)**3))

  do i = 1, N-1
     do j = i+1, N

        xi = x(i)
        yi = y(i)
        zi = z(i)
        xj = x(j)
        yj = y(j)
        zj = z(j)

        ! minimum image convention
        dx = abs(xi - xj)
        dy = abs(yi - yj)
        dz = abs(zi - zj)

        if (dx > halfL) then
           dx = L - dx
        end if

        if (dy > halfL) then
           dy = L - dy
        end if

        if (dz > halfL) then
           dz = L - dz
        end if

        r2 = dx**2 + dy**2 + dz**2

        if (r2 .le. r2c) then
           dU = 4.0d0*(((1.0d0/r2)**6) - ((1.0d0/r2)**3))
           U  = U + dU
        end if

     end do
  end do
  U = U + dble(N)*utail
  write(31,*) U/dble(N)
end subroutine energy

!     ==========================================
!     Attempt to move a randomly chosen particle
!     ==========================================
subroutine attempt_move(x,y,z,L,acc)
  use params
  implicit none

  real*8, intent(inout) :: x(1:N), y(1:N), z(1:N), L
  integer, intent(inout) :: acc
  real*8 :: q, x_old, y_old, z_old, qx, qy, qz, x_new, y_new, z_new, U_old, U_new
  real*8 :: halfL, dx_new, dy_new, dz_new, dx_old, dy_old, dz_old, r2_old, r2_new, r2c
  real*8 :: dU, w, v, xj, yj, zj, utail, pi
  integer :: k, j

  halfL = 0.5d0*L
  r2c = rc**2
  pi = 4.0d0*atan(1.d0)
  utail = (8.0d0/3.0d0)*pi*den*(((1.0d0/3.0d0)*((1.0d0/rc)**9)) - ((1.0d0/rc)**3))
  
  ! pick a random atom
  call random_number(q)
  q = q*dble(N)
  k = int(q) + 1
  x_old = x(k)
  y_old = y(k)
  z_old = z(k)

  ! make a move
  call random_number(qx)
  call random_number(qy)
  call random_number(qz)
  x_new = x_old + (2.0d0*qx - 1.0d0)*dr
  y_new = y_old + (2.0d0*qy - 1.0d0)*dr
  z_new = z_old + (2.0d0*qz - 1.0d0)*dr

  ! apply PBC
  if (x_new > L) then
     x_new = x_new - L
  end if

  if (x_new < 0.0d0) then
     x_new = x_new + L
  end if

  if (y_new > L) then
     y_new = y_new - L
  end if

  if (y_new < 0.0d0) then
     y_new = y_new + L
  end if

  if (z_new > L) then
     z_new = z_new - L
  end if

  if (z_new < 0.0d0) then
     z_new = z_new + L
  end if

  ! compute energy difference
  U_old = 0.0d0
  U_new = 0.0d0

  do j = 1, N
     if (j /= k) then

        xj = x(j)
        yj = y(j)
        zj = z(j)

        dx_old = abs(xj - x_old)
        dy_old = abs(yj - y_old)
        dz_old = abs(zj - z_old)
        dx_new = abs(xj - x_new)
        dy_new = abs(yj - y_new)
        dz_new = abs(zj - z_new)

        ! minimum image
        if (dx_old > halfL) then
           dx_old = L - dx_old
        end if

        if (dy_old > halfL) then
           dy_old = L - dy_old
        end if

        if (dz_old > halfL) then
           dz_old = L - dz_old
        end if

        if (dx_new > halfL) then
           dx_new = L - dx_new
        end if

        if (dy_new > halfL) then
           dy_new = L - dy_new
        end if

        if (dz_new > halfL) then
           dz_new = L - dz_new
        end if

        r2_old = dx_old**2 + dy_old**2 + dz_old**2
        r2_new = dx_new**2 + dy_new**2 + dz_new**2

        if (r2_old <= r2c) then
           U_old = U_old + 4.0d0*(((1.0d0/r2_old)**6) - ((1.0d0/r2_old)**3))
        end if

        if (r2_new <= r2c) then
           U_new = U_new + 4.0d0*(((1.0d0/r2_new)**6) - ((1.0d0/r2_new)**3))
        end if

     end if
  end do

  dU = U_new - U_old

  ! accept/reject translational move
  if (dU <= 0.0d0) then
     x(k) = x_new
     y(k) = y_new
     z(k) = z_new
     acc = acc + 1
  end if

  if (dU > 0.0d0) then
     w = exp(-dU/T)
     call random_number(v)

     if (v <= w) then
        x(k) = x_new
        y(k) = y_new
        z(k) = z_new
        acc = acc + 1
     end if
  end if

end subroutine attempt_move


  
     
