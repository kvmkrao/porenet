!**************************************************************
!Author: V Kotteda 
!Date :  Oct 22, 2021 
!**************************************************************
program main 
!implicit real(a-h,o-z) 
implicit none
integer :: nx, ny, nz 
double precision :: maxdis
double precision :: sigma, mu
double precision :: low, high
integer :: nodes, ic, il, i, j,k,ind = 0 

double precision, dimension (:), allocatable :: xc  
double precision, dimension (:), allocatable :: yc  
double precision, dimension (:), allocatable :: zc  
double precision, dimension (:), allocatable :: radn 

double precision :: calavg=0, calsd=0, radlnk, elmlen 
double precision :: r1, r2, rmin, rmax, radius
REAL snorm
double precision :: dx, dy, dz 
integer i1, i2 

open(20,file='input.dat', status='unknown') 
read(20,*) nx, ny, nz 
write(*,*) "nx ny nz", nx, ny, nz 
read(20,*) dx, dy, dz 
write(*,*) "dx dy dz", dx, dy, dz
read(20,*) sigma, mu 
write(*,*) "sigma mu", sigma, mu 
read(20,*) low, high 
write(*,*) "cutoff: low high", low, high
close(20) 
 
nodes = nx * ny * nz  
!allocate(xc(nx)) 
!allocate(yc(ny)) 
!allocate(zc(nz))
 
allocate(xc(nodes))
allocate(yc(nodes))
allocate(zc(nodes))
allocate(radn(nodes))

!call Normal_Distribution(1.0d0,1.0d0, 1000, pdist) 
!           avg, sd 
!http://www.netlib.org/random/  library 
!call gennor(5.0,1.0) 
open(unit=12,file='node_dia.dat', status='unknown') 

ind = 0 
do while(ind < nodes) 
  radius = sigma*snorm() + mu
  if(radius > low .and. radius < high) then 
      ind = ind+1
      radn(ind) = radius 
      write(12,*) ind, radn(ind) 
      calavg = calavg + radn(ind) 
  end if 
end do 


open(unit=10,file='node.dat',status='unknown') 
write(10,*) nodes 

ic = 0 
do i=1,nz 
   do j=1, ny  
      do k=1, nx
         ic = ic + 1 
         zc(ic) = float(i-1)*dz*rand_number()
         yc(ic) = float(j-1)*dy*rand_number() 
         xc(ic) = float(k-1)*dx*rand_number()
         write(*,*) i, j, k , dx, dy, dz 
         !write(10,*) ic, xc(i), yc(j),zc(k)
         write(10,*) xc(ic), yc(ic),zc(ic), radn(ic) 
      end do 
    end do 
end do 


rmin = 1e6 
rmax = 1e-6 
do i=1,nodes 
   if(radn(i) < rmin) rmin = radn(i)
   if(radn(i) > rmax) rmax = radn(i)
end do 

write(*,*) "min, max radius", rmin, rmax 

calavg = calavg/float(nodes) 

do i=1,nodes
  calsd = calsd + (radn(i) - calavg)**2.0 
end do
 
calsd = sqrt(calsd/float(nodes))
write(*,*) "avg, sd", calavg, calsd 

open(unit=11,file='link.dat', status='unknown') 
il=0

do i=1,nodes-1
   il   = il+1
   r1   = radn(i) 
   r2   = radn(i+1)
   i2   = i+1    
   !elmlen = sqrt( (xc(i)-xc(i2))**2.0 +(yc(i)-yc(i2))**2.0+(zc(i)-zc(i2))**2.0)
   elmlen = max(abs (xc(i)-xc(i2)), abs(yc(i)-yc(i2)) , abs(zc(i)-zc(i2)))
   call cal_link_rad(r1, r2, dy, radlnk) 
   write(11,*) i-1, i, radlnk, elmlen
end do 

do i=1,nodes-nx-1
  il = il +1 
  r1 = radn(i) 
  r2 = radn(i+nx)
  i2 = i+nx  
  !elmlen = sqrt( (xc(i)-xc(i2))**2.0 +(yc(i)-yc(i2))**2.0+(zc(i)-zc(i2))**2.0)
  elmlen = max(abs (xc(i)-xc(i2)), abs(yc(i)-yc(i2)) , abs(zc(i)-zc(i2)))
  call cal_link_rad(r1, r2,dx, radlnk) 
  write(11,*)  i-1,i+nx-1, radlnk, elmlen
end do 

do i=1,nodes-nx*ny-1
  r1 = radn(i) 
  r2 = radn(i+nx*ny) 
  i2  = i+nx*ny 
  call cal_link_rad(r1, r2,dz, radlnk) 
  !elmlen = sqrt( (xc(i)-xc(i2))**2.0 +(yc(i)-yc(i2))**2.0+(zc(i)-zc(i2))**2.0)
  elmlen = max(abs (xc(i)-xc(i2)), abs(yc(i)-yc(i2)) , abs(zc(i)-zc(i2)))
  write(11,*)  i-1,i+nx*ny-1, radlnk, elmlen
  il = il + 1 
end do 


!call Normal_Distribution(1.0d0,1.0d0, 1000, pdist) 
!           avg, sd 
!http://www.netlib.org/random/  library 
!call gennor(5.0,1.0) 


!deallocate(xc) 
!deallocate(yc) 
!deallocate(zc) 
!deallocate(radn) 


contains 
subroutine cal_link_rad(radi, radj, maxd, radlk)

double precision, intent(in) :: radi, radj, maxd
double precision :: qir, qjr, radin, radjn
double precision, intent(out) :: radlk
!refer to Joekar-Niasar, JFM, vol 655 (2010) Non-eqilibroum effects in capillarity and interfacial area in two-phase flow: 
double precision, parameter :: en=0.2 
double precision :: PI=4.D0*ATAN(1.D0) 

radin  = radi/maxd 
radjn  = radj/maxd
!write(*,*) radi, radj, radin, radjn, maxd 
qir   = radin*sin(PI/4.0)/(1.0-radin*cos(PI/4.0))**en
qjr   = radjn*sin(PI/4.0)/(1.0-radjn*cos(PI/4.0))**en
radlk = maxd*qir*qjr*( qir**(1.0/en)+qjr**(1.0/en) )**(-en) 

return 
end 

subroutine Normal_Distribution(sigma, mu, samples, x)
  implicit none
 
!!  integer, parameter :: i64 = selected_int_kind(18)
!!  integer, parameter :: r64 = selected_real_kind(15)
!!  integer(i64), parameter :: samples = 1000000_i64
!!  real(r64) :: mean, stddev
!!  real(r64) :: sumn = 0, sumnsq = 0
!!  integer(i64) :: n = 0 
!!  integer(i64) :: bin(-50:50) = 0
!!  integer :: i, ind
!!  real(r64) :: ur1, ur2, nr1, nr2, s

   integer :: samples
   double precision :: mean, stddev
   double precision :: sumn = 0, sumnsq = 0
   integer :: n = 0 
!   integer :: bin(-50:50) = 0
   integer :: i, ind
   double precision :: ur1, ur2, nr1, nr2, s
   double precision :: sigma, mu
   double precision :: x(samples) 

   do n = 1, samples/2
    call random_number(ur1)
    call random_number(ur2)
   
    ur1 = ur1 * 2.0 - sigma
    ur2 = ur2 * 2.0 - sigma
 
    s = ur1*ur1 + ur2*ur2  
    if(s >= (mu+sigma)) cycle
 
    nr1  = ur1 * sqrt(-2.0*log(s)/s)
    x(n) = nr1 
    !ind = floor(10.0*nr1)
    !bin(ind) = bin(ind) + 1
    sumn = sumn + nr1
    sumnsq = sumnsq + nr1*nr1
 
    nr2    = ur2 * sqrt(-2.0*log(s)/s)
    x(n+1) = nr2  
    !ind = floor(10.0*nr2)
    !bin(ind) = bin(ind) + 1
    sumn = sumn + nr2
    sumnsq = sumnsq + nr2*nr2
    !n = n + 2
    write(13,*) n, nr2
    write(13,*) n+1, nr1
    
  end do
 
  mean = sumn / n  ! =0 
  stddev = sqrt(sumnsq/n - mean*mean) ! =1 
 
  write(*, "(a, i0)") "sample size = ", samples
  write(*, "(a, f17.15)") "Mean :   ", mean
  write(*, "(a, f17.15)") "Stddev : ", stddev
 
!  do i = -15, 15 
!    write(*, "(f4.1, a, a)") real(i)/10.0, ": ", repeat("=", int(bin(i)*500/samples))
!  end do
 
return 
end 


FUNCTION rand_number()

IMPLICIT NONE
double precision :: rand_number

  rand_number  =  rand()
  do while (rand_number.le.0.7)
    rand_number = rand()
  end do 
  write(*,*) rand_number
  
end function rand_number

end program main 
