subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use results
!use chainsdat
use solver
use system
use molecules
use bulk
use const
use results
implicit none
 
integer ntot
real*8 x(3*dimz),f(3*dimz)
real*8 protemp, protemp1
integer i,j, k, ix, iy, iz, ii, ax, ay, az, temp, iiZ
real*8 xpotA(dimz),fdisbc,fdisAC
real*8 xpotB(dimz)
!real*16 psi2(0:dimz+1) ! psi plus boundaries at z=0 and dimz+1
!real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solvent
real*8 m, eta, Penality
! Kinsol
integer*4 ier
real*8 vectalpha(2),vectbeta(2)
real*8 x2alpha,x3alpha,x2beta,x3beta
real*8 mu2alpha,mu2beta,mu3alpha,mu3beta,fealpha,febeta
real*8 potquim2,elib,potquim3
! Recovers xh and psi from x
ntot = dimz 

if(VUELTA==1) then

do iz=1,dimz
   x3alpha=x(iz)
   x2beta=x(iz+ntot)	
   x3beta=x(iz+ntot+ntot)	!!!! 
enddo
x2alpha=x2alphafixed

endif

if (VUELTA==2) then
do iz=1,dimz
   x2alpha=x(iz)
   x2beta=x(iz+ntot)    
   x3beta=x(iz+ntot+ntot)       !!!! 
enddo
x3alpha=x3alphafixed

endif

vectalpha(1)=x2alpha
vectalpha(2)=x3alpha
vectbeta(1)=x2beta
vectbeta(2)=x3beta

! Pot quimico respecto de phiA 
Penality=((x2alpha-x2beta)**2+(x3alpha-x3beta)**2)
!Penality=1.
call mu2(vectalpha,potquim2)
mu2alpha = potquim2
!print*,mu2alpha

potquim2=0
call mu2(vectbeta,potquim2)
mu2beta=potquim2
potquim2=0
!print*,mu2beta

call mu3(vectalpha,potquim3)
mu3alpha = potquim3
potquim3=0
!print*,mu3alpha

call mu3(vectbeta,potquim3)
mu3beta=potquim3
potquim3=0
!print*,mu3beta
!stop  OK
 
call fe(vectalpha,elib)
fealpha=elib
elib=0 

!print*,fealpha

call fe(vectbeta,elib)
febeta=elib
elib=0

!print*,febeta
!stop
do iz=1,dimz
 f(iz)= mu2alpha-mu2beta
 f(iz)=f(iz)/Penality
enddo

! Pot quimico respecto de phiB

do iz=1,dimz
  f(iz+ntot)= mu3alpha-mu3beta
 f(iz+ntot)= f(iz+ntot)/Penality
enddo

! Recta tangente
 
do iz=1,dimz
	f(iz+ntot+ntot)= (fealpha-febeta&
-((x2alpha-x2beta)/(Ma*vp))*(mu2alpha+mu2beta )/2.&
-((x3alpha-x3beta)/(Mb*vp))*(mu3beta +mu3alpha)/2.)/Penality !	
enddo

iter = iter + 1
norma = 0.0

do i = 1, 3*ntot
norma = norma +(f(i))**2    
enddo

!if (norma.lt.1E-3) then
!   print*,x2alpha,x3alpha,x2beta,x3beta,norma
!   print*,f(1),f(2),f(3)
!endif
!print*, iter, norma

!print*,x2alpha,x3alpha,x2beta,x3beta
ier = 0.0
   
return
end subroutine
