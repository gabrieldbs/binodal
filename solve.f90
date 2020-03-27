subroutine solve(flagcrash)

use pks 
use system
use results
use solver
implicit none
real*16 x2betaOK
real*16 x3betaOK
real*8 phimin
real*8 phimax
integer nmax
real*16 xiter
real*8 x1(3*dimz)
real*8 x1g(3*dimz)
integer ier, i
integer flagcrash
real*8 potquim2,potquim3,mu2alpha,mu2beta
real*8 xatest(2),xbtest(2)
conteo=0
   x2betaOK=3.731051268530703E-01; !x2beta inicial 
   x3betaOK=3.731051268530703E-01; !x3beta inicial
   phimin=-8; !valor minimo del exponente
   phimax=-2; ! valor maximo del exponente
   nmax=npasos     ! npasos por consola i
VUELTA=1
   do i=1,nmax
      xiter=phimin+(i-1)*(phimax-phimin)/(nmax-1)
      x2alphafixed= 10**xiter ! x2phialpha

      x1(1)=x2alphafixed !x3phialpha inicial 
      x1g(1)=x1(1)
      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)

      call call_kinsol(x1, x1g, ier)

       if (norma.lt.1E-3) then
         conteo=conteo+1
         yes=yes+1
         print*,'Yes',yes
         print*,'x2alpha,x3alpha,x2beta,x3beta',x2alphafixed,x1(1),x1(2),x1(3)
         xatest(1)=x2alphafixed
         xatest(2)=x1(1)
         xbtest(1)=x1(2)
         xbtest(2)=x1(3)
         potquim2=0 
         call mu2(xatest,potquim2)
         mu2alpha=potquim2
         potquim2=0
         call mu2(xbtest,potquim2)
         mu2beta=potquim2
         potquim2=0
 !        print*,'mu2alpha',mu2alpha
!         print*,'mu2beta', mu2beta,ier
       ! stop
         
         x2betaOK = x1(2)
         x3betaOK = x1(3)
         arrayalpha(1,conteo)=x2alphafixed
         arrayalpha(2,conteo)=x1(1)
         arraybeta(1,conteo)=x1(2)
         arraybeta(2,conteo)=x1(3)

       endif
 
      ! if((ier.lt.0).or.(norma.gt.error)) then ! failed...
      !     print*, 'Error in solver: ', ier
      !     print*, 'norm ', norma
       !    print*, 'st', st
        !   print*, 'pH', pHbulk
           !call endall
       !endif
!
!     if(infile.ne.-1) then
!       if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! exploto...
!         print*, 'solve: Error in solver: ', ier
!         print*, 'solve: norma ', norma
!         flagcrash = 1
!         return
!       endif
!     endif

   enddo
    VUELTA=2
   x2betaOK=3.731051268530703E-01; !x2beta inicial 
   x3betaOK=3.731051268530703E-01; !x3beta inicial

   do i=1,nmax
      xiter=phimin+(i-1)*(phimax-phimin)/(nmax-1)
      x3alphafixed= 10**xiter ! x2phialpha

      x1(1)=x3alphafixed !x3phialpha inicial 
      x1g(1)=x1(1)
      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)

      call call_kinsol(x1, x1g, ier)

       if (norma.lt.1E-4) then
         conteo=conteo+1
         yes=yes+1
         print*,'Yes',yes
         print*,'x2alpha,x3alpha,x2beta,x3beta',x1(1),x3alphafixed,x1(2),x1(3)
         xatest(1)=x1(1)
         xatest(2)=x3alphafixed
         xbtest(1)=x1(2)
         xbtest(2)=x1(3)
         potquim2=0
         call mu2(xatest,potquim2)
         mu2alpha=potquim2
         potquim2=0
         call mu2(xbtest,potquim2)
         mu2beta=potquim2
         potquim2=0
      !   print*,'mu2alpha',mu2alpha
      !   print*,'mu2beta', mu2beta,ier
       ! stop


         x2betaOK = x1(2)
         x3betaOK = x1(3)
         arrayalpha(1,conteo)=x1(1)
         arrayalpha(2,conteo)=x3alphafixed
         arraybeta(1,conteo)=x1(2)
         arraybeta(2,conteo)=x1(3)

       endif

   enddo

return
end subroutine

