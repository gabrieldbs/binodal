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

   x2betaOK=3.731051268530703E-01; !x2beta inicial 
   x3betaOK=3.731051268530703E-01; !x3beta inicial
   phimin=-8; !valor minimo del exponente
   phimax=-2; ! valor maximo del exponente
   nmax=npasos     ! npasos por consola 
   do i=1,nmax
      xiter=phimin+(i-1)*(phimax-phimin)/(nmax-1)
      x2alphafixed= 10**xiter ! x2phialpha

      x1(1)=x2alphafixed           !x3phialpha inicial 
      x1g(1)=x1(1)
      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)

      call call_kinsol(x1, x1, ier)
 
       if((ier.lt.0).or.(norma.gt.error)) then ! failed...
           print*, 'Error in solver: ', ier
           print*, 'norm ', norma
           print*, 'st', st
           print*, 'pH', pHbulk
           !call endall
       endif

     if(infile.ne.-1) then
       if((ier.lt.0).or.(.not.((norma.gt.0).or.(norma.lt.0))).or.(norma.gt.error)) then ! exploto...
         print*, 'solve: Error in solver: ', ier
         print*, 'solve: norma ', norma
         flagcrash = 1
         return
       endif
     endif

   enddo

return
end subroutine

