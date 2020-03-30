subroutine solve(flagcrash)

use pks 
use system
use results
use solver
implicit none
real*8 x2betaOK,x2alphaOK
real*8 x3betaOK,x3alphaOK
real*8 phimin
real*8 phimax,criterio
integer nmax
real*8 xiter
real*8 x1(3*dimz)
real*8 x1g(3*dimz)
integer ier, i,newt
integer flagcrash
real*8 potquim2,potquim3,mu2alpha,mu2beta
real*8 xatest(2),xbtest(2),segsolalpha(2),segsolbeta(2)

conteo=0  !lo uso  para  guardar en el arrayalpha y arraybeta
criterio=1E-3!criterio pra la norma
!xo=	inicial!
!    x2betaOK=3.731051268530703E-01; !x2beta inicial 
!    x3betaOK=3.731051268530703E-01; !x3beta inicial
 !  x2betaOK= 3.9385204351639801E-02 !MA=10MB=10
 !  x3betaOK=7.9694435360414845E-02 !Ma=10Mb=10
 x2betaOK= 4.7669872426669833E-002 
 x3betaOK=  8.3887330029032467E-002  
!   x2betaOK=  1.7572211418669388E-002  !Ma100Mb=100
!   x3betaOK=   8.7546403827123531E-002!Mb=100Ma=100

   phimin=-12; !valor minimo del exponente
 phimax=-0.5! valor maximo del exponente
   nmax=npasos     ! npasos por consola i
   VUELTA=1 ! la idea es usarlo para que identifique si está fijo x2alpha o x3alpha
x3alphaOK=0! 
   do i=1,nmax
      xiter=phimin+(i-1)*(phimax-phimin)/(nmax-1)
      x2alphafixed= 10**xiter ! x2phialpha

      if (x3alphaOK.eq.0)then !si es 0   toma el inicial que habiamos seteado para matlab
      x1(1)=7.0569559047113746E-003  !x3phialpha inicial Ma=10=MB
    !  x1(1)=x2alphafixed
 !     x1(1)= 5.0187523024411020E-002 !MA=100=Mb
       x1(1)= 2.4791273427165763E-003   !x3phialpha inicial Ma=10=MB 10-12
      x1g(1)=x1(1)
      endif
      if (x3alphaOK.ne.0)then !cuando encuentre  una solucion la idea es que lo reutilice 
      x1(1)=x3alphaOk
      x1g(1)=x1(1)
      endif

      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)

      call call_kinsol(x1, x1g, ier)

       if (norma.lt.criterio) then !esto es para saber si encontró o no solución 
         conteo=conteo+1
         yes=yes+1
         print*,'Yes',yes
         print*,'x2alpha,x3alpha,x2beta,x3beta',x2alphafixed,x1(1),x1(2),x1(3) 
         x3alphaOK= x1(1)
         x2betaOK = x1(2)
         x3betaOK = x1(3)
         arrayalpha(1,conteo)=x2alphafixed
         arrayalpha(2,conteo)=x1(1)
         arraybeta(1,conteo)=x1(2)
         arraybeta(2,conteo)=x1(3)

!         segunda solución  intercambiando x2 con x3 chequeo si vale  en el caso de Ma=Mb deberia serlo
         xatest(1)=x1(1)
         xatest(2)=x2alphafixed
         xbtest(1)=x1(3)
         xbtest(2)=x1(2)
         nor=0
         newt=0
         call checksol(xatest,xbtest,nor)
         print*,nor
         if(nor.lt.criterio)then
         newt=newt+1
yes =yes+1
          conteo=conteo+1
          arrayalpha(1,conteo)=xatest(1)
          arrayalpha(2,conteo)=xatest(2)
          arraybeta(1,conteo)=xbtest(1)
          arraybeta(2,conteo)=xbtest(2)
         endif
       !  potquim2=0 
       !  call mu2(xatest,potquim2)
       !  mu2alpha=potquim2
       !  potquim2=0
       !  call mu2(xbtest,potquim2)
       !  mu2beta=potquim2
       !  potquim2=0
       !  print*,'mu2alpha',mu2alpha
       !  print*,'mu2beta', mu2beta,ier
       !  stop

      !  Tercer solución  intercambiando x2 con x3 chequeo si vale  en el caso de Ma=Mb deberia serlo
         xatest(1)=x1(2)
         xatest(2)=x1(3)
         xbtest(1)=x2alphafixed
         xbtest(2)=x1(1)
         nor=0
         newt=0
         call checksol(xatest,xbtest,nor)
         print*,nor
         if(nor.lt.criterio)then
         newt=newt+1
yes =yes+1
          conteo=conteo+1
          arrayalpha(1,conteo)=xatest(1)
          arrayalpha(2,conteo)=xatest(2)
          arraybeta(1,conteo)=xbtest(1)
          arraybeta(2,conteo)=xbtest(2)
         endif
       !  potquim2=0 
       !  call mu2(xatest,potquim2)
       !  mu2alpha=potquim2
       !  potquim2=0
       !  call mu2(xbtest,potquim2)
       !  mu2beta=potquim2
       !  potquim2=0
       !  print*,'mu2alpha',mu2alpha
       !  print*,'mu2beta', mu2beta,ier
       !  stop


   !  cuarta  solución  intercambiando x2 con x3 chequeo si vale  en el caso de Ma=Mb deberia serlo
         xatest(1)=x1(3)
         xatest(2)=x1(2)
         xbtest(1)=x1(1)
         xbtest(2)=x2alphafixed
         nor=0
         newt=0
         call checksol(xatest,xbtest,nor)
         print*,nor
         if(nor.lt.criterio)then
         newt=newt+1
yes=yes+1
          conteo=conteo+1
          arrayalpha(1,conteo)=xatest(1)
          arrayalpha(2,conteo)=xatest(2)
          arraybeta(1,conteo)=xbtest(1)
          arraybeta(2,conteo)=xbtest(2)
         endif
    !   !  potquim2=0 
       !  call mu2(xatest,potquim2)
       !  mu2alpha=potquim2
       !  potquim2=0
       !  call mu2(xbtest,potquim2)
       !  mu2beta=potquim2
       !  potquim2=0
       !  print*,'mu2alpha',mu2alpha
       !  print*,'mu2beta', mu2beta,ier
       !  stop


       endif

   enddo
!    VUELTA=2  !  la idea es que en la segunda vuelta recorra en x3 
!   x2betaOK= 3.9385204351639801E-002
!   x3betaOK=7.9694435360414845E-002

!   x2betaOK=3.731051268530703E-01; !x2beta inicial 
!   x3betaOK=3.731051268530703E-01; !x3beta inicial
!x2alphaOK=0
!   do i=1,nmax
!      xiter=phimax-(i-1)*(phimax-phimin)/(nmax-1)
!       xiter=phimin+(i-1)*(phimax-phimin)/(nmax-1)
!       x2alphafixed= 10**xiter ! x2phialpha
!      x3alphafixed= 10**xiter ! x2phialpha

!      if(x2alphaOK.eq.0)then
!        x1(1)=x3alphaOK !x2phialpha inicial 
     !   x1(1)=7.0569559047113746E-003  !x3phialpha inicial 
!        x1g(1)=x1(1)
!      endif
!      if(x2alphaOK.ne.0)then
!        x1(1)=x2alphaOK !x2phialpha inicial 
!       x1(1)=x3alphaOK !x2phialpha inicial 
 
!       x1g(1)=x1(1)
!      endif
!
!      x1(2)=x2betaOK     !x2phibeta inicial
!      x1g(2)=x1(2)
!      x1(3)=x3betaOK     !x3phibeta inicial
!      x1g(3)=x1(3)

!      call call_kinsol(x1, x1g, ier)

!       if (norma.lt.1E-4) then
!         conteo=conteo+1
!         yes=yes+1
!         print*,'Yes',yes
!         print*,'x2alpha,x3alpha,x2beta,x3beta',x1(1),x3alphafixed,x1(2),x1(3)
 !        xatest(1)=x1(1)
!         xatest(2)=x3alphafixed
 !        xbtest(1)=x1(2)
  !       xbtest(2)=x1(3)
         !potquim2=0
        ! call mu2(xatest,potquim2)
        ! mu2alpha=potquim2
        ! potquim2=0
        ! call mu2(xbtest,potquim2)
        ! mu2beta=potquim2
        ! potquim2=0
      !   print*,'mu2alpha',mu2alpha
      !   print*,'mu2beta', mu2beta,ier
       ! stop
!
!         x2alphaOK= x1(1)
!         x2betaOK = x1(2)
!         x3betaOK = x1(3)
!          arrayalpha(1,conteo)=x1(1)
!          arrayalpha(2,conteo)=x3alphafixed
!          arraybeta(1,conteo)=x1(2)
!          arraybeta(2,conteo)=x1(3)

!         xatest(1)=x1(1)
!         xatest(2)=x2alphafixed
!         xbtest(1)=x1(3)
!         xbtest(2)=x1(2)
!         nor=0
!         newt=0
!         call checksol(xatest,xbtest,nor)
 !        print*,nor
 !        if(nor.lt.criterio)then
 !        newt=newt+1
  !        conteo=conteo+1
  !        arrayalpha(1,conteo)=xatest(1)
  !        arrayalpha(2,conteo)=xatest(2)
  !        arraybeta(1,conteo)=xbtest(1)
  !        arraybeta(2,conteo)=xbtest(2)
 !        endif


  !   endif

 ! enddo

return
end subroutine

