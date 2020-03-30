subroutine checksol(xa,xb,nr)
use system 
use const
 
implicit none
real*8 f1,f2,f3,nr,Penalty,mu3alpha,mu3beta
real*8 elib,xatest(2),xbtest(2)
real*8 xa(2),xb(2),potquim,mu2alpha,mu2beta
real*8 fealpha, febeta
xatest(1)=xa(1)
xatest(2)=xa(2)
xbtest(1)=xb(1)
xbtest(2)=xb(2)
nor=nr
potquim=0
call mu2(xatest,potquim)
mu2alpha=potquim
potquim=0
call mu2(xbtest,potquim)
mu2beta=potquim
potquim=0
call mu3(xatest,potquim)
mu3alpha=potquim
potquim=0
call mu3(xbtest,potquim)
mu3beta=potquim
elib=0
call fe(xatest,elib)
fealpha=elib
elib=0
call fe(xbtest,elib)
febeta=elib
elib=0

Penalty=((xatest(1)-xbtest(1))**2+(xatest(2)-xbtest(2))**2)
!Penalty=1.
f1= (mu2alpha-mu2beta)/Penalty
nor=nor + f1**2
f2= (mu3alpha-mu3beta)/Penalty
nor=nor+f2**2
f3=(fealpha-febeta&
-((xatest(1)-xbtest(1))/(Ma*vp))*(mu2alpha+mu2beta )/2.&
-((xatest(2)-xbtest(2))/(Mb*vp))*(mu3beta +mu3alpha)/2.)/Penalty
nor=nor+f3**2
print*,f1,f2,f3,nor

end subroutine
