
subroutine mu2(x,potquim2)
use const
use system
implicit none
real*8 x(2)
real*8 potquim2,fa_A,fa_B
real*8 xphiA,xphiB,xmphiA,xmphiB,xsolv
real*8 frac(2)

xphiA=x(1)
xphiB=x(2)
xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)
xsolv=(1.0-xphiA-xphiB)
call fracasos(x,frac)
!frac=fracasos(x,xphiB)
fa_A=frac(1)
fa_B=frac(2)
potquim2= log(xmphiA*vs)-(Ma*vp/vs)*log(xsolv)+Ma*log(1.-fa_A)
!prpint* ,fa_A,fa_b,potquim2
!stop
end  subroutine 
