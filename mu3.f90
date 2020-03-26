
subroutine mu3(x,potquim3)
use const
use system

real*8 x(2)
real*8 potquim3,fa_A,fa_B
real*8 xphiA,xphiB,xmphiA,xmphiB,xsolv
real*8 frac(2)

xphiA=x(1)
xphiB=x(2)
xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)
xsolv=(1.0-xphiA-xphiB)
call fracasos(x,frac)
fa_A=frac(1)
fa_B=frac(2)
potquim3= log(xmphiB*vs)-(Mb*vp/vs)*log(xsolv)+Mb*log(1.-fa_B)
!return (potquim3)
end subroutine 
