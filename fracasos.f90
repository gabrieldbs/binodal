subroutine fracasos(x,frac)
use system
use const

real*8 frac(2),x(2)
real*8 xphiA,xphiB,xmphiA,xmphiB
real*8 eta,M ,etap
xphiA=x(1);
xphiB=x(2);
xmphiA=xphiA/(Ma*vp);
xmphiB=xphiB/(Mb*vp);

eta=xphiA/xphiB;
etap=xphiB/xphiA;
M=1.0/(xmphiA*Ma*Keo);
frac(1)=0.5*(1.0+etap+M )-( (0.5*(1.0+etap+M))**2-etap)**0.5;
frac(2)=frac(1)*eta;
!fa=y(1);
!fb=y(2);                                                    
end subroutine
