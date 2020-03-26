subroutine fe(x,elib)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine calculates the free energy of the system and chemical potential 
! of chains
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use results
!use chainsdat
use system
!use kai
use bulk
use molecules
use const
implicit none

integer cc, ccc
real*8 x(2)
real*8 elib
real*8 Free_energy, F_Mix_s, F_mix_A,F_mix_B,F_chem_A,F_chem_B,F_par_ion,F_par_ion_const
real*8 F_Mix_neg, F_Mix_Hplus
real*8 Free_energy2, sumpi, sumrho, sumel, sum, mupolA, mupolB, pilat, sumas,diffener,sumex!,hcapa,normhcapa
real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro
integer i, iz, iiz
!real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solventireal
real*8 xphiA,xphiB,xsolv,fa_A,fa_B,xmphiA,xmphiB
real*8 frac(2)

! Calculation of xsolvent
xphiA=x(1)
xphiB=x(2)
xsolv=(1.0-xphiA-xphiB)
 call fracasos(frac)
fa_A=frac(1)
fa_B=frac(2)
xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)

Free_Energy = 0.0

! 1. solvent entropy

  F_Mix_s = (xsolv/vs)*(log(xsolv)-1.0)   

! 2. pol A entropy

  F_mix_A = xmphiA*(log(xmphiA*vs)-1.0)

! 3. pol B e ntropy
 
  F_mix_B = xmphiB*(log(xmphiB*vs)-1.0)
 
Free_Energy = Free_Energy + F_Mix_s + F_mix_A +F_mix_B


! 4. Chemical Equilibrium
  F_chem_A= (xphiA/vp)*((1.-fa_A)*log(1.0-fa_A)+fa_A*log(fa_A))
  F_chem_B= (xphiB/vp)*((1.-fa_B)*log(1.0-fa_B)+fa_B*log(fa_B)) 

Free_Energy = Free_Energy + F_chem_A +F_chem_B

! 5 . Par ion.

  F_par_ion = (xphiA/vp)*fa_A*(log(xphiA*fa_A/vp)-1.0)
! 6. par ion constra.

  F_par_ion_const = fa_A*(xphiA/vp)*log(Keo)

  Free_Energy= Free_Energy-F_par_ion-F_par_ion_const
elib=Free_Energy
return 
end subroutine

