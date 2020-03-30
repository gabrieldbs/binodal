module system
real*8 delta   ! delta is the discretization lenght in z direction
integer  dimz  ! number of lattice sites in z direction
real*8 sigmaA
real*8 Ma
integer yes
real*8 Mb
integer npasos
real*8 sigmaB
real*8 csalt
real*8 nor
!real*16 pKaA
!real*16 pkaB
!real*16 pKaANa
!real*16 pkaBCl
real*16 pkEo
real*16 Keo
real*16 pHbulk
!real*8 st
integer VUELTA
endmodule

module results
real*8 x2alphafixed
real*8 x3alphafixed
real*8 arrayalpha(2,60000)
real*8 arraybeta(2,60000)
integer conteo
endmodule

module const
real*8, parameter :: pi = 3.14159 ! pi 
real*8, parameter :: Na = 6.02d23 ! Avogadro's number
real*8, parameter :: lb = 0.714   ! bjerrum lenght in water in nm
real*8, parameter :: vs = 0.03  ! bjerrum lenght in water in nm
real*8, parameter :: vp = 0.11   ! bjerrum lenght in water in nm

real*8 constq
real*8 pKw
endmodule


module solver
real*8, allocatable :: xflag(:)
integer infile
real*8, parameter :: error = 1.0d-6 ! maximum kinsol norm
real*8 norma
integer iter
endmodule


module molecules
real*8 zpos, zneg, zpolA, zpolB ! charges of cation, anions and polyelectrolyte segment
real*8 vsalt, vpol   ! volume of salt and polyelectrolyte segments in units of vsol
real*8 vsol             ! solvent volume 
real*16 K0A, K0B ,K0ANa,K0BCl, K0Eo!K0
endmodule

module bulk
real*16 xHplusbulk, xOHminbulk ! bulk volume fraction of H+ and OH-
real*16 xposbulk, xnegbulk     ! bulk volume fraction of cation and anion
real*16 xsolbulk               ! bulk volume fraction of solvent
real*16 expmupos, expmuneg, expmuHplus, expmuOHmin  ! exp(-beta*mu)*(bulk volume fraction), where mu is the chemical potential
endmodule

module rand
integer seed
endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule
