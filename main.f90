	!###############################################################################
	!     
	!     Simple brush: Standard Molecular Theory Program 
	!    
	!     Calculates a weak polyelectrolyte brush in poor sv conditions 
	!     Calculates free-energy
	!     Solve BINODAL METHOD
	!     MARCH 2020
	!         
	!###############################################################################
!use pks
use system
use const

	implicit none
	integer i, flagcrash
	!print*, 'Program Simple Brush'
	!print*, 'GIT Version: ', _VERSION
yes=0 ! es para  chequear si encuentra o no xalpha, xbeta
flagcrash=1
	call readinput  ! reads input variables from file
	call allocation ! allocates memory
	Keo=10**(-pKeo)
	call solve(flagcrash)
	!call fe(cc, ccc)         ! calculates and saves free energy to disk
	call salvar(flagcrash)
	print*, 'Save OK',yes
	call endall     ! clean up and terminate
end
subroutine endall
stop
end subroutine
