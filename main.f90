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
use pks
use system
use const

	implicit none
	integer i, flagcrash
	!print*, 'Program Simple Brush'
	!print*, 'GIT Version: ', _VERSION
yes=0
flagcrash=1
	call readinput  ! reads input variables from file
	!call init       ! initialize system dependent variables
	call allocation ! allocates memory
	!call kais       ! generates coefficients for poor-sv interactions
	!call creador    ! create chains conformationsi
	Keo=10**(-pKeo)
	call solve(flagcrash)
	!call fe(cc, ccc)         ! calculates and saves free energy to disk
	call salvar(flagcrash)
	print*, 'Save OK',yes
    ! call solve               ! solves the molecular theory
    ! save results to disk
    ! END HERE LOOP
call endall     ! clean up and terminate
end
subroutine endall
stop
end subroutine
