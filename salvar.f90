subroutine salvar(flagcrash)
use results

integer i,j
open (unit=1,file='alphaarray.txt',status='replace')

do i=1,conteo
   write (1,*) arrayalpha(1,i), arrayalpha(2,i)
end do

open (unit=1,file='betaarray.txt',status='replace')

do j=1,conteo
   write (1,*) arraybeta(1,j), arraybeta(2,j)
end do

end subroutine

