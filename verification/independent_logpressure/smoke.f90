program smoke
use defense_stellar, only: dp,star
implicit none
integer,parameter::n=2000
integer i
real(dp)::p(n),e(n),density,x(7),compact
do i=1,n
   density=exp(log(1.e-12_dp)+(i-1._dp)/(n-1)*log(2.e12_dp))
   p(i)=300*density**2
   e(i)=914*density+p(i)
enddo
call star(n,p,e,100._dp,0._dp,0._dp,.01_dp,x)
compact=x(1)*1.4766250385_dp/x(2)
if(.not.(x(1)>1 .and. x(2)>5 .and. x(3)>0 .and. x(3)<.2_dp))error stop 'stellar result outside expected range'
if(abs(x(4)*compact**5-2*x(3)/3)>1.e-12_dp)error stop 'tidal normalization'
print '(4es18.9)',x(1:4)
end program
