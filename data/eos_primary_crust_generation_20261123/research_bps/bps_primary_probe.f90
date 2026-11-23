program bps_primary_probe
 use compose_eos_types,only:wp,eos_data_t
 use compose_eos_bps_primary
 implicit none
 type(bps_catalogue_t)::cat
 type(bps_state_t)::s,eq
 type(eos_data_t)::eos
 character(1024)::source,out,table,arg,line
 character(256)::status
 real(wp)::rho,p,nb,z,a,gamma,erho,ep,egamma,maxrho,maxp,maxgamma
 integer::n,i,j,k,u,ios,count,mismatch,coexist,dups,phases(12)
 logical::inside
 call get_command_argument(1,source)
 call get_command_argument(2,out)
 call get_command_argument(3,arg)
 read(arg,*,iostat=ios)n
 if(ios/=0)error stop 'probe point count invalid'
 call generate_bps_primary_native(trim(source),n,eos,trim(out),status)
 write(*,'(a)')trim(status)
 if(status/='ok')stop 2
 write(*,'(a,i0)')'ROWS ',eos%n
 phases=0;dups=0
 do i=1,eos%n
  phases(eos%iphase(i)-100)=phases(eos%iphase(i)-100)+1
  if(i>1)then
   if(abs(eos%q(1,i)*eos%nb(i)/(eos%q(1,i-1)*eos%nb(i-1))-1)<1e-12_wp)dups=dups+1
  endif
 enddo
 if(any(phases<3).or.dups/=11)error stop 'missing phase or transition'
 call get_command_argument(4,table)
 if(len_trim(table)==0)stop
 call load_bps_primary_catalogue(trim(source),cat,status)
 if(status/='ok')error stop 'table load failed after generation'
 open(newunit=u,file=trim(table),status='old',iostat=ios)
 if(ios/=0)error stop 'reference file unavailable'
 read(u,'(a)')line
 count=0;mismatch=0;coexist=0;maxrho=0;maxp=0;maxgamma=0
 do
  read(u,*,iostat=ios)rho,p,nb,z,a,gamma
  if(ios<0)exit
  if(ios/=0)error stop 'reference row malformed'
  count=count+1;j=0
  do i=1,cat%n
   if(cat%z(i)==nint(z).and.cat%a(i)==nint(a))j=i
  enddo
  if(j==0)error stop 'reference isotope absent'
  call bps_isotope_state(nb*1e-39_wp,nint(z),nint(a),cat%binding(j),s,status)
  if(status/='ok')error stop 'reference isotope evaluation failed'
  erho=abs(s%epsilon*BPS_RHO_PER_MEVFM3/rho-1)
  ep=abs(s%pressure*BPS_P_CGS_PER_MEVFM3/p-1)
  egamma=abs(s%gamma/gamma-1)
  maxrho=max(maxrho,erho);maxp=max(maxp,ep);maxgamma=max(maxgamma,egamma)
  call bps_state_at_pressure(cat,p/BPS_P_CGS_PER_MEVFM3,eq,status)
  if(status/='ok')error stop 'reference equilibrium evaluation failed'
  if(eq%z/=nint(z).or.eq%a/=nint(a))then
   mismatch=mismatch+1;inside=.false.
   do k=2,eos%n
    if(eos%iphase(k)==eos%iphase(k-1))cycle
    if(abs(eos%q(1,k)*eos%nb(k)/(eos%q(1,k-1)*eos%nb(k-1))-1)>1e-12_wp)cycle
    if(nb*1e-39_wp>=eos%nb(k-1).and.nb*1e-39_wp<=eos%nb(k))then
     if((nint(eos%zbar(k))==nint(z).and.nint(eos%abar(k))==nint(a)).or. &
        (nint(eos%zbar(k-1))==nint(z).and.nint(eos%abar(k-1))==nint(a)))inside=.true.
    endif
   enddo
   if(inside)coexist=coexist+1
  endif
  write(*,'(a,i0,3(1x,es15.7),1x,i0,1x,i0)')'REFERENCE ',count,erho,ep,egamma,eq%z,eq%a
 enddo
 close(u)
 write(*,'(a,3(1x,es20.12),2(1x,i0))')'REFERENCE_MAX ',maxrho,maxp,maxgamma,mismatch,coexist
 if(count/=43)error stop 'reference row count'
 ! Independently predeclared allowances for printed 4-digit source n,P,rho and
 ! modern dimensional constants. They do not change matching/source gates.
 ! Four source pure-phase rows fall within exact Maxwell coexistence density
 ! gaps. This is a representation distinction verified against both limits,
 ! not a tolerance relaxation or an unclassified composition mismatch.
 if(maxrho>1e-3_wp.or.maxp>2e-3_wp.or.maxgamma>2e-3_wp.or.mismatch/=coexist.or.coexist/=4) &
  error stop 'BPS Table 5 independent reference mismatch'
 write(*,'(a)')'BPS_PRIMARY_REFERENCE_VERIFIED'
end program
