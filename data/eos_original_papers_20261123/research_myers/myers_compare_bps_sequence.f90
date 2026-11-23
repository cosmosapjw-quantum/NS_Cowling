! Comparison of a new formula model to the rounded BPS published sequence.
! This does not establish agreement with the inaccessible UCRL-11980 mass table.
program myers_compare_bps_sequence
 use compose_eos_types, only: wp
 use compose_eos_bps_primary, only: bps_catalogue_t,load_bps_primary_catalogue
 use compose_eos_myers1966, only: myers_mass_t,myers_mass
 implicit none
 type(bps_catalogue_t)::bps
 type(myers_mass_t)::m14,m20
 character(512)::path,status
 integer::i
 call get_command_argument(1,path)
 call load_bps_primary_catalogue(trim(path),bps,status)
 if(status/='ok')then
  print *,trim(status)
  stop 2
 endif
 print '(a)', 'Z,A,BPS1971_b_MeV_per_baryon,Myers1966_magic14_b,delta14,Myers1966_magic20_b,delta20'
 do i=1,bps%n
  call myers_mass(bps%a(i)-bps%z(i),bps%z(i),14,m14)
  call myers_mass(bps%a(i)-bps%z(i),bps%z(i),20,m20)
  if(.not.m14%valid.or..not.m20%valid)stop 2
  if(i>=4.and.abs(m14%binding-bps%binding(i))>5e-5_wp)stop 3
  print '(2(i0,","),4(es25.16e3,","),es25.16e3)',bps%z(i),bps%a(i),bps%binding(i),m14%binding, &
    m14%binding-bps%binding(i),m20%binding,m20%binding-bps%binding(i)
 enddo
end program
