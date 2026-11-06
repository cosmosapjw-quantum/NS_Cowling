program cell_checks
 use ns_precision, only: wp
 use ns_status
 use ns_eos_types
 use ns_matching_types, only: mu_bridge_t
 use ns_mu_bridge, only: construct_bridge_endpoints
 use ns_cell_eos, only: cell_state,cell_pressure_coordinate
 implicit none
 type(mu_bridge_t)::b
 type(eos_state_t)::s
 type(status_t)::st
 real(wp)::t,dn,p,h
 call construct_bridge_endpoints(1._wp,1._wp,2._wp,2._wp,4._wp,4._wp,b,st)
 call require(is_ok(st),'analytic construction')
 call cell_state(b,0.371_wp,s,dn,st)
 call require(is_ok(st),'analytic evaluation')
 call require(abs(s%pressure_pa-1.371_wp**2)<1e-14_wp,'analytic P')
 call require(abs(s%energy_j_m3-s%pressure_pa)<1e-14_wp,'analytic energy')
 call require(abs(s%chi-1._wp)<1e-14_wp,'analytic chi')
 call cell_pressure_coordinate(b,s%pressure_pa,t,st)
 call require(is_ok(st).and.abs(t-0.371_wp)<1e-14_wp,'local inverse')
 call cell_state(b,-0.01_wp,s,dn,st)
 call require(.not.is_ok(st),'out of domain')
 call cell_pressure_coordinate(b,0.5_wp,t,st)
 call require(.not.is_ok(st),'pressure out of domain')
 h=(1e-10_wp+1e-20_wp)-1e-10_wp
 call construct_bridge_endpoints(1e-10_wp,1e20_wp,1e44_wp,1e-10_wp+h, &
     1e20_wp+1.5_wp*h*1e44_wp,2e44_wp,b,st)
 call require(is_ok(st),'thin-cell construction')
 p=1e20_wp+b%a1*1e-8_wp+b%a2*1e-16_wp+b%a3*1e-24_wp
 call cell_pressure_coordinate(b,p,t,st)
 call require(is_ok(st),'thin-cell inverse status')
 call cell_state(b,t,s,dn,st)
 call require(abs(s%pressure_pa/p-1._wp)<5e-14_wp,'thin-cell pressure accuracy')
 call require(abs(t-1e-8_wp)<1e-19_wp,'thin-cell local coordinate')
 print *, 'cell_checks: 12 checks passed'
contains
 subroutine require(ok,label)
 logical,intent(in)::ok
 character(len=*),intent(in)::label
 if(.not.ok)then
  print *, 'FAIL ',label
  stop 1
 endif
 end subroutine
end program
