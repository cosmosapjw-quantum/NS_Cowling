program spg_finish_branch_probe
use compose_eos_types,only:wp
use compose_eos_gdfm_native
use compose_eos_spg_native
implicit none
type(gdfm_native_model_t)::m
type(spg_native_result_t)::cl,out
type(gdfm_native_beta_t)::un
type(spg_native_result_t),allocatable::a(:)
real(wp)::cross,bracket(2),s=.9202953869399835_wp,b=19.49999190059163_wp
integer::count,i
character(160)::status
call load_gdfm_native('SPG_M5',m,status)
call spg_native_transition(m,s,b,cross,cl,un,status,bracket,count)
write(*,*)'TRANSITION',trim(status),cross,cl%nb,bracket,count
call spg_native_solve(m,s,b,cl%nb,out,status,attempts=a)
write(*,*)'SOLVE',trim(status)
call gdfm_native_beta(m,cl%nb,un,status)
write(*,*)'UNIFORM',un%epsilon/cl%nb-gdfm_mass,un%yp
do i=1,size(a)
write(*,*)i,trim(a(i)%status),a(i)%stationary,a(i)%energy_per_baryon,a(i)%max_projected_gradient,a(i)%min_scaled_curvature
write(*,*)'coordinates',a(i)%z,'grad',a(i)%gradient
write(*,*)'phase',a(i)%volume_fraction,a(i)%gas_proton_margin,a(i)%muon_margin
enddo
end program
