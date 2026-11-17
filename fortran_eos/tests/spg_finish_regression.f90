program spg_finish_regression
use compose_eos_types,only:wp
use compose_eos_gdfm_native
use compose_eos_spg_native
implicit none
type(gdfm_native_model_t)::m
type(spg_native_result_t)::cl
type(gdfm_native_beta_t)::un
real(wp)::cross
character(160)::status
cl%evaluated=.true.;cl%volume_fraction=1.0127467909520987e-8_wp
cl%energy_per_baryon=8.0057667249900870_wp
if(.not.spg_native_dilute_uniform_limit(cl,8.0057667249895985_wp))stop 10
cl%volume_fraction=.1_wp
if(spg_native_dilute_uniform_limit(cl,8.0057667249895985_wp))stop 11
cl%volume_fraction=1e-8_wp;cl%energy_per_baryon=8.005_wp
if(spg_native_dilute_uniform_limit(cl,8.0057667249895985_wp))stop 12
call load_gdfm_native('SPG_M5',m,status)
if(status/='ok')stop 2
call spg_native_crossing(m,.9202953869399835_wp,19.49999190059163_wp, &
 [7.3000000000000051e-2_wp,7.3250000000000051e-2_wp],cross,cl,un,status)
write(*,*)trim(status),cross,cl%max_projected_gradient,cl%epsilon-un%epsilon
if(status/='ok')stop 3
if(.not.cl%stationary.or.cl%max_projected_gradient>1e-5_wp)stop 4
if(abs(cl%epsilon-un%epsilon)>1e-8_wp)stop 5
write(*,*)'PASS M5 finite cluster crossing with independent uniform-limit exclusion'
end program
