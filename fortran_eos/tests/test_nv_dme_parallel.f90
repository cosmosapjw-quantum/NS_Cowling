program test_nv_dme_parallel
 use compose_eos_types,only:wp
 use compose_eos_nv_original,only:nv_original_source_t,load_nv_original_source
 use compose_eos_nv_dme,only:nv_dme_coeff_t,nv_dme_coefficients,nv_dme_local
 use omp_lib,only:omp_get_max_threads
 implicit none
 integer,parameter::samples=1000
 type(nv_original_source_t)::src
 real(wp)::serial(8,samples),parallel(8,samples)
 character(1024)::path,status
 integer::i,failures
 call get_command_argument(1,path)
 call load_nv_original_source(trim(path),src,status)
 if(status/='ok')then;print*,trim(status);error stop 'source load';endif
 do i=1,samples
  call evaluate(i,serial(:,i),status)
  if(status/='ok')then;print*,trim(status);error stop 'serial functional evaluation';endif
 enddo
 failures=0
 !$omp parallel do schedule(dynamic,7) private(status) reduction(+:failures)
 do i=1,samples
  call evaluate(i,parallel(:,i),status)
  if(status/='ok')failures=failures+1
 enddo
 !$omp end parallel do
 if(failures/=0)error stop 'parallel functional errors'
 if(any(serial/=parallel))error stop 'parallel/serial functional discrepancy'
 print '(a,i0)','max_openmp_threads=',omp_get_max_threads()
 print '(a,i0)','points_with_alternating_modes=',samples
 print '(a)','NV_DME_PARALLEL_BITWISE_IDENTITY_PASS'
contains
 subroutine evaluate(i,out,status)
  integer,intent(in)::i
  real(wp),intent(out)::out(8)
  character(*),intent(out)::status
  type(nv_dme_coeff_t)::co
  real(wp)::n,p
  integer::mode
  n=.001_wp+.13_wp*real(mod(i,71),wp)/70
  p=.00001_wp+.09_wp*real(mod(i,53),wp)/52
  mode=mod(i,2)
  call nv_dme_coefficients(src,n,p,co,status,need_hessian=.false.,coefficient_mode=mode)
  if(status/='ok')return
  call nv_dme_local(src,n,p,1.07_wp*co%tau_uniform,[-.001_wp,-.0003_wp],.00002_wp, &
    out(1),out(2:3),out(4:5),out(6:7),out(8),status,coefficient_mode=mode)
 end subroutine
end program
