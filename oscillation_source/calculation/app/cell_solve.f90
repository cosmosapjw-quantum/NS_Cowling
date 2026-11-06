program cell_solve
 use ns_precision, only: wp
 use ns_status
 use ns_table_io, only: normalized_table_t,load_normalized_table
 use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos
 use ns_controlled_table_variant, only: controlled_variant_audit_t,build_controlled_variant
 use ns_stellar_methods, only: integrate_cells, match_cells, cell_run_t
 use ns_numerics_types, only: root_options_t,root_result_t
 use ns_bisection, only: bisect
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 type(normalized_table_t)::raw,variant
 type(table_eos_t)::model
 type(table_audit_t)::audit
 type(controlled_variant_audit_t)::vaudit
 type(status_t)::st
 type(root_options_t)::opts
 type(root_result_t)::root
 type(cell_run_t)::star,shot,checkstar
 character(len=2048)::input_path,output_prefix,table_path
 character(len=16)::task='spectrum',method='dopri',shooting='matched'
 integer::match_index=72,split_index=69,unit,ios,j,nmode=2,center_order=0
 real(wp)::alpha=0.1_wp,pc_pa=1.5e34_wp,p_cut_pa=1e20_wp,r0_m=1._wp
 real(wp)::target_mass_kg=2.644606257493485e30_wp,pc_lo_pa=1.2e34_wp,pc_hi_pa=1.8e34_wp
 real(wp)::rtol=1e-10_wp,mode_rtol=1e-10_wp,step_cap=0.05_wp,mass_rtol=1e-12_wp
 real(wp)::frequency_lo(2)=[2000._wp,6500._wp],frequency_hi(2)=[3000._wp,8000._wp]
 real(wp)::frequency_hz=2500._wp,root_atol_hz=1e-11_wp,residual_tol=1e-8_wp
 real(wp)::match_r(2),match_m(2),boundary_res(2),dummy_r,dummy_m
 real(wp)::f(2),res(2),mshot(2),rshot(2),width(2),lo,hi,mid,ml,mm,mh
 integer::steps(2),mcells(2)
 character(len=2048)::path
 namelist /precision/ table_path,task,method,shooting,match_index,split_index,alpha,pc_pa,p_cut_pa,r0_m, &
 target_mass_kg,pc_lo_pa,pc_hi_pa,rtol,mode_rtol,step_cap,mass_rtol,frequency_lo,frequency_hi, &
 frequency_hz,root_atol_hz,residual_tol,nmode,center_order
 if(command_argument_count()/=2)then
  print *, 'Usage: cell_solve input.nml output_prefix';stop 2
 endif
 call get_command_argument(1,input_path);call get_command_argument(2,output_prefix)
 table_path=''
 open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
 if(ios/=0)stop 2
 read(unit,nml=precision,iostat=ios);close(unit)
 if(ios/=0)then;print *,'invalid precision namelist';stop 2;endif
 if(.not.all(ieee_is_finite([alpha,pc_pa,p_cut_pa,r0_m,target_mass_kg,pc_lo_pa,pc_hi_pa, &
       rtol,mode_rtol,step_cap,mass_rtol,frequency_lo,frequency_hi,frequency_hz,root_atol_hz,residual_tol])))stop 2
 if(center_order/=0.and.center_order/=2)stop 2
 if(target_mass_kg<0._wp.or.root_atol_hz<=0._wp.or.residual_tol<=0._wp)stop 2
 call load_normalized_table(trim(table_path),raw,st);call must(st)
 call build_controlled_variant(raw,alpha,split_index,variant,vaudit,st);call must(st)
 call init_table_eos(variant,model,audit,st);call must(st)
 if(nmode<1.or.nmode>2.or.rtol<=0._wp.or.mode_rtol<=0._wp.or.mass_rtol<=0._wp)stop 2
 if(target_mass_kg>0._wp)then
  if(pc_lo_pa<=0._wp.or.pc_hi_pa<=pc_lo_pa)stop 2
  lo=log(pc_lo_pa);hi=log(pc_hi_pa);mid=(lo+hi)/2._wp
  call mass_residual(lo,ml,st);call must(st)
  call mass_residual(mid,mm,st);call must(st)
  call mass_residual(hi,mh,st);call must(st)
  if(ml>=mm.or.mm>=mh)then;print *,'nonmonotone sampled mass branch';stop 3;endif
  opts%xtol=2e-14_wp;opts%rtol=0._wp;opts%ftol=mass_rtol*0.01_wp;opts%max_iter=90
  call bisect(mass_residual,lo,hi,opts,root);call must(root%status)
  pc_pa=exp(root%root)
 endif
 call integrate_cells(model,pc_pa,p_cut_pa,r0_m,rtol,step_cap,method,0._wp,0._wp,star,center_order=center_order)
 call must(star%status)
 if(target_mass_kg>0._wp)then
  if(abs(star%mass_kg/target_mass_kg-1._wp)>mass_rtol)then
   print *, 'mass residual did not close',star%mass_kg/target_mass_kg-1._wp;stop 3
  endif
 endif
 call integrate_cells(model,pc_pa,p_cut_pa,r0_m,rtol,step_cap,method,0._wp,star%nu_surface, &
      checkstar,trim(output_prefix)//'_background.csv',center_order=center_order)
 call must(checkstar%status)
 match_r=0._wp;match_m=0._wp;boundary_res=0._wp;f=0._wp;res=0._wp;mshot=0._wp;rshot=0._wp;width=0._wp;steps=0;mcells=0
 if(task=='spectrum')then
  do j=1,nmode
   opts%xtol=root_atol_hz;opts%rtol=0._wp;opts%ftol=0._wp;opts%max_iter=100
   call bisect(mode_residual,frequency_lo(j),frequency_hi(j),opts,root);call must(root%status)
   f(j)=root%root;width(j)=root%hi-root%lo
   write(path,'(a,a,i0)')trim(output_prefix),'_mode',j
   if(shooting=='matched')then
    call match_cells(model,pc_pa,p_cut_pa,r0_m,mode_rtol,step_cap,method,f(j),star,match_index, &
         res(j),match_r(j),match_m(j),st,trim(path),steps(j),boundary_res(j),center_order=center_order)
    call must(st)
   else if(shooting=='outward')then
   path=trim(path)//'.csv'
   call integrate_cells(model,pc_pa,p_cut_pa,r0_m,mode_rtol,step_cap,method,f(j),star%nu_surface,shot,trim(path),center_order=center_order)
   call must(shot%status)
   res(j)=shot%residual;mshot(j)=shot%mass_kg;rshot(j)=shot%radius_m;steps(j)=shot%steps;mcells(j)=shot%cells
   else
    print *,'unknown shooting method';stop 2
   endif
   if(abs(res(j))>residual_tol)then
    print *,'mode root coordinate closed but residual failed',j,f(j),res(j)
    ! Still write complete numerical diagnostics below; overall status is failed.
   endif
  enddo
 else if(task=='shot')then
  if(shooting=='matched')then
   call match_cells(model,pc_pa,p_cut_pa,r0_m,mode_rtol,step_cap,method,frequency_hz,star,match_index, &
      res(1),match_r(1),match_m(1),st,trim(output_prefix)//'_shot',steps(1),boundary_res(1),center_order=center_order)
   call must(st);f(1)=frequency_hz
  else
  call integrate_cells(model,pc_pa,p_cut_pa,r0_m,mode_rtol,step_cap,method,frequency_hz,star%nu_surface, &
    shot,trim(output_prefix)//'_shot.csv',center_order=center_order)
  call must(shot%status);f(1)=frequency_hz;res(1)=shot%residual;mshot(1)=shot%mass_kg
  endif
 else if(task/='star')then
  print *,'unknown task';stop 2
 endif
 open(newunit=unit,file=trim(output_prefix)//'.json',status='replace',action='write')
 write(unit,'(a)')'{'
 write(unit,'(a,a,a)')'"method":"',trim(method),'",'
 write(unit,'(a,a,a)')'"shooting":"',trim(shooting),'",'
 write(unit,'(a,i0,a)')'"match_index":',match_index,','
 write(unit,'(a,i0,a)')'"center_order":',center_order,','
 write(unit,'(a,i0,a)')'"split_index":',split_index,','
 if(shooting=='matched')then
  write(unit,'(a)')'"residual_definition":"normalized_matching_determinant",'
 else
  write(unit,'(a)')'"residual_definition":"normalized_outward_surface_Q",'
 endif
 write(unit,'(a,a,a)')'"task":"',trim(task),'",'
 if(vaudit%core_function_unchanged)then
  write(unit,'(a)')'"core_function_unchanged": true,'
 else
  write(unit,'(a)')'"core_function_unchanged": false,'
 endif
 call scalar('max_nodal_pressure_fraction',vaudit%max_low_density_pressure_rel_change)
 call scalar('max_nodal_density_fraction',vaudit%max_low_density_density_rel_change)
 write(unit,'(a,i0,a)')'"admitted_eos_intervals":',audit%admitted_intervals,','
 call scalar('alpha',alpha);call scalar('rtol',rtol);call scalar('mode_rtol',mode_rtol)
 call scalar('step_cap',step_cap);call scalar('pc_pa',pc_pa);call scalar('p_cut_pa',p_cut_pa)
 call scalar('r0_m',r0_m);call scalar('target_mass_kg',target_mass_kg)
 call scalar('mass_kg',star%mass_kg);call scalar('radius_m',star%radius_m)
 call scalar('nu_surface',star%nu_surface);call scalar('surface_pressure_rel_error',star%surface_pressure_rel_error)
 call scalar('f_hz',f(1));call scalar('p1_hz',f(2));call scalar('f_residual',res(1));call scalar('p1_residual',res(2))
 if(shooting=='outward')then
  call scalar('f_mass_kg',mshot(1));call scalar('p1_mass_kg',mshot(2))
  call scalar('f_radius_m',rshot(1));call scalar('p1_radius_m',rshot(2))
 else
  call scalar('f_surface_residual',boundary_res(1));call scalar('p1_surface_residual',boundary_res(2))
 endif
 call scalar('f_match_radius_mismatch_m',match_r(1));call scalar('p1_match_radius_mismatch_m',match_r(2))
 call scalar('f_match_mass_mismatch_kg',match_m(1));call scalar('p1_match_mass_mismatch_kg',match_m(2))
 call scalar('f_bracket_hz',width(1));call scalar('p1_bracket_hz',width(2))
 write(unit,'(a,i0,a)')'"structure_steps":',star%steps,','
 write(unit,'(a,i0,a)')'"native_cells_traversed":',star%cells,','
 write(unit,'(a,i0,a)')'"f_steps":',steps(1),','
 write(unit,'(a,i0,a)')'"p1_steps":',steps(2),','
 if(task=='spectrum'.and.any(abs(res(:nmode))>residual_tol))then
  write(unit,'(a)')'"status":"residual_failed"}'
  close(unit);stop 4
 else
  write(unit,'(a)')'"status":"ok"}'
 endif
 close(unit)
 print *,'complete ',trim(output_prefix),f
contains
 subroutine must(status)
  type(status_t),intent(in)::status
  if(.not.is_ok(status))then
   print *,trim(status_name(status%code)),trim(status%category),trim(status%detail);stop 3
  endif
 end subroutine
 subroutine scalar(name,value)
  character(len=*),intent(in)::name
  real(wp),intent(in)::value
  write(unit,'(a,a,a,es25.16e3,a)')'"',name,'":',value,','
 end subroutine
 subroutine mass_residual(x,value,status)
  real(wp),intent(in)::x
  real(wp),intent(out)::value
  type(status_t),intent(out)::status
  type(cell_run_t)::run
  call integrate_cells(model,exp(x),p_cut_pa,r0_m,rtol,step_cap,method,0._wp,0._wp,run,center_order=center_order)
  status=run%status;value=0._wp
  if(is_ok(status))value=run%mass_kg/target_mass_kg-1._wp
 end subroutine
 subroutine mode_residual(freq,value,status)
  real(wp),intent(in)::freq
  real(wp),intent(out)::value
  type(status_t),intent(out)::status
  type(cell_run_t)::run
  if(shooting=='matched')then
   call match_cells(model,pc_pa,p_cut_pa,r0_m,mode_rtol,step_cap,method,freq,star,match_index, &
        value,dummy_r,dummy_m,status,center_order=center_order)
  else
   call integrate_cells(model,pc_pa,p_cut_pa,r0_m,mode_rtol,step_cap,method,freq,star%nu_surface,run,center_order=center_order)
   status=run%status;value=0._wp
   if(is_ok(status))value=run%residual
  endif
 end subroutine
end program
