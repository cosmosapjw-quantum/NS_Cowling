program mass_probe
 use ns_precision, only: wp
 use ns_status
 use ns_table_io, only: normalized_table_t,load_normalized_table
 use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos
 use ns_controlled_table_variant, only: controlled_variant_audit_t,build_controlled_variant
 use ns_stellar_methods, only: integrate_cells,cell_run_t
 implicit none
 type(normalized_table_t)::raw,variant
 type(table_eos_t)::model
 type(table_audit_t)::audit
 type(controlled_variant_audit_t)::va
 type(cell_run_t)::run
 type(status_t)::st
 character(len=2048)::path,arg
 real(wp)::pc,p_cut,rtol,step_cap,r0
 integer::ios
 if(command_argument_count()/=6) stop 2
 call get_command_argument(1,path)
 call get_command_argument(2,arg);read(arg,*,iostat=ios)pc;if(ios/=0)stop 2
 call get_command_argument(3,arg);read(arg,*,iostat=ios)p_cut;if(ios/=0)stop 2
 call get_command_argument(4,arg);read(arg,*,iostat=ios)rtol;if(ios/=0)stop 2
 call get_command_argument(5,arg);read(arg,*,iostat=ios)step_cap;if(ios/=0)stop 2
 call get_command_argument(6,arg);read(arg,*,iostat=ios)r0;if(ios/=0)stop 2
 call load_normalized_table(trim(path),raw,st);call must(st)
 call build_controlled_variant(raw,0._wp,69,variant,va,st);call must(st)
 call init_table_eos(variant,model,audit,st);call must(st)
 call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,'dopri',0._wp,0._wp,run,center_order=2)
 call must(run%status)
 write(*,'(a,es25.16e3,a,es25.16e3,a,i0)')'M=',run%mass_kg,' R=',run%radius_m,' cells=',run%cells
contains
 subroutine must(s)
  type(status_t),intent(in)::s
  if(.not.is_ok(s))then
   write(*,'(a,a,a)')trim(status_name(s%code)),' ',trim(s%detail)
   stop 3
  endif
 end subroutine
end program
