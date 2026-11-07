program audit_variant
 use ns_precision, only: wp
 use ns_status
 use ns_table_io, only: normalized_table_t,load_normalized_table
 use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos
 use ns_controlled_table_variant, only: controlled_variant_audit_t,build_controlled_variant
 implicit none
 type(normalized_table_t)::raw,variant
 type(table_eos_t)::model
 type(table_audit_t)::audit
 type(controlled_variant_audit_t)::va
 type(status_t)::st
 character(len=2048)::listpath,outpath,path,line
 integer::u,v,ios
 real(wp)::alpha
 integer::split_index
 if(command_argument_count()/=4) stop 2
 call get_command_argument(1,listpath);call get_command_argument(2,outpath)
 call get_command_argument(3,line);read(line,*)alpha
 call get_command_argument(4,line);read(line,*)split_index
 open(newunit=u,file=trim(listpath),status='old',action='read',iostat=ios);if(ios/=0)stop 2
 open(newunit=v,file=trim(outpath),status='replace',action='write',iostat=ios);if(ios/=0)stop 2
 write(v,'(a)') 'table_path,status,code,admitted,rejected,first_rejected,max_dP_frac,max_dn_frac,core_unchanged,detail'
 do
  read(u,'(a)',iostat=ios)line
  if(ios<0)exit
  if(ios/=0)stop 3
  path=trim(line)
  call load_normalized_table(trim(path),raw,st)
  if(.not.is_ok(st))then
   write(v,'(a)') '"'//trim(path)//'","LOAD_FAIL",'//trim(status_name(st%code))
   cycle
  endif
  call build_controlled_variant(raw,alpha,split_index,variant,va,st)
  if(.not.is_ok(st))then
   write(v,'(a)') '"'//trim(path)//'","VARIANT_FAIL",'//trim(status_name(st%code))
   cycle
  endif
  call init_table_eos(variant,model,audit,st)
  if(is_ok(st))then
   write(v,'(a,a,a,i0,a,i0,a,i0,a,i0,a,es15.7,a,es15.7,a,l1,a,a,a)')'"',trim(path),'","PASS",',st%code,',',audit%admitted_intervals,',',audit%rejected_intervals,',',audit%first_rejected_interval,',',va%max_low_density_pressure_rel_change,',',va%max_low_density_density_rel_change,',',va%core_function_unchanged,',"ok"'
  else
   write(v,'(a,a,a,i0,a,i0,a,i0,a,i0,a,es15.7,a,es15.7,a,l1,a,a,a)')'"',trim(path),'","REJECT",',st%code,',',audit%admitted_intervals,',',audit%rejected_intervals,',',audit%first_rejected_interval,',',va%max_low_density_pressure_rel_change,',',va%max_low_density_density_rel_change,',',va%core_function_unchanged,',"',trim(st%detail),'"'
  endif
 enddo
 close(u);close(v)
end program
