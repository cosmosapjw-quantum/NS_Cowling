program audit_many
 use ns_status
 use ns_table_io, only: normalized_table_t,load_normalized_table
 use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos
 implicit none
 type(normalized_table_t)::raw
 type(table_eos_t)::model
 type(table_audit_t)::audit
 type(status_t)::st
 character(len=2048)::listpath,outpath,path,line
 integer::u,v,ios,k
 if(command_argument_count()/=2) stop 2
 call get_command_argument(1,listpath);call get_command_argument(2,outpath)
 open(newunit=u,file=trim(listpath),status='old',action='read',iostat=ios)
 if(ios/=0) stop 2
 open(newunit=v,file=trim(outpath),status='replace',action='write',iostat=ios)
 if(ios/=0) stop 2
 write(v,'(a)') 'table_path,status,code,admitted,rejected,first_rejected,min_d2p,min_causal_margin,category,detail'
 do
  read(u,'(a)',iostat=ios)line
  if(ios<0)exit
  if(ios/=0)stop 3
  path=trim(line)
  call load_normalized_table(trim(path),raw,st)
  if(.not.is_ok(st))then
   write(v,'(a,a,a,i0,a,i0,a,a,a,a,a)')'"',trim(path),'","LOAD_FAIL",',st%code,',0,0,0,0,0,"',trim(st%category),'","',trim(st%detail),'"'
   cycle
  endif
  call init_table_eos(raw,model,audit,st)
  if(is_ok(st))then
   write(v,'(a,a,a,i0,a,i0,a,i0,a,i0,a,es24.15,a,es24.15,a,a,a,a,a)')'"',trim(path),'","PASS",',st%code,',',audit%admitted_intervals,',',audit%rejected_intervals,',',audit%first_rejected_interval,',',audit%first_bad_min_d2p,',',audit%first_bad_min_causal_margin,',"',trim(st%category),'","',trim(st%detail),'"'
  else
   write(v,'(a,a,a,i0,a,i0,a,i0,a,i0,a,es24.15,a,es24.15,a,a,a,a,a)')'"',trim(path),'","REJECT",',st%code,',',audit%admitted_intervals,',',audit%rejected_intervals,',',audit%first_rejected_interval,',',audit%first_bad_min_d2p,',',audit%first_bad_min_causal_margin,',"',trim(st%category),'","',trim(st%detail),'"'
  endif
 enddo
 close(u);close(v)
end program
