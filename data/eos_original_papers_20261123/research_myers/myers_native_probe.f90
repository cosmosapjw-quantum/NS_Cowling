program myers_native_probe
 use compose_eos_types, only: eos_data_t
 use compose_eos_myers1966, only: generate_myers_bps_native
 implicit none
 type(eos_data_t)::eos
 character(512)::source,out,arg,status
 integer::magic,points,ios
 call get_command_argument(1,source);call get_command_argument(2,out)
 call get_command_argument(3,arg);read(arg,*,iostat=ios)magic
 if(ios/=0)stop 2
 call get_command_argument(4,arg);read(arg,*,iostat=ios)points
 if(ios/=0)stop 2
 call generate_myers_bps_native(trim(source),points,eos,trim(out),status,magic)
 print '(a)',trim(status)
 if(status/='ok')stop 2
 print '(a,i0)', 'MYERS_PHASE_COMPLETE_ROWS=',eos%n
end program
