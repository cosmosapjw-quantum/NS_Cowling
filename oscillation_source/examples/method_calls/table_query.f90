program table_query
    use ns_table_io, only: normalized_table_t,load_normalized_table
    use ns_table_eos, only: table_eos_t,table_audit_t,init_table_eos
    use ns_status, only: status_t,is_ok,status_name
    implicit none
    type(normalized_table_t) :: raw
    type(table_eos_t) :: model
    type(table_audit_t) :: audit
    type(status_t) :: status
    character(len=4096) :: input_path,output_path
    integer :: u,ios
    if(command_argument_count()/=2) error stop 'Expected input CSV and output JSON'
    call get_command_argument(1,input_path)
    call get_command_argument(2,output_path)
    call load_normalized_table(trim(input_path),raw,status)
    if(.not.is_ok(status)) error stop 'Could not read normalized table'
    call init_table_eos(raw,model,audit,status)
    open(newunit=u,file=trim(output_path),status='new',action='write',iostat=ios)
    if(ios/=0) error stop 'Output exists or is not writable'
    write(u,'(A)') '{'
    write(u,'(A)') '"status":"'//trim(status_name(status%code))//'",'
    write(u,'(A,I0,A)') '"rows":',audit%row_count,','
    write(u,'(A,I0,A)') '"admitted_intervals":',audit%admitted_intervals,','
    write(u,'(A,I0,A)') '"rejected_intervals":',audit%rejected_intervals,','
    if(model%ready) then
        write(u,'(A)') '"model_ready":true'
    else
        write(u,'(A)') '"model_ready":false'
    end if
    write(u,'(A)') '}'
    close(u)
end program table_query
