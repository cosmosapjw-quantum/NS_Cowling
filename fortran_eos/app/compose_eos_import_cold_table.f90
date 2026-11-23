! Source-preserving cold n,P,energy projection. No interpolation or repair.
! n fm^-3; E/A offset and rest masses MeV; P and total epsilon MeV fm^-3.
! LAYOUT RAW_FILE NEW_OUTPUT_DIRECTORY MN_HEADER_MEV MP_HEADER_MEV [NMIN NMAX]
! Offset-energy layouts add MN; total-energy layouts preserve supplied epsilon.
! The caller owns the raw source/column/mass provenance and before/after hashes.
program compose_eos_import_cold_table
 use compose_eos_types,only:wp,eos_data_t,allocate_eos
 use compose_eos_io,only:write_compose,read_compose,write_internal_csv,write_metadata,write_normalized_csv
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 character(1024)::layout,source_path,output,arg
 character(65536)::line
 character(256)::status
 real(wp)::mn,mp,values(28),relative_readback,nlo,nhi
 real(wp)::density_scale,energy_scale,pressure_scale
 real(wp),parameter::c_si=299792458._wp,mev_fm3_to_j_m3=1.602176634e32_wp
 real(wp),allocatable::raw(:,:)
 type(eos_data_t)::eos,back
 integer::columns,density_column,energy_column,pressure_column,u,ios,rows,i,nt,code
 integer::nraw,first,last,excluded_lower,excluded_upper
 logical::exists,total_energy_density,select_density,cgs_units
 if(command_argument_count()/=5.and.command_argument_count()/=7) &
  call fail('usage: LAYOUT RAW_FILE NEW_OUTPUT_DIRECTORY MN_HEADER_MEV MP_HEADER_MEV [NMIN NMAX]')
 call get_command_argument(1,layout,status=ios);if(ios/=0)call fail('invalid layout argument')
 call get_command_argument(2,source_path,status=ios);if(ios/=0)call fail('invalid input path')
 call get_command_argument(3,output,status=ios);if(ios/=0)call fail('invalid output path')
 call get_command_argument(4,arg,status=ios);if(ios/=0)call fail('invalid rest-mass argument')
 call parse_mass(arg,mn,'rest-mass normalization')
 call get_command_argument(5,arg,status=ios);if(ios/=0)call fail('invalid proton-mass argument')
 call parse_mass(arg,mp,'proton header mass')
 select_density=command_argument_count()==7;nlo=0._wp;nhi=0._wp
 if(select_density)then
  call get_command_argument(6,arg,status=ios);if(ios/=0)call fail('invalid minimum density')
  call parse_scalar(arg,nlo,'minimum density')
  call get_command_argument(7,arg,status=ios);if(ios/=0)call fail('invalid maximum density')
  call parse_scalar(arg,nhi,'maximum density')
  if(nlo<0._wp.or.nhi<nlo)call fail('density window requires 0 <= NMIN <= NMAX')
 endif
 density_column=1;total_energy_density=.false.;cgs_units=.false.
 density_scale=1._wp;energy_scale=1._wp;pressure_scale=1._wp
 select case(trim(layout))
 case('gmrs_2022');columns=11;energy_column=8;pressure_column=9
 case('gogny_2020');columns=9;energy_column=5;pressure_column=6
 case('gogny_full_2020');columns=8;energy_column=2;pressure_column=3
 case('n_eps_pressure')
  columns=3;energy_column=2;pressure_column=3;total_energy_density=.true.
 case('bbp_cgs_2008')
  columns=3;density_column=3;energy_column=1;pressure_column=2
  total_energy_density=.true.;cgs_units=.true.
  ! rho[g cm^-3] c^2 -> J m^-3 -> MeV fm^-3; dyn cm^-2 = 0.1 Pa.
  ! Exact SI c and elementary charge define the conversion, no added rest mass.
  density_scale=1.e-39_wp
  energy_scale=(1000._wp*c_si**2)/mev_fm3_to_j_m3
  pressure_scale=0.1_wp/mev_fm3_to_j_m3
 case('bsk_2018');columns=28;energy_column=10;pressure_column=8
 case('user_core_energy_density')
  columns=7;energy_column=2;pressure_column=3;total_energy_density=.true.
 case('user_crust_energy_density')
  columns=5;density_column=2;energy_column=3;pressure_column=4;total_energy_density=.true.
 case default;call fail('unknown raw-table layout; columns are never inferred')
 end select
 inquire(file=trim(output),exist=exists)
 if(exists)call fail('output already exists; choose a new evidence directory')
 open(newunit=u,file=trim(source_path),status='old',action='read',iostat=ios)
 if(ios/=0)call fail('cannot open raw table')
 rows=0
 do
  read(u,'(a)',iostat=ios)line
  if(ios<0)exit
  if(ios/=0)call fail('raw table read error')
  call parse_row(line,values,nt)
  if(nt==0)cycle
  rows=rows+1
 enddo
 if(rows<3)call fail('at least three raw density rows are required')
 allocate(raw(3,rows));rewind(u);i=0
 do
  read(u,'(a)',iostat=ios)line
  if(ios<0)exit
  if(ios/=0)call fail('raw table read error')
  call parse_row(line,values,nt)
  if(nt==0)cycle
  i=i+1
  if(i>rows)call fail('source row count changed during read')
  raw(:,i)=[values(density_column)*density_scale,values(energy_column)*energy_scale, &
   values(pressure_column)*pressure_scale]
 enddo
 close(u)
 if(i/=rows)call fail('source row count changed during read')
 if(.not.all(ieee_is_finite(raw)))call fail('nonfinite converted source values')
 if(any(raw(1,:)<=0))call fail('density must be positive')
 if(any(raw(1,2:)<=raw(1,:rows-1)))call fail('density must be strictly increasing')
 ! Full source validation precedes the explicitly authorized inclusive selection.
 nraw=rows;excluded_lower=0;excluded_upper=0
 if(select_density)then
  excluded_lower=count(raw(1,:)<nlo);excluded_upper=count(raw(1,:)>nhi)
 endif
 first=1+excluded_lower;last=nraw-excluded_upper;rows=last-first+1
 if(rows<3)call fail('density window must retain at least three source rows')
 call allocate_eos(eos,rows)
 eos%mn=mn;eos%mp=mp;eos%il=1;eos%beta_equilibrated=.true.
 ! Total cold one-charge barotrope only: species and phase labels are unavailable.
 eos%nb=raw(1,first:last)
 eos%q(1,:)=raw(3,first:last)/eos%nb
 if(total_energy_density)then
  if(any(raw(2,first:last)<=0))call fail('nonpositive total energy density')
  ! Header MN is an encoding scale, never an added rest energy in this layout.
  eos%q(6,:)=(raw(2,first:last)/eos%nb)/mn-1._wp
  eos%q(3,:)=(raw(2,first:last)/eos%nb+eos%q(1,:))/mn-1._wp
 else
  if(any(mn+raw(2,first:last)<=0))call fail('nonpositive total energy per baryon')
  eos%q(3,:)=(raw(2,first:last)+eos%q(1,:))/mn
  eos%q(6,:)=raw(2,first:last)/mn
 endif
 eos%q(7,:)=eos%q(6,:)
 if(.not.all(ieee_is_finite(eos%q)))call fail('nonfinite normalized thermodynamics')
 ! Atomic directory creation prevents stale files and duplicate-worker overwrite.
 ! Its parent must exist; orchestration creates the campaign root.
 call execute_command_line('mkdir -- '//shellquote(trim(output)),exitstat=code,cmdstat=ios)
 if(ios/=0.or.code/=0)call fail('cannot create new output directory')
 call write_compose(trim(output),eos,status);call require_ok(status)
 call read_compose(trim(output),back,status);call require_ok(status)
 if(back%n/=rows)call fail('native readback row count differs')
 if(any(back%nb/=eos%nb))call fail('native readback density differs')
 relative_readback=maxval(abs(back%q-eos%q)/max(1._wp,abs(eos%q)))
 if(relative_readback>2e-14_wp)call fail('native readback thermodynamics differ')
 call write_internal_csv(trim(output)//'/eos_data.csv',back,status);call require_ok(status)
 call write_metadata(trim(output)//'/eos_data.nml',back,status);call require_ok(status)
 call write_normalized_csv(trim(output)//'/normalized_eos.csv',back,status);call require_ok(status)
 open(newunit=u,file=trim(output)//'/SOURCE_VALUES.csv',status='new',action='write',iostat=ios)
 if(ios/=0)call fail('cannot write source-values evidence')
 if(total_energy_density)then
  write(u,'(a)')'n_b_fm3,energy_mev_fm3,pressure_mev_fm3'
 else
  write(u,'(a)')'n_b_fm3,energy_per_baryon_offset_mev,pressure_mev_fm3'
 endif
 ! Evidence retains all input n/energy/P rows, including explicitly excluded rows.
 do i=1,nraw
  write(u,'(es25.16e3,2(",",es25.16e3))')raw(:,i)
 enddo
 close(u)
 open(newunit=u,file=trim(output)//'/IMPORT.json',status='new',action='write',iostat=ios)
 if(ios/=0)call fail('cannot write projection receipt')
 write(u,'(a)')'{"schema":"ns_cowling.cold_raw_projection.v2",'
 write(u,'(a)')'"status":"NATIVE_RAW_PROJECTION_READBACK_VERIFIED",'
 write(u,'(a)')'"layout":"'//trim(layout)//'",'
 write(u,'(a,i0,a)')'"rows":',rows,','
 write(u,'(a,i0,a)')'"raw_rows":',nraw,','
 write(u,'(a,i0,a)')'"selected_rows":',rows,','
 write(u,'(a,i0,a)')'"excluded_rows":',nraw-rows,','
 write(u,'(a,i0,a)')'"excluded_lower_rows":',excluded_lower,','
 write(u,'(a,i0,a)')'"excluded_upper_rows":',excluded_upper,','
 write(u,'(a,i0,a)')'"source_row_start_zero_based":',first-1,','
 write(u,'(a,i0,a)')'"source_row_stop_exclusive":',last,','
 write(u,'(a,es25.16e3,a,es25.16e3,a)')'"raw_density_domain_fm3":[',raw(1,1),',',raw(1,nraw),'],'
 write(u,'(a,es25.16e3,a,es25.16e3,a)')'"selected_density_domain_fm3":[',eos%nb(1),',',eos%nb(rows),'],'
 if(select_density)then
  write(u,'(a)')'"density_selection_applied":true,'
  write(u,'(a,es25.16e3,a,es25.16e3,a)')'"requested_density_domain_fm3":[',nlo,',',nhi,'],'
 else
  write(u,'(a)')'"density_selection_applied":false,"requested_density_domain_fm3":null,'
 endif
 write(u,'(a,i0,a,i0,a,i0,a)')'"source_columns_one_based":[',density_column,',',energy_column,',',pressure_column,'],'
 if(total_energy_density)then
  write(u,'(a)')'"energy_representation":"total_energy_density_mev_fm3","rest_mass_offset_mev":0,'
 else
  write(u,'(a)')'"energy_representation":"energy_per_baryon_offset_mev",'
  write(u,'(a,es25.16e3,a)')'"rest_mass_offset_mev":',mn,','
 endif
 write(u,'(a,es25.16e3,a)')'"neutron_header_mass_mev":',mn,','
 write(u,'(a,es25.16e3,a)')'"proton_header_mass_mev":',mp,','
 write(u,'(a,es25.16e3,a)')'"max_scaled_thermo_readback_difference":',relative_readback,','
 if(cgs_units)then
  write(u,'(a)')'"source_units":["cm^-3","g cm^-3","dyn cm^-2"],'
  write(u,'(a)')'"unit_conversion":"CGS_TO_MEV_FM3_EXACT_SI_CONSTANTS",'
 else
  if(total_energy_density)then
   write(u,'(a)')'"source_units":["fm^-3","MeV fm^-3","MeV fm^-3"],'
  else
   write(u,'(a)')'"source_units":["fm^-3","MeV","MeV fm^-3"],'
  endif
  write(u,'(a)')'"unit_conversion":"IDENTITY",'
 endif
 write(u,'(a,es25.16e3,a,es25.16e3,a,es25.16e3,a)') &
  '"source_unit_scale_to_internal":[',density_scale,',',energy_scale,',',pressure_scale,'],'
 write(u,'(a)')'"source_values_csv":"all_source_rows_in_declared_internal_units",'
 write(u,'(a)')'"chemical_potential":"Euler_from_supplied_pressure_and_energy",'
 write(u,'(a)')'"source_energy_pressure":"preserved_without_interpolation_or_repair",'
 write(u,'(a)')'"composition_available":false,"phase_labels_reconstructed":false,'
 write(u,'(a)')'"physical_acceptance_asserted":false}'
 close(u)
 write(*,'(a,i0)')'NATIVE_RAW_PROJECTION_READBACK_VERIFIED rows=',rows
contains
 subroutine parse_mass(text,value,label)
 character(*),intent(in)::text,label
 real(wp),intent(out)::value
 call parse_scalar(text,value,label)
 if(value<=0)call fail(label//' must be positive')
 end subroutine
 subroutine parse_scalar(text,value,label)
 character(*),intent(in)::text,label
 real(wp),intent(out)::value
 character(len(text))::token
 integer::io
 token=adjustl(text)
 if(len_trim(token)==0)call fail(label//' is empty')
 if(verify(trim(token),'0123456789+-.eEdD')/=0)call fail(label//' must be one numeric token')
 read(token,*,iostat=io)value
 if(io/=0)call fail(label//' is not numeric')
 if(.not.ieee_is_finite(value))call fail(label//' must be finite')
 end subroutine
 subroutine fail(message)
 character(*),intent(in)::message
 write(*,'(a)')'ERROR: '//trim(message)
 stop 2
 end subroutine
 subroutine require_ok(message)
 character(*),intent(in)::message
 if(message/='ok')call fail(message)
 end subroutine
 function shellquote(value)result(quoted)
 character(*),intent(in)::value
 character(:),allocatable::quoted
 integer::j
 quoted="'"
 do j=1,len_trim(value)
  if(value(j:j)=="'")then
   quoted=quoted//"'"//achar(34)//"'"//achar(34)//"'"
  else
   quoted=quoted//value(j:j)
  endif
 enddo
 quoted=quoted//"'"
 end function
 subroutine parse_row(text,nums,count)
 character(*),intent(in)::text
 real(wp),intent(out)::nums(28)
 integer,intent(out)::count
 character(len(text))::work
 character(96)::token
 integer::j,k,finish,io
 nums=0;count=0;work=adjustl(text);finish=len_trim(work)
 if(finish==0)return
 if(work(1:1)=='#')return
 if(finish==len(work))call fail('raw row exceeds parser line capacity')
 j=1
 do while(j<=finish)
  if(work(j:j)==' '.or.iachar(work(j:j))==9)then
   j=j+1;cycle
  endif
  k=j
  do while(k<=finish)
   if(work(k:k)==' '.or.iachar(work(k:k))==9)exit
   k=k+1
  enddo
  if(k-j>len(token))call fail('raw numeric token is too long')
  count=count+1
  if(count>columns)call fail('raw row has extra columns')
  token=work(j:k-1)
  if(verify(trim(token),'0123456789+-.eEdD')/=0)call fail('invalid or nonfinite numeric token')
  read(token,*,iostat=io)nums(count)
  if(io/=0)call fail('invalid numeric token')
  if(.not.ieee_is_finite(nums(count)))call fail('nonfinite raw table value')
  j=k
 enddo
 if(count/=columns)call fail('raw row column count differs from explicit layout')
 end subroutine
end program
