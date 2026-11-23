! C20 native finite contract: exact affine/shape examples, legacy ABI,
! phase-local batch interpolation, weak-field geometry and invalid inputs.
program cowling_c20_tov_pchip
  use iso_c_binding, only:c_double,c_int
  use, intrinsic :: ieee_arithmetic, only:ieee_value,ieee_quiet_nan
  use compose_eos_cowling_pchip
  use compose_eos_cowling_tov_geometry
  implicit none
  real(c_double)::s2(2),y2(2,2),d2(2,2),q2(3),v2(2,3)
  real(c_double)::s(4),y(2,4),d(2,4),dc(2,4),q(4),v(2,4),expected(2,4)
  real(c_double)::y5(5,4),d5(5,4),d5g(5,4),v5(5,4),v5g(5,4)
  real(c_double)::h(4),phase_v(2,4),rm(2,2),p(2),cs(2),hg(2),geo(5,2)
  real(c_double)::tol,nan
  integer(c_int)::status
  integer::j
  tol=2.e-13_c_double
  nan=ieee_value(1._c_double,ieee_quiet_nan)

  s2=[0._c_double,2._c_double]
  y2(:,1)=[1._c_double,4._c_double]
  y2(:,2)=[5._c_double,-2._c_double]
  call cowling_native_pchip_slopes_fields(2_c_int,2_c_int,s2,y2,d2,status)
  call require(status==0.and.all(d2(1,:)==2._c_double).and.all(d2(2,:)==-3._c_double),10)
  q2=[0._c_double,1._c_double,2._c_double]
  call cowling_native_pchip_values_fields(2_c_int,2_c_int,s2,y2,d2,3_c_int,q2,v2,status)
  call require(status==0.and.maxval(abs(v2(1,:)-[1._c_double,3._c_double,5._c_double]))<tol,11)
  call require(maxval(abs(v2(2,:)-[4._c_double,1._c_double,-2._c_double]))<tol,12)

  s=[0._c_double,.25_c_double,1._c_double,2._c_double]
  y(1,:)=1._c_double+2._c_double*s
  y(2,:)=[1._c_double,1._c_double,2._c_double,0._c_double]
  call cowling_native_pchip_slopes_fields(2_c_int,4_c_int,s,y,d,status)
  call require(status==0.and.maxval(abs(d(1,:)-2._c_double))<tol,13)
  call require(maxval(abs(d(2,:)-[0._c_double,0._c_double,0._c_double,-82._c_double/21]))<tol,14)
  call tov_native_pchip_slopes(4_c_int,s,y,dc,0_c_int,status)
  call require(status==0.and.maxval(abs(dc-d))<tol,35)
  call tov_native_pchip_slopes(4_c_int,s,y,dc,1_c_int,status)
  call require(status==0.and.maxval(abs(dc-d))<tol,36)
  call tov_native_pchip_slopes(4_c_int,s,y,dc,2_c_int,status)
  call require(status/=0.and.all(dc==0._c_double),37)
  q=[2._c_double,.625_c_double,0._c_double,1.5_c_double]
  expected(1,:)=1._c_double+2._c_double*q
  expected(2,:)=[0._c_double,1.5_c_double,1._c_double,125._c_double/84]
  call cowling_native_pchip_values_fields(2_c_int,4_c_int,s,y,d,4_c_int,q,v,status)
  call require(status==0.and.maxval(abs(v-expected))<tol,15)
  do j=1,5
    y5(j,:)=real(j,c_double)*y(1,:)+y(2,:)
  enddo
  call cowling_native_pchip_slopes(4_c_int,s,y5,d5,status)
  call require(status==0,16)
  call cowling_native_pchip_slopes_fields(5_c_int,4_c_int,s,y5,d5g,status)
  call require(status==0.and.all(d5==d5g),17)
  call cowling_native_pchip_values(4_c_int,s,y5,d5,4_c_int,q,v5,status)
  call require(status==0,18)
  call cowling_native_pchip_values_fields(5_c_int,4_c_int,s,y5,d5,4_c_int,q,v5g,status)
  call require(status==0.and.all(v5==v5g),19)
  h=5._c_double-q*q
  call tov_native_phase_values(4_c_int,5._c_double,s,y,d,4_c_int,h,phase_v,status)
  call require(status==0.and.maxval(abs(phase_v-expected))<tol,20)
  h(1)=6._c_double
  call tov_native_phase_values(4_c_int,5._c_double,s,y,d,4_c_int,h,phase_v,status)
  call require(status/=0.and.all(phase_v==0._c_double),34)
  h(1)=-1._c_double
  call tov_native_phase_values(4_c_int,5._c_double,s,y,d,4_c_int,h,phase_v,status)
  call require(status/=0.and.all(phase_v==0._c_double),21)

  q(1)=s(4)+.01_c_double
  call cowling_native_pchip_values_fields(2_c_int,4_c_int,s,y,d,4_c_int,q,v,status)
  call require(status/=0.and.all(v==0._c_double),22)
  s(2)=s(1)
  call cowling_native_pchip_slopes_fields(2_c_int,4_c_int,s,y,d,status)
  call require(status/=0.and.all(d==0._c_double),23)
  s(2)=-.1_c_double
  call cowling_native_pchip_slopes_fields(2_c_int,4_c_int,s,y,d,status)
  call require(status/=0.and.all(d==0._c_double),24)
  s(2)=.25_c_double;y(1,2)=nan
  call cowling_native_pchip_slopes_fields(2_c_int,4_c_int,s,y,d,status)
  call require(status/=0.and.all(d==0._c_double),25)

  rm(:,1)=[10._c_double,1.e-20_c_double]
  rm(:,2)=[20._c_double,1.e-20_c_double]
  p=0._c_double;cs=.5_c_double;hg=[.1_c_double,0._c_double]
  call tov_native_geometry(2_c_int,hg,0._c_double,1.e-20_c_double,20._c_double,rm,p,cs,geo,status)
  call require(status==0,26)
  call require(abs(geo(2,2)+5.e-22_c_double)<1.e-36_c_double,27)
  call require(abs(geo(3,1)-1.e-21_c_double)<1.e-36_c_double,28)
  call require(abs(geo(4,1)-1.e-22_c_double)<1.e-36_c_double,29)
  call require(all(geo(1,:)==rm(1,:)).and.all(geo(5,:)==cs),30)
  ! Frozen Python float/libm-power oracle for the isolated algebra kernel.
  ! This sampled pair is not asserted to solve a TOV equilibrium.
  rm(:,1)=[.3_c_double,0._c_double]
  rm(:,2)=[1.2345_c_double,0._c_double]
  p=1._c_double
  call tov_native_geometry(2_c_int,hg,0._c_double,0._c_double,2._c_double,rm,p,cs,geo,status)
  call require(status==0,38)
  call require(geo(4,1)==3.7699111843077513_c_double,39)
  call require(geo(4,2)==15.513184523426398_c_double,40)
  rm(2,1)=5._c_double
  call tov_native_geometry(2_c_int,hg,0._c_double,1.e-20_c_double,20._c_double,rm,p,cs,geo,status)
  call require(status/=0.and.all(geo==0._c_double),31)
  rm(2,1)=1.e-20_c_double;cs(1)=0._c_double
  call tov_native_geometry(2_c_int,hg,0._c_double,1.e-20_c_double,20._c_double,rm,p,cs,geo,status)
  call require(status/=0.and.all(geo==0._c_double),32)
  cs=.5_c_double;p(1)=nan
  call tov_native_geometry(2_c_int,hg,0._c_double,1.e-20_c_double,20._c_double,rm,p,cs,geo,status)
  call require(status/=0.and.all(geo==0._c_double),33)
  print '(A)', 'C20_NATIVE_TOV_PCHIP_AND_GEOMETRY_PASS'
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    integer,intent(in)::label
    if(.not.condition)then
      print '(A,I0)', 'C20_NATIVE_ASSERTION_FAILED ',label
      stop 1
    endif
  end subroutine
end program cowling_c20_tov_pchip
