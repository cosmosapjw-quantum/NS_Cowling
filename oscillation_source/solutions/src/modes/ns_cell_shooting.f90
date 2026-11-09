! Direct integration in native EOS cells. This adapter calls the selected
! workspace's TOV and Cowling equations; it does not read a saved spectrum.
module ns_cell_shooting
 use ns_precision, only: wp
 use ns_status
 use ns_physical_constants, only: c_si,g_si
 use ns_table_eos, only: table_eos_t, eos_table_at_pressure
 use ns_eos_types
 use ns_cell_eos
 use ns_tov, only: center_start,tov_rhs
 use ns_cowling, only: cowling_rhs_point
 use ns_center_state, only: cowling_center_state
 use ns_surface_residual, only: cowling_surface_residual
 use ns_numerics_types, only: ode_options_t,ode_result_t,step_result_t
 use ns_adaptive_rk, only: adaptive_integrate
 use ns_rk4, only: rk4_step
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 private
 real(wp),parameter::rscale=1e4_wp,mscale=1e30_wp,pi=acos(-1._wp)
 type,public::cell_run_t
  real(wp)::radius_m=0._wp,mass_kg=0._wp,nu_surface=0._wp,mu_surface=0._wp,residual=0._wp
  real(wp)::terminal_w=0._wp,terminal_v=0._wp,terminal_mu=0._wp
  real(wp)::surface_pressure_rel_error=0._wp,boundary_residual=0._wp
  integer::steps=0,cells=0
  type(status_t)::status
  logical::valid=.false.
 end type
 public::integrate_cells,match_cells
contains
 subroutine integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,nu_surface,result,profile_path,match_index,surface_start,r_limit_m)
  type(table_eos_t),intent(in)::model
  real(wp),intent(in)::pc,p_cut,r0,rtol,step_cap,f_hz,nu_surface
  character(len=*),intent(in)::method
  type(cell_run_t),intent(out)::result
  character(len=*),optional,intent(in)::profile_path
  integer,optional,intent(in)::match_index
  type(cell_run_t),optional,intent(in)::surface_start
  real(wp),optional,intent(in)::r_limit_m
  type(eos_state_t)::s0,ss,stage_state
  type(status_t)::st
  type(ode_options_t)::opts
  type(ode_result_t)::traj
  type(step_result_t)::step
  real(wp)::yc(3),t0,ts,tstart,tend,width,hcell,y(4),dy(4),dn,mu_s
  real(wp)::x,dx,phi,lambda,qterm,dytov(3),w0,v0,rr,mm
  real(wp)::yt(4)
  integer::i,ic,is,j,nd,u,ios,attempts,di,original_surface_cell
  real(wp)::tdir,b1,b2,normb,radius_limit
  logical::modes,save
  result=cell_run_t();modes=f_hz>0._wp;save=present(profile_path)
  radius_limit=1e6_wp
  if(present(r_limit_m))radius_limit=r_limit_m
  call fail_status(result%status,NS_INVALID_INPUT,'numerics','invalid direct-cell settings')
  if(.not.model%ready.or..not.all(ieee_is_finite([pc,p_cut,r0,rtol,step_cap,f_hz,nu_surface])))return
  if(pc<=p_cut.or.p_cut<=0._wp.or.r0<=0._wp.or.rtol<=0._wp.or.step_cap<=0._wp.or.f_hz<0._wp.or. &
     .not.ieee_is_finite(radius_limit).or.radius_limit<=r0)return
  if(method/='dopri'.and.method/='rk4')return
  nd=2;if(modes)nd=4
  call center_start(lookup,pc,r0,0._wp,yc,st)
  if(.not.is_ok(st))then;result%status=st;return;endif
  call table_pressure_coordinate(model,yc(2),ic,t0,s0,st)
  if(.not.is_ok(st))then;result%status=st;return;endif
  call table_pressure_coordinate(model,p_cut,is,ts,ss,st)
  if(.not.is_ok(st))then;result%status=st;return;endif
  mu_s=ss%mu_j;result%mu_surface=mu_s
  original_surface_cell=is;di=-1;tdir=-1._wp
  y=[r0/rscale,yc(1)/mscale,0._wp,0._wp]
  if(modes)then
   call cowling_center_state(2,r0,1._wp,w0,v0,st)
   if(.not.is_ok(st))then;result%status=st;return;endif
   y(3)=w0/r0**3;y(4)=v0/r0**2
  endif
  if(present(match_index))then
   if(match_index<=original_surface_cell.or.match_index>ic)then
    call fail_status(result%status,NS_OUT_OF_DOMAIN,'physics','matching knot outside stellar interior');return
   endif
   if(present(surface_start))then
    if(.not.surface_start%valid.or..not.modes)then
     call fail_status(result%status,NS_INVALID_INPUT,'implementation','inward solution needs valid surface and frequency');return
    endif
    ic=original_surface_cell;t0=ts;is=match_index-1;ts=1._wp;di=1;tdir=1._wp
    rr=surface_start%radius_m;mm=surface_start%mass_kg
    stage_state=ss
    call tov_rhs(stage_eos,rr,[mm,ss%pressure_pa,0._wp],dytov,st)
    if(.not.is_ok(st))then;result%status=st;return;endif
    lambda=-0.5_wp*log(1._wp-2._wp*g_si*mm/(rr*c_si*c_si))
    b1=rr*dytov(3);b2=(2._wp*pi*f_hz*rr/c_si)**2*exp(lambda-2._wp*nu_surface)
    normb=sqrt(b1*b1+b2*b2)
    y=[rr/rscale,mm/mscale,b2/normb,-b1/normb]
    call cowling_surface_residual(rr,2._wp*pi*f_hz,nu_surface,lambda,dytov(3), &
         rr**3*y(3),rr**2*y(4),result%boundary_residual,qterm,st)
    if(.not.is_ok(st))then;result%status=st;return;endif
   else
    is=match_index;ts=0._wp
   endif
  endif
  if(save)then
   open(newunit=u,file=profile_path,status='replace',action='write',iostat=ios)
   if(ios/=0)then
    call fail_status(result%status,NS_INVALID_INPUT,'environment','cannot open profile output');return
   endif
   write(u,'(a)')'cell,t,r_m,m_kg,p_pa,energy_j_m3,chi,nu,w_over_r3,v_over_r2'
  endif
  opts%rtol=rtol;allocate(opts%atol(nd))
  opts%atol(1:2)=[rtol*1e-7_wp,rtol*1e-10_wp]
  if(modes)opts%atol(3:4)=rtol*0.01_wp
  opts%h_min=1e-20_wp;opts%h_max=min(step_cap,1._wp);opts%max_steps=100000
  do i=ic,is,di
   tstart=1._wp;if(di==1)tstart=0._wp
   if(i==ic)tstart=t0
   tend=0._wp;if(di==1)tend=1._wp
   if(i==is)tend=ts
   width=abs(tstart-tend)
   if(width<=0._wp)cycle
   hcell=log(model%raw%pressure_pa(i+1))-log(model%raw%pressure_pa(i))
   call rhs(0._wp,y(:nd),dy(:nd),st)
   if(.not.is_ok(st))exit
   if(save)call write_row(0._wp,y)
   opts%h_initial=min(opts%h_max,0.01_wp*abs(y(1)/dy(1)),width)
   if(method=='dopri')then
    call adaptive_integrate(rhs,0._wp,width,y(:nd),opts,traj)
    st=traj%status
    if(.not.is_ok(st).or..not.traj%valid)exit
    y(:nd)=traj%y(:,size(traj%x))
    result%steps=result%steps+traj%accepted
    if(save)then
     do j=2,size(traj%x)
      yt=0._wp;yt(:nd)=traj%y(:,j);call write_row(traj%x(j),yt)
     enddo
    endif
   else
    x=0._wp;attempts=0
    do while(x<width)
     call rhs(x,y(:nd),dy(:nd),st);if(.not.is_ok(st))exit
     ! The additional radial cap regularizes the tiny first central substeps.
     ! Halving step_cap halves both caps; no error test is used in this RK4 comparison.
     dx=min(step_cap,width-x,step_cap*abs(y(1)/dy(1)))
     if(x+dx==x)then
      call fail_status(st,NS_STAGNATION,'numerics','RK4 cell coordinate stagnation');exit
     endif
     call rk4_step(rhs,x,y(:nd),dx,step);st=step%status
     if(.not.is_ok(st))exit
     x=x+dx;y(:nd)=step%y;attempts=attempts+1
     result%steps=result%steps+1
     if(save)call write_row(x,y)
     if(attempts>1000000)then
      call fail_status(st,NS_MAX_STEPS,'numerics','RK4 per-cell step budget');exit
     endif
    enddo
    if(.not.is_ok(st))exit
   endif
   result%cells=result%cells+1
  enddo
  if(save)close(u)
  result%status=st;if(.not.is_ok(st))return
  rr=rscale*y(1);mm=mscale*y(2)
  result%radius_m=rr;result%mass_kg=mm
  result%terminal_w=y(3);result%terminal_v=y(4)
  call table_cell_state(model,is,ts,s0,dn,st)
  result%terminal_mu=s0%mu_j
  if(2._wp*g_si*mm/(rr*c_si*c_si)>=1._wp)then
   call fail_status(result%status,NS_OUT_OF_DOMAIN,'physics','invalid surface compactness');return
  endif
  result%nu_surface=0.5_wp*log(1._wp-2._wp*g_si*mm/(rr*c_si*c_si))
  result%surface_pressure_rel_error=abs(ss%pressure_pa/p_cut-1._wp)
  if(modes.and..not.present(match_index))then
   stage_state=ss
   call tov_rhs(stage_eos,rr,[mm,ss%pressure_pa,0._wp],dytov,st)
   if(.not.is_ok(st))then;result%status=st;return;endif
   phi=dytov(3);lambda=-result%nu_surface
   call cowling_surface_residual(rr,2._wp*pi*f_hz,nu_surface,lambda,phi, &
       rr**3*y(3),rr**2*y(4),result%residual,qterm,st)
   if(.not.is_ok(st))then;result%status=st;return;endif
  endif
  result%valid=.true.;call clear_status(result%status)
 contains
  subroutine table_cell_state(eos,cell_index,t,state,dn_dmu,status)
   type(table_eos_t),intent(in)::eos
   integer,intent(in)::cell_index
   real(wp),intent(in)::t
   type(eos_state_t),intent(out)::state
   real(wp),intent(out)::dn_dmu
   type(status_t),intent(out)::status
   real(wp)::u0,u1,p
   dn_dmu=0._wp
   if(cell_index<1.or.cell_index>=eos%raw%n.or.t<0._wp.or.t>1._wp)then
    call clear_eos_state(state)
    call fail_status(status,NS_OUT_OF_DOMAIN,'physics','invalid log-pressure cell coordinate')
    return
   endif
   u0=log(eos%raw%pressure_pa(cell_index))
   u1=log(eos%raw%pressure_pa(cell_index+1))
   p=exp(u0+t*(u1-u0))
   call eos_table_at_pressure(eos,p,SIDE_REGULAR,state,status)
   if(.not.is_ok(status))return
   dn_dmu=state%n_b_m3**2/(state%chi*(state%energy_j_m3+state%pressure_pa))
  end subroutine table_cell_state
  subroutine lookup(p,side,s,status)
   real(wp),intent(in)::p
   integer,intent(in)::side
   type(eos_state_t),intent(out)::s
   type(status_t),intent(out)::status
   real(wp)::t
   integer::ii
   call table_pressure_coordinate(model,p,ii,t,s,status)
  end subroutine
  subroutine stage_eos(p,side,s,status)
   real(wp),intent(in)::p
   integer,intent(in)::side
   type(eos_state_t),intent(out)::s
   type(status_t),intent(out)::status
   s=stage_state
   if(abs(p-s%pressure_pa)>64._wp*epsilon(1._wp)*p)then
    call fail_status(status,NS_INVALID_INPUT,'implementation','stage EOS pressure mismatch')
   else
    call clear_status(status)
   endif
  end subroutine
  subroutine rhs(xx,z,dz,status)
   real(wp),intent(in)::xx,z(:)
   real(wp),intent(out)::dz(:)
   type(status_t),intent(out)::status
   real(wp)::t,r,m,phi,lam,nu,drdx,dn,dtov(3),dcow(2),qt
   dz=0._wp;t=tstart+tdir*xx
   ! Only snap roundoff from the exact cell endpoint; never cross an EOS interval.
   if(t<0._wp.and.t>=-16._wp*epsilon(1._wp))t=0._wp
   if(t>1._wp.and.t<=1._wp+16._wp*epsilon(1._wp))t=1._wp
   call table_cell_state(model,i,t,stage_state,dn,status);if(.not.is_ok(status))return
   r=z(1)*rscale;m=z(2)*mscale
   if(r<=0._wp.or.r>radius_limit.or.m<0._wp)then
    call fail_status(status,NS_OUT_OF_DOMAIN,'physics','invalid direct-cell stellar state');return
   endif
   call tov_rhs(stage_eos,r,[m,stage_state%pressure_pa,0._wp],dtov,status)
   if(.not.is_ok(status))return
   phi=dtov(3);drdx=-tdir*hcell*stage_state%pressure_pa/(stage_state%n_b_m3*stage_state%mu_j*phi)
   dz(1)=drdx/rscale;dz(2)=dtov(1)*drdx/mscale
   if(modes)then
    lam=-0.5_wp*log(1._wp-2._wp*g_si*m/(r*c_si*c_si))
    nu=nu_surface+log(mu_s/stage_state%mu_j)
    ! A common factor r^2 cancels from the linear Cowling RHS.
    call cowling_rhs_point(r,2._wp*pi*f_hz,stage_state%chi,nu,lam,phi,[r*z(3),z(4)],dcow,qt,status)
    if(.not.is_ok(status))return
    dz(3)=(dcow(1)/r-3._wp*z(3)/r)*drdx
    dz(4)=(dcow(2)-2._wp*z(4)/r)*drdx
   endif
   if(.not.all(ieee_is_finite(dz)))call fail_status(status,NS_NONFINITE_RESULT,'numerics','nonfinite cell derivative')
  end subroutine
  subroutine write_row(xx,z)
   real(wp),intent(in)::xx,z(4)
   type(eos_state_t)::ssave
   type(status_t)::stsave
   real(wp)::tt,dd,nn
   tt=max(0._wp,min(1._wp,tstart+tdir*xx))
   call table_cell_state(model,i,tt,ssave,dd,stsave)
   nn=nu_surface+log(mu_s/ssave%mu_j)
   write(u,'(i0,9(",",es25.16e3))')i,tt,z(1)*rscale,z(2)*mscale,ssave%pressure_pa, &
        ssave%energy_j_m3,ssave%chi,nn,z(3),z(4)
  end subroutine
 end subroutine
subroutine match_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star,match_index, &
                       residual,r_mismatch,m_mismatch,status,profile_prefix,step_count,surface_residual,r_limit_m)
  type(table_eos_t),intent(in)::model
  real(wp),intent(in)::pc,p_cut,r0,rtol,step_cap,f_hz
  character(len=*),intent(in)::method
  type(cell_run_t),intent(in)::star
  integer,intent(in)::match_index
  real(wp),intent(out)::residual,r_mismatch,m_mismatch
  type(status_t),intent(out)::status
  character(len=*),optional,intent(in)::profile_prefix
  integer,optional,intent(out)::step_count
  real(wp),optional,intent(out)::surface_residual
  real(wp),optional,intent(in)::r_limit_m
  type(cell_run_t)::a,b
  real(wp)::za(2),zb(2),na,nb
  residual=0._wp;r_mismatch=0._wp;m_mismatch=0._wp
  if(present(step_count))step_count=0
  if(present(surface_residual))surface_residual=0._wp
  if(present(profile_prefix))then
   if(present(r_limit_m))then
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,a, &
         trim(profile_prefix)//'_inner.csv',match_index=match_index,r_limit_m=r_limit_m)
   else
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,a, &
         trim(profile_prefix)//'_inner.csv',match_index=match_index)
   endif
   if(.not.is_ok(a%status))then;status=a%status;return;endif
   if(present(r_limit_m))then
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,b, &
         trim(profile_prefix)//'_outer.csv',match_index=match_index,surface_start=star,r_limit_m=r_limit_m)
   else
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,b, &
         trim(profile_prefix)//'_outer.csv',match_index=match_index,surface_start=star)
   endif
  else
   if(present(r_limit_m))then
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,a, &
         match_index=match_index,r_limit_m=r_limit_m)
   else
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,a,match_index=match_index)
   endif
   if(.not.is_ok(a%status))then;status=a%status;return;endif
   if(present(r_limit_m))then
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,b, &
         match_index=match_index,surface_start=star,r_limit_m=r_limit_m)
   else
    call integrate_cells(model,pc,p_cut,r0,rtol,step_cap,method,f_hz,star%nu_surface,b, &
         match_index=match_index,surface_start=star)
   endif
  endif
  status=b%status;if(.not.is_ok(status))return
  ! Dimensionless physical W,V, referred to one common surface length.
  za=[a%terminal_w*(a%radius_m/star%radius_m)**3,a%terminal_v*(a%radius_m/star%radius_m)**2]
  zb=[b%terminal_w*(b%radius_m/star%radius_m)**3,b%terminal_v*(b%radius_m/star%radius_m)**2]
  na=maxval(abs(za));nb=maxval(abs(zb))
  if(min(na,nb)<=0._wp.or..not.all(ieee_is_finite([za,zb])))then
   call fail_status(status,NS_NONFINITE_RESULT,'numerics','zero or nonfinite matching direction');return
  endif
  za=za/na;zb=zb/nb;za=za/sqrt(sum(za*za));zb=zb/sqrt(sum(zb*zb))
  residual=za(1)*zb(2)-za(2)*zb(1)
  r_mismatch=a%radius_m-b%radius_m;m_mismatch=a%mass_kg-b%mass_kg
  if(present(step_count))step_count=a%steps+b%steps
  if(present(surface_residual))surface_residual=b%boundary_residual
  call clear_status(status)
 end subroutine
end module
