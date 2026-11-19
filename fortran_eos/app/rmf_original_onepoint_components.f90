! Read-only independent native Fortran evaluation at frozen local mean fields.
! Does NOT receive or recreate inaccessible Pais source meson fields.
! Input:
!   model name (FSU2 or TM1e)
!   n_B, Y_p, n_e, n_mu [fm^-3], phi, omega, rho [MeV].
! Output: name/value lines with all fields, local EOM and physical totals.
program rmf_original_onepoint_components
 use compose_eos_types, only: wp
 use compose_eos_tf_models, only: tf_model_t,load_tf_model,tf_fermion, &
   tf_scalar_potential,tf_vector_potential,tf_hc3,tf_me,tf_mmu
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 type(tf_model_t)::model
 character(32)::name
 character(256)::status
 real(wp)::n,yp,ne,nmu,phi,w,r,nn,np,n3,mstar_n,mstar_p
 real(wp)::en,ep,ee,emu,sn,sp,sd,scalar_lepton,ds_lepton,efn,efp,efe,efmu
 real(wp)::U,du,ddu,Wpot,dw,dr,dww,drr,dwr
 real(wp)::pn,pp,pe,pm,eps_total,P_total
 real(wp)::sigma2,sigma3,sigma4,omega2,rho2,omega4,cross
 read(*,*,end=99,err=99)name
 read(*,*,end=99,err=99)n,yp,ne,nmu,phi,w,r
 if(trim(name)/='FSU2'.and.trim(name)/='TM1e')error stop 'NL only'
 if(.not.all([ieee_is_finite(n),ieee_is_finite(yp),ieee_is_finite(ne), &
      ieee_is_finite(nmu),ieee_is_finite(phi),ieee_is_finite(w), &
      ieee_is_finite(r)]))error stop 'nonfinite input'
 if(n<=0.or.yp<0.or.yp>1.or.ne<0.or.nmu<0)error stop 'bad domain'
 call load_tf_model(name,model,status)
 if(status/='ok')error stop 'unknown RMF source model'
 if(model%density_dependent)error stop 'nonlinear only'
 nn=n*(1-yp);np=n*yp;n3=np-nn
 mstar_n=model%mn-model%gs*phi
 mstar_p=model%mp-model%gs*phi
 if(min(mstar_n,mstar_p)<=0)error stop 'negative mass state'
 call tf_fermion(nn,mstar_n,en,sn,sd,efn)
 call tf_fermion(np,mstar_p,ep,sp,sd,efp)
 call tf_fermion(ne,tf_me,ee,scalar_lepton,ds_lepton,efe)
 call tf_fermion(nmu,tf_mmu,emu,scalar_lepton,ds_lepton,efmu)
 call tf_scalar_potential(model,phi,U,du,ddu)
 call tf_vector_potential(model,w,r,Wpot,dw,dr,dww,drr,dwr)
 pn=nn*efn-en;pp=np*efp-ep
 pe=ne*efe-ee;pm=nmu*efmu-emu
 eps_total=en+ep+ee+emu+(U-Wpot)/tf_hc3+ &
         model%gw*w*n+model%gr*r*n3/2
 P_total=pn+pp+pe+pm+(Wpot-U)/tf_hc3
 sigma2=.5_wp*model%ms**2*phi**2/tf_hc3
 sigma3=model%a3*phi**3/(3*tf_hc3)
 sigma4=model%a4*phi**4/(4*tf_hc3)
 omega2=.5_wp*model%mw**2*w**2/tf_hc3
 rho2=.5_wp*model%mr**2*r**2/tf_hc3
 omega4=model%b4*w**4/(4*tf_hc3)
 cross=model%cross*w*w*r*r/tf_hc3
 call val('epsilon_total',eps_total)
 call val('pressure_total',P_total)
 call val('Mstar_n',mstar_n)
 call val('local_phi',phi)
 call val('local_omega0',w)
 call val('local_rho30',r)
 call val('pressure_kinetic_neutron',pn)
 call val('pressure_kinetic_proton',pp)
 call val('pressure_kinetic_electron',pe)
 call val('pressure_kinetic_muon',pm)
 call val('epsilon_kinetic_neutron',en)
 call val('epsilon_kinetic_proton',ep)
 call val('epsilon_kinetic_electron',ee)
 call val('epsilon_kinetic_muon',emu)
 call val('epsilon_sigma_mass',sigma2)
 call val('epsilon_sigma_cubic',sigma3)
 call val('epsilon_sigma_quartic',sigma4)
 call val('epsilon_omega_mass',omega2)
 call val('epsilon_rho_mass',rho2)
 call val('epsilon_omega_quartic',3*omega4)
 call val('epsilon_omega_rho_mixed',3*cross)
 call val('pressure_sigma_mass',-sigma2)
 call val('pressure_sigma_cubic',-sigma3)
 call val('pressure_sigma_quartic',-sigma4)
 call val('pressure_omega_mass',omega2)
 call val('pressure_rho_mass',rho2)
 call val('pressure_omega_quartic',omega4)
 call val('pressure_omega_rho_mixed',cross)
 call val('scalar_EOM_fm3',du/tf_hc3-model%gs*(sn+sp))
 call val('omega_EOM_fm3',dw/tf_hc3-model%gw*n)
 call val('rho_EOM_fm3',dr/tf_hc3-model%gr*n3/2)
 stop
99 error stop 'invalid input'
contains
 subroutine val(k,x)
 character(*),intent(in)::k
 real(wp),intent(in)::x
 write(*,'(A,1X,ES25.16E3)')trim(k),x
 end subroutine val
end program rmf_original_onepoint_components
