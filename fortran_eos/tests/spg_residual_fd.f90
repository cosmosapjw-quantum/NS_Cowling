program spg_residual_fd
  use compose_eos_types, only: wp
  use compose_eos_gdfm_native, only: gdfm_native_model_t, load_gdfm_native, &
       gdfm_native_baryons, native_fermion
  implicit none
  real(wp), parameter :: pi=acos(-1._wp), e2=1.43996448_wp, me=.51099895_wp
  character(16) :: model_id
  character(160) :: status
  integer :: row, ios, k
  real(wp) :: nb, ye, aa, zz, ni, ng, f, radius, sigma0, bfit, y
  real(wp) :: eps_i, mn_i, mp_i, p_i, eps_g, mn_g, mp_g, p_g
  real(wp) :: ee, pe, mue, dval, dp, sigma, dsigma, phi, ccell
  real(wp) :: analytic, scaled, fd, step, ep, em
  real(wp), parameter :: steps(5)=[1.e-3_wp,3.e-4_wp,1.e-4_wp,3.e-5_wp,1.e-5_wp]
  type(gdfm_native_model_t) :: model
  write(*,'(a)') 'model,row,step,beta_fd_MeV,beta_full_MeV,beta_scaled_MeV,abs_fd_full_MeV'
  do
    read(*,*,iostat=ios) model_id,row,nb,ye,aa,zz,ni,ng,f,radius,sigma0,bfit
    if (ios < 0) exit
    if (ios > 0) stop 4
    call load_gdfm_native(trim(model_id),model,status)
    if (trim(status)/='ok') then
      write(0,'(a,1x,a,1x,i0,1x,a)') 'FAIL load',trim(model_id),row,trim(status)
      stop 20
    endif
    y=zz/aa
    call gdfm_native_baryons(model,ni,y,eps_i,mn_i,mp_i,p_i,status)
    if (trim(status)/='ok') then
      write(0,'(a,1x,a,1x,i0,1x,a)') 'FAIL dense reference',trim(model_id),row,trim(status)
      stop 21
    endif
    call gdfm_native_baryons(model,ng,0._wp,eps_g,mn_g,mp_g,p_g,status)
    if (trim(status)/='ok') then
      write(0,'(a,1x,a,1x,i0,1x,a)') 'FAIL gas reference',trim(model_id),row,trim(status)
      stop 22
    endif
    call native_fermion(nb*ye,me,ee,pe,mue)
    dval=y**(-3)+(1-y)**(-3)
    dp=-3*y**(-4)+3*(1-y)**(-4)
    sigma=sigma0*(bfit+16)/(bfit+dval)
    dsigma=-sigma*dp/(bfit+dval)
    phi=((2._wp-3._wp*f**(1._wp/3)+f)/5._wp)
    ccell=2*pi*e2*f*phi*radius**2*(ni*y)**2
    analytic=(mp_i-mn_i)+mue+2*ccell/(f*ni*y)+3*dsigma/(radius*ni)
    scaled=analytic-f*3*dsigma/(radius*ni)
    do k=1,size(steps)
      step=steps(k)
      ep=energy(y+step)
      em=energy(y-step)
      fd=(ep-em)/(2*step)/(f*ni)
      write(*,'(a,a,i0,a,5(es24.15e3,a))') trim(model_id),',',row,',',step,',',fd,',',analytic,',',scaled,',',abs(fd-analytic),''
    enddo
  enddo
contains
  function energy(yarg) result(e)
    real(wp), intent(in) :: yarg
    real(wp) :: e, ei, mni, mpi, pp, eg, mng, mpg, pg, el, pel, mul, d, s, ph, coul
    character(160) :: stat
    call gdfm_native_baryons(model,ni,yarg,ei,mni,mpi,pp,stat)
    if (trim(stat)/='ok') then
      write(0,'(a,1x,a,1x,i0,1x,a,1x,es24.15e3,1x,a)') 'FAIL dense perturbation',trim(model_id),row, &
        'y=',yarg,trim(stat)
      stop 23
    endif
    call gdfm_native_baryons(model,ng,0._wp,eg,mng,mpg,pg,stat)
    if (trim(stat)/='ok') then
      write(0,'(a,1x,a,1x,i0,1x,a)') 'FAIL gas perturbation',trim(model_id),row,trim(stat)
      stop 24
    endif
    call native_fermion(f*ni*yarg,me,el,pel,mul)
    d=yarg**(-3)+(1-yarg)**(-3)
    s=sigma0*(bfit+16)/(bfit+d)
    ph=(2._wp-3._wp*f**(1._wp/3)+f)/5._wp
    coul=2*pi*e2*f*ph*radius**2*(ni*yarg)**2
    e=f*ei+(1-f)*eg+el+coul+3*f*s/radius
  end function energy
end program
