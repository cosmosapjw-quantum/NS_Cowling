! Constrained semiclassical evaluation of the native Negele--Vautherin EDF.
! The density and kinetic-density closure are explicit reconstruction choices.
program nv_profile_generate
 use compose_eos_types, only: wp, eos_data_t, allocate_eos
 use compose_eos_nv_original, only: nv_original_source_t, load_nv_original_source
 use compose_eos_nv_uniform_branch, only: write_nv_uniform_branch
 use compose_eos_nv_hf, only: nv_hf_cell_energy, nv_hf_gradient, nv_hf_fields, nv_hf_electrons
 use compose_eos_io, only: write_compose, write_metadata, read_compose
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 real(wp), parameter :: pi=acos(-1._wp), tf=(3._wp/5)*(3*pi*pi)**(2._wp/3)
 integer, parameter :: rows=11
 type(nv_original_source_t) :: src
 type(eos_data_t) :: eos, back
 integer :: mesh, count_nb, i, j, mode, anchor_file, scan_file, last, ios, a, n, z
 real(wp) :: nb, eps, p, e, lam(2), perr, nrmax, offset, low, high, mul, mue, muerr, previous(2), pressure_fd_error
 real(wp), allocatable :: r(:), w(:), seed(:,:), jp(:)
 real(wp) :: h
 character(1024) :: path, out, arg, status
 character(8), parameter :: labels(0:1)=[character(8)::'printed','d03x100']
 character(*), parameter :: header='nb_fm3,N,Z,epsilon_MeV_fm3,pressure_MeV_fm3,E_minus_mn_MeV,'// &
  'neutron_contrast,proton_contrast,max_nucleon_density_fm3,max_particle_error'
 if(command_argument_count()/=4) error stop 'expected SOURCE OUT_DIRECTORY MESH N_ROWS'
 call get_command_argument(1,path)
 call get_command_argument(2,out)
 call get_command_argument(3,arg)
 read(arg,*,iostat=ios)mesh
 if(ios/=0.or.mesh<32) error stop 'need at least 32 radial cells'
 call get_command_argument(4,arg)
 read(arg,*,iostat=ios)count_nb
 if(ios/=0.or.count_nb<4) error stop 'need at least 4 density samples'
 allocate(r(mesh),w(mesh),seed(mesh,2),jp(mesh))
 jp=0
 call load_nv_original_source(trim(path),src,status)
 if(status/='ok')then
  print *,trim(status)
  error stop 1
 endif
 call write_nv_uniform_branch(src,trim(out),status)
 if(status/='ok')then
  print *,trim(status)
  error stop 1
 endif
 low=src%fit%table(1,rows-1)*1e-39_wp
 high=src%fit%table(1,rows)*1e-39_wp
 do mode=0,1
  open(newunit=anchor_file,file=trim(out)//'/'//trim(labels(mode))//'_anchors.csv',status='replace')
  write(anchor_file,'(a)')header//',source_ground_E_minus_mn_MeV'
  do j=1,rows
   nb=src%fit%table(1,j)*1e-39_wp
   call setup_cell(j,nb)
   call minimize(e,p,lam,perr,nrmax)
   offset=e/a
   eps=nb*(src%fit%mn+offset)
   write(anchor_file,'(es24.15e3,2(",",i0),*(:,",",es24.15e3))') &
    nb,n,z,eps,p,offset,lam,nrmax,perr,src%fit%table(8,j)
  enddo
  close(anchor_file)
  call allocate_eos(eos,count_nb)
  eos%mn=src%fit%mn
  eos%mp=src%fit%mp
  eos%charge_fraction_grid=real(z,wp)/a
  eos%beta_equilibrated=.false.
  eos%nT=1
  eos%nYq=1
  eos%t_first=1
  eos%t_last=1
  eos%yq_first=1
  eos%yq_last=1
  eos%selected_iT=1
  eos%selected_iYq=1
  eos%composition_file_present=.true.
  eos%composition_available=.true.
  eos%ye=eos%charge_fraction_grid
  eos%ye_present=.true.
  eos%row_status=1
  open(newunit=scan_file,file=trim(out)//'/'//trim(labels(mode))//'.csv',status='replace')
  write(scan_file,'(a)')header//',mu_L_MeV,mu_e_eff_MeV,mu_L_finite_difference_error_MeV,profile_sector,dilation_P_rel_error'
  previous=[1._wp,1._wp]
  do i=1,count_nb
   nb=low+(high-low)*(i-1)/(count_nb-1)
   call setup_cell(rows,nb)
   call minimize(e,p,lam,perr,nrmax,previous)
   previous=lam
   eos%iphase(i)=merge(1,2,maxval(lam)>3e-6_wp)
   offset=e/a
   eps=nb*(eos%mn+offset)
   call charge_response(lam,mul,mue,muerr)
   call check_dilation(lam,p,pressure_fd_error)
   write(scan_file,'(es24.15e3,2(",",i0),10(",",es24.15e3),",",i0,",",es24.15e3)') &
    nb,n,z,eps,p,offset,lam,nrmax,perr,mul,mue,muerr,eos%iphase(i),pressure_fd_error
   eos%nb(i)=nb
   eos%q(1,i)=p/nb
   ! Euler decomposition of the fixed-A cell family at n,Y_e=Z/A.
   eos%q(3,i)=(offset+p/nb-eos%charge_fraction_grid*mul)/eos%mn
   eos%q(4,i)=(mul-mue)/eos%mn
   eos%q(5,i)=mul/eos%mn
   eos%q(6,i)=offset/eos%mn
   eos%q(7,i)=offset/eos%mn
   if(i>1)then
    last=i-1
    if(p<=eos%q(1,last)*eos%nb(last))error stop 'unstable pressure branch'
    if(eps<=(1+eos%q(7,last))*eos%mn*eos%nb(last))error stop 'nonmonotone energy branch'
    if(p-eos%q(1,last)*eos%nb(last)>=eps-(1+eos%q(7,last))*eos%mn*eos%nb(last))error stop 'acausal branch'
   endif
  enddo
  close(scan_file)
  call write_compose(trim(out)//'/'//trim(labels(mode)),eos,status)
  if(status=='ok')call write_metadata(trim(out)//'/'//trim(labels(mode))//'/metadata.nml',eos,status)
  if(status=='ok')call read_compose(trim(out)//'/'//trim(labels(mode)),back,status)
  if(status=='ok')then
   if(back%beta_equilibrated.or.any(back%iphase/=eos%iphase).or.any(abs(back%q-eos%q)>1e-13_wp)) &
     status='profile round-trip mismatch'
  endif
  if(status/='ok')then
   print *,trim(status)
   error stop 1
  endif
 enddo
contains
 subroutine setup_cell(jref,density)
  integer,intent(in)::jref
  real(wp),intent(in)::density
  real(wp):: anchor_nb, gas, bound, v0, rad0, rc, s(mesh), volume, radius, smean
  integer::k
  n=nint(src%fit%table(2,jref))
  z=nint(src%fit%table(3,jref))
  a=n+z
  anchor_nb=src%fit%table(1,jref)*1e-39_wp
  gas=src%fit%table(6,jref)
  v0=a/anchor_nb
  bound=n-gas*v0
  if(bound<=0.or.bound>=n)error stop 'invalid Table 3 gas fraction for this profile'
  rad0=(3*v0/(4*pi))**(1._wp/3)
  ! Radius and diffuseness use the pre-existing solver's Fermi seed convention.
  rc=(3*(bound+z)/(4*pi*.14_wp))**(1._wp/3)
  volume=a/density
  radius=(3*volume/(4*pi))**(1._wp/3)
  h=radius/mesh
  do k=1,mesh
   r(k)=(k-.5_wp)*h
   w(k)=((real(k,wp)**3-real(k-1,wp)**3)/3)*h**3
   s(k)=1/(1+exp(((k-.5_wp)*rad0/mesh-rc)/.65_wp))
  enddo
  smean=sum(w*s)/sum(w)
  seed(:,1)=(n-bound+bound*s/smean)/a
  seed(:,2)=(real(z,wp)/a)*(s/smean)
 end subroutine

 subroutine evaluate(contrast,energy,pressure,particle_error,maxrho,valid)
  real(wp),intent(in)::contrast(2)
  real(wp),intent(out)::energy,pressure,particle_error,maxrho
  logical,intent(out)::valid
  real(wp)::rho(mesh,2),tau(mesh,2),gr(mesh),numbers(2)
  integer::s
  character(256)::msg
  numbers=real([n,z],wp)
  do s=1,2
   rho(:,s)=nb*((1-contrast(s))*numbers(s)/a+contrast(s)*seed(:,s))
   call nv_hf_gradient(rho(:,s),h,gr)
   tau(:,s)=tf*rho(:,s)**(5._wp/3)+gr**2/(36*rho(:,s))
  enddo
  maxrho=maxval(rho)
  valid=maxrho<=.15_wp
  if(.not.valid)return
  particle_error=maxval(abs(matmul(4*pi*w,rho)-numbers))
  if(particle_error>1e-8_wp)error stop 'particle-number closure failed'
  call nv_hf_cell_energy(src,r,w,h,rho,tau,jp,z,energy,pressure,msg,mode)
  if(msg/='ok')then
   print *,trim(msg)
   error stop 'invalid functional at an admissible trial'
  endif
  if(.not.all(ieee_is_finite([energy,pressure])))error stop 'nonfinite EDF'
 end subroutine

 subroutine check_dilation(contrast,pressure,rel_error)
  real(wp),intent(in)::contrast(2),pressure
  real(wp),intent(out)::rel_error
  real(wp),parameter::delta=1e-5_wp
  real(wp)::nb0,es(-1:1),pp,err,rhomax
  integer::k
  logical::valid
  nb0=nb
  do k=-1,1,2
   nb=nb0*(1+k*delta)
   call setup_cell(rows,nb)
   call evaluate(contrast,es(k),pp,err,rhomax,valid)
   if(.not.valid)error stop 'invalid compression probe'
  enddo
  nb=nb0
  call setup_cell(rows,nb)
  rel_error=abs(nb*(es(1)-es(-1))/(2*delta*a)-pressure)/pressure
  if(rel_error>2e-7_wp)error stop 'dilation pressure mismatch'
 end subroutine

 subroutine charge_response(contrast,mu_l,mu_e,fd_error)
  ! Extend each normalized species profile to continuous Z at fixed A,V and
  ! fixed shape. At a contrast extremum this gives its envelope derivative.
  real(wp),intent(in)::contrast(2)
  real(wp),intent(out)::mu_l,mu_e,fd_error
  real(wp),parameter::dz=1.e-3_wp
  real(wp)::rho(mesh,2),tau(mesh,2),tfpart(mesh,2),u(mesh,2),b(mesh,2),ss(mesh)
  real(wp)::gr(mesh),num(2),e(-1:1),ez,vol,en,ed,ex,ee,pe,mu,ec,p_here
  integer::k,s
  character(256)::msg
  vol=sum(4*pi*w)
  do k=-1,1
   ez=z+k*dz
   num=[a-ez,ez]
   do s=1,2
    rho(:,s)=nb*((1-contrast(s))*real(merge(n,z,s==1),wp)/a+contrast(s)*seed(:,s)) &
             *num(s)/real(merge(n,z,s==1),wp)
    call nv_hf_gradient(rho(:,s),h,gr)
    tfpart(:,s)=tf*rho(:,s)**(5._wp/3)
    tau(:,s)=tfpart(:,s)+gr*gr/(36*rho(:,s))
   enddo
   call nv_hf_fields(src,r,w,h,rho,tau,jp,ez/vol,u,b,ss,en,ed,ex,p_here,ec,msg,mode)
   if(msg/='ok')error stop 'invalid charge probe'
   call nv_hf_electrons(ez/vol,ee,pe,mu)
   e(k)=en+ed+ex+vol*ee+ez*(src%fit%mp-src%fit%mn)
   if(k==0)then
    mu_e=mu+ec
    mu_l=mu_e+src%fit%mp-src%fit%mn
    do s=1,2
     mu_l=mu_l+merge(-1,1,s==1)*sum(4*pi*w*(u(:,s)*rho(:,s)+ &
       b(:,s)*(tau(:,s)+(2._wp/3)*tfpart(:,s))))/num(s)
    enddo
   endif
  enddo
  fd_error=abs(mu_l-(e(1)-e(-1))/(2*dz))
  if(fd_error>2e-5_wp)error stop 'charge derivative mismatch'
 end subroutine

 subroutine minimize(best_energy,best_pressure,best,particle_error,maxrho,previous_best)
  real(wp),intent(out)::best_energy,best_pressure,best(2),particle_error,maxrho
  real(wp),intent(in),optional::previous_best(2)
  integer,parameter::cells=8, refinements=17
  integer::ix,iy,lev,pass,xi,yi,iseed,nseeds
  real(wp)::step,base(2),trial(2),energy,pressure,err,rhomax,grid(0:cells,0:cells)
  real(wp)::local_best(2),local_energy,starts(2,(cells+1)**2+1)
  logical::valid,improved
  best_energy=huge(1._wp)
  best=0
  grid=huge(1._wp)
  do ix=0,cells
   do iy=0,cells
    trial=real([ix,iy],wp)/cells
    call evaluate(trial,energy,pressure,err,rhomax,valid)
    if(valid)grid(ix,iy)=energy
   enddo
  enddo
  ! Refine every discrete local minimum, including competing uniform/clustered profiles.
  ! The two-shell stencil resolves an oblique valley leaving the uniform boundary.
  nseeds=0
  do ix=0,cells
   do iy=0,cells
    local_energy=grid(ix,iy)
    if(local_energy==huge(1._wp))cycle
    if(local_energy>minval(grid(max(ix-1,0):min(ix+1,cells),max(iy-1,0):min(iy+1,cells))))cycle
    nseeds=nseeds+1
    starts(:,nseeds)=real([ix,iy],wp)/cells
   enddo
  enddo
  if(present(previous_best))then
   nseeds=nseeds+1
   starts(:,nseeds)=previous_best
  endif
  do iseed=1,nseeds
    local_best=starts(:,iseed)
    call evaluate(local_best,local_energy,pressure,err,rhomax,valid)
    if(.not.valid)cycle
    step=1._wp/cells
    do lev=1,refinements
     step=step/2
     do pass=1,2*cells+1
      base=local_best
      improved=.false.
      do xi=-2,2
       do yi=-2,2
        if(xi==0.and.yi==0)cycle
        trial=base+step*[xi,yi]
        if(any(trial<0).or.any(trial>1))cycle
        call evaluate(trial,energy,pressure,err,rhomax,valid)
        if(.not.valid)cycle
        if(energy<local_energy)then
         local_energy=energy
         local_best=trial
         improved=.true.
        endif
       enddo
      enddo
      if(.not.improved)exit
     enddo
     if(pass>2*cells+1)error stop 'profile search stalled'
    enddo
    if(local_energy<best_energy)then
     best_energy=local_energy
     best=local_best
    endif
  enddo
  if(best_energy==huge(1._wp))error stop 'no admissible profile'
  call evaluate(best,best_energy,best_pressure,particle_error,maxrho,valid)
  if(.not.valid)error stop 'invalid minimum'
 end subroutine
end program
