! Equation-level audit of the printed NV1972 B/C/D tables.
! All arithmetic of the source-coefficient reconstruction is native Fortran.
program test_nv_dme_coefficient_consistency
 use compose_eos_types,only:wp
 use compose_eos_nv_dme,only:nv_dme_table1_coefficients
 implicit none
 real(wp)::b(0:7,0:7),c(0:6,0:7),d(0:6,0:6),dr(0:6,0:6)
 real(wp)::lo(0:13),hi(0:13),g(0:6),gc,rad,bpart,ctilde,scale,dlo,dhi,printed_radius,g0rad
 real(wp)::gg(0:6,0:6),gglo(0:6,0:6),gghi(0:6,0:6),gc2,rad2
 integer::i,j,k,bc,violations,interval_failures,outside,general_outside,symmetry_failures
 call nv_dme_table1_coefficients(b,c,d)
 lo=-huge(1._wp);hi=huge(1._wp);violations=0;interval_failures=0
 ! G(t)=V_D^U(t/2), F=x*L(x)+y*G(x+y)-B(x,y)/4.
 ! C_ij -(i+1)*B_(i+1,j)/4 = -j*choose(i+j,i)*g_(i+j), j>=1.
 do i=0,6
  do j=1,7
   k=i+j;scale=real(j*choose(k,i),wp)
   ctilde=c(i,j)-real(i+1,wp)*b(i+1,j)/4
   rad=(rounding(c(i,j))+real(i+1,wp)*rounding(b(i+1,j))/4)/scale
   gc=-ctilde/scale
   lo(k)=max(lo(k),gc-rad);hi(k)=min(hi(k),gc+rad)
  enddo
 enddo
 g(0)=(b(0,1)/2-d(0,0))/4
 g0rad=(rounding(b(0,1))/2+rounding(d(0,0)))/4
 lo(0)=g(0)-g0rad;hi(0)=g(0)+g0rad
 write(*,'(a)')'G_INTERVALS k,lower,upper,intersection_nonempty'
 do k=0,13
  write(*,'(i3,2(",",es25.16e3),",",l1)')k,lo(k),hi(k),lo(k)<=hi(k)
  if(lo(k)>hi(k))interval_failures=interval_failures+1
  if(k>=1.and.k<=6)g(k)=(lo(k)+hi(k))/2
  if(k>6)then
   if(lo(k)>0.or.hi(k)<0)violations=violations+1
  endif
 enddo
 ! General symmetric G(n,p) audit isolates functional identities separately from
 ! the stronger source assumption G(n,p)=U((n+p)/2).
 gg=0;gglo=0;gghi=0;general_outside=0;symmetry_failures=0;gg(0,0)=g(0);gglo(0,0)=lo(0);gghi(0,0)=hi(0)
 write(*,'(a)')'GENERAL_G i,j,Gij,symmetry_difference,rounding_sum'
 do i=1,6
  do j=0,6-i
   gg(i,j)=-(c(i-1,j+1)-real(i,wp)*b(i,j+1)/4)/real(i,wp)
   rad=rounding(c(i-1,j+1))/real(i,wp)+rounding(b(i,j+1))/4
   gglo(i,j)=gg(i,j)-rad;gghi(i,j)=gg(i,j)+rad
   if(j>0)then
    gc2=-(c(j-1,i+1)-real(j,wp)*b(j,i+1)/4)/real(j,wp)
    rad2=rounding(c(j-1,i+1))/real(j,wp)+rounding(b(j,i+1))/4
    write(*,'(i2,",",i2,3(",",es25.16e3))')i,j,gg(i,j),gg(i,j)-gc2,rad+rad2
    if(abs(gg(i,j)-gc2)>rad+rad2)symmetry_failures=symmetry_failures+1
   endif
  enddo
 enddo
 gg(0,1:6)=gg(1:6,0);gglo(0,1:6)=gglo(1:6,0);gghi(0,1:6)=gghi(1:6,0)
 write(*,'(a)')'GENERAL_G_D i,j,printed,reconstructed,discrepancy,rounding_bound'
 do i=0,6
  do j=i,6
   k=i+j
   bpart=(real(j+1,wp)*b(i,j+1)+real(i+1,wp)*b(j,i+1))/4
   gc=bpart-real(4+k,wp)*gg(i,j)
   rad=(real(j+1,wp)*rounding(b(i,j+1))+real(i+1,wp)*rounding(b(j,i+1)))/4
   rad=rad+real(4+k,wp)*(gghi(i,j)-gglo(i,j))/2+rounding(d(i,j))
   write(*,'(i2,",",i2,4(",",es25.16e3))')i,j,d(i,j),gc,d(i,j)-gc,rad
   if(abs(d(i,j)-gc)>rad)then
    general_outside=general_outside+1
    if(i/=0.or.j/=3)error stop 'unexpected D inconsistency beyond D03'
    if(abs(-6.62699e5_wp-gc)>rad-rounding(d(i,j))+rounding(-6.62699e5_wp))error stop 'D03x100 diagnostic outside derived interval'
   endif
  enddo
 enddo
 write(*,'(a,i0)')'general_G_symmetry_interval_failures=',symmetry_failures
 write(*,'(a,i0)')'general_G_D_rounding_failures=',general_outside
 if(interval_failures/=3.or.violations/=0)error stop 'unexpected average-density source inconsistency pattern'
 if(symmetry_failures/=0.or.general_outside/=1)error stop 'unexpected general-G source inconsistency pattern'
 write(*,'(a)')'NV_DME_EXPECTED_SOURCE_INCONSISTENCIES_CONFIRMED'
 write(*,'(a)')'D03x100_DIAGNOSTIC_IS_NOT_AN_AUTHENTICATED_ERRATUM'
contains
 pure integer function choose(n,k) result(v)
  integer,intent(in)::n,k
  integer::a
  v=1
  do a=1,min(k,n-k)
   v=v*(n-a+1)/a
  enddo
 end function
 pure real(wp) function rounding(x) result(v)
  real(wp),intent(in)::x
  integer::exponent
  if(x==0)then;v=0;return;endif
  ! Every nonzero B/C/D coefficient is printed with six significant digits.
  exponent=floor(log10(abs(x)))
  v=.5_wp*10._wp**(exponent-5)
 end function
end program test_nv_dme_coefficient_consistency
