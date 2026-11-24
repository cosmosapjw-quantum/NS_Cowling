! Outward TOV, l=2 static Love and Cowling in geometrized km.
! Input pressure and energy in MeV/fm^3, frequency at infinity in Hz.
module defense_stellar
use iso_c_binding, only: dp=>c_double,c_int
implicit none
real(dp),parameter:: pi=acos(-1._dp), conv=1.3238333135663825e-6_dp, msun=1.4766250385_dp
contains
subroutine star(n,p,e,pc,frequency,offset,limit,res) bind(c)
 integer(c_int),value:: n
 real(dp),intent(in)::p(n),e(n)
 real(dp),value::pc,frequency,offset,limit
 real(dp),intent(out)::res(7)
 real(dp):: y(6),a(6),b(6),c(6),d(6),s,stop,ds,ec,r,m,cslope,ecur,pcur,grad,A1,fact,omega2,el,eh
 real(dp):: compact,surface_y,coef,den,k,oldw,nodes,logbase,logpc
 integer i
 omega2=(2*pi*frequency/299792.458_dp)**2
 res=-1
 if(pc>=p(n).or.pc<=p(1))return
 i=n-1
 do while (i>1)
   if(p(i)<pc)exit
   i=i-1
 end do
 s=1.e-4_dp
 cslope=log(e(i+1)/e(i))/log(p(i+1)/p(i))
 ec=e(i)*(pc/p(i))**cslope*conv
 r=sqrt(2*pc*conv*s**2/(4*pi*(ec/3+pc*conv)*(ec+pc*conv)))
 y=[r,4*pi*ec*r**3/3,2._dp,offset,r**3,-r**2/2]
 nodes=0
 logpc=log(pc)
 do while(i>=1)
  if(p(i+1)==p(i))then
    r=y(1);m=y(2);pcur=p(i)*conv;eh=e(i+1)*conv;el=e(i)*conv
    A1=m+4*pi*r**3*pcur
    y(3)=y(3)-4*pi*r**3*(eh-el)/A1
    if(frequency>0)then
      fact=1-2*m/r;grad=A1/(r*r*fact)
      y(6)=(eh+pcur)/(el+pcur)*y(6)+(eh-el)/(el+pcur)*grad*sqrt(fact)*y(5)*exp(2*y(4))/(r*r*omega2)
    end if
    i=i-1
    cycle
  end if
  stop=sqrt(max(0._dp,logpc-log(p(i))))
  cslope=log(e(i+1)/e(i))/log(p(i+1)/p(i))
  logbase=log(e(i)*conv)+cslope*(logpc-log(p(i)))
  do while(s<stop-1.e-14_dp)
    ds=min(limit,0.09_dp*s,stop-s)
    call rhs(s,y,a)
    call rhs(s+ds/2,y+ds*a/2,b)
    call rhs(s+ds/2,y+ds*b/2,c)
    call rhs(s+ds,y+ds*c,d)
    oldw=y(5)
    y=y+ds*(a+2*b+2*c+d)/6
    if(oldw*y(5)<0)nodes=nodes+1
    s=s+ds
  end do
  i=i-1
 end do
 r=y(1);m=y(2);compact=m/r;surface_y=y(3)
 coef=2-surface_y+2*compact*(surface_y-1)
 den=2*compact*(6-3*surface_y+3*compact*(5*surface_y-8))+4*compact**3*(13-11*surface_y+compact*(3*surface_y-2)+2*compact**2*(1+surface_y))+3*(1-2*compact)**2*coef*log(1-2*compact)
 k=8*compact**5*(1-2*compact)**2*coef/(5*den)
 fact=1-2*compact;grad=(m+4*pi*r**3*p(1)*conv)/(r*r*fact)
 res=[m/msun,r,k,2*k/(3*compact**5),0.5_dp*log(fact)-y(4),omega2*exp(-2*y(4))*y(6)+grad*sqrt(fact)*y(5)/r**2,nodes]
contains
subroutine rhs(x,z,out)
 real(dp),intent(in)::x,z(6)
 real(dp),intent(out)::out(6)
 real(dp)::rad,mass,Pg,Eg,fac,g,dr,area,Fy,Qy,invc,met
 rad=z(1);mass=z(2);Pg=pc*conv*exp(-x*x);Eg=exp(logbase-cslope*x*x)
 fac=1-2*mass/rad;area=4*pi*rad**2;g=(mass+area*rad*Pg)/(rad**2*fac)
 invc=cslope*Eg/Pg;dr=2*x*Pg/((Eg+Pg)*g)
 Fy=(1+area*(Pg-Eg))/fac
 Qy=(area*(5*Eg+9*Pg+(Eg+Pg)*invc)-6)/fac-4*(rad*g)**2
 out(1:4)=dr*[1._dp,area*Eg,-(z(3)**2+Fy*z(3)+Qy)/rad,g]
 out(5:6)=0
 if(frequency>0)then
   met=1/sqrt(fac)
   out(5)=dr*(invc*(omega2*rad**2*met*exp(-2*z(4))*z(6)+g*z(5))-6*met*z(6))
   out(6)=dr*(2*g*z(6)-met*z(5)/rad**2)
 end if
end subroutine
end subroutine
end module

