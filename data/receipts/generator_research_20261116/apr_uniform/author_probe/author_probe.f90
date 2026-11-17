program author_probe
use Kind_Types_Mod, only: DP,I4B
use Physical_Constants_Mod
use Global_Variables_Mod, only: LDP_ONLY
use Skyrme_Coefficients_Mod
use Skyrme_Potential_Mod
use Skyrme_Effective_Mass_Mod
use Skyrme_Energy_Mod
use Skyrme_Force_Mod
implicit none
real(DP) :: nn,np,n,y,delta,g1,g2,f1,f2,u,un,up
real(DP) :: auxmn,auxmp,ln,lp,men,mep,dnn,dnp,dpn,dpp
real(DP) :: taun,taup,energy,vn,vp,mun,mup,pressure
integer(I4B) :: flag,ios
logical :: high
cp1=337.2_DP;cp2=-382._DP;cp3=89.8_DP;cp4=.457_DP;cp5=-59._DP
cp6=-19.1_DP;cp7=214.6_DP;cp8=-384._DP;cp9=6.4_DP;cp10=69._DP
cp11=-33._DP;cp12=.35_DP;cp13=0._DP;cp14=0._DP;cp15=287._DP
cp16=-1.54_DP;cp17=175._DP;cp18=-1.45_DP;cp19=.32_DP;cp20=.195_DP;cp21=0._DP
Mass_n=neutron_mass_MeV
Mass_p=proton_mass_MeV
Neut_Prot_Mass_Diff=Mass_n-Mass_p
LDP_ONLY=.false.
do
 read(*,*,iostat=ios) flag,nn,np
 if(ios/=0) exit
 n=nn+np;y=np/n;delta=(nn-np)/n
 call SKYRME_POTENTIAL(n,delta,nn,np,g1,g2,f1,f2,u,un,up,flag,high)
 call SKYRME_EFFECTIVE_MASS(nn,np,n,y,0._DP,auxmn,auxmp,ln,lp,men,mep,dnn,dnp,dpn,dpp,1)
 taun=C_tau*nn**(5._DP/3._DP);taup=C_tau*np**(5._DP/3._DP)
 call SKYRME_ENERGY(n,y,0._DP,taun,taup,ln,lp,u,energy)
 call SKYRME_FORCE(taun,taup,dnn,dnp,dpn,dpp,un,up,vn,vp)
 energy=energy+Mass_n*n
 mun=ln*(3._DP*PI*PI*nn)**(2._DP/3._DP)+vn+Mass_n
 mup=lp*(3._DP*PI*PI*np)**(2._DP/3._DP)+vp+Mass_n
 pressure=nn*mun+np*mup-energy
 write(*,'(I3,8ES26.17)') flag,nn,np,energy,pressure,mun,mup,men,mep
enddo
end program author_probe
