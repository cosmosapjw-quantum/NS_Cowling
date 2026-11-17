program reference
  implicit none
  double precision n, e, en, rho, plg, chi, gamma, ac, aw, zc, zb
  double precision mun,mup,mue,xrel,ap,an,cp,cn,xnuc,xnucn
  integer stat
  do
    read(*,*,iostat=stat) n
    if(stat/=0) exit
    call BSkEofN(22,n,e,en,rho)
    call BSKfit(22,log10(rho),plg,chi,gamma)
    call INCRST(22,n,ac,aw,zc,zb,mun,mup,mue,xrel,ap,an,cp,cn,xnuc,xnucn)
    write(*,'(12ES25.16)') n,en,rho,10d0**(plg-33.2047d0),chi,ac,aw,zc,zb,mun,mue,(aw-ac)/aw
  end do
end program reference
