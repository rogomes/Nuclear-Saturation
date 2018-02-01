! Programa desenvolvido por Rosana Gomes (janeiro/2014)
! O programa calcula a EoS e populaçao da materia de hadrons INCLUINDO O MESON DELTA sob efeito de campos magneticos e MMA

! modificado em 14/04/2014

!inclusao de pperp2  e mag2 (veronica), para particulas carregadas e nao carregadas
!para conferir valores com a pressao perperndicular total

!******* COMPILAR  ***********
!gfortran -I./ -c zen_coupling.f90
!gfortran zen_coupling.o -L./ -lnrf90


module global

implicit none
double precision :: mu_n, mu_e, Bmag, rhob, fs, sig, ome, r03, del3, gsn, grn, gwn, gdn
double precision, dimension(10), parameter :: &
gama=(/ 2./2., 2./2., 2./2., 2./2., 2./2., 2./2., 2./2., 2./2., 2./2., 2./2. /), & 
m=(/ 939.6/197.327, 938.3/197.327, 1116./197.327, 1189./197.327, 1193./197.327, &
1197./197.327, 1315./197.327, 1321./197.327, .511/197.327, 105.66/197.327 /), &
q= (/ 1., 0., 0., 1., 0., -1., 0., -1., -1., -1. /), & 
qb= (/ 1., 1., 1., 1., 1., 1., 1., 1., 0., 0. /), & 
qs= (/ 0., 0., -1., -1., -1., -1., -2., -2., 0., 0. /), & 
i3= (/ .5, -.5, 0., 1., 0., -1., .5, -.5, 0., 0. /), &

!Magneton em unidades (1/MeV)
mun=(/ 4.55213958d-5*197.327d0, 4.55213958d-5*197.327d0,&
4.55213958d-5*197.327d0,4.55213958d-5*197.327d0, &
4.55213958d-5*197.327d0,4.55213958d-5*197.327d0, &
4.55213958d-5*197.327d0,4.55213958d-5*197.327d0, &
8.3585617d-2*197.327d0, 0.d0 /)


double precision, parameter :: e = 0.0854245005776     !para Bmag em MeV**2
double precision, parameter :: ms = 550.d0/197.327d0   !massa sigma
double precision, parameter :: mw = 783.d0/197.327d0   !massa omega
double precision, parameter :: mr = 769.d0/197.327d0   !massa rho
double precision, parameter :: md = 980.d0/197.327d0   !massa delta
double precision, parameter :: mn = 939.d0/197.327d0   !massa nucleon
integer :: HYS, iB2, AMM
integer :: y, nu, nnu, counter

integer, DIMENSION(10,-1:1) :: iB 

integer, dimension(10), parameter :: spini = (/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  /), &
spinf = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), &
spins=(/ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2/)

double precision, DIMENSION (10) :: gr, gw, gs, gd, meff, mu, ef, kb
double precision, DIMENSION (10,-1:1) :: kf, rho, mspin, nu_max, argument
double precision, DIMENSION (10,0:200000,-1:1) :: kfmag, mmag, mmag2


!***************************************************************************************************************************************
!Constantes de acoplamento encontradas no outro programa, spin dos barions
!as = (gs/ms)**2,   aw = (gw/mw)**2,   ar = (gr/mr)**2 , ad=(gd/md)**2

!A1 (rho_0 = 0.153 fm^-3, E_L= -15.75 MeV)
double precision, parameter :: as = 13.04d0, aw = 7.27d0, ar = 5.31d0, ad = 2.55d0 
double precision, parameter :: lambda = 0.060d0, w= 0.695d0, rhosat = 0.153d0
!A5
!double precision, parameter :: w = 0.775, lambda = 0.131, as = 10.54, aw = 4.97, ar = 6.68, ad = 4.68

!A2 (rho_0 = 0.16 fm^-3, E_L= -16.25 MeV)
!double precision, parameter :: as = 12.56d0, aw = 6.87d0, ar = 5.03d0, ad = 2.58d0 
!double precision, parameter :: lambda = 0.063d0, w= 0.696d0, rhosat = 0.16d0
!A6
!double precision, parameter :: w = 0.775, lambda = 0.138, as = 10.17, aw = 4.69, ar = 6.30, ad = 4.59


!A3 (rho_0 = 0.16 fm^-3, E_L= -15.75 MeV)
!double precision, parameter :: as = 12.48d0, aw = 6.90d0, ar = 5.01d0, ad = 2.53d0 
!double precision, parameter :: lambda = 0.060d0, w= 0.695d0, rhosat = 0.16d0

!A4 (rho_0 = 0.153 fm^-3, E_L= -16.25 MeV)
!double precision, parameter :: as = 13.11d0, aw = 7.23d0, ar = 5.33d0, ad = 2.60d0 
!double precision, parameter :: lambda = 0.063d0, w= 0.696d0, rhosat = 0.153d0

end module global

!******************************************************************************************

program eos_hmma
USE nrtype; USE global
implicit none
logical(LGT) :: check
integer :: i
real(SP), allocatable :: x(:), xs(:), xw(:), xr(:), xd(:) 

double precision :: sigma1, omega1, m1
double precision :: elpn_lamb, xsl, meffn 
double precision :: potential_n, potential_lamb, potential_sig, potential_casc

double precision :: de, p, B_a, B_b, Bsurf, Bc
double precision :: p_perp, p_par, de_mag, pperp2, mag2,pmeson, mag, p_charged, p_perp0, p_par0
INTERFACE
	SUBROUTINE broydn(x,check)
	USE nrtype
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	END SUBROUTINE broydn
END INTERFACE


allocate(xs(8),xw(8),xr(8),xd(8))  
allocate(x(5)) 


write(*,*) 'iB= 0 (sem campo magnetico), 1 (com campo magnetico), 2 (campo variavel)'
!read(*,*) iB2
iB2 = 0
write(*,*) iB2
write(*,*) 'HYS = 0 (sem hyperons), 1 (Universal), 2 (Moszkowski), 3 (SU(6), 4 (hypernucleous)'
!read(*,*) HYS
HYS = 3
write(*,*) HYS
write(*,*) 'AMM= 0 (sem momentum magnetico anomalo), AMM= 1 (com momentum magnetico anomalo)'
!read(*,*) AMM
AMM = 0
write(*,*) AMM

!**********************************
!*    Acoplamento de hiperons      *    
!**********************************

!Acoplamento nucleons 
gsn = sqrt(as*((ms)**2)) 
gwn = sqrt(aw*((mw)**2))
grn = sqrt(ar*((mr)**2))
gdn = sqrt(ad*((md)**2))

do i=1,2
    xs(i)= 1.d0
    xw(i)= 1.d0
    xr(i)= 1.d0
    xd(i)= 1.d0
enddo


if(HYS.eq.0) then
do i=3,8
    xs(i)= 0.d0
    xw(i)= 0.d0
    xr(i)= 0.d0
    xd(i)= 0.d0	
enddo


else if(HYS.eq.1) then ! Universal 
do i=3,8
    xs(i)= 1.d0
    xw(i)= 1.d0
    xr(i)= 1.d0
    xd(i)= 1.d0
enddo


else if(HYS.eq.2) then  !Moszkowski 
do i=3,8
    xs(i)= sqrt(2.d0/3.d0)
    xw(i)= sqrt(2.d0/3.d0)
    xr(i)= sqrt(2.d0/3.d0)
    xd(i)= sqrt(2.d0/3.d0)
enddo

else if(HYS.eq.3) then !SU(6)  
do i=3,6  
    xw(i)= 2.d0/3.d0          !Omega coupling
enddo

do i=7,8
    xw(i)= 1.d0/3.d0               !Omega coupling
    xr(i)= 1.d0                    !Rho coupling
    xd(i)= 1.d0                    !Delta coupling
enddo

xr(3)= 0.d0     !Rho coupling
xd(3)= 0.d0     !Delta coupling
do i=4,6 
    xr(i)= 2.d0    !Rho coupling
    xd(i)= 2.d0    !Delta coupling
enddo

! Hyperon potentials (J. Schaffner-Bielich, Nuclear Physics A 881 (2012) 62-77)
potential_lamb = - 28.d0/197.327d0 
potential_sig = - 30.d0/197.327d0 
potential_casc = 40.d0/197.327d0 

meffn = w*mn

!Meson fields at saturation
sigma1 = (mn-meffn)/gsn
omega1 = (gwn/(mw**2))*rhosat     

! Nucleon potential
potential_n = (gwn*omega1 - gsn*sigma1)*197.327d0    !MeV

write(*,*) potential_n

! Lambda-sigma coupling
xs(3) = (xw(3)*gwn*omega1 - potential_lamb)/(gsn*sigma1)

!Sigma-sigma baryon coupling
do i=4,6           
xs(i) = (xw(i)*gwn*omega1 - potential_sig)/(gsn*sigma1)
enddo

! Cascade-sigma baryon coupling
do i=7,8           
xs(i) = (xw(i)*gwn*omega1 - potential_casc)/(gsn*sigma1)
enddo


endif !end of the HYS loop 


!Hyperons coupling
do i = 1, 8
  gw(i)= xw(i)*gwn
  gr(i)= xr(i)*grn
  gs(i)= xs(i)*gsn
  gd(i)= xd(i)*gdn
enddo


  gs(9)=0.d0
  gw(9)=0.d0
  gr(9)=0.d0
  gd(9)=0.d0


  gs(10)=0.d0
  gw(10)=0.d0
  gr(10)=0.d0
  gd(10)=0.d0


!Momentum Magnetico Anomalo (TESTE MMA - /1000.d0)
  kb(1)=1.79d0     ! proton 
  kb(2)=-1.91d0    ! neutron 
  kb(3)=-0.61d0    ! lambda
  kb(4)=1.67d0     ! sigma(+)
  kb(5)=1.61d0     ! sigma(0)
  kb(6)=-0.38d0    ! sigma(-)
  kb(7)=-1.25d0    ! cascade (0)
  kb(8)=0.06d0     ! cascade (-)
  kb(9)= 0.00116d0 ! eletron
  kb(10)=0.d0      ! muon 

!Teste de AMM
!do i=1,10
!kb(i)=kb(i)/1000.d0
!enddo

!********************
!* valores iniciais *
!********************

!potencial quimico neutron
mu_n = 940.d0/197.327d0
!mu_n = 1100.d0/197.327d0

x(1) = 45.d0/197.327d0    !potencial quimico eletron
!x(1) = 37.d0/197.327d0
x(2) = asin(sqrt(as*0.05d0/gsn)) 
x(3) = sqrt(aw*0.05d0/gwn) 
x(4) = -ar*0.05d0/(2*grn)   ! Rho03
x(5) = -ad*0.05d0/(2*gdn)   ! Delta

!*********************************** 
!*loop aumentando a densidade total*
!***********************************
open(71, file='eosz_swrd_B0_G1L93asym30_Ulm28_Usp30_Ucp40.dat')
!open(72, file='pop_zA2_B18_amm1.dat')
!open(73, file='g_zA2_B18_amm1.dat')
!open(74, file='anisotropiaB518AMM1.dat')

write(71,*) '#    de_mag,       p_perp,      p_par,        de,', &
           '           p,          rhob,     fs,      mu_n,        mu_e,        Bmag,          mag'

write(72,*) '#  rhob,  fs,   rho_p/rho, rho_n/rho, rho_l/rho, rhos+/rho,', &
            ' rhos0/rho, rhos-/rho, rhox0/rho, rhox-/rho, rho_e/rho, rho_mu/rho'

!write(73,*) '#     rhob,        g*p,         g*n,         m*p,         m*n'

!write(74,*) '#     rhob,        p_perp,         p_par,         mag,       p_perp/p_par,    Bmag'


!***********************************
!*   LOOP PARA POTENCIAL QUIMICO   * 
!***********************************

do while(mu_n<2000.d0/197.327d0)


!**********************
!*  CAMPO MAGNÉTICO   * 
!**********************

if(iB2.eq.1) then
!Bmag = 1.d4/197.327d0**2     !escala 10**6 (MeV**2)=1.444 X 10**19 G (B da superficie)
Bmag = 1.d-6/197.327d0**2    !B=0
!Bmag = 50.d0/197.327d0**2    !10^14 G
!Bmag = 1.d5/197.327d0**2     !10^18 G

endif


if(iB2.eq.2) then       !B_variavel - só usar para mu_n > 939
Bsurf = 69.252077562d0/197.327d0**2  
B_a = 2.5d0
B_b = 4.34708d-7
!Bc = 6.925d3/197.327d0**2  !B central (=1x10^17 G)
!Bc = 3.462d4/197.327d0**2  !B central (=5x10^17 G)
!Bc = 6.925d4/197.327d0**2  !B central (=1x10^18 G)
Bc = 3.462d5/197.327d0**2  !B central (=5x10^18 G)
!Bc = 6.925d5/197.327d0**2  !B central (=10^19 G)

Bmag = Bsurf + Bc*(1.d0 - exp(-B_b*(mu_n*197.327d0 - 938.d0)**B_a))
endif


!CONVERGENCIA

!write(*,*) 'ANTES DO BROYDN'!  

 call broydn(x,check)

 if (check) then
   counter=1
   do while (counter.lt.257)
     mu_n = mu_n + 1.1d0/197.327d0
     call broydn(x,check)
     if (check) then
       counter = counter + 1
     else
       counter = 258
     endif
   enddo
 endif

 if (check) STOP "no root in broydn"

!!write(*,*) 'voltou da brodyn'


!Calculo da EOS  
de = ((ms*sig)**2)/2.d0 +((mw*ome)**2)/2.d0 +((mr*r03)**2)/2.d0 +((md*del3)**2)/2.d0 

p = 0.d0 

pmeson = ((mw*ome)**2)/2.d0 +((mr*r03)**2)/2.d0 -((ms*sig)**2)/2.d0 -((md*del3)**2)/2.d0  !contribuicao dos mesons na pressao!
pperp2 = ((mw*ome)**2)/2.d0 +((mr*r03)**2)/2.d0 -((ms*sig)**2)/2.d0 -((md*del3)**2)/2.d0
p_perp = 0.d0
!ppar = 0.d0
p_charged = 0.d0

mag = 0.d0
mag2 = 0.d0

do i=1,10 !soma sobre barions
  do y = spini(i),spinf(i),spins(i) !soma sobre spin

if (kf(i,y).gt.1.d-6) then 
  
 if((dabs(q(i)).ge.1.d-4).and.(iB(i,y).ne.0)) then
!write(*,*) 'particulas CARREGADAS'
!write(*,*) i, y, kb(i), q(i)
      nnu=0

      do while ((nnu+1./2.-y/2.*q(i)/abs(q(i))).lt.nu_max(i,y))

      nu = nnu + 1./2.-y/2.*q(i)/abs(q(i))+0.1d0
      ef(i) = sqrt(mmag(i,nu,y)**2+ kfmag(i,nu,y)**2)        

	de = de + (gama(i)/(4.d0*pi**2))*abs(q(i))*e*Bmag*(kfmag(i,nu,y)*ef(i) &
+ mmag(i,nu,y)**2*dlog(dabs((kfmag(i,nu,y)+ef(i))/mmag(i,nu,y)))) 


!Novo - p_charged - guarda valor da pressao apenas para particulas carregadas

         p_charged = p_charged + (gama(i)/(4.d0*pi**2))*abs(q(i))*e*Bmag*(kfmag(i,nu,y)*ef(i) - &
mmag(i,nu,y)**2*dlog(dabs((kfmag(i,nu,y)+ef(i))/mmag(i,nu,y)))) 

         p = p + (gama(i)/(4.d0*pi**2))*abs(q(i))*e*Bmag*(kfmag(i,nu,y)*ef(i) - &
mmag(i,nu,y)**2*dlog(dabs((kfmag(i,nu,y)+ef(i))/mmag(i,nu,y)))) 
           
         mag = mag + gama(i)*dabs(e*q(i))*(kfmag(i,nu,y)*ef(i)-mmag2(i,nu,y)**2*(dlog(dabs((kfmag(i,nu,y)+ &
ef(i))/mmag(i,nu,y)))))/(4.d0*pi**2)   

!####################################################################################################################
!Charged particles (veronica)
        pperp2 = pperp2 + gama(i)*dabs(q(i))*e*(Bmag**2)*mmag(i,nu,y)* &
(dabs(q(i))*e*nu/(mmag(i,nu,y)  + y*mun(i)*kb(i)*Bmag) - y*mun(i)*kb(i)) &
*dlog(dabs((kfmag(i,nu,y)+ef(i))/mmag(i,nu,y)))/(2.d0*pi**2)


         mag2 = mag2 + gama(i)*dabs(q(i))*e*Bmag/(2.d0*pi**2)*mmag(i,nu,y)* &
(y*mun(i)*kb(i) - dabs(q(i))*e*nu/(mmag(i,nu,y) + y*mun(i)*kb(i)*Bmag))* &
dlog(dabs((kfmag(i,nu,y)+ef(i))/mmag(i,nu,y)))


!#####################################################################################################################

        nnu=nnu+1

    enddo !fim do loop do nivel de landau
   
  else ! EoS sem B  

!write(*,*) 'particulas NAO carregas'
!write(*,*) i,y, kb(i), q(i)

     ef(i) = sqrt(mspin(i,y)**2+ kf(i,y)**2)        

     de = de + (gama(i)/(2.d0*pi**2))*(.25d0*kf(i,y)*ef(i)**3-y*mun(i)*kb(i)* &
     Bmag/3.*ef(i)**3*(asin(argument(i,y))-pi/2.)-(y* &
     mun(i)*kb(i)*Bmag/6.+mspin(i,y)/8.)*(mspin(i,y) &
     *kf(i,y)*ef(i)+mspin(i,y)**3*dlog(dabs &
     ((kf(i,y)+ef(i))/mspin(i,y)))))
     

     p = p + (gama(i)/(2.d0*pi**2))*(kf(i,y)*(ef(i)**3)/12.-y*mun(i)*kb(i)* &
     Bmag/6.*ef(i)**3*(asin(argument(i,y))-pi/2.)-(y*  &
     mun(i)*kb(i)*Bmag/3.+mspin(i,y)*5./24.)*(mspin(i,y)*kf(i,y)* &
     ef(i))+(y*mun(i)*kb(i)*Bmag/6.+mspin(i,y)/8.)*mspin(i,y)**3* &
     dlog(dabs((kf(i,y)+ef(i))/mspin(i,y))))


!#####################################################################################################################
!Uncharged particles (veronica)

   pperp2 = pperp2 + gama(i)/(48.d0*pi**2)*(ef(i)*kf(i,y)*(-12.d0*(y*mun(i)*kb(i)*Bmag)**2 &
   + 2.d0*ef(i)**2 - 12.d0*mspin(i,y)*y*mun(i)*kb(i)*Bmag - 5.d0*mspin(i,y)**2) + 3.d0*mspin(i,y)**2* &
(mspin(i,y) + 2.d0*y*mun(i)*kb(i)*Bmag)**2 *dlog(dabs((kf(i,y)+ef(i))/mspin(i,y))))    


   mag2 = mag2 + gama(i)/(12.d0*pi**2)*y*mun(i)*kb(i)*(kf(i,y)*ef(i)*(3.d0*y*mun(i)*kb(i)*Bmag + mspin(i,y)) &
   - ef(i)**3*(asin(argument(i,y))-pi/2.) - (3.d0*y*mun(i)*kb(i)*Bmag + 2.d0*mspin(i,y))*mspin(i,y)**2*&
    dlog(dabs((kf(i,y)+ef(i))/mspin(i,y))))

!###################################################################################################################

endif !fim do loop de iB

!write(*,*) i, y, nu_max(i,y), de

else
!write(*,*) i,y, 'k=pequeno'
endif !fim do loop de kpequeno


enddo !fim do loop de spin
enddo !fim loop de particulas

   de_mag = de + (Bmag**2)/(8*pi)

   p_par = p + pmeson - (Bmag**2)/(8*pi)   !! ARRUMAR EXPRESSAO CORRETA???
    p_par0 = p
!   mag2 = mag2 + p/Bmag                      !p_par/Bmag (tirando pmeson do programa) 

!Nova expressao para magnetizacao - leva em conta p/B apenas para particulas carregadas
   mag2 = mag2 + p_charged/Bmag

   p_perp = p + pmeson + (Bmag**2)/(8*pi) - mag2*Bmag  !mag
    
   p_perp0 = p  - mag2*Bmag 
   pperp2 = pperp2 + (Bmag**2)/(8*pi)

!write(*,*) 'pperp=', p_perp, 'pperp2 =', pperp2 

de = de*197.327d0
p = p*197.327d0
pmeson= pmeson*197.327d0

de_mag = de_mag*197.327d0
p_perp = p_perp*197.327d0
p_par = p_par*197.327d0
p_perp0 = p_perp0*197.327d0
p_par0 = p_par0*197.327d0
!write(71,*) de, p+pmeson, rhob, mu_n*197.327d0, mu_e*197.327d0!, Bmag*197.327d0**2 
write(71,'(11F13.4)') de_mag, p_perp, p_par, de, p + pmeson, rhob, fs, mu_n*197.327d0, &
                      mu_e*197.327d0, Bmag*197.327d0**2, mag2*197.327d0**2
write(72,'(12F12.7)') rhob, fs, rho(1,1)/rhob + rho(1,-1)/rhob, & 
rho(2,1)/rhob + rho(2,-1)/rhob, rho(3,1)/rhob + rho(3,-1)/rhob, &
rho(4,1)/rhob + rho(4,-1)/rhob, rho(5,1)/rhob + rho(5,-1)/rhob, &
rho(6,1)/rhob + rho(6,-1)/rhob, rho(7,1)/rhob + rho(7,-1)/rhob, &
rho(8,1)/rhob + rho(8,-1)/rhob, rho(9,1)/rhob + rho(9,-1)/rhob, &
rho(10,1)/rhob + rho(10,-1)/rhob

!write(73,'(5F13.7)') rhob, as*m1(mn, sig, del3, i3(1), gsn, gdn)**2, &
!                           as*m1(mn, sig, del3, i3(2), gsn, gdn)**2, meff(1)*197.327d0, meff(2)*197.327d0

!write(74,*) rhob, p_perp0, p_par0, mag2, p_perp/p_par, Bmag



mu_n = mu_n + 1.d0/197.327d0
!mu_n = mu_n + .5d0/197.327d0

enddo


end program 


!*****************************************************
function funcv(x)
USE nrtype; USE global
implicit none

REAL(SP), DIMENSION(:), INTENT(IN) :: x
REAL(SP), DIMENSION(size(x)) :: funcv
double precision :: m1, qt
double precision :: a1,b1,c1, d1,e1,f1,cte
double precision :: sigma, omega, rho03, delta3
double precision, dimension(10) :: mueff
double precision, DIMENSION (8,-1:1) :: rs
integer :: i

!Mapping
mu_e = x(1) 	       !Potencial quimico do eletron
sig = sin(x(2))**2     !Sigma entre 0 e 1
ome = x(3)**2          !Omega positivo
r03 = x(4)             !Rho03 sem mapping
del3 = x(5)            !Delta3 sem mapping 

do i = 1, 10 !barions e leptons

!Massas efetivas 
  meff(i) = m(i) - gs(i)*sig*m1(m(i),sig,del3,i3(i),gs(i),gd(i)) &
          - gd(i)*i3(i)*del3*m1(m(i),sig,del3,i3(i),gs(i),gd(i))

!Potenciais quimicos 
  mu(i) = qb(i)*mu_n- q(i)*mu_e 

!Potencial quimico efetivo
  mueff(i) = mu(i) -gr(i)*r03*i3(i)-gw(i)*ome

enddo !fim loop particulas

	iB = float(iB2)

rhob = 0.d0
qt = 0.d0
fs = 0.d0
rho= 0.d0
!write(*,*) 'INICIO DO LOOP DE PARTICULAS'
do i= 1, 10 !loop de particulas

	ef(i) = mueff(i)

	if(AMM.eq.0)then  !caso sem MMA
	  kb(i)=0.
	endif
!write(*,*) 'INICIO DO LOOP DE spin'
      do y=spini(i), spinf(i),spins(i) !soma sobre spin: inicial, final, passo        

!Calculo de mspin (massa nao-magnetica com modificaçao de spin)

      mspin(i,y)=meff(i) -y*mun(i)*kb(i)*Bmag     


   argument(i,y)=mspin(i,y)/ef(i)                                       !expressao para rho
      if(argument(i,y).gt.1.)then
      !write(*,*)'????????????????????????????????????????argument>1'
      argument(i,y)=1.
      endif


!Limiar de mspin

  	if((ef(i)**2 - mspin(i,y)**2).gt.0.d0) then     
		kf(i,y) = sqrt(ef(i)**2 - mspin(i,y)**2)    
	else
		kf(i,y) = 0.d0
               	iB(i,y)=0
	endif


!****************************
! *    CASO MAGNETICO       *
!****************************

!Limiar para continuo
!write(*,*) i, y, q(i), iB(i,y)
if((dabs(q(i)).ge.1.d-4).and.(iB(i,y).ne.0)) then  

 	 nu_max(i,y) = ((ef(i)+y*mun(i)*kb(i)*Bmag)**2 - meff(i)**2)/  &
(2.d0*dabs(q(i))*e*Bmag)                                    !(caso com MMA)

!write(*,*) 'numax', nu_max(i,y), i, y, q(i), iB(i,y)

             if(nu_max(i,y).gt.15000.d0) then
                 iB(i,y) = 0                                     !desligar o campo magnetico
             endif
endif


!Calculo do nivel maximo de landau (nao muda com spin) 
 if((dabs(q(i)).ge.1.d-4).and.(iB(i,y).ne.0)) then      !inicio -if landau -     
	
!Calculo de mmag (massa magnetica)
      nnu=0
      do while ((nnu+1./2.-y/2.*q(i)/abs(q(i))).lt.nu_max(i,y))    !soma sobre niveis de landau

 !     !WRITE(*,*) 'ENTROU, nnu = ', nnu, 'nu_max = ', nu_max(i,y), iB(i,y), i, y

        nu = nnu + 1./2.-y/2.*q(i)/abs(q(i))+0.1d0

 	mmag(i,nu,y) = dsqrt(mspin(i,y)**2 + 2.d0*dabs(q(i))*e*Bmag*nu) -y*mun(i)*kb(i)*Bmag

        mmag2(i,nu,y)= dsqrt(mspin(i,y)**2 + 4.d0*nu*dabs(q(i))*Bmag*e) 


              !Limiar de mmag
	       if((ef(i)**2 - mmag(i,nu,y)**2).gt.0.d0) then
		kfmag(i,nu,y) = dsqrt(ef(i)**2 - mmag(i,nu,y)**2)
		  else
		kfmag(i,nu,y) = 0.d0

  		endif

!Densidade
        rho(i,y) = rho(i,y)+ (kfmag(i,nu,y)*abs(q(i))*e*Bmag/(2.d0*pi**2))
        
        nnu=nnu+1

      enddo        !fim da soma - niveis de landau -


          else  ! B=0 e/ou uncharged

!Densidade de particulas (B=0 e considerando efeitos MMA)

 rho(i,y)= rho(i,y) + gama(i)/(2.d0*pi**2)*(kf(i,y)**3/(3.d0) & 
     -y*mun(i)*kb(i)*Bmag/(2.d0)*(mspin(i,y)*kf(i,y)+ &
     ef(i)**2*(asin(argument(i,y))-pi/2.d0))) 

  endif !fim if iB


!Densidade barionica

  rhob= rhob + qb(i)*rho(i,y)
!WRITE(*,*) i, y, rhob

!Conservacao de carga
  qt= qt+ q(i)*rho(i,y)

!Strangeness
  fs= fs- qs(i)*rho(i,y)
   enddo! fim loop spin
!write(*,*) 'SAIU DO LOOP DE SPIN'
enddo ! fim loop particulas

!write(*,*) 'SAIU DO LOOP DE PARTICULAS'


fs = fs/rhob

!***************************************
!*   Calculo dos Campos Mesonicos      *
!***************************************

 a1 = 0.d0   ! Constantes para escrever sigma = a1 - b1*sigma - c1*delta
 b1 = 0.d0
 c1 = 0.d0

 d1 = 0.d0
 e1 = 0.d0   ! Constantes para escrever delta = d1 - e1*sigma - f1*delta
 f1 = 0.d0

 omega = 0.d0
 rho03 = 0.d0

!Campo Sigma
!write(*,*) 'ENTROU NO LOOP DE PARTICULAS'
do i=1,8 !loop de particulas
!write(*,*) 'ENTROU NO LOOP DE SPIN'
   do y=spini(i), spinf(i),spins(i) !soma sobre spin: inicial, final, passo        

    if((dabs(q(i)).ge.1.d-4).and.(iB(i,y).ne.0)) then      !inicio -if landau -     
!write(*,*) 'forbiden'
    nnu=0
    rs(i,y)=0.d0
        	do while ((nnu+1./2.-y/2.*q(i)/abs(q(i))).lt.nu_max(i,y))    !soma sobre niveis de landau
                !Densidade Escalar (com B) 
       !          write(*,*) 'ENTROU NO LOOP LANDAU' 
                nu = nnu + 1./2.-y/2.*q(i)/abs(q(i))+0.1d0                 

                ef(i) = sqrt(mmag(i,nu,y)**2+ kfmag(i,nu,y)**2)        

               rs(i,y)= rs(i,y) + (dabs(q(i))*e*Bmag*meff(i)*(mmag(i,nu,y)/(sqrt(meff(i)**2 + 2.d0*nu*dabs(q(i))*e*Bmag)))&
*log(dabs((kfmag(i,nu,y)+ef(i))/mmag(i,nu,y))))/(2.d0*pi**2)

		nnu=nnu+1

	        enddo        !fim da soma - niveis de landau -

     else
!Densidade Escalar (sem B e/ou uncharged) 
!write(*,*) 'ENTROU ELSE ############################'
     ef(i) = sqrt(mspin(i,y)**2+ kf(i,y)**2)        

     rs(i,y)=meff(i)*(kf(i,y)*ef(i) - (mspin(i,y)**2)*log(abs((kf(i,y)+ef(i))/mspin(i,y))))/(4.d0*pi**2)

     endif !fim do if de iB

 cte = (rs(i,y)/m(i))*m1(m(i),sig,del3,i3(i),gs(i),gd(i))**((1.d0 + lambda)/lambda)
 
 a1 = a1 + gs(i)*rs(i,y)*m1(m(i),sig,del3,i3(i),gs(i),gd(i))/(ms**2)
 b1 = b1 + cte*gsn*gs(i)/(ms**2)
 c1 = c1 + cte*gsn*gd(i)*i3(i)/(ms**2)

 d1 = d1 + gd(i)*rs(i,y)*i3(i)*m1(m(i),sig,del3,i3(i),gs(i),gd(i))/(md**2)
 e1 = e1 + cte*gdn*i3(i)*gs(i)/(md**2)
 f1 = f1 + cte*gdn*gd(i)*(i3(i)**2)/(md**2)

!Campo Omega e Rho_03 (formulas x2(2.371) e x3(2.372) generalizadas para outros acoplamentos)

  omega = omega + (gw(i)/(mw**2))*rho(i,y)       
  rho03 = rho03 + (gr(i)/(mr**2))*rho(i,y)*i3(i) 
!write(*,*) 'DENTRO DO LOOP DE spin', kf(i,y), i, y
    enddo !fim loop spin

enddo !fim loop particulas
!write(*,*) 'SAIU DO LOOP DE particulas'

!Campos Sigma e Delta (solucao do sistema de equacoes de campo escritas em termos das constantes)
sigma = (a1 - c1*d1/(1.d0+f1))/(1.d0 + b1 - c1*e1/(1.d0+f1))

delta3 = (d1 - e1*a1/(1.d0+b1))/(1.d0 + f1 - c1*e1/(1.d0+b1))

!Funcões a serem zeradas
funcv(1) = qt
funcv(2) = sig - sigma
funcv(3) = ome - omega
funcv(4) = r03 - rho03
funcv(5) = del3 - delta3

!write(*,*) funcv(1), funcv(2), funcv(3), funcv(4), funcv(5)
end function funcv

!********************************************

double precision FUNCTION m1(m, sig, del, i3, gs, gd)  !ARRUMAR - FUNCAO COM gsb e gdb
USE global, only: lambda
implicit none
!double precision :: m1
double precision, INTENT( IN ) :: m, sig, del, i3, gs, gd

m1 = ((lambda*m)/(lambda*m + sig*gs + del*gd*i3))**lambda

END FUNCTION m1

!********************************************************

include 'brodyn.f90'

