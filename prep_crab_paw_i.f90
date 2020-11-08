module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80

  integer,parameter :: Max_orbit=15
! constant
  real(8),parameter :: pi=3.1415926535897932d0
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
  real(8) gaunt(16,16,16),nYlm(0:3,-3:3)
! table
  integer,allocatable :: nR_phi(:),nL_phi(:)
    

! log-mash parameter
  integer Nx
  real(8) Dx,Rmax,Rp
  real(8),allocatable :: xL(:),rL(:),expXL(:)

! DFT parameter 
  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
  real(8),parameter :: CU=0.002d0,DU=-0.0116d0    

! GS
  real(8),allocatable :: upsi(:,:),Vh(:),Eexc(:),Vexc(:),Veff(:)
  real(8),allocatable :: VL(:),Vnucl(:)
  real(8),allocatable :: rho(:)
  real(8),allocatable :: occA(:),esp(:)
  integer Nupsi
  real(8),allocatable :: L_upsi(:)
  real(8),parameter :: rate=0.1d0
! atom date
  integer ZA
  character(5) atom_name
  real(8),allocatable :: occ_input(:)
  integer,allocatable :: orbit_L(:)
  character(10),allocatable :: orbit_state(:)
  
! PAW
  real(8) grid_size_3D
  integer NuPR_3D,NuPR
  integer Nunocc(0:3)
  integer,allocatable :: L_proj(:)
  real(8),allocatable :: Cr_phi_P(:,:)
  real(8),allocatable :: Rcut(:)
  integer,allocatable :: iRc(:)
  real(8) Rfilt,Rcomp,Rcore,RcutPR,softning
  integer iRcf,iRcomp,iRcore,iRcPR
  real(8),allocatable :: uAE(:,:),uPS(:,:),uPR(:,:)
  real(8),allocatable :: nc_A(:),nc_P(:),V_bar(:),gL(:,:)
! PAW constant
  real(8) A_c,Delta_a,F_c,core_kinetic_energy
  real(8),allocatable :: B_c(:,:),I_c(:,:),Delta_c(:,:,:),K_c(:)
  real(8),allocatable :: M_c(:,:,:),N_c(:,:)
  real(8),allocatable :: J_c(:,:,:,:),C_c(:,:,:,:),X_c(:,:)
  real(8),allocatable :: dip_c(:,:,:)
    


end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80
program main
  use global_variables
  implicit none
  integer iter,p,ix
  integer i,j

! calc  Ground state ------------------------------------------------!
  call prep_calc_parameter
  call init_wf
  call wf_normalize
  call upsi_rho
  call Hartree(0)
  call Exc_Cor


  do iter=1,200
    call solve_Kohn_Sham
    call wf_normalize
    call upsi_rho
    call Hartree(0)
    call Exc_Cor       
  end do

  do p=1,Nupsi
    write(*,*)p,esp(p)*2d0*Ry
  end do

  write(*,*)'End prepare'
  write(*,*)'Rmax=',Rmax
  write(*,*)'Dx=',Dx,1d0/32d0
  write(*,*)'Nx=',Nx
  write(*,*)'Rp=',Rp
    
! end calc  Ground state --------------------------------------------!
! preparation  PAW --------------------------------------------------!    

  write(*,*)'radial function'
  call make_radial_function

  write(*,*)'PAW const'
  call make_PAW_constant
  call output_constant

  call filtering_rojector
  call output_radial_function

! end preparation  PAW ----------------------------------------------!    

  rho=rL*exp(-rL**2)
  call hartree(1)
  open(20,file='vh.dat')
  do i=0,Nx
    write(20,'(100e20.10)')rL(i),Vh(i)
  end do
  close(20)
  
end program main
!-------10--------20--------30--------40--------50--------60----------72
subroutine prep_calc_parameter
  use global_variables
  implicit none
  integer i,j,ix,m
  
  call calc_gaunt_factor
  call prep_orbit_L
  
  allocate(occ_input(Max_orbit),orbit_state(Max_orbit))
    
  
! input file
  read(*,*)atom_name   ! atom name
  read(*,*)ZA          ! nuclear chage
  read(*,*)Nunocc(:)   !number unoccupied orbit
  do i=1,Max_orbit
    read(*,*)occ_input(i),orbit_state(i) ! occupation & state
  end do
  read(*,*)grid_size_3D! grid size used 3D calculation
  read(*,*)softning    ! soft compensation 
  read(*,*)Dx           ! grid size
  read(*,*)Nx          ! mesh points
  read(*,*)Rmax        ! radious simulation sphere

! count the number
  NuPR_3D=0
  NuPR=0
  Nupsi=0
  do i=1,Max_orbit
    if((orbit_state(i)=='v')) then
      NuPR=NuPR+1
      NuPR_3D=NuPR_3D+orbit_L(i)*2+1
      Nupsi=Nupsi+1
    end if
    if(orbit_state(i)=='c')Nupsi=Nupsi+1
  end do
  NuPR_3D=NuPR_3D &
    +Nunocc(0)*1 &
    +Nunocc(1)*(2*1+1) &
    +Nunocc(2)*(2*2+1) &
    +Nunocc(3)*(2*3+1)
  NuPR=NuPR+sum(Nunocc)
  Nupsi=Nupsi+sum(Nunocc)
  
  allocate(Rcut(NuPR),iRc(NuPR))
  read(*,*)Rcomp
  read(*,*)Rfilt
  read(*,*)Rcore
  do i=1,NuPR
    read(*,*)Rcut(i)
  end do

    
      
  Rp=Rmax/(exp(Dx*dble(Nx))-1d0)
  iRc=aint(log(Rcut/Rp+1d0)/Dx)
  iRcomp=aint(log(Rcomp/Rp+1d0)/Dx)
  iRcf=aint(log(Rfilt/Rp+1d0)/Dx)
  iRcore=aint(log(Rcore/Rp+1d0)/Dx)
  RcutPR=maxval(Rcut)
  iRcPR=aint(log(RcutPR/Rp+1d0)/Dx)
  
  
  allocate(xL(0:Nx),rL(0:Nx),expXL(0:Nx))
  do ix=0,Nx
    xL(ix)=Dx*dble(ix)
    rL(ix)=Rp*(exp(Dx*dble(ix))-1d0)
    expXL(ix)=exp(Dx*dble(ix))
  end do
  
    
  allocate(L_proj(NuPR),L_upsi(Nupsi))
  allocate(occA(Nupsi),esp(Nupsi))
  allocate(nR_phi(NuPR_3D),nL_phi(NuPR_3D))
  
  NuPR=0
  Nupsi=0
  do i=1,Max_orbit
    if((orbit_state(i)=='v'))then
      NuPR=NuPR+1
      Nupsi=Nupsi+1
      L_proj(NuPR)=orbit_L(i)
      L_upsi(Nupsi)=orbit_L(i)
      occA(Nupsi)=occ_input(i)
    end if
    if(orbit_state(i)=='c')then
      Nupsi=Nupsi+1
      L_upsi(Nupsi)=orbit_L(i)
      occA(Nupsi)=occ_input(i)
    end if
  end do
  
  do i=0,3
    do j=1,Nunocc(i)
      NuPR=NuPR+1
      Nupsi=Nupsi+1
      L_proj(NuPR)=i
      L_upsi(Nupsi)=i
      occA(Nupsi)=0d0
    end do
  end do

! make table
  j=0
  do i=1,NuPR
    do m=L_proj(i),-L_proj(i),-1
      j=j+1
      nR_phi(j)=i
      nL_phi(j)=nYlm(L_proj(i),m)
    end do
  end do
  

  allocate(upsi(0:Nx,Nupsi),rho(0:Nx))
  rho=0d0
  allocate(Vh(0:Nx),Vexc(0:Nx),Veff(0:Nx),Eexc(0:Nx))
  allocate(VL(0:Nx),Vnucl(0:Nx))
  
  Vnucl(1:Nx)=-dble(ZA)/rL(1:Nx);Vnucl(0)=0d0
  VL(1:Nx)=0.5d0/(rL(1:Nx)**2);VL(0)=0d0
  
  return
end subroutine prep_calc_parameter
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine prep_orbit_L
  use global_variables
  implicit none
  
  allocate(orbit_L(Max_orbit))

  orbit_L(1)=0   !1s
  orbit_L(2)=0   !2s
  orbit_L(3)=1   !2p
  orbit_L(4)=0   !3s
  orbit_L(5)=1   !3p
  orbit_L(6)=2   !3d
  orbit_L(7)=0   !4s
  orbit_L(8)=1   !4p
  orbit_L(9)=2   !4d
  orbit_L(10)=3  !4f
  orbit_L(11)=0  !5s
  orbit_L(12)=1  !5p
  orbit_L(13)=2  !5d
  orbit_L(14)=0  !6s
  orbit_L(15)=1  !6p
  return
end subroutine prep_orbit_L
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_gaunt_factor
  use global_variables
  implicit none
  integer,parameter :: Ntheta=15000+1,Nphi=15000+1
  integer itheta,iphi
  integer l1,l2,l3
  real(8) dtheta,dphi
  real(8) fco(16),fth(Ntheta,16),fph(Nphi,16)
  real(8) theta(Ntheta),phi(Nphi)
  real(8) Sco,Sth,Sph,Stot
  real(8) t1,t2,t3
! theta(i)=dtheta*(i-1) i=1,Ntheta
! phi(i)=dphi*(i-1) i=1,Nphi

  call cpu_time(t1)
  dtheta=pi/dble(Ntheta-1)
  dphi=2d0*pi/dble(Nphi-1)
  do itheta=1,Ntheta
    theta(itheta)=dtheta*dble(itheta-1)
  end do
  do iphi=1,Nphi
    phi(iphi)=dphi*dble(iphi-1)
  end do

!1 l=0 m=0
  nYlm(0,0)=1
  fco(1)=1d0/sqrt(4d0*pi)
  fth(:,1)=1d0
  fph(:,1)=1d0
!2 l=1 m=+1
  nYlm(1,1)=2
  fco(2)=sqrt(3d0/4d0/pi)
  fth(:,2)=sin(theta(:))
  fph(:,2)=cos(phi(:))
!3 l=1 m=0
  nYlm(1,0)=3
  fco(3)=sqrt(3d0/4d0/pi)
  fth(:,3)=cos(theta(:))
  fph(:,3)=1d0
!4 l=1 m=-1
  nYlm(1,-1)=4
  fco(4)=sqrt(3d0/4d0/pi)
  fth(:,4)=sin(theta(:))
  fph(:,4)=sin(phi(:))
!5 l=2 m=+2
  nYlm(2,2)=5
  fco(5)=sqrt(15d0/pi)/8d0
  fth(:,5)=1d0-cos(2d0*theta(:))
  fph(:,5)=cos(2d0*phi(:))
!6 l=2 m=+1
  nYlm(2,1)=6
  fco(6)=sqrt(15d0/pi)/4d0
  fth(:,6)=sin(2d0*theta(:))
  fph(:,6)=cos(phi(:))
!7 l=2 m=0
  nYlm(2,0)=7
  fco(7)=sqrt(5d0/pi)/8d0
  fth(:,7)=1d0+3d0*cos(2d0*theta(:))
  fph(:,7)=1d0
!8 l=2 m=-1
  nYlm(2,-1)=8
  fco(8)=sqrt(15d0/pi)/4d0
  fth(:,8)=sin(2d0*theta(:))
  fph(:,8)=sin(phi(:))
!9 l=2 m=-2
  nYlm(2,-2)=9
  fco(9)=sqrt(15d0/pi)/8d0
  fth(:,9)=1d0-cos(2d0*theta(:))
  fph(:,9)=sin(2d0*phi(:))
!10 l=3 m=+3
  nYlm(3,3)=10
  fco(10)=sqrt(70d0/pi)/32d0
  fth(:,10)=3d0*sin(theta(:))-sin(3d0*theta(:))
  fph(:,10)=cos(3d0*phi(:))
!11 l=3 m=+2
  nYlm(3,2)=11
  fco(11)=sqrt(105d0/pi)/16d0
  fth(:,11)=cos(theta(:))-cos(3d0*theta(:))
  fph(:,11)=cos(2d0*phi(:))
!12 l=3 m=+1
  nYlm(3,1)=12
  fco(12)=sqrt(42d0/pi)/32d0
  fth(:,12)=sin(theta(:))+5d0*sin(3d0*theta(:))
  fph(:,12)=cos(phi(:))
!13 l=3 m=0
  nYlm(3,0)=13
  fco(13)=sqrt(7d0/pi)/16d0
  fth(:,13)=3d0*cos(theta(:))+5d0*cos(3d0*theta(:))
  fph(:,13)=1d0
!14 l=3 m=-1
  nYlm(3,-1)=14
  fco(14)=sqrt(42d0/pi)/32d0
  fth(:,14)=sin(theta(:))+5d0*sin(3d0*theta(:))
  fph(:,14)=sin(phi(:))
!15 l=3 m=-2
  nYlm(3,-2)=15
  fco(15)=sqrt(105d0/pi)/16d0
  fth(:,15)=cos(theta(:))-cos(3d0*theta(:))
  fph(:,15)=sin(2d0*phi(:))
!16 l=3 m=-3
  nYlm(3,-3)=16
  fco(16)=sqrt(70d0/pi)/32d0
  fth(:,16)=3d0*sin(theta(:))-sin(3d0*theta(:))
  fph(:,16)=sin(3d0*phi(:))

  do l1=1,16
    do l2=1,16
      do l3=1,16
        Sco=fco(l1)*fco(l2)*fco(l3)
        Sth=sum(fth(2:Ntheta-1,l1)*fth(2:Ntheta-1,l2) &
          *fth(2:Ntheta-1,l3)*sin(theta(2:Ntheta-1)))
        Sph=sum(fph(2:Nphi-1,l1)*fph(2:Nphi-1,l2)*fph(2:Nphi-1,l3)) &
          +0.5d0*(fph(1,l1)*fph(1,l2)*fph(1,l3) &
          +fph(Nphi,l1)*fph(Nphi,l2)*fph(Nphi,l3))
        Stot=Sco*Sth*Sph*dtheta*dphi
        gaunt(l1,l2,l3)=Stot
        gaunt(l1,l3,l2)=Stot
        gaunt(l2,l1,l3)=Stot
        gaunt(l2,l3,l1)=Stot
        gaunt(l3,l2,l1)=Stot
        gaunt(l3,l1,l2)=Stot
      end do
    end do
  end do
  return
end subroutine calc_gaunt_factor
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine init_wf
  use global_variables
  implicit none
  real(8) r
  integer ix,p
  
  do p=1,Nupsi
    upsi(:,p)=rL(:)*exp(-5.125*rL(:)**2)
  end do
  
  return
end subroutine init_wf
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine wf_normalize
  use global_variables
  implicit none
  integer p
  real(8) S
  
  do p=1,Nupsi
    S=sum(upsi(:,p)**2*expXL(:))*Dx*Rp
    S=1d0/sqrt(S)
    upsi(:,p)=S*upsi(:,p)
  end do
  
  return
end subroutine wf_normalize
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine upsi_rho
  use global_variables
  implicit none
  real(8),allocatable :: rho_old(:)
  integer p
  
  allocate(rho_old(0:Nx))
  rho_old=rho
  rho=0d0
  do p=1,Nupsi
    rho(1:Nx)=rho(1:Nx)+occA(p)/(4d0*pi)*(upsi(1:Nx,p)/rL(1:Nx))**2
  end do
  rho(0)=2d0*rho(1)-rho(2)
  rho=rho*rate+(1d0-rate)*rho_old
  return
end subroutine upsi_rho
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine Hartree(l)
  use global_variables
  implicit none
  integer ix,l
  real(8),allocatable :: Sn(:),Vn(:),Kn(:)
  real(8) Qrho,Const
  real(8) fact1,fact2,fact3
  real(8) temp1,temp2
  
  allocate(Sn(0:Nx),Vn(0:Nx),Kn(0:Nx))
  
  Qrho=sum(expXL*rL**(2+l)*rho)*Rp*Dx*4d0*pi/dble(2*l+1)
  Sn=-4d0*pi*Rp**2*rL*exp(1.5d0*xL)*rho
  fact1=Dx**2/12d0
  fact2=5d0*Dx**2/12d0
  Vn(0:2)=rL(0:2)**(l+1)*exp(-0.5d0*xL(0:2))
  do ix=1,Nx
    Kn(ix)=-0.25d0-(Rp*expXL(ix)/rL(ix))**2*dble(l*(l+1))
  end do
  
  do ix=3,Nx
    temp1=2d0*(1d0-fact2*Kn(ix-1))*Vn(ix-1) &
      -(1d0+fact1*Kn(ix-2))*Vn(ix-2)
    temp2=fact1*(Sn(ix)+10d0*Sn(ix-1)+Sn(ix-2))
    Vn(ix)=(temp1+temp2)/(1d0+fact1*Kn(ix))
  end do
  
  Vh(1:Nx)=Vn(1:Nx)*exp(0.5d0*xL(1:Nx))/rL(1:Nx)
  Const=Qrho/(rL(Nx)**(2*l+1))-Vh(Nx)/(rL(Nx)**l)
  Vh(1:Nx)=Vh(1:Nx)+Const*rL(1:Nx)**l
  Vh(0)=2d0*Vh(1)-Vh(2)
  return
end subroutine Hartree
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine Exc_Cor
  use global_variables
  implicit none
  integer ix
  real(8) trho,rs,rssq,rsln,V_xc,E_xc
  
  do ix=0,Nx
    trho=rho(ix)+1d-20
    rs=(3d0/(4*Pi*trho))**(1d0/3d0)
    V_xc=-4d0/3d0*0.4582d0/rs
    E_xc=-.4582d0/rs
    if(rs>1d0) then
      rssq=sqrt(rs)
      V_xc=V_xc+gammaU*(1d0+7d0/6d0*beta1U*rssq+4d0/3d0*beta2U*rs) &
        /(1d0+beta1U*rssq+beta2U*rs)**2
      E_xc=E_xc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
    else
      rsln=log(rs)
      V_xc=V_xc+AU*rsln+(BU-AU/3d0) &
        +2d0/3d0*CU*rs*rsln+(2d0*DU-CU)/3d0*rs
      E_xc=E_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
    endif
    Vexc(ix)=V_xc
    Eexc(ix)=E_xc
  end do
  return
end subroutine Exc_Cor
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine solve_Kohn_Sham
  use global_variables
  implicit none
  integer Node_L(0:3),node,p,ix,n
  real(8) epsilon,ref_epsilon,epsilon_min,epsilon_max
  real(8) fact1,temp1,temp2
  real(8),allocatable :: v_upsi(:,:),K2(:)
  allocate(v_upsi(0:Nx,Nupsi),K2(0:Nx))

  Node_L=0
  ref_epsilon=1d-12
  fact1=Dx**2/12d0
  
  do p=1,Nupsi
    v_upsi(0:2,p)=rL(0:2)**(L_upsi(p)+1)*exp(-0.5d0*XL(0:2))
    Veff=Vh+Vexc+Vnucl+VL*dble(L_upsi(p)*(L_upsi(p)+1))
    epsilon_min=-0.5d0*dble(ZA)**2-10d0
    epsilon_max=3d0       
! search eigen value -------------------------------------------------!
    do
      epsilon=0.5d0*(epsilon_max+epsilon_min)
      K2=-0.25d0+2d0*((Rp*expXL)**2)*(epsilon-Veff)
      node=0
      do ix=3,Nx
        temp1=2d0*(1d0-5d0*fact1*K2(ix-1))*v_upsi(ix-1,p) &
          -(1d0+fact1*K2(ix-2))*v_upsi(ix-2,p)
        temp2=(1d0+fact1*K2(ix))
        v_upsi(ix,p)=temp1/temp2
        if(v_upsi(ix-1,p)*v_upsi(ix,p)<0d0)node=node+1
        if(abs(v_upsi(ix,p))>=1d2)exit
      end do
      
      if(node>node_L(L_upsi(p)))then
        epsilon_max=epsilon
      else
        epsilon_min=epsilon
      end if
      
      if(epsilon_max-epsilon_min<=ref_epsilon)exit
    end do
! end search eigen value ---------------------------------------------!
! calc wave function -------------------------------------------------!
    epsilon=epsilon_max
    esp(p)=epsilon
    K2=-0.25+2d0*(Rp*expXL)**2*(epsilon-Veff)
    node=0
    do ix=3,Nx
      temp1=2d0*(1d0-5d0*fact1*K2(ix-1))*v_upsi(ix-1,p) &
        -(1d0+fact1*K2(ix-2))*v_upsi(ix-2,p)
      temp2=(1d0+fact1*K2(ix))
      v_upsi(ix,p)=temp1/temp2
      if(v_upsi(ix,p)*v_upsi(ix-1,p)<0d0)node=node+1
      if(node>node_L(L_upsi(p)))exit
    end do
    v_upsi(ix:Nx,p)=0d0
! end calc wave function ---------------------------------------------!
    node_L(L_upsi(p))=node_L(L_upsi(p))+1
  end do
  
  do p=1,Nupsi
    upsi(:,p)=v_upsi(:,p)*exp(0.5d0*XL(:))
  end do
  return
end subroutine solve_Kohn_Sham
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_radial_function
  use global_variables
  implicit none
  
  allocate(uAE(0:Nx,NuPR))
  allocate(uPS(0:Nx,NuPR))
  allocate(uPR(0:Nx,NuPR))
  allocate(nc_A(0:Nx),nc_P(0:Nx),V_bar(0:Nx))
  allocate(gL(0:Nx,0:3))
  allocate(Cr_phi_P(4,NuPR))
  
  call make_gL
  call make_uAE
  call make_uPS
  call make_core_density
  call make_projector
  call make_V_bar
  call orthonormalize_phi_proj
    
  return
end subroutine make_radial_function
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_gL
  use global_variables
  implicit none
  integer l,l2
  real(8) alpha,fact1,fact2
  
  alpha=9d0/(Rcomp**2)
! l=0
  gL(:,0)=(4d0*alpha)**1.5d0/sqrt(4d0*pi)*exp(-alpha*rL(:)**2)
  do l=1,3
    l2=2*l+1
    call factorial(l,fact1)
    call factorial(l2,fact2)
    gL(:,l)=fact1/fact2*(4d0*alpha)**(dble(l)+1.5d0)/sqrt(4d0*pi) &
      *rL(:)**l*exp(-alpha*rL(:)**2)
  end do
  
  return
end subroutine make_gL
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_uAE
  use global_variables
  implicit none
  integer npsi,nphi,n
  nphi=0
  npsi=0
  do n=1,Max_orbit
    if((orbit_state(n)=='v'))then
      npsi=npsi+1
      nphi=nphi+1
      uAE(:,nphi)=upsi(:,npsi)
    else if(orbit_state(n)=='c')then
      npsi=npsi+1
    end if
  end do
  uAE(:,nphi+1:nphi+sum(Nunocc))=upsi(:,npsi+1:npsi+sum(Nunocc))
  return
end subroutine make_uAE
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_uPS
  use global_variables
  implicit none
  integer p,Nmax,i,j,l
  real(8) Amat(4,4),Bmat(4,4),bvec(4),cvec(4)
  real(8) rate1,rate2,temp1,temp2

  Nmax=4
  
  uPS=uAE
  do p=1,NuPR
    rate1=(Rcut(p)-rL(iRc(p)))/(rL(iRc(p)+1)-rL(iRc(p)))
    rate2=(rL(iRc(p)+1)-Rcut(p))/(rL(iRc(p)+1)-rL(iRc(p)))
    l=L_proj(p)
! make bvec
!! wave function
    temp1=uAE(iRc(p),p);temp2=uAE(iRc(p)+1,p)
    bvec(1)=temp1*rate2+temp2*rate1
!! first order
    temp1=0.5d0*exp(-xL(iRc(p)))/Rp*(uAE(iRc(p)+1,p)-uAE(iRc(p)-1,p))/Dx
    temp2=0.5d0*exp(-xL(iRc(p)+1)) &
      /Rp*(uAE(iRc(p)+2,p)-uAE(iRc(p),p))/Dx
    bvec(2)=temp1*rate2+temp2*rate1
!! second order
    temp1=exp(-2d0*xL(iRc(p)))/(Rp**2)*( &
      (uAE(iRc(p)+1,p)-2d0*uAE(iRc(p),p)+uAE(iRc(p)-1,p))/Dx**2 &
      -0.5d0*(uAE(iRc(p)+1,p)-uAE(iRc(p)-1,p))/Dx) 
    temp2=exp(-2d0*xL(iRc(p)+1))/(Rp**2)*( &
      (uAE(iRc(p)+2,p)-2d0*uAE(iRc(p)+1,p)+uAE(iRc(p),p))/Dx**2 &
      -0.5d0*(uAE(iRc(p)+2,p)-uAE(iRc(p),p))/Dx) 
    bvec(3)=temp1*rate2+temp2*rate1
!! third order
    temp1=exp(-3d0*xL(iRc(p)))/Rp**3*( &
      0.5d0*(uAE(iRc(p)+2,p)-2d0*uAE(iRc(p)+1,p) &
      +2d0*uAE(iRc(p)-1,p)-uAE(iRc(p)-2,p))/Dx**3 &
      -3d0*( &
      uAE(iRc(p)+1,p)-2d0*uAE(iRc(p),p)+uAE(iRc(p)-1,p))/Dx**2 &
      +(uAE(iRc(p)+1,p)-uAE(iRc(p)-1,p))/Dx)
    temp2=exp(-3d0*xL(iRc(p)+1))/Rp**3*( &
      0.5d0*(uAE(iRc(p)+3,p)-2d0*uAE(iRc(p)+2,p) &
      +2d0*uAE(iRc(p),p)-uAE(iRc(p)-1,p))/Dx**3 &
      -3d0*( &
      uAE(iRc(p)+2,p)-2d0*uAE(iRc(p)+1,p)+uAE(iRc(p),p))/Dx**2 &
      +(uAE(iRc(p)+2,p)-uAE(iRc(p),p))/Dx)
    bvec(4)=temp1*rate2+temp2*rate1
! make Amat
    do i=1,4
      Amat(1,i)=Rcut(p)**(2*i-1+l)
      Amat(2,i)=dble(2*i-1+l)*Rcut(p)**(2*i-2+l)
      Amat(3,i)=dble((2*i-1+l)*(2*i-2+l))*Rcut(p)**(2*i-3+l)
      Amat(4,i)=dble((2*i-1+l)*(2*i-2+l)*(2*i-3+l))*Rcut(p)**(2*i-4+l)
    end do
    call invers_mat(Nmax,Amat,Bmat)
    do i=1,Nmax
      Cvec(i)=sum(Bmat(i,:)*bvec(:))
    end do
    uPS(0:iRc(p),p)=Cvec(1)*rL(0:iRc(p))**(1+l) &
      +Cvec(2)*rL(0:iRc(p))**(3+l)+Cvec(3)*rL(0:iRc(p))**(5+l) &
      +Cvec(4)*rL(0:iRc(p))**(7+l)
    Cr_phi_P(:,p)=Cvec(:)
  end do
  return
end subroutine make_uPS
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_core_density
  use global_variables
  implicit none
  real(8) a,b
  integer npsi,n
  real(8) rate1,rate2,temp1,temp2

! make nc_A ---------------------------------------!
  nc_A=0d0
  npsi=0
  do n=1,Max_orbit
    if((orbit_state(n)=='v'))then
      npsi=npsi+1
    else if(orbit_state(n)=='c')then
      npsi=npsi+1
      nc_A(1:Nx)=nc_A(1:Nx) &
        +occA(npsi)/(4d0*pi)*(upsi(1:Nx,npsi)/rL(1:Nx))**2
    end if
  end do
  nc_A(0)=2d0*nc_A(1)-nc_A(2)
! end make nc_A -----------------------------------!
! make nc_P ---------------------------------------!
  nc_P=nc_A
  rate1=(Rcore-rL(iRcore))/(rL(iRcore+1)-rL(iRcore))
  rate2=(rL(iRcore+1)-Rcore)/(rL(iRcore+1)-rL(iRcore))
  
  temp1=0.5d0*exp(-xL(iRcore))/Rp*(nc_A(iRcore+1)-nc_A(iRcore-1))/Dx
  temp2=0.5d0*exp(-xL(iRcore+1))/Rp*(nc_A(iRcore+2)-nc_A(iRcore))/Dx
  b=0.5d0/Rcore*(temp1*rate2+temp2*rate1)

  temp1=nc_A(iRcore);temp2=nc_A(iRcore+1)
  a=(temp1*rate2+temp2*rate1)-b*Rcore**2

  nc_P(0:iRcore)=a+b*rL(0:iRcore)**2
! end make nc_P -----------------------------------!
  return
end subroutine make_core_density
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_projector
  use global_variables
  implicit none
  real(8) a,b
  integer npsi,nphi,p,n,ix,l,i
  real(8) rate1,rate2,temp1,temp2
  real(8),allocatable :: Vps(:),kphi(:,:)
  
  allocate(Vps(0:Nx),kphi(0:Nx,NuPR))
  
  rate1=(RcutPR-rL(iRcPR))/(rL(iRcPR+1)-rL(iRcPR))
  rate2=(rL(iRcPR+1)-RcutPR)/(rL(iRcPR+1)-rL(iRcPR))

! calc Vps ----------------------------------------!
  Vps=Vh+Vexc+Vnucl
  
  temp1=0.5d0*exp(-xL(iRcPR))/Rp*(Vps(iRcPR+1)-Vps(iRcPR-1))/Dx
  temp2=0.5d0*exp(-xL(iRcPR+1))/Rp*(Vps(iRcPR+2)-Vps(iRcPR))/Dx
  b=0.5d0/RcutPR*(temp1*rate2+temp2*rate1)
    
  temp1=Vps(iRcPR);temp2=Vps(iRcPR+1)
  a=(temp1*rate2+temp2*rate1)-b*RcutPR**2

  Vps(0:iRcPR)=a+b*rL(0:iRcPR)**2
! end calc Vps ------------------------------------!
! make projector-----------------------------------!
  do p=1,NuPR
    l=L_proj(p)
    Kphi(0:iRcPR,p)= -0.5d0*(&
      Cr_phi_P(1,p)*dble((l+1)*l)*rL(0:iRcPR)**(l-1) &
      +Cr_phi_P(2,p)*dble((3+l)*(2+l))*rL(0:iRcPR)**(1+l) &
      +Cr_phi_P(3,p)*dble((5+l)*(4+l))*rL(0:iRcPR)**(3+l) &
      +Cr_phi_P(4,p)*dble((7+l)*(6+l))*rL(0:iRcPR)**(5+l))
  end do

  nphi=0
  npsi=0
  uPR=0d0
  do n=1,Max_orbit
    if((orbit_state(n)=='v')) then
      npsi=npsi+1
      nphi=nphi+1
      
      uPR(0,nphi)=0d0
      uPR(1:iRcPR,nphi)=Kphi(1:iRcPR,nphi) &
        +(Vps(1:iRcPR) &
        +VL(1:iRcPR)*dble(L_proj(nphi)*(L_proj(nphi)+1)) &
        -esp(npsi))*uPS(1:iRcPR,nphi)
    else if(orbit_state(n)=='c')then
      npsi=npsi+1
    end if
  end do

  do n=0,3
    do i=1,Nunocc(n)
      npsi=npsi+1
      nphi=nphi+1
      uPR(0,nphi)=0d0
      uPR(1:iRcPR,nphi)=Kphi(1:iRcPR,nphi) &
        +(Vps(1:iRcPR) &
        +VL(1:iRcPR)*dble(L_proj(nphi)*(L_proj(nphi)+1)) &
        -esp(npsi))*uPS(1:iRcPR,nphi)
    end do
  end do


! end make projector-------------------------------!
  return
end subroutine make_projector
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_V_bar
  use global_variables
  implicit none
  integer npsi,nphi,n,ix,i
  real(8),allocatable :: nP_atom(:),nA_atom(:)
  real(8),allocatable :: rho_atom(:),Veff_atom(:)
  real(8),allocatable :: Vps(:)
  real(8) Q_comp,temp1,temp2,rate1,rate2,a,b
  
  allocate(nP_atom(0:Nx),nA_atom(0:Nx),rho_atom(0:Nx),Veff_atom(0:Nx))
  allocate(Vps(0:Nx))

  nP_atom=nc_P
  nA_atom=nc_A
  nphi=0
  npsi=0
  do n=1,Max_orbit
    if((orbit_state(n)=='v')) then
      npsi=npsi+1
      nphi=nphi+1
      nA_atom(1:Nx)=nA_atom(1:Nx) &
        +occA(npsi)/(4d0*pi)*(uAE(1:Nx,nphi)/rL(1:Nx))**2
      nA_atom(0)=2d0*nA_atom(1)-nA_atom(2)
      nP_atom(1:Nx)=nP_atom(1:Nx) &
        +occA(npsi)/(4d0*pi)*(uPS(1:Nx,nphi)/rL(1:Nx))**2
      nP_atom(0)=2d0*nP_atom(1)-nP_atom(2)
    else if(orbit_state(n)=='c')then
      npsi=npsi+1
    end if
  end do

  Q_comp=(sum(expXL*rL**2*(nA_atom-nP_atom))*Rp*Dx*4d0*pi-dble(ZA)) &
    /sqrt(4d0*pi)
  rho_atom(:)=nP_atom(:)+Q_comp*gL(:,0)/sqrt(4d0*pi)
  rho=nP_atom
  call Exc_Cor
  rho=rho_atom
  call Hartree(0)
  Veff=Vh+Vexc
  Vps=Veff

  rate1=(RcutPR-rL(iRcPR))/(rL(iRcPR+1)-rL(iRcPR))
  rate2=(rL(iRcPR+1)-RcutPR)/(rL(iRcPR+1)-rL(iRcPR))
  temp1=0.5d0*exp(-xL(iRcPR))/Rp*(Vps(iRcPR+1)-Vps(iRcPR-1))/Dx
  temp2=0.5d0*exp(-xL(iRcPR+1))/Rp*(Vps(iRcPR+2)-Vps(iRcPR))/Dx
  b=0.5d0/RcutPR*(temp1*rate2+temp2*rate1)

  temp1=Vps(iRcPR);temp2=Vps(iRcPR+1)
  a=(temp1*rate2+temp2*rate1)-b*RcutPR**2

  Vps(0:iRcPR)=a+b*rL(0:iRcPR)**2
  V_bar=Vps-Veff

  return
end subroutine make_V_bar
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine orthonormalize_phi_proj
  use global_variables
  implicit none
  integer i,j,li,lj,ix
  real(8) pp,pp1,pp2,max_amp,cc
  integer Node_L(0:3)

  max_amp=3.5d0
  Node_L=0
!orthonormalize
  do i=1,NuPR
    li=L_proj(i)
    Node_L(li)=Node_L(li)+1
    do j=1,i-1
      lj=L_proj(j)
      if(li==lj)then
        pp1=sum(expXL(:)*uPS(:,j)*uPR(:,i))*Rp*Dx
        pp2=sum(expXL(:)*uPS(:,i)*uPR(:,j))*Rp*Dx
        uPR(:,i)=uPR(:,i)-uPR(:,j)*pp1
        uAE(:,i)=uAE(:,i)-uAE(:,j)*pp2
        uPS(:,i)=uPS(:,i)-uPS(:,j)*pp2
      end if
    end do
    pp=sum(expXL(:)*uPS(:,i)*uPR(:,i))*Rp*Dx
    uPR(:,i)=uPR(:,i)/pp
  end do
  
!scaling
  do i=1,NuPR
    cc=0d0
    do ix=0,iRcPR
      if(abs(uPR(ix,i))>cc)cc=abs(uPR(ix,i))
    end do
    cc=max_amp/cc
    uPR(:,i)=uPR(:,i)*cc
    uAE(:,i)=uAE(:,i)/cc
    uPS(:,i)=uPS(:,i)/cc
  end do
  
  do i=1,NuPR
    li=L_proj(i)
    do j=1,NuPR
      lj=L_proj(j)
      pp=sum(expXL(:)*uPS(:,j)*uPR(:,i))*Rp*Dx
      write(*,'(4I3,e20.10)')i,j,li,lj,pp
    end do
  end do
       
  return
end subroutine orthonormalize_phi_proj
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine make_PAW_constant
  use global_variables
  implicit none
  integer Np
  
  Np=NuPR_3D
  allocate(B_c(Np,Np),I_c(Np,Np),Delta_c(Np,Np,9),K_c(9),M_c(Np,Np,9))
  allocate(N_c(9,9),J_c(Np,Np,Np,Np),C_c(Np,Np,Np,Np),X_c(Np,Np))
  allocate(dip_c(Np,Np,3))

  call calc_M_c
  call calc_N_c
  call calc_K_c
  call calc_J_c
  call calc_D_c_D_a
  call calc_I_c
  call calc_F_c
  
  call calc_A_c
  call calc_B_c
  call calc_C_c
  call calc_X_c
  return
end subroutine make_PAW_constant
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_M_c
  use global_variables
  implicit none
  integer l,m,iProj,jProj,i,j,nYlm_i,nYlm_g,nYlm_j
  real(8) s

  do l=0,2
    rho(:)=gL(:,l)
    call Hartree(l)
    do m=-l,l
      nYlm_g=nYlm(l,m)
      do i=1,NuPR_3D
        do j=1,NuPR_3D
          nYlm_i=nL_phi(i)
          nYlm_j=nL_phi(j)
          iProj=nR_phi(i)
          jProj=nR_phi(j)
          s=sum(expXL(0:iRcPR)*Vh(0:iRcPR)*uPS(0:iRcPR,iProj) &
            *uPS(0:iRcPR,jProj))*Rp*Dx
          M_c(i,j,nYlm_g)=-gaunt(nYlm_i,nYlm_j,nYlm_g)*s
        end do
      end do
    end do
  end do
  return
end subroutine calc_M_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_N_c
  use global_variables
  implicit none
  real(8) alpha,const,temp

  alpha=9d0/(Rcomp**2)
  const=0.5d0/alpha
  N_c=0d0

! l=0
  temp=(4.d0*alpha)**1.5d0/sqrt(4.d0*pi)
  temp=temp**2
  N_c(1,1)=-0.5d0*2d0*pi**1.5d0*const**2.5d0*temp
! l=1
  temp=(4d0*alpha)**2.5d0/sqrt(4d0*pi)/6.d0
  temp=temp**2
  N_c(2,2)=-0.5d0*pi**1.5d0*const**3.5d0*temp
  N_c(3,3)=-0.5d0*pi**1.5d0*const**3.5d0*temp
  N_c(4,4)=-0.5d0*pi**1.5d0*const**3.5d0*temp
! l=2
  temp=(4.d0*alpha)**3.5d0/sqrt(4.d0*pi)/60.d0
  temp=temp**2
  N_c(5,5)=-0.5d0*1.5d0*pi**1.5d0*const**4.5d0*temp
  N_c(6,6)=-0.5d0*1.5d0*pi**1.5d0*const**4.5d0*temp
  N_c(7,7)=-0.5d0*1.5d0*pi**1.5d0*const**4.5d0*temp
  N_c(8,8)=-0.5d0*1.5d0*pi**1.5d0*const**4.5d0*temp
  N_c(9,9)=-0.5d0*1.5d0*pi**1.5d0*const**4.5d0*temp
  return
end subroutine calc_N_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_K_c
  use global_variables
  implicit none
    
  K_c=0d0
  rho(:)=gL(:,0)
  call Hartree(0)
  Vh=Vh/sqrt(4d0*pi)
  K_c(1)=-sum(expXL(0:iRcore)*Vh(0:iRcore)*nc_P(0:iRcore) &
    *rL(0:iRcore)**2)*Rp*Dx*4d0*pi
  return
end subroutine calc_K_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_J_c
  use global_variables
  implicit none
  real(8),allocatable :: trho(:),srho(:)
  integer i1,i2,i3,i4,l,m
  real(8) S0,S1,S2

  allocate(trho(0:Nx),srho(0:Nx))
  
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      srho(1:Nx)=(uAE(1:Nx,nR_phi(i1))*uAE(1:Nx,nR_phi(i2)) &
        -uPS(1:Nx,nR_phi(i1))*uPS(1:Nx,nR_phi(i2)))/(rL(1:Nx)**2)
      srho(0)=2d0*srho(1)-srho(2)
      do i3=1,NuPR_3D
        do i4=1,NuPR_3D
          trho(1:Nx)=(uAE(1:Nx,nR_phi(i3))*uAE(1:Nx,nR_phi(i4)) &
            +uPS(1:Nx,nR_phi(i3))*uPS(1:Nx,nR_phi(i4)))/(rL(1:Nx)**2)
          trho(0)=2d0*trho(1)-trho(2)
          S0=0d0
          do l=0,2
            do m=-l,l
              rho=gaunt(nYlm(l,m),nL_phi(i1),nL_phi(i2))*srho
              call hartree(l)
              S1=sum(expXL(0:iRcPR)*Vh(0:iRcPR)*trho(0:iRcPR) &
                *rL(0:iRcPR)**2)*Rp*Dx*gaunt(nYlm(l,m),nL_phi(i3),nL_phi(i4))
              S0=S0+S1
            end do
          end do
          J_c(i1,i2,i3,i4)=0.5d0*S0
        end do
      end do
    end do
  end do
  return
end subroutine calc_J_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_D_c_D_a
  use global_variables
  implicit none
  real(8),allocatable :: trho(:)
  real(8) S
  integer i1,i2,l,m
  
  allocate(trho(0:Nx))
  
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      trho(1:Nx)=(uAE(1:Nx,nR_phi(i1))*uAE(1:Nx,nR_phi(i2)) &
        -uPS(1:Nx,nR_phi(i1))*uPS(1:Nx,nR_phi(i2)))/(rL(1:Nx)**2)
      trho(0)=2d0*trho(1)-trho(2)
      do l=0,2
        S=sum(expXL(:)*trho(:)*rL(:)**(l+2))*Rp*Dx
        do m=-l,l
          Delta_c(i1,i2,nYlm(l,m))=S &
            *gaunt(nL_phi(i1),nL_phi(i2),nYlm(l,m))
        end do
      end do
    end do
  end do
  
  Delta_a=(sum(expXL*(nc_A-nc_P)*rL**2)*4d0*pi*Rp*Dx-dble(ZA))/sqrt(4d0*pi)
  
  return
end subroutine calc_D_c_D_a
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_I_c
  use global_variables
  implicit none
  integer i,j,ip,jp,ix
  real(8),allocatable :: kphi_Ai(:),kphi_Aj(:)
  real(8),allocatable :: kphi_Pi(:),kphi_Pj(:),KLphi(:)
  real(8) S1,S2

  allocate(Kphi_Ai(0:Nx),Kphi_Aj(0:Nx),KLphi(0:Nx))
  allocate(Kphi_Pi(0:Nx),Kphi_Pj(0:Nx))
  
  I_c=0d0
  do i=1,NuPR_3D
    ip=nR_phi(i)
    do ix=1,Nx-1
      Kphi_Ai(ix)=0.5d0*exp(-xL(ix)) &
        /Rp*(uAE(ix+1,ip)-uAE(ix-1,ip))/Dx
      Kphi_Pi(ix)=0.5d0*exp(-xL(ix)) &
        /Rp*(uPS(ix+1,ip)-uPS(ix-1,ip))/Dx
    end do
    Kphi_Ai(0)=2d0*Kphi_Ai(1)-Kphi_Ai(2);Kphi_Ai(Nx)=0d0
    Kphi_Pi(0)=2d0*Kphi_Pi(1)-Kphi_Pi(2);Kphi_Pi(Nx)=0d0
    do j=1,NuPR_3D
      if(nL_phi(i)/=nL_phi(j))cycle
      jp=nR_phi(j)
      do ix=1,Nx-1
        Kphi_Aj(ix)=0.5d0*exp(-xL(ix)) &
          /Rp*(uAE(ix+1,jp)-uAE(ix-1,jp))/Dx
        Kphi_Pj(ix)=0.5d0*exp(-xL(ix)) &
          /Rp*(uPS(ix+1,jp)-uPS(ix-1,jp))/Dx
      end do
      Kphi_Aj(0)=2d0*Kphi_Aj(1)-Kphi_Aj(2);Kphi_Aj(Nx)=0d0
      Kphi_Pj(0)=2d0*Kphi_Pj(1)-Kphi_Pj(2);Kphi_Pj(Nx)=0d0
      
      KLphi(1:Nx)=VL(1:Nx)*(uAE(1:Nx,ip)*uAE(1:Nx,jp) &
        -uPS(1:Nx,ip)*uPS(1:Nx,jp))
      KLphi(0)=2d0*KLphi(1)-KLphi(2)
      
      S1=.5d0*sum(expXL(0:iRcPR)*(Kphi_Ai(0:iRcPR) &
        *Kphi_Aj(0:iRcPR)-Kphi_Pi(0:iRcPR)*Kphi_Pj(0:iRcPR)))*Rp*Dx
      S2=dble(L_proj(ip)*(L_proj(ip)+1))*sum(expXL(0:iRcPR) &
        *KLphi(0:iRcPR))*Rp*Dx
      I_c(i,j)=(S1+S2)
    end do
  end do

  do i=1,NuPR_3D
    do j=1,NuPR_3D
      rho(1:Nx)=(uAE(1:Nx,nR_phi(i))*uAE(1:Nx,nR_phi(j)) &
        -uPS(1:Nx,nR_phi(i))*uPS(1:Nx,nR_phi(j)))/(rL(1:Nx)**2)
      rho(0)=2d0*rho(1)-rho(2)       
      call Hartree(0)
      Vh=Vh*gaunt(1,nL_phi(i),nL_phi(j))/sqrt(4d0*pi)
      S1=sum(expXL(0:iRcPR)*rL(0:iRcPR)**2*Vh(0:iRcPR)*nc_A(0:iRcPR))*Rp*Dx*4d0*pi
      
      rho=nc_A-nc_P
      call Hartree(0)
      Vh=Vh*sqrt(4d0*pi)
      S2=sum(expXL(0:iRcPR)*Vh(0:iRcPR)*uPS(0:iRcPR,nR_phi(i)) &
        *uPS(0:iRcPR,nR_phi(j)))*Rp*Dx*gaunt(1,nL_phi(i),nL_phi(j))
      I_c(i,j)=I_c(i,j)+S1+S2
    end do
  end do

  do i=1,NuPR_3D
    do j=1,NuPR_3D
      if(nL_phi(i)/=nL_phi(j))cycle
      
      I_c(i,j)=I_c(i,j) &
        -sum(expXL(1:iRcPR)*uAE(1:iRcPR,nR_phi(i))*uAE(1:iRcPR,nR_phi(j)) &
        /rL(1:iRcPR))*Rp*Dx*dble(ZA)
      
      I_c(i,j)=I_c(i,j) &
        -sum(expXL(0:iRcPR)*uPS(0:iRcPR,nR_phi(i)) &
        *uPS(0:iRcPR,nR_phi(j))*V_bar(0:iRcPR))*Rp*Dx
    end do
  end do
  
  return
end subroutine calc_I_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_F_c
  use global_variables
  implicit none
  real(8),allocatable :: Kphi(:),KLphi(:)
  real(8) S1,S2
  integer i,p,ix

  allocate(kphi(0:Nx),KLphi(0:Nx))
  
  F_c=0d0
  p=0
  do i=1,Max_orbit
    if((orbit_state(i)=='v'))then
      p=p+1
    end if
    if(orbit_state(i)=='c')then
      p=p+1
      do ix=1,Nx-1
        Kphi(ix)=0.5d0*exp(-xL(ix)) &
          /Rp*(upsi(ix+1,p)-upsi(ix-1,p))/Dx
      end do
      Kphi(0)=2d0*Kphi(1)-Kphi(2);Kphi(Nx)=0d0
      
      KLphi(1:Nx)=VL(1:Nx)*upsi(1:Nx,p)**2
      KLphi(0)=2d0*KLphi(1)-KLphi(2)
      
      S1=.5d0*sum(expXL*Kphi**2)*Rp*Dx
      S2=dble(L_upsi(p)*(L_upsi(p)+1))*sum(expXL*KLphi)*Rp*Dx
      F_c=F_c+occ_input(i)*(S1+S2)
    end if
  end do
  
  F_c=F_c-sum(expXL*nc_A*rL)*Rp*Dx*4d0*pi*dble(ZA)
  
  rho=nc_A-nc_P
  call hartree(0)
  F_c=F_c+0.5d0*sum(expXL*Vh*(nc_A+nc_P)*rL**2)*Rp*Dx*4d0*pi
  
  F_c=F_c-sum(expXL*nc_P*v_bar*rL**2)*Rp*Dx*4d0*pi
  
  return
end subroutine calc_F_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_A_c
  use global_variables
  implicit none
  A_c=F_c+Delta_a*K_c(1)+N_c(1,1)*Delta_a**2
  return
end subroutine calc_A_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_B_c
  use global_variables
  implicit none
  integer L
  
  B_c(:,:)=I_c(:,:)+Delta_a*M_c(:,:,1)

  do L=1,9
    B_c(:,:)=B_c(:,:)+Delta_c(:,:,L)*K_c(L) &
      +2d0*Delta_c(:,:,L)*N_c(L,1)*Delta_a
  end do
  return
end subroutine calc_B_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_C_c
  use global_variables
  implicit none
  integer L1,L2,i1,i2,i3,i4,Np
  real(8),allocatable :: tC_c(:,:,:,:)
  
  Np=NuPR_3D
  allocate(tC_c(Np,Np,Np,Np))
  
  C_c=J_c
  
  do i1=1,Np
    do i2=1,Np
      do i3=1,Np
        do i4=1,Np
          do L1=1,9
            C_C(i1,i2,i3,i4)=C_c(i1,i2,i3,i4) &
              +M_c(i1,i2,L1)*Delta_c(i3,i4,L1)
          end do
        end do
      end do
    end do
  end do
  
  do i1=1,Np
    do i2=1,Np
      do i3=1,Np
        do i4=1,Np
          do L1=1,9
            do L2=1,9
              C_c(i1,i2,i3,i4)=C_c(i1,i2,i3,i4) &
                +Delta_c(i1,i2,L1)*N_c(L1,L2)*Delta_c(i3,i4,L2)
            end do
          end do
        end do
      end do
    end do
  end do
  
  tC_c=C_c
  
  do i1=1,Np
    do i2=1,Np
      do i3=1,Np
        do i4=1,Np
          C_c(i1,i2,i3,i4)=0.5d0*(tC_c(i1,i2,i3,i4)+tC_c(i3,i4,i1,i2))
        end do
      end do
    end do
  end do
    
  return
end subroutine calc_C_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_X_c
  use global_variables
  implicit none
  real(8),allocatable :: S_c(:,:),P_c(:,:)
  real(8),allocatable :: Tmat(:,:),INV(:,:)
  integer Np,i1,i2
  real(8) s
  integer i3,i4
  
  Np=NuPR_3D
  allocate(S_c(Np,Np),P_c(Np,Np))
  allocate(Tmat(Np,Np),INV(Np,Np))
  
  S_c=0d0
  P_c=0d0
    
  do i1=1,NP
    do i2=1,NP
      if(nL_phi(i1)/=nL_phi(i2))cycle
      P_c(i1,i2)=sum(expXL(:)*uPR(:,nR_phi(i1)) &
        *uPR(:,nR_phi(i2)))*Rp*Dx
      S_c(i1,i2)=sqrt(4d0*pi)*delta_c(i1,i2,1)
    end do
  end do

  do i1=1,Np
    do i2=1,Np
      Tmat(i1,i2)=sum(P_c(i1,:)*S_c(:,i2))
      if(i1==i2)Tmat(i1,i2)=Tmat(i1,i2)+1d0
    end do
  end do
  call invers_mat(Np,Tmat,INV)
  
  do i1=1,Np
    do i2=1,Np
      X_c(i1,i2)=-sum(S_c(i1,:)*INV(:,i2))
    end do
  end do
  return
end subroutine calc_X_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine calc_dip_c
  use global_variables
  implicit none
  real(8),allocatable :: trho(:)
  real(8) S
  integer i1,i2,ip1,ip2,il1,il2

  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      ip1=nR_phi(i1)
      ip2=nR_phi(i2)
      il1=nL_phi(i1)
      il2=nL_phi(i2)
      
      trho(:)=(uAE(:,ip1)*uAE(:,ip2) &
        -uPS(:,ip1)*uPS(:,ip2))
      
      S=sum(expXL*trho*rL)*Rp*Dx
      dip_c(i1,i2,1)=S*2d0*dsqrt(pi/3d0)*gaunt(2,il1,il2)
      dip_c(i1,i2,2)=S*2d0*dsqrt(pi/3d0)*gaunt(4,il1,il2)
      dip_c(i1,i2,3)=S*2d0*dsqrt(pi/3d0)*gaunt(3,il1,il2)
    end do
  end do
  return
end subroutine calc_dip_c
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine output_constant
  use global_variables
  implicit none
  integer i1,i2,i3,i4,L1,L2,p
  character(50) filename
  
  
    
! output parameter
  filename=trim(atom_name)//'_parameter.dat'
  open(10,file=filename)
  write(10,'(100I7)')NuPR,NuPR_3D,Nx,iRcPR,iRcomp,iRcf,iRcore
  write(10,'(100e25.15)')Rp,Dx,RcutPR,Rcomp,Rfilt,Rcore,softning
  do p=1,NuPR
    write(10,'(I7,2X,e25.15)')iRc(p),Rcut(p)
  end do
  do p=1,NuPR
    write(10,'(100I7)')L_proj(p)
  end do
  close(10)

! output ylm
  filename=trim(atom_name)//'_ylm.dat'
  open(10,file=filename)
  do p=1,NuPR
    write(10,'(100I7)')L_proj(p)
  end do
  close(10)

! output A_cn
  filename=trim(atom_name)//'_Aatom.dat'
  open(10,file=filename)
  write(10,'(100e25.15)')A_c
  close(10)

! output B_c
  filename=trim(atom_name)//'_Batom.dat'
  open(10,file=filename)
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      write(10,'(100e25.15)')B_c(i1,i2)
    end do
  end do
  close(10)

! output C_c
  filename=trim(atom_name)//'_Catom.dat'
  open(10,file=filename)
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      do i3=1,NuPR_3D
        do i4=1,NuPR_3D
          write(10,'(100e25.15)')C_c(i1,i2,i3,i4) 
        end do
      end do
    end do
  end do
  close(10)

! output Delta_c
  filename=trim(atom_name)//'_Delta_c.dat'
  open(10,file=filename)
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      do L1=1,9
        write(10,'(100e25.15)')Delta_c(i1,i2,L1)
      end do
    end do
  end do
  close(10)

! output Delta_a
  filename=trim(atom_name)//'_Delta_a.dat'
  open(10,file=filename)
  write(10,'(100e25.15)')Delta_a
  close(10)

! output X_c
  filename=trim(atom_name)//'_Xatom.dat'
  open(10,file=filename)
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      write(10,'(100e25.15)')X_c(i1,i2)
    end do
  end do
  close(10)

! output X_c
  filename=trim(atom_name)//'_dip.dat'
  open(10,file=filename)
  do i1=1,NuPR_3D
    do i2=1,NuPR_3D
      write(10,'(100e25.15)')dip_c(i1,i2,1),dip_c(i1,i2,2) &
        ,dip_c(i1,i2,3)
    end do
  end do
  close(10)
  
  return
end subroutine output_constant
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine output_radial_function
  use global_variables
  implicit none
  character(50) filename
  real(8),allocatable :: phi_A(:,:),phi_P(:,:),proj(:,:)
  integer Np,p,ix
  
  Np=NuPR
  allocate(phi_A(0:Nx,Np),phi_P(0:Nx,Np),proj(0:Nx,Np))
  do p=1,Np
    phi_A(1:Nx,p)=uAE(1:Nx,p)/rL(1:Nx)
    phi_A(0,p)=2d0*phi_A(1,p)-phi_A(2,p)
    phi_P(1:Nx,p)=uPS(1:Nx,p)/rL(1:Nx)
    phi_P(0,p)=2d0*phi_P(1,p)-phi_P(2,p)
    proj(1:Nx,p)=uPR(1:Nx,p)/rL(1:Nx)
    proj(0,p)=2d0*proj(1,p)-proj(2,p)
  end do

! output weight function
  filename=trim(atom_name)//'_weight.dat'
  open(10,file=filename)
  do ix=0,Nx
    write(10,'(100e25.15)')rL(ix),expXL(ix)*rL(ix)**2
  end do
  close(10)

! output partial wave
  filename=trim(atom_name)//'_phi.dat'
  open(10,file=filename)
  do ix=0,Nx
    write(10,'(100e25.15)')rL(ix),(phi_A(ix,p),phi_P(ix,p),p=1,Np)
  end do
  close(10)

! output projector
  filename=trim(atom_name)//'_proj.dat'
  open(10,file=filename)
  do ix=0,Nx
    write(10,'(100e25.15)')rL(ix),(proj(ix,p),p=1,Np)
  end do
  close(10)
  
! output core density
  filename=trim(atom_name)//'_core_density.dat'
  open(10,file=filename)
  do ix=0,Nx
    write(10,'(100e25.15)')rL(ix),nc_A(ix),nc_P(ix)
  end do
  close(10)

! output local potential
  filename=trim(atom_name)//'_Vloc.dat'
  open(10,file=filename)
  do ix=0,Nx
    write(10,'(100e25.15)')rL(ix),V_bar(ix)
  end do
  close(10)

  return
end subroutine output_radial_function
!-------10--------20--------30--------40--------50--------60----------72
subroutine filtering_rojector
  use global_variables
  implicit none
  real(8),allocatable :: Fq(:),qL(:)
  real(8) S
  real(8) qmax,dq,qc,alpha,beta,q,qr,sph_bess
  real(8) gamma
  integer Nq,p,l,ix,iq
  integer iRc_f
  
  alpha=1.1d0
  Nq=10000
  allocate(Fq(0:Nq),qL(0:Nq))
  qmax=pi/grid_size_3D
  qc=qmax/alpha
  dq=qmax/dble(Nq)
  beta=5d0/(qmax-qc)*log(10d0)
  gamma=Rfilt/RcutPR
  iRc_f=aint(log(Rfilt/Rp+1d0)/Dx)
  
  do iq=0,Nq
    qL(iq)=dq*dble(iq)
  end do

! q=dq*i i=0,Nq
! fourier transform
  do p=1,NuPR
    select case(L_proj(p))
    case(0)
      do iq=0,Nq
        S=0d0
        do ix=1,iRcPR
          qr=qL(iq)*rL(ix)
          if(qr/=0d0)then
            sph_bess=sin(qr)/qr
          else if(qr==0d0)then
            sph_bess=1d0
          end if
          S=S+expXL(ix)*uPR(ix,p)*rL(ix)*sph_bess &
            *exp(gamma*(rL(ix)/Rfilt)**2)
        end do
        Fq(iq)=S*Rp*Dx
      end do
    case(1)
      do iq=0,Nq
        S=0d0
        do ix=1,iRcPR
          qr=qL(iq)*rL(ix)
          if(qr>=0.9760d-3)then
            sph_bess=(sin(qr)-qr*cos(qr))/qr**2
          else 
            sph_bess=qr/3d0
          end if
          S=S+expXL(ix)*uPR(ix,p)*rL(ix)*sph_bess &
            *exp(gamma*(rL(ix)/Rfilt)**2)
        end do
        Fq(iq)=S*Rp*Dx
      end do
    case(2)
      do iq=0,Nq
        S=0d0
        do ix=1,iRcPR
          qr=qL(iq)*rL(ix)
          if(qr>=0.7823d-2)then
            sph_bess=((3d0-qr**2)*sin(qr)-3d0*qr*cos(qr))/qr**3
          else 
            sph_bess=qr**2/15d0
          end if
          S=S+expXL(ix)*uPR(ix,p)*rL(ix)*sph_bess &
            *exp(gamma*(rL(ix)/Rfilt)**2)
        end do
        Fq(iq)=S*Rp*Dx
      end do
    case(3)
      do iq=0,Nq
        S=0d0
        do ix=1,iRcPR
          qr=qL(iq)*rL(ix)
          if(qr>=0.3125d-1)then
            sph_bess=((15d0-6d0*qr**2)*sin(qr) &
              -qr*(15d0-qr**2)*cos(qr))/qr**4
          else 
            sph_bess=qr**3/105d0
          end if
          S=S+expXL(ix)*uPR(ix,p)*rL(ix)*sph_bess &
            *exp(gamma*(rL(ix)/Rfilt)**2)
        end do
        Fq(iq)=S*Rp*Dx
      end do
    end select
! filtering
    do iq=0,Nq
      if(qL(iq)>=qc)Fq(iq)=Fq(iq)*exp(-beta*(qL(iq)/qc-1d0)**2)
    end do
! fourier transform inverse
    select case(L_proj(p))
    case(0)
      do ix=0,iRc_f
        S=0d0
        do iq=0,Nq
          qr=qL(iq)*rL(ix)
          if(qr/=0d0)then
            sph_bess=sin(qr)/qr
          else if(qr==0d0)then
            sph_bess=1d0
          end if
          S=S+Fq(iq)*sph_bess*qL(iq)**2
        end do
        uPR(ix,p)=S*dq*2d0/pi*rL(ix)*exp(-gamma*(rL(ix)/Rfilt)**2)
      end do
    case(1)
      do ix=0,iRc_f
        S=0d0
        do iq=0,Nq
          qr=qL(iq)*rL(ix)
          if(qr>=0.9760d-3)then
            sph_bess=(sin(qr)-qr*cos(qr))/qr**2
          else 
            sph_bess=qr/3d0
          end if
          S=S+Fq(iq)*sph_bess*qL(iq)**2
        end do
        uPR(ix,p)=S*dq*2d0/pi*rL(ix)*exp(-gamma*(rL(ix)/Rfilt)**2)
      end do
    case(2)
      do ix=0,iRc_f
        S=0d0
        do iq=0,Nq
          qr=qL(iq)*rL(ix)
          if(qr>=0.7823d-2)then
            sph_bess=((3d0-qr**2)*sin(qr)-3d0*qr*cos(qr))/qr**3
          else 
            sph_bess=qr**2/15d0
          end if
          S=S+Fq(iq)*sph_bess*qL(iq)**2
        end do
        uPR(ix,p)=S*dq*2d0/pi*rL(ix)*exp(-gamma*(rL(ix)/Rfilt)**2)
      end do
    case(3)
      do ix=0,iRc_f
        S=0d0
        do iq=0,Nq
          qr=qL(iq)*rL(ix)
          if(qr>=0.3125d-1)then
            sph_bess=((15d0-6d0*qr**2)*sin(qr) &
              -qr*(15d0-qr**2)*cos(qr))/qr**4
          else 
            sph_bess=qr**3/105d0
          end if
          S=S+Fq(iq)*sph_bess*qL(iq)**2
        end do
        uPR(ix,p)=S*dq*2d0/pi*rL(ix)*exp(-gamma*(rL(ix)/Rfilt)**2)
      end do
    end select
  end do
  return
end subroutine filtering_rojector
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine factorial(input,output)
  integer input
  real(8) output
  integer i

  output=1d0
  do i=1,input
    output=output*dble(i)
  end do
  return
end subroutine factorial
!-------10--------20--------30--------40--------50--------60--------70--------80
subroutine invers_mat(Nmax,Amat,INV)
  implicit none
! Amat^-1=INV
  integer Nmax
  real(8) Amat(Nmax,Nmax)
  real(8) Hmat(Nmax,Nmax),INV(Nmax,Nmax)
  real(8) Cvec(Nmax),Cvec2(Nmax)
  real(8) d1,d2,f1,f2
  integer nd1,nd2,i,j,k
  
  INV=0d0
  Hmat=Amat
  do i=1,Nmax
    INV(i,i)=1.d0
  end do
    
  do i=1,Nmax
 
! pivot     
    nd1=i
    d1=dabs(Hmat(i,i))
    do j=i+1,Nmax
      if(abs(Hmat(i,j))>d1)then
        d1=dabs(Hmat(i,j))
        nd1=j
      end if
    end do
    Cvec(:)=Hmat(:,i)
    Cvec2(:)=INV(:,i)
    
    Hmat(:,i)=Hmat(:,nd1)
    INV(:,i)=INV(:,nd1)
    
    Hmat(:,nd1)=Cvec(:)
    INV(:,nd1)=Cvec2(:)
      
! deleat
    do j=1,Nmax
      if(i/=j)then
        f1=Hmat(i,j)/Hmat(i,i)
        Hmat(:,j)=Hmat(:,j)-f1*Hmat(:,i)
        INV(:,j)=INV(:,j)-f1*INV(:,i)
      end if
    end do
    
    d1=1.d0/Hmat(i,i)
    Hmat(:,i)=Hmat(:,i)*d1
    INV(:,i)=INV(:,i)*d1
    
  end do
  return
end subroutine invers_mat
!-------10--------20--------30--------40--------50--------60----------72
