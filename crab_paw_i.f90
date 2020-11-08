 module global_variables
    include 'mpif.h'

! constants
    real(8),parameter :: pi=3.1415926535897932d0
    real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
    real(8),parameter :: sq4p=dsqrt(4d0*pi)
    complex(8),parameter :: zI=(.0d0,1d0)

! DFT parameter 
    real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
    real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
    real(8),parameter :: CU=0.002d0,DU=-0.0116d0

! parallelization
    integer Npx,Npy,Npz,Nprocs,myrank,my_px,my_py,my_pz
    integer rsize_x,rsize_y,rsize_z
    integer,allocatable :: Itable(:,:,:)
    real(8),allocatable :: sendxp(:,:,:),sendxm(:,:,:),sendyp(:,:,:)
    real(8),allocatable :: sendym(:,:,:),sendzp(:,:,:),sendzm(:,:,:)
    real(8),allocatable :: recvxp(:,:,:),recvxm(:,:,:),recvyp(:,:,:)
    real(8),allocatable :: recvym(:,:,:),recvzp(:,:,:),recvzm(:,:,:)
    complex(8),allocatable :: zsendxp(:,:,:),zsendxm(:,:,:)
    complex(8),allocatable :: zsendyp(:,:,:),zsendym(:,:,:)
    complex(8),allocatable :: zsendzp(:,:,:),zsendzm(:,:,:)
    complex(8),allocatable :: zrecvxp(:,:,:),zrecvxm(:,:,:)
    complex(8),allocatable :: zrecvyp(:,:,:),zrecvym(:,:,:)
    complex(8),allocatable :: zrecvzp(:,:,:),zrecvzm(:,:,:)

! grid
    integer NLx,NLy,NLz,NLpx,NLpy,NLpz,Nd
    integer Nsx,Nsy,Nsz,Nex,Ney,Nez
    real(8) H
    real(8),allocatable :: xL(:),yL(:),zL(:)
    real(8) length_x,length_y,length_z

! finite difference
    real(8),parameter :: cN0=-205d0/72d0,cN1=8d0/5d0
    real(8),parameter :: cN2=-1d0/5d0,cN3=8d0/315d0
    real(8),parameter :: cN4=-1d0/560d0    

! Gauss-Lgendre integral      
    real(8) GLz(20),GLome(20),GLYlm(20,20,9)
! table
    real(8),allocatable :: l_azim(:),m_magn(:)
! PAW parameter
! paw constant
    integer maxproj,maxprojR,irmax
    real(8),allocatable :: softning(:)
    real(8),allocatable :: Ac(:),Bc(:,:,:),Cc(:,:,:,:,:)
    real(8),allocatable :: Vc(:,:,:,:),Dc(:,:,:),Qc(:,:)
    real(8),allocatable :: Delta_c(:,:,:,:),Delta_a(:)
    real(8),allocatable :: Hc(:,:,:),Sc(:,:,:),Xc(:,:,:)
    real(8),allocatable :: Dmu_c(:,:,:),Wc(:,:),Exc_c(:)
    real(8),allocatable :: rk_c(:,:,:)
    complex(8),allocatable :: dip_c(:,:,:,:)
!! radial function
    real(8),allocatable :: phiA(:,:,:),phiP(:,:,:),projR(:,:,:)      
    real(8),allocatable :: ncA(:,:),ncP(:,:),Vbar_a(:,:)
    real(8),allocatable :: nacA(:,:,:,:),nacP(:,:,:,:)
    real(8),allocatable :: rL_atom(:,:),expXr2(:,:)
!! parameter
    real(8),allocatable :: Rcut(:),H_atom(:),Rfilt(:)
    real(8),allocatable :: Rcore(:),Rcomp(:)
    integer,allocatable :: iRcut(:),iRfilt(:),iRcore(:)
    integer,allocatable :: Npw(:),lpw(:,:),Nproj(:)
    integer,allocatable :: NLproj(:,:),NRproj(:,:),NL_R(:)
    real(8),allocatable :: Rp_atom(:)
    real(8) faclog(500),clebma
    integer,allocatable :: iproj(:,:,:)


! material
    integer NI,NE,NST,NSTocc
    integer,allocatable :: Kion(:)
    real(8) Etot
    real(8),allocatable :: Rion(:,:)
    real(8),allocatable :: occ(:),esp(:)

! common
    real(8),allocatable :: nP(:,:,:),rho(:,:,:)
    real(8),allocatable :: nP_core(:,:,:)
    real(8),allocatable :: proj_P(:,:,:,:,:),g_hat(:,:,:,:,:)
    real(8),allocatable :: V_hat(:,:,:,:,:),Vbar(:,:,:)
    real(8),allocatable :: Veff(:,:,:),Vhat_sum(:,:,:)
    real(8),allocatable :: Ylm(:,:,:,:,:)

! gs
    integer Ncg,Niter,NCiter
    real(8),allocatable :: psi(:,:,:,:),Pc(:,:,:)
    real(8),parameter :: mixingrate=0.2d0                 !simple mixing
    real(8),allocatable :: Veff_old(:,:,:),Hc_old(:,:,:)  !simple mixing

! rt
    integer NTiter,iter_rt
    real(8) f0,omega,tt
    complex(8),allocatable :: zpsi(:,:,:,:),zPc(:,:,:)
    real(8) DT,ex,ey,ez,rk_mom,pol_theta,pol_phi
    real(8) W0,deltaW,Wcm2,pulse_time_fs,wave_length,ft
    real(8),allocatable :: Vext(:,:,:),Wabs(:,:,:)
    real(8) dipole_td
    real(8) length_R

! Hartree
    real(8),allocatable :: Vh(:,:,:)

! Exc_Corr potential
    real(8),allocatable :: Vexc(:,:,:),Eexc(:,:,:)

! temporary
    integer iterVh
    real(8),allocatable :: wk(:,:,:),Lwk(:,:,:)
    complex(8),allocatable :: zwk(:,:,:),zLwk(:,:,:)
    real(8),allocatable :: tpsi(:,:,:),htpsi(:,:,:),Stpsi(:,:,:)
    complex(8),allocatable :: ztpsi(:,:,:),zXHtpsi(:,:,:),zStpsi(:,:,:)
    
! job control,I/O,files
    character(10) rt_option
    character(10) absorber,Bcondition
    character(50) SYSname,file_laser,file_dipole
    character(50) file_elec_all,file_energy
    character(10) predict_corr
    character(50),allocatable :: atom_name(:)

  end module global_variables
!-------10--------20--------30--------40--------50--------60----------72
  program main
    use global_variables
    implicit none
    integer iter,ierr


    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,Myrank,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!-----preparation for calculation


    call prep_calc_parameter
    call init_wf
    call Gram_Schmidt
    call psi_rho
    call local_potential
    call nonlocal_potential
    Veff_old=Veff;Hc_old=Hc
! GS---------------!
    do iter=1,Niter
    if(myrank==0)write(*,*)'iter=',iter
    call simple_mixing
    call DTcg
    call Gram_Schmidt
    call diag
    call psi_rho
    call local_potential
    call nonlocal_potential
    call total_energy('gs')
    end do
! end GS------------!

! end rt------------!
    call prep_rt
    call particle_number_projection
    call zpsi_rho
    iter_rt=0
    tt=Dt*dble(iter_rt)
    if(Myrank == 0)then
       write(14,*)iter_rt*DT*0.02419d0,ft
    end if
    call calc_dipole_moment
    call calc_particle_number
    call total_energy('rt')
    do iter_rt=1,NTiter
       call RTI
       tt=Dt*dble(iter_rt)
       if(Myrank == 0)then
          write(14,*)iter_rt*DT*0.02419d0,ft
       end if
       call calc_dipole_moment
       call calc_particle_number
       call total_energy('rt')
    end do

    call particle_number_projection

    if(Myrank == 0)then
       close(22)
       close(8)
       close(12)
       close(14)
    end if
    
    if(myrank==0)write(*,*)'END PROGRAM'
    call MPI_FINALIZE(ierr)

  end program main

!-------10--------20--------30--------40--------50--------60----------72
  subroutine prep_calc_parameter
    use global_variables
    implicit none
    integer ierr,a,j


    if(myrank==0)then
       read(*,*)SYSname
       read(*,*)Ncg
       read(*,*)Niter
       read(*,*)rt_option
       read(*,*)NTiter,DT
       read(*,*)rk_mom
       read(*,*)pol_theta
       read(*,*)pol_phi
       read(*,*)Bcondition
       read(*,*)absorber
       read(*,*)predict_corr
       read(*,*)W0,deltaW
       read(*,*)Wcm2,pulse_time_fs,wave_length
       read(*,*)H
       read(*,*)NLx,NLy,NLz
       read(*,*)Npx,Npy,Npz
       read(*,*)Nd
       read(*,*)NST,NSTocc
       read(*,*)NI,NE
    end if
    call MPI_BCAST(NI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    allocate(atom_name(NE))
    allocate(Rion(3,NI),Kion(NI))
    
    if(myrank==0)then
       do a=1,NE
          read(*,*)j,atom_name(a)
       end do
       do a=1,NI
          read(*,*)j,Rion(1,a),Rion(2,a),Rion(3,a),Kion(a)
       end do
    end if

!-------Broadcast input parameter-----------------------
    call MPI_BCAST(SYSname,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Ncg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Niter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rt_option,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NTiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(DT,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rk_mom,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(pol_theta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(pol_phi,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Bcondition,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(absorber,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(predict_corr,10,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(W0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(deltaW,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Wcm2,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(pulse_time_fs,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(wave_length,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(H,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NLx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NLy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NLz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Npx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Npy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Npz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NSTocc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rion,3*NI,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Kion,NI,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      if(Myrank == 0)then
        file_dipole  = trim(SYSname)//'_Dt_p.out'
        file_laser = trim(SYSname)//'_laser_p.out'
        file_energy= trim(SYSname)//'_energy_p.out'
        file_elec_all=trim(SYSname)//'_elec_all.out'
      endif
    
      call MPI_BCAST(file_dipole,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(file_energy,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(file_laser,50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(file_elec_all,50,MPI_CHARACTER &
           ,0,MPI_COMM_WORLD,ierr)

    call init_parallel

    call input_paw_parameter

    allocate(psi(Nsx:Nex,Nsy:Ney,Nsz:Nez,1:NST))
    allocate(nP(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(rho(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(nP_core(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(proj_P(Nsx:Nex,Nsy:Ney,Nsz:Nez,maxproj,NI))
    allocate(g_hat(Nsx:Nex,Nsy:Ney,Nsz:Nez,9,NI))
    allocate(V_hat(Nsx:Nex,Nsy:Ney,Nsz:Nez,9,NI))
    allocate(Vbar(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Vh(Nsx:Nex,Nsy:Ney,Nsz:Nez),Veff(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Veff_old(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Hc_old(maxproj,maxproj,NI))
    allocate(Vexc(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Eexc(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Vhat_sum(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Ylm(Nsx-Nd:Nex+Nd,Nsy-Nd:Ney+Nd,Nsz-Nd:Nez+nd,0:4,-4:4))
    allocate (tpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Stpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (htpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (wk(Nsx-Nd:Nex+Nd,Nsy-Nd:Ney+Nd,Nsz-Nd:Nez+Nd), &
         Lwk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (zwk(Nsx-Nd:Nex+Nd,Nsy-Nd:Ney+Nd,Nsz-Nd:Nez+Nd), &
         zLwk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    wk=0d0
    zwk=0d0
    Vh=0d0
    call Make_Ylm
    call prep_paw

    ex=sin(pol_theta/180d0*Pi)*cos(pol_phi/180d0*Pi)
    ey=sin(pol_theta/180d0*Pi)*sin(pol_phi/180d0*Pi)
    ez=cos(pol_theta/180d0*Pi)

    allocate(occ(NST),esp(NST))
    occ=0d0
    occ(1:NSTocc)=2d0
    return
  end subroutine prep_calc_parameter
!-------10--------20--------30--------40--------50--------60----------72
  subroutine init_parallel
    use global_variables
    implicit none
    integer px,py,pz,irank,ierr
    integer ix,iy,iz
    NLpx=NLx/Npx
    NLpy=NLy/Npy
    NLpz=NLz/Npz
    if(Myrank == 0) write(*,*) 'NLpx,NLpy,NLpz',NLpx,NLpy,NLpz
    if((NLx /= NLpx*Npx).or.(NLy /= NLpy*Npy) &
         .or.(NLz /= NLpz*Npz))then
       call MPI_FINALIZE(ierr)
       stop 'NLx /= NLpx*Npx'
    end if

    allocate(Itable(-1:Npx,-1:Npy,-1:Npz))
    do px=-1,Npx
    do py=-1,Npy
    do pz=-1,Npz
       Itable(px,py,pz)=MPI_PROC_NULL
    enddo
    enddo
    enddo

    irank=0
    do px=0,Npx-1
    do py=0,Npy-1
    do pz=0,Npz-1
       Itable(px,py,pz)=irank
       if(Myrank == irank) then
          My_px=px
          My_py=py
          My_pz=pz
       endif
       irank=irank+1
    enddo
    enddo
    enddo

    Nsx=My_px*NLpx+1
    Nsy=My_py*NLpy+1
    Nsz=My_pz*NLpz+1
    Nex=(My_px+1)*NLpx
    Ney=(My_py+1)*NLpy
    Nez=(My_pz+1)*NLpz

!--- Set size for zLaplacan, Laplacian
    rsize_x=4*(Ney-Nsy+1)*(Nez-Nsz+1)
    rsize_y=(Nex-Nsx+1)*4*(Nez-Nsz+1)
    rsize_z=(Nex-Nsx+1)*(Ney-Nsy+1)*4

    length_x=(NLx+1)*dble(H)
    length_y=(NLy+1)*dble(H)
    length_z=(NLz+1)*dble(H)

    allocate(xL(Nsx-Nd:Nex+Nd),yL(Nsy-Nd:Ney+Nd),zL(Nsz-Nd:Nez+Nd))

    do ix=Nsx-Nd,Nex+Nd
       xL(ix)=dble(ix)*H-length_x/2d0
    enddo
    do iy=Nsy-Nd,Ney+Nd
       yL(iy)=iy*dble(H)-length_y/2d0
    enddo
    do iz=Nsz-Nd,Nez+Nd
       zL(iz)=iz*dble(H)-length_z/2d0
    enddo

    allocate(sendxp(1:Nd,    Nsy:Ney, Nsz:Nez),&
         sendxm(1:Nd,    Nsy:Ney, Nsz:Nez),&
         sendyp(Nsx:Nex, 1:Nd,    Nsz:Nez),&
         sendym(Nsx:Nex, 1:Nd,    Nsz:Nez),&
         sendzp(Nsx:Nex, Nsy:Ney, 1:Nd   ),&
         sendzm(Nsx:Nex, Nsy:Ney, 1:Nd   ),&
         recvxp(1:Nd,    Nsy:Ney, Nsz:Nez),&
         recvxm(1:Nd,    Nsy:Ney, Nsz:Nez),&
         recvyp(Nsx:Nex, 1:Nd,    Nsz:Nez),&
         recvym(Nsx:Nex, 1:Nd,    Nsz:Nez),&
         recvzp(Nsx:Nex, Nsy:Ney, 1:Nd   ),&
         recvzm(Nsx:Nex, Nsy:Ney, 1:Nd   ) )

    allocate(zsendxp(1:Nd,    Nsy:Ney, Nsz:Nez)&
         ,zsendxm(1:Nd,    Nsy:Ney, Nsz:Nez)&
         ,zsendyp(Nsx:Nex, 1:Nd,    Nsz:Nez)&
         ,zsendym(Nsx:Nex, 1:Nd,    Nsz:Nez)&
         ,zsendzp(Nsx:Nex, Nsy:Ney, 1:Nd  )&
         ,zsendzm(Nsx:Nex, Nsy:Ney, 1:Nd  )&
         ,zrecvxp(1:Nd,    Nsy:Ney, Nsz:Nez)&
         ,zrecvxm(1:Nd,    Nsy:Ney, Nsz:Nez)&
         ,zrecvyp(Nsx:Nex, 1:Nd,    Nsz:Nez)&
         ,zrecvym(Nsx:Nex, 1:Nd,    Nsz:Nez)&
         ,zrecvzp(Nsx:Nex, Nsy:Ney, 1:Nd  )&
         ,zrecvzm(Nsx:Nex, Nsy:Ney, 1:Nd  ) )

    return
  end subroutine init_parallel
!-------10--------20--------30--------40--------50--------60----------72
  subroutine input_paw_parameter
    use global_variables
    implicit none
    real(8) r
    integer ierr,a,p,i,j,i1,i2,j1,j2,ir,l,m
    character(50) filename


! allocate
!!   parameter
    allocate (Rcut(NE),iRcut(NE),iRfilt(NE),NL_R(NE),softning(NE))
    allocate (Rcore(NE),iRcore(NE),Rcomp(NE))
    allocate (Npw(NE),H_atom(NE),Nproj(NE),Rp_atom(NE),Rfilt(NE))


    if(myrank == 0)then
       do a=1,NE
         filename=trim(atom_name(a))//'_parameter.dat'
         open(20,file=filename)
         read(20,*)Npw(a),Nproj(a),NL_R(a),iRcut(a),iRfilt(a),iRcore(a)
         read(20,*)Rp_atom(a),H_atom(a),Rcut(a),Rcomp(a),Rfilt(a),Rcore(a),softning(a)
       end do
       maxproj=maxval(Nproj)
       maxprojR=maxval(Npw)
       irmax=maxval(NL_R)
     end if


      call MPI_BCAST(maxproj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(maxprojR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        
! allocate
!!   paw constant
      allocate (Ac(NE),Bc(maxproj,maxproj,NE))
      allocate (Cc(maxproj,maxproj,maxproj,maxproj,NE))
      allocate (Vc(9,9,NI,NI),Dc(maxproj,maxproj,NI),Qc(9,NI))
      allocate (Delta_c(maxproj,maxproj,9,NE),Delta_a(NE))
      allocate (Sc(maxproj,NST,NE),Hc(maxproj,maxproj,NI))
      allocate (Dmu_c(maxproj,maxproj,NI),Wc(9,NI),Exc_c(NI))
      allocate (Xc(maxproj,maxproj,NE),dip_c(maxproj,maxproj,3,NE))
      allocate (Pc(maxproj,NST,NI))
!!   radial function
      allocate (phiA(0:irmax,maxprojR,NE),phiP(0:irmax,maxprojR,NE))
      allocate (projR(0:irmax,maxprojR,NE))
      allocate (ncA(0:irmax,NE),ncP(0:irmax,NE),Vbar_a(0:irmax,NE))
      allocate (nacA(20,20,0:irmax,NI),nacP(20,20,0:irmax,NI))
      allocate (rL_atom(0:irmax,NE),expXr2(0:irmax,NE))
!!   parameter
      allocate (NLproj(maxproj,NE),NRproj(maxproj,NE))
      allocate (iproj(2,3,NI),lpw(maxprojR,NE))
      Ac=0d0;Bc=0d0;Cc=0d0;Vc=0d0;Delta_c=0d0
      Delta_a=0d0;Sc=0d0;Hc=0d0;Dmu_c=0d0;Wc=0d0;Xc=0d0;dip_c=0d0
      Pc=0d0;Dc=0d0

    if(myrank == 0)then
      do a=1,NE
         filename=trim(atom_name(a))//'_ylm.dat'
         open(20,file=filename)
         do j=1,Npw(a)
           read(20,*)lpw(j,a)
         end do
         close(20)

        filename=trim(atom_name(a))//'_weight.dat'
        open(20,file=filename)
        do ir=0,irmax
          read(20,*)rL_atom(ir,a),expXr2(ir,a)
        end do
        close(20)

        filename=trim(atom_name(a))//'_phi.dat'
        open(20,file=filename)
        do ir=0,irmax
          read(20,*)rL_atom(ir,a), &
            (phiA(ir,p,a),phiP(ir,p,a),p=1,Npw(a))
        end do
        close(20)

        filename=trim(atom_name(a))//'_proj.dat'
        open(20,file=filename)
        do ir=0,irmax
          read(20,*)r,(projR(ir,p,a),p=1,Npw(a))
        end do
        close(20)
        filename=trim(atom_name(a))//'_core_density.dat'      
        open(20,file=filename)
        do ir=0,irmax
          read(20,*)r,ncA(ir,a),ncP(ir,a)
        end do
        close(20)
        filename=trim(atom_name(a))//'_Vloc.dat'
        open(20,file=filename)
        do ir=0,irmax
          read(20,*)r,Vbar_a(ir,a)
        end do
        close(20)

        filename=trim(atom_name(a))//'_Aatom.dat'
        open(10,file=filename)
        read(10,*)Ac(a)
        close(10)
        filename=trim(atom_name(a))//'_Catom.dat'
        open(10,file=filename)
        do i1=1,Nproj(a)
          do j1=1,Nproj(a)
            do i2=1,Nproj(a)
              do j2=1,Nproj(a)
                read(10,*)Cc(i1,j1,i2,j2,a)
              end do
            end do
          end do
        end do
        close(10)
        filename=trim(atom_name(a))//'_Batom.dat'       
        open(10,file=filename)
        do i1=1,Nproj(a)
          do j1=1,Nproj(a)
            read(10,*)Bc(i1,j1,a)
          end do
        end do
        close(10)
        filename=trim(atom_name(a))//'_Delta_c.dat'       
        open(10,file=filename)
        do i1=1,Nproj(a)
          do j1=1,Nproj(a)
            do l=1,9
              read(10,*)Delta_c(i1,j1,l,a)
            end do
          end do
        end do
        close(10)
        filename=trim(atom_name(a))//'_Delta_a.dat'                
        open(10,file=filename)
        read(10,*)Delta_a(a)
        close(10)
        filename=trim(atom_name(a))//'_Xatom.dat'       
        open(10,file=filename)
        do i1=1,Nproj(a)
          do j1=1,Nproj(a)
            read(10,*)Xc(i1,j1,a)
          end do
        end do
        close(10)
        filename=trim(atom_name(a))//'_dip.dat'       
        open(10,file=filename)
        do i1=1,Nproj(a)
          do j1=1,Nproj(a)
            read(10,*)dip_c(i1,j1,1,a),dip_c(i1,j1,2,a),dip_c(i1,j1,3,a)
          end do
        end do
        close(10)
      end do
    end if

    call MPI_BCAST(Npw,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nproj,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NL_R,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iRcut,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iRfilt,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iRcore,NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rp_atom,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(H_atom,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rcut,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rcomp,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rfilt,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Rcore,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(softning,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(lpw,maxprojR*NE,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phiA,(irmax+1)*maxprojR*NE,MPI_REAL8,0, &
      MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phiP,(irmax+1)*maxprojR*NE,MPI_REAL8,0, &
      MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rL_atom,(irmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(expXr2,(irmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phiP,(irmax+1)*maxprojR*NE,MPI_REAL8,0, &
      MPI_COMM_WORLD,ierr)
    call MPI_BCAST(projR,(irmax+1)*maxprojR*NE,MPI_REAL8,0, &
      MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ncA,(irmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ncP,(irmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Vbar_a,(irmax+1)*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Ac,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Cc,maxproj**4*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Bc,maxproj**2*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Delta_c,maxproj**2*NE*9,MPI_REAL8,0, &
      MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Delta_a,NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Xc,maxproj**2*NE,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dip_c,maxproj**2*NE*3,MPI_DOUBLE_COMPLEX,0, &
      MPI_COMM_WORLD,ierr)
    do a=1,NE
       i=0
       do j=1,Npw(a)
          do m=lpw(j,a),-lpw(j,a),-1
             i=i+1
             if(lpw(j,a)==0)NLproj(i,a)=1
             if(lpw(j,a)==1 .and. m==1)NLproj(i,a)=2
             if(lpw(j,a)==1 .and. m==0)NLproj(i,a)=3
             if(lpw(j,a)==1 .and. m==-1)NLproj(i,a)=4
             if(lpw(j,a)==2 .and. m==2)NLproj(i,a)=5
             if(lpw(j,a)==2 .and. m==1)NLproj(i,a)=6
             if(lpw(j,a)==2 .and. m==0)NLproj(i,a)=7
             if(lpw(j,a)==2 .and. m==-1)NLproj(i,a)=8
             if(lpw(j,a)==2 .and. m==-2)NLproj(i,a)=9
             NRproj(i,a)=j
          end do
       end do
    end do
    return
  end subroutine input_paw_parameter
!-------10--------20--------30--------40--------50--------60----------72
!=======================================================================
! This function Ylm is related to the real spherical harmonics Ylm0 by
! Ylm=sqrt(4*pi/(2l+1))*r**l*Ylm0 and is a monomial of x,y,z

  SUBROUTINE Make_Ylm
    use Global_Variables
    implicit none
    real(8) x,y,z
    real(8) r2
    integer ix,iy,iz
    
    Ylm=0d0
    do iz = Nsz-Nd, Nez+Nd
       do iy = Nsy-Nd, Ney+Nd
          do ix = Nsx-Nd, Nex+Nd
             x = xL(ix)
             y = yL(iy)
             z = zL(iz)
             r2=x*x+y*y+z*z
             Ylm(ix,iy,iz,0,0)=1.d0
             Ylm(ix,iy,iz,1,-1)=-y
             Ylm(ix,iy,iz,1,0)=z
             Ylm(ix,iy,iz,1,1)=-x
             Ylm(ix,iy,iz,2,-2)=sqrt(3.d0)*x*y
             Ylm(ix,iy,iz,2,-1)=-sqrt(3.d0)*y*z
             Ylm(ix,iy,iz,2,0)=(3*z*z-r2)/2.d0
             Ylm(ix,iy,iz,2,1)=-sqrt(3.d0)*x*z
             Ylm(ix,iy,iz,2,2)=sqrt(3.d0/4.d0)*(x*x-y*y)
             Ylm(ix,iy,iz,3,-3)=-sqrt(5.d0/8.d0)*y*(3*x*x-y*y)
             Ylm(ix,iy,iz,3,-2)=sqrt(15.d0)*x*y*z
             Ylm(ix,iy,iz,3,-1)=-sqrt(3.d0/8.d0)*y*(5*z*z-r2)
             Ylm(ix,iy,iz,3,0)=z*(5*z*z-3*r2)/2.d0
             Ylm(ix,iy,iz,3,1)=-sqrt(3.d0/8.d0)*x*(5*z*z-r2)
             Ylm(ix,iy,iz,3,2)=sqrt(15.d0/4.d0)*z*(x*x-y*y)
             Ylm(ix,iy,iz,3,3)=-sqrt(5.d0/8.d0)*x*(x*x-3*y*y)
             Ylm(ix,iy,iz,4,-4)=sqrt(35.d0)/2.d0*x*y*(x*x-y*y)
             Ylm(ix,iy,iz,4,-3)=-sqrt(35.d0/8.d0)*y*z*(3*x*x-y*y)
             Ylm(ix,iy,iz,4,-2)=sqrt(5.d0)/2.d0*x*y*(7*z*z-r2)
             Ylm(ix,iy,iz,4,-1)=-sqrt(5.d0/8.d0)*y*z*(7*z*z-3*r2)
             Ylm(ix,iy,iz,4,0)=(35*z**4-30*z*z*r2+3.d0*r2*r2)/8.d0
             Ylm(ix,iy,iz,4,1)=-sqrt(5.d0/8.d0)*x*z*(7*z**2-3*r2)
             Ylm(ix,iy,iz,4,2)=sqrt(5.d0)/4.d0*(7*z*z-r2)*(x*x-y*y)
             Ylm(ix,iy,iz,4,3)=-sqrt(35.d0/8.d0)*x*z*(x*x-3*y*y)
             Ylm(ix,iy,iz,4,4)=sqrt(35.d0)/8.d0*(x**4+y**4-6*x*x*y*y)
          end do
       end do
    end do
    return
  end SUBROUTINE Make_Ylm
!-------10--------20--------30--------40--------50--------60----------72
  subroutine prep_paw
    use global_variables
    implicit none
    integer ix,iy,iz,ir,ierr
    integer nr,na,l1,m1,L,a,iRc
    integer irmax_g
    real(8) x,y,z,r,theta,phi,rate,alpha,alpha2,const,const2,dr,f1
    real(8) rate1,rate2,dr_g
    real(8) Rkaizyo,Yfunc
    real(8),allocatable :: gdens(:),pot(:)

    irmax_g=2000
    dr_g=0.01d0
    allocate(gdens(irmax_g),pot(irmax_g))


    call set_angular_momentum_table
    call gauss_LUgendre(GLz,GLome)
    call set_yfunc_d
    call calc_Vc

! set 3D projector       
  do a=1,NI
    na=Kion(a)
    L=0
    do nr=1,Npw(na)
      l1=lpw(nr,na)
      do m1=l1,-l1,-1
        L=L+1
        do iz=Nsz,Nez
        z=zL(iz)-Rion(3,a)
        do iy=Nsy,Ney
        y=yL(iy)-Rion(2,a)
        do ix=Nsx,Nex
        x=xL(ix)-Rion(1,a)
        r=sqrt(x*x+y*y+z*z)
        ir=aint(log(r/Rp_atom(na)+1d0)/H_atom(na))

        if(ir < NL_R(na))then
          rate1=(r-rL_atom(ir,na))/(rL_atom(ir+1,na)-rL_atom(ir,na))
          rate2=(rL_atom(ir+1,na)-r)/(rL_atom(ir+1,na)-rL_atom(ir,na))
          theta=acos(z/(r+1d-30))
          if(y >= 0d0)then
            phi=acos(x/(sqrt(x**2+y**2)+1d-30))
          else
            phi=2.d0*pi-acos(x/(sqrt(x**2+y**2)+1d-30))
          end if
          proj_P(ix,iy,iz,L,a)=(projR(ir,nr,na)*rate2 &
            &+projR(ir+1,nr,na)*rate1)*Yfunc(theta,phi,l1,m1)
        else
          proj_P(ix,iy,iz,L,a)=0d0          
        end if
        end do
        end do
        end do
      end do
    end do
  end do
! set g_hat
     do a=1,NI
        na=Kion(a)
        alpha=9.d0/softning(na)/(Rcut(na)**2)
        do L=1,9
           l1=l_azim(L)
           m1=m_magn(L)
           const=Rkaizyo(l1)/Rkaizyo(2*l1+1)&
                /sqrt(4.d0*pi)*(4.d0*alpha)**(dble(l1)+1.5d0)
           do iz=Nsz,Nez
           z=zL(iz)-Rion(3,a)
           do iy=Nsy,Ney
           y=yL(iy)-Rion(2,a)
           do ix=Nsx,Nex
           x=xL(ix)-Rion(1,a)
           r=sqrt(x*x+y*y+z*z)
           theta=acos(z/(r+1d-30))
           if(y >= (.0d0))then
              phi=acos(x/(sqrt(x*x+y*y)+1d-30))
           else
              phi=2.d0*pi-acos(x/(sqrt(x*x+y*y)+1d-30))
           end if
           g_hat(ix,iy,iz,L,a)=const*r**l1*exp(-alpha*r*r)&
                *Yfunc(theta,phi,l1,m1)
           end do
           end do
           end do
        end do
     end do
! set V_hat
     do a=1,NI
        na=Kion(a)
        alpha=9.d0/(Rcut(na)**2)
        alpha2=9.d0/softning(na)/(Rcut(na)**2)
        do L=1,9
           l1=l_azim(L)
           m1=m_magn(L)
           const=Rkaizyo(l1)/Rkaizyo(2*l1+1)&
                /sqrt(4.d0*pi)*(4.d0*alpha)**(dble(l1)+1.5d0) 
           const2=Rkaizyo(l1)/Rkaizyo(2*l1+1)&
                /sqrt(4.d0*pi)*(4.d0*alpha2)**(dble(l1)+1.5d0)
           do ir=1,irmax_g
              r=dr_g*dble(ir-1)
              gdens(ir)=(const*exp(-alpha*r*r) &
                   -const2*exp(-alpha2*r*r))*r**l1 
           end do
           call Rhartree(gdens,pot,l1,dr_g,irmax_g)
           do iz=Nsz,Nez
           z=zL(iz)-Rion(3,a)
           do iy=Nsy,Ney
           y=yL(iy)-Rion(2,a)
           do ix=Nsx,Nex
           x=xL(ix)-Rion(1,a)
           r=sqrt(x*x+y*y+z*z)
           ir=1+aint(r/dr_g)
           if(ir<irmax_g)then
              rate=r/dr_g-dble(aint(r/dr_g))
              theta=acos(z/(r+1d-30))
              if(y>=0d0)then
                 phi=acos(x/(sqrt(x*x+y*y)+1.d-20))
              else
                 phi=2.d0*pi-acos(x/(sqrt(x*x+y*y)+1d-30))
              end if
              V_hat(ix,iy,iz,L,a)=(pot(ir)*(1.d0-rate) &
                   &+rate*pot(ir+1))*Yfunc(theta,phi,l1,m1)
           else
              V_hat(ix,iy,iz,L,a)=.0d0
           end if
           end do
           end do
           end do
        end do
     end do
! set Dens_core
     nP_core=.0d0
     do a=1,NI
        na=Kion(a)
        do iz=Nsz,Nez
        z=zL(iz)-Rion(3,a)
        do iy=Nsy,Ney
        y=yL(iy)-Rion(2,a)
        do ix=Nsx,Nex
        x=xL(ix)-Rion(1,a)
        r=sqrt(x*x+y*y+z*z)
        ir=aint(log(r/Rp_atom(na)+1d0)/H_atom(na))
        if(ir<NL_R(na)-1)then
          rate1=(r-rL_atom(ir,na))/(rL_atom(ir+1,na)-rL_atom(ir,na))
          rate2=(rL_atom(ir+1,na)-r)/(rL_atom(ir+1,na)-rL_atom(ir,na))
           nP_core(ix,iy,iz)=nP_core(ix,iy,iz)+(ncP(ir,na)*rate2 &
             +rate1*ncP(ir+1,na))
        end if
        end do
        end do
        end do
     end do
! set localized potential vbar
     Vbar=.0d0
     do a=1,NI
        na=Kion(a)
        do iz=Nsz,Nez
        z=zL(iz)-Rion(3,a)
        do iy=Nsy,Ney
        y=yL(iy)-Rion(2,a)
        do ix=Nsx,Nex
        x=xL(ix)-Rion(1,a)
        r=sqrt(x*x+y*y+z*z)
        ir=aint(log(r/Rp_atom(na)+1d0)/H_atom(na))
        if(ir<iRfilt(na))then
          rate1=(r-rL_atom(ir,na))/(rL_atom(ir+1,na)-rL_atom(ir,na))
          rate2=(rL_atom(ir+1,na)-r)/(rL_atom(ir+1,na)-rL_atom(ir,na))
           f1=(Vbar_a(ir,na)*rate2 &
                &+rate1*Vbar_a(ir+1,na))
           Vbar(ix,iy,iz)=Vbar(ix,iy,iz)+f1
        end if

        end do
        end do
        end do
     end do
     call do_loop_finit
    return
  end subroutine prep_paw
!-------10--------20--------30--------40--------50--------60----------72
  subroutine gauss_LUgendre(x,omega)
    implicit none
    real(8) x(20),omega(20)
    x(1)=-0.993128599185094924786122388471d0 ; x(20)=-x(1)
    x(2)=-0.963971927277913791267666131197d0 ; x(19)=-x(2)
    x(3)=-0.912234428251325905867752441203d0 ; x(18)=-x(3)
    x(4)=-0.839116971822218823394529061701d0 ; x(17)=-x(4)
    x(5)=-0.746331906460150792614305070355d0 ; x(16)=-x(5)
    x(6)=-0.636053680726515025452836696226d0 ; x(15)=-x(6)
    x(7)=-0.510867001950827098004364050955d0 ; x(14)=-x(7)
    x(8)=-0.373706088715419560672548177024d0 ; x(13)=-x(8)
    x(9)=-0.227785851141645078080496195368d0 ; x(12)=-x(9)
    x(10)=-0.076526521133497333754640409398d0; x(11)=-x(10)
    omega(1)=0.017614007139152118311861962351d0  ; omega(20)=omega(1)
    omega(2)=0.040601429800386941331039952274d0  ; omega(19)=omega(2)
    omega(3)=0.062672048334109063569506535187d0  ; omega(18)=omega(3)
    omega(4)=0.083276741576704748724758143222d0  ; omega(17)=omega(4)
    omega(5)=0.101930119817240435036750135480d0  ; omega(16)=omega(5)
    omega(6)=0.118194531961518417312377377711d0  ; omega(15)=omega(6)
    omega(7)=0.131688638449176626898494499748d0  ; omega(14)=omega(7)
    omega(8)=0.142096109318382051329298325067d0  ; omega(13)=omega(8)
    omega(9)=0.149172986472603746787828737001d0  ; omega(12)=omega(9)
    omega(10)=0.152753387130725850698084331955d0 ; omega(11)=omega(10)
    return
  end subroutine gauss_LUgendre
!-------10--------20--------30--------40--------50--------60----------72
  subroutine set_yfunc_d
    use global_variables
    implicit none
    integer iz,iy,L,l1,m1
    real(8) theta,phi,dp,Yfunc
    
    dp=2.d0*pi/20.d0
    do L=1,9
      l1=l_azim(L)
      m1=m_magn(L)
      do iy=1,20
        phi=dp*dble(iy)
        do iz=1,20
          theta=acos(GLz(iz)+1d-30)
          GLYlm(iy,iz,L)=Yfunc(theta,phi,l1,m1)
        end do
      end do
    end do
  end subroutine set_yfunc_d
!-------10--------20--------30--------40--------50--------60----------72
  subroutine calc_Vc
    use global_variables
    implicit none
    complex(8) fcomp1,fcomp2,fcomp3,fcomp4
    real(8) f1,f2,f3,f4,a1,a2,a3,a4,x,y,z,r,theta,phi,Rab
    complex(8)Y_Y_Gauss
    integer nuc1,nuc2,no1,no2,l1,l2,m1,m2,max,min,na1,na2
    
    Vc=.0d0
    call inirac 
    do nuc1=1,NI
      do nuc2=1,NI
        if(nuc1==nuc2)cycle
        na1=Kion(nuc1)
        na2=Kion(nuc1)
        do no1=1,9
          do no2=1,9
            l1=l_azim(no1)
            m1=m_magn(no1)
            l2=l_azim(no2)
            m2=m_magn(no2)
            if(m1>=m2)then
              max=m1
              min=m2
            else
              max=m2
              min=m1
            end if
            x=Rion(1,nuc1)-Rion(1,nuc2)
            y=Rion(2,nuc1)-Rion(2,nuc2)
            z=Rion(3,nuc1)-Rion(3,nuc2)
            Rab=sqrt(x*x+y*y+z*z)
            theta=acos(z/(Rab+1d-30))
            if(y>=(.0d0))then
              phi=acos(x/(sqrt(x*x+y*y)+1d-30))
            else
              phi=2d0*pi-acos(x/(sqrt(x*x+y*y)+1d-30))
            end if
            a1=9.d0/(rcut(na1)**2)
            a2=9.d0/(rcut(na2)**2)
            a3=9.d0/softning(na1)/(Rcut(na1)**2)
            a4=9.d0/softning(na2)/(Rcut(na2)**2)
            if((max>0).and.(min>0))then
              fcomp1=(-1.d0)**max* &
                &Y_Y_gauss(l1,l2,-max,min,a1,a2,Rab,theta,phi)
              fcomp2=Y_Y_gauss(l1,l2,max,min,a1,a2,Rab,theta,phi)
              f1=(-1.d0)**(m1+m2)*real(fcomp1+fcomp2)
              fcomp3=(-1.d0)**max* &
                &Y_Y_gauss(l1,l2,-max,min,a3,a4,Rab,theta,phi)
              fcomp4=Y_Y_gauss(l1,l2,max,min,a3,a4,Rab,theta,phi)
              f2=(-1.d0)**(m1+m2)*real(fcomp3+fcomp4)
              Vc(no1,no2,nuc1,nuc2)=f1-f2
            else if((max>0).and.(min==0))then
              fcomp1=Y_Y_gauss(l1,l2,max,min,a1,a2,Rab,theta,phi)
              f1=dsqrt(2.d0)*(1.d0)**max*real(fcomp1)
              fcomp3=Y_Y_gauss(l1,l2,max,min,a3,a4,Rab,theta,phi)
              f2=dsqrt(2.d0)*(1.d0)**max*real(fcomp3)
              Vc(no1,no2,nuc1,nuc2)=f1-f2
            else if((max>0).and.(min<0))then
              fcomp1=(-1.d0)**max* &
                &Y_Y_gauss(l1,l2,-max,-min,a1,a2,Rab,theta,phi)
              fcomp2=Y_Y_gauss(l1,l2,max,-min,a1,a2,Rab,theta,phi)
              f1=(-1.d0)**(max-min)*aimag(fcomp1+fcomp2)
              fcomp3=(-1.d0)**max* &
                &Y_Y_gauss(l1,l2,-max,-min,a3,a4,Rab,theta,phi)
              fcomp4=Y_Y_gauss(l1,l2,max,-min,a3,a4,Rab,theta,phi)
              f2=(-1.d0)**(max-min)*aimag(fcomp3+fcomp4)
              Vc(no1,no2,nuc1,nuc2)=f1-f2
            else if((max==0).and.(min==0))then
              f1=real(Y_Y_gauss(l1,l2,max,min,a1,a2,Rab,theta,phi))
              f2=real(Y_Y_gauss(l1,l2,max,min,a3,a4,Rab,theta,phi))
              Vc(no1,no2,nuc1,nuc2)=f1-f2
            else if((max==0).and.(min<0))then
              fcomp1=Y_Y_gauss(l1,l2,max,min,a1,a2,Rab,theta,phi)
              f1=-sqrt(2d0)*aimag(fcomp1)
              fcomp2=Y_Y_gauss(l1,l2,max,min,a3,a4,Rab,theta,phi)
              f2=-sqrt(2.d0)*aimag(fcomp2)
              Vc(no1,no2,nuc1,nuc2)=f1-f2
            else if((max<0).and.(min<0))then
              fcomp1=(-1.d0)**abs(max)* &
                &Y_Y_gauss(l1,l2,max,-min,a1,a2,Rab,theta,phi)
              fcomp2=Y_Y_gauss(l1,l2,abs(max), &
                abs(min),a1,a2,Rab,theta,phi)
              f1=(-1.d0)*(-1.d0)**abs(m1+m2)*real(fcomp1-fcomp2)
              fcomp3=(-1.d0)**abs(max)* &
                &Y_Y_gauss(l1,l2,max,-min,a3,a4,Rab,theta,phi)
              fcomp4=Y_Y_gauss(l1,l2,abs(max), &
                abs(min),a3,a4,Rab,theta,phi)
              f2=(-1.d0)*(-1.d0)**abs(m1+m2)*real(fcomp3-fcomp4)
              Vc(no1,no2,nuc1,nuc2)=f1-f2
            else
              write(*,*)'Error matrix element Vc',max,min
              stop
            end if
          end do
        end do
      end do
    end do
    do nuc1=1,NI
      na1=Kion(nuc1)
      a1=9.d0/(Rcut(na1)**2)
      a2=9.d0/softning(na1)/(Rcut(na1)**2)
      do no1=1,9
        do no2=1,9
          l1=L_azim(no1)
          m1=M_magn(no1)
          l2=L_azim(no2)
          m2=M_magn(no2)
          if(no1.eq.no2)then
            if(l1.eq.0)then
              ! l=0
              f1=(4.d0*a1)**1.5d0/sqrt(4.d0*pi)
              f1=f1*f1*2.d0*pi**1.5d0*(.5d0/a1)**2.5d0
              f2=(4.d0*a2)**1.5d0/sqrt(4.d0*pi)
              f2=f2*f2*2.d0*pi**1.5d0*(.5d0/a2)**2.5d0
              Vc(no1,no2,nuc1,nuc1)=f1-f2
            else if(l1.eq.1)then
              ! l=1
              f1=(4.d0*a1)**2.5d0/dsqrt(4.d0*pi)/6.d0
              f1=f1*f1*pi**1.5d0*(.5d0/a1)**3.5d0
              f2=(4.d0*a2)**2.5d0/dsqrt(4.d0*pi)/6.d0
              f2=f2*f2*pi**1.5d0*(.5d0/a2)**3.5d0
              Vc(no1,no2,nuc1,nuc1)=f1-f2      
            else if(l1.eq.2)then
              ! l=2
              f1=(4.d0*a1)**3.5d0/dsqrt(4.d0*pi)/60.d0
              f1=f1*f1*pi**1.5d0*1.5d0*(.5d0/a1)**4.5d0
              f2=(4.d0*a2)**3.5d0/dsqrt(4.d0*pi)/60.d0
              f2=f2*f2*pi**1.5d0*1.5d0*(.5d0/a2)**4.5d0
              Vc(no1,no2,nuc1,nuc1)=f1-f2
            else
              write(*,*)'Error matrix element Vc'
            end if
          else
            Vc(no1,no2,nuc1,nuc1)=.0d0
          end if
        end do
      end do
    end do
    return
  end subroutine calc_Vc
!-------10--------20--------30--------40--------50--------60----------72
  subroutine set_angular_momentum_table
    use global_variables
      
    allocate(l_azim(9),m_magn(9))
    l_azim(1)=0;m_magn(1)=0
    l_azim(2)=1;m_magn(2)=1
    l_azim(3)=1;m_magn(3)=0
    l_azim(4)=1;m_magn(4)=-1
    l_azim(5)=2;m_magn(5)=2
    l_azim(6)=2;m_magn(6)=1
    l_azim(7)=2;m_magn(7)=0
    l_azim(8)=2;m_magn(8)=-1
    l_azim(9)=2;m_magn(9)=-2
    return
  end subroutine set_angular_momentum_table
!-------10--------20--------30--------40--------50--------60----------72
  subroutine Rhartree(gdens,pot,l,dr,irmax_g)
!    use global_variables
    implicit none
    integer l,ir,irmax_g
    real(8),parameter :: pi=3.1415926535897932d0    
    real(8) Rmax,q,g
    real(8) r,coeff1,coeff2
    real(8) f1,f2,f3,const,dr
    real(8) gdens(irmax_g),pot(irmax_g)
    real(8),allocatable :: chi(:)
! l:angular momentim number
! q:multipole moment
! g:gaunt factor
    allocate (chi(irmax_g))
    coeff1=-4.d0*pi
    coeff2=dble(l*(l+1))
    Rmax=dr*dble(irmax_g-1)
      
    chi(1)=.0d0
    chi(2)=dr
    do ir=2,irmax_g-1
       r=dr*dble(ir-1)
       f1=coeff2*chi(ir)/r/r+coeff1*r*gdens(ir)
       chi(ir+1)=2.d0*chi(ir)-chi(ir-1)+f1*dr*dr
    end do
    q=.0d0
    do ir=1,irmax_g
       r=dr*dble(ir-1)
       q=q+gdens(ir)*r**(l+2)
    end do
    q=q*dr
    const=q/(Rmax**(l+1))*4.d0*pi/dble(2*l+1)-chi(irmax_g)/Rmax
    const=const/(Rmax**l)
    do ir=2,irmax_g
       r=dr*dble(ir-1)
       pot(ir)=chi(ir)/r+const*r**l
    end do
    
    pot(1)=2.d0*pot(2)-pot(3)
    return
    
  end subroutine Rhartree
!-------10--------20--------30--------40--------50--------60----------72
  subroutine do_loop_finit
    use global_variables
    implicit none
    real(8) x
    integer a,na

! calc cost cut      
    do a=1,NI
      na=Kion(a)
! a x
      x=Rion(1,a)-Rfilt(na)
      iproj(1,1,a)=aint((x+0.5d0*length_x)/H)
      x=Rion(1,a)+Rfilt(na)
      iproj(2,1,a)=1+aint((x+0.5d0*length_x)/H)
! a y
      x=Rion(2,a)-Rfilt(na)
      iproj(1,2,a)=aint((x+0.5d0*length_y)/H)
      x=Rion(2,a)+Rfilt(na)
      iproj(2,2,a)=1+aint((x+0.5d0*length_y)/H)
! a z
      x=Rion(3,a)-Rfilt(na)
      iproj(1,3,a)=aint((x+0.5d0*length_z)/H)
      x=Rion(3,a)+Rfilt(na)
      iproj(2,3,a)=1+aint((x+0.5d0*length_z)/H)
    end do
    return
  end subroutine do_loop_finit
!-------10--------20--------30--------40--------50--------60----------72
  subroutine init_wf
    use global_variables
    implicit none
    integer ix,iy,iz,p,iseed,ierr
    real(8) x1,y1,z1,r2,rnd,Rx_max,Ry_max,Rz_max
    real(8) x2,y2,z2,r22
    integer a,na,ir,l1,m1
    real(8) phi,theta,rate1,rate2,r
    real(8) Yfunc
    real(8) s0,s
    
    
    Rx_max=maxval(abs(Rion(1,:)))+4d0
    Ry_max=maxval(abs(Rion(2,:)))+4d0
    Rz_max=maxval(abs(Rion(3,:)))+4d0
    
    iseed=123
    do p=1,NST
       call quickrnd(iseed,rnd)
       x1=Rx_max*(2*rnd-1)
       call quickrnd(iseed,rnd)
       y1=Ry_max*(2*rnd-1)
       call quickrnd(iseed,rnd)
       z1=Rz_max*(2*rnd-1)
       call quickrnd(iseed,rnd)
       x2=Rx_max*(2*rnd-1)
       call quickrnd(iseed,rnd)
       y2=Ry_max*(2*rnd-1)
       call quickrnd(iseed,rnd)
       z2=Rz_max*(2*rnd-1)
       do iz=Nsz,Nez
          do iy=Nsy,Ney
             do ix=Nsx,Nex
                r2=(xL(ix)-x1)**2+(yL(iy)-y1)**2+(zL(iz)-z1)**2
                r22=(xL(ix)-x2)**2+(yL(iy)-y2)**2+(zL(iz)-z2)**2
                psi(ix,iy,iz,p)=exp(-0.5d0*r2)+exp(-0.25*r22)
             enddo
          enddo
       enddo
    end do

    return
  end subroutine init_wf
!-------10--------20--------30--------40--------50--------60----------72
  subroutine quickrnd(iseed,rnd)
    implicit none
    integer,parameter :: im=6075,ia=106,ic=1283
    integer iseed
    real(8) rnd
    iseed=mod(iseed*ia+ic,im)
    rnd=dble(iseed)/dble(im)
    return
  end subroutine quickrnd
!-------10--------20--------30--------40--------50--------60----------72
  subroutine Gram_Schmidt
    use global_variables
    implicit none
    integer p,q,ierr
    real(8) s0,s
    real(8), allocatable :: tov0(:),tov(:)

    allocate(tov0(NST),tov(NST))

    do p=1,NST
      tpsi(:,:,:)=psi(:,:,:,p)
      call Spsi
      do q=1,p-1
        tov0(q)=sum(Stpsi(:,:,:)*psi(:,:,:,q))*H**3
      end do
      call MPI_ALLREDUCE(tov0(1),tov(1),p-1, &
        MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      
      do q=1,p-1
        psi(:,:,:,p)=psi(:,:,:,p)-tov(q)*psi(:,:,:,q)
      end do
      
      tpsi(:,:,:)=psi(:,:,:,p)
      call Spsi
      s0=sum(psi(:,:,:,p)*Stpsi(:,:,:))*H**3
      call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      psi(:,:,:,p)=psi(:,:,:,p)/sqrt(s)
    end do
    return
  end subroutine Gram_Schmidt
!-------10--------20--------30--------40--------50--------60----------72
  subroutine diag
    use global_variables
    implicit none
    integer p,q,ierr
    real(8),allocatable :: psi1(:,:,:,:)
    real(8),allocatable :: a(:,:),e(:),t(:,:),aik(:)
    real(8),allocatable :: a0(:,:)
    
    allocate(psi1(Nsx:Nex,Nsy:Ney,Nsz:Nez,NST))
    allocate(a(NST,NST),e(NST),t(NST,NST),aik(NST))
    allocate(a0(NST,NST))
    
    do p=1,NST
       tpsi(:,:,:)=psi(:,:,:,p)
       call hpsi
       do q=1,p
          a0(p,q)=sum(psi(:,:,:,q)*htpsi(:,:,:))*H**3
          a0(q,p)=a0(p,q)
       end do
    end do
    
    call MPI_ALLREDUCE(a0,a,NST**2,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    call JACOBI(NST,a,e,t,aik)
    psi1=0.d0
    do p=1,NST
       do q=1,NST
          psi1(:,:,:,p)=psi1(:,:,:,p)+psi(:,:,:,q)*t(q,p)
       enddo
    enddo

    esp=e
    if(myrank==0)then
       do p=1,NST
          write(*,*)'p=',p,'esp=',esp(p)*2d0*Ry,'eV'
       end do
    end if
    deallocate(a,e,t,aik,a0,psi1)         
    return
  end subroutine diag
!-------10--------20--------30--------40--------50--------60----------72
!  N: DIMENSION OF MATRIX
!  A: MATRIX TO BE DIAGONALIZED
!  E: EIGENVALUE
!  T: EIGENVECTOR C(J,I) FOR EIGENVALUE E(I)
!  AIK: WORKING AREA
  
  SUBROUTINE JACOBI(N,A,E,T,AIK)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),E(N),T(N,N),AIK(N)
      ITMAX=10**4
      EPS1=1D-10
      EPS2=1D-10
      EPS3=1D-8
      DO 2 I=1,N
      DO 2 J=1,N
    2 T(I,J)=0D0
      NM1=N-1
      SIGMA1=0D0
      OFFDSQ=0D0
      DO 5 I=1,N
      SIGMA1=SIGMA1+A(I,I)**2
      T(I,I)=1D0
      IP1=I+1
      IF(I.GE.N) GO TO 6
      DO 5 J=IP1,N
    5 OFFDSQ=OFFDSQ+A(I,J)**2
    6 S=2D0*OFFDSQ+SIGMA1
      DO 26 ITER=1,ITMAX
      DO 20 I=1,NM1
      IP1=I+1
      DO 20 J=IP1,N
      Q=ABS(A(I,I)-A(J,J))
      IF(Q.LE.EPS1) GO TO 9
      IF(ABS(A(I,J)).LE.EPS2) GO TO 20
      P=2D0*A(I,J)*Q/(A(I,I)-A(J,J))
      SPQ=SQRT(P*P+Q*Q)
      CSA=SQRT((1D0+Q/SPQ)/2D0)
      SNA=P/(2D0*CSA*SPQ)
      GO TO 10 
    9 CSA=1D0/SQRT(2D0)
      SNA=CSA
   10 CONTINUE
      DO 11 K=1,N
      HOLDKI=T(K,I)
      T(K,I)=HOLDKI*CSA+T(K,J)*SNA
   11 T(K,J)=HOLDKI*SNA-T(K,J)*CSA
      DO 16 K=1,N
      IF(K.GT.J) GO TO 15
      AIK(K)=A(I,K)
      A(I,K)=CSA*AIK(K)+SNA*A(K,J)
      IF(K.NE.J) GO TO 14
      A(J,K)=SNA*AIK(K)-CSA*A(J,K)
   14 GO TO 16
   15 HOLDIK=A(I,K)
      A(I,K)=CSA*HOLDIK+SNA*A(J,K)
      A(J,K)=SNA*HOLDIK-CSA*A(J,K)
   16 CONTINUE
      AIK(J)=SNA*AIK(I)-CSA*AIK(J)
      DO 19 K=1,J
      IF(K.LE.I) GO TO 18
      A(K,J)=SNA*AIK(K)-CSA*A(K,J)
      GO TO 19
   18 HOLDKI=A(K,I)
      A(K,I)=CSA*HOLDKI+SNA*A(K,J)
      A(K,J)=SNA*HOLDKI-CSA*A(K,J)
   19 CONTINUE
   20 A(I,J)=0D0
      SIGMA2=0D0
      DO 21 I=1,N
      E(I)=A(I,I)
   21 SIGMA2=SIGMA2+E(I)**2
      IF(ABS(1D0-SIGMA1/SIGMA2).GE.EPS3) GO TO 25
      GO TO 100
   25 SIGMA1=SIGMA2
   26 CONTINUE
      WRITE(*,*) 'ITER EXCEEDS ITMAX, NO CONVERGENCE IN JACOBI'
 100  itemp=1
!      write(*,*) 'iter=',iter
 110  is=1
      es=e(1)
      do 120 j=2,n
      if(e(j).lt.e(is)) then
      is=j
      es=e(j)
      end if
 120  continue
      aik(itemp)=e(is)
      do 130 j=1,n
      a(j,itemp)=t(j,is)
 130  continue
      itemp=itemp+1
      e(is)=1D+10
      if(itemp.gt.n) go to 140
      go to 110
 140  do 150 i=1,n
      e(i)=aik(i)
      do 160 j=1,n
      t(j,i)=a(j,i)
 160  continue
 150  continue
      RETURN
  END 
!-------10--------20--------30--------40--------50--------60----------72
  subroutine Spsi
    use global_variables
    implicit none
    integer ix,iy,iz,a,i,ierr,na
    integer isx,iex,isy,iey,isz,iez
    real(8) s0
    real(8),allocatable :: Pc_sub(:,:),Pc_sub0(:,:)
    allocate(Pc_sub(maxproj,NI),Pc_sub0(maxproj,NI))
    
    Pc_sub0=0d0
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=0d0
          do iz=isz,iez
             do iy=isy,iey
                do ix=isx,iex
                   s0=s0+tpsi(ix,iy,iz)*proj_P(ix,iy,iz,i,a)
                end do
             end do
          end do
          Pc_sub0(i,a)=s0*H**3
       end do
    end do
    call MPI_ALLREDUCE(Pc_sub0,Pc_sub,maxproj*NI,MPI_REAL8&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)


    Stpsi=tpsi

    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=sq4p*sum(Delta_c(i,:,1,na)*Pc_sub(:,a))
          do iz=isz,iez
          do iy=isy,iey
          do ix=isx,iex
             Stpsi(ix,iy,iz)=Stpsi(ix,iy,iz)+s0*proj_P(ix,iy,iz,i,a)
          end do
          end do
          end do
       end do
    end do
    return
  end subroutine Spsi
!-------10--------20--------30--------40--------50--------60----------72
  subroutine hpsi
    use global_variables
    implicit none
    integer ix,iy,iz,a,i,ierr,na
    integer isx,iex,isy,iey,isz,iez
    real(8) s0
    real(8),allocatable :: Pc_sub(:,:),Pc_sub0(:,:)
    allocate(Pc_sub(maxproj,NI),Pc_sub0(maxproj,NI))

    Pc_sub0=0d0
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=0d0
          do iz=isz,iez
             do iy=isy,iey
                do ix=isx,iex
                   s0=s0+tpsi(ix,iy,iz)*proj_P(ix,iy,iz,i,a)
                end do
             end do
          end do
          Pc_sub0(i,a)=s0*H**3
       end do
    end do
    call MPI_ALLREDUCE(Pc_sub0,Pc_sub,maxproj*NI,MPI_REAL8&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    wk(Nsx:Nex,Nsy:Ney,Nsz:Nez)=tpsi(:,:,:)
    call laplacian
    htpsi=-0.5d0*Lwk+Veff*tpsi

    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       do i=1,Nproj(na)
          s0=sum(Hc(i,:,a)*Pc_sub(:,a))
          do iz=isz,iez
          do iy=isy,iey
          do ix=isx,iex
             htpsi(ix,iy,iz)=htpsi(ix,iy,iz)+s0*proj_P(ix,iy,iz,i,a)
          end do
          end do
          end do
       end do
    end do
    return
  end subroutine hpsi
!-------10--------20--------30--------40--------50--------60----------72
  subroutine zSpsi
    use global_variables
    implicit none
    integer ix,iy,iz,a,i,ierr,na
    integer isx,iex,isy,iey,isz,iez
    complex(8) s0
    complex(8),allocatable :: zPc_sub(:,:),zPc_sub0(:,:)
    allocate(zPc_sub(maxproj,NI),zPc_sub0(maxproj,NI))
    
    zPc_sub0=0d0
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=0d0
          do iz=isz,iez
             do iy=isy,iey
                do ix=isx,iex
                   s0=s0+ztpsi(ix,iy,iz)*proj_P(ix,iy,iz,i,a)
                end do
             end do
          end do
          zPc_sub0(i,a)=s0*H**3
       end do
    end do
    call MPI_ALLREDUCE(zPc_sub0,zPc_sub,maxproj*NI,MPI_DOUBLE_COMPLEX&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    zStpsi=ztpsi

    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=sq4p*sum(Delta_c(i,:,1,na)*zPc_sub(:,a))
          do iz=isz,iez
          do iy=isy,iey
          do ix=isx,iex
             zStpsi(ix,iy,iz)=zStpsi(ix,iy,iz)+s0*proj_P(ix,iy,iz,i,a)
          end do
          end do
          end do
       end do
    end do
    
    return
  end subroutine zSpsi
!-------10--------20--------30--------40--------50--------60----------72
  subroutine psi_rho
    use global_variables
    implicit none
    integer ix,iy,iz,a,na,i,j,p,L,ierr
    integer isx,iex,isy,iey,isz,iez
    real(8) s0
    real(8),allocatable :: Pc0(:,:,:)
    allocate(Pc0(maxproj,NST,NI))
    Pc0=0d0
! calc <psi|projctor>
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          do p=1,NST
             s0=0d0
             do iz=isz,iez
             do iy=isy,iey
             do ix=isx,iex
                s0=s0+psi(ix,iy,iz,p)*proj_P(ix,iy,iz,i,a)
             end do
             end do
             end do
             Pc0(i,p,a)=s0*H**3
          end do
       end do
    end do
    call MPI_ALLREDUCE(Pc0,Pc,maxproj*NST*NI,MPI_REAL8&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
! calc D_c
    do a=1,NI
       na=Kion(a)
       do i=1,Nproj(na)
          do j=1,Nproj(na)
             s0=.0d0
             do p=1,NST
                s0=s0+occ(p)*Pc(i,p,a)*Pc(j,p,a)
             end do
             Dc(i,j,a)=s0
          end do
       end do
    end do
! calc density
    nP=0d0
    do p=1,NST
       nP(:,:,:)=nP(:,:,:)+occ(p)*psi(:,:,:,p)**2
    end do
    nP=nP+nP_core
! initial charge density
    do a=1,NI
       na=Kion(a)
       do L=1,9
          if(L.eq.1)then
             s0=Delta_a(na)
          else
             s0=.0d0
          end if
          do i=1,Nproj(na)
             do j=1,Nproj(na)
                s0=s0+Delta_c(i,j,L,na)*Dc(i,j,a)
             end do
          end do
          Qc(L,a)=s0
       end do
    end do
    rho=nP
    do a=1,NI
       do L=1,9
          rho(:,:,:)=rho(:,:,:)+Qc(L,a)*g_hat(:,:,:,L,a)
       end do
    end do
    return
  end subroutine psi_rho
!-------10--------20--------30--------40--------50--------60----------72
  subroutine Hartree         ! rho -> Vh
    use global_variables
    implicit none
    integer iter,ierr,lm,l,m,icount,ix,iy,iz,i
    real(8) rhoLM0(25),rhoLM(25)
    real(8) s10,s1,s20,s2,s30,s3,ak,ck
    real(8) ,allocatable :: zk(:,:,:)
    real(8) ,allocatable :: pk(:,:,:)
    real(8) ,allocatable :: tk(:,:,:)
    
    allocate(zk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(pk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(tk(Nsx:Nex,Nsy:Ney,Nsz:Nez))

    wk=0.d0
    if(Bcondition=='yes') then
       rhoLM0=0
       lm=0
       do l=0,4
          do m=-l,l
             lm=lm+1
             do iz=Nsz,Nez
             do iy=Nsy,Ney
             do ix=Nsx,Nex
                rhoLM0(lm)=rhoLM0(lm)+rho(ix,iy,iz)*Ylm(ix,iy,iz,l,m)
             end do
             end do
             end do
          end do
       end do
       rhoLM0=rhoLM0*H**3
       call MPI_ALLREDUCE(rhoLM0,rhoLM,25,MPI_REAL8,MPI_SUM,&
            MPI_COMM_WORLD,ierr)

       lm=0
       do l=0,4
          do m=-l,l
             lm=lm+1
             do iz=Nsz-Nd,Nez+Nd
                do iy=Nsy-Nd,Ney+Nd
                   do ix=Nsx-Nd,Nex+Nd
                      icount=0
                      if(ix<1.or.ix>NLx) icount=icount+1
                      if(iy<1.or.iy>NLy) icount=icount+1
                      if(iz<1.or.iz>NLz) icount=icount+1
                      if(icount /= 1) cycle
                      wk(ix,iy,iz)=wk(ix,iy,iz)&
                           +rhoLM(lm)*Ylm(ix,iy,iz,l,m)&
                           /(xL(ix)**2+yL(iy)**2+zL(iz)**2)**(l+0.5d0)
                   end do
                end do
             end do
          end do
       end do
    end if
    
    zk=-4*Pi*rho

    do i=1,2
       if(i == 2)Vh=0d0
       wk(Nsx:Nex,Nsy:Ney,Nsz:Nez)=Vh(:,:,:)
       call Laplacian
       zk=zk-Lwk
       pk=zk
       s10=sum(zk**2)*H**3
       call MPI_ALLREDUCE(s10,s1,1,MPI_REAL8,MPI_SUM &
            ,MPI_COMM_WORLD,ierr)
       
       do iter=1,1000
          wk=0.d0
          wk(Nsx:Nex,Nsy:Ney,Nsz:Nez)=pk(:,:,:)
          call Laplacian
          tk=Lwk
          s20=sum(zk*tk)*H**3
          call MPI_ALLREDUCE(s20,s2,1,MPI_REAL8,MPI_SUM, &
               MPI_COMM_WORLD,ierr)
          ak=s1/s2
          Vh=Vh+ak*pk
          zk=zk-ak*tk
          s30=sum(zk**2)*H**3
          call MPI_ALLREDUCE(s30,s3,1,MPI_REAL8,MPI_SUM, &
               MPI_COMM_WORLD,ierr)
          if(abs(s3)<1d-11) then
             iterVh=iter
             deallocate(zk,pk,tk)
             if(myrank==0)then
                write(*,*)'monopole'
                write(*,*)rhoLM(1)
                write(*,*)'dipole'
                write(*,*)-rhoLM(2),-rhoLM(4),rhoLM(3)
                write(*,*)'end hartree=',iterVh
             end if
             return
          endif
          ck=s3/s1
          s1=s3
          pk=zk+ck*pk
       end do
    end do

    if(Myrank == 0) then
       write(*,*) 'Warning:Vh iteration not converged, s3=',s3
    endif
    call MPI_FINALIZE(ierr)
    stop
    iterVh=iter
    return
    
    deallocate(zk)
    deallocate(pk)
    deallocate(tk)
      
  end subroutine Hartree
!-------10--------20--------30--------40--------50--------60----------72
  subroutine Exc_Cor
    use global_variables
    implicit none
    integer ix,iy,iz
    real(8) trho,rs,rssq,rsln,V_xc,E_xc
    real(8) trho_sigma,x_sigma
    
    do iz=Nsz,Nez
    do iy=Nsy,Ney
    do ix=Nsx,Nex
       trho=nP(ix,iy,iz)+1d-20
       rs=(3d0/(4*Pi*trho))**(1d0/3d0)
       V_xc=-4d0/3d0*0.4582d0/rs
       E_xc=-.4582d0/rs
       if(rs>1d0) then
          rssq=sqrt(rs)
          V_xc=V_xc+gammaU*(1d0+7d0/6d0*beta1U*rssq &
               +4d0/3d0*beta2U*rs) &
               /(1d0+beta1U*rssq+beta2U*rs)**2
          E_xc=E_xc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
       else
          rsln=log(rs)
          V_xc=V_xc+AU*rsln+(BU-AU/3d0) &
               +2d0/3d0*CU*rs*rsln+(2d0*DU-CU)/3d0*rs
          E_xc=E_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
       endif
       Vexc(ix,iy,iz)=V_xc
       Eexc(ix,iy,iz)=E_xc
    enddo
    enddo
    enddo
    return
  end subroutine Exc_Cor
!-------10--------20--------30--------40--------50--------60----------72
  subroutine local_potential
    use global_variables
    implicit none
    integer a,L
      
    call Hartree
    call Exc_Cor

    Vhat_sum=.0d0
    do a=1,NI
       do L=1,9
          Vhat_sum(:,:,:)=Vhat_sum(:,:,:)+Qc(L,a)*V_hat(:,:,:,L,a)
       end do
    end do
    
    Veff=Vh+Vexc+Vhat_sum+Vbar
    
    return
  end subroutine local_potential
!-------10--------20--------30--------40--------50--------60----------72
  subroutine nonlocal_potential
    use global_variables
    implicit none
    integer a1,a2,i,j,l1,l2,no1,no2,na,na2,L,ierr
    integer ir,itheta,iphi
    real(8) rs1,rs2,rssq,rsln
    real(8) r,s,s0
    real(8) f1,f2,f3,f4,f5,f6
    real(8),allocatable ::muA(:,:,:,:),muP(:,:,:,:)
    allocate(muA(20,20,0:irmax,NI),muP(20,20,0:irmax,NI))
    if(myrank==0)write(*,*)'make na'
    Dmu_c=0d0
    do a1=1,NI
       na=Kion(a1)
       do ir=0,iRcut(na)
       do itheta=1,20
       do iphi=1,20
          f1=ncA(ir,na)
          f4=ncP(ir,na)
          do i=1,Nproj(na)
             f3=0d0
             f6=0d0
             l1=NLproj(i,na)
             no1=NRproj(i,na)
             do j=i+1,Nproj(na)
                l2=NLproj(j,na)
                no2=NRproj(j,na)
                f2=Dc(i,j,a1) &
                     *GLYlm(iphi,itheta,l1)*GLYlm(iphi,itheta,l2)
                f3=f3+f2*phiA(ir,no1,na)*phiA(ir,no2,na)
                f6=f6+f2*phiP(ir,no1,na)*phiP(ir,no2,na)
             end do
             f1=f1+f3*2d0 &
                  +Dc(i,i,a1) &
                  *(phiA(ir,no1,na)*GLYlm(iphi,itheta,l1))**2
             f4=f4+f6*2d0 &
                  +Dc(i,i,a1) &
                  *(phiP(ir,no1,na)*GLYlm(iphi,itheta,l1))**2
          end do
          nacA(iphi,itheta,ir,a1)=f1
          nacP(iphi,itheta,ir,a1)=f4
       end do
       end do
       end do
    end do
! atomcentered exchange-correlation 1,2
    do a1=1,NI
       na=Kion(a1)
       do ir=0,ircut(na)
       do itheta=1,20
       do iphi=1,20
! exchange potentual
          rs1=nacA(iphi,itheta,ir,a1)+1d-20
          rs2=nacP(iphi,itheta,ir,a1)+1d-20
          rs1=(3.d0/(4.d0*pi*rs1))**(1d0/3d0)
          rs2=(3.d0/(4.d0*pi*rs2))**(1d0/3d0)
          f1=-4d0/3d0*0.4582d0/rs1
          f2=-4d0/3d0*0.4582d0/rs2
          if(rs1>=1d0) then
             rssq=sqrt(rs1)
             f1=f1+gammaU*(1d0+7d0/6d0*beta1U*rssq &
                  +4d0/3d0*beta2U*rs1) &
                  /(1d0+beta1U*rssq+beta2U*rs1)**2
          else
             rsln=log(rs1)
             f1=f1+AU*rsln+(BU-AU/3d0) &
               +2d0/3d0*CU*rs1*rsln+(2d0*DU-CU)/3d0*rs1
          end if
          if(rs2>=1d0) then
             rssq=sqrt(rs2)
             f2=f2+gammaU*(1d0+7d0/6d0*beta1U*rssq &
                  +4d0/3d0*beta2U*rs2) &
                  /(1d0+beta1U*rssq+beta2U*rs2)**2
          else
             rsln=log(rs2)
             f2=f2+AU*rsln+(BU-AU/3d0) &
               +2d0/3d0*CU*rs2*rsln+(2d0*DU-CU)/3d0*rs2
          end if
          muA(iphi,itheta,ir,a1)=f1
          muP(iphi,itheta,ir,a1)=f2
       end do
       end do
       end do
    end do
! calc Exc_c
    do a1=1,NI
       na=Kion(a1)
       do i=1,Nproj(na)
       do j=1,Nproj(na)
          no1=NRproj(i,na)
          no2=NRproj(j,na)
          l1=NLproj(i,na)
          l2=NLproj(j,na)
          f1=.0d0
          f2=.0d0
          do ir=0,ircut(na)
             r=rL_atom(ir,na)
             f3=.0d0
             f4=.0d0
             do itheta=1,20
                f3=f3+sum(muA(:,itheta,ir,a1)*GLYlm(:,itheta,l1) &
                     *GLYlm(:,itheta,l2))*GLome(itheta) 
                f4=f4+sum(muP(:,itheta,ir,a1)*GLYlm(:,itheta,l1) &
                     *GLYlm(:,itheta,l2))*GLome(itheta) 
             end do
             
             f3=f3*2.d0*pi/20.d0
             f4=f4*2.d0*pi/20.d0
             f1=f1+f3*phiA(ir,no1,na)*phiA(ir,no2,na)*expXr2(ir,na)
             f2=f2+f4*phiP(ir,no1,na)*phiP(ir,no2,na)*expXr2(ir,na)
          end do
          Dmu_c(i,j,a1)=(f1-f2)*H_atom(na)*Rp_atom(na)
       end do
       end do
    end do
! calc W_c
    do a1=1,NI
       na=Kion(a1)
       do L=1,9
          s0=sum(g_hat(:,:,:,L,a1)*Vh(:,:,:))*H**3
          call MPI_ALLREDUCE(s0,f1,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          s0=sum(nP(:,:,:)*v_hat(:,:,:,L,a1))*H**3
          call MPI_ALLREDUCE(s0,f2,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          f3=.0d0
          do a2=1,NI
             na2=Kion(a2)
             do L2=1,9
                f3=f3+Vc(L,L2,a1,a2)*Qc(L2,a2)
             end do
          end do
          Wc(L,a1)=f1+f2+f3
       end do
    end do
    Hc=0d0
    do a1=1,NI
       na=Kion(a1)
       do i=1,Nproj(na)
          do j=1,Nproj(na)
             f1=sum(Delta_c(i,j,:,na)*Wc(:,a1))
             f2=2.d0*sum(Cc(i,j,:,:,na)*Dc(:,:,a1))
             f1=f1+Dmu_c(i,j,a1)+Bc(i,j,na)+f2
             Hc(i,j,a1)=f1
          end do
       end do
    end do
    return
  end subroutine nonlocal_potential
!-------10--------20--------30--------40--------50--------60----------72
  subroutine DTcg
    use global_variables
    implicit none
    integer iter,p,q,ierr
    real(8) s0,s
    real(8), allocatable :: sov0(:),sov(:)
    real(8) xkHxk,xkxk,Rk,gkgk,xkpk,pkpk,pkHxk,pkHpk
    real(8) uk,alpha,Ak,Bk,Ck
    real(8) ,allocatable :: xk(:,:,:),hxk(:,:,:),Sxk(:,:,:)
    real(8) ,allocatable :: gk(:,:,:),pk(:,:,:),hpk(:,:,:),Spk(:,:,:)
      
    allocate(sov0(NST),sov(NST))
    allocate(xk(Nsx:Nex,Nsy:Ney,Nsz:Nez),hxk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(pk(Nsx:Nex,Nsy:Ney,Nsz:Nez),hpk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Sxk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(Spk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(gk(Nsx:Nex,Nsy:Ney,Nsz:Nez))


    do p=1,NST
       tpsi(:,:,:)=psi(:,:,:,p)
       xk=tpsi
       call hpsi
       call Spsi
       hxk=htpsi
       Sxk=Stpsi
       s0=sum(xk*hxk)*H**3
       call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM, MPI_COMM_WORLD,ierr)
       xkHxk=s
       xKxK=1d0
       Rk=xkHxk/xkxk
       do iter=1,Ncg
          gk=2d0*(hxk-Rk*Sxk)/xkxk
          tpsi(:,:,:)=gk(:,:,:)
          call Spsi
          do q=1,p-1
             sov0(q)=sum(psi(:,:,:,q)*Stpsi)*H**3
          end do
          call MPI_ALLREDUCE(sov0(1),sov(1),p-1,MPI_REAL8,MPI_SUM, &
               MPI_COMM_WORLD,ierr)
          do q=1,p-1
            gk=gk-sov(q)*psi(:,:,:,q)
          enddo
          s0=sum(gk**2)*H**3             
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          select case (iter)
            case(1)
              pk=-gk
            case default
              uk=s/gkgk
              pk=-gk+uk*pk
          end select
          gkgk=s
          tpsi(:,:,:)=pk(:,:,:)
          call hpsi;call Spsi
          hpk=htpsi;Spk=Stpsi
          s0=sum(Sxk*pk)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          xkpk=s
          s0=sum(Spk*pk)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          pkpk=s
          s0=sum(hxk*pk)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          pkHxk=s
          s0=sum(hpk*pk)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          pkHpk=s

          Ak=pkHpk*xkpk-pkHxk*pkpk
          Bk=pkHpk*xkxk-xkHxk*pkpk
          Ck=pkHxk*xkxk-xkHxk*xkpk
          alpha=(-Bk+sqrt(Bk*Bk-4d0*Ak*Ck))/(2d0*Ak)
          xk=xk+alpha*pk
          hxk=hxk+alpha*hpk
          Sxk=Sxk+alpha*Spk
          s0=sum(xk*Sxk)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          xkxk=s
          s0=sum(xk*hxk)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
               ,MPI_COMM_WORLD,ierr)
          xkHxk=s
          Rk=xkHxk/xkxk
       end do
       s0=sum(xk*Sxk)*H**3
       call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
            ,MPI_COMM_WORLD,ierr)
       psi(:,:,:,p)=xk/sqrt(s)
    end do
    deallocate(sov0,sov)
    deallocate(xk)
    deallocate(hxk)
    deallocate(gk)
    deallocate(pk)
    return       
  end subroutine DTcg
!-------10--------20--------30--------40--------50--------60----------72
  subroutine total_energy(type)
    use global_variables
    implicit none
    character(2) type
    real(8),allocatable ::Exc(:,:,:),Kpsi(:,:,:)
    real(8),allocatable ::excA(:,:,:,:),excP(:,:,:,:)
    real(8),allocatable ::Delta_E(:)
    real(8),allocatable :: u(:,:,:),hu(:,:,:)
    real(8) rs,rs1,rs2,f1,f2,f3,f4,f5,f6,r
    integer a1,a2,ir,ix,iy,iz,na,i1,i2,no1,ierr
    integer l1,l2,l3,m1,m2,m3
    real(8) s1,s2,s3,s4,s5

    allocate(Exc(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(kpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate(excA(20,20,0:irmax,NI))
    allocate(excP(20,20,0:irmax,NI),Delta_E(NI))
    allocate(u(Nsx:Nex,Nsy:Ney,Nsz:Nez))
! atomcentered exchange-correlation 1,2
    do a1=1,NI
       na=Kion(a1)
       do ir=0,ircut(na)
       do iz=1,20
       do iy=1,20
! exchange energy
          f1=-0.75d0*(3.d0*nacA(iy,iz,ir,a1)/pi)**(1.d0/3.d0)
          f2=-0.75d0*(3.d0*nacP(iy,iz,ir,a1)/pi)**(1.d0/3.d0)
          rs1=nacA(iy,iz,ir,a1)
          rs2=nacP(iy,iz,ir,a1)
          if((rs1.le.(.0d0)).or.(rs2.le.(.0d0)))then
             write(*,*)'Error! density is negative'
             stop
          end if
          rs1=(3.d0/(4.d0*pi*rs1+1.d-15))**(1.d0/3.d0)
          rs2=(3.d0/(4.d0*pi*rs2+1.d-15))**(1.d0/3.d0)
          if(rs1.ge.(1.d0))then
             f1=f1+gammaU/(1d0+beta1U*dsqrt(rs1)+beta2U*rs1)
          else
             f1=f1+AU*dlog(rs1)+BU+CU*rs1*dlog(rs1)+DU*rs1
          end if
          if(rs2.ge.(1.d0))then
             f2=f2+gammaU/(1d0+beta1U*dsqrt(rs2)+beta2U*rs2)
          else
             f2=f2+AU*dlog(rs2)+BU+CU*rs2*dlog(rs2)+DU*rs2
          end if
          excA(iy,iz,ir,a1)=f1
          excP(iy,iz,ir,a1)=f2
       end do
       end do
       end do
    end do
! calc Ex_c
    do a1=1,NI
       na=Kion(a1)
       f1=.0d0
       f2=.0d0
       do ir=0,ircut(na)
          r=rL_atom(ir,na)
          f3=.0d0
          f4=.0d0
          do iz=1,20
             f3=f3+sum(excA(:,iz,ir,a1)*nacA(:,iz,ir,a1))*GLome(iz) 
             f4=f4+sum(excP(:,iz,ir,a1)*nacP(:,iz,ir,a1))*GLome(iz)
          end do
          f3=f3*2.d0*pi/20.d0
          f4=f4*2.d0*pi/20.d0
          f1=f1+f3*expXr2(ir,na)
          f2=f2+f4*expXr2(ir,na)
       end do
       Exc_c(a1)=(f1-f2)*H_atom(na)*Rp_atom(na)
    end do
! calc Delta_E
    do a1=1,NI
       na=Kion(a1)
       f1=sum(Bc(:,:,na)*Dc(:,:,a1))
       f2=.0d0
       do i1=1,Nproj(na)
          do i2=1,Nproj(na)
             f2=f2+Dc(i1,i2,a1)*sum(Cc(i1,i2,:,:,na)*Dc(:,:,a1))
          end do
       end do
       Delta_E(a1)=Ac(na)+f1+f2+Exc_c(a1)
    end do

! calc Exc
    do iz=Nsz,Nez
       do iy=Nsy,Ney
          do ix=Nsx,Nex
             f1=-0.75d0*(3.d0*nP(ix,iy,iz)/pi)**(1.d0/3.d0)
             rs=(3.d0/(4.d0*pi*nP(ix,iy,iz)+1.d-15))**(1d0/3d0)
             if(rs.ge.(1.d0))then
                f1=f1+gammaU/(1d0+beta1U*dsqrt(rs)+beta2U*rs)
             else
                f1=f1+AU*dlog(rs)+BU+CU*rs*dlog(rs)+DU*rs
             end if
             Exc(ix,iy,iz)=f1
          end do
       end do
    end do
    

    s1=.0d0
    select case(type)
    case('gs')
       do no1=1,NST
          wk(Nsx:Nex,Nsy:Ney,Nsz:Nez)=psi(:,:,:,no1)
          call laplacian
          s1=s1-0.5d0*occ(no1)*sum(Lwk(:,:,:)*psi(:,:,:,no1))*H**3
       end do
    case('rt')
       do no1=1,NST
          zwk(Nsx:Nex,Nsy:Ney,Nsz:Nez)=zpsi(:,:,:,no1)
          call zlaplacian
          s1=s1-0.5d0*occ(no1) &
               *sum(zLwk(:,:,:)*conjg(zpsi(:,:,:,no1)))*H**3
       end do
    end select
    s2=0.5d0*sum(rho(:,:,:)*Vh(:,:,:))*H**3
    s3=sum(Exc(:,:,:)*nP(:,:,:))*H**3
    s4=sum(nP(:,:,:)*Vbar(:,:,:))*H**3
    s5=sum(nP(:,:,:)*Vhat_sum(:,:,:))*H**3
    f6=.0d0
    do a1=1,NI
       do a2=1,NI
          do l1=1,9
             do l2=1,9
                f6=f6+Qc(l1,a1)*Qc(l2,a2)*Vc(l1,l2,a1,a2)
             end do
          end do
       end do
    end do
    f6=f6*0.5d0
    call MPI_ALLREDUCE(s1,f1,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(s2,f2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(s3,f3,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(s4,f4,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(s5,f5,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    Etot=f1+f2+f3+f4+f5+f6+sum(Delta_E(:))
    if(myrank == 0)then
       write(*,*)'Total Energy'
       write(*,'(2x,e24.16,"eV",3x,e24.16,"a.u.")')Etot*2d0*Ry,Etot
    end if

      select case(type)
      case('rt')
      if(Myrank == 0)then
         write(12,*)iter_rt*DT*0.02419d0,Etot
      end if
      end select
  end subroutine total_energy
!-------10--------20--------30--------40--------50--------60----------72
  subroutine calc_dipole_moment
    use global_variables
    implicit none
    integer ix,iy,iz,a,p,ierr
    integer i,j,na
    real(8) s0,s
    complex(8) co1

    dipole_td=0d0
    do p=1,NST
       s0=0d0
       do iz=Nsz,Nez
       do iy=Nsy,Ney
       do ix=Nsx,Nex
       s0=s0+(ex*xL(ix)+ey*yL(iy)+ez*zL(iz))*abs(zpsi(ix,iy,iz,p))**2
       end do
       end do
       end do
       s0=s0*H**3
       call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)       
       co1=0d0
       do a=1,NI
          na=Kion(a)
          do i=1,Nproj(na)
          do j=1,Nproj(na)
             co1=co1+conjg(zPc(i,p,a))*zPc(j,p,a)*rk_c(i,j,a)
          end do
          end do
       end do
       dipole_td=dipole_td+occ(p)*(s+real(co1))
    end do
    if(Myrank == 0)then
       write(8,*)iter_rt*DT*0.02419d0,dipole_td
    end if

    return
  end subroutine calc_dipole_moment
!-------10--------20--------30--------40--------50--------60----------72
  subroutine calc_particle_number
    use global_variables
    implicit none
    integer p,ierr
    real(8) norm_rt0,norm_rt,particle_no0,particle_no

    particle_no=0d0
    do p=1,NST
       ztpsi(:,:,:)=zpsi(:,:,:,p)
       call zSpsi
       norm_rt0=sum(conjg(zpsi(:,:,:,p))*zStpsi(:,:,:))*H**3
       call MPI_ALLREDUCE(norm_rt0,norm_rt,1,MPI_REAL8 &
            ,MPI_SUM,MPI_COMM_WORLD,ierr)
       particle_no=particle_no+norm_rt*occ(p)
    end do

    if(Myrank == 0)then
       write(22,*)iter_rt*DT*0.02419d0,particle_no
    end if
    return
  end subroutine calc_particle_number
!-------10--------20--------30--------40--------50--------60----------72
  subroutine prep_rt
    use global_variables
    implicit none
    integer a,ix,iy,iz,i,j,na
! absorber    
    integer min_M
    real(8) R_abs


    allocate (zpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez,1:NST),zPc(maxproj,NST,NI))
    allocate (ztpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (zStpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (zXHtpsi(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (Vext(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (Wabs(Nsx:Nex,Nsy:Ney,Nsz:Nez))
    allocate (rk_c(maxproj,maxproj,NI))



    if(Myrank == 0)then
       open(22,file=file_elec_all)
       open(8,file=file_dipole)
       open(12,file=file_energy)
       open(14,file=file_laser)
    end if


    do a=1,NI
       na=Kion(a)
       do i=1,Nproj(na)
          do j=1,Nproj(na)
             rk_c(i,j,a)= &
              ex*(dip_c(i,j,1,na)+sq4p*Delta_c(i,j,1,na)*Rion(1,a))&
              +ey*(dip_c(i,j,2,na)+sq4p*Delta_c(i,j,1,na)*Rion(2,a))&
              +ez*(dip_c(i,j,3,na)+sq4p*Delta_c(i,j,1,na)*Rion(3,a))
          end do
       end do
    end do

! set absorber
    Wabs=0d0
    min_M=NLx
    if(min_M > NLy)min_M=NLy
    if(min_M > NLz)min_M=NLz
    length_R=min_M*dble(H)/2d0-deltaW
    if(absorber=='yes') then
       do iz=Nsz,Nez
          do iy=Nsy,Ney
             do ix=Nsx,Nex
                R_abs=sqrt(xL(ix)**2+yL(iy)**2+zL(iz)**2)
                if(R_abs > length_R)then
                   Wabs(ix,iy,iz)=-W0*(R_abs-length_R)/deltaW
                end if
             end do
          end do
       end do
    endif

! rt_option
    select case(rt_option)
    case('dipole')
       Vext=0.d0
       ft=rk_mom
    case('laser1')
       f0=5.338d-9*sqrt(Wcm2)      ! electric field in a.u.
       omega=45.569d0/wave_length  ! frequency in a.u.
       ft=0d0
       do iz=Nsz,Nez
       do iy=Nsy,Ney
       do ix=Nsx,Nex
          Vext(ix,iy,iz)=xL(ix)*ex+yL(iy)*ey+zL(iz)*ez
       end do
       end do
       end do

    case('laser2')
       f0=5.338d-9*sqrt(Wcm2)      ! electric field in a.u.
       omega=45.569d0/wave_length  ! frequency in a.u.
       ft=0d0
       do iz=Nsz,Nez
       do iy=Nsy,Ney
       do ix=Nsx,Nex
          Vext(ix,iy,iz)=xL(ix)*ex+yL(iy)*ey+zL(iz)*ez
       enddo
       enddo
       enddo

    end select



!--- wave function at t=0 ---


      select case(rt_option)
      case('dipole')
!         zpsi=psi
         call impulse 
      case('laser1')
         zpsi=psi
      case('laser2')
         zpsi=psi
      end select

      deallocate(psi)
      return
    end subroutine prep_rt
!-------10--------20--------30--------40--------50--------60----------72
    subroutine impulse
      use global_variables
      implicit none
      integer a,p,ix,iy,iz,na,i,j,no1,no2
      integer l1,l2,ir,iphi,itheta,iter,ierr
      integer Nimpulse
      complex(8),allocatable :: zPc_sub(:,:),zPc_sub0(:,:)
      complex(8),allocatable :: impulse_c(:,:,:)
      complex(8),allocatable :: xk(:,:,:),rk(:,:,:),bk(:,:,:)
      complex(8),allocatable :: qk(:,:,:),pk(:,:,:),Sxk(:,:,:)
      complex(8) alpha,Stheta,Sphi,pkqk,co0,co,co1,co2
      real(8) beta,kx,ky,kz,rkrk,s0,s,r,kr
      integer isx,iex,isy,iey,isz,iez

      allocate(zPc_sub(maxproj,NI),zPc_sub0(maxproj,NI))
      allocate(impulse_c(maxproj,maxproj,NI))
      allocate(xk(Nsx:Nex,Nsy:Ney,Nsz:Nez),rk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
      allocate(bk(Nsx:Nex,Nsy:Ney,Nsz:Nez),qk(Nsx:Nex,Nsy:Ney,Nsz:Nez))
      allocate(pk(Nsx:Nex,Nsy:Ney,Nsz:Nez),Sxk(Nsx:Nex,Nsy:Ney,Nsz:Nez))

      Nimpulse=50
! calc impulse_c
      do a=1,NI
        na=Kion(a)
        do i=1,Nproj(na)
        do j=1,Nproj(na)
          no1=NRproj(i,na)
          no2=NRproj(j,na)
          l1=NLproj(i,na)
          l2=NLproj(j,na)
          co1=.0d0
          co2=.0d0
            do ir=0,ircut(na)
              r=rL_atom(ir,na)
              Stheta=0d0
              do itheta=1,20
                Sphi=0d0
                do iphi=1,20
                  kx=ex*(sqrt(4d0*pi/3d0)*r*GLYlm(iphi,itheta,2)+Rion(1,a))
                  ky=ey*(sqrt(4d0*pi/3d0)*r*GLYlm(iphi,itheta,4)+Rion(2,a))
                  kz=ez*(sqrt(4d0*pi/3d0)*r*GLYlm(iphi,itheta,3)+Rion(3,a))
                  Sphi=Sphi+exp(zI*rk_mom*(kx+ky+kz)) &
                    *GLYlm(iphi,itheta,l1)*GLYlm(iphi,itheta,l2)
                end do
                Stheta=Stheta+Sphi*2.d0*pi/20.d0*GLome(itheta)
              end do
              co1=co1+Stheta*phiA(ir,no1,na)*phiA(ir,no2,na)*expXr2(ir,na)
              co2=co2+Stheta*phiP(ir,no1,na)*phiP(ir,no2,na)*expXr2(ir,na)
            end do
            impulse_c(i,j,a)=(co1-co2)*H_atom(na)*Rp_atom(na)
          end do
        end do
      end do
      
      do p=1,NST
        xk(:,:,:)=psi(:,:,:,p)

! calc bk -------!
        zPc_sub0=0d0
        do a=1,NI
          isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
          if(isx>iex)cycle
          isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
          if(isy>iey)cycle
          isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
          if(isz>iez)cycle
          na=Kion(a)
          do i=1,Nproj(na)
            co0=0d0
            do iz=isz,iez
              do iy=isy,iey
                do ix=isx,iex
                  co0=co0+xk(ix,iy,iz)*proj_P(ix,iy,iz,i,a)
                end do
              end do
            end do
            zPc_sub0(i,a)=co0*H**3
          end do
        end do
        call MPI_ALLREDUCE(zPc_sub0,zPc_sub,maxproj*NI, &
              MPI_double_complex,MPI_SUM,MPI_COMM_WORLD,ierr)


        do iz=Nsz,Nez
          do iy=Nsy,Ney
            do ix=Nsx,Nex
              kr=ex*xL(ix)+ey*yL(iy)+ez*zL(iz)
              bk(ix,iy,iz)=exp(zI*rk_mom*kr)*xk(ix,iy,iz)
            end do
          end do
        end do


        do a=1,NI
          isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
          if(isx>iex)cycle
          isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
          if(isy>iey)cycle
          isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
          if(isz>iez)cycle
          na=Kion(a)
          do i=1,Nproj(na)
            co0=sum(impulse_c(i,:,a)*zPc_sub(:,a))
            do iz=isz,iez
              do iy=isy,iey
                do ix=isx,iex
                  bk(ix,iy,iz)=bk(ix,iy,iz)+co0*proj_P(ix,iy,iz,i,a)
                end do
              end do
            end do
          end do
        end do
! end calc bk -------!
        ztpsi(:,:,:)=xk(:,:,:)
        call zSpsi
        rk=bk-zStpsi
        pk=rk
        s0=sum(abs(rk)**2)*H**3
        call MPI_ALLREDUCE(s0,rkrk,1,MPI_REAL8,MPI_SUM, &
          MPI_COMM_WORLD,ierr)

! CG loop -----------!
        do iter=1,Nimpulse
          ztpsi(:,:,:)=pk(:,:,:)
          call zSpsi
          qk=zStpsi
          co0=sum(conjg(pk)*qk)*H**3
          call MPI_ALLREDUCE(co0,pkqk,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            MPI_COMM_WORLD,ierr)
          
          alpha=rkrk/pkqk
          xk=xk+alpha*pk
          rk=rk-alpha*qk
          
          s0=sum(abs(rk)**2)*H**3
          call MPI_ALLREDUCE(s0,s,1,MPI_REAL8,MPI_SUM &
            ,MPI_COMM_WORLD,ierr)
          
          beta=s/rkrk
          rkrk=s
          pk=rk+beta*pk
          if(myrank==0)write(*,*)p,iter,'impulse res=',rkrk
        end do
! end CG loop -------!
        zpsi(:,:,:,p)=xk(:,:,:)
      end do
      return
    end subroutine impulse
!-------10--------20--------30--------40--------50--------60----------72
    subroutine RTI
      use global_variables
      implicit none
      real(8),allocatable :: Veff0(:,:,:),Hc0(:,:,:)
      complex(8),allocatable :: zpsi0(:,:,:,:)
      real(8) :: pulse_time

      pulse_time=pulse_time_fs/0.02418d0

      select case(rt_option)
      case('dipole')
         ft=0d0   
      case('laser1')   
         if(tt < pulse_time) then
            ft=f0*sin(Pi*tt/pulse_time)**2*sin(omega*(tt-pulse_time*0.5d0))
         else
            ft=0.d0
         endif
      case('laser2')   
         if(tt < pulse_time/2d0) then
            ft=f0*sin(Pi*tt/pulse_time)**2*sin(omega*tt)
         else
            ft=f0*sin(omega*tt)
         endif
      end select

      if(predict_corr == 'yes')then
         allocate(Veff0(Nsx:Nex,Nsy:Ney,Nsz:Nez))
         allocate(Hc0(maxproj,maxproj,NI))
         allocate(zpsi0(Nsx:Nex,Nsy:Ney,Nsz:Nez,1:NST))
         Veff0=Veff
         Hc0=Hc
         zpsi0=zpsi
         call zNDTexp
         call zpsi_rho
         call local_potential
         call nonlocal_potential
         Veff=0.5d0*(Veff+Veff0)
         Hc=0.5d0*(Hc+Hc0)
!-- set ft(dt/2) -------
         select case(rt_option)
         case('dipole')
            continue
         case('laser1')   
            if(tt+0.5d0*DT < pulse_time) then
               ft=f0*sin(Pi*(tt+0.5d0*DT)/pulse_time)**2&
                    *sin(omega*(tt+0.5d0*DT-pulse_time*0.5d0))
            else
               ft=0.d0
            endif
         case('laser2')   
            if((iter_rt+0.5d0)*DT < pulse_time/2d0) then
               ft=f0*sin(Pi*(tt+0.5d0*DT)/pulse_time)**2&
                    *sin(omega*(tt+0.5d0*DT))
            else
               ft=f0*sin(omega*(tt+0.5d0*DT))
            endif
         end select
         zpsi=zpsi0
!-- end set ft --------         
      end if
      call zNDTexp
      call zpsi_rho
      call local_potential
      call nonlocal_potential

      return
    end subroutine RTI
!-------10--------20--------30--------40--------50--------60----------72
  subroutine zpsi_rho
    use global_variables
    implicit none
    integer ix,iy,iz,a,na,i,j,p,L,ierr
    integer isx,iex,isy,iey,isz,iez
    complex(8) zs0
    real(8) s0
    complex(8),allocatable :: zPc0(:,:,:)
    allocate(zPc0(maxproj,NST,NI))
    zPc0=0d0
! calc <psi|projctorx

    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          do p=1,NST
             zs0=0d0
             do iz=isz,iez
             do iy=isy,iey
             do ix=isx,iex
                zs0=zs0+zpsi(ix,iy,iz,p)*proj_P(ix,iy,iz,i,a)
             end do
             end do
             end do
             zPc0(i,p,a)=zs0*H**3
          end do
       end do
    end do
    call MPI_ALLREDUCE(zPc0,zPc,maxproj*NST*NI,MPI_DOUBLE_COMPLEX&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
! calc D_c
    do a=1,NI
       na=Kion(a)
       do i=1,Nproj(na)
          do j=1,Nproj(na)
             zs0=.0d0
             do p=1,NST
                zs0=zs0+occ(p)*conjg(zPc(i,p,a))*zPc(j,p,a)
             end do
             Dc(i,j,a)=zs0
          end do
       end do
    end do
! calc density
    nP=0d0
    do p=1,NST
       nP(:,:,:)=nP(:,:,:)+occ(p)*abs(zpsi(:,:,:,p))**2
    end do
    nP=nP+nP_core
! initial charge density
    do a=1,NI
       na=Kion(a)
       do L=1,9
          if(L.eq.1)then
             s0=Delta_a(na)
          else
             s0=.0d0
          end if
          do i=1,Nproj(na)
             do j=1,Nproj(na)
                s0=s0+Delta_c(i,j,L,na)*Dc(i,j,a)
             end do
          end do
          Qc(L,a)=s0
       end do
    end do
    rho=nP
    do a=1,NI
       do L=1,9
          rho(:,:,:)=rho(:,:,:)+Qc(L,a)*g_hat(:,:,:,L,a)
       end do
    end do
    return
  end subroutine zpsi_rho
!-------10--------20--------30--------40--------50--------60----------72
  subroutine zNDTexp
    use global_variables
    implicit none
    integer iexp,p
    complex(8) zfac

    do p=1,NST
       zfac=1d0
       ztpsi(:,:,:)=zpsi(:,:,:,p)
       do iexp=1,4
          zfac=zfac*(-zI*DT)/dble(iexp)
          call zXHpsi
          zpsi(:,:,:,p)=zpsi(:,:,:,p)+zfac*zXHtpsi(:,:,:)
          ztpsi=zXHtpsi
       end do
    end do
    return
  end subroutine zNDTexp
!-------10--------20--------30--------40--------50--------60----------72
  subroutine zXHpsi
    use global_variables
    implicit none
    integer ix,iy,iz,a,i,ierr,na
    integer isx,iex,isy,iey,isz,iez
    complex(8) s0
    complex(8),allocatable :: zPc_sub(:,:),zPc_sub0(:,:)
    allocate(zPc_sub(maxproj,NI),zPc_sub0(maxproj,NI))


! operate H -----------1
    zPc_sub0=0d0
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=0d0
          do iz=isz,iez
             do iy=isy,iey
                do ix=isx,iex
                   s0=s0+ztpsi(ix,iy,iz)*proj_P(ix,iy,iz,i,a)
                end do
             end do
          end do
          zPc_sub0(i,a)=s0*H**3
       end do
    end do
    call MPI_ALLREDUCE(zPc_sub0,zPc_sub,maxproj*NI,MPI_DOUBLE_COMPLEX&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    zwk(Nsx:Nex,Nsy:Ney,Nsz:Nez)=ztpsi(:,:,:)
    call zLaplacian
    zXHtpsi=-0.5d0*zLwk+(Veff+ft*Vext+zI*Wabs)*ztpsi

    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=sum((Hc(i,:,a)+ft*rk_c(i,:,a))*zPc_sub(:,a))
          do iz=isz,iez
          do iy=isy,iey
          do ix=isx,iex
             zXHtpsi(ix,iy,iz)=zXHtpsi(ix,iy,iz)+s0*proj_P(ix,iy,iz,i,a)
          end do
          end do
          end do
       end do
    end do
! end operate H -------1
! operate X -----------1
    zPc_sub0=0d0
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=0d0
          do iz=isz,iez
             do iy=isy,iey
                do ix=isx,iex
                   s0=s0+zXHtpsi(ix,iy,iz)*proj_P(ix,iy,iz,i,a)
                end do
             end do
          end do
          zPc_sub0(i,a)=s0*H**3
       end do
    end do
    call MPI_ALLREDUCE(zPc_sub0,zPc_sub,maxproj*NI,MPI_DOUBLE_COMPLEX&
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    do a=1,NI
       isx=max(Nsx,iproj(1,1,a));iex=min(Nex,iproj(2,1,a))
       if(isx>iex)cycle
       isy=max(Nsy,iproj(1,2,a));iey=min(Ney,iproj(2,2,a))
       if(isy>iey)cycle
       isz=max(Nsz,iproj(1,3,a));iez=min(Nez,iproj(2,3,a))
       if(isz>iez)cycle
       na=Kion(a)
       do i=1,Nproj(na)
          s0=sum(Xc(i,:,na)*zPc_sub(:,a))
          do iz=isz,iez
          do iy=isy,iey
          do ix=isx,iex
             zXHtpsi(ix,iy,iz)=zXHtpsi(ix,iy,iz)+s0*proj_P(ix,iy,iz,i,a)
          end do
          end do
          end do
       end do
    end do
! end operate X -------1
    return    
  end subroutine zXHpsi

!-------10--------20--------30--------40--------50--------60----------72
!end function
    SUBROUTINE Laplacian       ! wk -> Lwk
      use global_variables
      implicit none
      integer ix,iy,iz,ierr
      integer isend_xp,isend_xm,isend_yp,isend_ym,isend_zp,isend_zm
      integer irecv_xp,irecv_xm,irecv_yp,irecv_ym,irecv_zp,irecv_zm
      integer istatus(MPI_STATUS_SIZE)

! boundary transfer

      if(Itable(My_px+1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  sendxp=wk(Nex-Nd+1:Nex,Nsy:Ney,Nsz:Nez)
      if(Itable(My_px-1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  sendxm=wk(Nsx:Nsx+Nd-1,Nsy:Ney,Nsz:Nez)
      if(Itable(My_px,My_py+1,My_pz)/=MPI_PROC_NULL)&
     &  sendyp=wk(Nsx:Nex,Ney-Nd+1:Ney,Nsz:Nez)
      if(Itable(My_px,My_py-1,My_pz)/=MPI_PROC_NULL)&
     &  sendym=wk(Nsx:Nex,Nsy:Nsy+Nd-1,Nsz:Nez)
      if(Itable(My_px,My_py,My_pz+1)/=MPI_PROC_NULL)&
     &  sendzp=wk(Nsx:Nex,Nsy:Ney,Nez-Nd+1:Nez)
      if(Itable(My_px,My_py,My_pz-1)/=MPI_PROC_NULL)&
     &  sendzm=wk(Nsx:Nex,Nsy:Ney,Nsz:Nsz+Nd-1)


      call MPI_SEND(sendxp,rsize_x,MPI_REAL8,&
     &               Itable(My_px+1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(recvxm,rsize_x,MPI_REAL8,&
     &               Itable(My_px-1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(sendxm,rsize_x,MPI_REAL8,&
     &               Itable(My_px-1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(recvxp,rsize_x,MPI_REAL8,&
     &               Itable(My_px+1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(sendyp,rsize_y,MPI_REAL8,&
     &               Itable(My_px,My_py+1,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(recvym,rsize_y,MPI_REAL8,&
     &               Itable(My_px,My_py-1,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(sendym,rsize_y,MPI_REAL8,&
     &               Itable(My_px,My_py-1,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(recvyp,rsize_y,MPI_REAL8,&
     &               Itable(My_px,My_py+1,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(sendzp,rsize_z,MPI_REAL8,&
     &               Itable(My_px,My_py,My_pz+1),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(recvzm,rsize_z,MPI_REAL8,&
     &               Itable(My_px,My_py,My_pz-1),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(sendzm,rsize_z,MPI_REAL8,&
     &               Itable(My_px,My_py,My_pz-1),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(recvzp,rsize_z,MPI_REAL8,&
     &               Itable(My_px,My_py,My_pz+1),1,&
     &               MPI_COMM_WORLD,istatus,ierr)



      if(Itable(My_px-1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  wk(Nsx-Nd:Nsx-1,Nsy:Ney,Nsz:Nez)=recvxm

      if(Itable(My_px+1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  wk(Nex+1:Nex+Nd,Nsy:Ney,Nsz:Nez)=recvxp

      if(Itable(My_px,My_py-1,My_pz)/=MPI_PROC_NULL)&
     &  wk(Nsx:Nex,Nsy-Nd:Nsy-1,Nsz:Nez)=recvym

      if(Itable(My_px,My_py+1,My_pz)/=MPI_PROC_NULL)&
     &  wk(Nsx:Nex,Ney+1:Ney+Nd,Nsz:Nez)=recvyp

      if(Itable(My_px,My_py,My_pz-1)/=MPI_PROC_NULL)&
     &  wk(Nsx:Nex,Nsy:Ney,Nsz-Nd:Nsz-1)=recvzm

      if(Itable(My_px,My_py,My_pz+1)/=MPI_PROC_NULL)&
     &  wk(Nsx:Nex,Nsy:Ney,Nez+1:Nez+Nd)=recvzp

      do iz=Nsz,Nez
      do iy=Nsy,Ney
      do ix=Nsx,Nex
        Lwk(ix,iy,iz)=&
     &       3*cN0*wk(ix,iy,iz)&
     &      +cN1*(wk(ix+1,iy,iz)+wk(ix-1,iy,iz)+wk(ix,iy+1,iz)&
     &           +wk(ix,iy-1,iz)+wk(ix,iy,iz+1)+wk(ix,iy,iz-1))&
     &      +cN2*(wk(ix+2,iy,iz)+wk(ix-2,iy,iz)+wk(ix,iy+2,iz)&
     &           +wk(ix,iy-2,iz)+wk(ix,iy,iz+2)+wk(ix,iy,iz-2))&
     &      +cN3*(wk(ix+3,iy,iz)+wk(ix-3,iy,iz)+wk(ix,iy+3,iz)&
     &           +wk(ix,iy-3,iz)+wk(ix,iy,iz+3)+wk(ix,iy,iz-3))&
     &      +cN4*(wk(ix+4,iy,iz)+wk(ix-4,iy,iz)+wk(ix,iy+4,iz)&
     &           +wk(ix,iy-4,iz)+wk(ix,iy,iz+4)+wk(ix,iy,iz-4))
      enddo
      enddo
      enddo

      Lwk=Lwk/H**2
      return
    end SUBROUTINE Laplacian
!-------10--------20--------30--------40--------50--------60----------72
    SUBROUTINE zLaplacian       ! zwk -> zLwk
      use global_variables
      implicit none
      integer ix,iy,iz,ierr
      integer isend_xp,isend_xm,isend_yp,isend_ym,isend_zp,isend_zm
      integer irecv_xp,irecv_xm,irecv_yp,irecv_ym,irecv_zp,irecv_zm
      integer istatus(MPI_STATUS_SIZE)

! boundary transfer

      if(Itable(My_px+1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  zsendxp=zwk(Nex-Nd+1:Nex,Nsy:Ney,Nsz:Nez)
      if(Itable(My_px-1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  zsendxm=zwk(Nsx:Nsx+Nd-1,Nsy:Ney,Nsz:Nez)
      if(Itable(My_px,My_py+1,My_pz)/=MPI_PROC_NULL)&
     &  zsendyp=zwk(Nsx:Nex,Ney-Nd+1:Ney,Nsz:Nez)
      if(Itable(My_px,My_py-1,My_pz)/=MPI_PROC_NULL)&
     &  zsendym=zwk(Nsx:Nex,Nsy:Nsy+Nd-1,Nsz:Nez)
      if(Itable(My_px,My_py,My_pz+1)/=MPI_PROC_NULL)&
     &  zsendzp=zwk(Nsx:Nex,Nsy:Ney,Nez-Nd+1:Nez)
      if(Itable(My_px,My_py,My_pz-1)/=MPI_PROC_NULL)&
     &  zsendzm=zwk(Nsx:Nex,Nsy:Ney,Nsz:Nsz+Nd-1)


      call MPI_SEND(zsendxp,rsize_x,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px+1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(zrecvxm,rsize_x,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px-1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(zsendxm,rsize_x,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px-1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(zrecvxp,rsize_x,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px+1,My_py,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(zsendyp,rsize_y,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py+1,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(zrecvym,rsize_y,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py-1,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(zsendym,rsize_y,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py-1,My_pz),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(zrecvyp,rsize_y,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py+1,My_pz),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(zsendzp,rsize_z,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py,My_pz+1),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(zrecvzm,rsize_z,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py,My_pz-1),1,&
     &               MPI_COMM_WORLD,istatus,ierr)

      call MPI_SEND(zsendzm,rsize_z,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py,My_pz-1),1,&
     &               MPI_COMM_WORLD,ierr)

      call MPI_RECV(zrecvzp,rsize_z,MPI_DOUBLE_COMPLEX,&
     &               Itable(My_px,My_py,My_pz+1),1,&
     &               MPI_COMM_WORLD,istatus,ierr)



      if(Itable(My_px-1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  zwk(Nsx-Nd:Nsx-1,Nsy:Ney,Nsz:Nez)=zrecvxm

      if(Itable(My_px+1,My_py,My_pz)/=MPI_PROC_NULL)&
     &  zwk(Nex+1:Nex+Nd,Nsy:Ney,Nsz:Nez)=zrecvxp

      if(Itable(My_px,My_py-1,My_pz)/=MPI_PROC_NULL)&
     &  zwk(Nsx:Nex,Nsy-Nd:Nsy-1,Nsz:Nez)=zrecvym

      if(Itable(My_px,My_py+1,My_pz)/=MPI_PROC_NULL)&
     &  zwk(Nsx:Nex,Ney+1:Ney+Nd,Nsz:Nez)=zrecvyp

      if(Itable(My_px,My_py,My_pz-1)/=MPI_PROC_NULL)&
     &  zwk(Nsx:Nex,Nsy:Ney,Nsz-Nd:Nsz-1)=zrecvzm

      if(Itable(My_px,My_py,My_pz+1)/=MPI_PROC_NULL)&
     &  zwk(Nsx:Nex,Nsy:Ney,Nez+1:Nez+Nd)=zrecvzp

      do iz=Nsz,Nez
      do iy=Nsy,Ney
      do ix=Nsx,Nex
        zLwk(ix,iy,iz)=&
     &       3*cN0*zwk(ix,iy,iz)&
     &      +cN1*(zwk(ix+1,iy,iz)+zwk(ix-1,iy,iz)+zwk(ix,iy+1,iz)&
     &           +zwk(ix,iy-1,iz)+zwk(ix,iy,iz+1)+zwk(ix,iy,iz-1))&
     &      +cN2*(zwk(ix+2,iy,iz)+zwk(ix-2,iy,iz)+zwk(ix,iy+2,iz)&
     &           +zwk(ix,iy-2,iz)+zwk(ix,iy,iz+2)+zwk(ix,iy,iz-2))&
     &      +cN3*(zwk(ix+3,iy,iz)+zwk(ix-3,iy,iz)+zwk(ix,iy+3,iz)&
     &           +zwk(ix,iy-3,iz)+zwk(ix,iy,iz+3)+zwk(ix,iy,iz-3))&
     &      +cN4*(zwk(ix+4,iy,iz)+zwk(ix-4,iy,iz)+zwk(ix,iy+4,iz)&
     &           +zwk(ix,iy-4,iz)+zwk(ix,iy,iz+4)+zwk(ix,iy,iz-4))
      enddo
      enddo
      enddo

      zLwk=zLwk/H**2
      return
    end SUBROUTINE zLaplacian      
!-------10--------20--------30--------40--------50--------60----------72
  real(8) function Yfunc(theta,phi,L,M)
    implicit none
! Yfunc:real spherical harmonics function
    integer L,M
    real(8) theta,phi
    real(8),parameter ::pi=3.1415926535897932d0
            
    select case(l)
    case(0)
       Yfunc=1.d0/dsqrt(4.d0*pi)
    case(1)
       if(M.eq.1)then
          Yfunc=dsqrt(3.d0/4.d0/pi)*dsin(theta)*dcos(phi)
       else if(M.eq.-1)then
          Yfunc=dsqrt(3.d0/4.d0/pi)*dsin(theta)*dsin(phi)
       else
          Yfunc=dsqrt(3.d0/4.d0/pi)*dcos(theta)
       end if
    case(2)
       if(M.eq.2)then
          Yfunc=dsqrt(15.d0/16.d0/pi)*dsin(theta)**2*dcos(2.d0*phi)
       else if(M.eq.-2)then
          Yfunc=dsqrt(15.d0/16.d0/pi)*dsin(theta)**2*dsin(2.d0*phi)
       else if(M.eq.1)then
          Yfunc=dsqrt(15.d0/4.d0/pi)*dsin(theta)*dcos(theta)*dcos(phi)
       else if(M.eq.-1)then
          Yfunc=dsqrt(15.d0/4.d0/pi)*dsin(theta)*dcos(theta)*dsin(phi)
       else
          Yfunc=dsqrt(5.d0/16.d0/pi)*(3.d0*dcos(theta)**2-1.d0)
       end if
    case default
       write(*,*)'error Yfunc l',l,m
       stop
    end select
    return
  end function Yfunc
!-------10--------20--------30--------40--------50--------60----------72
  subroutine inirac
    use global_variables
    integer i
    real(8) fn
    faclog(1)=0.0d0
    faclog(2)=0.0d0
    fn=1.d0
    do 200 i=3,500
       fn=fn+1.0d0
200    faclog(i)=faclog(i-1)+dlog(fn)
    return
  end subroutine inirac 
!-------10--------20--------30--------40--------50--------60----------72
  real(8) function cleb(ia,id,ib,ie,ic,if)
    use global_variables
    implicit real(8) (a-h,o-z)
    implicit integer (i-n)
    integer ia,id,ib,ie,ic,if
    integer nzt2,nzm1,ibme,iabcp,ibpe,ibca,iabc
    integer iapd,iamd,icmf

10  cleb=0.0d0
    if(id+ie-if) 7000,105,7000
105 k1=ia+ib+ic
    if((-1.)**k1) 7000,110,110
110 k1=ia+ib-ic
    k2=ic-iabs(ia-ib)
    k3=min0(k1,k2)
    if(k3) 7000,130,130
130 if((-1.)**(ib+ie)) 7000,7000,140
140 if((-1.)**(ic+if)) 7000,7000,150
150 if(ia-iabs(id)) 7000,152,152
152 if(ib-iabs(ie)) 7000,154,154
154 if(ic-iabs(if)) 7000,160,160
160 if(ia) 7000,175,165
165 if(ib) 7000,175,170
170 if(ic) 7000,180,250
175 cleb=1.0d0
    go to 7000
180 fb=ib+1
    cleb=((-1.0d0)**((ia-id)/2))/dsqrt(fb)
    go to 7000
250 fc2=ic+1
    iabcp=(ia+ib+ic)/2+1
    iabc=iabcp-ic
    icab=iabcp-ib
    ibca=iabcp-ia
    iapd=(ia+id)/2+1
    iamd=iapd-id
    ibpe=(ib+ie)/2+1
    ibme=ibpe-ie
    icpf=(ic+if)/2+1
    icmf=icpf-if
    sqfclg=0.5d0*(dlog(fc2)-faclog(iabcp+1) &
         &      +faclog(iabc)+faclog(icab)+faclog(ibca) &
         &      +faclog(iapd)+faclog(iamd)+faclog(ibpe) &
         +faclog(ibme)+faclog(icpf)+faclog(icmf))
    nzmic2=(ib-ic-id)/2
    nzmic3=(ia-ic+ie)/2
    nzmi= max0 (0,nzmic2,nzmic3)+1
    nzmx= min0 (iabc,iamd,ibpe)
    if(nzmx.lt.nzmi) go to 7000
    s1=(-1.0d0)**(nzmi-1)
    do 400 nz=nzmi,nzmx
       nzm1=nz-1
       nzt1=iabc-nzm1
       nzt2=iamd-nzm1
       nzt3=ibpe-nzm1
       nzt4=nz-nzmic2
       nzt5=nz-nzmic3
       termlg=sqfclg-faclog(nz)-faclog(nzt1)-faclog(nzt2) &
            &           -faclog(nzt3)-faclog(nzt4)-faclog(nzt5)
       ssterm=s1*dexp(termlg)
       cleb=cleb+ssterm
400    s1=-s1
         
7000 return
  end function cleb
!-------10--------20--------30--------40--------50--------60----------72
  complex(8) function Y_Y_gauss(l1,l2,m1,m2,a1,a2,Rab,theta,phi)
    use global_variables
    implicit none
    real(8) a1,a2,Rab,theta,phi,const,c1,c2
    integer l1,l2,m1,m2,l0,m0,L,kaizyo
    complex(8) a,s,f1
    real(8) gamma
    complex(8) f2,f3
    complex(8) yf,Y_func_im
    real(8) Gauss_Bessel,cleb
      
    a=(.0d0,1.d0)
    S=(.0d0,.0d0)
    gamma=(a1+a2)/(4.d0*a1*a2)
    l0=0
    m0=0
    const=dble(kaizyo(l1)*kaizyo(l2))*(4.d0*a1)**(dble(l1)+1.5d0)* &
         (4.d0*a2)**(dble(l2)+1.5d0)/dble(kaizyo(2*l1+1) &
         *kaizyo(2*l2+1))/4.d0/pi 
    do L=0,4
       if(abs(m1-m2).le.L)then
          c1=cleb(2*l2,m0,2*L,m0,2*l1,m0)
          c2=cleb(2*l2,2*m2,2*L,2*(m1-m2),2*l1,2*m1)
          f1=sqrt(dble((2*L+1)*(2*l2+1))/dble(2*l1+1)/pi*0.25d0)*c1*c2
          yf=Y_func_im(L,m1-m2,theta,phi)
          f1=f1*conjg(yf)
          f1=f1*GAuss_Bessel(Rab,L,l1+l2,Gamma)*a**L
          S=S+f1
       end if
    end do
    S=S*(.5d0/a1)**(dble(l1)+1.5d0)*(.5d0/a2)**(dble(l2)+1.5d0)
    Y_Y_gauss=S*16.d0*pi*pi*a**(l1-l2)*const
    return
  end function Y_Y_gauss
!-------10--------20--------30--------40--------50--------60----------72
  real(8) function Gauss_Bessel(Rab,L1,L2,Gamma)
    implicit none
    real(8) Rab,k,Gamma
    integer L1,L2
    real(8) alpha,x,xmin,xmax,xf,dx,epsilon,S,Sbes
    integer N,i
    
    alpha=gamma/(Rab**2)
    xmin=1d0
    xmax=10000d0
    xf=xmax
    dx=0.05d0
    epsilon=1d-12

    if(xmax**L2*exp(-alpha*xmax**2)<=epsilon)then
       xf=0.5d0*(xmax+xmin)
       do 
          if(xf**L2*exp(-alpha*xf**2)<=epsilon)then
             xmax=xf
          else
             xmin=xf
          end if
          xf=0.5d0*(xmax+xmin)
          if(xmax-xmin<=dx)exit
       end do
    end if
    
    N=xf/dx+1
    xf=dx*dble(N)
    if(L2==0)then
       S=0.5d0*(Sbes(0d0,L1)+Sbes(xf,L1)*exp(-alpha*xf**2))
    else
       S=0.5d0*(Sbes(xf,L1)*xf**L2*exp(-alpha*xf**2))         
    end if
    do i=1,N-1
       x=dx*dble(i)
       S=S+Sbes(x,L1)*x**L2*exp(-alpha*x**2)
    end do
    Gauss_Bessel=S*dx/(Rab**(L2+1))
    return
  end function Gauss_Bessel
!-------10--------20--------30--------40--------50--------60----------72
  real(8) function Sbes(x,l)
    implicit none
    ! Spherical Bessel function
    integer l
    real(8) x
    
    if(x.le.(1.0d-6))then
       if(L.eq.0)then
          Sbes=1.d0
       else
          Sbes=.0d0
       end if
       return
    end if
    
    if(l.eq.1)goto 10
    if(l.eq.2)goto 20
    if(l.eq.3)goto 30
    if(l.eq.4)goto 40
    
    sbes=dsin(x)/x
    return
10  sbes=(dsin(x)-x*dcos(x))/(x*x)
    return
20  sbes=((3.d0-x*x)*dsin(x)-3.d0*x*dcos(x))/(x**3)
    return
30  sbes=((15.d0-6.d0*x*x)*dsin(x)-(15.d0*x-x**3)*dcos(x)) &
         &/(x**4)
    return
40  sbes=((105.d0-45.d0*x*x+x**4)*dsin(x) &
         &-(105.d0*x-10.d0*x**3)*dcos(x))/(x**5)
    return
      
  end function Sbes
!-------10--------20--------30--------40--------50--------60----------72
  integer function kaizyo(N)
    implicit none
    integer N,i
    if(N.eq.0)then
       kaizyo=1
       return
    end if
    kaizyo=1
    do i=1,N
       kaizyo=kaizyo*i
    end do
    return
  end function kaizyo
!-------10--------20--------30--------40--------50--------60----------72
  real(8) function  Rkaizyo(N)
    implicit none
    integer N,ido,i
    if(N.eq.0)then
       ido=1
       Rkaizyo=dble(ido)
       return
    end if
    ido=1
    do i=1,N
       ido=ido*i
    end do
    Rkaizyo=dble(ido)
    return
  end function Rkaizyo
!-------10--------20--------30--------40--------50--------60----------72
  complex(8) function Y_func_im(l,m,theta,phi)
    implicit none
    real(8),parameter ::pi=3.1415926535897932d0
    complex(8),parameter :: zI=(0d0,1d0)
    integer l,m
    real(8) theta,phi
    complex(8) a
    a=(.0d0,1.d0)
      
    if(l==0)then
       if(m==0)then
          Y_func_im=.5d0*sqrt(1.d0/pi)
          return
       end if
       write(*,*)'error spherical harmonics',l,m
       stop
    else if(l==1)then
       if(m==(-1))then
          Y_func_im=.5d0*sqrt(3.d0/(2.d0*pi))*exp(-zI*phi)*dsin(theta)
          return
       else if(m==1)then
          Y_func_im=-.5d0*sqrt(3.d0/(2.d0*pi))*exp(zI*phi)*dsin(theta)
          return
       else if(m==0)then
          Y_func_im=.5d0*sqrt(3.d0/pi)*dcos(theta)
          return
       end if
       write(*,*)'error spherical harmonics',l,m
    else if(l==2)then
       if(m==-2)then
          Y_func_im=0.25d0*sqrt(15.d0/(2.d0*pi)) &
               *exp(-zI*2.d0*phi)*sin(theta)**2
          return
       else if(m==2)then
          Y_func_im=0.25d0*sqrt(15.d0/(2.d0*pi))*exp(zI*2.d0*phi)* &
               sin(theta)**2.d0
          return
       else if(m==-1)then
          Y_func_im=0.5d0*sqrt(15.d0/(2.d0*pi))*exp(-zI*phi)* &
               & sin(theta)*cos(theta)
          return
       else if(m.eq.1)then
          Y_func_im=-0.5d0*sqrt(15.d0/(2.d0*pi))*exp(zI*phi)* &
               & sin(theta)*cos(theta)
          return
       else if(m.eq.0)then
          Y_func_im=0.25d0*sqrt(5.d0/pi)* &
               & (3.d0*cos(theta)**2-1.d0)
          return
       end if
       write(*,*)'error spherical harmonics',l,m
       stop
    else if(l==3)then
       if(m==-3)then
          Y_func_im=0.125d0*sqrt(35.d0/pi)*exp(-zI*3.d0*phi)* &
               & sin(theta)**3
          return
       else if(m==3)then
          Y_func_im=-0.125d0*sqrt(35.d0/pi)*exp(zI*3.d0*phi)* &
               & sin(theta)**3
          return
       else if(m==-2)then
          Y_func_im=0.25d0*sqrt(105.d0/(2.d0*pi))*exp(-zI*2.d0*phi)* &
               & sin(theta)**2*cos(theta)
          return
       else if(m==2)then
          Y_func_im=0.25d0*sqrt(105.d0/(2.d0*pi))*exp(zI*2.d0*phi)* &
               & sin(theta)**2*cos(theta)
          return
       else if(m==-1)then
          Y_func_im=0.125d0*sqrt(21.d0/pi)*exp(-zI*phi)* &
               & sin(theta)*(5.d0*cos(theta)**2-1.d0)
          return
       else if(m==1)then
          Y_func_im=-0.125d0*sqrt(21.d0/pi)*exp(zI*phi)* &
               & sin(theta)*(5.d0*cos(theta)**2-1.d0)
          return
       else if(m==0)then
          Y_func_im=0.25d0*sqrt(7.d0/pi)* &
               & (5.d0*cos(theta)**3-3.d0*cos(theta))
          return
       end if
       write(*,*)'error spherical harmonics',l,m
       stop
    else if(l==4)then
       if(m==-4)then
          Y_func_im=(3.d0/16.d0)*sqrt(35.d0/(2.d0*pi)) &
               *exp(-zI*4.d0*phi)*sin(theta)**4
          return
       else if(m==4)then
          Y_func_im=(3.d0/16.d0)*sqrt(35.d0/(2.d0*pi)) &
               *exp(zI*4.d0*phi)*sin(theta)**4
          return
       else if(m==-3)then
          Y_func_im=(3.d0/8.d0)*sqrt(35.d0/pi)*exp(-zI*3.d0*phi)* &
               & sin(theta)**3*cos(theta)
          return
       else if(m==3)then
          Y_func_im=-(3.d0/8.d0)*sqrt(35.d0/pi)*exp(zI*3.d0*phi)* &
               & sin(theta)**3*cos(theta)
          return
       else if(m==-2)then
          Y_func_im=(3.d0/8.d0)*sqrt(5.d0/(2.d0*pi)) &
               *exp(-zI*2d0*phi)*sin(theta)**2*(7d0*cos(theta)**2-1d0)
          return
       else if(m==2)then
          Y_func_im=(3d0/8d0)*sqrt(5d0/(2d0*pi))*exp(zI*2d0*phi)* &
                sin(theta)**2*(7d0*cos(theta)**2-1d0)
          return
       else if(m==-1)then
          Y_func_im=(3.d0/8.d0)*sqrt(5.d0/pi)*exp(-zI*phi)* &
               & sin(theta)*(7.d0*cos(theta)**3-3.d0*cos(theta))
          return
       else if(m==1)then
          Y_func_im=-(3.d0/8.d0)*sqrt(5.d0/pi)*exp(zI*phi)* &
               & sin(theta)*(7.d0*cos(theta)**3-3.d0*cos(theta))
          return
       else if(m==0)then
          Y_func_im=(3.d0/16.d0)*dsqrt(1.d0/pi)* &
               & (35.d0*cos(theta)**4-30.d0*cos(theta)**2+3.d0)
          return
       end if
       write(*,*)'error spherical harmonics',l,m
       stop
       
    end if
    write(*,*)'error spherical harmonics',l,m
    stop
    return
  end function Y_func_im
!-------10--------20--------30--------40--------50--------60----------72
  subroutine simple_mixing
    use global_variables
    implicit none
    Veff=mixingrate*Veff+(1d0-mixingrate)*Veff_old
    Hc=mixingrate*Hc+(1d0-mixingrate)*Hc_old
    Veff_old=Veff;Hc_old=Hc
    return
  end subroutine simple_mixing
!-------10--------20--------30--------40--------50--------60----------72
  subroutine particle_number_projection
    use global_variables
    implicit none
    integer,allocatable :: int_state(:)
    complex(8),allocatable :: zpsipsi(:,:,:),za(:,:)
    real(8),allocatable :: Prob(:)
    complex(8) zs0,zs,zdet1,zdet2
    integer ip,jp,inout,n,ierr

!0=in 1=out
    allocate(zpsipsi(NST,NST,0:1),za(NST,NST))
    allocate(int_state(2*NST),Prob(0:2*NST))
    int_state=0
    Prob=0d0

    do ip=1,NST
      ztpsi(:,:,:)=zpsi(:,:,:,ip)
      call zSpsi
      do jp=ip,NST
        zs0=sum(conjg(zpsi(:,:,:,jp))*zStpsi(:,:,:))*H**3
        call MPI_ALLREDUCE(zs0,zs,1,MPI_DOUBLE_COMPLEX &
          ,MPI_SUM,MPI_COMM_WORLD,ierr)
        zpsipsi(ip,jp,0)=zs
        zpsipsi(jp,ip,0)=conjg(zs)
      end do
    end do

    do ip=1,NST
      do jp=1,NST
        zpsipsi(ip,jp,1)=-zpsipsi(ip,jp,0)
        if(ip == jp) zpsipsi(ip,jp,1)=zpsipsi(ip,jp,1)+1d0
      end do
    end do

    inout=0
    do 
      inout=inout+1
! upspin determinant
      do ip=1,NST
        do jp=1,NST
          za(ip,jp)=zpsipsi(ip,jp,int_state(ip))
        end do
      end do
      call zdeterminant(za,NST,NST,zdet1)
! down determinant
      do ip=1,NST
        do jp=1,NST
          za(ip,jp)=zpsipsi(ip,jp,int_state(ip+NST))
        end do
      end do
      call zdeterminant(za,NST,NST,zdet2)

      Prob(sum(int_state))=Prob(sum(int_state))+real(zdet1*zdet2)
      if(inout == 2**(NST*2))exit
      n=1
      call add_intstate(int_state,n,NST*2)
    end do

      if(Myrank == 0)then
        write(*,*)'particle number projection'
        write(*,'(100e26.16E3)')iter_rt*DT*0.02419d0,Prob
      end if
    
    return
  end subroutine particle_number_projection
!-------10--------20--------30--------40--------50--------60----------72
  recursive subroutine add_intstate(int_state,n,Nsize)
    implicit none
    integer :: Nsize,n
    integer :: int_state(Nsize)

    if(int_state(n) == 0)then
      int_state(n)=1
    else
      int_state(n)=0
      call add_intstate(int_state,n+1,Nsize)
    end if
    return
  end subroutine add_intstate
!-------10--------20--------30--------40--------50--------60----------72
  subroutine zdeterminant(za,n,n0,zdet)
    IMPLICIT NONE
    INTEGER :: n,k,i,j,n0
    COMPLEX(8) :: za(n0,n0),zdet,zda,zdi
    zdet=1.d0
    DO k=1,n
      zda=za(k,k)
      zdet=zdet*zda
      IF(ABS(zda) < 1d-10) RETURN
      DO i=k,n
        za(k,i)=za(k,i)/zda
      ENDDO
      DO j=1,n
        IF(j == k) CYCLE
        zdi=za(j,k)
        DO i=k,n
          za(j,i)=za(j,i)-zdi*za(k,i)
        ENDDO
      ENDDO
    ENDDO
    RETURN
  end subroutine zdeterminant
