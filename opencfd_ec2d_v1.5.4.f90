!----OpenCFD-EC 2d--------------------------------------------------------------
!  Copyright by Li Xinliang, LHD, Institute of Mechanics, lixl@imech.ac.cn                 
!  A multi-block Finite Volume Method Navier-Stokes solver
!  Ref. J. Blazek's Book: Computational Fluid Dynamics: principle and application 
!----------------------------------------------------------------------------------
!  Ver 0.91: 2010-8-19:  Goemetry message is not save
!  Ver 0.92: 2010-8-20:  New feature of Steger-Warming spliting    
!  Ver 0.94: 2010-8-23:  Support Multi-Block mesh  (Cell-Vetex finite volume method)
!  Ver 0.95: 2010-9-4:   Cell-Centred Finite Volume method
!  Ver 0.96: 2010-9-6:   Reconstruction by using original value 
!  Ver 0.97: 2010-9-17:  BL Turbulent model is adding
!  Ver 0.971: 2010-9-23: Bug in viscous flux (j) is removed; User can choose using the original variables or Conservative variables
!  Ver 0.98: User can choose using characteristic variables in the reconstruction
!  Ver 0.99: New Residual, MUSCL 3
!  Ver 1.0 : 2010-11-30:   Boundary conditions (Far field) are modified
!  Ver 1.0a: 2011-3-21:    A bug in Roe splitting is removed.
!  Ver 1.0b: 2011-4-13:    A bug in subroutine boundary_Farfield( ) is removed
!  Ver 1.1:  2011-3-20:    Multigrid method
!  Ver 1.1.1: 2011-4-14:   Outlet boudary condition is modified, P_outlet is input (if <0 P_outlet=p1)
!  Ver 1.1.2: 2011-5-6:    A bug in Roe splitting is removed
!  Ver 1.1.3: 2011-5-22:   OMUSCL2 is adding
!  Ver 1.2.0: 2011-5-22:   Ghost Cell with LAP cells is used (for original version only one Ghost Cell is used)
!  Ver 1.2.0a: 2011-7-6:   �� comput_new_U( )�еĴ�������޸ģ��������Intel Fortran 32λ�����������ж������� �����ܻ�����������!!!��
!  Ver 1.2.1: 2011-11-23:  Symmetry boundary condition is adding
!  Ver 1.2.2: 2011-12-9:   isothormal wall boundary is adding
!  Ver 1.3.0: 2012-4-29:   LU-SGS is supported
!  Ver 1.4.0: 2012-4-29:   k-w SST model is supported
!  Ver 1.4.1: 2012-5-11:   Bug in boundary-condition for SST model is removed
!  Ver 1.4.2: 2012-5-20:   Bug in update_buffer_onemesh() is removed
!  Ver 1.4.3: 2012-5-22:   Support SA model
!  Ver 1.5.0: 2012-5-26:   Hybrid Finite-Volume Finite-Difference Method
!  Ver 1.5.1: 2012-10-8:   Dual-time LU-SGS method
!  Ver 1.5.2: 2012-12-17:  ��Զ���߽����������˲��䣻 �ָ����������������
!  Ver 1.5.3: 2013-4-7:   ������������ڱ߽�������ԭ��ֻ��Զ���߽�������,�ָ��������ٳ��ڱ߽�, ȡ����KS_Farfield
!  Ver 1.5.4: 2013-6-3:  Bugs in kw-SST model are removed;  ԭ�Ȱ汾�������ٻ��������� ǿ��K, W �Ǹ� 
!-----------------------------------------------------------------------------------
 
! Consts  ����   
  module const_var
  implicit none
  real*8,parameter::  PI=3.1415926535897932d0, P_rato_limit=2.d0
  integer,parameter:: Scheme_UD1=0,Scheme_NND2=1, Scheme_UD3=2,Scheme_WENO3=3,Scheme_MUSCL2=4,Scheme_MUSCL3=5, Scheme_OMUSCL2=6
  integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_VanLeer=5,Flux_Ausm=6
  integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2
  integer,parameter:: BC_Wall=-10, BC_Farfield=-20, BC_Symmetry_or_slidewall=-30, BC_Inlet=-40, BC_Outlet=-50
  integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0,Time_dual_LU_SGS=-1
  integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2,Turbulence_SST=3
  integer,parameter:: LAP=2                   ! ��Ϳ�֮��Ľ�����(overlap)��� ��LAP=2���֧��4�ף�LAP=3���֧��6�ף�LAP=4���֧��8�׾��ȣ�
  real*8,parameter::  PrT=0.9d0                ! PrT ����Plandtl��
  integer,parameter:: Method_FVM=0, Method_FDM=1   ! �������������ַ�
  integer,parameter:: LFDM=4    ! ��ַ���ı߽����񣨿��ǵ���-�����Ӵ�������Ĳ��⻬�ԣ����Ե��4������ʹ�ò�ַ���

  end module const_var
 
! ȫ�ֱ�����������������������������
!========================================================================================
  module Global_Var    
  use const_var
  implicit none
  
  TYPE BC_MSG_TYPE              ! �߽�������Ϣ
   integer::  f_no, face, ist, iend, jst, jend,  neighb, subface, orient
  END TYPE BC_MSG_TYPE

!------------------------------------�����--------------------------------------
   TYPE Block_TYPE                               ! ���ݽṹ������� ���������α����������������Ϣ 
     integer :: Block_no,nx,ny,subface           ! ��ţ�������nx,ny��������
	 integer:: FVM_FDM                           ! ���������/��ַ�
     real*8,pointer,dimension(:,:):: x,y,x1,y1   ! (x,y) : coordinate of vortex; (x1,y1): coordinate of cell center  
	 real*8,pointer,dimension(:,:):: Akx,Aky,Aix,Aiy,AJac   ! Jocabian for Finite-Difference method; 
	 real*8,pointer,dimension(:,:,:) :: U,Un,Un1     ! �غ���� (��ʱ�䲽��ǰ1,2��ʱ�䲽��ֵ), conversation variables 
     real*8,pointer,dimension(:,:,:) :: Res      ! �в� ����ͨ����
     real*8,pointer,dimension(:,:)::    dt       ! (�ֲ�)ʱ�䲽��
     real*8,pointer,dimension(:,:,:) :: QF       ! ǿ�Ⱥ��� (���������д�����ʹ��)
     real*8,pointer,dimension(:,:,:) :: deltU    ! �غ�����Ĳ�ֵ, dU=U(n+1)-U(n)  ��������ʹ��
	 real*8,pointer,dimension(:,:,:)::  dU       ! �غ�����Ĳ�ֵdU=U(n+1)-U(n), LU-SGS������ʹ��
	 real*8,pointer,dimension(:,:):: Amu,Amu_t ! ����ճ��ϵ��������ճ��ϵ��
	        ! ������: vol���������; si,sj i,j-���������߽糤��; (ni1,ni2) i-���������߽編����; (nj1,nj2) j-���������߽編����;
     real*8,pointer,dimension(:,:):: vol,si,sj,ni1,ni2,nj1,nj2        
     real*8,pointer,dimension(:,:):: Lci,Lcj,Lvi,Lvj       ! �װ뾶 (Lci,Lcj i-,j-������ճ���װ뾶; Lvi,Lvj i-, j-����ճ�����װ뾶)
     real*8,pointer,dimension(:,:):: dw                    ! ������ľ��� ��k-w SSTģ����ʹ�ã�
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg    ! �߽�������Ϣ
   End TYPE Block_TYPE  

!---------------------------���� -------------------------------------------------------- 
                                      ! (�絥������ֻ��1�ף���������񣬿����ж���) 
   TYPE Mesh_TYPE                     ! ���ݽṹ�����񡱣� �������α��������������Ϣ
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep          ! ������ (1��Ϊ��ϸ����2��Ϊ������ 3��Ϊ��������...)�����������������Ŀ,ʱ�䲽 
     integer:: Nvar  ! ���������̣���Ŀ����ʹ��k-wģ�ͣ�����6������ ��ϡ����ʹ����ģ�ͣ������������4����
     real*8::  Res_max(6),Res_rms(6),tt                 ! ���в�������в�, �ƽ���ʱ��
	 TYPE (Block_TYPE),pointer,dimension(:):: Block     ! ������顱  �������ڡ����񡱣�

!                                                       ���Ʋ��������ڿ�����ֵ������ͨ������������ģ�͵�    
!             ��Щ���Ʋ��������ڡ����񡱣���ͬ�����񡱿��Բ��ò�ͬ�ļ��㷽��������ģ�͵ȡ�	 �����磬�������õ;��ȷ�����������ʹ������ģ��,...��
	integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction
   End TYPE Mesh_TYPE
!---------------------------------------------------------------------------------------------


! global variables                                       ���ӳ�����ɼ���ȫ�ֱ���
!----------------------------------------------------------------------------
   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh       ! ������ ������
   integer,save:: Num_Mesh                            ! ���������  
   integer,save:: Nvar   
   integer,save::  Kstep_save,If_viscous, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method,Kstep_show,Num_Threads
   real*8,save:: Ma,Re,gamma,Pr,AoA,Cp,Cv,t_end,p_outlet,T_inf,Twall,vt_inf,Kt_inf,Wt_inf
   real*8,save :: dt_global,CFL,dtmax,dtmin                                        ! ��ʱ�䲽���йص���
   real*8,save::  Res_Inner_Limit     ! �ڵ����в�����
   integer,save:: Pre_Step_Mesh(3), Nstep_Inner_Limit                 ! ������ֵʱ��������Ԥ���������� �ڵ�����������

!-----------------------------------------------
 !   ȫ�ֿ��Ʋ�����������ֵ������ͨ������������ģ�͵� ����Щֻ����ϸ������Ч��
 ! Nvar���̣�����������Ŀ�� ��ʹ��BLģ��Nvar=4;  ʹ��SAģ�� Nvar=5, ��ʹ��K-W SSTģ�� Nvar=6 (4����������+k����+w���̣�                                                        
 ! ���Ʊ���  If_viscous=0 Euler���̣�1 N_S���̣�
 ! Iflag_turbulence_model ����ģ�ͣ�BL,SA,SST��;
 ! Iflag_Scheme ��ֵ��ʽ��
 ! Iflag_flux ͨ�������� 
 ! Iflag_local_dt �Ƿ���þֲ�ʱ�䲽�� ; 
 !Num_Threads OpenMP���õ��߳���
 ! global parameter (for all Meshes )                     ��������, ��ȫ�塰���񡱶�����
 ! Ma: Mach��;  Re: Reynolds��; Pr: Prandtl��; Cp,Cv: ��ѹ�����ݱ���;
 ! t_end: End time ���������ʱ��; P_outlet: ����ѹ������������������������,����������;  Twall: ���£�������������λK��
 ! T_inf: �����¶� (������ֵ����λK),��surthland��ʽ��ʹ��; 
 ! Kt_inf, Wt_inf: �����Ķ��ܡ����ܱȺ�ɢ��  (Amut_inf=Kt_inf/Wt_inf), SSTģ�͵ĳ�ֵ����ڱ߽�����ʹ��
 ! vt_inf: ����������ճ��ϵ���������ճ��ϵ��֮�ȣ��� SAģ��ʹ��
  end module Global_Var
   
!----------------------------------------------------------------------------
! ���������� ������ÿ��ʱ�����ڴ棬�ÿ����������ͷţ�������ʱ������ 
! ��Щ��������ֵ���豣��;  
module Flow_Var
   real*8, save,pointer,dimension(:,:)::  d,uu,v,T,p,cc  ! �ܶȡ�x-�ٶȡ�y-�ٶȡ�ѹ�������٣�
   real*8, save,pointer,dimension(:,:,:):: Fluxi,Fluxj         ! i- ��j-�����ͨ��
end module Flow_Var



!------------------------------------------------------------------------------------------
! ������ ������ ������
!-----------------------------------------------------------------------------------------
  program main
  use Global_Var
  implicit none
  integer  nMesh
    print*,  "----------------- OpenCFD-EC2D ver 1.4.3-------------------------------"
	print*,  "        Copyright by Li Xinliang, lixl@imech.ac.cn                     "
	print*,  "        Programming by Li Xinliang  2012-5-22                          "
    print*,  "-----------------------------------------------------------------------" 

    call read_parameter                     ! ��ȡ����������������Ϣ

!$ call omp_set_num_threads(NUM_THREADS)   ! ����OpenMP�������߳��� ��������Ŀ���� ������openmp����������ע��!
!$OMP Parallel
   print*, "run ..."                        ! ����һ�����е��߳��������OpenMP���������ӡ��NUM_THREADS�� "run ..."��
!$OMP END parallel

    call check_mesh_multigrid               ! �������������������������,���趨�������������
    call Init                               ! ��ʼ�����������ݽṹ; ���뼸�μ�������Ϣ
    call set_control_para                   ! �趨���������ϵĿ�����Ϣ����ֵ������ͨ������������ģ�͡�ʱ���ƽ���ʽ��
	call Update_coordinate_buffer           ! ����������Ϣ�����������������
    call Init_FiniteDifference              ! �趨���޲�ַ������򣬼���Jocabian�任ϵ��
	
	if (Iflag_turbulence_model .eq. Turbulence_SST .or. Iflag_turbulence_model .eq. Turbulence_SA )  call comput_dw        ! ������������ģ��㵽����ľ��� ������SA��SSTģ��ʱ��Ҫ�˼���)
    call Init_flow

   
   print*, " Start ......"

!------------------------------------------------------------------------
! ʱ���ƽ������õ������񡢶���������������� ����1��Euler,3��RK��LU-SGS
  do while(Mesh(1)%tt .lt. t_end )
   
    if(Num_Mesh .eq. 1)  then                          ! ��������ʱ���ƽ� (1��Euler,3��RK, LU-SGS)
       call NS_Time_advance(1)
	else  if(Num_Mesh .eq. 2)  then                    ! ʱ���ƽ�������2������
  	  call NS_2stge_multigrid
    else                                               ! ʱ���ƽ�������3������
  	  call NS_3stge_multigrid 
    endif	  
  

   if(mod(Mesh(1)%Kstep, Kstep_show).eq.0)    call output_Res(1)         ! ��ӡ�в� (��������)
   if(mod(Mesh(1)%Kstep, Kstep_save) .eq. 0)   call output (1)           ! ������� (��������)
  
  enddo
    call output (1)



  end
!----------------------------------------------------------------------------------------------



!  -----------------------��ȡ�������������Ʊ���----------------------
  subroutine read_parameter
  use Global_var
  implicit none
  integer:: k
  open(99,file="control.in")
    read(99,*)
    read(99,*) Ma, Re, gamma, AoA,Pr,t_end,Kstep_save,If_viscous,Iflag_turbulence_model,Iflag_init
	read(99,*)
	read(99,*) Iflag_local_dt,dt_global,CFL,dtmax,dtmin,Time_Method,P_outlet,T_inf,Twall,vt_inf,Kt_inf,Wt_inf
    read(99,*)
    read(99,*) Iflag_Scheme,Iflag_Flux,IFlag_Reconstruction,Kstep_show
    read(99,*)
	read(99,*) Num_Mesh,Num_Threads,Nstep_Inner_Limit, Res_Inner_Limit
	read(99,*)
	read(99,*) (Pre_Step_Mesh(k), k=1,Num_Mesh)
    close(99)
!-----------------------------------------------------------------------------
!  �ڱ��汾��ver 1.3����LU_SGS������֧�ֶ��������û�����ͬʱʹ�������ּ�����������
!  �����汾�У����߽����иĽ�    
	if( (Time_Method .eq. Time_LU_SGS .or. Time_Method .eq. Time_dual_LU_SGS)  .and. Num_Mesh .ne. 1) then 

	   print*, " In this version (ver 1.5.1 ), LU_SGS method Do Not support Multigrid !!! "
	   print*, "Please modify 'control.in' to choose single-grid or other time method"
	   stop
	endif
    if(Iflag_turbulence_model .eq. Turbulence_SST ) then
	  Nvar=6           ! 6������ (4����������+�Ķ���k+�Ⱥ�ɢ��w)
	else if(Iflag_turbulence_model .eq. Turbulence_SA ) then
	  Nvar=5
    else
	  Nvar=4
	endif
!-----------------------------------------------
     AoA=AoA*PI/180.d0
     Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma) 
	 Cp=Cv*gamma 
	 Twall=Twall/T_inf

   end subroutine read_parameter

!------------------------------------------------------------------    
!  �趨���������ϵĿ�����Ϣ
   subroutine set_control_para
   use Global_var
   implicit none
   integer nMesh
   TYPE (Mesh_TYPE),pointer:: MP
     MP=>Mesh(1)            ! ��ϸ������
 !   ��ϸ�����ϵĿ��Ʋ����������Ʋ�����ͬ
     MP%Iflag_turbulence_model=Iflag_turbulence_model
     MP%Iflag_Scheme=Iflag_Scheme
     MP%IFlag_flux=IFlag_flux
     MP%IFlag_Reconstruction=IFlag_Reconstruction
     MP%Nvar=Nvar   ! ���������̣���Ŀ����������ʹ������ģ�ͣ���Ŀ��Nvar��ͬ��
!   �趨�������ϵĿ��Ʋ���
   do nMesh=2,Num_Mesh
     MP=>Mesh(nMesh)
     MP%Iflag_turbulence_model=Turbulence_NONE    ! ������ʹ������ģ��
     MP%Iflag_Scheme=Scheme_UD1                   ! ������ʹ��1��ӭ���ʽ
     MP%IFlag_flux=IFlag_flux                     ! �������ͨ�����Ѽ�����ʱ���ƽ����Ƽ��ع�������ϸ������ͬ
     MP%IFlag_Reconstruction=IFlag_Reconstruction
     MP%Nvar=4          ! ���������̣���Ŀ ��������ʹ������ģ�ͣ�������ĿΪ4��
   enddo
  end subroutine set_control_para



!---------------------------------------
!         ��ӡ�в���в�;������в
    subroutine output_Res(nMesh)
    use Global_var
	implicit none
	integer:: nMesh
      print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
	  print*, "----------The Max Residuals are-------- ", " ---Mesh---",nMesh
      write(*, "(4E20.10)") Mesh(nMesh)%Res_max(1:NVAR)
      print*, "  The R.M.S Residuals are "
      write(*, "(4E20.10)") Mesh(nMesh)%Res_rms(1:NVAR)
      open(99,file="Residual.dat",position="append")
      write(99,"(I8,13E24.10)") Mesh(nMesh)%Kstep, Mesh(nMesh)%Res_max(:),Mesh(nMesh)%Res_rms(1:NVAR)
      close(99) 

    end

!----------------------------------------------------------------------
!  ������μ������� ��tecplot��ʽ��, ��ϸ����flow2d.dat; ������ flow2d-2.dat ; ������� flow2d-3.dat
   
      subroutine output (nMesh)
       use Global_Var
       implicit none
       TYPE (Mesh_TYPE),pointer:: MP
       Type (Block_TYPE),pointer:: B
       integer  nMesh,i,j,m
       real*8:: d1,u1,v1,T1,p1,s1
       character(len=50):: filename,filename1
       
	   if(nMesh .eq. 1) then
          filename="flow2d.dat"
	   else
          write(filename,"('flow2d-'I1.1'.dat')") nMesh    ! flow2d-2.dat ; flow2d-3.dat 
       endif

       if(nMesh .eq. 1) then
         write(filename1,"('flow2d-'I7.7'.dat')") Mesh(1)%Kstep
	     MP=>Mesh(1)
         open(99,file=filename1)
         write(99,*) "variables=x,y,d,u,v,T,p"
         do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx+1, " j= ", B%ny+1
  	     do j=0,B%ny
	     do i=0,B%nx
           d1=B%U(1,i,j)
           u1=B%U(2,i,j)/d1
           v1=B%U(3,i,j)/d1
           T1=(B%U(4,i,j)-0.5d0*d1*(u1*u1+v1*v1))/(Cv*d1)
		   p1=d1*T1/(gamma*Ma*Ma)
           write(99,"(7f20.10)") B%x1(i,j),B%y1(i,j), d1,u1,v1,T1,p1
         enddo
         enddo
         enddo
		close(99)  
       endif

	   print*, "write data file ...", filename
       
	   MP=>Mesh(nMesh)
       open(99,file=filename)
       write(99,*) "variables=x,y,d,u,v,T,p,Amut"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx+1, " j= ", B%ny+1
  	     do j=0,B%ny
	     do i=0,B%nx
           d1=B%U(1,i,j)
           u1=B%U(2,i,j)/d1
           v1=B%U(3,i,j)/d1
           T1=(B%U(4,i,j)-0.5d0*d1*(u1*u1+v1*v1))/(Cv*d1)
         write(99,"(8f20.10)") B%x1(i,j),B%y1(i,j), d1,u1,v1,T1,d1*T1/(gamma*Ma*Ma),B%Amu_t(i,j)*Re
         enddo
         enddo
        enddo
       close(99)

      if(MP%Nvar .eq. 5) then
      open(99,file="SA2d.dat")
      write(99,*) "variables=x,y,vt"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx+1, " j= ", B%ny+1
  	     do j=0,B%ny
	     do i=0,B%nx
         write(99,"(3f20.10)") B%x1(i,j),B%y1(i,j),B%U(5,i,j) 
         enddo
         enddo
        enddo
       close(99)
      endif

      if(MP%Nvar .eq. 6) then
      open(99,file="SST2d.dat")
       write(99,*) "variables=x,y,Kt,Wt"
        do m=1,MP%Num_Block
         B => Mesh(nMesh)%Block(m)
         write(99,*)  " zone ", "i= ", B%nx+1, " j= ", B%ny+1
  	     do j=0,B%ny
	     do i=0,B%nx
         write(99,"(4f20.10)") B%x1(i,j),B%y1(i,j),B%U(5,i,j),B%U(6,i,j) 
         enddo
         enddo
        enddo
       close(99)
       endif

     end subroutine output

!----------------------------------------------------------------------------
   include "sub_Residual.f90"
   include "sub_flux_split.f90"
   include "sub_boundary.f90"
   include "sub_init.f90"
   include "sub_NS_singlegrid.f90"
   include "sub_NS_multigrid.f90"
   include "sub_LU_SGS.f90"  
   include "sub_turbulence_BL.f90"
   include "sub_turbulence_SST.f90"
   include "sub_turbulence_SA.f90"
   include "sub_Finite_Difference.f90"
