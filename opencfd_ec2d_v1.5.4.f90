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
!  Ver 1.2.0a: 2011-7-6:   对 comput_new_U( )中的代码进行修改，程序可在Intel Fortran 32位编译器下运行而不出错 （可能还有其他问题!!!）
!  Ver 1.2.1: 2011-11-23:  Symmetry boundary condition is adding
!  Ver 1.2.2: 2011-12-9:   isothormal wall boundary is adding
!  Ver 1.3.0: 2012-4-29:   LU-SGS is supported
!  Ver 1.4.0: 2012-4-29:   k-w SST model is supported
!  Ver 1.4.1: 2012-5-11:   Bug in boundary-condition for SST model is removed
!  Ver 1.4.2: 2012-5-20:   Bug in update_buffer_onemesh() is removed
!  Ver 1.4.3: 2012-5-22:   Support SA model
!  Ver 1.5.0: 2012-5-26:   Hybrid Finite-Volume Finite-Difference Method
!  Ver 1.5.1: 2012-10-8:   Dual-time LU-SGS method
!  Ver 1.5.2: 2012-12-17:  对远场边界条件进行了补充； 恢复了亚声速入口条件
!  Ver 1.5.3: 2013-4-7:   添加了入口与出口边界条件（原先只有远场边界条件）,恢复了亚声速出口边界, 取消了KS_Farfield
!  Ver 1.5.4: 2013-6-3:  Bugs in kw-SST model are removed;  原先版本在无量纲化方面有误； 强制K, W 非负 
!-----------------------------------------------------------------------------------
 
! Consts  常量   
  module const_var
  implicit none
  real*8,parameter::  PI=3.1415926535897932d0, P_rato_limit=2.d0
  integer,parameter:: Scheme_UD1=0,Scheme_NND2=1, Scheme_UD3=2,Scheme_WENO3=3,Scheme_MUSCL2=4,Scheme_MUSCL3=5, Scheme_OMUSCL2=6
  integer,parameter:: Flux_Steger_Warming=1, Flux_HLL=2, Flux_HLLC=3,Flux_Roe=4,Flux_VanLeer=5,Flux_Ausm=6
  integer,parameter:: Reconst_Original=0,Reconst_Conservative=1,Reconst_Characteristic=2
  integer,parameter:: BC_Wall=-10, BC_Farfield=-20, BC_Symmetry_or_slidewall=-30, BC_Inlet=-40, BC_Outlet=-50
  integer,parameter:: Time_Euler1=1,Time_RK3=3,Time_LU_SGS=0,Time_dual_LU_SGS=-1
  integer,parameter:: Turbulence_NONE=0, Turbulence_BL=1, Turbulence_SA=2,Turbulence_SST=3
  integer,parameter:: LAP=2                   ! 块和块之间的交叠区(overlap)宽度 （LAP=2最高支持4阶，LAP=3最高支持6阶，LAP=4最高支持8阶精度）
  real*8,parameter::  PrT=0.9d0                ! PrT 湍流Plandtl数
  integer,parameter:: Method_FVM=0, Method_FDM=1   ! 有限体积法，差分法
  integer,parameter:: LFDM=4    ! 差分法块的边界网格（考虑到块-块连接处的网格的不光滑性，块边缘的4层网格不使用差分法）

  end module const_var
 
! 全局变量，包括参数、几何量及物理量
!========================================================================================
  module Global_Var    
  use const_var
  implicit none
  
  TYPE BC_MSG_TYPE              ! 边界链接信息
   integer::  f_no, face, ist, iend, jst, jend,  neighb, subface, orient
  END TYPE BC_MSG_TYPE

!------------------------------------网格块--------------------------------------
   TYPE Block_TYPE                               ! 数据结构：网格块 ；包含几何变量及物理变量的信息 
     integer :: Block_no,nx,ny,subface           ! 块号；网格数nx,ny；子面数
	 integer:: FVM_FDM                           ! 有限体积法/差分法
     real*8,pointer,dimension(:,:):: x,y,x1,y1   ! (x,y) : coordinate of vortex; (x1,y1): coordinate of cell center  
	 real*8,pointer,dimension(:,:):: Akx,Aky,Aix,Aiy,AJac   ! Jocabian for Finite-Difference method; 
	 real*8,pointer,dimension(:,:,:) :: U,Un,Un1     ! 守恒变量 (本时间步及前1,2个时间步的值), conversation variables 
     real*8,pointer,dimension(:,:,:) :: Res      ! 残差 （净通量）
     real*8,pointer,dimension(:,:)::    dt       ! (局部)时间步长
     real*8,pointer,dimension(:,:,:) :: QF       ! 强迫函数 (多重网格法中粗网格使用)
     real*8,pointer,dimension(:,:,:) :: deltU    ! 守恒变量的差值, dU=U(n+1)-U(n)  多重网格使用
	 real*8,pointer,dimension(:,:,:)::  dU       ! 守恒变量的插值dU=U(n+1)-U(n), LU-SGS方法中使用
	 real*8,pointer,dimension(:,:):: Amu,Amu_t ! 层流粘性系数、湍流粘性系数
	        ! 几何量: vol控制体面积; si,sj i,j-方向控制体边界长度; (ni1,ni2) i-方向控制体边界法方向; (nj1,nj2) j-方向控制体边界法方向;
     real*8,pointer,dimension(:,:):: vol,si,sj,ni1,ni2,nj1,nj2        
     real*8,pointer,dimension(:,:):: Lci,Lcj,Lvi,Lvj       ! 谱半径 (Lci,Lcj i-,j-方向无粘项谱半径; Lvi,Lvj i-, j-方向粘性项谱半径)
     real*8,pointer,dimension(:,:):: dw                    ! 到壁面的距离 （k-w SST模型中使用）
	 TYPE(BC_MSG_TYPE),pointer,dimension(:)::bc_msg    ! 边界链接信息
   End TYPE Block_TYPE  

!---------------------------网格 -------------------------------------------------------- 
                                      ! (如单重网格，只有1套；如多重网格，可以有多套) 
   TYPE Mesh_TYPE                     ! 数据结构“网格”； 包含几何变量及物理变量信息
     integer:: Mesh_no,Num_Block, Num_Cell,Kstep          ! 网格编号 (1号为最细网格，2号为粗网格， 3号为更粗网格...)，网格块数，网格数目,时间步 
     integer:: Nvar  ! 变量（方程）数目，如使用k-w模型，则有6个变量 （稀网格不使用湍模型，因而变量仍是4个）
     real*8::  Res_max(6),Res_rms(6),tt                 ! 最大残差，均方根残差, 推进的时间
	 TYPE (Block_TYPE),pointer,dimension(:):: Block     ! “网格块”  （从属于“网格”）

!                                                       控制参数，用于控制数值方法、通量技术、湍流模型等    
!             这些控制参数从属于“网格”，不同“网格”可以采用不同的计算方法、湍流模型等。	 （例如，粗网格用低精度方法，粗网格不使用湍流模型,...）
	integer::   Iflag_turbulence_model,  Iflag_Scheme,IFlag_flux,IFlag_Reconstruction
   End TYPE Mesh_TYPE
!---------------------------------------------------------------------------------------------


! global variables                                       各子程序均可见的全局变量
!----------------------------------------------------------------------------
   TYPE (Mesh_TYPE),pointer,dimension(:):: Mesh       ! 主数据 “网格”
   integer,save:: Num_Mesh                            ! 网格的套数  
   integer,save:: Nvar   
   integer,save::  Kstep_save,If_viscous, Iflag_turbulence_model,Iflag_init,  &
      Iflag_Scheme,IFlag_flux,Iflag_local_dt,IFlag_Reconstruction,Time_Method,Kstep_show,Num_Threads
   real*8,save:: Ma,Re,gamma,Pr,AoA,Cp,Cv,t_end,p_outlet,T_inf,Twall,vt_inf,Kt_inf,Wt_inf
   real*8,save :: dt_global,CFL,dtmax,dtmin                                        ! 与时间步长有关的量
   real*8,save::  Res_Inner_Limit     ! 内迭代残差下限
   integer,save:: Pre_Step_Mesh(3), Nstep_Inner_Limit                 ! 构建初值时，粗网格预迭代步数； 内迭代步数限制

!-----------------------------------------------
 !   全局控制参数，控制数值方法、通量技术及湍流模型等 （有些只对最细网格有效）
 ! Nvar方程（变量）的数目， 如使用BL模型Nvar=4;  使用SA模型 Nvar=5, 如使用K-W SST模型 Nvar=6 (4个基本方程+k方程+w方程）                                                        
 ! 控制变量  If_viscous=0 Euler方程，1 N_S方程；
 ! Iflag_turbulence_model 湍流模型（BL,SA,SST）;
 ! Iflag_Scheme 数值格式；
 ! Iflag_flux 通量技术； 
 ! Iflag_local_dt 是否采用局部时间步长 ; 
 !Num_Threads OpenMP采用的线程数
 ! global parameter (for all Meshes )                     流动参数, 对全体“网格”都适用
 ! Ma: Mach数;  Re: Reynolds数; Pr: Prandtl数; Cp,Cv: 定压、定容比热;
 ! t_end: End time 计算结束的时间; P_outlet: 出口压力（亚声速内流计算必须给定,无量纲量）;  Twall: 壁温（有量纲量，单位K）
 ! T_inf: 来流温度 (有量纲值，单位K),在surthland公式中使用; 
 ! Kt_inf, Wt_inf: 来流湍动能、湍能比耗散率  (Amut_inf=Kt_inf/Wt_inf), SST模型的初值及入口边界条件使用
 ! vt_inf: 来流的湍流粘性系数（与层流粘性系数之比）， SA模型使用
  end module Global_Var
   
!----------------------------------------------------------------------------
! 流场物理量 （计算每块时申请内存，该块计算结束后释放；属于临时变量） 
! 这些物理量的值无需保留;  
module Flow_Var
   real*8, save,pointer,dimension(:,:)::  d,uu,v,T,p,cc  ! 密度、x-速度、y-速度、压力、声速；
   real*8, save,pointer,dimension(:,:,:):: Fluxi,Fluxj         ! i- 及j-方向的通量
end module Flow_Var



!------------------------------------------------------------------------------------------
! 主程序 主程序 主程序
!-----------------------------------------------------------------------------------------
  program main
  use Global_Var
  implicit none
  integer  nMesh
    print*,  "----------------- OpenCFD-EC2D ver 1.4.3-------------------------------"
	print*,  "        Copyright by Li Xinliang, lixl@imech.ac.cn                     "
	print*,  "        Programming by Li Xinliang  2012-5-22                          "
    print*,  "-----------------------------------------------------------------------" 

    call read_parameter                     ! 读取流动参数及控制信息

!$ call omp_set_num_threads(NUM_THREADS)   ! 设置OpenMP的运行线程数 （并行数目）， 本语句对openmp编译器不是注释!
!$OMP Parallel
   print*, "run ..."                        ! 测试一下运行的线程数（如果OpenMP正常，则打印出NUM_THREADS个 "run ..."）
!$OMP END parallel

    call check_mesh_multigrid               ! 检查网格配置所允许的最大重数,并设定多重网格的重数
    call Init                               ! 初始化，创建数据结构; 读入几何及物理信息
    call set_control_para                   ! 设定各重网格上的控制信息（数值方法、通量技术、湍流模型、时间推进方式）
	call Update_coordinate_buffer           ! 利用连接信息，给出虚网格的坐标
    call Init_FiniteDifference              ! 设定有限差分法的区域，计算Jocabian变换系数
	
	if (Iflag_turbulence_model .eq. Turbulence_SST .or. Iflag_turbulence_model .eq. Turbulence_SA )  call comput_dw        ! 计算各网格（中心）点到壁面的距离 （采用SA或SST模型时需要此计算)
    call Init_flow

   
   print*, " Start ......"

!------------------------------------------------------------------------
! 时间推进，采用单重网格、二重网格或三重网格； 采用1阶Euler,3阶RK或LU-SGS
  do while(Mesh(1)%tt .lt. t_end )
   
    if(Num_Mesh .eq. 1)  then                          ! 单重网格时间推进 (1阶Euler,3阶RK, LU-SGS)
       call NS_Time_advance(1)
	else  if(Num_Mesh .eq. 2)  then                    ! 时间推进，采用2重网格法
  	  call NS_2stge_multigrid
    else                                               ! 时间推进，采用3重网格法
  	  call NS_3stge_multigrid 
    endif	  
  

   if(mod(Mesh(1)%Kstep, Kstep_show).eq.0)    call output_Res(1)         ! 打印残差 (最密网格)
   if(mod(Mesh(1)%Kstep, Kstep_save) .eq. 0)   call output (1)           ! 输出流场 (最密网格)
  
  enddo
    call output (1)



  end
!----------------------------------------------------------------------------------------------



!  -----------------------读取流动参数及控制变量----------------------
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
!  在本版本（ver 1.3），LU_SGS方法不支持多重网格，用户不能同时使用这两种加速收敛技术
!  后续版本中，作者将进行改进    
	if( (Time_Method .eq. Time_LU_SGS .or. Time_Method .eq. Time_dual_LU_SGS)  .and. Num_Mesh .ne. 1) then 

	   print*, " In this version (ver 1.5.1 ), LU_SGS method Do Not support Multigrid !!! "
	   print*, "Please modify 'control.in' to choose single-grid or other time method"
	   stop
	endif
    if(Iflag_turbulence_model .eq. Turbulence_SST ) then
	  Nvar=6           ! 6个变量 (4个流场变量+湍动能k+比耗散率w)
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
!  设定各重网格上的控制信息
   subroutine set_control_para
   use Global_var
   implicit none
   integer nMesh
   TYPE (Mesh_TYPE),pointer:: MP
     MP=>Mesh(1)            ! 最细的网格
 !   最细网格上的控制参数与主控制参数相同
     MP%Iflag_turbulence_model=Iflag_turbulence_model
     MP%Iflag_Scheme=Iflag_Scheme
     MP%IFlag_flux=IFlag_flux
     MP%IFlag_Reconstruction=IFlag_Reconstruction
     MP%Nvar=Nvar   ! 变量（方程）数目（最密网格使用湍流模型，数目与Nvar相同）
!   设定粗网格上的控制参数
   do nMesh=2,Num_Mesh
     MP=>Mesh(nMesh)
     MP%Iflag_turbulence_model=Turbulence_NONE    ! 粗网格不使用湍流模型
     MP%Iflag_Scheme=Scheme_UD1                   ! 粗网格使用1阶迎风格式
     MP%IFlag_flux=IFlag_flux                     ! 粗网格的通量分裂技术、时间推进近似及重构技术与细网格相同
     MP%IFlag_Reconstruction=IFlag_Reconstruction
     MP%Nvar=4          ! 变量（方程）数目 （粗网格不使用湍流模型，方程数目为4）
   enddo
  end subroutine set_control_para



!---------------------------------------
!         打印残差（最大残差和均方根残差）
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
!  输出几何及物理量 （tecplot格式）, 最细网格flow2d.dat; 粗网格 flow2d-2.dat ; 最粗网格 flow2d-3.dat
   
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
