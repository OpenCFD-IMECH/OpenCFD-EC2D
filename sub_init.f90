! Copyright By Li Xinliang, Code by Li Xinliang
! 初始化：包括创建数据结构及赋初值
! 对于多重网格，根据上级网格的信息，创建各级网格
!---------------------------------------------------------------------------------
!------------------------------------------------------------------------------     
   subroutine init
   use Global_var
   implicit none
   integer :: i,j,k,m,nx1,ny1,Num_Block1,ksub,Kmax_grid
   real*8,allocatable,dimension(:,:):: x1,y1
   integer,allocatable,dimension(:):: NI,NJ
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
 !--------------------------------------------------------------------
 ! initial of const variables
     Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma)

!--------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
   allocate( Mesh(Num_Mesh) )             ! 主数据结构： “网格” （其成员是“网格块”）
   call Creat_Mesh1            ! 创建最细的网格 (从网格文件Mesh2d.dat)
   call read_bcin              ! 读网格连接信息 (bc2d.in)
   if(Num_Mesh .ge. 2) then
     call Creat_Mesh(1,2)    ! 根据1号网格（最细网格）信息，创建2号网格（粗网格）
   endif
   if(Num_Mesh .ge. 3) then
     call Creat_Mesh(2,3)    ! 根据2号网格信息（粗网格）， 创建3号网格（最粗网格）
   endif


 end   
!---------------------------------------------------------



!  检查网格是否适用于多重网格 
!  单方向网格数= 2*K+1 可用2重网格，=4*K+1 可用3重网格，=8*K+1 可用4重网格 ...
  
     subroutine check_mesh_multigrid 
	 use Global_var
	 implicit none
     integer,allocatable,dimension(:):: NI,NJ
     integer:: NB,NST,m,k,NN,Km,Km_grid,N_Cell,Ntmp,Bsub,ksub
     integer:: f_no, face, ist, iend, jst, jend,  neighb, subface, orient
	 print*, "Read Mesh2d.dat, Check if Multi-Grid can be used ..."
	 open(99,file="Mesh2d.dat")
     read(99,*) NB
	 allocate(NI(NB),NJ(NB))
     read(99,*) (NI(m), NJ(m), m=1,NB)
     close(99)
	 N_Cell=0
	 Km_grid=NI(1)  ! 初始值
     
	 do m=1,NB 
	   N_Cell=N_Cell+(NI(m)-1)*(NJ(m)-1)  ! 统计总网格单元数
 
 !    判断可使用的网格重数      
 	    Km=1
	    NN=2
!                                判断准则： 网格数-1 能被2**km 整除， 且最稀的网格单元数不小于2
       do while( mod((NI(m)-1),NN) .eq. 0 .and. (NI(m)-1)/NN .ge. 2     &
		     .and. mod((NJ(m)-1),NN) .eq. 0 .and. (NJ(m)-1)/NN .ge. 2 ) 
        Km=Km+1              ! 所允许的网格重数
	    NN=NN*2
       enddo
       Km_grid=min(Km_grid,Km)
    enddo
    Print*, " Finished check Mesh2d.dat,  Most stage is ", Km_grid
    print*,  "Check bc2d.in ..." 
    open(88,file="bc2d_1.in")
    read(88,*)
    read(88,*)
    read(88,*) Ntmp

    do m=1,NB
     read(88,*)
     read(88,*)
     read(88,*) Bsub    !number of the subface in the Block m
     read(88,*)
     do ksub=1, Bsub
      read(88,*)  f_no, face, ist, iend, jst, jend,  neighb, subface, orient
	     NN=1
		 Km=1
         do while( mod((ist-1),NN) .eq. 0 .and. mod((iend-1),NN) .eq.0       &
		     .and. mod((jst-1),NN) .eq. 0 .and. mod((jend-1),NN) .eq. 0  ) 
          NN=NN*2
		  Km=Km+1
		 enddo
        Km_grid=min(Km_grid,Km)
     enddo
     enddo


     close(88)


	 print*, "Check multigrid OK"
	 print*, "Total Block number is ", NB, "Total Cell number is " , N_Cell
	 print*, "Most stage number of multi-grid is ", Km_grid
	 print*, "--------------------------------------------------------------"
!-------------------------------------------------------
      Num_Mesh=min(3,Km_grid,Num_Mesh)    ! 设定网格重数，本版本最多允许3重网格
      print*, Num_Mesh, "  stage grids is used !"
!-------------------------------------------------------
	     
	 deallocate(NI,NJ)
    end subroutine check_mesh_multigrid

!--------------------------------------------------------------------------------------
!   创建数据结构： 最细网格 （储存几何量及守恒变量）
   subroutine Creat_Mesh1
   use Global_var
   implicit none
   integer,allocatable,dimension(:):: NI,NJ
   integer:: NB,m,nx,ny,i,j,k
   Type (Block_TYPE),pointer:: B
   real*8:: dx,dy

   print*, "-------------------------------------"
   print*, "read Mesh2d.dat"
   open(99,file="Mesh2d.dat")
   read(99,*) NB   ! Block number
     Mesh(1)%Num_Block=NB
     Mesh(1)%Num_Cell=0
     allocate(Mesh(1)%Block(NB))
     allocate( NI(NB),NJ(NB) )
     read(99,*) (NI(k), NJ(k), k=1,NB)
    do m=1,NB
     B => Mesh(1)%Block(m)
     B%Block_no=m
	 B%FVM_FDM=METHOD_FVM   ! 数值方法，默认为有限体积法
     B%nx=NI(m); B%ny=NJ(m) 
     nx=B%nx ; ny= B%ny
     Mesh(1)%Num_Cell=Mesh(1)%Num_Cell+(nx-1)*(ny-1)

!   申请内存   (x,y) 节点坐标； (x1,y1)网格中心坐标; s0 控制体体积； U, Un 守恒变量
   
     allocate(B%x(1-LAP:nx+LAP,1-LAP:ny+LAP), B%y(1-LAP:nx+LAP,1-LAP:ny+LAP))   
     allocate(B%x1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%y1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
     allocate(B%U(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))    !守恒变量 (4个流体变量；或6个变量：4个流体变量+k+w) 
     allocate(B%deltU(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))   !两步之差，由粗网格插值过来,多重网格使用；4个变量，k,w方程不采用 

	 allocate(B%Res(NVAR,1:nx-1,1:ny-1))        ! 残差
	 allocate(B%dt(1:nx-1,1:ny-1))              ! 时间步长
     
	 allocate(B%dU(NVAR,nx,ny))                 ! =U(n+1)-U(n), LU-SGS中使用
     allocate(B%Amu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), B%Amu_t(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))

!  几何量
     allocate(B%vol(nx,ny),B%si(nx,ny),B%sj(nx,ny),B%ni1(nx,ny),B%ni2(nx,ny),B%nj1(nx,ny),B%nj2(nx,ny)) 
     allocate(B%Lci(nx,ny),B%Lcj(nx,ny),B%Lvi(nx,ny),B%Lvj(nx,ny))        ! 谱半径 

     allocate(B%Un(NVAR,0:nx,0:ny))             ! 上一时间步的值
     if( Time_Method .eq. Time_dual_LU_SGS ) then
	  allocate(B%Un1(NVAR,0:nx,0:ny))    ! n-1时间步的值， 双时间步LU_SGS方法中采用
     endif


!   初始化	   
       B%x1(:,:)=0.d0; B%y1(:,:)=0.d0; B%vol(:,:)=0.d0                     
       B%U(1,:,:)=1.d0; B%U(2,:,:)=0.d0; B%U(3,:,:)=0.d0; B%U(4,:,:)=1.d0
       B%Un(1,:,:)=1.d0; B%Un(2,:,:)=0.d0; B%Un(3,:,:)=0.d0; B%Un(4,:,:)=0.d0
       B%Lci(:,:)=0.d0; B%Lcj(:,:)=0.d0; B%Lvi(:,:)=0.d0; B%Lvj(:,:)=0.d0
       B%dU(:,:,:)=0.d0
	   B%Amu(:,:)=0.d0;  B%Amu_t(:,:)=0.d0
       B%DeltU=0.d0
       B%Res=0.d0
       B%dt=0.d0
       if(Nvar .eq. 5) then
	     B%U(5,:,:)=1.d0/Re
       else if(NVar .eq. 6) then
	    B%U(5,:,:)=0.d0          ! 湍动能
        B%U(6,:,:)=1.d0          ! 湍能比耗散率
       endif
       
     read(99,*) ((B%x(i,j),i=1,nx),j=1,ny),  ((B%y(i,j),i=1,nx),j=1,ny)

!   计算控制体的体积 
!   控制体中心坐标需要等获得Ghost Cell 区信息后计算       
	   do j=1,ny-1
       do i=1,nx-1
         B%vol(i,j) =abs((B%x(i,j)-B%x(i+1,j+1))*(B%y(i+1,j)-B%y(i,j+1)) -  &
           (B%x(i+1,j)-B%x(i,j+1))*(B%y(i,j)-B%y(i+1,j+1)) )*0.5d0           ! 控制体体积（面积）
       enddo
	   enddo

!     几何量 （边长，法方向）      
	  do j=1,ny-1
      do i=1,nx
       dx=B%x(i,j+1)-B%x(i,j)
       dy=B%y(i,j+1)-B%y(i,j)
       B%si(i,j)=sqrt(dx*dx+dy*dy)
       B%ni1(i,j)=dy/B%si(i,j); B%ni2(i,j)=-dx/B%si(i,j)   ! normal vector at (i,j) or (I-1/2,J) 
	  enddo
	  enddo

      do j=1,ny
      do i=1,nx-1
       dx=B%x(i+1,j)-B%x(i,j)
       dy=B%y(i+1,j)-B%y(i,j)
       B%sj(i,j)=sqrt(dx*dx+dy*dy)     ! length 
       B%nj1(i,j)=-dy/B%sj(i,j); B%nj2(i,j)=dx/B%sj(i,j)     ! normal vector at i, j+1/2
      enddo
	  enddo
 
	
	 enddo
     close(99) 
	 deallocate(NI,NJ)

!    时间步、时间
       Mesh(1)%Kstep=0
       Mesh(1)%tt=0.d0

	 print*, "read Mesh2d.dat OK"
     end subroutine Creat_Mesh1


!----Mesh control message (bc2d.in)------------------------------------------
    subroutine read_bcin 
    use Global_Var
    implicit none
    integer:: NB,m,ksub
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc

    print*, "read bc2d.in ......"
    open(88,file="bc2d_1.in")
    read(88,*)
    read(88,*)
    read(88,*) NB
    if(NB .ne. Mesh(1)%Num_Block) then
      print*, "Error!  Block number in bc2d.in is not equal to that in Mesh2d.dat !"
      stop
    endif

    do m=1,NB
     B => Mesh(1)%Block(m)
     read(88,*)
     read(88,*)
     read(88,*) B%subface   !number of the subface in the Block m
     read(88,*)
     allocate(B%bc_msg(B%subface))
     do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     read(88,*)  Bc%f_no, Bc%face, Bc%ist, Bc%iend, Bc%jst, Bc%jend,  Bc%neighb, Bc%subface, Bc%orient
     enddo
     enddo
     close(88)
	 print*, "read bc2d.in OK"
   end  subroutine read_bcin

!-----------------------------------------------------------------------------------
! 根据上级网格信息，创建新网格m2 (稀疏网格)
! 稀疏网格不使用湍流模型，变量数为4 
   subroutine Creat_Mesh(m1,m2)
    use Global_Var
    implicit none
    integer:: NB,m,m1,m2,ksub,nx,ny,i,j,i1,j1,Bsub
    Type (Block_TYPE),pointer:: B1,B2
    TYPE (BC_MSG_TYPE),pointer:: Bc1,Bc2
    Type (Mesh_TYPE),pointer:: MP1,MP2
	real*8:: dx,dy
    print*,  "Creat Mesh ......", m2

    MP1=>Mesh(m1)             ! 上一级网格 （细网格）
	Mp2=>Mesh(m2)             ! 本级网格   （粗网格）

    NB=MP1%Num_Block
	MP2%Num_Block=NB      !  网格m2与m1 块数相同
    MP2%Mesh_no=m2        ! 网格号
    MP2%Num_Cell=0      

	allocate(MP2%Block(NB))   ! 在MP2中创建数据结构：“块”

	do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)
	 B2%Block_no=m
 	 B2%FVM_FDM=METHOD_FVM   ! 数值方法，粗网格只能使用有限体积法（高精度方法只在密网格上运行）
 
     nx=(B1%nx-1)/2+1          ! 粗网格的点数
	 ny=(B1%ny-1)/2+1
	 B2%nx=nx
	 B2%ny=ny
     
	 MP2%Num_Cell=MP2%Num_Cell+(nx-1)*(ny-1)          ! 统计MP2的总网格单元数

!    创建几何量及物理量
!--------------------------------------------------
     allocate(B2%x(1-LAP:nx+LAP,1-LAP:ny+LAP), B2%y(1-LAP:nx+LAP,1-LAP:ny+LAP))   
     allocate(B2%x1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B2%y1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
     allocate(B2%U(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))   !守恒变量  （不使用湍流模型，不含k,w, 总共4个变量）
     allocate(B2%deltU(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))  ! 两时间步U的差值 （多重网格使用，从粗网格插值而来）；
	 allocate(B2%Un(4,0:nx,0:ny))   ! 上一时间步的值
	 allocate(B2%Res(4,1:nx-1,1:ny-1))        ! 残差
	 allocate(B2%dt(1:nx-1,1:ny-1))           ! 时间步长
     allocate(B2%QF(4,1:nx-1,1:ny-1))         ! 强迫函数
	 allocate(B2%dU(4,nx,ny))   ! LU-SGS中使用
     allocate(B2%Amu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), B2%Amu_t(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
!    几何量
     allocate(B2%vol(nx,ny),B2%si(nx,ny),B2%sj(nx,ny),B2%ni1(nx,ny),B2%ni2(nx,ny),B2%nj1(nx,ny),B2%nj2(nx,ny)) 
     allocate(B2%Lci(nx,ny),B2%Lcj(nx,ny),B2%Lvi(nx,ny),B2%Lvj(nx,ny))        ! 谱半径 
    
	 if( Time_Method .eq. Time_dual_LU_SGS ) then
	  allocate(B2%Un1(NVAR,0:nx,0:ny))    ! n-1时间步的值， 双时间步LU_SGS方法中采用
     endif

!    初始化
       B2%x1(:,:)=0.d0; B2%y1(:,:)=0.d0; B2%vol(:,:)=0.d0                     
       B2%U(1,:,:)=1.d0; B2%U(2,:,:)=0.d0; B2%U(3,:,:)=0.d0; B2%U(4,:,:)=1.d0
       B2%Un(1,:,:)=1.d0; B2%Un(2,:,:)=0.d0; B2%Un(3,:,:)=0.d0; B2%Un(4,:,:)=1.d0
       B2%Res(:,:,:)=0.d0
	   B2%dt(:,:)=0.d0
	   B2%QF(:,:,:)=0.d0                       ! 强迫函数初始化为0
       B2%Lci(:,:)=0.d0; B2%Lcj(:,:)=0.d0; B2%Lvi(:,:)=0.d0; B2%Lvj(:,:)=0.d0
       B2%dU(:,:,:)=0.d0
	   B2%Amu(:,:)=0.d0;  B2%Amu_t(:,:)=0.d0
       
!   设定坐标信息（根据粗、细网格的对应关系）     
	  do j=1,ny
	  do i=1,nx
	    i1=2*i-1 ; j1=2*j-1
	    B2%x(i,j)=B1%x(i1,j1)         !粗网格与细网格的对应关系 （隔一个点设置一个粗网格点）
        B2%y(i,j)=B1%y(i1,j1)
	  enddo
	  enddo
!  计算网格面积
	   do j=1,ny-1
       do i=1,nx-1
         B2%vol(i,j) =abs((B2%x(i,j)-B2%x(i+1,j+1))*(B2%y(i+1,j)-B2%y(i,j+1)) -  &
           (B2%x(i+1,j)-B2%x(i,j+1))*(B2%y(i,j)-B2%y(i+1,j+1)) )*0.5d0           ! 控制体体积（面积）
       enddo
	   enddo

!     几何量 （边长，法方向）      
	  do j=1,ny-1
      do i=1,nx
       dx=B2%x(i,j+1)-B2%x(i,j)
       dy=B2%y(i,j+1)-B2%y(i,j)
       B2%si(i,j)=sqrt(dx*dx+dy*dy)
       B2%ni1(i,j)=dy/B2%si(i,j); B2%ni2(i,j)=-dx/B2%si(i,j)   ! normal vector at (i,j) or (I-1/2,J) 
	  enddo
	  enddo

      do j=1,ny
      do i=1,nx-1
       dx=B2%x(i+1,j)-B2%x(i,j)
       dy=B2%y(i+1,j)-B2%y(i,j)
       B2%sj(i,j)=sqrt(dx*dx+dy*dy)     ! length 
       B2%nj1(i,j)=-dy/B2%sj(i,j); B2%nj2(i,j)=dx/B2%sj(i,j)     ! normal vector at i, j+1/2
      enddo
	  enddo


!-----------------------------------------------  
!    创建连接信息
      Bsub=B1%subface        ! 子面数
      B2%subface=Bsub
      allocate(B2%bc_msg(Bsub))
      do ksub=1, Bsub
	   Bc1=> B1%bc_msg(ksub)    ! 上一级网格的连接信息
	   Bc2=> B2%bc_msg(ksub)    ! 本级网格的连接信息
       Bc2%f_no=Bc1%f_no
	   Bc2%face=Bc1%face
	   Bc2%ist=(Bc1%ist-1)/2+1   ! 粗、细网格下标的对应关系
	   Bc2%iend=(Bc1%iend-1)/2+1
	   Bc2%jst=(Bc1%jst-1)/2+1
	   Bc2%jend=(Bc1%jend-1)/2+1
       Bc2%neighb= Bc1%neighb
       Bc2%subface=Bc1%subface
       Bc2%orient=Bc1%orient
      enddo

  
     enddo
       Mesh(m2)%Kstep=0
       Mesh(m2)%tt=0.d0
	  Print*, "Creat Mesh ",m2," OK ", "Total Cell number is", MP2%Num_Cell 
   end  subroutine Creat_Mesh
!--------------------------------------------------------------------------------
! 初始化流场 
    subroutine Init_flow 
     use Global_var
	 implicit none
     integer:: i,j,m,m1
     Type (Block_TYPE),pointer:: B
     Type (Mesh_TYPE),pointer:: MP
    
	 if(Iflag_init .eq. 0) then
	   call init_flow_zero                      ! 从零流场算起 （先从粗网格计算，再插值到细网格）
     else
	   call init_flow_read                      ! 从文件读取流场
     endif
    
	 if(Time_Method .eq. Time_Dual_LU_SGS) then
	    MP=>Mesh(1)             
       do m=1,MP%Num_Block
     
		B => MP%Block(m)    ! Mesh(Num_Mesh) 是最粗的网格
         do j=1,B%ny-1
         do i=1,B%nx-1
         do m1=1,NVAR
		  B%Un(m1,i,j)=B%U(m1,i,j)
		  B%Un1(m1,i,j)=B%U(m1,i,j) 
         enddo
         enddo
         enddo
		enddo
     endif 
	
	
	print*, " Initialize OK ......"
   end subroutine Init_flow







! 用来流初始化； 多重网格情况下，从最粗网格开始计算（然后插值到细网格）  
     subroutine init_flow_zero
     use Global_var
	 implicit none
     real*8:: d0,u0,v0,p0,T0,tmp
     integer:: i,j,m,step,nMesh
     Type (Block_TYPE),pointer:: B
     Type (Mesh_TYPE),pointer:: MP
	 
!-------------------------------------------------------------------
	  d0=1.d0; u0=1.d0*cos(AoA); v0=1.d0*sin(AoA) ; p0=1.d0/(gamma*Ma*Ma)
      if(Nvar .eq. 5) then
        MP=>Mesh(1)             ! 密网格
        do m=1,MP%Num_Block
         B => MP%Block(m)    ! Mesh(Num_Mesh) 是最粗的网格
         do j=1,B%ny-1
         do i=1,B%nx-1
          B%U(5,i,j)=vt_inf/Re   ! Vt值 （通常设定为层流粘性系数的3-5倍）
         enddo
         enddo
        enddo
      endif
      
!  k和w的初值 （0初值，仅在最密网格上有值）
      if(Nvar .eq. 6) then
        MP=>Mesh(1)             ! 密网格
        do m=1,MP%Num_Block
         B => MP%Block(m)    ! Mesh(Num_Mesh) 是最粗的网格
         do j=1,B%ny-1
         do i=1,B%nx-1
          B%U(5,i,j)=10.d0*Kt_inf     ! 湍动能， 初值为来流值的10倍（湍流易于发展）
          B%U(6,i,j)=Wt_inf    ! 湍能比耗散率，初值为来流值 
         enddo
         enddo
        enddo
      endif
!---------------------------------------------------------------------
  !  流体量的初值，由粗网格逐级计算插值而来
      
      MP=>Mesh(Num_Mesh)             ! 最稀疏的网格
      do m=1,MP%Num_Block
      B => MP%Block(m)    ! Mesh(Num_Mesh) 是最粗的网格
      do j=1,B%ny-1
      do i=1,B%nx-1
       B%U(1,i,j)=d0
       B%U(2,i,j)=d0*u0
       B%U(3,i,j)=d0*v0
       B%U(4,i,j)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0)
      enddo
      enddo
      enddo
    call Boundary_condition_onemesh(Num_Mesh)                   ! 边界条件 （设定Ghost Cell的值）
    call update_buffer_onemesh(Num_Mesh)                        ! 同步各块的交界区
!-----------------------------------------------------------------
!------------------------------------------------------
!   准备初值的过程
!   从最粗网格计算，逐级插值到细网格
  

  do nMesh=Num_Mesh,1,-1 
   do step=1, Pre_Step_Mesh(nMesh)      
    call NS_Time_advance(nMesh)
    if(mod(step,Kstep_show) .eq. 0) call output_Res(nMesh)
   enddo
!    call output (nMesh)
    if(nMesh .gt. 1) then
    call prolong_U(nMesh,nMesh-1,1)                            ! 把nMesh重网格上的物理量插值到上一重网格; flag=1 插值U本身
    call Boundary_condition_onemesh(nMesh-1)                   ! 边界条件 （设定Ghost Cell的值）
    call update_buffer_onemesh(nMesh-1)                        ! 同步各块的交界区
   
    print*, " Prolong  to mesh ", nMesh-1, "   OK"           
    endif
  enddo
  end subroutine init_flow_zero 
 
!-----------------------------------------------
!  从flow2d.dat 文件读取初值 （最密的网格）， 数据为tecplot格式
     subroutine init_flow_read
     use Global_var
     real*8:: d0,u0,v0,p0,T0,kt0,kw0,tmp
     Type (Mesh_TYPE),pointer:: MP

     integer:: i,j,step,nMesh
     Type (Block_TYPE),pointer:: B
      MP=>Mesh(1)             ! 最密的网格
      print*, "Init from 'flow2d.dat' ......"

       open(99,file="flow2d.dat")
	   read(99,*)
       do m=1,MP%Num_Block
       B => MP%Block(m)                 ! 网格块
         read(99,*)
         do j=0,B%ny
	     do i=0,B%nx
	     read(99,*) tmp,tmp,d0,u0,v0,T0
         B%U(1,i,j)=d0                    ! 密度
         B%U(2,i,j)=d0*u0                 ! x-方向动量密度
         B%U(3,i,j)=d0*v0                 ! y-方向动量密度
         B%U(4,i,j)=d0*Cv*T0+0.5d0*d0*(u0*u0+v0*v0)  ! 总能量密度
         enddo
         enddo
       enddo
	   close(99)
       
	   if(MP%NVAR .eq. 5) then
       open(99,file="SA2d.dat")   ! 读取vt
	   read(99,*)
       do m=1,MP%Num_Block
       B => MP%Block(m)                 
         read(99,*)
         do j=0,B%ny
	     do i=0,B%nx
	     read(99,*) tmp,tmp,B%U(5,i,j)
         enddo
         enddo
       enddo
	   close(99)
       endif

	   if(MP%NVAR .eq. 5) then
       open(99,file="SST2d.dat")   ! 读取Kt,Wt
	   read(99,*)
       do m=1,MP%Num_Block
       B => MP%Block(m)                 
         read(99,*)
         do j=0,B%ny
	     do i=0,B%nx
	     read(99,*) tmp,tmp,B%U(5,i,j),B%U(6,i,j)
         enddo
         enddo
       enddo
	   close(99)
       endif



      call Boundary_condition_onemesh(1)                   ! 物理边界条件 
      call update_buffer_onemesh(1)                        ! 内边界条件
    end  subroutine init_flow_read