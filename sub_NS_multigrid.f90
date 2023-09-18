!----------------------------------------------------------------------
! Copyright by Li Xinliang, Code by Li Xinliang

!----------------------------------------------------------------------
! 多重网格求解N-S方程 （推进1个时间步）
!  nMesh=1,2,3 分别对应用细网格、粗网格、更粗网格
! 包括2重网格和3重网格两个子程序；
!---------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! 两重网格上推进1个时间步 (3阶RK or 1th Euler)
  subroutine NS_2stge_multigrid
   use Global_var
   implicit none
   integer::nMesh,m
   Type (Block_TYPE),pointer:: B
   integer,parameter:: Time_step_coarse_mesh=3     ! 粗网格迭代步数

!---------------------------------------------------
! ----  网格1 -----------------
      if(Time_Method .eq. Time_Euler1) then
	     call  NS_Time_advance_1Euler(1)                 ! 细网格，1阶Euler方法推进1步 -> U(n+1)
      else
	     call  NS_Time_advance_RK3(1)                    ! 细网格，RK方法推进1步 -> U(n+1)
      endif

! ???????????     
	  call  Comput_Residual_one_mesh(1)           ! 计算网格1的残差 R(n+1)   ! ???????? 似乎可以取消该步
	  
	  call  interpolation2h(1,2,2)                ! 把残差插值到网格2 (储存在QF里面)
	  call  interpolation2h(1,2,1)                ! 把守恒变量从网格1插值到网格2   （flag=1 插值守恒变量，=2 插值残差）

!------------------------------
      call  Boundary_condition_onemesh(2)         ! 物理边界条件
      call  update_buffer_onemesh(2)              ! 内边界条件
      call  Comput_Residual_one_mesh(2)           ! 计算网格2的残差
      call  comput_force_function(2)              ! 计算强迫函数QF
	  
      if(Time_Method .eq. Time_Euler1) then
	    call Set_Un(2)                               ! 记录初始值  （RK方法中已经包含了该步） 
        do m=1, Time_step_coarse_mesh
	     call  NS_Time_advance_1Euler(2)             ! 1阶Euler迭代若干步
		enddo
	   else 
		 call  NS_Time_advance_RK3(2)                ! RK方法推进1步 （网格2）
       endif

	  call  comput_delt_U(2)                      ! 计算修正量deltU （储存在Un里面）
      call  prolong_U(2,1,2)                      ! 把修正量插值到细网格 (储存在Un里面); flag=2 插值deltU (储存在Un里)
!------------------------------------	 
	  call  comput_new_U(1)                       ! 计算新的U  (U=U+deltU)
      call  Boundary_condition_onemesh(1)         ! 物理边界条件
      call  update_buffer_onemesh(1)              ! 内边界条件

!     print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
  end subroutine NS_2stge_multigrid

!------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! 三重网格上迭代1个时间步 （V-型迭代） 3阶RK or 1阶Euler
  subroutine NS_3stge_multigrid
   use Global_var
   implicit none
   integer::nMesh,m
   integer,parameter:: Time_step_coarse_mesh=3     ! 粗网格迭代步数 (对1阶Euler有效)

!------------------
       Type (Block_TYPE),pointer:: B
       integer:: i,j
!------------------------

!---------------------------------------------------
! ----  网格1 -----------------
      if(Time_Method .eq. Time_Euler1) then
	     call  NS_Time_advance_1Euler(1)                 ! 细网格，1阶Euler方法推进1步 -> U(n+1)
      else
	     call  NS_Time_advance_RK3(1)                    ! 细网格，RK方法推进1步 -> U(n+1)
      endif

! ?????????????????????????????
	  call  Comput_Residual_one_mesh(1)            ! 计算网格1的残差 R(n+1)   !    
! -----------------------------  
	  
	  call  interpolation2h(1,2,2)                 ! 把残差插值到网格2 (储存在网格2的QF里面)
	  call  interpolation2h(1,2,1)                 ! 把守恒变量从网格1插值到网格2 （储存到U里面）  （flag=1 插值守恒变量，=2 插值残差）
      
!-------网格2 --------------------
      call  Boundary_condition_onemesh(2)          ! 物理边界条件
      call  update_buffer_onemesh(2)               ! 内边界条件
      call  Comput_Residual_one_mesh(2)            ! 计算网格2的残差         Res_2h(0)
      call  comput_force_function(2)               ! 计算强迫函数QF （网格2）QF_2h=QF_2h-Res_2h(0)
	  
      if(Time_Method .eq. Time_Euler1) then
	    call Set_Un(2)                               ! 记录初始值  （RK方法中已经包含了该步） 
        do m=1, Time_step_coarse_mesh
	     call  NS_Time_advance_1Euler(2)             ! 1阶Euler迭代若干步
		enddo
	   else 
		 call  NS_Time_advance_RK3(2)                ! RK方法推进1步 （网格2）
       endif

      call  Comput_Residual_one_mesh(2)             ! 计算网格2的残差 R_2h(n+1)  
      call  Add_force_function(2)                   ! 添加上强迫残差(储存在Res里面)  RF_2h(n+1)=R_2h(n+1)+QF_2h  ;  目的：插值到网格3上
	  call  interpolation2h(2,3,2)                  ! 把残差插值到网格3 (储存在网格3的QF里面)

	  call  interpolation2h(2,3,1)                 ! 把守恒变量从网格2插值到网格3 （储存到U里面）  （flag=1 插值守恒变量，=2 插值残差）

!------网格3----------------------	  
      call  Boundary_condition_onemesh(3)          ! 边界条件: 物理边界 
      call  update_buffer_onemesh(3)               ! 内边界
	  call  Comput_Residual_one_mesh(3)           ! 计算网格3的残差
	
      call  comput_force_function(3)              ! 计算强迫函数QF （网格3）
	  
      if(Time_Method .eq. Time_Euler1) then
	    call Set_Un(3)
        do m=1, Time_step_coarse_mesh
	     call  NS_Time_advance_1Euler(3)              ! 1阶Euler迭代若干步
		enddo
	   else 
		 call  NS_Time_advance_RK3(3)                 ! RK方法推进1步 （网格3）
       endif
	  
	  call  comput_delt_U(3)                      ! 计算修正量deltU (=U-Un)
	  call  prolong_U(3,2,2)                      ! 把修正量插值到网格2 (储存在deltU里面); flag=2 插值deltU 

!------网格2------------------------      
	  call  comput_new_U(2)                       ! 网格2计算新的U  (U=U+deltU)
	  call  comput_delt_U(2)                      ! 计算修正量deltU =U-Un
      call  prolong_U(2,1,2)                      ! 把修正量插值到细网格 (储存在deltU里面); flag=2 插值deltU 
!------网格1------------------------------
	  call  comput_new_U(1)                       ! 计算新的U  (U=U+deltU)
      call Boundary_condition_onemesh(1)          ! 物理边界条件
      call update_buffer_onemesh(1)               ! 内边界条件

!     print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
  end subroutine NS_3stge_multigrid

!------------------------------------------------------------------------------------------








!------------------------------------------------------------------------------------------
  
!  计算强迫函数 QF=Ih_to_2h Res(n-1) - Res(n)      ! QF中储存着细网格插值过来的残差
  subroutine comput_force_function(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m
   Type (Block_TYPE),pointer:: B
       do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
       do j=1,B%ny-1
	   do i=1,B%nx-1
	   do m=1,4
	   B%QF(m,i,j)=B%QF(m,i,j)-B%Res(m,i,j)            ! QF原先储存着从细网格插值过来的残差
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO
	   enddo
  end  subroutine comput_force_function

!------------------------------------------------------------
!  把强迫函数添加到残差中 RF=R+QF        
  subroutine Add_force_function(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m
   Type (Block_TYPE),pointer:: B
       do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
       do j=1,B%ny-1
	   do i=1,B%nx-1
	   do m=1,4
	   B%Res(m,i,j)=B%Res(m,i,j)+B%QF(m,i,j)            ! 添加强迫函数后的残差仍储存在B%Res里面 （节省内存）
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO

	   enddo
  end  subroutine Add_force_function

!----------------------------------------------------------------------  
!  计算修正量 deltU=U-Un
  subroutine comput_delt_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m
   Type (Block_TYPE),pointer:: B
       do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
       do j=0,B%ny
	   do i=0,B%nx
	   do m=1,4
	   B%deltU(m,i,j)=B%U(m,i,j)-B%Un(m,i,j)
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO
	   enddo
  end  subroutine comput_delt_U
!-----------------------------------------------------------------------
! 设定Un=U
  subroutine Set_Un(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m
   Type (Block_TYPE),pointer:: B
       do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
       do j=0,B%ny
	   do i=0,B%nx
	   do m=1,4
	   B%Un(m,i,j)=B%U(m,i,j)
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO

	   enddo
  end  subroutine Set_Un

!-------------------------------修正U --------------------------------------
  subroutine comput_new_U(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m
   Type (Block_TYPE),pointer:: B,B1
       do mBlock=1,Mesh(nMesh)%Num_Block
       B => Mesh(nMesh)%Block(mBlock)
       B1 => Mesh(nMesh)%Block(mBlock)
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
	   do j=1,B%ny-1
	   do i=1,B%nx-1
	   do m=1,4
!	   B%U(m,i,j)=B%U(m,i,j)+B %deltU(m,i,j)         !deltU里面储存的是U的修正量 （从粗网格插值而来）
 	   B%U(m,i,j)=B%U(m,i,j)+B1%deltU(m,i,j)           ! 修改2011-7-6 ,修改后使用Intel Fortran 32位编译器不发散 （尚未找到原因，可能还有其他问题!!!）
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO

	   enddo
  end  subroutine comput_new_U
!---------------------------------------------------------------------------


!----------------------------------------------------------------------
! 粗网格向细网格的插值(Prolong) 及 细网格向粗网格上插值 (interpolation)
!----------------------------------------------------------------------
! 将网格m1的守恒变量(U) 或U的差插值到网格m2 (上一级细网格)
! flag=1时，将U插值到上一级网格；  (准备初值时使用)
! flag=2时，将deltU插值到上一级网格 （deltU储存着本时间步与上个时间步U的差） 

   Subroutine prolong_U(m1,m2,flag)
   use Global_Var
   implicit none
   integer:: m1,m2,mb,flag
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   integer,allocatable,dimension(:,:):: ia,ja
   real*8:: a1=9.d0/16.d0,a2=3.d0/16.d0,a3=1.d0/16.d0   ! 插值系数
   integer:: i,j,m,nx1,ny1,nx2,ny2

    if(m1 .le. 1 .or. m1-m2 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2)
   do mb=1,MP1%Num_Block
     B1=>MP1%Block(mb)
	 B2=>Mp2%Block(mb)
     nx1=B1%nx; ny1=B1%ny
     nx2=B2%nx; ny2=B2%ny
	 allocate( ia(2,0:nx2),ja(2,0:ny2))

!    寻找插值基架点的下标 
!    ia(1,i) 是距离i点最近的粗网格点的下标；ia(2,i)是次近点的下标	 
	 do i=0,nx2
	  if(mod(i,2).eq.0) then
	   ia(1,i)=i/2                    !最近点
	   ia(2,i)=i/2+1                  !次近点
	  else  
	   ia(1,i)=i/2+1                  !最近点
	   ia(2,i)=i/2                    !次近点
	  endif
     enddo

    do j=0,ny2
	 if( mod(j,2).eq. 0) then
	  ja(1,j)=j/2
	  ja(2,j)=j/2+1
	 else
	  ja(1,j)=j/2+1
	  ja(2,j)=j/2
	 endif
	enddo

	if(flag .eq. 1) then
!     插值守恒变量U

!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
 	 do j=0,ny2
	 do i=0,nx2
	 do m=1,4
!               插值，最近点的权重a1, 次近点的权重a2, 最远点的权重a3	 
	  B2%U(m,i,j)=a1*B1%U(m,ia(1,i),ja(1,j))+a2*(B1%U(m,ia(2,i),ja(1,j))+B1%U(m,ia(1,i),ja(2,j)) ) &
	                                   +a3*B1%U(m,ia(2,i),ja(2,j))
     enddo
	 enddo
	 enddo
!$OMP END PARALLEL DO
	else
!  插值deltU
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
 	 do j=0,ny2
	 do i=0,nx2
	 do m=1,4
!               插值，最近点的权重a1, 次近点的权重a2, 最远点的权重a3	 
	  B2%deltU(m,i,j)=a1*B1%deltU(m,ia(1,i),ja(1,j))+a2*(B1%deltU(m,ia(2,i),ja(1,j))+B1%deltU(m,ia(1,i),ja(2,j)) ) &
	                                   +a3*B1%deltU(m,ia(2,i),ja(2,j))
     enddo
	 enddo
	 enddo
!$OMP END PARALLEL DO
    endif
    deallocate(ia,ja)
   enddo
   end Subroutine prolong_U

!-------------------------------------------------------------

!-----
! 将网格m1的守恒变量U插值到网格m2 (细网格->粗网格) 
   Subroutine interpolation2h(m1,m2,flag)
   use Global_Var
   implicit none
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   real*8,dimension(:,:,:),pointer:: P1,P2
   integer:: flag,m1,m2,mb,i,j,m,i1,i2,j1,j2
!   flag==1 插值守恒变量； flag==2 插值残差
	 if( m2-m1 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2)
   
   do mb=1,MP1%Num_Block
     B1=>MP1%Block(mb)
  	 B2=>Mp2%Block(mb)
    if(flag .eq. 1) then  ! 插值守恒变量
	  P1=>B1%U
	  P2=>B2%U
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,i1,j1,i2,j2,m)
      do j=1,B2%ny-1
	  do i=1,B2%nx-1
	   i1=2*i-1 ; i2=2*i
	   j1=2*j-1 ; j2=2*j
      do m=1,4
!     以控制体体积为权重的加权平均 	 
	    P2(m,i,j)=(P1(m,i1,j1)*B1%vol(i1,j1)+P1(m,i1,j2)*B1%vol(i1,j2)   &
	      +P1(m,i2,j1)*B1%vol(i2,j1)+P1(m,i2,j2)*B1%vol(i2,j2))/B2%vol(i,j)
  	  enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO

    else     ! 插值残差  （把m1网格上的残差B%Res 插值到m2网格上B%QF (然后减去本m2网格上的残差，形成强迫函数)）
   	  P1=>B1%Res
	  P2=>B2%QF
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,i1,j1,i2,j2,m)
	  do j=1,B2%ny-1
	  do i=1,B2%nx-1
	   i1=2*i-1 ; i2=2*i
	   j1=2*j-1 ; j2=2*j
       do m=1,4
	     P2(m,i,j)=P1(m,i1,j1)+P1(m,i1,j2) +P1(m,i2,j1)+P1(m,i2,j2)    ! 残差的插值： 简单相加
  	   enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO
 	endif
   enddo

  end Subroutine interpolation2h
