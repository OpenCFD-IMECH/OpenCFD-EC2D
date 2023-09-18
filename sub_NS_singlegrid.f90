!----------------------------------------------------------------------
! Copyright by Li Xinliang,  Code by Li Xinliang
!
!----------------------------------------------------------------------
! 在给定的网格上求解N-S方程 （推进1个时间步）
! 对于单重网格，nMesh=1;  对于多重网格，nMesh=1,2,3, ... 分别对应用细网格、粗网格、更粗网格 ...
 
  subroutine NS_Time_advance(nMesh)
   use Global_var
   implicit none
   integer:: nMesh
   if(Time_Method .eq. Time_Euler1) then
      call NS_Time_advance_1Euler(nMesh)                   ! 1阶Euler
   else if (Time_Method .eq. Time_RK3) then
     call NS_Time_advance_RK3(nMesh)                       ! 3阶RK
   else if (Time_Method .eq. Time_LU_SGS) then             ! LU-SGS隐格式
     call NS_time_advance_LU_SGS(nMesh)
   else if(Time_Method .eq. Time_Dual_LU_SGS) then         ! 双时间步LU-SGS (不支持多重网格)
     call Dual_time_LU_SGS
   else
      print*, "This time advance method is not supported!!!"
   endif
   !  强制 k,w, vt 非负
    call force_vt_kw(nMesh)

   end subroutine NS_Time_advance

!---------------------------------------------------------------------------------------------
! 强制vt, k,w非负   
   subroutine force_vt_kw(nMesh)
    use Global_var
    implicit none
    integer:: nMesh,mBlock,nx,ny,i,j
    Type (Block_TYPE),pointer:: B
    Type (Mesh_TYPE),pointer:: MP
  
    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
    if(MP%NVAR == 5) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m)
 	  do j=1,ny-1
      do i=1,nx-1
       if(B%U(5,i,j) < 0)  B%U(5,i,j)=1.d-10
      enddo
      enddo
!$OMP END PARALLEL DO
    else if (MP%NVAR == 6) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m)
 	  do j=1,ny-1
      do i=1,nx-1
       if(B%U(5,i,j) < 0)  B%U(5,i,j)=1.d-10
       if(B%U(6,i,j) < 0)  B%U(6,i,j)=1.d-10
	  enddo
      enddo
!$OMP END PARALLEL DO
    endif
   enddo    
  end
!--------------------------------------------------------------------------------------




! 采用1阶Euler法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
 subroutine NS_Time_advance_1Euler(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Block_TYPE),pointer:: B
   Type (Mesh_TYPE),pointer:: MP
   real*8:: du
 
    MP=>Mesh(nMesh)
    call Comput_Residual_one_mesh(nMesh)     ! 单重网格上计算残差
    if(nMesh .ne. 1) call Add_force_function(nMesh)   !  添加强迫函数（多重网格的粗网格使用）

    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   时间推进 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m,du)
      do j=1,ny-1
      do i=1,nx-1
      do m=1,MP%NVAR
        du=B%Res(m,i,j)/B%vol(i,j)    
        B%U(m,i,j)=B%U(m,i,j)+B%dt(i,j)*du
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO 

    enddo    

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------  
      call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
      call update_buffer_onemesh(nMesh)              ! 同步各块的交界区


     Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! 时间 （使用全局时间步长法时有意义）
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! 计算步数
!    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 end subroutine NS_Time_advance_1Euler






!----------------------------------------------------------------------------------------
! 采用3阶RK方法推进1个时间步 （第nMesh重网格 的单重网格）
 subroutine NS_Time_advance_RK3(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,KRK,i,j,m,nx,ny
   Type (Block_TYPE),pointer:: B
   Type (Mesh_TYPE),pointer:: MP
   real*8:: du
   real*8:: Ralfa(3), Rbeta(3) , Rgamma(3)  

	 Ralfa(1)=1.d0 ;  Ralfa(2)=3.d0/4.d0 ; Ralfa(3)=1.d0/3.d0
     Rbeta(1)=1.d0 ;  Rbeta(2)=1.d0/4.d0 ; Rbeta(3)=2.d0/3.d0
     Rgamma(1)=0.d0;  Rgamma(2)=1.d0/4.d0; Rgamma(3)=2.d0/3.d0

      MP=>Mesh(nMesh)
      do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m)
       do j=0,B%ny
	   do i=0,B%nx
	   do m=1,MP%NVAR
	   B%Un(m,i,j)=B%U(m,i,j)
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO
	  enddo

 
   do KRK=1,3                                         ! 3-step Runge-Kutta Method
      
	  call Comput_Residual_one_mesh(nMesh)              ! 计算残差
      if(nMesh .ne. 1) call Add_force_function(nMesh)   ! 添加强迫函数（多重网格的粗网格使用）
   
	  do mBlock=1,Mesh(nMesh)%Num_Block
      B => Mesh(nMesh)%Block(mBlock)                  ! 第nMesh 重网格的第mBlock块
      nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!    时间推进 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m,du)
      do j=1,ny-1
      do i=1,nx-1
      do m=1,MP%NVAR
		du=B%Res(m,i,j)/B%vol(i,j)  
        B%U(m,i,j)=Ralfa(KRK)*B%Un(m,i,j)+Rgamma(KRK)*B%U(m,i,j)+B%dt(i,j)*Rbeta(KRK)*du        ! 3阶RK
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO

    enddo    
!---------------------------------------------------------------------------------------

      call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
      call update_buffer_onemesh(nMesh)              ! 同步各块的交界区

  enddo
     
	 Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! 时间 （使用全局时间步长法时有意义）
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1      ! 计算步数
!    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 end subroutine NS_Time_advance_RK3






!---------------------------------------------------------------------------------------------
! 采用LU_SGS法进行时间推进一个时间步 （第nMesh重网格 的单重网格）
 subroutine NS_Time_advance_LU_SGS(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real*8:: alfa1
   
    call Comput_Residual_one_mesh(nMesh)     ! 单重网格上计算残差
    alfa1=0.d0                               
    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock,alfa1)                          ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   时间推进 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m)
      do j=1,ny-1
      do i=1,nx-1
      do m=1,MP%NVAR
        B%U(m,i,j)=B%U(m,i,j)+B%du(m,i,j)                      ! U(n+1)=U(n)+dU
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

    enddo    

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------  
      call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
      call update_buffer_onemesh(nMesh)              ! 同步各块的交界区


     Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! 时间 （使用全局时间步长法时有意义）
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! 计算步数
!    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 end subroutine NS_Time_advance_LU_SGS


!------------------------------------------------------------------------------------
! 双时间步 LU_SGS方法
! 采用Dual time step LU_SGS法进行时间推进一个时间步 
! 目前只支持单重网格；
 subroutine Dual_time_LU_SGS
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny,kt_in
   real*8:: max_res
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real*8:: alfa1
   integer,save:: Iflag=0

       if(Iflag .eq. 0) then
	    Iflag=1                     
		alfa1=1.d0/dt_global         ! 仅针对第1个时间步
	   else
	    alfa1=3.d0/(2.d0*dt_global)
	   endif
   
   nMesh=1
  
  do kt_in=1, Nstep_inner_Limit    ! 内循环迭代
   
    call Comput_Residual_one_mesh(nMesh)     ! 单重网格上计算残差
    MP=>Mesh(nMesh)

    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock,alfa1)                          ! 采用LU_SGS方法计算DU=U(n+1)-U(n)
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m)
      do j=1,ny-1
      do i=1,nx-1
      do m=1,MP%NVAR
        B%U(m,i,j)=B%U(m,i,j)+B%du(m,i,j)                      ! U(n+1)=U(n)+dU
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo    

      call Boundary_condition_onemesh(nMesh)         ! 边界条件 （设定Ghost Cell的值）
      call update_buffer_onemesh(nMesh)              ! 同步各块的交界区

     max_res=MP%Res_rms(1)       ! 最大均方根残差
     do m=1,NVAR
	 max_res=max(max_res,MP%Res_rms(m))
	 enddo
     if( max_res .le. Res_Inner_Limit) exit   ! 达到残差标准，跳出内迭代
  enddo
 	
	 print*, "Inner step ... ", kt_in
	 print*, "rms residual eq =", MP%Res_rms(1:NVAR)


! 内迭代结束，更新U(n),U(n-1)    
    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m)
      do j=1,ny-1
      do i=1,nx-1
      do m=1,MP%NVAR
        B%Un1(m,i,j)=B%Un(m,i,j)                      ! U(n-1) <= U(n)
        B%Un(m,i,j)=B%U(m,i,j)                     
	  enddo
      enddo
      enddo
!$OMP END PARALLEL DO
    enddo    



     Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! 时间 （使用全局时间步长法时有意义）
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! 计算步数

     print*, "----------------------------------------------------"
     print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 
 end subroutine Dual_time_LU_SGS












! 计算（当地）时间步长
   subroutine comput_dt(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real*8,parameter:: C0=1
     Type (Block_TYPE),pointer:: B
    
     B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
     do j=1,ny-1
     do i=1,nx-1
      if(Iflag_local_dt .eq. 0) then   ! 全局时间步长
        B%dt(i,j)=dt_global                                 
      else                             ! 当地时间步长
         B%dt(i,j)=CFL*B%vol(i,j)/(B%Lci(i,j)+B%Lcj(i,j) +C0*(B%Lvi(i,j)+B%Lvj(i,j)) )     ! 局部时间步长
         if(B%dt(i,j) .gt. dtmax) B%dt(i,j)=dtmax
         if(B%dt(i,j) .lt. dtmin) B%dt(i,j)=dtmin
      endif
	 enddo
	 enddo
!$OMP END PARALLEL DO 
  end subroutine comput_dt

!----------------------------------------------------------
! 计算无粘性及粘性项的谱半径，在加速收敛技术（局部时间步长，隐格式，残差光顺等）中使用
   subroutine comput_Lvc(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real*8:: C0,si,sj,s0,ni1,ni2,nj1,nj2,uni,unj,tmp1
     Type (Block_TYPE),pointer:: B
    
	 C0=1.d0
     B => Mesh(nMesh)%Block(mBlock)                 !第nMesh 重网格的第mBlock块
     nx=B%nx; ny=B%ny

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,si,sj,s0,ni1,ni2,nj1,nj2,uni,unj,tmp1)

     do j=1,ny-1
     do i=1,nx-1
 		 s0 =B%vol(i,j)         ! 面积
		 si=0.5d0*(B%si(i,j)+B%si(i+1,j))         ! 控制体边长（两侧平均）
		 ni1=0.5d0*(B%ni1(i,j)+B%ni1(i+1,j))      ! 界面法方向（两侧平均）
		 ni2=0.5d0*(B%ni2(i,j)+B%ni2(i+1,j))
         sj=0.5d0*(B%sj(i,j)+B%sj(i,j+1))
		 nj1=0.5d0*(B%nj1(i,j)+B%nj1(i,j+1))
		 nj2=0.5d0*(B%nj2(i,j)+B%nj2(i,j+1))
         uni=uu(i,j)*ni1+v(i,j)*ni2             !法向速度
		 unj=uu(i,j)*nj1+v(i,j)*nj2

!        谱半径
		 B%Lci(i,j)=(abs(uni)+cc(i,j))*si          ! 无粘项Jocabian矩阵的谱半径 （Blazek's Book 6.1.4节）
		 B%Lcj(i,j)=(abs(unj)+cc(i,j))*sj

		 tmp1=gamma/d(i,j)*(B%Amu(i,j)/Pr+B%Amu_t(i,j)/PrT)
		 B%Lvi(i,j)=tmp1*si*si/s0         ! 粘性项Jocabian矩阵谱半径
         B%Lvj(i,j)=tmp1*sj*sj/s0
 	 enddo
	 enddo
!$OMP END PARALLEL DO
  end subroutine comput_Lvc
