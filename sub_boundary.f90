! modules for Boundary layer condition
! For OpenCFD-EC 2D ver 1.1
! Copyright by Li Xinliang, lixl@imech.ac.cn
! Code by Li Xinliang, 2011-3-27 
! -------modified-----------------------------------------------------
! 2011-11-23: Symmetry boundary condition is adding
! 2011-12-9:  Isothormal wall boundary condition is adding
!---------------------------------------------------------------------
! 处理边界条件（非内边界） （处理一套网格）
     subroutine Boundary_condition_onemesh(nMesh)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     integer:: nMesh,mBlock,ksub

 ! ------------------------------------------------------ 
 do mBlock=1,Mesh(nMesh)%Num_Block
  B => Mesh(nMesh)%Block(mBlock)
  do  ksub=1,B%subface
  Bc=> B%bc_msg(ksub)

   if(Bc%neighb .lt. 0 ) then               ! 非内边界
      if( Bc%neighb .eq. BC_Wall ) then
          call boundary_wall(nMesh,mBlock,ksub)         ! 固壁边界
      else if( Bc%neighb .eq. BC_Farfield  .or. Bc%neighb .eq. BC_Inlet .or. Bc%neighb .eq. BC_Outlet) then
          call boundary_Farfield_inlet_outlet(nMesh,mBlock,ksub, Bc%neighb)     ! 远场边界 or 入口 or 出口
      else if( Bc%neighb .eq. BC_Symmetry_or_slidewall ) then
          call Symmetry_or_slidewall(nMesh,mBlock,ksub)   ! 对称（或滑移）边界
      endif
   endif
   enddo
   enddo
  end subroutine Boundary_condition_onemesh

!------------------------------------------------------------



!---------------------------------------------------------------------------
! inner boundary condition  内边界
! 根据网格连接关系，更新缓冲区内物理量的信息 (处理一套网格)
! 缓冲区为LAP层网格 (目前版本设定LAP=2)

     subroutine update_buffer_onemesh(nMesh)
     use Global_Var
     implicit none
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc,Bc1
     real*8,allocatable:: Ua(:,:,:)
     integer:: nMesh,Nvar1,i,j,k,k1,m,mBlock,ksub,n1,m_neighbour,msub,orient,Kflag_initial,nx1,ny1
     integer:: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2

!----------------------------------------------------------
 MP=>Mesh(nMesh)   ! 网格 （nMesh=3,2,1代表 粗、中、密网格） 
 Nvar1=MP%Nvar     ! 变量数目 （4或6个），粗网格只有4个变量； 密网格可以有4或6个变量（包括SST模型的两个变量）
 do mBlock=1,MP%Num_Block
  B => MP%Block(mBlock)
  do  ksub=1,B%subface
  
    Bc=> B%bc_msg(ksub)
    if( Bc%neighb .ge. 0 ) then   ! inner boundary
     ibegin2=Bc%ist; iend2=Bc%iend; jbegin2=Bc%jst; jend2=Bc%jend       ! write
     n1=iend2-ibegin2+jend2-jbegin2     ! +1 (number of cells)
  
  !                                     把连接点的信息读入临时数组Ua
     allocate(Ua(Nvar1,n1,LAP))             !存储交换信息的临时数组

       m_neighbour=Bc%neighb;  msub= Bc%subface   
       B1 =>  Mesh(nMesh)%Block(m_neighbour)
       Bc1 => B1%bc_msg(msub)

       ibegin1=Bc1%ist; iend1=Bc1%iend; jbegin1=Bc1%jst; jend1=Bc1%jend   ! read
 
       if(Bc1%face .eq. 1 ) then                     !  boundary  i-
         do j=jbegin1,jend1-1
		 do i=1,LAP   
!		   Ua(:,j,i)=B1%U(:,ibegin1+i-1,j) 
		   Ua(:,j-jbegin1+1,i)=B1%U(:,ibegin1+i-1,j) 
         enddo
	     enddo

	   else if ( Bc1%face .eq. 2) then  !  boundary  j-
         do  j=1,LAP
		 do  i=ibegin1,iend1-1  
!		   Ua(:,i,j)=B1%U(:,i,jbegin1+j-1)
		   Ua(:,i-ibegin1+1,j)=B1%U(:,i,jbegin1+j-1)
         enddo
		 enddo

	   else if( Bc1%face .eq. 3) then    !  boundary  i+
         do j=jbegin1, jend1-1
		 do i=1,LAP  
!		   Ua(:,j,i)=B1%U(:,iend1-LAP-1+i,j)  
		   Ua(:,j-jbegin1+1,i)=B1%U(:,iend1-LAP-1+i,j)  
         enddo
		 enddo

	   else                                !  boundary  j+
         do j=1,LAP
	     do i=ibegin1,iend1-1
!	 	   Ua(:,i,j)=B1%U(:,i,jend1-LAP-1+j)
	 	   Ua(:,i-ibegin1+1,j)=B1%U(:,i,jend1-LAP-1+j)
         enddo
		 enddo

	   endif

!                                    把临时数组Ua中的信息写入缓冲区
       orient=Bc%orient              ! orient==2 顺序;   其他值为 逆序
       do k=1,n1
        if(orient .eq. 2) then
          k1=k
        else
          k1=n1+1-k
        endif
       if(Bc%face .eq. 1 ) then             !  boundary  i-
         do i=1,LAP  
		   B%U(:,ibegin2-LAP-1+i,jbegin2+k-1)=Ua(:,k1,i)
         enddo

	   else if(Bc%face .eq. 2 ) then        !  boundary  j-
         do j=1,LAP   
		  B%U(:,ibegin2+k-1,jbegin2-LAP-1+j)=Ua(:,k1,j) 
         enddo
	   else if (Bc%face .eq. 3 ) then       !  boundary  i+
         do i=1,LAP
	       B%U(:,iend2+i-1,jbegin2+k-1)=Ua(:,k1,i)
         enddo
	   else                                          !  boundary  j+
         do j=1,LAP
	       B%U(:,ibegin2+k-1,jend2+j-1)=Ua(:,k1,j)
         enddo
      endif

      enddo
      deallocate(Ua)
     endif
    enddo
 
  !  处理缓冲区的角点，用插值的方法赋值     
        nx1=B%nx; ny1=B%ny
        call U_average_conner(Nvar1,B%U(:,1,0),B%U(:,1,1),B%U(:,0,1),B%U(:,0,0),Cv)
        call U_average_conner(Nvar1,B%U(:,1,ny1),B%U(:,1,ny1-1),B%U(:,0,ny1-1),B%U(:,0,ny1),Cv)
        call U_average_conner(Nvar1,B%U(:,nx1-1,0),B%U(:,nx1-1,1),B%U(:,nx1,1),B%U(:,nx1,0),Cv)
        call U_average_conner(Nvar1,B%U(:,nx1-1,ny1),B%U(:,nx1-1,ny1-1),B%U(:,nx1,ny1-1),B%U(:,nx1,ny1),Cv)
  
 
    enddo


   end  subroutine update_buffer_onemesh



!-----------------------------------------------------------------------------   
!    根据连接关系，将对应点的坐标写入缓冲区，并计算中心点坐标      
   subroutine Update_coordinate_buffer
    use Global_Var
    implicit none
    integer nMesh
	do nMesh=1,Num_Mesh
  	  call Update_coordinate_buffer_onemesh(nMesh)
	enddo
   end subroutine Update_coordinate_buffer

!---------Continue boundary (inner boundary) -------------------------------
!    根据连接关系，将对应点的坐标写入缓冲区，并计算中心点坐标 (处理1套网格) 
!    两块网格交界区布置 LAP层虚网格
!    交换虚网格的几何信息
!    对于物理边界，采用外插的方法构造虚网格的几何信息 （利用高阶外插）
     
	 subroutine Update_coordinate_buffer_onemesh(nMesh)
     use Global_Var
     implicit none
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B,B1
     Type (BC_MSG_TYPE),pointer:: Bc,Bc1
     real*8,allocatable:: Ux(:,:,:)
     integer:: nMesh,i,j,k,k1,m,mBlock,ksub,n1,m_neighbour,msub,orient,Kflag_initial
     integer:: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2
	 integer:: nx,ny,i1,i2,j1,j2
     character(len=30)::filename
 ! ------------------------------------------------------ 
  MP=>Mesh(nMesh)   ! 网格 （nMesh=3,2,1代表 粗、中、密网格） 

 do mBlock=1,MP%Num_Block
  B => MP%Block(mBlock)
  do  ksub=1,B%subface
   Bc=> B%bc_msg(ksub)
    ibegin2=Bc%ist; iend2=Bc%iend; jbegin2=Bc%jst; jend2=Bc%jend     ! write
    if(Bc%neighb .gt. 0 ) then   ! inner boundary
       n1=Bc%iend-Bc%ist+Bc%jend-Bc%jst+1  ! i_end-i_begin+j_end-j_begin+1
       allocate(Ux(2,n1,LAP))
       m_neighbour=Bc%neighb;  msub= Bc%subface   
       B1 => Mesh(nMesh)%Block(m_neighbour)
       Bc1=> B1%bc_msg(msub)
       ibegin1=Bc1%ist; iend1=Bc1%iend; jbegin1=Bc1%jst; jend1=Bc1%jend   ! read
 
       do i=1,LAP
         if(Bc1%face .eq. 1 ) then                     !  boundary  i-
  	      Ux(1,:,i)=B1%x(ibegin1+i,jbegin1:jend1) ;  Ux(2,:,i)=B1%y(ibegin1+i,jbegin1:jend1); 
  	     else if ( Bc1%face .eq. 2) then               !  boundary  j-
          Ux(1,:,i)=B1%x(ibegin1:iend1,jbegin1+i);   Ux(2,:,i)=B1%y(ibegin1:iend1,jbegin1+i)
         else if( Bc1%face .eq. 3) then                !  boundary  i+
          Ux(1,:,i)=B1%x(iend1-LAP-1+i,jbegin1:jend1) ;    Ux(2,:,i)=B1%y(iend1-LAP-1+i,jbegin1:jend1); 
         else                                                   !  boundary  j+
          Ux(1,:,i)=B1%x(ibegin1:iend1,jend1-LAP-1+i);     Ux(2,:,i)=B1%y(ibegin1:iend1,jend1-LAP-1+i)
         endif
       enddo

       orient=Bc%orient
       do k=1,n1
       if(orient .eq. 2) then
         k1=k
       else
         k1=n1+1-k
       endif
       do i=1,LAP
        if(Bc%face .eq. 1 ) then             !  boundary  i-
          B%x(ibegin2-LAP-1+i,jbegin2+k-1)=Ux(1,k1,i) ; B%y(ibegin2-LAP-1+i,jbegin2+k-1)=Ux(2,k1,i)
        else if(Bc%face .eq. 2 ) then        !  boundary  j-
          B%x(ibegin2+k-1,jbegin2-LAP-1+i)=Ux(1,k1,i) ; B%y(ibegin2+k-1,jbegin2-LAP-1+i)=Ux(2,k1,i)
        else if (Bc%face .eq. 3 ) then       !  boundary  i+
          B%x(iend2+i,jbegin2+k-1)=Ux(1,k1,i) ; B%y(iend2+i,jbegin2+k-1)=Ux(2,k1,i)
        else                                          !  boundary  j+
          B%x(ibegin2+k-1,jend2+i)=Ux(1,k1,i) ; B%y(ibegin2+k-1,jend2+i)=Ux(2,k1,i)
        endif
       enddo

	  enddo
      
	  deallocate(Ux)

    else  ! not inner boundary
	      ! 非内边界点，Ghost Cell的坐标采用外推方法获得 （有待改进）
	  do i=1,LAP
       if(Bc%face .eq. 1 ) then             !  boundary  i-
          B%x(ibegin2-i,jbegin2:jend2)=2.d0*B%x(ibegin2,jbegin2:jend2)-B%x(ibegin2+i,jbegin2:jend2) 
          B%y(ibegin2-i,jbegin2:jend2)=2.d0*B%y(ibegin2,jbegin2:jend2)-B%y(ibegin2+i,jbegin2:jend2)   
       else if(Bc%face .eq. 2 ) then        !  boundary  j-
          B%x(ibegin2:iend2,jbegin2-i)=2.d0*B%x(ibegin2:iend2,jbegin2)-B%x(ibegin2:iend2,jbegin2+i) 
          B%y(ibegin2:iend2,jbegin2-i)=2.d0*B%y(ibegin2:iend2,jbegin2)-B%y(ibegin2:iend2,jbegin2+i)  
       else if (Bc%face .eq. 3 ) then       !  boundary  i+
          B%x(iend2+i,jbegin2:jend2)=2.d0*B%x(iend2,jbegin2:jend2)-B%x(iend2-i,jbegin2:jend2) 
          B%y(iend2+i,jbegin2:jend2)=2.d0*B%y(iend2,jbegin2:jend2)-B%y(iend2-i,jbegin2:jend2) 
       else                                          !  boundary  j+
          B%x(ibegin2:iend2,jend2+i)=2.d0*B%x(ibegin2:iend2,jend2)-B%x(ibegin2:iend2,jend2-i)
          B%y(ibegin2:iend2,jend2+i)=2.d0*B%y(ibegin2:iend2,jend2)-B%y(ibegin2:iend2,jend2-i)
       endif
	   enddo
     endif
    enddo
  !  Conner point     
       nx=B%nx; ny=B%ny
!       B%x(0,0)=B%x(1,0)+B%x(0,1)-B%x(1,1);   B%y(0,0)=B%y(1,0)+B%y(0,1)-B%y(1,1) 
!       B%x(0,ny1+1)=B%x(1,ny1+1)+B%x(0,ny1)-B%x(1,ny1);  B%y(0,ny1+1)=B%y(1,ny1+1)+B%y(0,ny1)-B%y(1,ny1) 
!       B%x(nx1+1,0)=B%x(nx1,0)+B%x(nx1+1,1)-B%x(nx1,1);  B%y(nx1+1,0)=B%y(nx1,0)+B%y(nx1+1,1)-B%y(nx1,1) 
!       B%x(nx1+1,ny1+1)=B%x(nx1,ny1+1)+B%x(nx1+1,ny1)-B%x(nx1,ny1); B%y(nx1+1,ny1+1)=B%y(nx1,ny1+1)+B%y(nx1+1,ny1)-B%y(nx1,ny1) 
  ! 区域的角点（共LAP*LAP个）
      do j=1,LAP
	  do i=1,LAP
	   i1=1-i; j1=1-j
	   i2=nx+i;j2=ny+j 
	   B%x(i1,j1)=B%x(i1,1)+B%x(1,j1)-B%x(1,1); B%y(i1,j1)=B%y(i1,1)+B%y(1,j1)-B%y(1,1) 
  	   B%x(i2,j2)=B%x(i2,ny)+B%x(nx,j2)-B%x(nx,ny); B%y(i2,j2)=B%y(i2,ny)+B%y(nx,j2)-B%y(nx,ny) 
 	   B%x(i1,j2)=B%x(i1,ny)+B%x(1,j2)-B%x(1,ny); B%y(i1,j2)=B%y(i1,ny)+B%y(1,j2)-B%y(1,ny)
 	   B%x(i2,j1)=B%x(i2,1)+B%x(nx,j1)-B%x(nx,1); B%y(i2,j1)=B%y(i2,1)+B%y(nx,j1)-B%y(nx,1)
     enddo
	 enddo
	    

  ! Coordinates of the Cell center 
     do j=1-LAP,ny+LAP-1
     do i=1-LAP,nx+LAP-1
	  B%x1(i,j)=(B%x(i,j)+B%x(i+1,j)+B%x(i,j+1)+B%x(i+1,j+1))*0.25
	  B%y1(i,j)=(B%y(i,j)+B%y(i+1,j)+B%y(i,j+1)+B%y(i+1,j+1))*0.25
     enddo
	 enddo

   enddo
 
! 输出含缓冲区的网格信息（tecplot格式），便于调试
    write(filename,"('mesh-test'I1.1'.dat')") nMesh
	 open(99,file=filename)           ! mesh-test1.dat 细网格， mesh-test2.dat 粗网格， mesh-test3.dat 更粗网格
        write(99,*) "variables=x,y"
        do m=1,Mesh(nMesh)%Num_Block
          B => Mesh(nMesh)%Block(m)
!          write(99,*)  " zone ", "i= ", B%nx+2, " j= ", B%ny+2
!  	      do j=0,B%ny+1
!         do i=0,B%nx+1
           write(99,*)  " zone ", "i= ", B%nx+4, " j= ", B%ny+4
  	      do j=-1,B%ny+2
          do i=-1,B%nx+2
          write(99,*) B%x(i,j),B%y(i,j)
          enddo
          enddo
         enddo
      close(99)
    end  subroutine Update_coordinate_buffer_onemesh

!-------------------------------------------------------------------
!  处理壁面边界条件
!  使用LAP层虚网格
    subroutine boundary_wall(nMesh,mBlock,ksub)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP

     integer:: nMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i2,j2,i1,j1,k
     real*8:: d1,p1,T1,d2,p2,T2,u1,v1,u2,v2,Tsb,Amu2,wt
     real*8,parameter:: beta1_SST=0.075d0
     
     MP => Mesh(nMesh)
     B  => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)

     Tsb=110.4d0/T_inf       ! 参考温度 （Surthland公式中使用)
     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      

       do j=jbegin,jend
       do i=ibegin,iend
       do k=1,LAP
    !      (i1,j1)  内点；  (i2,j2) 为对应的buffer区的点   
		 if(Bc%face .eq. 1) then
           i1=i+k-1; j1=j; i2=i-k ; j2=j  
         else if(Bc%face .eq. 2) then
           i1=i; j1=j+k-1; i2=i; j2=j-k
         else if(Bc%face .eq. 3) then
           i1=i-k; j1=j;  i2=i+k-1; j2=j
         else
           i1=i; j1=j-k; i2=i; j2=j+k-1
         endif
        
		if(Twall .le. 0) then   ! 绝热壁

 		   B%U(1,i2,j2)= B%U(1,i1,j1)       ! d(0)=d(1)   对称
           B%U(2,i2,j2)=-B%U(2,i1,j1)       ! u(0)=-u(1)  -> d*u 反对称
           B%U(3,i2,j2)=-B%U(3,i1,j1)       ! v(0)=-v(1)  -> d*v 反对称
           B%U(4,i2,j2)= B%U(4,i1,j1)       ! E(0)=E(1)   对称
        
		else   ! 等温壁
            
			d1=B%U(1,i1,j1)    ! 内点处的密度、压力、温度、速度
			u1=B%U(2,i1,j1)/d1 
			v1=B%U(3,i1,j1)/d1
		    p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)             
            T1=gamma*Ma*Ma*p1/d1     
			 
            p2=p1               ! 边界层假设，壁面处法向压力梯度为0
            T2=2.d0*Twall-T1    ! 等温壁  0.5*(T1+T2)=Twall
			u2=-u1              ! 无滑移壁
			v2=-v1
			d2=gamma*Ma*Ma*p2/T2 

             B%U(1,i2,j2)=d2
             B%U(2,i2,j2)=d2*u2
             B%U(3,i2,j2)=d2*v2
  			 B%U(4,i2,j2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2)
       endif
          if(MP%Nvar .eq. 5) then
            B%U(5,i2,j2)=0.d0
		  endif

	      if(MP%Nvar .eq. 6) then
             d2=B%U(1,i2,j2)
             T2=(B%U(4,i2,j2)-0.5d0*(B%U(2,i2,j2)**2+B%U(3,i2,j2)**2)/d2 )/(Cv*d2)
             Amu2=(1.d0+Tsb)*sqrt(T2**3)/(Tsb+T2)   !层流粘性系数，sutherland equation

!              B%U(5,i2,j2)=0.d0                               ! 湍动能 （固壁上为0）
!              B%U(6,i2,j2)=10.d0*6.d0*Amu2/(d2*beta1_SST*B%dw(i1,j1)**2*Re*Re)   ! 湍能比耗散率 , Bug removed 2012-5-10                  
         
		 	   wt=60.d0*Amu2/(d2*beta1_SST*B%dw(i1,j1)**2*Re*Re)
 			   B%U(5,i2,j2)=- B%U(5,i1,j1)           ! 壁面镜像点， 使得壁面上k=0
               B%U(6,i2,j2)=(d1+d2)*wt-B%U(6,i1,j1)
 
           endif
       enddo
	   enddo
       enddo 
      
      end

!-------------------------------------------------------------------
!  对称边界条件(或滑移固壁)
!  使用LAP层虚网格
    subroutine Symmetry_or_slidewall(nMesh,mBlock,ksub)
     use Global_Var
     implicit none
     
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP

     integer:: nMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i2,j2,i1,j1,k
     real*8:: p2,dx,dy,si,n1,n2,Vn
     
     MP=>Mesh(nMesh)
     B => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)
     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      

       do j=jbegin,jend
       do i=ibegin,iend
!      计算边界法方向    
	     if(Bc%face .eq. 1 .or. Bc%face .eq. 3) then    ! i+ or i- boundary
	      dx=B%x(i,j+1)-B%x(i,j) ;      dy=B%y(i,j+1)-B%y(i,j)
	      si=sqrt(dx*dx+dy*dy)
          n1=dy/si; n2=-dx/si   ! normal vector at (i,j) or (I-1/2,J) 
      
		 else   ! j+ or j-
	      dx=B%x(i+1,j)-B%x(i,j) ;      dy=B%y(i+1,j)-B%y(i,j)
  	      si=sqrt(dx*dx+dy*dy)
          n1=-dy/si; n2=dx/si   ! normal vector at i, j+1/2
 		 endif


       do k=1,LAP
         if(Bc%face .eq. 1) then
           i1=i+k-1; j1=j; i2=i-k ; j2=j  
         else if(Bc%face .eq. 2) then
           i1=i; j1=j+k-1; i2=i; j2=j-k
         else if(Bc%face .eq. 3) then
           i1=i-k; j1=j;  i2=i+k-1; j2=j
         else
           i1=i; j1=j-k; i2=i; j2=j+k-1
         endif
           Vn=B%U(2,i1,j1)*n1+B%U(3,i1,j1)*n2   ! 法向动量

         B%U(1,i2,j2)= B%U(1,i1,j1)       ! d(0)=d(1)   对称
         B%U(2,i2,j2)= B%U(2,i1,j1)-2.d0*Vn*n1       ! 法向动量相反，切向动量不变
         B%U(3,i2,j2)= B%U(3,i1,j1)-2.d0*Vn*n2       ! 法向动量相反，切向动量不变
         B%U(4,i2,j2)= B%U(4,i1,j1)       ! E(0)=E(1)   对称
	   
	    if(MP%Nvar .eq. 5) then
	      B%U(5,i2,j2)=B%U(5,i1,j1)     ! 标量，对称
	    elseif(MP%Nvar .eq. 6) then         ! k,w
	      B%U(5,i2,j2)=B%U(5,i1,j1)     ! 标量，对称
	      B%U(6,i2,j2)=B%U(6,i1,j1)
	    endif
	    
	   enddo
	   enddo
       enddo 
      
      end


!-------------------------------------------------------------------------

! 可处理远场边界条件 以及入口/出口边界条件
   subroutine boundary_Farfield_inlet_outlet(nMesh,mBlock,ksub,Bctype)
     use Global_Var
     implicit none
     Type (Block_TYPE),pointer:: B
     Type (BC_MSG_TYPE),pointer:: Bc
     Type (Mesh_TYPE),pointer:: MP

     integer:: Bctype,NMesh,mBlock,ksub,ibegin,iend,jbegin,jend,i,j,i1,j1,i2,j2,k
     real*8:: d_inf,u_inf,v_inf,p_inf,p_out,n1,n2 
     real*8:: d1,u1,v1,p1,c1,d2,u2,v2,p2,pb,db,ub,vb
     real*8:: dx,dy,s,Ma_n

!  本软件目前用来计算内流，给定无穷远条件
     d_inf=1.d0; u_inf=1.d0*cos(AoA); v_inf=1.d0*sin(AoA) ; p_inf=1.d0/(gamma*Ma*Ma)

     MP=>Mesh(nMesh)
     B => MP%Block(mBlock)
     Bc => B%bc_msg(ksub)

     ibegin=Bc%ist; iend=Bc%iend; jbegin=Bc%jst; jend=Bc%jend      
      
       do j=jbegin,jend
       do i=ibegin,iend
!-----------------------------------------------------------------
 !  (i1,j1) 是靠近边界的内点， (i2,j2) 是边界外的Ghost Cell点    
         if(Bc%face .eq. 1) then       ! i- 面， 
           i1=i; j1=j  
         else if(Bc%face .eq. 2) then  ! j- 面
           i1=i; j1=j
         else if(Bc%face .eq. 3) then  ! i+ 面 (i=ibegin=iend=nx), i1=i-1 是内点, i2=i=nx是Ghost Cell
           i1=i-1; j1=j
         else
           i1=i; j1=j-1
         endif
            d1=B%U(1,i1,j1) ;   u1=B%U(2,i1,j1)/d1 ;     v1=B%U(3,i1,j1)/d1
            p1=(B%U(4,i1,j1)-0.5d0*d1*(u1*u1+v1*v1))*(gamma-1.d0)              ! 内点处的值
            c1=sqrt(gamma*p1/d1) 

!   计算边界 外 法方向       
          if( Bc%face .eq. 1 ) then 
             dx=B%x(i,j+1)-B%x(i,j) ;    dy=B%y(i,j+1)-B%y(i,j);  s=sqrt(dx*dx+dy*dy)
             n1=-dy/s; n2= dx/s    ! 外法线  
          else if(Bc%face .eq. 3) then 
             dx=B%x(i,j+1)-B%x(i,j) ;    dy=B%y(i,j+1)-B%y(i,j);  s=sqrt(dx*dx+dy*dy)
             n1=dy/s; n2= -dx/s    ! 外法线  
          else if(Bc%face .eq. 2) then
             dx=B%x(i+1,j)-B%x(i,j) ;   dy=B%y(i+1,j)-B%y(i,j) ;  s=sqrt(dx*dx+dy*dy)   
             n1=dy/s; n2=-dx/s      ! 外法线 
          else
             dx=B%x(i+1,j)-B%x(i,j) ;   dy=B%y(i+1,j)-B%y(i,j) ;  s=sqrt(dx*dx+dy*dy)   
             n1=-dy/s; n2=dx/s        ! 外法线 
          endif
            
!			if(KS_Farfield == 0) then       ! 远场处理方式0 （均按照超声速情况处理）
!			 Ma_n=(u_inf*n1+v_inf*n2)*Ma    ! 以无穷远来流量计算 （避免了出口回流区的误算，2012-5-23）
!            else
			 Ma_n=(u1*n1+v1*n2)/c1    ! 法向Mach数
!            endif

!-------------------------------------------------------------------------------------------
          do k=1,LAP
	       if(Bc%face .eq. 1) then       ! i- 面， 
             i2=i-k ; j2=j  
           else if(Bc%face .eq. 2) then  ! j- 面
             i2=i; j2=j-k
           else if(Bc%face .eq. 3) then  ! i+ 面 (i=ibegin=iend=nx), i1=i-1 是内点, i2=i=nx是Ghost Cell
             i2=i+k-1; j2=j
           else
            i2=i; j2=j+k-1
          endif
	   
	      !
		   if(Bctype .eq. BC_outlet .or. (Bctype .eq. BC_Farfield .and.  Ma_n >= 0.d0 ) ) then   ! 出口 
              if(Ma_n > 1.d0 .or. P_outlet < 0 .or. Bctype== BC_Farfield) then 
			    B%U(:,i2,j2)=B%U(:,i1,j1)          ! 超声速出口，外推 (BC_Farfield 均按照外推处理)
		      else   ! 亚声速出口，给定背压
  			   pb=p_outlet 
               db=d1+(p_outlet-p1)/(c1*c1)
               ub=u1+(p1-p_outlet)/(d1*c1)*n1
               vb=v1+(p1-p_outlet)/(d1*c1)*n2
               p2=2.d0*pb-p1 ; d2=2.d0*db-d1 ; u2=2.d0*ub-u1 ; v2=2.d0*vb-v1
               B%U(1,i2,j2)=d2
               B%U(2,i2,j2)=d2*u2
               B%U(3,i2,j2)=d2*v2
               B%U(4,i2,j2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2)
			    if(Mp%Nvar >= 5) then             ! 湍流模型有关量
                   B%U(5:MP%Nvar,i2,j2)=B%U(5:MP%Nvar,i1,j1)
		        endif
			  endif

		   
		   else        ! 入口

             if(P_outlet<0 .or. Ma_n <= -1.d0 ) then   ! 超声速情况 (设定P_outlet<0 默认按超声速情况处理入口和出口边界)
		      B%U(1,i2,j2)=d_inf
              B%U(2,i2,j2)=d_inf*u_inf
              B%U(3,i2,j2)=d_inf*v_inf
              B%U(4,i2,j2)=p_inf/(gamma-1.d0)+0.5d0*d_inf*(u_inf*u_inf+v_inf*v_inf)
             else
		      pb=0.5d0*(p1+p_inf-d1*c1*((u_inf-u1)*n1+(v_inf-v1)*n2) )
              db=d_inf+(pb-p_inf)/(c1*c1)
              ub=u_inf-(p_inf-pb)/(d1*c1)*n1
              vb=v_inf-(p_inf-pb)/(d1*c1)*n2
              p2=2.d0*pb-p1 ; d2=2.d0*db-d1 ; u2=2.d0*ub-u1 ; v2=2.d0*vb-v1
              B%U(1,i2,j2)=d2
              B%U(2,i2,j2)=d2*u2
              B%U(3,i2,j2)=d2*v2
              B%U(4,i2,j2)=p2/(gamma-1.d0)+0.5d0*d2*(u2*u2+v2*v2)
             endif
			   
			  
			  if(Mp%Nvar .eq. 5) then
               B%U(5,i2,j2)=vt_inf/Re          ! 来流vt值（通常设定层流粘性系数的3-5倍）
			  else if(MP%Nvar .eq. 6) then
               B%U(5,i2,j2)=d_inf*Kt_inf      ! 来流值
               B%U(6,i2,j2)=d_inf*Wt_inf
             endif
          endif  



        enddo
	   enddo
      enddo 
    end subroutine boundary_Farfield_inlet_outlet

!-------------------------------------------------------------------------









!-------------------------------
! 计算缓冲区角点处的值 （例如U1(0,0)点的值）， 用外插方法计算
  subroutine U_average_conner(Nvar1,U1,U2,U3,U4,Cv)
  implicit none
  integer::Nvar1
  real*8,dimension(Nvar1):: U1,U2,U3,U4
  real*8:: Cv,d1,uu1,v1,T1,d2,uu2,v2,T2,d3,uu3,v3,T3,d4,uu4,v4,T4
    d1=U1(1); uu1=U1(2)/d1; v1=U1(3)/d1; T1=(U1(4)-(uu1*U1(2)+v1*U1(3))*0.5d0)/(d1*Cv)  ! density, velocity, Temperature 
    d2=U2(1); uu2=U2(2)/d2; v2=U2(3)/d2; T2=(U2(4)-(uu2*U2(2)+v2*U2(3))*0.5d0)/(d2*Cv)   
    d3=U3(1); uu3=U3(2)/d3; v3=U3(3)/d3; T3=(U3(4)-(uu3*U3(2)+v3*U3(3))*0.5d0)/(d3*Cv) 
    d4=d1+d3-d2; uu4=uu1+uu3-uu2; v4=v1+v3-v2; T4=T1+T3-T2
    U4(1)=d4; U4(2)=d4*uu4; U4(3)=d4*v4; U4(4)=d4*(Cv*T4+(uu4*uu4+v4*v4)*0.5d0)
    if(Nvar1 .eq. 5) then
      U4(5)=U1(5)+U3(5)-U2(5)
	endif
    if(Nvar1 .eq. 6) then
      U4(5)=U1(5)+U3(5)-U2(5)
      U4(6)=U1(6)+U3(6)-U2(6)
    endif
  end subroutine U_average_conner
