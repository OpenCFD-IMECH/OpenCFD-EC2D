!------------------------------------------------------------------------------
! Source term of SA model  (See: J. Blazek's book, P240-243) 
! Do not consider transition (full turbulence)
  subroutine  Turbulence_model_SA(nMesh,mBlock)
   use Global_Var
   Use Flow_Var
   implicit none
   integer:: nMesh,mBlock,nx,ny,i,j

   real*8:: S, S1,X,fv1,fv2,ft2,r,g,fw,Q_SA,vn1,vn2,vfi,s0,n1,n2,v0
   real*8:: Dix,Diy,Djx,Djy,Ds,Div,Djv,vtx,vty,Diu,Dju,ux,vx,uy,vy
   real*8,parameter:: SA_sigma=2.d0/3.d0,Cv1=7.1d0,Cv2=5.d0,Cb1=0.1355d0,Cb2=0.622d0,SA_k=0.41d0
   real*8,parameter:: Cw1=Cb1/(SA_k*SA_k)+(1.d0+Cb2)/SA_sigma, Cw2=0.3d0, Cw3=2.d0, Ct1=1.d0, Ct2=2.d0, Ct3=1.3d0, Ct4=0.5d0
   real*8,Pointer,dimension(:,:):: vt,Fluxv,fluxv2
   Type (Block_TYPE),pointer:: B
 
  
   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny
 ! 计算湍流粘性系数
   allocate(vt(0:nx,0:ny),Fluxv(nx,ny),fluxv2(nx,ny))
	
! OpenMP的编译指示符（不是注释）， 指定Do 循环并行执行； 指定一些各进程私有的变量
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx,ny,B,vt,d)
    do j=0,ny
    do i=0,nx
     vt(i,j)=B%U(5,i,j)
	 X=d(i,j)*vt(i,j)/B%Amu(i,j)   ! 湍流粘性系数与层流粘性系数之比
     fv1=X**3/(X**3+Cv1**3)
     B%Amu_t(i,j)=fv1*d(i,j)*vt(i,j)
	enddo
    enddo
!$OMP END PARALLEL DO

! 设定湍流粘性系数虚网格的值
 call  Amut_boundary(nMesh,mBlock)

! 计算vt方程的残差 B%Res(6,:,:,:)

!$OMP PARALLEL  DEFAULT(PRIVATE) SHARED(nx,ny,B,vt,d,uu,v,Fluxv,fluxv2)
!------i- direcion ---------------------------
!$OMP DO
     do j=1,ny-1
     do i=1,nx
     
	 s0=B%si(i,j)   ! 控制体边界长度
	 n1=B%ni1(i,j); n2=B%ni2(i,j)   ! 法方向

 !   对流项，采用1阶迎风格式      
     vn1=uu(i-1,j)*n1+v(i-1,j)*n2         ! 法向速度分量
	 vn2=uu(i,j)*n1+v(i,j)*n2
!   1阶L-F 格式
	 vfi=-0.5d0*( (vn1+abs(vn1))*vt(i-1,j) + (vn2-abs(vn2))*vt(i,j) )*s0   ! 无粘通量
     v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j)+B%Amu(i,j)/d(i,j) + vt(i-1,j)+B%Amu(i-1,j)/d(i-1,j))  ! (I-1/2,J) 点的动力学粘性系数
     
! 计算物理量（k,w）对坐标x,y的导数 （采用Jocabian变换）
!----Jocabian系数 （物理坐标对计算坐标的导数, 用于计算物理量的导数）
      Dix=B%x1(i,j)-B%x1(i-1,j)
      Diy=B%y1(i,j)-B%y1(i-1,j)
      Djx=(B%x1(i-1,j+1)+B%x1(i,j+1)-B%x1(i-1,j-1)-B%x1(i,j-1))*0.25d0
      Djy=(B%y1(i-1,j+1)+B%y1(i,j+1)-B%y1(i-1,j-1)-B%y1(i,j-1))*0.25d0
      Ds=1.d0/(Dix*Djy-Djx*Diy)
! 物理量对计算坐标的导数    
      Div=vt(i,j)-vt(i-1,j)
      Djv=(vt(i-1,j+1)+vt(i,j+1)-vt(i-1,j-1)-vt(i,j-1))*0.25d0
! 物理量对x,y坐标的导数
      vtx=(Div*Djy-Djv*Diy)*Ds
      vty=(-Div*Djx+Djv*Dix)*Ds
!  粘性应力及能量通量
      fluxv(i,j)=vfi+v0*(vtx*n1+vty*n2)*s0
      fluxv2(i,j)=(vtx*n1+vty*n2)*s0
	  enddo
     enddo
!$OMP END DO

!$OMP DO
      do j=1,ny-1
      do i=1,nx-1
	     v0=(vt(i,j)+B%Amu(i,j)/d(i,j))*Cb2/SA_sigma
         B%Res(5,i,j)=Fluxv(i+1,j)-Fluxv(i,j) -v0*(Fluxv2(i+1,j)-Fluxv2(i,j))       
      enddo
      enddo
!$OMP END DO
!----- j-direction-----------------------

!$OMP DO
    do j=1,ny
    do i=1,nx-1
 	 s0=B%sj(i,j)
 	 n1=B%nj1(i,j); n2=B%nj2(i,j)   ! 法方向

!   对流项，采用1阶迎风格式 （L-F分裂）
     
	 vn1=uu(i,j-1)*n1+v(i,j-1)*n2   ! 法向速度
	 vn2=uu(i,j)*n1+v(i,j)*n2
	 vfi=-0.5d0*( (vn1+abs(vn1))*vt(i,j-1) + (vn2-abs(vn2))*vt(i,j) )*s0   ! 无粘通量

!  粘性项
!---------Viscous term -----------------------------------------------------------------------------
    v0=0.5d0*(1.d0+Cb2)/SA_sigma*(vt(i,j)+B%Amu(i,j)/d(i,j) + vt(i,j-1)+B%Amu(i,j-1)/d(i,j-1))  ! (I,J-1/2) 点的动力学粘性系数  
! 计算物理量（k,w）对坐标x,y的导数 （采用Jocabian变换）
    Dix=(B%x1(i+1,j-1)+B%x1(i+1,j)-B%x1(i-1,j-1)-B%x1(i-1,j))*0.25d0
    Diy=(B%y1(i+1,j-1)+B%y1(i+1,j)-B%y1(i-1,j-1)-B%y1(i-1,j))*0.25d0
    Djx=B%x1(i,j)-B%x1(i,j-1)
    Djy=B%y1(i,j)-B%y1(i,j-1)
    Ds=1.d0/(Dix*Djy-Djx*Diy)
    Div=(vt(i+1,j-1)+vt(i+1,j)-vt(i-1,j-1)-vt(i-1,j))*0.25d0
    Djv=vt(i,j)-vt(i,j-1)
    vtx=(Div*Djy-Djv*Diy)*Ds
    vty=(-Div*Djx+Djv*Dix)*Ds
    fluxv(i,j)=vfi+v0*(vtx*n1+vty*n2)*s0
    fluxv2(i,j)=(vtx*n1+vty*n2)*s0
   enddo
   enddo
!$OMP END DO

!$OMP DO
      do j=1,ny-1
      do i=1,nx-1
	     v0=(vt(i,j)+B%Amu(i,j)/d(i,j))*Cb2/SA_sigma
         B%Res(5,i,j)=B%Res(5,i,j)+Fluxv(i,j+1)-Fluxv(i,j) -v0*(Fluxv2(i,j+1)-Fluxv2(i,j))       
      enddo
      enddo
!$OMP END DO


!--------------源项---------------------------------------------------------
!$OMP DO
   do j=1,ny-1
   do i=1,nx-1
!----- get S (normal of vorticity)  S=sqrt(0.5*Omiga_ij*Omiga_ij) at the cell's center ------
!  计算涡量
! 计算物理量对坐标x,y的导数，使用Jocabian变换
     Dix=(B%x1(i+1,j)-B%x1(i-1,j))*0.5d0
     Diy=(B%y1(i+1,j)-B%y1(i-1,j))*0.5d0
     Djx=(B%x1(i,j+1)-B%x1(i,j-1))*0.5d0
     Djy=(B%y1(i,j+1)-B%y1(i,j-1))*0.5d0
     Ds=1.d0/(Dix*Djy-Djx*Diy)

     Diu=(uu(i+1,j)-uu(i-1,j))*0.5d0
     Div=(v(i+1,j)-v(i-1,j))*0.5d0
     Dju=(uu(i,j+1)-uu(i,j-1))*0.5d0
     Djv=(v(i,j+1)-v(i,j-1))*0.5d0

! 导数值    
     ux=(Diu*Djy-Dju*Diy)*Ds
     vx=(Div*Djy-Djv*Diy)*Ds
     uy=(-Diu*Djx+Dju*Dix)*Ds
     vy=(-Div*Djx+Djv*Dix)*Ds
     S=abs(vx-uy)      ! 涡量的模

!--------------------------------------------------------------------------   
! Source term, original form; See: http://turbmodels.larc.nasa.gov/spalart.html

   X=d(i,j)*vt(i,j)/B%Amu(i,j)   ! 湍流粘性系数与层流粘性系数之比
   fv1=X**3/(X**3+Cv1**3)

   fv2=1.d0-X/(1.d0+X*fv1)
   S1=max(S+fv2*vt(i,j)/(SA_k*B%dw(i,j))**2,0.3*S)
   r=min(vt(i,j)/(S1*SA_K*SA_K*B%dw(i,j)*B%dw(i,j)),10.0)

   g=r+Cw2*(r**6-r)
   fw=g*((1.d0+Cw3**6)/(g**6+Cw3**6))**(1.d0/6.d0)
   ft2=Ct3*exp(-Ct4*X*X)            
   Q_SA=Cb1*(1.d0-ft2)*S1*vt(i,j) -(Cw1*fw-Cb1*ft2/(SA_k*SA_k))*(vt(i,j)/B%dw(i,j))**2
!---------------------------------------------------------------------------
   B%Res(5,i,j)=B%Res(5,i,j)+Q_SA*B%vol(i,j)     
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

  deallocate(vt,Fluxv,fluxv2)

end  subroutine  Turbulence_model_SA

!-----------------------------------------------------------------------------------

! 粘性系数虚网格上的值 （固壁边界采用反值，以保证固壁上的平均湍流粘性系数为0）
subroutine Amut_boundary(nMesh,mBlock)
   use Global_Var
   use Flow_Var
   implicit none
   integer:: mBlock,nx,ny,i,j,m,ksub,nMesh
   integer:: ib,ie,jb,je
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc

   B => Mesh(nMesh)%Block(mBlock)
   nx=B%nx ; ny=B%ny

! Ghost Cell 点的 mut值为 内点mut值*（-1） (这样可以使壁面上mut=0)

  do  ksub=1,B%subface
    Bc=> B%bc_msg(ksub)
      if( Bc%neighb .eq. BC_Wall  ) then   ! (粘性) 壁面边界条件 
      Bc => B%bc_msg(ksub)
      ib=Bc%ist; ie=Bc%iend; jb=Bc%jst; je=Bc%jend       

     if(Bc%face .eq. 1 ) then   ! i- 
       do j=jb,je-1
         B%Amu_t(0,j)= -B%Amu_t(1,j)           
       enddo
     else if (Bc%face .eq. 3 ) then   ! i+
       do j=jb,je-1
         B%Amu_t(nx,j)= -B%Amu_t(nx-1,j)       ! mut
       enddo
     else if(Bc%face .eq. 2 ) then   !j-
       do i=ib,ie-1
         B%Amu_t(i,0)= -B%Amu_t(i,1)       ! mut
      enddo
     else if(Bc%face .eq. 4 ) then   !j+
       do i=ib,ie-1
         B%Amu_t(i,ny)= -B%Amu_t(i,ny-1)       ! mut
       enddo
     endif
    endif
   enddo

end  subroutine Amut_boundary
