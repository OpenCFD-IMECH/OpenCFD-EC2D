! 用差分法计算通量 （相当于补丁程序）
! 多重网格情况下，仅在最密的网格采用 （稀网格通常不采用高精度方法）

   subroutine Residual_FDM(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
!   integer,parameter:: LFDM=4    ! 差分法块的边界网格（块边缘的4层网格不使用差分法）
   integer,parameter:: NMax=10000    ! 一维数组的长度
   integer:: nMesh,mBlock
   integer:: nx,ny,i,j,m
   real*8:: A1,A2,ui,vi,Ti,uj,vj,Tj,ux,vx,Tx,uy,vy,Ty,t11,t12,t22,E1,E2,mu0,Ak0
   real*8,allocatable,dimension(:,:,:):: EV1,EV2
   real*8:: fp(4),fm(4),fx1(NMax,4),fx2(NMax,4),hi1(NMax),hi2(NMax)    ! 使用静态数组，便于利用OpenMP
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
 
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   
   nx=B%nx-1    !  差分法的网格（由原网格中心点构成的新网格）
   ny=B%ny-1
   
   if(max(nx,ny) .gt. NMax) then
     print*, "Error in Residual_FDM() !!! nx or ny > NMAX !"
	 stop
   endif

   allocate(EV1(nx,ny,3),Ev2(nx,ny,3))

! 无粘通量的计算
!--------i- -----------------------------------
! LFDM 差分块临近边界区的网格层数（通常为4，这些层网格不使用差分法计算）

!$OMP PARALLEL DEFAULT (PRIVATE) SHARED(nx,ny,MP,B,d,uu,v,p,cc,T,Fluxi,Fluxj,Ev1,Ev2,If_viscous,gamma,Cp,Pr)
!$OMP DO
        do j=LFDM,ny-LFDM
        do i=1,nx
	     A1=B%Akx(i,j)/B%Ajac(i,j)
		 A2=B%Aky(i,j)/B%Ajac(i,j)
         call split_Steger_Warming(d(i,j),uu(i,j),v(i,j),cc(i,j),fP,fm,A1,A2,gamma)   ! Steger-Warming分裂
        do m=1,4
	     fx1(i,m)=fp(m)
	     fx2(i,m)=fm(m)
	    enddo
        enddo
        do m=1,4
	    call fp_weno5(nx,fx1(1,m),hi1)                         ! 调用差分格式
	    call fm_weno5(nx,fx2(1,m),hi2)
	    do i=LFDM,nx-LFDM
	     Fluxi(m,i,j)=-(hi1(i)+hi2(i))
	    enddo
	    enddo
       enddo
!$OMP END DO

!-------j- -------------------------------------
!$OMP DO
        do i=LFDM,nx-LFDM
        do j=1,ny
		 A1=B%Aix(i,j)/B%Ajac(i,j)
         A2=B%Aiy(i,j)/B%Ajac(i,j)
      	 call split_Steger_Warming(d(i,j),uu(i,j),v(i,j),cc(i,j),fP,fm,A1,A2,gamma)
	     do m=1,4
	      fx1(j,m)=fp(m)
	      fx2(j,m)=fm(m)
         enddo
        enddo
		do m=1,4
		 call fp_weno5(ny,fx1(1,m),hi1)
	     call fm_weno5(ny,fx2(1,m),hi2)
	     do j=LFDM,ny-LFDM
	      Fluxj(m,i,j)=-(hi1(j)+hi2(j))
	     enddo
        enddo
       enddo
!$OMP END DO

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!   粘性通量的计算

!  计算应力张量与热流项  (2阶中心差分)
  if(If_viscous .eq. 1) then

!$OMP DO
   do j=2,ny-1
   do i=2,nx-1
     mu0=B%Amu(i,j)+B%Amu_t(i,j)
	 Ak0=Cp*(B%Amu(i,j)/Pr+B%Amu_t(i,j)/Prt)
     ui=0.5d0*(uu(i+1,j)-uu(i-1,j))
     vi=0.5d0*(v(i+1,j)-v(i-1,j))
     Ti=0.5d0*(T(i+1,j)-T(i-1,j))
     uj=0.5d0*(uu(i,j+1)-uu(i,j-1))
     vj=0.5d0*(v(i,j+1)-v(i,j-1))
     Tj=0.5d0*(T(i,j+1)-T(i,j-1))
     ux=ui*B%Akx(i,j)+uj*B%Aix(i,j)
     vx=vi*B%Akx(i,j)+vj*B%Aix(i,j)
     Tx=Ti*B%Akx(i,j)+Tj*B%Aix(i,j)
     uy=ui*B%Aky(i,j)+uj*B%Aiy(i,j)
     vy=vi*B%Aky(i,j)+vj*B%Aiy(i,j)
     Ty=Ti*B%Aky(i,j)+Tj*B%Aiy(i,j)
     t11=(4.d0/3.d0*ux-2.d0/3.d0*vy)*mu0
     t22=(4.d0/3.d0*vy-2.d0/3.d0*ux)*mu0
     t12=(uy+vx)*mu0
     E1=uu(i,j)*t11+v(i,j)*t12+Ak0*Tx
     E2=uu(i,j)*t12+v(i,j)*t22+Ak0*Ty

 	 Ev1(i,j,1)=(B%Akx(i,j)*t11+B%Aky(i,j)*t12)/B%Ajac(i,j)
     Ev1(i,j,2)=(B%Akx(i,j)*t12+B%Aky(i,j)*t22)/B%Ajac(i,j)
	 Ev1(i,j,3)=(B%Akx(i,j)*E1+B%Aky(i,j)*E2  )/B%Ajac(i,j)
 
 	 Ev2(i,j,1)=(B%Aix(i,j)*t11+B%Aiy(i,j)*t12)/B%Ajac(i,j)
     Ev2(i,j,2)=(B%Aix(i,j)*t12+B%Aiy(i,j)*t22)/B%Ajac(i,j)
	 Ev2(i,j,3)=(B%Aix(i,j)*E1+ B%Aiy(i,j)*E2 )/B%Ajac(i,j)
   enddo
   enddo
!$OMP END DO

!  计算粘性通量

!$OMP DO
   do j=3,ny-2
   do i=3,nx-2
    Fluxi(2,i,j)=Fluxi(2,i,j)   +0.5*(Ev1(i,j,1)+Ev1(i+1,j,1))   
    Fluxi(3,i,j)=Fluxi(3,i,j)   +0.5*(Ev1(i,j,2)+Ev1(i+1,j,2))
    Fluxi(4,i,j)=Fluxi(4,i,j)   +0.5*(Ev1(i,j,3)+Ev1(i+1,j,3))

    Fluxj(2,i,j)=Fluxj(2,i,j)   +0.5*(Ev2(i,j,1)+Ev2(i,j+1,1))   
    Fluxj(3,i,j)=Fluxj(3,i,j)   +0.5*(Ev2(i,j,2)+Ev2(i,j+1,2))
    Fluxj(4,i,j)=Fluxj(4,i,j)   +0.5*(Ev2(i,j,3)+Ev2(i,j+1,3))
   enddo
   enddo
!$OMP END DO
  endif

!$OMP DO
    do j=LFDM+1,ny-LFDM
    do i=LFDM+1,nx-LFDM
    do m=1,4
      B%Res(m,i,j)= (Fluxi(m,i,j)-Fluxi(m,i-1,j)+Fluxj(m,i,j)-Fluxj(m,i,j-1))
	enddo
	enddo
	enddo
!$OMP END DO
!$OMP END PARALLEL
    deallocate(EV1,EV2)


  end  

!---------------------------------------------------------------------------------------

    subroutine split_Steger_Warming(d,u,v,cc,fP,fm,A1,A2,gamma)
	  implicit none
	  integer::i,j,i1,i2,j1,j2
	  
      real*8:: d,u,v,cc,A1,A2
      real*8,dimension (4):: fp,fm
	  real*8:: gamma,Ak1,AK2,tmp0,tmp1,tmp2,tmp3,ss,E1,E2,E3,E1P,E2P,E3P,E1M,E2M,E3M,vs,uc1,uc2,vc1,vc2,vvc1,vvc2,vv,W2 

!c El 为特征值,其中El(:,1)为x方向的特征值 (4个， u, u, u+c, u-c)
	  tmp1=2.d0*(gamma-1.d0)
	  tmp2=1.d0/(2.d0*gamma)
	  tmp3=(3.d0-gamma)/(2.d0*(gamma-1.d0)) 
	  
        vs=A1*u+A2*v  
        ss=sqrt(A1*A1+A2*A2)
        ak1=A1/ss
	    ak2=A2/ss
!c       E1 is lamda1, lamda2 , E2 is lamda3; E3 is lamda4 
	    E1=vs
        E2=vs-cc*ss
        E3=vs+cc*ss
	    E1P=(E1+abs(E1))*0.5d0
	    E2P=(E2+abs(E2))*0.5d0
	    E3P=(E3+abs(E3))*0.5d0
	    E1M=E1-E1P
	    E2M=E2-E2P
	    E3M=E3-E3P
        tmp0=d/(2.d0*gamma) 
	    uc1=u-cc*ak1
	    uc2=u+cc*ak1
	    vc1=v-cc*ak2
	    vc2=v+cc*ak2
	    vvc1=(uc1*uc1+vc1*vc1)/2.d0
	    vvc2=(uc2*uc2+vc2*vc2)/2.d0
	    vv=(gamma-1.d0)*(u*u+v*v)
        W2=tmp3*cc*cc

        fp(1)=tmp0*(tmp1*E1P+E2P+E3P)
        fp(2)=tmp0*(tmp1*E1P*u+E2P*uc1+E3P*uc2)
	    fp(3)=tmp0*(tmp1*E1P*v+E2P*vc1+E3P*vc2)
        fp(4)=tmp0*(E1P*vv+E2p*vvc1+E3P*vvc2+W2*(E2P+E3P))
        
        fm(1)=tmp0*(tmp1*E1M+E2M+E3M)
        fm(2)=tmp0*(tmp1*E1M*u+E2M*uc1+E3M*uc2)
	    fm(3)=tmp0*(tmp1*E1M*v+E2M*vc1+E3M*vc2)
        fm(4)=tmp0*(E1M*vv+E2M*vvc1+E3M*vvc2+W2*(E2M+E3M))

     return 
     end


!c---------------------------------------------------------------
   
       subroutine fp_weno5(nx,v,hj)
       implicit none
       integer:: nx,k,LAP
	   real*8:: v(nx),hj(nx)
       real*8:: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23
   	   ep=1.d-6
 	   C03=3.d0/10.d0
   	   C13=3.d0/5.d0
	   C23=1.d0/10.d0

      do k=3,nx-2
         S0=13.d0/12.d0*(v(k)-2.d0*v(k+1)+v(k+2))**2+  1.d0/4.d0*(3.d0*v(k)-4.d0*v(k+1)+v(k+2))**2
         S1=13.d0/12.d0*(v(k-1)-2.d0*v(k)+v(k+1))**2+  1.d0/4.d0*(v(k-1)-v(k+1))**2
         S2=13.d0/12.d0*(v(k-2)-2.d0*v(k-1)+v(k))**2+  1.d0/4.d0*(v(k-2)-4.d0*v(k-1)+3.*v(k))**2

      a0=C03/((ep+S0)**2)
	  a1=C13/((ep+S1)**2)
	  a2=C23/((ep+S2)**2)

	  W0=a0/(a0+a1+a2)
      W1=a1/(a0+a1+a2)
	  W2=a2/(a0+a1+a2)

	  q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k+1)-1.d0/6.d0*v(k+2)
	  q13=-1.d0/6.d0*v(k-1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k+1)
	  q23=1.d0/3.d0*v(k-2)-7.d0/6.d0*v(k-1)+11.d0/6.d0*v(k)
	  hj(k)=W0*q03+W1*q13+W2*q23
     enddo

	 return
	 end
!-------------------------------------------------
!----- 对于负通量：
       subroutine fm_weno5(nx,v,hj)
       implicit none
       integer:: nx,k,LAP
	   real*8:: v(nx),hj(nx)
       real*8:: ep,C03,C13,C23,S0,S1,S2,a0,a1,a2,W0,W1,W2,q03,q13,q23

	 ep=1.d-6
	 
	 C03=3.d0/10.d0
   	 C13=3.d0/5.d0
	 C23=1.d0/10.d0

       do k=3,nx-2

         S0=13.d0/12.d0*(v(k)-2.d0*v(k-1)+v(k-2))**2+ 1.d0/4.d0*(3.d0*v(k)-4.d0*v(k-1)+v(k-2))**2
         S1=13.d0/12.d0*(v(k+1)-2.d0*v(k)+v(k-1))**2+ 1.d0/4.d0*(v(k+1)-v(k-1))**2
         S2=13.d0/12.d0*(v(k+2)-2.d0*v(k+1)+v(k))**2+ 1.d0/4.d0*(v(k+2)-4.d0*v(k+1)+3.d0*v(k))**2

       a0=C03/((ep+S0)**2)
	   a1=C13/((ep+S1)**2)
	   a2=C23/((ep+S2)**2)

	  W0=a0/(a0+a1+a2)
      W1=a1/(a0+a1+a2)
	  W2=a2/(a0+a1+a2)

	 q03=1.d0/3.d0*v(k)+5.d0/6.d0*v(k-1)-1.d0/6.d0*v(k-2)
	 q13=-1.d0/6.d0*v(k+1)+5.d0/6.d0*v(k)+1.d0/3.d0*v(k-1)
	 q23=1.d0/3.d0*v(k+2)-7.d0/6.d0*v(k+1)+11.d0/6.d0*v(k)

	 hj(k-1)=W0*q03+W1*q13+W2*q23
	 enddo

     return

	 end	  

!  有限差分模块的初始化
   subroutine Init_FiniteDifference
   use Global_var
   implicit none
   integer:: NB_FDM,m,md
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   logical EX
    MP=> Mesh(1)   ! 最密的网格

  Inquire(file="FDM.in",exist=EX)
  if(EX) then
    print*, "Find 'FDM.in',  Blocks in FDM.in using Finite Difference Method "
	open(99,file="FDM.in")
	read(99,*)
	read(99,*)
	read(99,*)
	read(99,*)
    read(99,*) NB_FDM
	read(99,*) 
	do m=1,NB_FDM
	read(99,*) md
	if(md .gt. MP%Num_Block) then
	  print*, "error in 'FDM.in' "
	else
	 MP%Block(md)%FVM_FDM=Method_FDM
	endif
	enddo
    close(99)
   else
     print*, "Not find 'FDM.in', all blocks using Finite Volume Method"
   endif

  call Comput_Jacobian
  end
      
!c---------------------------------------------------------------
 subroutine comput_Jacobian
  use Global_var
  implicit none
  Type (Block_TYPE),pointer:: B
  real*8,allocatable,dimension(:,:):: xx,yy,xk,xi, yk,yi
  integer:: i,j,m,nMesh,NB,nx,ny
  real*8:: hx,hy
 
    open(99,file="Jocabian.dat")
    write(99,*) "variables=x,y,kx,ky,ix,iy,Jac,s1" 
	do m=1,Mesh(1)%Num_Block
     B => Mesh(1)%Block(m)
     if(B%FVM_FDM .ne. METHOD_FDM) cycle  
     print*, " Comput Jocabian Coefficient, Block No.", m
     
	 nx=B%nx-1
	 ny=B%ny-1

!   申请内存  (x1,y1)网格中心坐标
!   差分法使用的Jocabian系数
    allocate(B%Akx(nx,ny),B%Aky(nx,ny), B%Aix(nx,ny),B%Aiy(nx,ny),B%Ajac(nx,ny))
    allocate(xx(nx,ny),yy(nx,ny),xk(nx,ny), xi(nx,ny),  yk(nx,ny),yi(nx,ny))
    
	 hx=1.d0   ! 
	 hy=1.d0   !
   
     do j=1,ny
  	 do i=1,nx
      xx(i,j)=B%x1(i,j)
 	  yy(i,j)=B%y1(i,j)
	 enddo
	 enddo
	  

	 call dx0(xx,xk,nx,ny,hx)
	 call dx0(yy,yk,nx,ny,hx)
	 call dy0(xx,xi,nx,ny,hy)
	 call dy0(yy,yi,nx,ny,hy)
	 
	 do j=1,ny
	 do i=1,nx
     B%Ajac(i,j)=xk(i,j)*yi(i,j)-xi(i,j)*yk(i,j)
     B%Ajac(i,j)=1.d0/B%Ajac(i,j)
	 B%Akx(i,j)=B%Ajac(i,j)*yi(i,j)
	 B%Aky(i,j)=-B%Ajac(i,j)*xi(i,j)
	 B%Aix(i,j)=-B%Ajac(i,j)*yk(i,j)
	 B%Aiy(i,j)=B%Ajac(i,j)*xk(i,j)
     enddo
	 enddo

	 do j=1,ny
	 do i=1,nx
	   if(B%Ajac(i,j).lt. 0 ) print*, i,j,B%Ajac(i,j)
	 enddo
	 enddo


	 write(99,*) "zone i= ", nx, " j= ", ny
	 do j=1,ny
	 do i=1,nx
	    write(99,"(8f20.10)") xx(i,j),yy(i,j),B%Akx(i,j),B%Aky(i,j),B%Aix(i,j),B%Aiy(i,j),B%Ajac(i,j),B%vol(i,j)*B%Ajac(i,j)
	 enddo
	 enddo
     deallocate(xk,xi,yk,yi,xx,yy)
	enddo
	 print*, "comput Jocabian OK ..."
     close(99)

	end   subroutine comput_Jacobian



!c==================================================================
       subroutine dx0(f,fx,nx,ny,hx)
        implicit none
		integer:: nx,ny,i,j
		real*8:: hx,b1,b2,a1,a2,a3
        real*8:: f(nx,ny),fx(nx,ny)
         b1=8.d0/(12.d0*hx)
         b2=1.d0/(12.d0*hx)
         a1=1.d0/(60.d0*hx)
         a2=-3.d0/(20.d0*hx)
         a3=3.d0/(4.d0*hx)

         do j=1,ny
         do i=4,nx-3
          fx(i,j)  =a1*(f(i+3,j)-f(i-3,j)) +a2*(f(i+2,j)-f(i-2,j)) +a3*(f(i+1,j)-f(i-1,j))
         enddo
         enddo

	    do j=1,ny 
          fx(1,j)=(-3.d0*f(1,j)+4.d0*f(2,j)-f(3,j))/(2.d0*hx)            
          fx(2,j)=(-2.d0*f(1,j)-3.d0*f(2,j)+6.d0*f(3,j)-f(4,j)) /(6.d0*hx)  
	      fx(3,j)=b1*(f(4,j)-f(2,j)) -b2*(f(5,j)-f(1,j))
          fx(nx-2,j)=b1*(f(nx-1,j)-f(nx-3,j)) -b2*(f(nx,j)-f(nx-4,j))
          fx(nx-1,j)=(f(nx-3,j)-6.d0*f(nx-2,j)+3.d0*f(nx-1,j) +2.d0*f(nx,j))/(6.d0*hx)
          fx(nx,j)=(f(nx-2,j)-4.d0*f(nx-1,j)+3.d0*f(nx,j))/(2.d0*hx)
         enddo
       end

!----------------------------------------------
       subroutine dy0(f,fy,nx,ny,hy)
       implicit none
	   integer:: nx,ny,i,j
	   real*8:: hy,b1,b2,a1,a2,a3
       real*8:: f(nx,ny),fy(nx,ny)
	     b1=8.d0/(12.d0*hy)
	     b2=1.d0/(12.d0*hy)
         a1=1.d0/(60.d0*hy)
         a2=-3.d0/(20.d0*hy)
         a3=3.d0/(4.d0*hy)

         do j=4,ny-3
         do i=1,nx
          fy(i,j)=a1*(f(i,j+3)-f(i,j-3)) +a2*(f(i,j+2)-f(i,j-2)) +a3*(f(i,j+1)-f(i,j-1))
         enddo
         enddo
	     do i=1,nx 
          fy(i,1)=(-3.d0*f(i,1)+4.d0*f(i,2)-f(i,3))/(2.d0*hy)            
          fy(i,2)=(-2.d0*f(i,1)-3.d0*f(i,2)+6.d0*f(i,3)-f(i,4)) /(6.d0*hy)  
	      fy(i,3)=b1*(f(i,4)-f(i,2)) -b2*(f(i,5)-f(i,1))
          fy(i,ny-2)=b1*(f(i,ny-1)-f(i,ny-3))  -b2*(f(i,ny)-f(i,ny-4))
          fy(i,ny-1)=(f(i,ny-3)-6.d0*f(i,ny-2)+3.d0*f(i,ny-1)  +2.d0*f(i,ny))/(6.d0*hy)
          fy(i,ny)=(f(i,ny-2)-4.d0*f(i,ny-1)+3.d0*f(i,ny))/(2.d0*hy)
         enddo
        end


