!---------------------------------------------------------------------------------------
! 湍流模型模块，For OpenCFD-EC2D ver 1.1; Copyright by Li Xinliang, lixl@imech.ac.cn
! Code by Li Xinliang
!--------------------------------------------------------------------------------------
!  BL模型
   subroutine turbulence_model_BL(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   integer:: nMesh,mBlock,nx1,ny1,ksub,i,j,kflag,i1,j1
   real*8:: ui,vi,uj,vj,ux,vx,uy,vy
   real*8:: xi,yi,xj,yj,Ds,Jac,ix,iy,jx,jy,x0,y0
   real*8,allocatable,dimension(:,:):: omiga
   real*8,allocatable,dimension(:):: Amu1d,Amut1d,d1d,u1d,v1d,yy,omiga1d
   integer,allocatable:: flag1(:,:) 
   Type (Block_TYPE),pointer:: B
   Type (BC_MSG_TYPE),pointer:: Bc
   
   B => Mesh(nMesh)%Block(mBlock)
   nx1=B%nx ; ny1=B%ny
  

! test if the block contains wall    
   kflag=0
   do ksub=1, B%subface
     if(B%bc_msg(ksub)%neighb .eq. BC_WALL) kflag=1
   enddo  
   if(Kflag .eq. 0) return    ! No wall in this block

!  This Block Contains Wall
   allocate(omiga(0:nx1,0:ny1),flag1(0:nx1,0:ny1))
   B%Amu_t(:,:)=0.d0
   omiga=0.d0
   flag1=0

!----- get Omiga (vorticity)  omiga=vx-uy at the cell's center ------
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nx1,ny1,B,uu,v,omiga)
   do j=1,ny1-1
   do i=1,nx1-1
    xi=B%x1(i+1,j)-B%x1(i-1,j)
    yi=B%y1(i+1,j)-B%y1(i-1,j)
    xj=B%x1(i,j+1)-B%x1(i,j-1)
    yj=B%y1(i,j+1)-B%y1(i,j-1)
    Jac=1.d0/(xi*yj-xj*yi)
    ix=Jac*yj; iy=-Jac*xj; jx=-Jac*yi; jy=Jac*xi   ! Jocabian:  ksix,ksiy,itax,itay

    ui=uu(i+1,j)-uu(i-1,j)
    vi=v(i+1,j)-v(i-1,j)
    uj=uu(i,j+1)-uu(i,j-1)
    vj=v(i,j+1)-v(i,j-1)

    vx=vi*ix+vj*jx
    uy=ui*iy+uj*jy
    omiga(i,j)=abs(vx-uy)
   enddo
   enddo
!$OMP END PARALLEL DO
!----------------------------------------------------------------------
 do ksub=1, B%subface
  Bc=> B%bc_msg(ksub)
  if(Bc%neighb .eq. BC_WALL) then   ! Wall boundary
  
  if(Bc%face .eq. 1) then   ! face of i+ 
    allocate(yy(nx1),Amu1d(nx1),Amut1d(nx1),d1d(nx1),u1d(nx1),v1d(nx1),omiga1d(nx1)) 
	Amu1d=0.d0; Amut1d=0.d0
	 
    do j=Bc%jst,Bc%jend-1
     x0=(B%x1(1,j)+B%x1(0,j))*0.5d0 ; y0=(B%y1(1,j)+B%y1(0,j))*0.5d0
     do i=1,nx1-1
       yy(i)=sqrt((B%x1(i,j)-x0)**2+(B%y1(i,j)-y0)**2)
       u1d(i)=uu(i,j) ;    v1d(i)=v(i,j) ;   d1d(i)=d(i,j)  ;   omiga1d(i)=omiga(i,j)
       Amu1d(i)=B%Amu(i,j)
     enddo
   
     call BL_model_1d(nx1,yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
   
     do i=1,nx1-1
     if(flag1(i,j) .eq. 0) then
      flag1(i,j)=1
      B%Amu_t(i,j)=Amut1d(i)
     else
      B%Amu_t(i,j)=min(B%Amu_t(i,j),Amut1d(i))
     endif
     enddo
    enddo

     B%Amu_t(:,Bc%jst-1)=B%Amu_t(:,Bc%jst) ; B%Amu_t(:,Bc%jend)=B%Amu_t(:,Bc%jend-1)
     B%Amu_t(nx1,Bc%jst:Bc%jend)=B%Amu_t(nx1-1,Bc%jst:Bc%jend)
     B%Amu_t(0,Bc%jst:Bc%jend)=-B%Amu_t(1,Bc%jst:Bc%jend)   ! To set Amu_t=0 in the wall
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
  
  else if (Bc%face .eq. 3) then  ! face i-
    allocate(yy(nx1),Amu1d(nx1),Amut1d(nx1),d1d(nx1),u1d(nx1),v1d(nx1),omiga1d(nx1))  
   	Amu1d=0.d0; Amut1d=0.d0

	do j=Bc%jst,Bc%jend-1
     x0=(B%x1(nx1-1,j)+B%x1(nx1,j))*0.5d0 ; y0=(B%y1(nx1-1,j)+B%y1(nx1,j))*0.5d0
     do i=nx1-1,1,-1
       i1=nx1-i
       yy(i1)=sqrt((B%x1(i,j)-x0)**2+(B%y1(i,j)-y0)**2)
       u1d(i1)=uu(i,j) ;    v1d(i1)=v(i,j) ;   d1d(i1)=d(i,j)  ;   omiga1d(i1)=omiga(i,j)
       Amu1d(i1)=B%Amu(i,j)
     enddo
     call BL_model_1d(nx1,yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)

     do i=nx1-1,1,-1
       i1=nx1-i
       if(flag1(i,j) .eq. 0) then
         flag1(i,j)=1
         B%Amu_t(i,j)=Amut1d(i1)
       else
        B%Amu_t(i,j)=min(B%Amu_t(i,j),Amut1d(i1))
      endif
     enddo
    enddo

     B%Amu_t(:,Bc%jst-1)=B%Amu_t(:,Bc%jst) ; B%Amu_t(:,Bc%jend)=B%Amu_t(:,Bc%jend-1)
     B%Amu_t(0,Bc%jst:Bc%jend)=B%Amu_t(1,Bc%jst:Bc%jend)
     B%Amu_t(nx1,Bc%jst:Bc%jend)=-B%Amu_t(nx1-1,Bc%jst:Bc%jend)   ! To set B%Amu_t=0 in the wall
    deallocate(yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
 
  else if(Bc%face .eq. 2) then   ! face of j- 
  
    allocate(yy(ny1),Amu1d(ny1),Amut1d(ny1),d1d(ny1),u1d(ny1),v1d(ny1),omiga1d(ny1))  
  	Amu1d=0.d0; Amut1d=0.d0

    do i=Bc%ist,Bc%iend-1
     x0=(B%x1(i,1)+B%x1(i,0))*0.5d0 ; y0=(B%y1(i,1)+B%y1(i,0))*0.5d0
     do j=1,ny1-1
       yy(j)=sqrt((B%x1(i,j)-x0)**2+(B%y1(i,j)-y0)**2)
       u1d(j)=uu(i,j) ;    v1d(j)=v(i,j) ;   d1d(j)=d(i,j)  ;   omiga1d(j)=omiga(i,j)
       Amu1d(j)=B%Amu(i,j)
     enddo
   
     call BL_model_1d(ny1,yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
   
     do j=1,ny1-1
     if(flag1(i,j) .eq. 0) then
      flag1(i,j)=1
      B%Amu_t(i,j)=Amut1d(j)
     else
      B%Amu_t(i,j)=min(B%Amu_t(i,j),Amut1d(j))
     endif
     enddo
    enddo

     B%Amu_t(Bc%ist:Bc%iend,0)=-B%Amu_t(Bc%ist:Bc%iend,1) ! To set B%Amu_t=0 in the wall
     B%Amu_t(Bc%ist:Bc%iend,ny1)=B%Amu_t(Bc%ist:Bc%iend,ny1-1)
     B%Amu_t(Bc%ist-1,:)=B%Amu_t(Bc%ist,:)
     B%Amu_t(Bc%iend,:)=B%Amu_t(Bc%iend-1,:)   
     deallocate(yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
  else   ! face j-
    allocate(yy(ny1),Amu1d(ny1),Amut1d(ny1),d1d(ny1),u1d(ny1),v1d(ny1),omiga1d(ny1))  
   	Amu1d=0.d0; Amut1d=0.d0

	do i=Bc%ist,Bc%iend-1
     x0=(B%x1(i,ny1-1)+B%x1(i,ny1))*0.5d0 ; y0=(B%y1(i,ny1-1)+B%y1(i,ny1))*0.5d0
     do j=ny1-1,1,-1
       j1=ny1-j
       yy(j1)=sqrt((B%x1(i,j)-x0)**2+(B%y1(i,j)-y0)**2)
       u1d(j1)=uu(i,j) ;    v1d(j1)=v(i,j) ;   d1d(j1)=d(i,j)  ;   omiga1d(j1)=omiga(i,j)
       Amu1d(j1)=B%Amu(i,j)
     enddo
   
     call BL_model_1d(ny1,yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
   
     do j=ny1-1,1,-1
     j1=ny1-j
     if(flag1(i,j) .eq. 0) then
      flag1(i,j)=1
      B%Amu_t(i,j)=Amut1d(j1)
     else
      B%Amu_t(i,j)=min(B%Amu_t(i,j),Amut1d(j1))
     endif
     enddo
    enddo

     B%Amu_t(Bc%ist:Bc%iend,ny1)=-B%Amu_t(Bc%ist:Bc%iend,ny1-1) ! To set B%Amu_t=0 in the wall
     B%Amu_t(Bc%ist:Bc%iend,0)=B%Amu_t(Bc%ist:Bc%iend,1)
     B%Amu_t(Bc%ist-1,:)=B%Amu_t(Bc%ist,:)
     B%Amu_t(Bc%iend,:)=B%Amu_t(Bc%iend-1,:)   
     deallocate(yy,Amu1d,Amut1d,d1d,u1d,v1d,omiga1d)
   endif
  endif
 enddo

! 取消限制条件(2010-10-12). 很多工况出现粘性不足，产生虚假分离，无法收敛现象，需要加大粘性。
!   do j=0,ny1
!   do i=0,nx1
!        限制条件，避免计算出的湍流粘性系数过大
!   if(Amu_t(i,j) .gt. 1.d0/sqrt(Re)) Amu_t(i,j)=1.d0/sqrt(Re)   
!   enddo
!   enddo

 
 !  write(99,*) "zone i=", nx1+1, " j= ", ny1+1
 !  do j=0, ny1
 !  do i=0, nx1
 !  write(99,"(4f20.10)") B%x1(i,j),B%y1(i,j), omiga(i,j),Amu_t(i,j)*Re
 !  enddo
 !  enddo

   deallocate(omiga,flag1)

end


!c------------------------------------------------------------------------
! B-L model of turbulence
! Ref:  Wilox DC. Turbulence Modeling for CFD (2nd Edition), p77
   subroutine BL_model_1d(ny,yy,Amu,Amu_t,d,u,v,omiga)
      implicit none
      integer ny,j,Iflag,Ny_boundary
      real*8,dimension(ny) :: yy, Amu,Amu_t,d,u,v,omiga
      real*8,parameter::  AP=26.,Ccp=1.6,Ckleb=0.3,Cwk=0.25d0,AKT=0.4,AK=0.0168  ! Cwk= 1.d0
      real*8:: Tw,Ret,Fmax,etamax,Udif,etap,FF,Fwak,bl,Fkleb,Visti,Visto,uu
      
           TW=abs(Amu(1)*omiga(1))
            do j=1,Ny-1
             if(abs(Amu(j)*omiga(j)) .gt. TW) Tw=abs(Amu(j)*omiga(j))  
            enddo
             Ret=sqrt(d(1)*TW)/Amu(1)

           Fmax=0.d0 ;   etamax=0.d0 ;     Udif=0.d0

!            Ny_boundary=Ny
             Ny_boundary=Ny-1   !!! 修改 2011-10-18

 ! 临时修改，防止搜索范围过大，造成ymax,Fmax 出现过大值          
 !          if(Ny .gt. 10) then
 !           Ny_boundary=Ny-5
 !          else
 !           Ny_boundary=Ny
 !          endif   
!----------------------------------
           do j=1,Ny_boundary
            uu=sqrt(u(j)**2+v(j)**2)
            if(uu .gt. Udif) Udif=uu
             etap=yy(j)*Ret
             FF=yy(j)*abs(omiga(j))*(1.d0-exp(-etap/AP))

            if(FF.gt.Fmax) then
             Fmax=FF
             etamax=yy(j)    ! 某些位置算出的值偏大
            endif

!            if(FF.gt.Fmax) then
!             Fmax=FF
!             etamax=yy(j)    
!           else
!            goto 100    ! Find the first peak of F(y)   ! 重要的修改 !!!         
!           endif

           enddo
100        continue

           IFlag=0
           Fwak=dmin1(etamax*Fmax,Cwk*etamax*Udif*Udif/Fmax)
           do j=1,Ny-1
            etap=Ret*yy(j)
            bl=AKT*yy(j)*(1.d0-exp(-etap/AP))
            visti=d(j)*bl*bl*abs(omiga(j))
            Fkleb=1.d0/(1.d0+5.5d0*(Ckleb*yy(j)/etamax)**6)
            visto=AK*Ccp*d(j)*Fwak*Fkleb
            if(abs(visto).lt.abs(visti)) IFlag=1
            if(Iflag.eq.0) then
             Amu_t(j)=visti
            else
             Amu_t(j)=visto
            endif
           enddo 

  end


