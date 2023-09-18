
! SST k-wģ��
! ��������ճ��ϵ��Amu_t;
! ����k���̼�w���̵Ĳв�Ҷ��; 

 subroutine SST_kw(nMesh,mBlock)  
   Use Global_Var
   Use Flow_Var
   implicit none
   
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   integer:: nMesh,mBlock,nx,ny,i,j
   real*8:: si,sj,ni1,ni2,nj1,nj2,un,un1,un2,an1,an2,Amu1,Amu2
   real*8:: Dix,Diy,Djx,Djy,Ds,Dik,Diw,Djk,Djw,Diu,Div,Dju,Djv,Kx,Ky,Wx,Wy,ux,vx,uy,vy, &
             omega,Kws,CD_kw,arg1,arg2,arg3,f2,t11,t22,t12,Pk,Pk0
             
   real*8,allocatable,dimension(:,:):: Kt,Wt,f1,Qk,Qw,Fluxk,Fluxw
   
!  ģ��ϵ�� (1Ϊk-wģ�͵�ϵ�����ڽ�����ʹ�ã� 2Ϊk-epslģ�͵�ϵ����Զ����ʹ�ã�
!  SST��һ��k-w��k-epslģ�͵Ļ��ģ�ͣ�ͨ�����غ���f���л�
   real*8,parameter:: sigma_k1_SST=0.85d0,sigma_w1_SST=0.5d0,beta1_SST=0.075d0,Cw1_SST=0.533d0, &
                        sigma_k2_SST=1.d0,  sigma_w2_SST=0.856d0,beta2_SST=0.0828d0,Cw2_SST=0.440d0
   real*8,parameter::a1_SST=0.31d0, betas_SST=0.09d0
   real*8:: sigma_k_SST,sigma_w_SST,beta_SST,Cw_SST 
    
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx ; ny=B%ny
   allocate(f1(nx,ny),Qk(nx,ny),Qw(nx,ny))
   allocate(Kt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),Wt(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
   allocate(Fluxk(nx,ny),Fluxw(nx,ny))

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nx,ny,B,Kt,Wt,f1,Qk,Qw,Fluxk,Fluxw,d,uu,v,Re)

!$OMP DO   	
    do j=1-LAP,ny+LAP-1
    do i=1-LAP,nx+LAP-1
      Kt(i,j)=B%U(5,i,j)/d(i,j)
      Wt(i,j)=B%U(6,i,j)/d(i,j)
    enddo
    enddo
    
!$OMP END DO
f1=0.
!  ��������ճ��ϵ��Amu_t;   
!  ���㺯��f1 �����ֽ�������Զ������ 
!  ����Դ���������2�����Ĳ�֣�����Jocabian�任��
!  Դ����k��w���̵�������ȶ�������ɢ����ؼ�)
! See Blazek's Book, section 7.2.3

!$OMP DO   	
    do j=1,ny-1
    do i=1,nx-1

! ����������������x,y�ĵ�����ʹ��Jocabian�任
     Dix=(B%x1(i+1,j)-B%x1(i-1,j))*0.5d0
     Diy=(B%y1(i+1,j)-B%y1(i-1,j))*0.5d0
     Djx=(B%x1(i,j+1)-B%x1(i,j-1))*0.5d0
     Djy=(B%y1(i,j+1)-B%y1(i,j-1))*0.5d0
     Ds=1.d0/(Dix*Djy-Djx*Diy)

     Diu=(uu(i+1,j)-uu(i-1,j))*0.5d0
     Div=(v(i+1,j)-v(i-1,j))*0.5d0
     Dik=(Kt(i+1,j)-Kt(i-1,j))*0.5d0
     Diw=(Wt(i+1,j)-Wt(i-1,j))*0.5d0

     Dju=(uu(i,j+1)-uu(i,j-1))*0.5d0
     Djv=(v(i,j+1)-v(i,j-1))*0.5d0
     Djk=(Kt(i,j+1)-Kt(i,j-1))*0.5d0
     Djw=(Wt(i,j+1)-Wt(i,j-1))*0.5d0

! ����ֵ    
     ux=(Diu*Djy-Dju*Diy)*Ds
     vx=(Div*Djy-Djv*Diy)*Ds
     kx=(Dik*Djy-Djk*Diy)*Ds
     wx=(Diw*Djy-Djw*Diy)*Ds

     uy=(-Diu*Djx+Dju*Dix)*Ds
     vy=(-Div*Djx+Djv*Dix)*Ds
     ky=(-Dik*Djx+Djk*Dix)*Ds
     wy=(-Diw*Djx+Djw*Dix)*Ds
 
 !  ��������ճ��ϵ�� �� Blazek's Book Eq. (7.66)    
     B%Amu(i,j)=B%Amu(i,j)*Re
     omega=vx-uy      ! ����
	 arg2=max( 2.d0* sqrt(abs(Kt(i,j)))/(0.09*Wt(i,j)*B%dw(i,j)*Re) , 500.d0*B%Amu(i,j)/(d(i,j)*Wt(i,j)*B%dw(i,j)**2*Re**2) )
   
     f2=tanh(arg2**2)

!!!  Revised by Wang XiangYu
     B%Amu_t(i,j)=min(min(d(i,j)*Kt(i,j)/Wt(i,j),a1_SST*d(i,j)*Kt(i,j)*Re/(f2*abs(omega))),100000.)
 ! ����f1 (ʶ���Ƿ�Ϊ��������������������1��      
     Kws=2.d0*(kx*wx+ky*wy)*d(i,j)*sigma_w2_SST/(Wt(i,j)+1.d-20)      ! ����������
     CD_kw=max(Kws,1.d-20)
     arg3=max(sqrt(abs(Kt(i,j)))/(0.09*Wt(i,j)*B%dw(i,j)*Re)  , 500.d0*B%Amu(i,j)/(d(i,j)*Wt(i,j)*B%dw(i,j)**2*Re**2) )
    
	 arg1=min(arg3,4.d0*d(i,j)*sigma_w2_SST*Kt(i,j)/(CD_kw*B%dw(i,j)**2 ))
     f1(i,j)=tanh(arg1**4)             ! ���غ�����������������1��Զ����������0  �������л�k-w��k-epsl����)
     
     
!    ��Ӧ�� ��ʹ������ճģ�ͣ�     
     t11=((4.d0/3.d0)*ux-(2.d0/3.d0)*vy)*B%Amu_t(i,j) -(2.d0/3.d0)*d(i,j)*Kt(i,j)*Re   ! Blazek's Book, Eq. (7.25)
     t22=((4.d0/3.d0)*vy-(2.d0/3.d0)*ux)*B%Amu_t(i,j) -(2.d0/3.d0)*d(i,j)*Kt(i,j)*Re
     t12=(uy+vx)*B%Amu_t(i,j)

!    ���ܷ��̵�Դ�����-��ɢ)     
!     Pk=t11*ux+t22*vy+t12*(uy+vx)  ! ���������� ����Ӧ������Ӧ���ʣ�       
      Pk=B%Amu_t(i,j)*omega*omega        ! ��   
	
	 Pk0=min(Pk,20.d0*betas_SST*Kt(i,j)*Wt(i,j)*Re**2)    ! ������������������ƣ���ֹ���ܹ���

	 Qk(i,j)=Pk0/Re-betas_SST*d(i,j)*Wt(i,j)*Kt(i,j)*Re    ! k���̵�Դ��  ��������-��ɢ�

     Cw_SST=f1(i,j)*Cw1_SST+(1.d0-f1(i,j))*Cw2_SST    ! ģ��ϵ��������f1���������л�
     beta_SST=f1(i,j)*beta1_SST+(1.d0-f1(i,j))*beta2_SST    ! ģ��ϵ��������f1���������л�
!     Qw(i,j)= Cw_SST*d(i,j)*Pk/(B%Amu_t(i,j)+1.d-20)/Re-beta_SST*d(i,j)*Wt(i,j)**2*Re+(1.d0-f1(i,j))*Kws/Re     ! W���̵�Դ��    
      Qw(i,j)= Cw_SST*d(i,j)*omega*omega/Re-beta_SST*d(i,j)*Wt(i,j)**2*Re+(1.d0-f1(i,j))*Kws/Re     ! W���̵�Դ��    
	 
	enddo
	enddo
!$OMP END DO
 

 
! ���������ɢ��
!----- i- direction ----------------------------------------------------------------------------------
!$OMP  DO
   do j=1,ny-1
   do i=1,nx
	 si=B%si(i,j)   ! ������߽糤��
	 ni1=B%ni1(i,j); ni2=B%ni2(i,j)   ! ������

 !   ���������1��ӭ���ʽ      
     un1=uu(i-1,j)*ni1+v(i-1,j)*ni2
	 un2=uu(i,j)*ni1+v(i,j)*ni2

!   1��L-F ��ʽ
	 fluxk(i,j)=-0.5d0*( (un1+abs(un1))*Kt(i-1,j) + (un2-abs(un2))*Kt(i,j) )*si
	 fluxw(i,j)=-0.5d0*( (un1+abs(un1))*Wt(i-1,j) + (un2-abs(un2))*Wt(i,j) )*si


 !   ճ�����ɢ�������2�����ĸ�ʽ

 !   ��ʽϵ����k-w��k-epsl��ʽϵ��֮��ѡ�� (f1��Ϊ�л����غ���)
     sigma_K_SST=f1(i,j)*sigma_k1_SST+(1.d0-f1(i,j))*sigma_k2_SST
     sigma_W_SST=f1(i,j)*sigma_w1_SST+(1.d0-f1(i,j))*sigma_w2_SST
     
 ! �����ϵ�ֵ=����ֵ��ƽ��, �߽��ϵ���ɢϵ��=�ڲ��ֵ   
     if(i .eq. 1) then
       Amu1=B%Amu(i,j) + sigma_K_SST*B%Amu_t(i,j)         ! ��ɢϵ�� (k����)
       Amu2=B%Amu(i,j) + sigma_W_SST*B%Amu_t(i,j)         ! ��ɢϵ�� (w����)
     else if (i .eq. nx ) then  
       Amu1=B%Amu(i-1,j) + sigma_K_SST*B%Amu_t(i-1,j)        ! ��ɢϵ�� (k����)
       Amu2=B%Amu(i-1,j) + sigma_W_SST*B%Amu_t(i-1,j)        ! ��ɢϵ�� (w����)
     else
       Amu1=(B%Amu(i-1,j)+B%Amu(i,j) + sigma_K_SST*(B%Amu_t(i-1,j)+B%Amu_t(i,j)) )*0.5d0        ! ��ɢϵ�� (k����), �����ϵ�ֵ=����ֵ��ƽ��
       Amu2=(B%Amu(i-1,j)+B%Amu(i,j) + sigma_W_SST*(B%Amu_t(i-1,j)+B%Amu_t(i,j)) )*0.5d0        ! ��ɢϵ�� (w����)
     endif
     
! ������������k,w��������x,y�ĵ��� ������Jocabian�任��
!----Jocabianϵ�� ����������Լ�������ĵ���, ���ڼ����������ĵ�����
      Dix=B%x1(i,j)-B%x1(i-1,j)
      Diy=B%y1(i,j)-B%y1(i-1,j)
      Djx=(B%x1(i-1,j+1)+B%x1(i,j+1)-B%x1(i-1,j-1)-B%x1(i,j-1))*0.25d0
      Djy=(B%y1(i-1,j+1)+B%y1(i,j+1)-B%y1(i-1,j-1)-B%y1(i,j-1))*0.25d0
      Ds=1.d0/(Dix*Djy-Djx*Diy)
! �������Լ�������ĵ���    
      Dik=Kt(i,j)-Kt(i-1,j)
      Diw=Wt(i,j)-Wt(i-1,j)
      Djk=(Kt(i-1,j+1)+Kt(i,j+1)-Kt(i-1,j-1)-Kt(i,j-1))*0.25d0
      Djw=(Wt(i-1,j+1)+Wt(i,j+1)-Wt(i-1,j-1)-Wt(i,j-1))*0.25d0
! ��������x,y����ĵ���
      Kx=(Dik*Djy-Djk*Diy)*Ds
      Wx=(Diw*Djy-Djw*Diy)*Ds
      Ky=(-Dik*Djx+Djk*Dix)*Ds
      Wy=(-Diw*Djx+Djw*Dix)*Ds
!  ճ��Ӧ��������ͨ��
      fluxk(i,j)=fluxk(i,j)+Amu1*(Kx*ni1+Ky*ni2)*si/Re
      fluxw(i,j)=fluxw(i,j)+Amu2*(Wx*ni1+Wy*ni2)*si/Re
	enddo
    enddo
!$OMP END DO

!$OMP  DO
   do j=1,ny-1
   do i=1,nx-1
      B%Res(5,i,j)= Fluxk(i+1,j)-Fluxk(i,j)
      B%Res(6,i,j)= Fluxw(i+1,j)-Fluxw(i,j)
   enddo
   enddo
!$OMP END DO


!----- j- direction ----------------------------------------------------------------------------------

!$OMP  DO
   do j=1,ny
   do i=1,nx-1
   ! �߳���������   
 	 sj=B%sj(i,j)
 	 nj1=B%nj1(i,j); nj2=B%nj2(i,j)

!   ���������1��ӭ���ʽ ��L-F���ѣ�
     
	 un1=uu(i,j-1)*nj1+v(i,j-1)*nj2
	 un2=uu(i,j)*nj1+v(i,j)*nj2
	 fluxk(i,j)=-0.5d0*( (un1+abs(un1))*Kt(i,j-1) + (un2-abs(un2))*Kt(i,j) )*sj
	 fluxw(i,j)=-0.5d0*( (un1+abs(un1))*Wt(i,j-1) + (un2-abs(un2))*Wt(i,j) )*sj

!  ճ����
 !---------Viscous term -----------------------------------------------------------------------------
     sigma_K_SST=f1(i,j)*sigma_k1_SST+(1.d0-f1(i,j))*sigma_k2_SST
     sigma_W_SST=f1(i,j)*sigma_w1_SST+(1.d0-f1(i,j))*sigma_w2_SST
   
    if( j .eq. 1) then
    Amu1=B%Amu(i,j) +sigma_K_SST*B%Amu_t(i,j)
    Amu2=B%Amu(i,j) +sigma_W_SST*B%Amu_t(i,j)
    else if (j .eq. ny) then
    Amu1=B%Amu(i,j-1) +sigma_K_SST*B%Amu_t(i,j-1)
    Amu2=B%Amu(i,j-1) +sigma_W_SST*B%Amu_t(i,j-1)
    else
    Amu1=(B%Amu(i,j)+B%Amu(i,j-1) +sigma_K_SST*(B%Amu_t(i,j)+B%Amu_t(i,j-1)))*0.5d0
    Amu2=(B%Amu(i,j)+B%Amu(i,j-1) +sigma_W_SST*(B%Amu_t(i,j)+B%Amu_t(i,j-1)))*0.5d0
   endif
  
! ������������k,w��������x,y�ĵ��� ������Jocabian�任��
    Dix=(B%x1(i+1,j-1)+B%x1(i+1,j)-B%x1(i-1,j-1)-B%x1(i-1,j))*0.25d0
    Diy=(B%y1(i+1,j-1)+B%y1(i+1,j)-B%y1(i-1,j-1)-B%y1(i-1,j))*0.25d0
    Djx=B%x1(i,j)-B%x1(i,j-1)
    Djy=B%y1(i,j)-B%y1(i,j-1)
    Ds=1.d0/(Dix*Djy-Djx*Diy)

    Dik=(Kt(i+1,j-1)+Kt(i+1,j)-Kt(i-1,j-1)-Kt(i-1,j))*0.25d0
    Diw=(Wt(i+1,j-1)+Wt(i+1,j)-Wt(i-1,j-1)-Wt(i-1,j))*0.25d0
    Djk=Kt(i,j)-Kt(i,j-1)
    Djw=Wt(i,j)-Wt(i,j-1)
!     
    Kx=(Dik*Djy-Djk*Diy)*Ds
    Wx=(Diw*Djy-Djw*Diy)*Ds
    Ky=(-Dik*Djx+Djk*Dix)*Ds
    Wy=(-Diw*Djx+Djw*Dix)*Ds

     fluxk(i,j)=fluxk(i,j)+Amu1*(Kx*nj1+Ky*nj2)*sj/Re
     fluxw(i,j)=fluxw(i,j)+Amu2*(Wx*nj1+Wy*nj2)*sj/Re
 
   enddo
   enddo
!$OMP END DO

!  �в�=������+ճ����+Դ��  

!$OMP  DO
   do j=1,ny-1
   do i=1,nx-1
      B%Res(5,i,j)= B%Res(5,i,j)+Fluxk(i,j+1)-Fluxk(i,j)+Qk(i,j)*B%vol(i,j)
      B%Res(6,i,j)= B%Res(6,i,j)+Fluxw(i,j+1)-Fluxw(i,j)+QW(i,j)*B%vol(i,j)
      
      B%Amu(i,j)=B%Amu(i,j)/Re
      B%Amu_t(i,j)=B%Amu_t(i,j)/Re
   enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
   
  deallocate(Kt,Wt,f1,Qk,Qw,Fluxk,Fluxw)
  
   ! �趨����ճ��ϵ���������ֵ
   call  Amut_boundary(nMesh,mBlock)

 end  subroutine SST_kw
!---------------------------------------------------------------------------
! ���������㵽����ľ��룬SSTģ����Ҫʹ�ø�ֵ
 subroutine comput_dw
   Use Global_Var
   Use Flow_Var
   implicit none
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B,B1
   TYPE (BC_MSG_TYPE),pointer:: Bc1
   integer:: nx,ny,i,j,k,m,mBlock,NB,i1,j1,m1,n1,nx1,ny1
   logical ex
   real*8,allocatable,dimension(:,:):: dw
   real*8:: d1
   
    MP=>Mesh(1)  ! ���ܵ�����
    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny   
       allocate (B%dw(0:nx,0:ny))
       B%dw(:,:)=1.d0
    enddo
    
 ! ������ľ��룬�������ļ����ڣ����ȡ   
  inquire(file="dist_wall.dat",exist=ex)
   if(ex) then
     open(99,file="dist_wall.dat")
     read(99,*) NB
     read(99,*) (nx1,ny1, k=1,NB)
     do k=1,NB
     B => MP%Block(k)
     read(99,*) ((B%dw(i,j),i=1,B%nx-1),j=1,B%ny-1)
     enddo
     close(99)
     return
  endif
  
  ! �������ļ������ڣ��򴴽�
  ! ����ÿ�������ĵ㣬����������㣬�ҳ���С����
   do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny   
       allocate(dw(nx,ny))
       do j=1,ny
       do i=1,nx
         dw(i,j)=1.d10
          do m1=1,MP%Num_Block
           B1=>MP%Block(m1)
            do n1=1,B1%subface
              Bc1=>B1%bc_msg(n1)
              if(Bc1%neighb .ne. BC_Wall) cycle
              do j1=Bc1%jst,Bc1%jend
              do i1=Bc1%ist,Bc1%iend
                d1=sqrt((B%x(i,j)-B1%x(i1,j1))**2+(B%y(i,j)-B1%y(i1,j1))**2)
               if(d1 .lt. dw(i,j)) dw(i,j)=d1                      ! ��Ѱ���������С�ľ���
               enddo
               enddo
             enddo
            enddo
       enddo
       enddo
       
       do j=1,ny-1
       do i=1,nx-1
         B%dw(i,j)=(dw(i,j)+dw(i+1,j)+dw(i,j+1)+dw(i+1,j+1))*0.25d0       ! ���ĵ� ��������ľ��룩=���ܵ��ƽ��ֵ
       enddo
       enddo
      deallocate(dw)
   enddo   
  
  !  �洢�ļ�
   open(99,file="dist_wall.dat")
   write(99,*) MP%Num_Block
   write(99,*) (MP%Block(k)%nx,MP%Block(k)%ny,k=1,MP%Num_Block)
   do m=1,MP%Num_Block
   B=>MP%Block(m)
   write(99,*) ((B%dw(i,j),i=1,B%nx-1),j=1,B%ny-1)
   enddo
   close(99)
 ! tecplot��ʽ�ļ�
   open(100,file="dw.dat")
   write(100,*) "variables=x,y,dw"
   do m=1,MP%Num_Block
   B=>MP%Block(m)
   nx=B%nx; ny=B%ny
   write(100,*) "zone i= ", nx-1, " j= ", ny-1     
   do j=1,ny-1
   do i=1,nx-1
   write(100,*) B%x1(i,j),B%y1(i,j),B%dw(i,j)
   enddo
   enddo
   enddo
   close(100)
   
  end  subroutine  comput_dw 
