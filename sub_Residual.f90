 !-----The core subroutines: Comput Residual: includes the inviscous and viscous flux ---------
 !     Copyright by Li Xinliang, LHD, Institute of Mechanics, CAS. lixl@imech.ac.cn
 !     Ver 1.1
 !--------------------------------------------------------------------------------------------
 ! Since V0.96, Original value is used in Reconstruction
 ! Since V0.971, User can choose using original variables or conservative variables
 ! Since V0.99, User can choose using characteristic variables
 ! Ver 1.1 , 2011-3-28: multi-grid is supported, Code by Li Xinliang
 ! Ver 1.13, 2011-5-22: OMUSCL2 is adding
 !--------------------------------------------------------------------------------------------
 !-------------------------------------------------------------------------------------------------------   
! ����в�����ȫ���飩
   Subroutine Comput_Residual_one_mesh(nMesh)
   use Global_Var
   use Flow_Var 
   implicit none
   real*8:: Res,Res_max(6),Res_rms(6),Sfac
   integer:: nMesh,mBlock,NV,nx,ny,i,j,m,Nvar1
   integer,save:: Iflag1=0
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
!---------------------------------------------  
       MP=>Mesh(nMesh)               ! ��������
       NV=MP%NVAR
       MP%Res_max(:)=0.d0
	   MP%Res_rms(:)=0.d0
!-------Dual Time LU-SGS ----------------------
!   ����1��Sfac=1/3 (�����U(n-1)=U(n)), �����Ϊ1/2     
     if(Iflag1 .eq. 0) then
	   Iflag1=1
	   Sfac=1.d0/(3.d0*dt_global)
	 else
	   Sfac=1.d0/(2.d0*dt_global)
	 endif
!----------------------------------------------
   do mBlock=1,MP%Num_Block
     B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny
!    ������ʱ����
!     allocate(d(0:nx,0:ny),uu(0:nx,0:ny),v(0:nx,0:ny),T(0:nx,0:ny),cc(0:nx,0:ny),p(0:nx,0:ny))   ! Bug found, 2012-5-1
     allocate(d(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),uu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), &
              v(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),T(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),  &
              p(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), cc(nx-1,ny-1))
     allocate(Fluxi(4,nx,ny),Fluxj(4,nx,ny))


!------------------------------------------------------------------------------------     
	 call comput_duvtpckw(nMesh,mBlock)    ! ��������� d,u,v,T,p,cc,Kt,Wt
!---------------------------------------------------------------------------------------
!  ���N-S���̵������ģ��  
!  ����һ�������Ĳв�Ҷ�� ; ��nMesh ������ĵ�mBlock��  ������ģ��Ҳ�ڸÿ��м��㣩   
     call Residual (nMesh,mBlock)       
 
!  ˫ʱ�䲽 LU-SGS�����������     
   if( Time_Method .eq. Time_Dual_LU_SGS ) then 

     do j=1,ny-1
     do i=1,nx-1
     do m=1,NVAR
       B%Res(m,i,j)=B%Res(m,i,j)-(3.d0*B%U(m,i,j)-4.d0*B%Un(m,i,j)+B%Un1(m,i,j))*B%vol(i,j)*Sfac
     enddo
	 enddo
	 enddo

	endif

 
    	   
!--------------------------------------------------------------------------------------
!   ͳ�����в�;������в� (OpenMP��ʹ���˹�Լ����)
  
       Res_max(:)=0.d0   ! ���в�
 	   Res_rms(:)=0.d0   ! �������в� 

! $OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(MP,B,NV) REDUCTION(MAX: Res_max) REDUCTION(+: Res_rms)   
    do j=1,ny-1
    do i=1,nx-1
! -------------------------------------------------------------------------------------------
!    ʱ���ƽ�
       do m=1,NV
        Res=B%Res(m,i,j)
!--------------------------------------------------------------------------------------------------
        if(abs(Res) .gt. Res_max(m))  Res_max(m)=abs(Res)          ! ���в�
        Res_rms(m)=Res_rms(m)+Res*Res                              ! �������в�      
!--------------------------------------------------------------------------------------------------       
	   enddo
    enddo
    enddo
! $OMP END PARALLEL DO
     do m=1,NV 
	  MP%Res_max(m)=max(MP%Res_max(m), Res_max(m))          ! ���в�
      MP%Res_rms(m)=MP%Res_rms(m)+Res_rms(m)                !    ȫ���������ܾ������в�
     enddo
   
!---------------------------------------------------------------------------------------- 
     call comput_Lvc(nMesh,mBlock)
	 call comput_dt(nMesh,mBlock)        ! ����(����) ʱ�䲽��
 
 !   ɾ������ʱ����    
	 deallocate(d,uu,v,T,cc,p,Fluxi,Fluxj)
 
   enddo    
    
	 do m=1,NV 
	   MP%Res_rms(m)=sqrt(MP%Res_rms(m)/(1.d0*MP%Num_Cell))   !    ȫ���������ܾ������в�
     enddo
  end Subroutine Comput_Residual_one_mesh



!---------------------------------------------------------
!  �����غ��������������� (d,u,v,T,p,c) 
!----------------------------------------------------------
    subroutine comput_duvtpckw(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     implicit none
 
     Type (Mesh_TYPE),pointer:: MP
     Type (Block_TYPE),pointer:: B
	 integer nMesh,mBlock,nx,ny,i,j
     real*8 p00
	 p00=1.d0/(gamma*Ma*Ma)
     MP=> mesh(nMesh)
     B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny
 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
	  do j=1-LAP,ny+LAP-1
      do i=1-LAP,nx+LAP-1
        d(i,j)=B%U(1,i,j)
        uu(i,j)=B%U(2,i,j)/d(i,j)
        v(i,j)=B%U(3,i,j)/d(i,j)
        T(i,j)=(B%U(4,i,j)-0.5d0*d(i,j)*(uu(i,j)*uu(i,j)+v(i,j)*v(i,j)))/(Cv*d(i,j))
        p(i,j)=p00*d(i,j)*T(i,j)
      enddo
      enddo
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do j=1,ny-1
      do i=1,nx-1
        if(T(i,j) .lt. 1.d-5) then
          print*, "T <1.d-5 ! ", i,j,mBlock, T(i,j)
          stop
        endif
        cc(i,j)=sqrt(T(i,j))/Ma
      enddo
      enddo
!$OMP END PARALLEL DO      
    
    end subroutine comput_duvtpckw
!---------------------------------------------------------------------------------

 
 
 ! �����nMesh������ĵ�mBlock��Ĳв� ���Ҷ�� 
   subroutine Residual(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   integer:: nMesh,mBlock
 
     MP=> Mesh(nMesh)
     B => MP%Block(mBlock)

! ����ճ��ϵ��(������������ģ�ͼ�������ճ��ϵ��)  
   if(If_viscous .eq. 1) then   
    call get_viscous(nMesh,mBlock)         ! �������ճ��ϵ��
    B%Amu_t(:,:)=0.d0                       ! ����ճ��ϵ������
!   ����ģ�ͣ���������ճ��ϵ��Amu_t  ��BL,SA,SSTģ�ͣ�
    if(MP%Iflag_turbulence_model .eq. Turbulence_BL) then
     call  turbulence_model_BL(nMesh,mBlock) 
    else if(MP%Iflag_turbulence_model .eq. Turbulence_SA) then
	 call turbulence_model_SA(nMesh,mBlock)
	else if (MP%Iflag_turbulence_model .eq. Turbulence_SST) then
     call SST_kw(nMesh,mBlock)
    endif
   endif
 
 ! �������������в� (��ÿ�ʹ�ò�ַ�����ֻ�����ٽ���߽��4������
   call Residual_FVM(nMesh,mBlock)

! �������޲�ַ�����в� (�ٽ���߽��4���������)
   if(B%FVM_FDM .eq. Method_FDM) then
       call Residual_FDM(nMesh,mBlock)
   endif

   end subroutine Residual

 

!------------------------------------------------------------------------------
! ʹ���������������в� (�в�=�Ҷ���=������); �����N-S���̵ĺ���ģ��; 

! Flux (= inviscous flux + viscous flux )

   subroutine Residual_FVM(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   real*8,dimension(4):: UL,UR,UL1,UR1,UL2,UR2,UL3,UR3,QL,QR,Flux0
   real*8:: U0(4,4)
   integer:: nMesh,mBlock,Scheme,IFlux,Reconstruction
   integer:: nx,ny,i,j,m
   real*8:: minmod,dx,dy,si,ni1,ni2,sj,nj1,nj2
   real*8:: Diu,Div,Dju,Djv,DiT,DjT,ux,vx,Tx,uy,vy,Ty,t11,t12,t22,E1,E2
   real*8:: Dix,Diy,Djx,Djy,Ds,Amu1,u1,v1,Amk1
   real*8:: dl,uul,vl,ppl,dr,uur,vr,ppr,pr1
   TYPE (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   MP=> Mesh(nMesh)
   B => MP%Block(mBlock)
   nx=B%nx ; ny=B%ny

   Scheme=MP%Iflag_Scheme         ! ��ֵ��ʽ ����ͬ�����ϲ��ò�ͬ��ʽ��
   IFlux=MP%Iflag_Flux            ! ͨ������
   Reconstruction=MP%IFlag_Reconstruction  ! �ع���ʽ
  
!$OMP PARALLEL DEFAULT (PRIVATE) SHARED(nx,ny,MP,B,d,uu,v,p,T,Fluxi,Fluxj,Scheme,IFlux,Reconstruction,If_viscous,gamma,Cp,Pr)

!$OMP DO
   do j=1,ny
   do i=1,nx
   do m=1,4
   Fluxi(m,i,j)=0.d0  ! ����
   Fluxj(m,i,j)=0.d0
   enddo
   enddo
   enddo
!$OMP ENDDO


!----- i- direction ----------------------------------------------------------------------------------

!$OMP DO
   do j=1,ny-1
   do i=1,nx

!  ����ò�ַ����㣬����������ʹ��������������� 
    if(B%FVM_FDM .eq. Method_FDM  .and. (i .ge. LFDM+2 .and. i .le. nx-LFDM-1) &
	   .and.  (j .ge. LFDM+2 .and. j .le. ny-LFDM-1) ) cycle    ! �������� �����޲�ַ��ļ�������

	 si=B%si(i,j)   ! ������߽糤��
	 ni1=B%ni1(i,j); ni2=B%ni2(i,j)   ! ������
!-----��ճ�� -----------Inviscous term-------------------------------------------------------------

!-----ѡ���ع�����ֵ������ֵ�� ��ʽ  ��ʹ��ԭʼ�������غ���������������ع���

     if(Reconstruction .eq. Reconst_Original) then
       U0(:,1)=d(i-2:i+1,j) ; U0(:,2)=uu(i-2:i+1,j) ; U0(:,3)=v(i-2:i+1,j); U0(:,4)=p(i-2:i+1,j)
       call Reconstuction_original(U0,UL,UR,gamma,Scheme)
     else if (Reconstruction .eq. Reconst_Conservative) then
      do m=1,4
       U0(:,m)=B%U(m,i-2:i+1,j)
      enddo
      call Reconstuction_conservative(U0,UL,UR,gamma,Scheme)
     else
      do m=1,4
       U0(:,m)=B%U(m,i-2:i+1,j)
      enddo
      call Reconstuction_Characteristic(U0,UL,UR,gamma,Scheme)
     endif
   

!-------������ת����ֱ�ڽ��� ������-���� ����ϵ��
      QL(1)=UL(1); QL(2)=UL(2)*ni1+UL(3)*ni2 ; QL(3)= -UL(2)*ni2+UL(3)*ni1 ; QL(4)=UL(4) ! �ܶȡ�ѹ���������ٶȡ������ٶ� ����ֵ��
      QR(1)=UR(1); QR(2)=UR(2)*ni1+UR(3)*ni2 ; QR(3)= -UR(2)*ni2+UR(3)*ni1 ; QR(4)=UR(4) ! �ܶȡ�ѹ���������ٶȡ������ٶ� ����ֵ��

!-------ѡ��ͨ�����ѷ���----(Steger_Warming,Van Leer, Roe, AUSM, HLL, HLLC----
      if(IFlux .eq. Flux_Steger_Warming ) then
         call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)   
      else if (IFlux .eq. Flux_HLL ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,0)
      else if (IFlux .eq. Flux_HLLC ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,1)
      else if  (IFlux .eq. Flux_VanLeer ) then
         call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_AUSM ) then
         call Flux_Ausm_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_Roe ) then
         call Flux_Roe_1D(QL,QR,Flux0,gamma)   
      endif


!------------------------------------------------	  
!   ���ͨ�� ���任��x-y����ϵ��  
        Fluxi(1,i,j)=-Flux0(1)*si                            ! ����ͨ��
        Fluxi(2,i,j)=-(Flux0(2)*ni1-Flux0(3)*ni2)*si         ! x-������ͨ��
        Fluxi(3,i,j)=-(Flux0(2)*ni2+Flux0(3)*ni1)*si         ! y-������ͨ��
        Fluxi(4,i,j)=-Flux0(4)*si                            ! ����ͨ��
      
	  if (If_viscous .eq. 0) cycle   ! �����Euler���̣���ճ����������ճ�������

!--i-������ճͨ���������������i-����ճ��ͨ��
!--------------------------------------------------------------------------------------------------------- 
!----------- Viscous term---------------------------------------------------------------------------------  
! ��ɢϵ����ճ��ϵ�����ȴ���ϵ����= ��������ֵ��ƽ�� ���߽紦���õ���ֵ)    
    if( i .eq. 1) then 
     Amu1=B%Amu(i,j) + B%Amu_t(i,j)                  ! ճ��ϵ�� (����+����), �����ϵ�ֵ=����ֵ��ƽ��
     Amk1=Cp*(B%Amu(i,j)/Pr + B%Amu_t(i,j)/PrT )     ! �ȴ���ϵ��
    else if (i .eq. nx ) then
      Amu1=B%Amu(i-1,j) + B%Amu_t(i-1,j)                  ! ճ��ϵ�� (����+����), �����ϵ�ֵ=����ֵ��ƽ��
      Amk1=Cp*(B%Amu(i-1,j)/Pr + B%Amu_t(i-1,j)/PrT )     ! �ȴ���ϵ��
    else
      Amu1=(B%Amu(i-1,j)+B%Amu(i,j) + B%Amu_t(i-1,j)+B%Amu_t(i,j) )*0.5d0                 ! ճ��ϵ�� (����+����), �����ϵ�ֵ=����ֵ��ƽ��
      Amk1=Cp*((B%Amu(i-1,j)+B%Amu(i,j))/Pr + (B%Amu_t(i-1,j)+B%Amu_t(i,j))/PrT )*0.5d0   ! �ȴ���ϵ��
    endif

!----Jocabianϵ�� ����������Լ�������ĵ���, ���ڼ����������ĵ�����
    Dix=B%x1(i,j)-B%x1(i-1,j)
    Diy=B%y1(i,j)-B%y1(i-1,j)
    Djx=(B%x1(i-1,j+1)+B%x1(i,j+1)-B%x1(i-1,j-1)-B%x1(i,j-1))*0.25d0
    Djy=(B%y1(i-1,j+1)+B%y1(i,j+1)-B%y1(i-1,j-1)-B%y1(i,j-1))*0.25d0
    Ds=1.d0/(Dix*Djy-Djx*Diy)
! �������Լ�������ĵ���    
    Diu=uu(i,j)-uu(i-1,j)
    Div=v(i,j)-v(i-1,j)
    DiT=T(i,j)-T(i-1,j)
    Dju=(uu(i-1,j+1)+uu(i,j+1)-uu(i-1,j-1)-uu(i,j-1))*0.25d0
    Djv=(v(i-1,j+1)+v(i,j+1)-v(i-1,j-1)-v(i,j-1))*0.25d0
    DjT=(T(i-1,j+1)+T(i,j+1)-T(i-1,j-1)-T(i,j-1))*0.25d0
! ��������x,y����ĵ���
    ux=(Diu*Djy-Dju*Diy)*Ds
    vx=(Div*Djy-Djv*Diy)*Ds
    Tx=(DiT*Djy-DjT*Diy)*Ds
    uy=(-Diu*Djx+Dju*Dix)*Ds
    vy=(-Div*Djx+Djv*Dix)*Ds
    Ty=(-DiT*Djx+DjT*Dix)*Ds
!  ճ��Ӧ��������ͨ��
    u1=(uu(i,j)+uu(i-1,j))*0.5d0
    v1=(v(i,j)+v(i-1,j))*0.5d0
    t11=((4.d0/3.d0)*ux-(2.d0/3.d0)*vy)*Amu1
    t22=((4.d0/3.d0)*vy-(2.d0/3.d0)*ux)*Amu1
    t12=(uy+vx)*Amu1
    E1=u1*t11+v1*t12+Amk1*Tx
    E2=u1*t12+v1*t22+Amk1*Ty
! ���ճ��ͨ��
    Fluxi(2,i,j)=Fluxi(2,i,j)   +(t11*ni1+t12*ni2)*si   
    Fluxi(3,i,j)=Fluxi(3,i,j)   +(t12*ni1+t22*ni2)*si
    Fluxi(4,i,j)=Fluxi(4,i,j)   +(E1*ni1+ E2*ni2)*si
  
   enddo
   enddo
!$OMP ENDDO
  
!==================================================================================================================
!                                         j-�������ճ��ճ��ͨ��    
!-----------------------------------------j- direction -------------------------------------------------------------	 
!$OMP DO
   do j=1,ny
   do i=1,nx-1
    
	 if(B%FVM_FDM .eq. Method_FDM  .and. (i .ge. LFDM+2 .and. i .le. nx-LFDM-1) &
	   .and.  (j .ge. LFDM+2 .and. j .le. ny-LFDM-1) ) cycle    ! �������� �����޲�ַ��ļ�������
 
   ! �߳���������   
	sj=B%sj(i,j)
	nj1=B%nj1(i,j); nj2=B%nj2(i,j)

!-----Inviscous term --------------------------------------------
 
    if(Reconstruction .eq. Reconst_Original) then
       U0(:,1)=d(i,j-2:j+1) ; U0(:,2)=uu(i,j-2:j+1) ; U0(:,3)=v(i,j-2:j+1); U0(:,4)=p(i,j-2:j+1)
       call Reconstuction_original(U0,UL,UR,gamma,Iflag_Scheme)
    else if (Reconstruction .eq. Reconst_Conservative) then
      do m=1,4
       U0(:,m)=B%U(m,i,j-2:j+1)
      enddo
       call Reconstuction_conservative(U0,UL,UR,gamma,Scheme)
    else
      do m=1,4
       U0(:,m)=B%U(m,i,j-2:j+1)
      enddo
      call Reconstuction_characteristic(U0,UL,UR,gamma,Scheme)
    endif	   
    

       QL(1)=UL(1); QL(2)=UL(2)*nj1+UL(3)*nj2 ; QL(3)= -UL(2)*nj2+UL(3)*nj1 ; QL(4)=UL(4) !density, normal/tangitial velocity and pressure
       QR(1)=UR(1); QR(2)=UR(2)*nj1+UR(3)*nj2 ; QR(3)= -UR(2)*nj2+UR(3)*nj1 ; QR(4)=UR(4)
 
!-------ѡ��ͨ�����ѷ���----(Steger_Warming,Van Leer, Roe, AUSM, HLL, HLLC)----
      if(IFlux .eq. Flux_Steger_Warming ) then
         call Flux_steger_warming_1Da(QL,QR,Flux0,gamma)   
      else if (IFlux .eq. Flux_HLL ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,0)
      else if (IFlux .eq. Flux_HLLC ) then
         call Flux_HLL_HLLC_1D(QL,QR,Flux0,gamma,1)
      else if  (IFlux .eq. Flux_VanLeer ) then
         call Flux_Van_Leer_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_AUSM ) then
         call Flux_Ausm_1Da(QL,QR,Flux0,gamma)   
      else if  (IFlux .eq. Flux_Roe ) then
         call Flux_Roe_1D(QL,QR,Flux0,gamma)   
      endif
     
!   ��ճͨ��   	   
    Fluxj(1,i,j)=-Flux0(1)*sj 
    Fluxj(2,i,j)=-(Flux0(2)*nj1-Flux0(3)*nj2)*sj
    Fluxj(3,i,j)=-(Flux0(2)*nj2+Flux0(3)*nj1)*sj
    Fluxj(4,i,j)=-Flux0(4)*sj

    if (If_viscous .eq. 0) cycle  ! �����Euler���̣���ճ����������ճ�������
 !---------Viscous term -----------------------------------------------------------------------------
 ! ��ɢϵ����ճ��ϵ�����ȴ���ϵ����= ��������ֵ��ƽ�� ���߽紦���õ���ֵ)    
    if( j .eq. 1) then
     Amu1=B%Amu(i,j)+B%Amu_t(i,j)
     Amk1=Cp*(B%Amu(i,j)/Pr +B%Amu_t(i,j)/PrT )   ! �ȴ���ϵ��
    else if (j .eq. ny) then
      Amu1=B%Amu(i,j-1)+B%Amu_t(i,j-1)
      Amk1=Cp*(B%Amu(i,j-1)/Pr +B%Amu_t(i,j-1)/PrT )   ! �ȴ���ϵ��
    else
     Amu1=(B%Amu(i,j)+B%Amu(i,j-1) +B%Amu_t(i,j)+B%Amu_t(i,j-1))*0.5d0
     Amk1=Cp*((B%Amu(i,j-1)+B%Amu(i,j))/Pr + (B%Amu_t(i,j-1)+B%Amu_t(i,j))/PrT )*0.5d0   ! �ȴ���ϵ��
    endif
    
    Dix=(B%x1(i+1,j-1)+B%x1(i+1,j)-B%x1(i-1,j-1)-B%x1(i-1,j))*0.25d0
    Diy=(B%y1(i+1,j-1)+B%y1(i+1,j)-B%y1(i-1,j-1)-B%y1(i-1,j))*0.25d0
    Djx=B%x1(i,j)-B%x1(i,j-1)
    Djy=B%y1(i,j)-B%y1(i,j-1)
   
    Ds=1.d0/(Dix*Djy-Djx*Diy)

    Diu=(uu(i+1,j-1)+uu(i+1,j)-uu(i-1,j-1)-uu(i-1,j))*0.25d0
    Div=(v(i+1,j-1)+v(i+1,j)-v(i-1,j-1)-v(i-1,j))*0.25d0
    DiT=(T(i+1,j-1)+T(i+1,j)-T(i-1,j-1)-T(i-1,j))*0.25d0
    Dju=uu(i,j)-uu(i,j-1)
    Djv=v(i,j)-v(i,j-1)
    DjT=T(i,j)-T(i,j-1)
!     
    ux=(Diu*Djy-Dju*Diy)*Ds
    vx=(Div*Djy-Djv*Diy)*Ds
    Tx=(DiT*Djy-DjT*Diy)*Ds
    uy=(-Diu*Djx+Dju*Dix)*Ds
    vy=(-Div*Djx+Djv*Dix)*Ds
    Ty=(-DiT*Djx+DjT*Dix)*Ds
    t11=((4.d0/3.d0)*ux-(2.d0/3.d0)*vy)*Amu1
    t22=((4.d0/3.d0)*vy-(2.d0/3.d0)*ux)*Amu1
    t12=(uy+vx)*Amu1
    u1=(uu(i,j)+uu(i,j-1))*0.5d0
    v1=(v(i,j)+v(i,j-1))*0.5d0

    E1=u1*t11+v1*t12+Amk1*Tx
    E2=u1*t12+v1*t22+Amk1*Ty

    Fluxj(2,i,j)=Fluxj(2,i,j) +(t11*nj1+t12*nj2)*sj
    Fluxj(3,i,j)=Fluxj(3,i,j) +(t12*nj1+t22*nj2)*sj
    Fluxj(4,i,j)=Fluxj(4,i,j) +(E1*nj1+ E2*nj2)*sj

 
   enddo
   enddo
!$OMP ENDDO
  
 !-------����в� ����������----------------------------------------------------  
 !                                   ����ǿ�Ⱥ��� ����������ʱ�������ǿ�Ⱥ�����
!$OMP DO 
    do j=1,ny-1
    do i=1,nx-1
	 if(B%FVM_FDM .eq. Method_FDM  .and. (i .ge. LFDM+1 .and. i .le. nx-LFDM-1) &
	   .and.  (j .ge. LFDM+1 .and. j .le. ny-LFDM-1) ) cycle    ! �������� �����޲�ַ��ļ�������
    do m=1,4
      B%Res(m,i,j)= Fluxi(m,i+1,j)-Fluxi(m,i,j)+Fluxj(m,i,j+1)-Fluxj(m,i,j)
	enddo
	enddo
	enddo
!$OMP ENDDO

!$OMP END PARALLEL 
  end  

!---------------------------------------------------------------------------------------


      function minmod(a,b)
      implicit none
      real*8 a,b,minmod
      if(a*b .le. 0.d0) then
       minmod=0.d0
      else if (abs(a) .le. abs(b)) then
       minmod=a
      else
       minmod=b
      endif
      end
!--------------------------------------------------------------------------------


!-------------------------------------------------------------------
! ����ԭʼ��������ƽ�� �����ڼ���ǵ��ֵ��
  subroutine U_average(U1,U2,U0,gamma)
  implicit none
  real*8,dimension(4):: U1,U2,U0
  real*8:: gamma,d1,uu1,v1,p1,d2,uu2,v2,p2,d0,uu0,v0,p0
    d1=U1(1); uu1=U1(2)/d1; v1=U1(3)/d1; p1=(U1(4)-(uu1*U1(2)+v1*U1(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure 
    d2=U2(1); uu2=U2(2)/d2; v2=U2(3)/d2; p2=(U2(4)-(uu2*U2(2)+v2*U2(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure 
    d0=(d1+d2)*0.5d0; uu0=(uu1+uu2)*0.5d0; v0=(v1+v2)*0.5d0; p0=(p1+p2)*0.5d0
    U0(1)=d0; U0(2)=d0*uu0; U0(3)=d0*v0; U0(4)=p0/(gamma-1.d0)+(d0*uu0*uu0+d0*v0*v0)*0.5d0
  end

!--------------------------------------------------------------------------
! ����ԭʼ�����ع� 
   subroutine Reconstuction_original(U0,UL,UR,gamma,Iflag_Scheme)
     implicit none
     real*8:: U0(4,4),UL(4),UR(4),gamma
     integer:: Iflag_Scheme,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1    ;   m=1,4 for d,u,v,p
     do m=1,4
       call scheme_fP(UL(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
       call scheme_fm(UR(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
     enddo
   end

! �����غ�����ع�
   subroutine Reconstuction_conservative(U0,UL,UR,gamma,Iflag_Scheme)
     implicit none
     real*8:: U0(4,4),UL(4),UR(4),QL(4),QR(4),gamma
     integer:: Iflag_Scheme,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1   ; m for the conservative variables U0(1,m)=d, U0(2,m)=d*u, ....
     do m=1,4
        call scheme_fP(QL(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
        call scheme_fm(QR(m),U0(1,m),U0(2,m),U0(3,m),U0(4,m),Iflag_Scheme) 
      enddo
       UL(1)=QL(1); UL(2)=QL(2)/UL(1); UL(3)=QL(3)/UL(1); UL(4)=(QL(4)-(UL(2)*QL(2)+UL(3)*QL(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
       UR(1)=QR(1); UR(2)=QR(2)/UR(1); UR(3)=QR(3)/UR(1); UR(4)=(QR(4)-(UR(2)*QR(2)+UR(3)*QR(3))*0.5d0)*(gamma-1.d0)  ! find a bug, removed
   end

!  �������������ع�
   subroutine Reconstuction_Characteristic(U0,UL,UR,gamma,Iflag_Scheme)
     implicit none
     real*8:: U0(4,4),UL(4),UR(4),gamma
     real*8:: Uh(4),S(4,4),S1(4,4),V0(4,4),VL(4),VR(4),QL(4),QR(4)
     real*8:: V2,d1,u1,v1,p1,c1,tmp0,tmp1,tmp3,tmp5
     integer:: Iflag_Scheme,i,j,k,m
!    U0(k,m) : k=1,4 for  i-2,i-1,i,i+1   ; m for the conservative variables U0(1,m)=d, U0(2,m)=d*u, ....
     Uh(:)=0.5d0*(U0(2,:)+U0(3,:))           ! conservative variables in the point I-1/2  (or i)
     d1=Uh(1); u1=Uh(2)/d1; v1=Uh(3)/d1; p1=(Uh(4)-(Uh(2)*u1+Uh(3)*v1)*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
     c1=sqrt(gamma*p1/d1)

	  V2=(u1*u1+v1*v1)*0.5d0
	  tmp1=(gamma-1.d0)/c1
      tmp3=(gamma-1.d0)/(c1*c1)
	  tmp5=1.d0/(2.d0*c1)
	  tmp0=1.d0/tmp3

!  A=S(-1)*LAMDA*S    see �������������ѧ�� 158-159ҳ   (with alfa=1, beta=0)

          S(1,1)=V2-tmp0;       S(1,2)=-u1 ;         S(1,3)=-v1 ;      S(1,4)=1.d0
          S(2,1)=-v1 ;          S(2,2)=0.d0 ;        S(2,3)=1.d0 ;     S(2,4)=0.d0 
          S(3,1)=-u1-V2*tmp1;   S(3,2)=1.d0+tmp1*u1; S(3,3)=tmp1*v1;   S(3,4)=-tmp1
          S(4,1)=-u1+V2*tmp1;   S(4,2)=1.d0-tmp1*u1; S(4,3)=-tmp1*v1;  S(4,4)=tmp1 
 
          S1(1,1)=-tmp3;    S1(1,2)=0.d0;   S1(1,3)=-tmp5 ;         S1(1,4)=tmp5
          S1(2,1)=-tmp3*u1; S1(2,2)=0.d0;   S1(2,3)=0.5d0-u1*tmp5 ; S1(2,4)=0.5d0+u1*tmp5
          S1(3,1)=-tmp3*v1; S1(3,2)=1.d0;   S1(3,3)=-v1*tmp5;       S1(3,4)=v1*tmp5
          S1(4,1)=-tmp3*V2; S1(4,2)=v1;     S1(4,3)=(c1*u1-V2-tmp0)*tmp5; S1(4,4)=(c1*u1+V2+tmp0)*tmp5
     
! V=SU      V(k)=S*U(k)
    do k=1,4
     
     do m=1,4
     V0(k,m)=0.d0
      do j=1,4
      V0(k,m)=V0(k,m)+S(m,j)*U0(k,j)
      enddo
     enddo
    enddo
     
     do m=1,4
       call scheme_fP(VL(m),V0(1,m),V0(2,m),V0(3,m),V0(4,m),Iflag_Scheme) 
       call scheme_fm(VR(m),V0(1,m),V0(2,m),V0(3,m),V0(4,m),Iflag_Scheme) 
     enddo

    do m=1,4
    QL(m)=0.d0; QR(m)=0.d0
    do j=1,4
    QL(m)=QL(m)+S1(m,j)*VL(j)
    QR(m)=QR(m)+S1(m,j)*VR(j)
    enddo
    enddo
     UL(1)=QL(1); UL(2)=QL(2)/UL(1); UL(3)=QL(3)/UL(1); UL(4)=(QL(4)-(UL(2)*QL(2)+UL(3)*QL(3))*0.5d0)*(gamma-1.d0)  ! density, velocity, pressure and sound speed
     UR(1)=QR(1); UR(2)=QR(2)/UR(1); UR(3)=QR(3)/UR(1); UR(4)=(QR(4)-(UR(2)*QR(2)+UR(3)*QR(3))*0.5d0)*(gamma-1.d0)  ! find a bug, removed

   end


!-----------------------------------------------------------
! ��ֵ��ʽ������UL=U(j+1/2,L); u1=u(j-1),u2=u(j),u3=u(j+1),u4=u(j+2)
    subroutine scheme_fP(uL,u1,u2,u3,u4,Iflag_Scheme)
     use   const_var
     implicit none
     real*8:: uL,u1,u2,u3,u4,minmod
     real*8:: IS1,IS2,a1,a2,w1,w2,up,um,s
     real*8,parameter:: k=1.d0,b=2.d0
	 real*8:: r1,r2,f1,f
     real*8,parameter::k3=1.d0/3.d0,ep=1.d-6
     integer:: Iflag_Scheme
     if(Iflag_Scheme .eq. Scheme_UD1 ) then
       UL=u2                            ! 1��ӭ���ʽ
	 else if(Iflag_Scheme .eq. Scheme_UD3 ) then
       UL=(-u1+5.d0*u2+2.d0*u3  )/6.d0   !UD 3nd order  3��ӭ���ʽ
     else if(Iflag_Scheme .eq. Scheme_NND2 ) then
       UL=u2+0.5d0*minmod(u2-u1,u3-u2)   ! NND 2nd order  2��NND
     else if(Iflag_Scheme .eq. Scheme_WENO3) then  ! WENO3 scheme 
       IS1=(u2-u1)**2 ;    IS2=(u3-u2)**2
       a1=1.d0/(3.d0*(1.d-6+IS1)**2) ;   a2=2.d0/(3.d0*(1.d-6+IS2)**2)
       w1=a1/(a1+a2) ;     w2=a2/(a1+a2) 
       UL=w1*(-u1/2.d0+3.d0*u2/2.d0)+w2*(u2/2.d0+u3/2.d0)   !3��WENO
     else if(Iflag_Scheme .eq. Scheme_MUSCL2 ) then         !2��MUSCL(Minmod��������
       UL=u2+0.25d0*((1.d0-k)*minmod(u2-u1,b*(u3-u2))+(1.d0+k)*minmod(u3-u2,b*(u2-u1))) ! MUSCL
     else if(Iflag_Scheme .eq. Scheme_MUSCL3 ) then          ! 3��MUSCL (Van Albada������)
        up=u3-u2; um=u2-u1                                   ! 1�� ǰ����
        s=(2.d0*up*um+ep)/(up*up+um*um+ep)                   ! Van Albada������ ���⻬����ǰ������ӽ�����ֵ�ӽ�1��
        UL=u2+0.25d0*s*((1.d0-k3*s)*um+(1.d0+k3*s)*up)       ! 3��MUSCL (�⻬���ƽ�3��ӭ��)
     else if(Iflag_Scheme .eq. Scheme_OMUSCL2 ) then         ! 2���Ż���MUSCL ��ʽ (Developed by Leng Yan)
        r1=(u2-u1+ep)/(u3-u2+ep); r2=(u3-u2+ep)/(u4-u3+ep)
		f1=0.8d0-0.175d0/r2+0.375d0*r1
		f=max(0.d0,min(2.d0,f1,2.d0*r1))
        UL=u2+0.5d0*f*(u3-u2)                               ! 2���Ż���MUSCL: OMUSCL2
	 endif
     end

! ��ֵ��ʽ������UR=U(j+1/2,R) ; u1=u(j-1),u2=u(j),u3=u(j+1),u4=u(j+2)
     subroutine scheme_fm(uR,u1,u2,u3,u4,Iflag_Scheme)
     use   const_var
     implicit none
     real*8:: uR,u1,u2,u3,u4,minmod,up,um,s
     real*8:: IS1,IS2,a1,a2,w1,w2
  	 real*8:: r1,r2,f1,f
     real*8,parameter:: k=1.d0/3.0,b=1.d0
     real*8,parameter::k3=1.d0/3.d0,ep=1.d-6
     integer:: Iflag_Scheme
     if(Iflag_Scheme .eq. Scheme_UD1) then
	   UR=u3                              ! 1��ӭ��
	 else if(Iflag_Scheme .eq. Scheme_UD3 ) then
       UR=(2.d0*u2+5.d0*u3-u4  )/6.d0     !UD 3nd order
     else if(Iflag_Scheme .eq.Scheme_NND2 ) then
        UR=u3-0.5d0*minmod(u3-u2,u4-u3)    ! NND 2nd order
     else if(Iflag_Scheme .eq. Scheme_WENO3 ) then   ! WENO3   
        IS1=(u3-u4)**2;   IS2=(u2-u3)**2
        a1=1.d0/(3.d0*(1.d-6+IS1)**2) ;   a2=2.d0/(3.d0*(1.d-6+IS2)**2)
        w1=a1/(a1+a2) ;     w2=a2/(a1+a2) 
        UR=w1*(-u4/2.d0+3.d0*u3/2.d0)+w2*(u3/2.d0+u2/2.d0)       !WENO 3nd order
     else if(Iflag_Scheme .eq. Scheme_MUSCL2 ) then
        uR=u3-0.25d0*((1.d0-k)*minmod(u4-u3,b*(u3-u2))+(1.d0+k)*minmod(u3-u2,b*(u4-u3)))   ! MUSCL
     else if(Iflag_Scheme .eq. Scheme_MUSCL3 ) then      ! 3��MUSCL (Van Albada������)
         up=u4-u3 ; um=u3-u2                             !ǰ����
         s=(2.d0*up*um+ep)/(up*up+um*um+ep)
         UR=u3-0.25d0*s*((1.d0-k3*s)*up+(1.d0+k3*s)*um)
     else if(Iflag_Scheme .eq. Scheme_OMUSCL2 ) then         ! 2���Ż���MUSCL ��ʽ (Developed by Leng Yan)
        r1=(u4-u3+ep)/(u3-u2+ep); r2=(u3-u2+ep)/(u2-u1+ep)
		f1=0.8d0-0.175d0/r2+0.375d0*r1
		f=max(0.d0,min(2.d0,f1,2.d0*r1))
        UR=u3-0.5d0*f*(u3-u2)                               ! 2���Ż���MUSCL: OMUSCL2
    
     endif
     end
!------------------------------------------------------------------------------
 






!----------------------------------------------------------------
! ����ճ��ϵ���ļ���
   subroutine get_viscous(nMesh,mBlock)
   Use Global_Var
   Use Flow_Var
   implicit none
   real*8:: Tsb
   integer:: nMesh,mBlock,i,j,nx,ny
   Type (Block_TYPE),pointer:: B

   B => Mesh(nMesh)%Block(mBlock)
!  Ref_Amu_T0=288.15d0
!  T_inf ���������ο��¶�
     Tsb=110.4d0/T_inf
     nx=B%nx; ny=B%ny
     do j=0,ny
     do i=0,nx
       B%Amu(i,j)=1.d0/Re*(1.d0+Tsb)*sqrt(T(i,j)**3)/(Tsb+T(i,j))   !sutherland equation
 !     Amu(i,j)=1.d0/Re
     enddo
     enddo
  end

