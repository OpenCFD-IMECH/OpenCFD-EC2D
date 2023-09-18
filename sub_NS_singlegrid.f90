!----------------------------------------------------------------------
! Copyright by Li Xinliang,  Code by Li Xinliang
!
!----------------------------------------------------------------------
! �ڸ��������������N-S���� ���ƽ�1��ʱ�䲽��
! ���ڵ�������nMesh=1;  ���ڶ�������nMesh=1,2,3, ... �ֱ��Ӧ��ϸ���񡢴����񡢸������� ...
 
  subroutine NS_Time_advance(nMesh)
   use Global_var
   implicit none
   integer:: nMesh
   if(Time_Method .eq. Time_Euler1) then
      call NS_Time_advance_1Euler(nMesh)                   ! 1��Euler
   else if (Time_Method .eq. Time_RK3) then
     call NS_Time_advance_RK3(nMesh)                       ! 3��RK
   else if (Time_Method .eq. Time_LU_SGS) then             ! LU-SGS����ʽ
     call NS_time_advance_LU_SGS(nMesh)
   else if(Time_Method .eq. Time_Dual_LU_SGS) then         ! ˫ʱ�䲽LU-SGS (��֧�ֶ�������)
     call Dual_time_LU_SGS
   else
      print*, "This time advance method is not supported!!!"
   endif
   !  ǿ�� k,w, vt �Ǹ�
    call force_vt_kw(nMesh)

   end subroutine NS_Time_advance

!---------------------------------------------------------------------------------------------
! ǿ��vt, k,w�Ǹ�   
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




! ����1��Euler������ʱ���ƽ�һ��ʱ�䲽 ����nMesh������ �ĵ�������
 subroutine NS_Time_advance_1Euler(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Block_TYPE),pointer:: B
   Type (Mesh_TYPE),pointer:: MP
   real*8:: du
 
    MP=>Mesh(nMesh)
    call Comput_Residual_one_mesh(nMesh)     ! ���������ϼ���в�
    if(nMesh .ne. 1) call Add_force_function(nMesh)   !  ���ǿ�Ⱥ�������������Ĵ�����ʹ�ã�

    do mBlock=1,MP%Num_Block
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   ʱ���ƽ� 
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
      call Boundary_condition_onemesh(nMesh)         ! �߽����� ���趨Ghost Cell��ֵ��
      call update_buffer_onemesh(nMesh)              ! ͬ������Ľ�����


     Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! ʱ�� ��ʹ��ȫ��ʱ�䲽����ʱ�����壩
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! ���㲽��
!    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 end subroutine NS_Time_advance_1Euler






!----------------------------------------------------------------------------------------
! ����3��RK�����ƽ�1��ʱ�䲽 ����nMesh������ �ĵ�������
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
      
	  call Comput_Residual_one_mesh(nMesh)              ! ����в�
      if(nMesh .ne. 1) call Add_force_function(nMesh)   ! ���ǿ�Ⱥ�������������Ĵ�����ʹ�ã�
   
	  do mBlock=1,Mesh(nMesh)%Num_Block
      B => Mesh(nMesh)%Block(mBlock)                  ! ��nMesh ������ĵ�mBlock��
      nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!    ʱ���ƽ� 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m,du)
      do j=1,ny-1
      do i=1,nx-1
      do m=1,MP%NVAR
		du=B%Res(m,i,j)/B%vol(i,j)  
        B%U(m,i,j)=Ralfa(KRK)*B%Un(m,i,j)+Rgamma(KRK)*B%U(m,i,j)+B%dt(i,j)*Rbeta(KRK)*du        ! 3��RK
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO

    enddo    
!---------------------------------------------------------------------------------------

      call Boundary_condition_onemesh(nMesh)         ! �߽����� ���趨Ghost Cell��ֵ��
      call update_buffer_onemesh(nMesh)              ! ͬ������Ľ�����

  enddo
     
	 Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! ʱ�� ��ʹ��ȫ��ʱ�䲽����ʱ�����壩
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1      ! ���㲽��
!    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 end subroutine NS_Time_advance_RK3






!---------------------------------------------------------------------------------------------
! ����LU_SGS������ʱ���ƽ�һ��ʱ�䲽 ����nMesh������ �ĵ�������
 subroutine NS_Time_advance_LU_SGS(nMesh)
   use Global_var
   implicit none
   integer::nMesh,mBlock,i,j,m,nx,ny
   Type (Mesh_TYPE),pointer:: MP
   Type (Block_TYPE),pointer:: B
   real*8:: alfa1
   
    call Comput_Residual_one_mesh(nMesh)     ! ���������ϼ���в�
    alfa1=0.d0                               
    MP=>Mesh(nMesh)
    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock,alfa1)                          ! ����LU_SGS��������DU=U(n+1)-U(n)
       B => MP%Block(mBlock)
       nx=B%nx; ny=B%ny
!--------------------------------------------------------------------------------------
!   ʱ���ƽ� 

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
      call Boundary_condition_onemesh(nMesh)         ! �߽����� ���趨Ghost Cell��ֵ��
      call update_buffer_onemesh(nMesh)              ! ͬ������Ľ�����


     Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! ʱ�� ��ʹ��ȫ��ʱ�䲽����ʱ�����壩
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! ���㲽��
!    print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 end subroutine NS_Time_advance_LU_SGS


!------------------------------------------------------------------------------------
! ˫ʱ�䲽 LU_SGS����
! ����Dual time step LU_SGS������ʱ���ƽ�һ��ʱ�䲽 
! Ŀǰֻ֧�ֵ�������
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
		alfa1=1.d0/dt_global         ! ����Ե�1��ʱ�䲽
	   else
	    alfa1=3.d0/(2.d0*dt_global)
	   endif
   
   nMesh=1
  
  do kt_in=1, Nstep_inner_Limit    ! ��ѭ������
   
    call Comput_Residual_one_mesh(nMesh)     ! ���������ϼ���в�
    MP=>Mesh(nMesh)

    do mBlock=1,MP%Num_Block
      call  du_LU_SGS_2D(nMesh,mBlock,alfa1)                          ! ����LU_SGS��������DU=U(n+1)-U(n)
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

      call Boundary_condition_onemesh(nMesh)         ! �߽����� ���趨Ghost Cell��ֵ��
      call update_buffer_onemesh(nMesh)              ! ͬ������Ľ�����

     max_res=MP%Res_rms(1)       ! ���������в�
     do m=1,NVAR
	 max_res=max(max_res,MP%Res_rms(m))
	 enddo
     if( max_res .le. Res_Inner_Limit) exit   ! �ﵽ�в��׼�������ڵ���
  enddo
 	
	 print*, "Inner step ... ", kt_in
	 print*, "rms residual eq =", MP%Res_rms(1:NVAR)


! �ڵ�������������U(n),U(n-1)    
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



     Mesh(nMesh)%tt=Mesh(nMesh)%tt+dt_global      ! ʱ�� ��ʹ��ȫ��ʱ�䲽����ʱ�����壩
     Mesh(nMesh)%Kstep=Mesh(nMesh)%Kstep+1        ! ���㲽��

     print*, "----------------------------------------------------"
     print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
 
 end subroutine Dual_time_LU_SGS












! ���㣨���أ�ʱ�䲽��
   subroutine comput_dt(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real*8,parameter:: C0=1
     Type (Block_TYPE),pointer:: B
    
     B => Mesh(nMesh)%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
     do j=1,ny-1
     do i=1,nx-1
      if(Iflag_local_dt .eq. 0) then   ! ȫ��ʱ�䲽��
        B%dt(i,j)=dt_global                                 
      else                             ! ����ʱ�䲽��
         B%dt(i,j)=CFL*B%vol(i,j)/(B%Lci(i,j)+B%Lcj(i,j) +C0*(B%Lvi(i,j)+B%Lvj(i,j)) )     ! �ֲ�ʱ�䲽��
         if(B%dt(i,j) .gt. dtmax) B%dt(i,j)=dtmax
         if(B%dt(i,j) .lt. dtmin) B%dt(i,j)=dtmin
      endif
	 enddo
	 enddo
!$OMP END PARALLEL DO 
  end subroutine comput_dt

!----------------------------------------------------------
! ������ճ�Լ�ճ������װ뾶���ڼ��������������ֲ�ʱ�䲽��������ʽ���в��˳�ȣ���ʹ��
   subroutine comput_Lvc(nMesh,mBlock)
     use Global_Var
     use Flow_Var 
     implicit none
	 integer  nMesh,mBlock,nx,ny,i,j
     real*8:: C0,si,sj,s0,ni1,ni2,nj1,nj2,uni,unj,tmp1
     Type (Block_TYPE),pointer:: B
    
	 C0=1.d0
     B => Mesh(nMesh)%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,si,sj,s0,ni1,ni2,nj1,nj2,uni,unj,tmp1)

     do j=1,ny-1
     do i=1,nx-1
 		 s0 =B%vol(i,j)         ! ���
		 si=0.5d0*(B%si(i,j)+B%si(i+1,j))         ! ������߳�������ƽ����
		 ni1=0.5d0*(B%ni1(i,j)+B%ni1(i+1,j))      ! ���淨��������ƽ����
		 ni2=0.5d0*(B%ni2(i,j)+B%ni2(i+1,j))
         sj=0.5d0*(B%sj(i,j)+B%sj(i,j+1))
		 nj1=0.5d0*(B%nj1(i,j)+B%nj1(i,j+1))
		 nj2=0.5d0*(B%nj2(i,j)+B%nj2(i,j+1))
         uni=uu(i,j)*ni1+v(i,j)*ni2             !�����ٶ�
		 unj=uu(i,j)*nj1+v(i,j)*nj2

!        �װ뾶
		 B%Lci(i,j)=(abs(uni)+cc(i,j))*si          ! ��ճ��Jocabian������װ뾶 ��Blazek's Book 6.1.4�ڣ�
		 B%Lcj(i,j)=(abs(unj)+cc(i,j))*sj

		 tmp1=gamma/d(i,j)*(B%Amu(i,j)/Pr+B%Amu_t(i,j)/PrT)
		 B%Lvi(i,j)=tmp1*si*si/s0         ! ճ����Jocabian�����װ뾶
         B%Lvj(i,j)=tmp1*sj*sj/s0
 	 enddo
	 enddo
!$OMP END PARALLEL DO
  end subroutine comput_Lvc
