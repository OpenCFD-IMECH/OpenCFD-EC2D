!----------------------------------------------------------------------
! Copyright by Li Xinliang, Code by Li Xinliang

!----------------------------------------------------------------------
! �����������N-S���� ���ƽ�1��ʱ�䲽��
!  nMesh=1,2,3 �ֱ��Ӧ��ϸ���񡢴����񡢸�������
! ����2�������3�����������ӳ���
!---------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! �����������ƽ�1��ʱ�䲽 (3��RK or 1th Euler)
  subroutine NS_2stge_multigrid
   use Global_var
   implicit none
   integer::nMesh,m
   Type (Block_TYPE),pointer:: B
   integer,parameter:: Time_step_coarse_mesh=3     ! �������������

!---------------------------------------------------
! ----  ����1 -----------------
      if(Time_Method .eq. Time_Euler1) then
	     call  NS_Time_advance_1Euler(1)                 ! ϸ����1��Euler�����ƽ�1�� -> U(n+1)
      else
	     call  NS_Time_advance_RK3(1)                    ! ϸ����RK�����ƽ�1�� -> U(n+1)
      endif

! ???????????     
	  call  Comput_Residual_one_mesh(1)           ! ��������1�Ĳв� R(n+1)   ! ???????? �ƺ�����ȡ���ò�
	  
	  call  interpolation2h(1,2,2)                ! �Ѳв��ֵ������2 (������QF����)
	  call  interpolation2h(1,2,1)                ! ���غ����������1��ֵ������2   ��flag=1 ��ֵ�غ������=2 ��ֵ�в

!------------------------------
      call  Boundary_condition_onemesh(2)         ! ����߽�����
      call  update_buffer_onemesh(2)              ! �ڱ߽�����
      call  Comput_Residual_one_mesh(2)           ! ��������2�Ĳв�
      call  comput_force_function(2)              ! ����ǿ�Ⱥ���QF
	  
      if(Time_Method .eq. Time_Euler1) then
	    call Set_Un(2)                               ! ��¼��ʼֵ  ��RK�������Ѿ������˸ò��� 
        do m=1, Time_step_coarse_mesh
	     call  NS_Time_advance_1Euler(2)             ! 1��Euler�������ɲ�
		enddo
	   else 
		 call  NS_Time_advance_RK3(2)                ! RK�����ƽ�1�� ������2��
       endif

	  call  comput_delt_U(2)                      ! ����������deltU ��������Un���棩
      call  prolong_U(2,1,2)                      ! ����������ֵ��ϸ���� (������Un����); flag=2 ��ֵdeltU (������Un��)
!------------------------------------	 
	  call  comput_new_U(1)                       ! �����µ�U  (U=U+deltU)
      call  Boundary_condition_onemesh(1)         ! ����߽�����
      call  update_buffer_onemesh(1)              ! �ڱ߽�����

!     print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
  end subroutine NS_2stge_multigrid

!------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! ���������ϵ���1��ʱ�䲽 ��V-�͵����� 3��RK or 1��Euler
  subroutine NS_3stge_multigrid
   use Global_var
   implicit none
   integer::nMesh,m
   integer,parameter:: Time_step_coarse_mesh=3     ! ������������� (��1��Euler��Ч)

!------------------
       Type (Block_TYPE),pointer:: B
       integer:: i,j
!------------------------

!---------------------------------------------------
! ----  ����1 -----------------
      if(Time_Method .eq. Time_Euler1) then
	     call  NS_Time_advance_1Euler(1)                 ! ϸ����1��Euler�����ƽ�1�� -> U(n+1)
      else
	     call  NS_Time_advance_RK3(1)                    ! ϸ����RK�����ƽ�1�� -> U(n+1)
      endif

! ?????????????????????????????
	  call  Comput_Residual_one_mesh(1)            ! ��������1�Ĳв� R(n+1)   !    
! -----------------------------  
	  
	  call  interpolation2h(1,2,2)                 ! �Ѳв��ֵ������2 (����������2��QF����)
	  call  interpolation2h(1,2,1)                 ! ���غ����������1��ֵ������2 �����浽U���棩  ��flag=1 ��ֵ�غ������=2 ��ֵ�в
      
!-------����2 --------------------
      call  Boundary_condition_onemesh(2)          ! ����߽�����
      call  update_buffer_onemesh(2)               ! �ڱ߽�����
      call  Comput_Residual_one_mesh(2)            ! ��������2�Ĳв�         Res_2h(0)
      call  comput_force_function(2)               ! ����ǿ�Ⱥ���QF ������2��QF_2h=QF_2h-Res_2h(0)
	  
      if(Time_Method .eq. Time_Euler1) then
	    call Set_Un(2)                               ! ��¼��ʼֵ  ��RK�������Ѿ������˸ò��� 
        do m=1, Time_step_coarse_mesh
	     call  NS_Time_advance_1Euler(2)             ! 1��Euler�������ɲ�
		enddo
	   else 
		 call  NS_Time_advance_RK3(2)                ! RK�����ƽ�1�� ������2��
       endif

      call  Comput_Residual_one_mesh(2)             ! ��������2�Ĳв� R_2h(n+1)  
      call  Add_force_function(2)                   ! �����ǿ�Ȳв�(������Res����)  RF_2h(n+1)=R_2h(n+1)+QF_2h  ;  Ŀ�ģ���ֵ������3��
	  call  interpolation2h(2,3,2)                  ! �Ѳв��ֵ������3 (����������3��QF����)

	  call  interpolation2h(2,3,1)                 ! ���غ����������2��ֵ������3 �����浽U���棩  ��flag=1 ��ֵ�غ������=2 ��ֵ�в

!------����3----------------------	  
      call  Boundary_condition_onemesh(3)          ! �߽�����: ����߽� 
      call  update_buffer_onemesh(3)               ! �ڱ߽�
	  call  Comput_Residual_one_mesh(3)           ! ��������3�Ĳв�
	
      call  comput_force_function(3)              ! ����ǿ�Ⱥ���QF ������3��
	  
      if(Time_Method .eq. Time_Euler1) then
	    call Set_Un(3)
        do m=1, Time_step_coarse_mesh
	     call  NS_Time_advance_1Euler(3)              ! 1��Euler�������ɲ�
		enddo
	   else 
		 call  NS_Time_advance_RK3(3)                 ! RK�����ƽ�1�� ������3��
       endif
	  
	  call  comput_delt_U(3)                      ! ����������deltU (=U-Un)
	  call  prolong_U(3,2,2)                      ! ����������ֵ������2 (������deltU����); flag=2 ��ֵdeltU 

!------����2------------------------      
	  call  comput_new_U(2)                       ! ����2�����µ�U  (U=U+deltU)
	  call  comput_delt_U(2)                      ! ����������deltU =U-Un
      call  prolong_U(2,1,2)                      ! ����������ֵ��ϸ���� (������deltU����); flag=2 ��ֵdeltU 
!------����1------------------------------
	  call  comput_new_U(1)                       ! �����µ�U  (U=U+deltU)
      call Boundary_condition_onemesh(1)          ! ����߽�����
      call update_buffer_onemesh(1)               ! �ڱ߽�����

!     print*, "Kstep, t=", Mesh(nMesh)%Kstep, Mesh(nMesh)%tt
  end subroutine NS_3stge_multigrid

!------------------------------------------------------------------------------------------








!------------------------------------------------------------------------------------------
  
!  ����ǿ�Ⱥ��� QF=Ih_to_2h Res(n-1) - Res(n)      ! QF�д�����ϸ�����ֵ�����Ĳв�
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
	   B%QF(m,i,j)=B%QF(m,i,j)-B%Res(m,i,j)            ! QFԭ�ȴ����Ŵ�ϸ�����ֵ�����Ĳв�
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO
	   enddo
  end  subroutine comput_force_function

!------------------------------------------------------------
!  ��ǿ�Ⱥ�����ӵ��в��� RF=R+QF        
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
	   B%Res(m,i,j)=B%Res(m,i,j)+B%QF(m,i,j)            ! ���ǿ�Ⱥ�����Ĳв��Դ�����B%Res���� ����ʡ�ڴ棩
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO

	   enddo
  end  subroutine Add_force_function

!----------------------------------------------------------------------  
!  ���������� deltU=U-Un
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
! �趨Un=U
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

!-------------------------------����U --------------------------------------
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
!	   B%U(m,i,j)=B%U(m,i,j)+B %deltU(m,i,j)         !deltU���洢�����U�������� ���Ӵ������ֵ������
 	   B%U(m,i,j)=B%U(m,i,j)+B1%deltU(m,i,j)           ! �޸�2011-7-6 ,�޸ĺ�ʹ��Intel Fortran 32λ����������ɢ ����δ�ҵ�ԭ�򣬿��ܻ�����������!!!��
       enddo
	   enddo
       enddo
!$OMP END PARALLEL DO

	   enddo
  end  subroutine comput_new_U
!---------------------------------------------------------------------------


!----------------------------------------------------------------------
! ��������ϸ����Ĳ�ֵ(Prolong) �� ϸ������������ϲ�ֵ (interpolation)
!----------------------------------------------------------------------
! ������m1���غ����(U) ��U�Ĳ��ֵ������m2 (��һ��ϸ����)
! flag=1ʱ����U��ֵ����һ������  (׼����ֵʱʹ��)
! flag=2ʱ����deltU��ֵ����һ������ ��deltU�����ű�ʱ�䲽���ϸ�ʱ�䲽U�Ĳ 

   Subroutine prolong_U(m1,m2,flag)
   use Global_Var
   implicit none
   integer:: m1,m2,mb,flag
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   integer,allocatable,dimension(:,:):: ia,ja
   real*8:: a1=9.d0/16.d0,a2=3.d0/16.d0,a3=1.d0/16.d0   ! ��ֵϵ��
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

!    Ѱ�Ҳ�ֵ���ܵ���±� 
!    ia(1,i) �Ǿ���i������Ĵ��������±ꣻia(2,i)�Ǵν�����±�	 
	 do i=0,nx2
	  if(mod(i,2).eq.0) then
	   ia(1,i)=i/2                    !�����
	   ia(2,i)=i/2+1                  !�ν���
	  else  
	   ia(1,i)=i/2+1                  !�����
	   ia(2,i)=i/2                    !�ν���
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
!     ��ֵ�غ����U

!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
 	 do j=0,ny2
	 do i=0,nx2
	 do m=1,4
!               ��ֵ��������Ȩ��a1, �ν����Ȩ��a2, ��Զ���Ȩ��a3	 
	  B2%U(m,i,j)=a1*B1%U(m,ia(1,i),ja(1,j))+a2*(B1%U(m,ia(2,i),ja(1,j))+B1%U(m,ia(1,i),ja(2,j)) ) &
	                                   +a3*B1%U(m,ia(2,i),ja(2,j))
     enddo
	 enddo
	 enddo
!$OMP END PARALLEL DO
	else
!  ��ֵdeltU
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,m)
 	 do j=0,ny2
	 do i=0,nx2
	 do m=1,4
!               ��ֵ��������Ȩ��a1, �ν����Ȩ��a2, ��Զ���Ȩ��a3	 
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
! ������m1���غ����U��ֵ������m2 (ϸ����->������) 
   Subroutine interpolation2h(m1,m2,flag)
   use Global_Var
   implicit none
   Type (Mesh_TYPE),pointer:: MP1,MP2
   Type (Block_TYPE),pointer:: B1,B2
   real*8,dimension(:,:,:),pointer:: P1,P2
   integer:: flag,m1,m2,mb,i,j,m,i1,i2,j1,j2
!   flag==1 ��ֵ�غ������ flag==2 ��ֵ�в�
	 if( m2-m1 .ne. 1) print*, "Error !!!!"
     MP1=>Mesh(m1)
     MP2=>Mesh(m2)
   
   do mb=1,MP1%Num_Block
     B1=>MP1%Block(mb)
  	 B2=>Mp2%Block(mb)
    if(flag .eq. 1) then  ! ��ֵ�غ����
	  P1=>B1%U
	  P2=>B2%U
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,i1,j1,i2,j2,m)
      do j=1,B2%ny-1
	  do i=1,B2%nx-1
	   i1=2*i-1 ; i2=2*i
	   j1=2*j-1 ; j2=2*j
      do m=1,4
!     �Կ��������ΪȨ�صļ�Ȩƽ�� 	 
	    P2(m,i,j)=(P1(m,i1,j1)*B1%vol(i1,j1)+P1(m,i1,j2)*B1%vol(i1,j2)   &
	      +P1(m,i2,j1)*B1%vol(i2,j1)+P1(m,i2,j2)*B1%vol(i2,j2))/B2%vol(i,j)
  	  enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO

    else     ! ��ֵ�в�  ����m1�����ϵĲв�B%Res ��ֵ��m2������B%QF (Ȼ���ȥ��m2�����ϵĲв�γ�ǿ�Ⱥ���)��
   	  P1=>B1%Res
	  P2=>B2%QF
!$OMP PARALLEL DO DEFAULT (SHARED) PRIVATE (i,j,i1,j1,i2,j2,m)
	  do j=1,B2%ny-1
	  do i=1,B2%nx-1
	   i1=2*i-1 ; i2=2*i
	   j1=2*j-1 ; j2=2*j
       do m=1,4
	     P2(m,i,j)=P1(m,i1,j1)+P1(m,i1,j2) +P1(m,i2,j1)+P1(m,i2,j2)    ! �в�Ĳ�ֵ�� �����
  	   enddo
	  enddo
	  enddo
!$OMP END PARALLEL DO
 	endif
   enddo

  end Subroutine interpolation2h
