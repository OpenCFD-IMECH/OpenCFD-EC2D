! Copyright By Li Xinliang, Code by Li Xinliang
! ��ʼ���������������ݽṹ������ֵ
! ���ڶ������񣬸����ϼ��������Ϣ��������������
!---------------------------------------------------------------------------------
!------------------------------------------------------------------------------     
   subroutine init
   use Global_var
   implicit none
   integer :: i,j,k,m,nx1,ny1,Num_Block1,ksub,Kmax_grid
   real*8,allocatable,dimension(:,:):: x1,y1
   integer,allocatable,dimension(:):: NI,NJ
   Type (Block_TYPE),pointer:: B
   TYPE (BC_MSG_TYPE),pointer:: Bc
 !--------------------------------------------------------------------
 ! initial of const variables
     Cv=1.d0/(gamma*(gamma-1.d0)*Ma*Ma)

!--------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
   allocate( Mesh(Num_Mesh) )             ! �����ݽṹ�� ������ �����Ա�ǡ�����顱��
   call Creat_Mesh1            ! ������ϸ������ (�������ļ�Mesh2d.dat)
   call read_bcin              ! ������������Ϣ (bc2d.in)
   if(Num_Mesh .ge. 2) then
     call Creat_Mesh(1,2)    ! ����1��������ϸ������Ϣ������2�����񣨴�����
   endif
   if(Num_Mesh .ge. 3) then
     call Creat_Mesh(2,3)    ! ����2��������Ϣ�������񣩣� ����3�������������
   endif


 end   
!---------------------------------------------------------



!  ��������Ƿ������ڶ������� 
!  ������������= 2*K+1 ����2������=4*K+1 ����3������=8*K+1 ����4������ ...
  
     subroutine check_mesh_multigrid 
	 use Global_var
	 implicit none
     integer,allocatable,dimension(:):: NI,NJ
     integer:: NB,NST,m,k,NN,Km,Km_grid,N_Cell,Ntmp,Bsub,ksub
     integer:: f_no, face, ist, iend, jst, jend,  neighb, subface, orient
	 print*, "Read Mesh2d.dat, Check if Multi-Grid can be used ..."
	 open(99,file="Mesh2d.dat")
     read(99,*) NB
	 allocate(NI(NB),NJ(NB))
     read(99,*) (NI(m), NJ(m), m=1,NB)
     close(99)
	 N_Cell=0
	 Km_grid=NI(1)  ! ��ʼֵ
     
	 do m=1,NB 
	   N_Cell=N_Cell+(NI(m)-1)*(NJ(m)-1)  ! ͳ��������Ԫ��
 
 !    �жϿ�ʹ�õ���������      
 	    Km=1
	    NN=2
!                                �ж�׼�� ������-1 �ܱ�2**km ������ ����ϡ������Ԫ����С��2
       do while( mod((NI(m)-1),NN) .eq. 0 .and. (NI(m)-1)/NN .ge. 2     &
		     .and. mod((NJ(m)-1),NN) .eq. 0 .and. (NJ(m)-1)/NN .ge. 2 ) 
        Km=Km+1              ! ���������������
	    NN=NN*2
       enddo
       Km_grid=min(Km_grid,Km)
    enddo
    Print*, " Finished check Mesh2d.dat,  Most stage is ", Km_grid
    print*,  "Check bc2d.in ..." 
    open(88,file="bc2d_1.in")
    read(88,*)
    read(88,*)
    read(88,*) Ntmp

    do m=1,NB
     read(88,*)
     read(88,*)
     read(88,*) Bsub    !number of the subface in the Block m
     read(88,*)
     do ksub=1, Bsub
      read(88,*)  f_no, face, ist, iend, jst, jend,  neighb, subface, orient
	     NN=1
		 Km=1
         do while( mod((ist-1),NN) .eq. 0 .and. mod((iend-1),NN) .eq.0       &
		     .and. mod((jst-1),NN) .eq. 0 .and. mod((jend-1),NN) .eq. 0  ) 
          NN=NN*2
		  Km=Km+1
		 enddo
        Km_grid=min(Km_grid,Km)
     enddo
     enddo


     close(88)


	 print*, "Check multigrid OK"
	 print*, "Total Block number is ", NB, "Total Cell number is " , N_Cell
	 print*, "Most stage number of multi-grid is ", Km_grid
	 print*, "--------------------------------------------------------------"
!-------------------------------------------------------
      Num_Mesh=min(3,Km_grid,Num_Mesh)    ! �趨�������������汾�������3������
      print*, Num_Mesh, "  stage grids is used !"
!-------------------------------------------------------
	     
	 deallocate(NI,NJ)
    end subroutine check_mesh_multigrid

!--------------------------------------------------------------------------------------
!   �������ݽṹ�� ��ϸ���� �����漸�������غ������
   subroutine Creat_Mesh1
   use Global_var
   implicit none
   integer,allocatable,dimension(:):: NI,NJ
   integer:: NB,m,nx,ny,i,j,k
   Type (Block_TYPE),pointer:: B
   real*8:: dx,dy

   print*, "-------------------------------------"
   print*, "read Mesh2d.dat"
   open(99,file="Mesh2d.dat")
   read(99,*) NB   ! Block number
     Mesh(1)%Num_Block=NB
     Mesh(1)%Num_Cell=0
     allocate(Mesh(1)%Block(NB))
     allocate( NI(NB),NJ(NB) )
     read(99,*) (NI(k), NJ(k), k=1,NB)
    do m=1,NB
     B => Mesh(1)%Block(m)
     B%Block_no=m
	 B%FVM_FDM=METHOD_FVM   ! ��ֵ������Ĭ��Ϊ���������
     B%nx=NI(m); B%ny=NJ(m) 
     nx=B%nx ; ny= B%ny
     Mesh(1)%Num_Cell=Mesh(1)%Num_Cell+(nx-1)*(ny-1)

!   �����ڴ�   (x,y) �ڵ����ꣻ (x1,y1)������������; s0 ����������� U, Un �غ����
   
     allocate(B%x(1-LAP:nx+LAP,1-LAP:ny+LAP), B%y(1-LAP:nx+LAP,1-LAP:ny+LAP))   
     allocate(B%x1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B%y1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
     allocate(B%U(NVAR,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))    !�غ���� (4�������������6��������4���������+k+w) 
     allocate(B%deltU(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))   !����֮��ɴ������ֵ����,��������ʹ�ã�4��������k,w���̲����� 

	 allocate(B%Res(NVAR,1:nx-1,1:ny-1))        ! �в�
	 allocate(B%dt(1:nx-1,1:ny-1))              ! ʱ�䲽��
     
	 allocate(B%dU(NVAR,nx,ny))                 ! =U(n+1)-U(n), LU-SGS��ʹ��
     allocate(B%Amu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), B%Amu_t(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))

!  ������
     allocate(B%vol(nx,ny),B%si(nx,ny),B%sj(nx,ny),B%ni1(nx,ny),B%ni2(nx,ny),B%nj1(nx,ny),B%nj2(nx,ny)) 
     allocate(B%Lci(nx,ny),B%Lcj(nx,ny),B%Lvi(nx,ny),B%Lvj(nx,ny))        ! �װ뾶 

     allocate(B%Un(NVAR,0:nx,0:ny))             ! ��һʱ�䲽��ֵ
     if( Time_Method .eq. Time_dual_LU_SGS ) then
	  allocate(B%Un1(NVAR,0:nx,0:ny))    ! n-1ʱ�䲽��ֵ�� ˫ʱ�䲽LU_SGS�����в���
     endif


!   ��ʼ��	   
       B%x1(:,:)=0.d0; B%y1(:,:)=0.d0; B%vol(:,:)=0.d0                     
       B%U(1,:,:)=1.d0; B%U(2,:,:)=0.d0; B%U(3,:,:)=0.d0; B%U(4,:,:)=1.d0
       B%Un(1,:,:)=1.d0; B%Un(2,:,:)=0.d0; B%Un(3,:,:)=0.d0; B%Un(4,:,:)=0.d0
       B%Lci(:,:)=0.d0; B%Lcj(:,:)=0.d0; B%Lvi(:,:)=0.d0; B%Lvj(:,:)=0.d0
       B%dU(:,:,:)=0.d0
	   B%Amu(:,:)=0.d0;  B%Amu_t(:,:)=0.d0
       B%DeltU=0.d0
       B%Res=0.d0
       B%dt=0.d0
       if(Nvar .eq. 5) then
	     B%U(5,:,:)=1.d0/Re
       else if(NVar .eq. 6) then
	    B%U(5,:,:)=0.d0          ! �Ķ���
        B%U(6,:,:)=1.d0          ! ���ܱȺ�ɢ��
       endif
       
     read(99,*) ((B%x(i,j),i=1,nx),j=1,ny),  ((B%y(i,j),i=1,nx),j=1,ny)

!   ������������� 
!   ����������������Ҫ�Ȼ��Ghost Cell ����Ϣ�����       
	   do j=1,ny-1
       do i=1,nx-1
         B%vol(i,j) =abs((B%x(i,j)-B%x(i+1,j+1))*(B%y(i+1,j)-B%y(i,j+1)) -  &
           (B%x(i+1,j)-B%x(i,j+1))*(B%y(i,j)-B%y(i+1,j+1)) )*0.5d0           ! ����������������
       enddo
	   enddo

!     ������ ���߳���������      
	  do j=1,ny-1
      do i=1,nx
       dx=B%x(i,j+1)-B%x(i,j)
       dy=B%y(i,j+1)-B%y(i,j)
       B%si(i,j)=sqrt(dx*dx+dy*dy)
       B%ni1(i,j)=dy/B%si(i,j); B%ni2(i,j)=-dx/B%si(i,j)   ! normal vector at (i,j) or (I-1/2,J) 
	  enddo
	  enddo

      do j=1,ny
      do i=1,nx-1
       dx=B%x(i+1,j)-B%x(i,j)
       dy=B%y(i+1,j)-B%y(i,j)
       B%sj(i,j)=sqrt(dx*dx+dy*dy)     ! length 
       B%nj1(i,j)=-dy/B%sj(i,j); B%nj2(i,j)=dx/B%sj(i,j)     ! normal vector at i, j+1/2
      enddo
	  enddo
 
	
	 enddo
     close(99) 
	 deallocate(NI,NJ)

!    ʱ�䲽��ʱ��
       Mesh(1)%Kstep=0
       Mesh(1)%tt=0.d0

	 print*, "read Mesh2d.dat OK"
     end subroutine Creat_Mesh1


!----Mesh control message (bc2d.in)------------------------------------------
    subroutine read_bcin 
    use Global_Var
    implicit none
    integer:: NB,m,ksub
    Type (Block_TYPE),pointer:: B
    TYPE (BC_MSG_TYPE),pointer:: Bc

    print*, "read bc2d.in ......"
    open(88,file="bc2d_1.in")
    read(88,*)
    read(88,*)
    read(88,*) NB
    if(NB .ne. Mesh(1)%Num_Block) then
      print*, "Error!  Block number in bc2d.in is not equal to that in Mesh2d.dat !"
      stop
    endif

    do m=1,NB
     B => Mesh(1)%Block(m)
     read(88,*)
     read(88,*)
     read(88,*) B%subface   !number of the subface in the Block m
     read(88,*)
     allocate(B%bc_msg(B%subface))
     do ksub=1, B%subface
     Bc => B%bc_msg(ksub)
     read(88,*)  Bc%f_no, Bc%face, Bc%ist, Bc%iend, Bc%jst, Bc%jend,  Bc%neighb, Bc%subface, Bc%orient
     enddo
     enddo
     close(88)
	 print*, "read bc2d.in OK"
   end  subroutine read_bcin

!-----------------------------------------------------------------------------------
! �����ϼ�������Ϣ������������m2 (ϡ������)
! ϡ������ʹ������ģ�ͣ�������Ϊ4 
   subroutine Creat_Mesh(m1,m2)
    use Global_Var
    implicit none
    integer:: NB,m,m1,m2,ksub,nx,ny,i,j,i1,j1,Bsub
    Type (Block_TYPE),pointer:: B1,B2
    TYPE (BC_MSG_TYPE),pointer:: Bc1,Bc2
    Type (Mesh_TYPE),pointer:: MP1,MP2
	real*8:: dx,dy
    print*,  "Creat Mesh ......", m2

    MP1=>Mesh(m1)             ! ��һ������ ��ϸ����
	Mp2=>Mesh(m2)             ! ��������   ��������

    NB=MP1%Num_Block
	MP2%Num_Block=NB      !  ����m2��m1 ������ͬ
    MP2%Mesh_no=m2        ! �����
    MP2%Num_Cell=0      

	allocate(MP2%Block(NB))   ! ��MP2�д������ݽṹ�����顱

	do m=1,NB
     B1=>MP1%Block(m)
     B2=>MP2%Block(m)
	 B2%Block_no=m
 	 B2%FVM_FDM=METHOD_FVM   ! ��ֵ������������ֻ��ʹ��������������߾��ȷ���ֻ�������������У�
 
     nx=(B1%nx-1)/2+1          ! ������ĵ���
	 ny=(B1%ny-1)/2+1
	 B2%nx=nx
	 B2%ny=ny
     
	 MP2%Num_Cell=MP2%Num_Cell+(nx-1)*(ny-1)          ! ͳ��MP2��������Ԫ��

!    ������������������
!--------------------------------------------------
     allocate(B2%x(1-LAP:nx+LAP,1-LAP:ny+LAP), B2%y(1-LAP:nx+LAP,1-LAP:ny+LAP))   
     allocate(B2%x1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1),B2%y1(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1)) 
     allocate(B2%U(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))   !�غ����  ����ʹ������ģ�ͣ�����k,w, �ܹ�4��������
     allocate(B2%deltU(4,1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))  ! ��ʱ�䲽U�Ĳ�ֵ ����������ʹ�ã��Ӵ������ֵ��������
	 allocate(B2%Un(4,0:nx,0:ny))   ! ��һʱ�䲽��ֵ
	 allocate(B2%Res(4,1:nx-1,1:ny-1))        ! �в�
	 allocate(B2%dt(1:nx-1,1:ny-1))           ! ʱ�䲽��
     allocate(B2%QF(4,1:nx-1,1:ny-1))         ! ǿ�Ⱥ���
	 allocate(B2%dU(4,nx,ny))   ! LU-SGS��ʹ��
     allocate(B2%Amu(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1), B2%Amu_t(1-LAP:nx+LAP-1,1-LAP:ny+LAP-1))
!    ������
     allocate(B2%vol(nx,ny),B2%si(nx,ny),B2%sj(nx,ny),B2%ni1(nx,ny),B2%ni2(nx,ny),B2%nj1(nx,ny),B2%nj2(nx,ny)) 
     allocate(B2%Lci(nx,ny),B2%Lcj(nx,ny),B2%Lvi(nx,ny),B2%Lvj(nx,ny))        ! �װ뾶 
    
	 if( Time_Method .eq. Time_dual_LU_SGS ) then
	  allocate(B2%Un1(NVAR,0:nx,0:ny))    ! n-1ʱ�䲽��ֵ�� ˫ʱ�䲽LU_SGS�����в���
     endif

!    ��ʼ��
       B2%x1(:,:)=0.d0; B2%y1(:,:)=0.d0; B2%vol(:,:)=0.d0                     
       B2%U(1,:,:)=1.d0; B2%U(2,:,:)=0.d0; B2%U(3,:,:)=0.d0; B2%U(4,:,:)=1.d0
       B2%Un(1,:,:)=1.d0; B2%Un(2,:,:)=0.d0; B2%Un(3,:,:)=0.d0; B2%Un(4,:,:)=1.d0
       B2%Res(:,:,:)=0.d0
	   B2%dt(:,:)=0.d0
	   B2%QF(:,:,:)=0.d0                       ! ǿ�Ⱥ�����ʼ��Ϊ0
       B2%Lci(:,:)=0.d0; B2%Lcj(:,:)=0.d0; B2%Lvi(:,:)=0.d0; B2%Lvj(:,:)=0.d0
       B2%dU(:,:,:)=0.d0
	   B2%Amu(:,:)=0.d0;  B2%Amu_t(:,:)=0.d0
       
!   �趨������Ϣ�����ݴ֡�ϸ����Ķ�Ӧ��ϵ��     
	  do j=1,ny
	  do i=1,nx
	    i1=2*i-1 ; j1=2*j-1
	    B2%x(i,j)=B1%x(i1,j1)         !��������ϸ����Ķ�Ӧ��ϵ ����һ��������һ��������㣩
        B2%y(i,j)=B1%y(i1,j1)
	  enddo
	  enddo
!  �����������
	   do j=1,ny-1
       do i=1,nx-1
         B2%vol(i,j) =abs((B2%x(i,j)-B2%x(i+1,j+1))*(B2%y(i+1,j)-B2%y(i,j+1)) -  &
           (B2%x(i+1,j)-B2%x(i,j+1))*(B2%y(i,j)-B2%y(i+1,j+1)) )*0.5d0           ! ����������������
       enddo
	   enddo

!     ������ ���߳���������      
	  do j=1,ny-1
      do i=1,nx
       dx=B2%x(i,j+1)-B2%x(i,j)
       dy=B2%y(i,j+1)-B2%y(i,j)
       B2%si(i,j)=sqrt(dx*dx+dy*dy)
       B2%ni1(i,j)=dy/B2%si(i,j); B2%ni2(i,j)=-dx/B2%si(i,j)   ! normal vector at (i,j) or (I-1/2,J) 
	  enddo
	  enddo

      do j=1,ny
      do i=1,nx-1
       dx=B2%x(i+1,j)-B2%x(i,j)
       dy=B2%y(i+1,j)-B2%y(i,j)
       B2%sj(i,j)=sqrt(dx*dx+dy*dy)     ! length 
       B2%nj1(i,j)=-dy/B2%sj(i,j); B2%nj2(i,j)=dx/B2%sj(i,j)     ! normal vector at i, j+1/2
      enddo
	  enddo


!-----------------------------------------------  
!    ����������Ϣ
      Bsub=B1%subface        ! ������
      B2%subface=Bsub
      allocate(B2%bc_msg(Bsub))
      do ksub=1, Bsub
	   Bc1=> B1%bc_msg(ksub)    ! ��һ�������������Ϣ
	   Bc2=> B2%bc_msg(ksub)    ! ���������������Ϣ
       Bc2%f_no=Bc1%f_no
	   Bc2%face=Bc1%face
	   Bc2%ist=(Bc1%ist-1)/2+1   ! �֡�ϸ�����±�Ķ�Ӧ��ϵ
	   Bc2%iend=(Bc1%iend-1)/2+1
	   Bc2%jst=(Bc1%jst-1)/2+1
	   Bc2%jend=(Bc1%jend-1)/2+1
       Bc2%neighb= Bc1%neighb
       Bc2%subface=Bc1%subface
       Bc2%orient=Bc1%orient
      enddo

  
     enddo
       Mesh(m2)%Kstep=0
       Mesh(m2)%tt=0.d0
	  Print*, "Creat Mesh ",m2," OK ", "Total Cell number is", MP2%Num_Cell 
   end  subroutine Creat_Mesh
!--------------------------------------------------------------------------------
! ��ʼ������ 
    subroutine Init_flow 
     use Global_var
	 implicit none
     integer:: i,j,m,m1
     Type (Block_TYPE),pointer:: B
     Type (Mesh_TYPE),pointer:: MP
    
	 if(Iflag_init .eq. 0) then
	   call init_flow_zero                      ! ������������ ���ȴӴ�������㣬�ٲ�ֵ��ϸ����
     else
	   call init_flow_read                      ! ���ļ���ȡ����
     endif
    
	 if(Time_Method .eq. Time_Dual_LU_SGS) then
	    MP=>Mesh(1)             
       do m=1,MP%Num_Block
     
		B => MP%Block(m)    ! Mesh(Num_Mesh) ����ֵ�����
         do j=1,B%ny-1
         do i=1,B%nx-1
         do m1=1,NVAR
		  B%Un(m1,i,j)=B%U(m1,i,j)
		  B%Un1(m1,i,j)=B%U(m1,i,j) 
         enddo
         enddo
         enddo
		enddo
     endif 
	
	
	print*, " Initialize OK ......"
   end subroutine Init_flow







! ��������ʼ���� ������������£����������ʼ���㣨Ȼ���ֵ��ϸ����  
     subroutine init_flow_zero
     use Global_var
	 implicit none
     real*8:: d0,u0,v0,p0,T0,tmp
     integer:: i,j,m,step,nMesh
     Type (Block_TYPE),pointer:: B
     Type (Mesh_TYPE),pointer:: MP
	 
!-------------------------------------------------------------------
	  d0=1.d0; u0=1.d0*cos(AoA); v0=1.d0*sin(AoA) ; p0=1.d0/(gamma*Ma*Ma)
      if(Nvar .eq. 5) then
        MP=>Mesh(1)             ! ������
        do m=1,MP%Num_Block
         B => MP%Block(m)    ! Mesh(Num_Mesh) ����ֵ�����
         do j=1,B%ny-1
         do i=1,B%nx-1
          B%U(5,i,j)=vt_inf/Re   ! Vtֵ ��ͨ���趨Ϊ����ճ��ϵ����3-5����
         enddo
         enddo
        enddo
      endif
      
!  k��w�ĳ�ֵ ��0��ֵ������������������ֵ��
      if(Nvar .eq. 6) then
        MP=>Mesh(1)             ! ������
        do m=1,MP%Num_Block
         B => MP%Block(m)    ! Mesh(Num_Mesh) ����ֵ�����
         do j=1,B%ny-1
         do i=1,B%nx-1
          B%U(5,i,j)=10.d0*Kt_inf     ! �Ķ��ܣ� ��ֵΪ����ֵ��10�����������ڷ�չ��
          B%U(6,i,j)=Wt_inf    ! ���ܱȺ�ɢ�ʣ���ֵΪ����ֵ 
         enddo
         enddo
        enddo
      endif
!---------------------------------------------------------------------
  !  �������ĳ�ֵ���ɴ������𼶼����ֵ����
      
      MP=>Mesh(Num_Mesh)             ! ��ϡ�������
      do m=1,MP%Num_Block
      B => MP%Block(m)    ! Mesh(Num_Mesh) ����ֵ�����
      do j=1,B%ny-1
      do i=1,B%nx-1
       B%U(1,i,j)=d0
       B%U(2,i,j)=d0*u0
       B%U(3,i,j)=d0*v0
       B%U(4,i,j)=p0/(gamma-1.d0)+0.5d0*d0*(u0*u0+v0*v0)
      enddo
      enddo
      enddo
    call Boundary_condition_onemesh(Num_Mesh)                   ! �߽����� ���趨Ghost Cell��ֵ��
    call update_buffer_onemesh(Num_Mesh)                        ! ͬ������Ľ�����
!-----------------------------------------------------------------
!------------------------------------------------------
!   ׼����ֵ�Ĺ���
!   �����������㣬�𼶲�ֵ��ϸ����
  

  do nMesh=Num_Mesh,1,-1 
   do step=1, Pre_Step_Mesh(nMesh)      
    call NS_Time_advance(nMesh)
    if(mod(step,Kstep_show) .eq. 0) call output_Res(nMesh)
   enddo
!    call output (nMesh)
    if(nMesh .gt. 1) then
    call prolong_U(nMesh,nMesh-1,1)                            ! ��nMesh�������ϵ���������ֵ����һ������; flag=1 ��ֵU����
    call Boundary_condition_onemesh(nMesh-1)                   ! �߽����� ���趨Ghost Cell��ֵ��
    call update_buffer_onemesh(nMesh-1)                        ! ͬ������Ľ�����
   
    print*, " Prolong  to mesh ", nMesh-1, "   OK"           
    endif
  enddo
  end subroutine init_flow_zero 
 
!-----------------------------------------------
!  ��flow2d.dat �ļ���ȡ��ֵ �����ܵ����񣩣� ����Ϊtecplot��ʽ
     subroutine init_flow_read
     use Global_var
     real*8:: d0,u0,v0,p0,T0,kt0,kw0,tmp
     Type (Mesh_TYPE),pointer:: MP

     integer:: i,j,step,nMesh
     Type (Block_TYPE),pointer:: B
      MP=>Mesh(1)             ! ���ܵ�����
      print*, "Init from 'flow2d.dat' ......"

       open(99,file="flow2d.dat")
	   read(99,*)
       do m=1,MP%Num_Block
       B => MP%Block(m)                 ! �����
         read(99,*)
         do j=0,B%ny
	     do i=0,B%nx
	     read(99,*) tmp,tmp,d0,u0,v0,T0
         B%U(1,i,j)=d0                    ! �ܶ�
         B%U(2,i,j)=d0*u0                 ! x-�������ܶ�
         B%U(3,i,j)=d0*v0                 ! y-�������ܶ�
         B%U(4,i,j)=d0*Cv*T0+0.5d0*d0*(u0*u0+v0*v0)  ! �������ܶ�
         enddo
         enddo
       enddo
	   close(99)
       
	   if(MP%NVAR .eq. 5) then
       open(99,file="SA2d.dat")   ! ��ȡvt
	   read(99,*)
       do m=1,MP%Num_Block
       B => MP%Block(m)                 
         read(99,*)
         do j=0,B%ny
	     do i=0,B%nx
	     read(99,*) tmp,tmp,B%U(5,i,j)
         enddo
         enddo
       enddo
	   close(99)
       endif

	   if(MP%NVAR .eq. 5) then
       open(99,file="SST2d.dat")   ! ��ȡKt,Wt
	   read(99,*)
       do m=1,MP%Num_Block
       B => MP%Block(m)                 
         read(99,*)
         do j=0,B%ny
	     do i=0,B%nx
	     read(99,*) tmp,tmp,B%U(5,i,j),B%U(6,i,j)
         enddo
         enddo
       enddo
	   close(99)
       endif



      call Boundary_condition_onemesh(1)                   ! ����߽����� 
      call update_buffer_onemesh(1)                        ! �ڱ߽�����
    end  subroutine init_flow_read