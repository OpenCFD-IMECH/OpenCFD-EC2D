!  ����LU-SGS����������DU=U(n+1)-U(n)
!  2D Code by Li Xinliang, 2012-4-29
!-----------------------------------------------------------------------------------------
    subroutine  du_LU_SGS_2D(nMesh,mBlock,alfa1)                          ! ����LU_SGS��������DU=U(n+1)-U(n)
    use Global_Var
    use Flow_Var 
    implicit none
	integer:: nMesh,mBlock,NV,nx,ny,nz,plane,i,j,m

    real*8,dimension(6):: alfa,dui,duj,DF   ! alfa,�Խ���Ԫ��ֵ;  dui,duj,DF����ͨ��
	real*8,parameter:: w_LU=1.d0   !  LU_SGS���ɳ�����(1-2֮��), ����w_LU������ȶ��ԣ����ή�������ٶ�
	real*8:: alfa1                 ! ˫ʱ�䲽LU_SGS���ӶԽ���ֵ
	Type (Block_TYPE),pointer:: B
    TYPE (Mesh_TYPE),pointer:: MP
     MP=> Mesh(nMesh)
	 NV=MP%NVAR
	 B => MP%Block(mBlock)                 !��nMesh ������ĵ�mBlock��
     nx=B%nx; ny=B%ny

! LU-SGS������ɨ��
!----------------------------------
!   ��i=1,j=1 ��i=nx-1,j=ny-1��ɨ�����  (����ɨ�����)
!   ɨ�� i+j=plane ��ƽ��
!   w_LU���ɳ����ӣ�1��2֮�䣩������w_LU������ȶ��ԣ����ή�������ٶ�
   do plane=2,nx+ny-2            
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(plane,nx,ny,NV,B,gamma,If_viscous,alfa1)
   do j=1,ny-1
   i=plane-j
	 if( i .lt. 1 .or. i .gt. nx-1) cycle    ! ���������ƽ��
!     ����Խ���Ԫ�� ���翼��ճ�����ټ���ճ���װ뾶�� ��������ģ���еķ��̣�������Դ���װ뾶��
	  alfa(:)=B%vol(i,j)/B%dt(i,j)+w_LU*(B%Lci(i,j)+B%Lcj(i,j)) + B%vol(i,j)*alfa1                    ! �Խ�Ԫ
	  
	  if(If_viscous .eq. 1)  alfa(:)=alfa(:)+2.d0*(B%Lvi(i,j)+B%Lvj(i,j))           ! ����ճ����
       if(NV .eq. 5) then
	     alfa(5)=alfa(5)+B%Lvi(i,j)+B%Lvj(i,j)                !SAģ��
!       elseif(NV .eq. 6) then                                 !SSTģ��
!	      alfa(5)=alfa(5)+0.09*B%U(6,i,j)/B%U(1,i,j)
!         alfa(6)=alfa(6)+2.d0*0.0828*B%U(6,i,j)/B%U(1,i,j)
	   endif

	  		   
	 if(i.ne. 1) then
!                                                      ͨ���Ĳ������������Ƽ���A*W (See Blazek's book, page 208)
       call comput_DFn(NV,DF(1:NV),B%U(:,i-1,j),B%DU(:,i-1,j),B%ni1(i,j),B%ni2(i,j),gamma)  
       dui(1:NV)=0.5d0*(DF(1:NV)*B%si(i,j)+w_LU*B%Lci(i-1,j)*B%dU(:,i-1,j))
       if(If_viscous .eq. 1)    dui(1:NV)=dui(1:NV)+B%Lvi(i-1,j)*B%dU(:,i-1,j)         ! 2012-2-29, �ȶ��Ը���Щ
     else
	   dui=0.d0                             ! ���û�е�
     endif
	 
	 if(j.ne.1) then
       call comput_DFn(NV,DF(1:NV),B%U(:,i,j-1),B%DU(:,i,j-1),B%nj1(i,j),B%nj2(i,j),gamma)  
       duj(1:NV)=0.5d0*(DF(1:NV)*B%sj(i,j)+w_LU*B%Lcj(i,j-1)*B%dU(:,i,j-1))
       if(If_viscous .eq. 1)    duj(1:NV)=duj(1:NV)+B%Lvj(i,j-1)*B%dU(:,i,j-1)    ! 2012-2-29

	 else
	   duj=0.d0
	 endif

	do m=1,NV
      B%dU(m,i,j)=(B%Res(m,i,j)+dui(m)+duj(m))/alfa(m)
	enddo

   enddo
! $OMP END PARALLEL DO 
   enddo
!----------------------------------------------------------
!  �� (nx-1,ny-1)��(1,1)��ɨ����� ������ɨ����̣�
!  plane=i+j+k
   do plane=nx+ny-2,2,-1   

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(plane,nx,ny,NV,B,gamma,If_viscous,alfa1)
	 do j=ny-1,1,-1
	 i=plane-j
	 if( i .lt. 1 .or. i .gt. nx-1) cycle            ! ���������ƽ��
      alfa(:)=B%vol(i,j)/B%dt(i,j)+w_LU*(B%Lci(i,j)+B%Lcj(i,j)) +B%vol(i,j)*alfa1
      if(If_viscous .eq. 1)  alfa(:)=alfa(:)+2.d0*(B%Lvi(i,j)+B%Lvj(i,j) )         
       if(NV .eq. 5) then
         alfa(5)=alfa(5)+B%Lvi(i,j)+B%Lvj(i,j)
!       elseif(NV .eq. 6) then
!  	      alfa(5)=alfa(5)+0.09*B%U(6,i,j)/B%U(1,i,j)
!         alfa(6)=alfa(6)+2.d0*0.0828*B%U(6,i,j)/B%U(1,i,j)
	   endif

	 if(i.ne. nx-1) then
!                              ͨ���Ĳ������������Ƽ���A*W (See Blazek's book, page 208)
       call comput_DFn(NV,DF(1:NV),B%U(:,i+1,j),B%DU(:,i+1,j),B%ni1(i,j),B%ni2(i,j),gamma)  
       dui(1:NV)=-0.5d0*(DF(1:NV)*B%si(i+1,j)-w_LU*B%Lci(i+1,j)*B%dU(:,i+1,j))
       if(If_viscous .eq. 1)    dui(1:NV)=dui(1:NV)+B%Lvi(i+1,j)*B%dU(:,i+1,j)
	   else
	   dui=0.d0
       endif
	 
	 if(j.ne. ny-1) then
       call comput_DFn(NV,DF(1:NV),B%U(:,i,j+1),B%DU(:,i,j+1),B%nj1(i,j),B%nj2(i,j),gamma)  
       duj(1:NV)=-0.5d0*(DF(1:NV)*B%sj(i,j+1)-w_LU*B%Lcj(i,j+1)*B%dU(:,i,j+1))
       if(If_viscous .eq. 1)    duj(1:NV)=duj(1:NV)+B%Lvj(i,j+1)*B%dU(:,i,j+1)
	 else
	   duj=0.d0
	 endif

	do m=1,NV
      B%dU(m,i,j)=B%dU(m,i,j)+(dui(m)+duj(m))/alfa(m)
	enddo
   enddo
! $OMP END PARALLEL DO 
   enddo
   end subroutine  du_LU_SGS_2D







!  ����ͨ���Ĳ��� DF=F(U+DU)-F(U),  LU-SGS������ʹ�ã���������A*DU  
    subroutine comput_DFn(NV,DF,U,DU,n1,n2,gamma)
	implicit none
    integer:: NV
	real*8,dimension(NV):: DF,U,DU,U2
	real*8:: n1,n2,un1,un2,gamma,p1,p2
    U2=U+DU                                 ! �µ��غ����
	un1=(U(2)*n1+U(3)*n2)/U(1)              !un �����ٶ�
	p1=(gamma-1.d0)*(U(4)-0.5d0*(U(2)**2+U(3)**2)/U(1))   ! ѹ��
	un2=(U2(2)*n1+U2(3)*n2)/U2(1)          !�����ٶ�un
	p2=(gamma-1.d0)*(U2(4)-0.5d0*(U2(2)**2+U2(3)**2)/U2(1))    !ѹ��
!  ͨ��֮�� DF=F(U+DU)-F(U)
    DF(1)=U2(1)*un2-U(1)*un1                  ! d*un
	DF(2)=(U2(2)*un2+p2*n1)-(U(2)*un1+p1*n1)
	DF(3)=(U2(3)*un2+p2*n2)-(U(3)*un1+p1*n2)
	DF(4)=(U2(4)+p2)*un2-(U(4)+p1)*un1
     if(NV .eq. 5) then
      DF(5)=U2(5)*un2-U(5)*un1
	 elseif (NV .eq. 6) then       
	  DF(5)=U2(5)*un2-U(5)*un1   ! k���̵Ķ���ͨ��
	  DF(6)=U2(6)*un2-U(6)*un1   ! w���̵Ķ���ͨ��
	 endif
	end subroutine comput_dFn