!----OpenCFD-EC 2d--------------------------------------------------------------
!  Copyright by Li Xinliang, Institute of Mechanics, lixl@imech.ac.cn                 
!  read flow and plot the wall pressure coefficient
!  Ver 1.2 , read opencfd-ec ver 1.4.3 
!==============================================
  module Goem_Var
  implicit none
 ! global parameter (for all blocks) 
   integer:: Num_Block
   real*8,parameter:: gamma=1.4d0
   TYPE Block_TYPE           ! local variables (for each block) 
     integer :: Block_no,nx,ny,subface
     real*8,pointer,dimension(:,:):: x,y,x1,y1  ! (x,y) : coordinate of vortex; (x1,y1): coordinate of cell center  
     real*8,pointer,dimension(:,:) :: d,u,v,T,p
     integer,pointer,dimension(:,:):: bc_msg    ! see bc2d.in
   End TYPE Block_TYPE  
  
   TYPE (Block_TYPE), save,dimension(:),allocatable,target:: Block
   real*8,save:: Re,Ma,P00
  end module Goem_Var
!------------------------------------------------------------------


  program main
  use Goem_Var
  implicit none
  Type (Block_TYPE),pointer:: B
  integer  mBlock,i,j,k,m,ibegin,iend,jbegin,jend,ksub

  open(99,file="control.in")
    read(99,*)
    read(99,*) Ma, Re
  close(99)
   p00=1.d0/(gamma*Ma*Ma)
   call init
  open(99,file="Cp.dat")
  write(99,*) "variable=x,Cp"
  do m=1, Num_Block 
   B=> Block(m)
   do  ksub=1,B%subface
   if(B%bc_msg(7,ksub) .eq. -10) then
   ibegin=B%bc_msg(3,ksub); iend=B%bc_msg(4,ksub); jbegin=B%bc_msg(5,ksub); jend=B%bc_msg(6,ksub)
   write(99,*) "zone i=", (jend-jbegin+1)*(iend-ibegin+1)
   do j=jbegin,jend
   do i=ibegin,iend
   write(99,*) B%x1(i,j),-(B%p(i,j)-p00)*2.d0
   enddo
   enddo
   endif
   enddo
  enddo

  end

   

 
!------------------------------------------------------------------------------     
   subroutine init
   use Goem_Var
   implicit none
   integer :: i,j,k,m,nx1,ny1,Num_Block1
   integer,allocatable,dimension(:):: NI,NJ
   Type (Block_TYPE),pointer:: B
   integer,parameter::LAP=3
 !--------------------------------------------------------
! ---------node Coordinates----------------------------------------  
   print*, "read Mesh2d.dat"
   open(99,file="Mesh2d.dat")
    read(99,*) Num_Block
    allocate(Block(Num_Block))
    allocate( NI(Num_Block),NJ(Num_Block) )
    read(99,*) (NI(k), NJ(k), k=1,Num_Block)
    do m=1,Num_Block
     B => Block(m)
     B%Block_no=m
     B%nx=NI(m); B%ny=NJ(m) 
     nx1=B%nx ; ny1= B%ny
     allocate(B%x(0:nx1+1,0:ny1+1), B%y(0:nx1+1,0:ny1+1))
     allocate(B%x1(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1),B%y1(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1))
     allocate(B%d(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1),B%u(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1),B%v(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1),B%T(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1),B%p(1-LAP:nx1+LAP-1,1-LAP:ny1+LAP-1))
     read(99,*) ((B%x(i,j),i=1,nx1),j=1,ny1),  ((B%y(i,j),i=1,nx1),j=1,ny1)
     enddo
     close(99) 
!----Mesh control message (bc2d.in)------------------------------------------
    print*, "read bc2d.in"
    open(88,file="bc2d_1.in")
    read(88,*)
    read(88,*)
    read(88,*) Num_Block1
    if(Num_Block1 .ne. Num_Block) then
      print*, "Error!  Block number in bc2d.in is not equal to that in Mesh2d.dat !"
      stop
    endif
    do m=1,Num_Block
     B => Block(m)
     read(88,*)
     read(88,*)
     read(88,*) B%subface   !number of the subface in the Block m
     allocate(B%bc_msg(9,B%subface))   ! Boundary control message
     read(88,*)
     read(88,*) ((B%bc_msg(i,j),i=1,9),j=1,B%subface)
     enddo
     close(88)
!-------Initial of U ------------------------------------------------------------
        open(99,file="flow2d.dat")
         print*, "Init from 'flow2d.dat' "
         read(99,*)
        do m=1,Num_Block
          B => Block(m)
          read(99,*)
	      do j=0,B%ny
	      do i=0,B%nx
	        read(99,*) B%x1(i,j),B%y1(i,j),B%d(i,j),B%u(i,j),B%v(i,j),B%T(i,j),B%p(i,j)
          enddo
          enddo
        enddo
     close(99)

 end   




!---------Continue boundary (inner boundary) -------------------------------
     subroutine Update_coordinate_buffer
     use Goem_Var
     implicit none
     Type (Block_TYPE),pointer:: B,B1
     real*8,allocatable:: Ux(:,:)
     integer:: i,j,k,k1,m,mBlock,ksub,n1,m_neighbour,msub,orient,Kflag_initial,nx1,ny1
     integer:: ibegin1,iend1,jbegin1,jend1,ibegin2,iend2,jbegin2,jend2

 ! ------------------------------------------------------ 
 do mBlock=1,Num_Block
  B => Block(mBlock)
  do  ksub=1,B%subface
  ibegin2=B%bc_msg(3,ksub); iend2=B%bc_msg(4,ksub); jbegin2=B%bc_msg(5,ksub); jend2=B%bc_msg(6,ksub)       ! write
  if(B%bc_msg(7, ksub) .gt. 0 ) then   ! inner boundary
       n1=B%bc_msg(4,ksub)-B%bc_msg(3,ksub)+B%bc_msg(6,ksub)-B%bc_msg(5,ksub)+1  ! i_end-i_begin+j_end-j_begin+1
       allocate(Ux(2,n1))
       m_neighbour=B%bc_msg(7,ksub);  msub= B%bc_msg(8,ksub)   
       B1 => Block(m_neighbour)
       ibegin1=B1%bc_msg(3,msub); iend1=B1%bc_msg(4,msub); jbegin1=B1%bc_msg(5,msub); jend1=B1%bc_msg(6,msub)   ! read
 
       if(B1%bc_msg(2,msub) .eq. 1 ) then                     !  boundary  i-
          Ux(1,:)=B1%x(ibegin1+1,jbegin1:jend1) ;  Ux(2,:)=B1%y(ibegin1+1,jbegin1:jend1); 
       else if ( B1%bc_msg(2,msub) .eq. 2) then               !  boundary  j-
          Ux(1,:)=B1%x(ibegin1:iend1,jbegin1+1);   Ux(2,:)=B1%y(ibegin1:iend1,jbegin1+1)
       else if( B1%bc_msg(2,msub) .eq. 3) then                !  boundary  i+
          Ux(1,:)=B1%x(iend1-1,jbegin1:jend1) ;  Ux(2,:)=B1%y(iend1-1,jbegin1:jend1); 
       else                                                   !  boundary  j+
          Ux(1,:)=B1%x(ibegin1:iend1,jend1-1);   Ux(2,:)=B1%y(ibegin1:iend1,jend1-1)
       endif

       orient=B%bc_msg(9,ksub)
       do k=1,n1
       if(orient .eq. 2) then
         k1=k
       else
         k1=n1+1-k
       endif
       if(B%bc_msg(2,ksub) .eq. 1 ) then             !  boundary  i-
          B%x(ibegin2-1,jbegin2+k-1)=Ux(1,k1) ; B%y(ibegin2-1,jbegin2+k-1)=Ux(2,k1)
       else if(B%bc_msg(2,ksub) .eq. 2 ) then        !  boundary  j-
          B%x(ibegin2+k-1,jbegin2-1)=Ux(1,k1) ; B%y(ibegin2+k-1,jbegin2-1)=Ux(2,k1)
       else if (B%bc_msg(2,ksub) .eq. 3 ) then       !  boundary  i+
          B%x(iend2+1,jbegin2+k-1)=Ux(1,k1) ; B%y(iend2+1,jbegin2+k-1)=Ux(2,k1)
       else                                          !  boundary  j+
          B%x(ibegin2+k-1,jend2+1)=Ux(1,k1) ; B%y(ibegin2+k-1,jend2+1)=Ux(2,k1)
       endif
      enddo
      deallocate(Ux)

    else  ! not inner boundary
       if(B%bc_msg(2,ksub) .eq. 1 ) then             !  boundary  i-
          B%x(ibegin2-1,jbegin2:jend2)=2.d0*B%x(ibegin2,jbegin2:jend2)-B%x(ibegin2+1,jbegin2:jend2) 
          B%y(ibegin2-1,jbegin2:jend2)=2.d0*B%y(ibegin2,jbegin2:jend2)-B%y(ibegin2+1,jbegin2:jend2)   
       else if(B%bc_msg(2,ksub) .eq. 2 ) then        !  boundary  j-
          B%x(ibegin2:iend2,jbegin2-1)=2.d0*B%x(ibegin2:iend2,jbegin2)-B%x(ibegin2:iend2,jbegin2+1) 
          B%y(ibegin2:iend2,jbegin2-1)=2.d0*B%y(ibegin2:iend2,jbegin2)-B%y(ibegin2:iend2,jbegin2+1)  
       else if (B%bc_msg(2,ksub) .eq. 3 ) then       !  boundary  i+
          B%x(iend2+1,jbegin2:jend2)=2.d0*B%x(iend2,jbegin2:jend2)-B%x(iend2-1,jbegin2:jend2) 
          B%y(iend2+1,jbegin2:jend2)=2.d0*B%y(iend2,jbegin2:jend2)-B%y(iend2-1,jbegin2:jend2) 
       else                                          !  boundary  j+
          B%x(ibegin2:iend2,jend2+1)=2.d0*B%x(ibegin2:iend2,jend2)-B%x(ibegin2:iend2,jend2-1)
          B%y(ibegin2:iend2,jend2+1)=2.d0*B%y(ibegin2:iend2,jend2)-B%y(ibegin2:iend2,jend2-1)
       endif
     endif
    enddo
  !  Conner point     
       nx1=B%nx; ny1=B%ny
       B%x(0,0)=B%x(1,0)+B%x(0,1)-B%x(1,1);   B%y(0,0)=B%y(1,0)+B%y(0,1)-B%y(1,1) 
       B%x(0,ny1+1)=B%x(1,ny1+1)+B%x(0,ny1)-B%x(1,ny1);  B%y(0,ny1+1)=B%y(1,ny1+1)+B%y(0,ny1)-B%y(1,ny1) 
       B%x(nx1+1,0)=B%x(nx1,0)+B%x(nx1+1,1)-B%x(nx1,1);  B%y(nx1+1,0)=B%y(nx1,0)+B%y(nx1+1,1)-B%y(nx1,1) 
       B%x(nx1+1,ny1+1)=B%x(nx1,ny1+1)+B%x(nx1+1,ny1)-B%x(nx1,ny1); B%y(nx1+1,ny1+1)=B%y(nx1,ny1+1)+B%y(nx1+1,ny1)-B%y(nx1,ny1) 
  
  ! Coordinates of the Cell center 
     do j=0,ny1
     do i=0,nx1
	  B%x1(i,j)=(B%x(i,j)+B%x(i+1,j)+B%x(i,j+1)+B%x(i+1,j+1))*0.25
	  B%y1(i,j)=(B%y(i,j)+B%y(i+1,j)+B%y(i,j+1)+B%y(i+1,j+1))*0.25
     enddo
	 enddo

   enddo
 
   end
