!----通量分裂-----------------------------------------------
! 包含：Steger-Warming, Van Leer, Roe, AUSM, HLL, HLLC  (用户可选择其中一种)
! 程序编制： 李新亮，冷岩
!c========================================================
! Van Leer 流通矢量分裂 
!-------Extend 1D Van Leer FLux for finite volume method
! F=F+ (uL) + F- (uR)	

!-------------------------------------------------------------

!  Flux by using  HLL or  HLLC approximate Riemann slover
!  Copyright by Li Xinliang, Institute of Mechanics, CAS; lixl@imech.ac.cn
!  Revised by Leng Yan (In HLLC, the second kind method in Toro's book is used; the first method is unstable)
!  Ref. Toro : Riemann Solvers and Numerical Methods for Fluid Dynamics,  p330       
!  Iflag_HLL_HLLC==0 : HLL  ;   ==1 : HLLC       
       subroutine Flux_HLL_HLLC_1D(QL,QR,Flux,gamma,Iflag_HLL_HLLC)   
       use const_var
       implicit none
       integer:: Iflag_HLL_HLLC,IFlag_Reconstruction
       real*8:: QL(4),QR(4),UL(4),UR(4),Flux(4),gamma
       real*8:: dl,uul,vvl,pl,al, dr,uur,vvr,pr,ar    ! uu velocity
       real*8:: p_pvrs,p_star,qql,qqr,Sl,Sr,Fl(4),Fr(4),S_star,tmpl,tmpr,tmp,F_starl(4),F_starr(4)
       
        dl=QL(1); uul=QL(2); vvl=QL(3); pl=QL(4)
        dr=QR(1); uur=QR(2); vvr=QR(3); pr=QR(4)
        UL(1)=dl; UL(2)=dl*uul; UL(3)=dl*vvl; UL(4)=pl/(gamma-1.d0)+dl*(uul*uul+vvl*vvl)*0.5d0 
        UR(1)=dr; UR(2)=dr*uur; UR(3)=dr*vvr; UR(4)=pr/(gamma-1.d0)+dr*(uur*uur+vvr*vvr)*0.5d0 
    	al=sqrt(gamma*pl/dl); ar=sqrt(gamma*pr/dr)
 	    p_pvrs=0.5d0*(pl+pr)-(uur-uul)*(dl+dr)*(al+ar)*0.125d0
        p_star=max(0.d0,p_pvrs)  
    
       if(p_star .le. pl) then
         qql=1
       else
         qql=sqrt(1.d0+ ((gamma+1.d0)/(2.d0*gamma)) * (p_star/pl-1.d0) )
       endif
       if(p_star .le. pr) then
         qqr=1
       else
         qqr=sqrt(1.d0+ ((gamma+1.d0)/(2.d0*gamma)) * (p_star/pr-1.d0) )
       endif
       Sl=uul-al*qql; Sr=uur+ar*qqr    ! speed of the lift and the right shockwaves
       
       Fl(1)=UL(2); Fl(2)=UL(2)*uul+pl; Fl(3)=UL(3)*uul; Fl(4)=uul*(UL(4)+pl)
       Fr(1)=UR(2); Fr(2)=UR(2)*uur+pr; Fr(3)=UR(3)*uur; Fr(4)=uur*(UR(4)+pr)
!----HLL---------------------------------
     
    if(Iflag_HLL_HLLC .eq. 0) then  ! HLL Flux
       if( Sl .ge. 0 ) then
         Flux=Fl
       else if (Sr .le. 0) then
         Flux=Fr
       else
         Flux=(Sr*Fl-Sl*Fr+Sl*Sr*(UR-UL))/(Sr-Sl)
       endif
    else  
!    HLLC Flux
     S_star=(pr-pl+dl*uul*(Sl-uul)-dr*uur*(Sr-uur))/(dl*(Sl-uul)-dr*(Sr-uur))
     tmpl=dl*(Sl-uul)/(Sl-S_star) ; tmpr=dr*(Sr-uur)/(Sr-S_star)
      F_starl(1)=Fl(1)+Sl*(tmpl-UL(1))
      F_starl(2)=Fl(2)+Sl*(tmpl*S_star-UL(2))
      F_starl(3)=Fl(3)+Sl*(tmpl*vvl-UL(3))
      tmp=UL(4)/dl+(S_star-uul)*(S_star+pl/(dl*(Sl-uul)))
      F_starl(4)=Fl(4)+Sl*(tmpl*tmp-UL(4))
      F_starr(1)=Fr(1)+Sr*(tmpr-UR(1))
      F_starr(2)=Fr(2)+Sr*(tmpr*S_star-UR(2))
      F_starr(3)=Fr(3)+Sr*(tmpr*vvr-UR(3))
      tmp=UR(4)/dr+(S_star-uur)*(S_star+pr/(dr*(Sr-uur)))
      F_starr(4)=Fr(4)+Sr*(tmpr*tmp-UR(4))

! Revised by Leng Yan  
!     S_star=(pr-pl+dl*uul*(Sl-uul)-dr*uur*(Sr-uur))/(dl*(Sl-uul)-dr*(Sr-uur))
!     F_starl(1)=S_star*(Sl*UL(1)-Fl(1))/(Sl-S_star)
!     F_starl(2)=(S_star*(Sl*UL(2)-Fl(2))+Sl*(pl+dl*(Sl-uul)*(S_star-uul)))/(Sl-S_star)
!     F_starl(3)=S_star*(Sl*UL(3)-Fl(3))/(Sl-S_star)
!     F_starl(4)=(S_star*(Sl*UL(4)-Fl(4))+Sl*(pl+dl*(Sl-uul)*(S_star-uul))*S_star)/(Sl-S_star)
!     F_starr(1)=S_star*(Sr*UR(1)-Fr(1))/(Sr-S_star)
!     F_starr(2)=(S_star*(Sr*UR(2)-Fr(2))+Sr*(pr+dr*(Sr-uur)*(S_star-uur)))/(Sr-S_star)
!     F_starr(3)=S_star*(Sr*UR(3)-Fr(3))/(Sr-S_star)
!     F_starr(4)=(S_star*(Sr*UR(4)-Fr(4))+Sr*(pr+dr*(Sr-uur)*(S_star-uur))*S_star)/(Sr-S_star)
 
  
      if( Sl .ge. 0 ) then
        Flux=Fl
      else if (Sr .le. 0) then
        Flux=Fr
      else if (S_star .ge. 0) then
        Flux=F_starl
      else
        Flux=F_starr
      endif
   endif
   end


!c========================================================
 !-------Extend 1D steger_warming FLux for finite volume method
! F=F+ (uL) + F- (uR)	
     subroutine Flux_steger_warming_1Da(QL,QR,Flux,gamma)   
     use const_var
     implicit none
    integer:: IFlag_Reconstruction
     real*8:: QL(4),QR(4),Flux(4),gamma
     real*8:: dl,uul,vvl,pl,al, dr,uur,vvr,pr,ar ,pr1   ! uu velocity
     real*8:: tmp0,tmp1,tmp2,tmp3,E1P,E2P,E3P,E1M,E2M,E3M,fp(4),fm(4)
 
        dl=QL(1); uul=QL(2); vvl=QL(3); pl=QL(4) 
        dr=QR(1); uur=QR(2); vvr=QR(3); pr=QR(4)
	
	    al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
        ar=sqrt(gamma*pr/dr)  ! find a bug, removed
      
	    tmp1=2.d0*(gamma-1.d0)
        tmp3=(3.d0-gamma)/(2.d0*(gamma-1.d0)) 
! eigenvalues---------	      
        E1P=(uul+abs(uul))*0.5d0
        E2P=(uul-al+abs(uul-al))*0.5d0
        E3P=(uul+al+abs(uul+al))*0.5d0
        tmp0=dl/(2.d0*gamma) 
        fp(1)=tmp0*(tmp1*E1P+E2P+E3P)
        fp(2)=tmp0*(tmp1*E1P*uul+E2P*(uul-al)+E3P*(uul+al))
        fp(3)=tmp0*(tmp1*E1P*vvl+E2P*vvl+E3P*vvl)
        fp(4)=tmp0*(E1P*(gamma-1.d0)*(uul*uul+vvl*vvl)+E2p*((uul-al)**2+vvl*vvl)*0.5d0    &
	             +E3P*((uul+al)**2+vvl*vvl)*0.5d0+tmp3*al*al*(E2P+E3P))
     
        E1M=(uur-abs(uur))*0.5d0
        E2M=(uur-ar-abs(uur-ar))*0.5d0
        E3M=(uur+ar-abs(uur+ar))*0.5d0
        tmp0=dr/(2.d0*gamma) 
 	    fm(1)=tmp0*(tmp1*E1M+E2M+E3M)
        fm(2)=tmp0*(tmp1*E1M*uur+E2M*(uur-ar)+E3M*(uur+ar))
        fm(3)=tmp0*(tmp1*E1M*vvr+E2M*vvr+E3M*vvr)
        fm(4)=tmp0*(E1M*(gamma-1.d0)*(uur*uur+vvr*vvr)+E2M*((uur-ar)**2+vvr*vvr)*0.5d0    &
              +E3M*((uur+ar)**2+vvr*vvr)*0.5d0+tmp3*ar*ar*(E2M+E3M))
        Flux=fp+fm
      end
!-----------------------------------------------------------
!  Code by Leng Yan
     subroutine Flux_Van_Leer_1Da(QL,QR,Flux,gamma)   
     use Const_Var
     implicit none
     real*8:: QL(4),QR(4),Flux(4),gamma
     real*8:: dl,uul,vvl,pl,al, dr,uur,vvr,pr,ar ,Ml,Mr,Mp,Mm  ! uu velocity
     real*8:: tmp0,fp(4),fm(4)

        dl=QL(1); uul=QL(2); vvl=QL(3);  pl=QL(4) 
        dr=QR(1); uur=QR(2); vvr=QR(3);  pr=QR(4)
        al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
        ar=sqrt(gamma*pr/dr)  
        Ml=uul/(al); Mr=uur/(ar)
        if(Ml>=1.d0) then
         fp(1)=dl*uul
         fp(2)=dl*uul*uul+pl
         fp(3)=dl*uul*vvl
         fp(4)=uul*(gamma*pl/(gamma-1.d0)+0.5d0*dl*(uul*uul+vvl*vvl))
        else if(abs(Ml)<1.d0) then 
         Mp=0.25d0*(1.d0+Ml)*(1.d0+Ml)
         tmp0=dl*al*Mp
         fp(1)=tmp0
         fp(2)=tmp0*((gamma-1.d0)*uul+2.d0*al)/gamma
         fp(3)=tmp0*vvl
         fp(4)=tmp0*(((gamma-1.d0)*uul+2.d0*al)*((gamma-1.d0)*uul+2.d0*al)*0.5d0/(gamma*gamma-1.d0)+0.5d0*(vvl*vvl))
        else if(Ml<=-1.d0)   then
         fp(1)=0.d0
         fp(2)=0.d0
         fp(3)=0.d0
         fp(4)=0.d0
        end if	
		
        if(Mr>= 1.d0) then
         fm(1)=0.d0
         fm(2)=0.d0
         fm(3)=0.d0
         fm(4)=0.d0
        else if(abs(Mr) < 1.d0) then 
         Mm=-0.25d0*(Mr-1.d0)*(Mr-1.d0)
         tmp0=dr*ar*Mm
         fm(1)=tmp0
         fm(2)=tmp0*((gamma-1.d0)*uur-2.d0*ar)/gamma
         fm(3)=tmp0*vvr
         fm(4)=tmp0*(((gamma-1.d0)*uur-2.d0*ar)*((gamma-1.d0)*uur-2.d0*ar)*0.5d0/(gamma*gamma-1.d0)+0.5d0*(vvr*vvr))   
        else if(Mr<=-1.d0)    then
         fm(1)=dr*uur
         fm(2)=dr*uur*uur+pr
         fm(3)=dr*uur*vvr
         fm(4)=uur*(gamma*pr/(gamma-1.d0)+0.5d0*dr*(uur*uur+vvr*vvr))
        end if
	  Flux=fp+fm
      end
!------------------------------------------------------------------------------------------------
! AUSM+ 通量分裂 (by Leng Yan)
!-------Extend 1D Ausm FLux for finite volume method
! F=F+ (uL) + F- (uR)	

     subroutine Flux_Ausm_1Da(QL,QR,Flux,gamma)   
     use Const_Var
     implicit none
     real*8:: QL(4),QR(4),Flux(4),gamma
     real*8:: dl,uul,vvl,pl,al, hl,dr,uur,vvr,pr,ar,hr  ! uu velocity
     real*8:: a,Ml,Mr,Mp4,Mm4,Pp5,Pm5,M,mp,mm,p,dm,xm
     real*8:: fp(4),fm(4)

 
        dl=QL(1); uul=QL(2); vvl=QL(3);  pl=QL(4)
        dr=QR(1); uur=QR(2); vvr=QR(3);  pr=QR(4)
        hl=gamma*pl/((gamma-1.d0)*dl)+0.5*(uul*uul+vvl*vvl)
        hr=gamma*pr/((gamma-1.d0)*dr)+0.5*(uur*uur+vvr*vvr)
	    
    	al=sqrt(gamma*pl/dl)  ! density, velocity, pressure and sound speed
        ar=sqrt(gamma*pr/dr) 
        a=0.5d0*(al+ar)
        Ml=uul/a; Mr=uur/a
 
        if(abs(Ml)>=1.d0) then
          Mp4=0.5d0*(Ml+abs(Ml))
        else
          Mp4=0.25*(Ml+1.d0)*(Ml+1.d0)+0.125d0*(Ml*Ml-1.d0)*(Ml*Ml-1.d0)
        end if
        if(abs(Mr)>=1.d0) then
           Mm4=0.5d0*(Mr-abs(Mr))
         else
           Mm4=-0.25*(Mr-1.d0)*(Mr-1.d0)-0.125d0*(Mr*Mr-1.d0)*(Mr*Mr-1.d0)
         end if
        
        M=Mp4+Mm4
        mp=dl*a*max(0.d0,M)
        mm=dr*a*min(0.d0,M)
        xm=mp+mm
        if(M>=0.d0) then
         dm=a*abs(M)*dl
        else
         dm=a*abs(M)*dr
        end if

        if(abs(Ml)>=1.d0) then
         Pp5=0.5d0*(Ml+abs(Ml))/Ml
        else
         Pp5=0.25*(Ml+1.d0)*(Ml+1.d0)*(2.d0-Ml)+3.d0*Ml*(Ml*Ml-1.d0)*(Ml*Ml-1.d0)/16.d0
        end if
        if(abs(Mr)>=1.d0) then
         Pm5=0.5d0*(Mr-abs(Mr))/Mr
        else
         Pm5=0.25*(Mr-1.d0)*(Mr-1.d0)*(2.d0+Mr)-3.d0*Mr*(Mr*Mr-1.d0)*(Mr*Mr-1.d0)/16.d0
        end if
         p=pl*Pp5+pr*Pm5
    
        fp(1)=1.d0
        fp(2)=uul
        fp(3)=vvl
        fp(4)=hl
        
        fm(1)=1.d0
        fm(2)=uur
        fm(3)=vvr
        fm(4)=hr
        
        Flux(1)=0.5d0*xm*(fp(1)+fm(1))-0.5d0*dm*(fm(1)-fp(1))
        Flux(2)=0.5d0*xm*(fp(2)+fm(2))-0.5d0*dm*(fm(2)-fp(2))+p
        Flux(3)=0.5d0*xm*(fp(3)+fm(3))-0.5d0*dm*(fm(3)-fp(3))
        Flux(4)=0.5d0*xm*(fp(4)+fm(4))-0.5d0*dm*(fm(4)-fp(4))
      end
!-----------------------------------------------------------
!c==================================================================================
! Roe 型FVS分裂  (Code by LengYan, Revised by Li Xinliang)
   subroutine Flux_Roe_1D(QL,QR,Flux,gamma)   
   implicit none
   integer:: flag
   real*8:: QL(4),QR(4),Flux(4),gamma
   real*8:: dl,uul,vvl,pl,al, dr,uur,vvr,pr,ar    ! uu velocity
   real*8:: Fl(4),Fr(4),  lamda(4),arr(4), avd,avu,avv,avh,ava, hl,hr,tmpp
   real*8,parameter:: delt=0.1d0   !Harten熵修正的系数 （0至0.125，系数越大粘性越大)

  
   dl=QL(1); uul=QL(2); vvl=QL(3); pl=QL(4)
   dr=QR(1); uur=QR(2); vvr=QR(3); pr=QR(4)

   hl=gamma*pl/((gamma-1.d0)*dl)+0.5d0*(uul*uul+vvl*vvl)
   hr=gamma*pr/((gamma-1.d0)*dr)+0.5d0*(uur*uur+vvr*vvr)  ! find a bug "dl" should be "dr"

! 通量
   Fl(1)=dl*uul; Fl(2)=dl*uul*uul+pl; Fl(3)=dl*uul*vvl; Fl(4)=uul*dl*hl
   Fr(1)=dr*uur; Fr(2)=dr*uur*uur+pr; Fr(3)=dr*uur*vvr; Fr(4)=uur*dr*hr
!---Roe 平均 ------------------  
   tmpp=sqrt(dr/dl)
   avd=0.25d0*dl*(1.d0+tmpp)*(1.d0+tmpp)  !!! ??? avd=dl*tmpp
   avu=(uul+tmpp*uur)/(1.d0+tmpp)
   avv=(vvl+tmpp*vvr)/(1.d0+tmpp)
   avh=(hl+tmpp*hr)/(1.d0+tmpp)                         ! 平均总焓
!   ava=sqrt((gamma-1.d0)*(avh-0.5d0*(avu*avu+avv+avv)))   ! 平均声速  find a bug  (in Ver 1.0) !!!
    ava=sqrt((gamma-1.d0)*(avh-0.5d0*(avu*avu+avv*avv)))    ! 平均声速

    
!---Harten 型熵修正---------
!  A bug is reomved, 2011-5-5
       lamda(1)=abs(avu-ava)
       lamda(2)=abs(avu)
       lamda(3)=abs(avu)
       lamda(4)=abs(avu+ava)


   if(lamda(1) < delt) then
     lamda(1)=(lamda(1)**2 +delt*delt)/(2.d0*delt)
   end if
   if(lamda(2) < delt) then
     lamda(2)=(lamda(2)**2 +delt*delt)/(2.d0*delt)
   end if
   if(lamda(3) < delt) then
     lamda(3)=(lamda(3)**2 +delt*delt)/(2.d0*delt)
   end if
   if(lamda(4) < delt) then
     lamda(4)=(lamda(4)**2 +delt*delt)/(2.d0*delt)
   end if

!---------------------------
   arr(1)=(avd*ava*(-uur+uul)+pr-pl)/(2.d0*ava*ava)
   arr(2)=dr-dl-(pr-pl)/(ava*ava)
   arr(3)=avd*(vvr-vvl)/ava
   arr(4)=(avd*ava*(uur-uul)+pr-pl)/(2.d0*ava*ava)
   flux(1)=0.5d0*(fl(1)+fr(1))-0.5d0*(lamda(1)*arr(1)+lamda(2)*arr(2)+lamda(4)*arr(4))
   flux(2)=0.5d0*(fl(2)+fr(2))-0.5d0*(lamda(1)*arr(1)*(avu-ava)+lamda(2)*arr(2)*avu+lamda(4)*arr(4)*(avu+ava))
   flux(3)=0.5d0*(fl(3)+fr(3))-0.5d0*(lamda(1)*arr(1)*avv+lamda(2)*arr(2)*avv+lamda(3)*arr(3)+lamda(4)*arr(4)*avv)
   flux(4)=0.5d0*(fl(4)+fr(4))-0.5d0*(lamda(1)*arr(1)*(avh-avu*ava) &
           +lamda(2)*arr(2)*0.5d0*(avu*avu+avv*avv)+lamda(3)*arr(3)*avv+lamda(4)*arr(4)*(avh+avu*ava))
   end
!--------------------------------------------------------------------------------------------