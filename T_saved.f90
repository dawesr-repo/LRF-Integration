
    !********************************************************
module Tensors_constant  
   implicit none 
   integer,parameter::L=13
   real*8 :: T_Tensor(2*L+1,2*L+1,2*L+1,2*L+1),T_Changed(2*L+1,2*L+1,2*L+1,2*L+1)

  contains

  subroutine init_Tensors()  
     implicit none 

     T_Tensor  = 0d0
     T_Changed = 0d0


  end subroutine init_Tensors

  subroutine change_Tensors(la_,ka1_,ka2_,lb_,kb1_,kb2_,val)  
    integer ,INTENT(IN)::la_,ka1_,lb_,kb1_
      Character(len = 1) , INTENT(IN)::ka2_,kb2_
      real*8 ,INTENT(In)::val
      integer :: ka,kb,la,ka1,lb,kb1
      Character(len = 1) ::ka2,kb2
      la = la_
      ka1 = ka1_
      ka2 = ka2_
      lb = lb_
      kb1 = kb1_
      kb2 = kb2_ 

      IF (ka1<0 .or. kb1<0 .or. la<0 .or. lb<0 .or. ka1>la .or. kb1>lb) THEN
        RETURN
      ELSE IF (la==0 .and. lb==0) THEN  
        RETURN
      ELSE
        call T_component(la,ka1,ka2,ka) 
        call T_component(lb,kb1,kb2,kb)  

        T_Tensor(la+1,ka,lb+1,kb)=val
        T_Changed(la+1,ka,lb+1,kb)=1d0

        Return

    endif
  end subroutine change_Tensors
  
  subroutine get_Tensor(la_,ka1_,ka2_,lb_,kb1_,kb2_,T,chang)  
      
      integer ,INTENT(IN)::la_,ka1_,lb_,kb1_
      Character(len = 1) , INTENT(IN)::ka2_,kb2_
      real*8 ,INTENT(OUT)::T,chang
      integer :: ka,kb,la,ka1,lb,kb1
      Character(len = 1) ::ka2,kb2
      la = la_
      ka1 = ka1_
      ka2 = ka2_
      lb = lb_
      kb1 = kb1_
      kb2 = kb2_ 

      IF (ka1<0 .or. kb1<0 .or. la<0 .or. lb<0 .or. ka1>la .or. kb1>lb) THEN
        T=0d0
        chang =0d0
        RETURN
      ELSE IF (la==0 .and. lb==0) THEN  
        T = 1d0
        chang=1d0
        RETURN
      ELSE
        call T_component(la,ka1,ka2,ka) 
        call T_component(lb,kb1,kb2,kb)  
        if (ka>0 .and. kb>0)Then
          T     =  T_Tensor (la+1,ka,lb+1,kb)
          chang =  T_Changed(la+1,ka,lb+1,kb)
          RETURN
        else
          T=0d0
          chang =0d0
          RETURN

        end if
        
       endif
  end subroutine  get_Tensor


  subroutine T_component(l_,k1_,k2_,k)  
      
      integer ,INTENT(IN)::k1_,l_
      Character(len = 1) , INTENT(IN)::k2_
      integer , INTENT(OUT)::k
      integer :: k1,l
      Character(len = 1) ::k2

      k1 = k1_
      k2 = k2_
      l = l_

      IF (k1<0  .or. l<0 .or. k1>l ) THEN
         k = 0;
      ELSE IF (l==0 .and. k1==0) THEN  
        if (k2=="0" )then
          k = 1;
        else
          k = 0;
        end if
      ELSE
        if (k2=="s")then
          k=2*k1+1
        elseif (k2=="c")then
            k=2*k1
        else
           k = 1
        endif
      
    endif
      
  end subroutine T_component
      
end module Tensors_constant
  
