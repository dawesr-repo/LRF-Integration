
    !********************************************************
module Geometry_Constant  
   implicit none
   
    integer,parameter::Lev = 15
   public
    real*8 :: T_Tensor(2*Lev+1,2*Lev+1,2*Lev+1,2*Lev+1),T_Changed(2*Lev+1,2*Lev+1,2*Lev+1,2*Lev+1)
    real*8 , dimension(3):: Ar 
    real*8 , dimension(3):: Br
    real*8 , dimension(9):: CC
    real*8 , dimension(11):: cal_coord
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

  RECURSIVE SUBROUTINE T_lk(la,ka1,ka2,lb,kb1,kb2,Tlk)
    
    IMPLICIT NONE
    

    real*8 , dimension(9):: Cab
    Integer , INTENT(IN):: la,ka1,lb,kb1
    Character(len = 1) , INTENT(IN)::ka2,kb2
    Character(len = 20), dimension(3):: coordinates = ["z","y","x"]
    Integer :: i ,j   ,k_res1,ka_res1,kb_res1
    real*8 , INTENT(INOUT):: Tlk
    real*8 :: M,M1,M2,NN_,Tlk_new,NN_1,NN_2,Tlk_1,Tlk_2
    Character (len = 1)::k_res2,ka_res2,kb_res2


    real*8:: eps = 1.0D-8 ,T,chang




    Cab(1) = CC(1) !zz
    Cab(2) = CC(7) !zy
    Cab(3) = CC(4) !zx
    Cab(4) = CC(3) !yz
    Cab(5) = CC(9) !yy
    Cab(6) = CC(6) !yx
    Cab(7) = CC(2) !xz
    Cab(8) = CC(8) !xy
    Cab(9) = CC(5) !xx


    IF (ka1<0 .or. kb1<0 .or. la<0 .or. lb<0 .or. ka1>la .or. kb1>lb) THEN
        Tlk = 0d0
        RETURN
    ELSE IF (la==0 .and. lb==0) THEN  
        Tlk = 1d0
        RETURN
    ELSE

    
      call get_Tensor(la,ka1,ka2,lb,kb1,kb2,T,chang) 


      if (DABS(chang-1d0)<eps)Then
          Tlk =T
        RETURN
      else

        IF (lb==0) THEN

          Tlk = 0d0;
          
        
          do i = 1,3
            call MT(coordinates(i),ka1,ka2,M)
            call N_eta(coordinates(i),ka1,ka2,k_res1,k_res2)
            call NN(la-1,k_res1,0,0,NN_)

            if (abs(M)> eps .and. abs(2.d0*la-1d0) > eps .and. abs(Ar(i))>eps  ) Then
              call T_lk(la-1,k_res1,k_res2,0,0,"0",Tlk_new)
              Tlk = Tlk + ((2.d0*la-1d0)/(la*1.d0))*M*Ar(i)*Tlk_new*NN_
            end if
          end do
        
          call NN(la-2,ka1 ,0,0,NN_)

          if (abs(la-1d0) > eps ) Then
              call T_lk(la-2,ka1,ka2,0,0,"0",Tlk_new)
              Tlk = Tlk  - ((la-1d0)/(la*1d0))*Tlk_new*NN_;
          end if
          call NN(la,ka1 ,lb,kb1,NN_)
          Tlk = Tlk /NN_;
          call change_Tensors(la,ka1,ka2,lb,kb1,kb2,Tlk) 
          RETURN
      ELSE IF (la==0 ) THEN  

        Tlk = 0d0;
          
        
        do i = 1,3
          call MT(coordinates(i),kb1,kb2,M)
          call N_eta(coordinates(i),kb1,kb2,k_res1,k_res2)
          call NN(lb-1,k_res1,0,0,NN_)

          if (abs(M)> eps .and. abs(2.d0*lb-1d0) > eps .and. abs(Br(i))>eps ) Then

            call T_lk(0,0,"0",lb-1,k_res1,k_res2,Tlk_new)
            Tlk = Tlk + ((2.d0*lb-1d0)/(lb*1.d0))*M*Br(i)*Tlk_new*NN_

          end if 
        end do
      
        if (abs(1d0*lb-1d0) > eps ) Then 
          call NN(lb-2,kb1 ,0,0,NN_)
          call T_lk(0,0,"0",lb-2,kb1,kb2,Tlk_new)
          
          Tlk = Tlk  - ((1d0*lb-1d0)/(lb*1d0))*Tlk_new*NN_;
        end if 
        call NN(la,ka1 ,lb,kb1,NN_)
        Tlk = Tlk /NN_;
        call change_Tensors(la,ka1,ka2,lb,kb1,kb2,Tlk) 
        RETURN
      ELSE
        

        Tlk = 0d0;

        Call NN(la-2,ka1,lb,kb1,NN_1)
        Call NN(la,ka1,lb-2,kb1,NN_2)
        Call T_lk(la-2,ka1,ka2,lb,kb1,kb2,Tlk_1)
        

        Tlk  =Tlk  + Tlk_1*NN_1

        if (abs(2d0*la+1.d0*lb-1d0) > eps ) Then 

          Call T_lk(la,ka1,ka2,lb-2,kb1,kb2,Tlk_2)
          Tlk  =Tlk  - ((2d0*la+1.d0*lb-1d0)/lb*1d0)*Tlk_2*NN_2;
        end if

        do i = 1,3

            Call N_eta(coordinates(i),kb1,kb2,k_res1,k_res2);
            call MT(coordinates(i),kb1,kb2,M)
            Call NN(la,ka1,lb-1,k_res1,NN_);

            if (abs(M) > eps .and. abs(2d0*la+2d0*lb-1d0) > eps .and. abs(Br(i))>eps) Then

              Call T_lk(la,ka1,ka2,lb-1,k_res1,k_res2,Tlk_new)
            
              Tlk  =Tlk + ((2d0*la+2d0*lb-1d0)/(lb*1d0))*M*Br(i)*Tlk_new*NN_;

            end if 

          
          
        end do


        do i = 1,3
          do j = 1,3

                Call N_eta(coordinates(i),ka1,ka2,ka_res1,ka_res2);
                Call N_eta(coordinates(j),kb1,kb2,kb_res1,kb_res2);
                call MT(coordinates(i),ka1,ka2,M1)
                call MT(coordinates(j),kb1,kb2,M2)

                if (abs(M1) > eps .and. abs(M2) > eps .and.  abs(2d0*la-1d0) > eps .and. abs(Cab((i-1)*3+j))>eps) Then

                  Call NN(la-1,ka_res1,lb-1,kb_res1,NN_);
                
                  Call T_lk(la-1,ka_res1,ka_res2,lb-1,kb_res1,kb_res2,Tlk_new)

                  Tlk  =Tlk +  ((2d0*la-1d0)/(lb*1d0))*M1*M2*Cab((i-1)*3+j)*Tlk_new*NN_;
                end if 

                
            
          end do
        end do

        call NN(la,ka1 ,lb,kb1,NN_)
        Tlk = Tlk /NN_;
        call change_Tensors(la,ka1,ka2,lb,kb1,kb2,Tlk) 
        RETURN
      END IF 
      endif

    END IF 

  End   SUBROUTINE T_lk
  !********************************************************

  



  !********************************************************
  SUBROUTINE T_l0(q,M2,IND,lk , result)
          IMPLICIT NONE

          INTEGER ::  i
          real*8, INTENT(INOUT) ::result 
          Integer, INTENT(IN):: ind,lk
          real*8 :: q,t_k1,t_k2
          real*8 , dimension(2*lk+1):: M2
          real*8:: eps=EPSILON(result)



          result = 0d0 

          do i = 0,lk
          

              if (ind==0) then

                  if (i==0) then
          
                      if (DABS(q)>eps .and. DABS(M2(1))>eps) Then
                          Call T_lk(lk,0,"0",0,0,"0",t_k1)
                          
                          result = result + q*M2(1)*t_k1 
                      end if 
                  else

                      if (DABS(q)>eps .and. DABS(M2(2*i))>eps) Then
                          Call T_lk(lk,i,"c",0,0,"0",t_k1)
                          result = result + q*M2(2*i)*t_k1 
                      end if 
                      if ( DABS(q)>eps .and. DABS(M2(2*i+1))>eps) Then
                          Call T_lk(lk,i,"s",0,0,"0",t_k2)
                          result = result + q*M2(2*i+1)*t_k2
                      end if 
                      
                  end if
          ELSE
              if (i==0) then

                  if (DABS(q)>eps .and. DABS(M2(1))>eps) Then
                      Call T_lk(0,0,"0",lk,0,"0",t_k1)
                      result = result + q*M2(1)*t_k1 
                  end if 
              else

                      if (DABS(q)>eps .and. DABS(M2(2*i))>eps) Then
                          Call T_lk(0,0,"0",lk,i,"c",t_k1)
                          result = result + q*M2(2*i)*t_k1 
                      end if 
                      if ( DABS(q)>eps .and. DABS(M2(2*i+1))>eps) Then
                          Call T_lk(0,0,"0",lk,i,"s",t_k2)
                          result = result + q*M2(2*i+1)*t_k2
                      end if 
                      
                  end if
                  
              end if 
          end do

          RETURN
  END SUBROUTINE T_l0

    !********************************************************
  SUBROUTINE T_ll(M1,M2,lk1,lk2 , result)
      IMPLICIT NONE

      INTEGER ::  i,j
      real*8, INTENT(INOUT) ::result 
      Integer, INTENT(IN):: lk1,lk2
      real*8 :: t_k1
      real*8 , dimension(2*lk1+1):: M1
      real*8 , dimension(2*lk2+1):: M2
      real*8:: eps=EPSILON(result)
      Character :: cp1,cp2


          result = 0d0 

          do i = 1,2*lk1+1
                  do j = 1,2*lk2+1

                      if (DABS(M1(i))>eps .and. DABS(M2(j))>eps) Then

                          call Get_Comp(i,cp1)
                          call Get_Comp(j,cp2)

                          Call T_lk(lk1,Floor((i*1d0)/2d0),cp1,lk2,Floor((j*1d0)/2d0),cp2,t_k1)
                          result = result + M1(i)*M2(j)*t_k1 
                      end if 
                  
                  end do
          end do

      RETURN
  END SUBROUTINE T_ll




  SUBROUTINE MT(mu,k1,k2,res)
    IMPLICIT NONE

    Integer , INTENT(IN):: k1
    Character(len = 1) , INTENT(IN)::mu,k2
    real*8 ,INTENT(INOUT) :: res
        
    if (mu=="x") Then
        if (k1==0) then
            res =  0d0 ;
        else if (k1==1 .and. k2=="c") then
                res = sqrt(2d0);
        else if (k1==1 .and. k2=="s") then
                res = 0d0;    
        else
                res = k1*1d0;  
        end if 
      
        elseif (mu=="y")then
          if (k1==0) then
            res =  0d0 ;
          else if (k1==1 .and. k2=="s") then
                  res = sqrt(2d0);
          else if (k1==1 .and. k2=="c") then
                  res = 0d0;    
          elseif (k1>1 .and. k2=="c") then
                    res = -k1*1d0;
          elseif (k1>1 .and. k2=="s") then
                      res = k1*1d0;   
          end if
        else
          res=1d0  
      end if
  End   SUBROUTINE MT


  SUBROUTINE N_eta(mu,k1,k2,k_res1,k_res2)
    IMPLICIT NONE

    Integer , INTENT(IN):: k1
    Character(len = 1) , INTENT(IN)::mu,k2
    Integer ,INTENT(INOUT) :: k_res1
    Character(len = 1) , INTENT(INOUT)::k_res2

    if (k1==0) then
      k_res1 = 0
      k_res2 ="0"
    else

      if (mu=="x") Then
        if (k1==1) then
            k_res1 = 0
            k_res2 ="0"
        else if (k1>1 .and. k2=="c") then
            k_res1 = k1-1
            k_res2 ="c"  
        else
            k_res1 = k1-1
            k_res2 ="s"  
        end if 
      
        elseif (mu=="y")then
          if (k1==1) then
            k_res1 = 0
            k_res2 ="0"
          else if (k1>1 .and. k2=="c") then
              k_res1 = k1-1
              k_res2 ="s"  
          else
              k_res1 = k1-1
              k_res2 ="c"  
          end if 
        else
              k_res1 = k1
              k_res2 = k2  
      end if

    end if
        

  End   SUBROUTINE N_eta


  SUBROUTINE NN(la,ka1,lb,kb1,NN_fact)
    IMPLICIT NONE

    Integer , INTENT(IN):: la,ka1,lb,kb1
    real*8,INTENT(INOUT) :: NN_fact
    real*8 :: nf1,nf2,df1,df2

    if (la<0 .or. lb<0 .or. ka1<0 .or. kb1<0 .or. ka1 > la .or. kb1>lb) Then
      NN_fact  = 0d0
    else
      call Factorial(la+ka1,nf1)
      call Factorial(lb+kb1,nf2)
      call Factorial(la-ka1,df1)
      call Factorial(lb-kb1,df2)

      NN_fact  = DSQRT( (nf1/df1))*DSQRT((nf2/df2));
    end if
    
    

  End   SUBROUTINE NN


  SUBROUTINE Factorial(n,fact)
    IMPLICIT NONE

    integer, intent(in) :: n
    real*8, intent(inout) :: fact
    integer :: i


    fact = 1.0d0
    do i = 2, n
      fact = fact * (i*1d0)
    end do
    
    

  End   SUBROUTINE Factorial
  

  SUBROUTINE Get_Comp(i,cp)
              !! DECLARACTIONS
              IMPLICIT NONE

              INTEGER, INTENT(IN):: i
              Character, INTENT(INOUT):: cp
              !! EXECUTABLES
              if (i==1) then 
                  cp = "0"
              elseif(mod(i,2)==0)then
                  cp = "c"
              else
                  cp = "s"
              end if
          
              
          RETURN	!!ONLY NEEDED IF WE PLAN TO REACH END FUNCTION ALL OF THE TIME
  END SUBROUTINE Get_Comp







  SUBROUTINE Generate_Coordenates(coordenates)

 
      real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
      real*8, parameter :: pii = DACOS(-1.d0)  


      real*8 :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

          cal_coord(1) = coordenates(1)
            

          cal_coord(2)  = DCOS(coordenates(2)*pii/180d0)
          cal_coord(3) = DSIN(coordenates(2)*pii/180d0)

          cal_coord(4) = DCOS(coordenates(3)*pii/180d0)
          cal_coord(5) = DSIN(coordenates(3)*pii/180d0)
      
          cal_coord(6) = DCOS(coordenates(4)*pii/180d0)
          cal_coord(7) = DSIN(coordenates(4)*pii/180d0)

          cal_coord(8) = DCOS(coordenates(5)*pii/180d0)
          cal_coord(9) = DSIN(coordenates(5)*pii/180d0)

          cal_coord(10) = DCOS(coordenates(6)*pii/180d0)
          cal_coord(11) = DSIN(coordenates(6)*pii/180d0)


          cos_b1 =    cal_coord(2)
          sin_b1 =    cal_coord(3)
          cos_b2 =    cal_coord(4)
          sin_b2 =    cal_coord(5)
          cos_phi =   cal_coord(6)
          sin_phi =   cal_coord(7)
          cos_c1 =    cal_coord(8)
          sin_c1 =    cal_coord(9)
          cos_c2 =    cal_coord(10)
          sin_c2 =    cal_coord(11)


          

          Ar(1) = cos_b1           !Az
          Ar(2) = cos_c1*sin_b1    !Ay
          Ar(3) = sin_b1*sin_c1    !Ax
          
          
          Br(1) = -cos_b2           !Bz
          Br(2) = -cos_c2*sin_b2    !By
          Br(3) = -sin_b2*sin_c2    !Bx
          
          CC(1)= cos_b1*cos_b2 + cos_phi*sin_b1*sin_b2 !Czz
          CC(2)= cos_b2*sin_b1*sin_c1 - sin_b2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1)  !Cxz
          CC(3)= cos_b2*cos_c1*sin_b1 + sin_b2 *(-cos_phi*cos_b1*cos_c1 + sin_phi*sin_c1) !Cyz

          CC(4)= cos_c2*sin_phi*sin_b1 + (-cos_phi*cos_b2*sin_b1 + cos_b1*sin_b2)*sin_c2; !Czx
          CC(5)= -cos_b1*cos_c2*sin_phi*sin_c1 + (cos_b2*cos_c1*sin_phi + sin_b1*sin_b2*sin_c1)*sin_c2 &
          + cos_phi *(cos_c1*cos_c2 + cos_b1*cos_b2*sin_c1*sin_c2) !Cxx
          CC(6)= cos_c1*sin_b1*sin_b2*sin_c2 + cos_b1*cos_c1 *(-cos_c2*sin_phi + cos_phi*cos_b2*sin_c2) - &
          sin_c1 *(cos_phi*cos_c2 + cos_b2*sin_phi*sin_c2)   !Cyx

          CC(7)= -cos_phi*cos_b2*cos_c2*sin_b1 + cos_b1*cos_c2*sin_b2 - sin_phi*sin_b1*sin_c2 !Czy
          CC(8)= cos_c2*sin_b1*sin_b2*sin_c1 + cos_b2*cos_c2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1) &
          + (-cos_phi*cos_c1 + cos_b1*sin_phi*sin_c1)*sin_c2 !Cxy
          CC(9)= -cos_b2*cos_c2*sin_phi*sin_c1 + cos_c1 *(cos_c2*sin_b1*sin_b2 + cos_b1*sin_phi*sin_c2) &
          + cos_phi*(cos_b1*cos_b2*cos_c1*cos_c2 + sin_c1*sin_c2) !Cyy
        

  End   SUBROUTINE Generate_Coordenates

end module Geometry_Constant
  
