
!********************************************************
RECURSIVE SUBROUTINE T_lk(Ar,Br,C,la,ka1,ka2,lb,kb1,kb2,Tlk)
    use Tensors_constant
    IMPLICIT NONE
    
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 , dimension(9):: Cab
    Integer , INTENT(IN):: la,ka1,lb,kb1
    Character(len = 1) , INTENT(IN)::ka2,kb2
    Character(len = 20), dimension(3):: coordinates = ["z","y","x"]
    Integer :: i ,j   ,k_res1,ka_res1,kb_res1
    real*8 , INTENT(INOUT):: Tlk
    real*8 :: M,M1,M2,NN_,Tlk_new,NN_1,NN_2,Tlk_1,Tlk_2
    Character (len = 1)::k_res2,ka_res2,kb_res2


    real*8:: eps = 1.0D-8 ,T,chang




    Cab(1) = C(1) !zz
    Cab(2) = C(7) !zy
    Cab(3) = C(4) !zx
    Cab(4) = C(3) !yz
    Cab(5) = C(9) !yy
    Cab(6) = C(6) !yx
    Cab(7) = C(2) !xz
    Cab(8) = C(8) !xy
    Cab(9) = C(5) !xx


    IF (ka1<0 .or. kb1<0 .or. la<0 .or. lb<0 .or. ka1>la .or. kb1>lb) THEN
        Tlk = 0d0
        RETURN
    ELSE IF (la==0 .and. lb==0) THEN  
        Tlk = 1d0
        RETURN
    ELSE

      !call change_Tensors(la,ka1,ka2,lb,kb1,kb2,3d0) 
      call get_Tensor(la,ka1,ka2,lb,kb1,kb2,T,chang) 

      ! if (la==0 .and. ka1==0 .and. ka2=="0" .and. &
      !     lb==1 .and. kb1==0 .and. kb2=="0")then
      !   write(*,*)"found it",T,chang
      ! end if 

      if (DABS(chang-1d0)<eps)Then
          !write(*,*)"T from constant!",T,chang,la,ka1,ka2,lb,kb1,kb2
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
              call T_lk(Ar,Br,C,la-1,k_res1,k_res2,0,0,"0",Tlk_new)
              Tlk = Tlk + ((2.d0*la-1d0)/(la*1.d0))*M*Ar(i)*Tlk_new*NN_
            end if
          end do
        
          call NN(la-2,ka1 ,0,0,NN_)

          if (abs(la-1d0) > eps ) Then
              call T_lk(Ar,Br,C,la-2,ka1,ka2,0,0,"0",Tlk_new)
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

            call T_lk(Ar,Br,C,0,0,"0",lb-1,k_res1,k_res2,Tlk_new)
            Tlk = Tlk + ((2.d0*lb-1d0)/(lb*1.d0))*M*Br(i)*Tlk_new*NN_

          end if 
        end do
      
        if (abs(1d0*lb-1d0) > eps ) Then 
          call NN(lb-2,kb1 ,0,0,NN_)
          call T_lk(Ar,Br,C,0,0,"0",lb-2,kb1,kb2,Tlk_new)
          
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
        Call T_lk(Ar,Br,C,la-2,ka1,ka2,lb,kb1,kb2,Tlk_1)
        

        Tlk  =Tlk  + Tlk_1*NN_1

        if (abs(2d0*la+1.d0*lb-1d0) > eps ) Then 

          Call T_lk(Ar,Br,C,la,ka1,ka2,lb-2,kb1,kb2,Tlk_2)
          Tlk  =Tlk  - ((2d0*la+1.d0*lb-1d0)/lb*1d0)*Tlk_2*NN_2;
        end if

        do i = 1,3

            Call N_eta(coordinates(i),kb1,kb2,k_res1,k_res2);
            call MT(coordinates(i),kb1,kb2,M)
            Call NN(la,ka1,lb-1,k_res1,NN_);

            if (abs(M) > eps .and. abs(2d0*la+2d0*lb-1d0) > eps .and. abs(Br(i))>eps) Then

              Call T_lk(Ar,Br,C,la,ka1,ka2,lb-1,k_res1,k_res2,Tlk_new)
            
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
                
                  Call T_lk(Ar,Br,C,la-1,ka_res1,ka_res2,lb-1,kb_res1,kb_res2,Tlk_new)

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

    NN_fact  = DSQRT( (nf1/df1)*(nf2/df2));
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

