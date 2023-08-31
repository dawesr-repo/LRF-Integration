!FittingConstant is the module in charge of all constants related with the fit 
!It can handle several coefficient files at the same time


MODULE FitConstants
    save
    public
    real*8, parameter ::  C1=627.5095d0
    real*8, parameter ::  C2=0.529177249d0
    real*8, parameter ::  C3=349.757d0

    TYPE Fit_Contant
     character(:), allocatable :: filename
     real*8 :: Zero
     Integer::initflag
     Integer, dimension (8) :: M_Fit     
     Integer, dimension (3):: D_Fit     
     Integer, dimension (5):: I_Fit     
     Integer, dimension (2):: H_Fit     
     
     !Multipoles !

     real*8 , dimension(64)   :: A_Mult,B_Mult 

     !Polarizability!

     real*8 , dimension(57)   :: A_Pol,B_Pol 

     !Hyperpolarizability!
    
     real*8 , dimension(40)   :: A_HPol,B_HPol 

     !Dispersion!

     real*8 , dimension(873)   :: Disp

     CONTAINS
        PROCEDURE, PASS :: Initializer
        PROCEDURE, PASS :: Read_Parameters
    END TYPE
    Integer,parameter :: NArray = 5
    TYPE(Fit_Contant) :: Coeff(NArray)
    
   
CONTAINS
  SUBROUTINE Initializer(this,filename)
    IMPLICIT NONE
    CLASS(Fit_Contant), INTENT(OUT) :: this
    Character(len=*), INTENT(IN) ::filename

        this%initflag = 1
        this%filename = filename

  END SUBROUTINE Initializer

  SUBROUTINE Read_Parameters(this)
    IMPLICIT NONE
    CLASS(Fit_Contant), INTENT(InOut) :: this

    Character(len = 200) :: row

    

    if (this%initflag==1)Then
       
       
        this%initflag = 2
        !write(*,*)"initflag",this%initflag
        Open( 10, file = this%filename )

        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row

        read(10, *)  this%M_Fit
        read(10, *)  this%D_Fit
        read(10, *)  this%I_Fit
        read(10, *)  this%H_Fit

        Read( 10, *) row
        read(10, *)  this%Zero

        Read( 10, *) row
        read(10, *) this%A_Mult(1)          !q
        read(10, *) this%A_Mult(2:4)        !m
        read(10, *) this%A_Mult(5:9)        !Qd
        read(10, *) this%A_Mult(10:16)      !O
        read(10, *) this%A_Mult(17:25)      !Phi
        read(10, *) this%A_Mult(26:36)      !M5
        read(10, *) this%A_Mult(37:49)      !M6
        read(10, *) this%A_Mult(50:64)      !M7

        Read( 10, *) row
        read(10, *) this%B_Mult(1)          !q
        read(10, *) this%B_Mult(2:4)        !m
        read(10, *) this%B_Mult(5:9)        !Qd
        read(10, *) this%B_Mult(10:16)      !O
        read(10, *) this%B_Mult(17:25)      !Phi
        read(10, *) this%B_Mult(26:36)      !M5
        read(10, *) this%B_Mult(37:49)      !M6
        read(10, *) this%B_Mult(50:64)      !M7

        Read( 10, *) row
        read(10, *) this%A_Pol(1:6)         !mm
        read(10, *) this%A_Pol(7:21)        !Qdm
        read(10, *) this%A_Pol(22:36)       !QdQd
        read(10, *) this%A_Pol(37:57)       !Om

        Read( 10, *) row
        read(10, *)  this%B_Pol(1:6)        !mm
        read(10, *)  this%B_Pol(7:21)       !Qdm
        read(10, *)  this%B_Pol(22:36)      !QdQd
        read(10, *)  this%B_Pol(37:57)      !Om

        Read( 10, *) row
        read(10, *)  this%A_HPol(1:10)      !mmm
        read(10, *)  this%A_HPol(11:40)     !Qdmm

        
        Read( 10, *) row
        read(10, *)  this%B_HPol(1:10)      !mmm
        read(10, *)  this%B_HPol(11:40)     !Qdmm

        Read( 10, *) row
        read(10, *) this%Disp(1:36)         !mm_mm      
        read(10, *) this%Disp(37:126)       !Qdm_mm 
        read(10, *) this%Disp(127:216)      !mm_Qdm 
        read(10, *) this%Disp(217:342)      !Om_mm 
        read(10, *) this%Disp(343:468)      !mm_Om 
        read(10, *) this%Disp(469:558)      !QdQd_mm
        read(10, *) this%Disp(559:648)      !mm_QdQd
        read(10, *) this%Disp(649:873)      !Qdm_Qdm


        close(10)



    else
        
    end if

  END SUBROUTINE Read_Parameters

  SUBROUTINE Find_Coeff_Set(filename,ind)
    IMPLICIT NONE
    Character(*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: ind 
    integer:: i
    ind = -1

    if (NArray<1)Then
        Return
    else
        do i=1,NArray
            if (Coeff(i)%filename == filename)then
                ind = i 
                return
            end if
        end do 

    end if

  END SUBROUTINE Find_Coeff_Set

  SUBROUTINE Last_Coeff_Set(lastIndex)
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: lastIndex
        integer:: i
        lastIndex = 0

        if (NArray<1)Then
            Return
        else
            
            do i=1,NArray
                if (LEN(Coeff(i)%filename)<1)then
                    lastIndex = i-1 
                    return
                end if
            end do 

            if (lastIndex==NArray)then
            write(*,*)"The maximun number of coefficients sets is :",NArray,&
                        "to change the maximun go to module MODULE Fit and change NARRAY" 
                  lastIndex=-10      
            end if
        end if


  END SUBROUTINE Last_Coeff_Set

  SUBROUTINE Get_Coeff_Index(filename,indx)
       IMPLICIT NONE
       Character(*), INTENT(IN) :: filename
       Integer, INTENT(OUT) :: indx
       integer::ind,lastIndex
       
       Call Find_Coeff_Set(filename,ind)

       
       if (ind < 1)then

        Call Last_Coeff_Set(lastIndex)

        indx = lastIndex + 1
        CALL Coeff(indx)%Initializer(filename)
        Call Coeff(indx)%Read_Parameters()

       else
        indx = ind
        return 
       end if 
  
  END SUBROUTINE Get_Coeff_Index

END module FitConstants