!FittingConstant is the module in charge of all constants related with the fit 
!It can handle several coefficient files at the same time


MODULE FitConstants
    save
    public
    real*8, parameter ::  C1=627.5095d0
    real*8, parameter ::  C2=0.529177249d0
    real*8, parameter ::  C3=349.757d0

    INTEGER :: max_T = 0


    TYPE Fit_Contant
     character(:), allocatable :: filename
     real*8 :: Zero
     Integer::initflag

     Integer, dimension (15) :: M_Fit     
     Integer, dimension (15) :: D_Fit     
     Integer, dimension (15) :: I_Fit   


     !Multipoles !

     real*8 , dimension(225)   :: A_Mult,B_Mult 

     !Polarizability!

     real*8 , dimension(6,12,195)  :: A_Pol,B_Pol 


     !Dispersion!

     real*8 , dimension(5,10,5,10,3087)   :: Disp

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
    Integer::iord,mord,dord
    Integer::i,j,l1,l2,t1,t2,ln
    real*8 :: polarr(195)

    if (this%initflag==1)Then
       
       
        this%initflag = 2
      
        Open( 10, file = this%filename )

        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row
        Read( 10, *) row


        Read(10, *)  row,this%Zero
        Read( 10, *) row

        read(10, *)  row,this%M_Fit
        read(10, *)  row,this%I_Fit
        read(10, *)  row,this%D_Fit

        iord = MAXVAL(this%I_Fit);
        mord = MAXVAL(this%M_Fit);
        mord = MAXVAL([mord,iord-3]);
        dord = MAXVAL(this%D_Fit);

        max_T = MAXVAL([max_T,iord-2,mord,dord-3])

        
        read(10, *)row, this%A_Mult(1:mord**2)  !A_Mult        
        read(10, *)row, this%B_Mult(1:mord**2)  !B_Mult 



        if (iord>=4) then
            do i=1,iord-3
                do j=i,iord-3
                    if (i+j<=iord-2) then
                        ln = (2*i+1)*(2*j+1);
                        read(10, *)row, this%A_Pol(i,j,1:ln)  !PA_i-j
                        read(10, *)row,  this%B_Pol(i,j,1:ln) !PB_i-j
                    end if
                end do
            end do
        end if

        if (dord>=6) then
            do l1=1,dord-5
                do l2=l1,dord-5
                    do t1=1,dord-5
                        do t2=t1,dord-5
                    
                            if (l1+l2+t1+t2<=dord-2) then
                                ln = (2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1);
                                read(10, *)row, this%Disp(l1,l2,t1,t2,1:ln)  !Dispersion l1,l2, t1,t2                                                           disp_kk_vv_coeff{l1,l2,t1,t2});
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        
        close(10)

        
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