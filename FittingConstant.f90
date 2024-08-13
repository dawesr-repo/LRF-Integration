!FittingConstant is the module in charge of all constants related with the fit 
!It can handle several coefficient files at the same time


MODULE FitConstants
    save
    public
    real*8, parameter ::  C1=627.5095d0
    real*8, parameter ::  C2=0.529177249d0
    real*8, parameter ::  C3=349.757d0

    Integer, parameter::  polArr_ij(5,5) = transpose(reshape((/ 1, 7, 37, 0, 0&
                                                              , 0, 22, 0, 0, 0&
                                                              , 0, 0, 0, 0,  0&
                                                              , 0, 0, 0, 0,  0&
                                                              , 0, 0, 0, 0,  0/), shape(polArr_ij))) 

    Integer, parameter::  hpolArr_ijk(3,3,3) = reshape((/ 1, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                                   11, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                                   0, 0, 0, 0, 0, 0, 0, 0, 0&
                                                                   /), shape(hpolArr_ijk))                                                           

    Integer, parameter::  disp_ijkl(3,3,3,3) = reshape((/ 1, 0, 0, 37, 469, 0, 217, 0, 0,&
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                          127, 0, 0, 649, 0, 0, 0, 0, 0,&
                                                          559, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                          343, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0,&
                                                          0, 0, 0, 0, 0, 0, 0, 0, 0&
                                                                   /), shape(disp_ijkl))  

    !polArr_ij = 0
    ! polArr_ij(1,:)=(/1,7,37,0,0/)
    ! polArr_ij(2,:)=(/0,22,0,0,0/)

    TYPE Fit_Contant
     character(:), allocatable :: filename
     real*8 :: Zero
     Integer::initflag

     Integer, dimension (15) :: M_Fit     
     Integer, dimension (15) :: D_Fit     
     Integer, dimension (15) :: I_Fit   

     INTEGER :: max_T

     !Multipoles !

     real*8 , dimension(225)   :: A_Mult,B_Mult 

     !Polarizability!

     real*8 , dimension(6,7,195)   :: A_Pol,B_Pol 


     !Dispersion!

     real*8 , dimension(5,10,5,10,3087)   :: Disp

     CONTAINS
        PROCEDURE, PASS :: Initializer
        PROCEDURE, PASS :: Read_Parameters
        ! PROCEDURE, PASS :: Get_Mult_Comp
        ! PROCEDURE, PASS :: Get_Pol_Comp
        ! PROCEDURE, PASS :: Get_Disp_Comp
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

        this%max_T = MAXVAL([iord-2,mord,dord-3])

        
        read(10, *)row, this%A_Mult(1:mord**2)  !A_Mult        
        read(10, *)row, this%B_Mult(1:mord**2)  !B_Mult   

        if (iord>=4) then
            do i=1,iord-3
                do j=i,iord-3
                    if (i+j<=iord-2) then
                        ln = (2*i+1)*(2*j+1);
                        read(10, *)row, this%A_Pol(i,j,1:ln)  !PA_i-j
                        read(10, *)row, this%B_Pol(i,j,1:ln)  !PB_i-j
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

!   SUBROUTINE Get_Mult_Comp(this,Molec,i,arr)

!     IMPLICIT NONE
!     CLASS(Fit_Contant), INTENT(IN) :: this
!     Character(len=*), INTENT(IN) ::Molec
!     Integer, INTENT(IN) ::i
!     real*8, INTENT(OUT)  :: arr(2*i+1) 
!     real*8 , dimension(64)   :: Mult 
!     integer::n,init
    
!     init = i**2+1
!     n = 2*i+1
    
!     if (Molec=="A")then
!         Mult = this%A_Mult
!     else
!         Mult = this%B_Mult    
!     end if

!     arr = Mult(init:init+n-1)

!   end SUBROUTINE Get_Mult_Comp

!   SUBROUTINE Get_Pol_Comp(this,Molec,i,j,arr,n)

!     IMPLICIT NONE
!     CLASS(Fit_Contant), INTENT(IN) :: this
!     Character(len=*), INTENT(IN) ::Molec
!     Integer, INTENT(IN) ::i,j
!     real*8, allocatable, INTENT(OUT)  :: arr(:) 
!     real*8 , dimension(57)   :: Pol 
!     integer::init
!     integer, INTENT(OUT)::n


!     init = polArr_ij(i,j)
   
!     if (Molec=="A")then
!         Pol = this%A_Pol
!     else
!         Pol = this%B_Pol    
!     end if

!     if (i==j)then
!         n = (i+1)*(2*i+1)
!     else
!         n = (2*j+1)*(2*i+1)
!     end if 
    
!     Allocate(arr(n)) 
!     arr = Pol(init: n + init-1)

!   end SUBROUTINE Get_Pol_Comp


!   SUBROUTINE Get_Disp_Comp(this,i,j,k,l,arr,n)

!     IMPLICIT NONE
!     CLASS(Fit_Contant), INTENT(IN) :: this
!     Integer, INTENT(IN) ::i,j,k,l
!     real*8, allocatable, INTENT(OUT)  :: arr(:)    
!     integer, INTENT(OUT)::n

!     integer::init

!     init = disp_ijkl(i,j,k,l)


!     if (i==j .and. k==l)then
!         n = (i+1)*(2*i+1)*(k+1)*(2*k+1)
!     elseif (i.NE.j .and. k==l)then
!         n = (2*i+1)*(2*j+1)*(k+1)*(2*k+1)
!     elseif (i==j .and. k.NE.l)then
!         n = (i+1)*(2*i+1)*(2*k+1)*(2*l+1)   
!     else
!         n = (2*i+1)*(2*j+1)*(2*k+1)*(2*l+1)  
!     end if 


    
!     Allocate(arr(n))
   
!     arr = this%Disp(init: n + init-1)

!   end SUBROUTINE Get_Disp_Comp


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