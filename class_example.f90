PROGRAM test
    USE FitConstants
  IMPLICIT NONE
  integer::i

    Call Evaluate()
    ! DO i=1,5
    !     write(*,*)"filename  ",C(i)%filename
        
    ! end do

    
    ! Call Find_Coeff_Set('./files/coefficients_Adp.txt',5,C,ind(1))
    ! Call Find_Coeff_Set('./files/coefficients.txt',5,C,ind(2))
    ! Call Find_Coeff_Set('./files/coefficients_Ap.txt',5,C,ind(3))

    ! Call Last_Coeff_Set(5,C,lastIndex)
    ! write(*,*)"lastIndex  ",lastIndex

    ! write(*,*)ind
  
END PROGRAM


SUBROUTINE Evaluate()
    USE FitConstants
    IMPLICIT NONE
    integer::N,i,indx1,indx2,indx3

    
    ! CALL Coeff_(1)%Initializer('./files/coefficients_Adp.txt')
    ! CALL Coeff_(2)%Initializer('./files/coefficients_Ap.txt')
    ! CALL Coeff_(3)%Initializer('./files/coefficients.txt')

    ! N=10


    DO i=1,10
        
        Call Get_Coeff_Index('./files/coefficients_Adp.txt',indx1)
        Call Get_Coeff_Index('./files/coefficients_Ap.txt',indx2)
        Call Get_Coeff_Index('./files/coefficients.txt',indx3)

        write(*,*)"Indexes: ",indx1,indx2,indx3
    END DO




END SUBROUTINE Evaluate





