module Search
 

 character(len=2),private,parameter ::strNum(40) = (/"01","02","03","04","05","06","07","08","09","10",&
                                        "11","12","13","14","15","16","17","18","19","20",&
                                        "21","22","23","24","25","26","27","28","29","30",&
                                        "31","32","33","34","35","36","37","38","39","40"/)
 
 save
    Integer::initflag

    type strArrInd
        character(len=5), allocatable :: arr(:)
    end type

    type strArrDisp
        character(len=11), allocatable :: arr(:)
    end type


    type(strArrInd),public::IndArr(5,5)
    type(strArrDisp),public::DispArr(5,5,5,5)
  
 

 contains
   
    Subroutine Initialize_Search()
        IMPLICIT NONE
        Integer ::c
   
        if (initflag .NE. 2)then
            initflag = 2 
            call Init_Arr_Ind()
            call Init_Arr_Disp()
        end if
        ! write(*,*)"DispArr(1,2,1,1)"
        ! do c = 0,14
        !     write(*,*)DispArr(1,2,1,1)%arr(6*c+1),"    ",DispArr(1,2,1,1)%arr(6*c+2),"    ",DispArr(1,2,1,1)%arr(6*c+3),&
        !      "    ",DispArr(1,2,1,1)%arr(6*c+4),"    ",DispArr(1,2,1,1)%arr(6*c+5),"    ",DispArr(1,2,1,1)%arr(6*c+6)
        ! end do 
    end Subroutine Initialize_Search

    SUBROUTINE Init_Arr_Ind()
        
        IMPLICIT NONE
        integer i,j

        do i=1,5
            do j=i,5
                call Gen_Arr_Ind(i,j,IndArr(i,j)%arr)
            end do
        end do


    END SUBROUTINE Init_Arr_Ind

    SUBROUTINE Init_Arr_Disp()
        
        IMPLICIT NONE
        integer i,j,k,l

        do i=1,5
            do j=i,5
                do k=1,5
                    do l=k,5
                            Call Gen_Arr_Disp(i,j,k,l,DispArr(i,j,k,l)%arr) 
                    end do
                end do
            end do
        end do
        
        

    END SUBROUTINE Init_Arr_Disp

    SUBROUTINE Gen_Arr_Ind(i,j,arr)
        
        IMPLICIT NONE
        
        Integer , INTENT(IN) :: i,j

        integer :: n,ci,cj,counter
        character(len=5),allocatable, INTENT(OUT)::arr(:)


        if (i==j)then
            n = (i+1)*(2*i+1)
        else
            n = (2*j+1)*(2*i+1)
        end if 

        Allocate(arr(n))
        

        if (i==j)then
                counter = 0
                do ci=1,2*i+1
                    do cj=ci,2*i+1
                        counter = counter+1
                        arr(counter) = strNum(ci)//"_"//strNum(cj)
                    end do
                end do
        else
                counter = 0
                do ci=1,2*i+1
                    do cj=1,2*j+1
                        counter = counter+1
                        arr(counter) = strNum(ci)//"_"//strNum(cj)
                    end do
                end do

        end if

    END SUBROUTINE Gen_Arr_Ind

    SUBROUTINE Gen_Arr_Disp(i,j,k,l,arr)
        
        IMPLICIT NONE
        
        Integer , INTENT(IN) :: i,j,k,l

        integer :: n,nA,nB,Nmax,Nmin,counter,cl,cm
        character(len=11),allocatable, INTENT(OUT)::arr(:)


        if (i==j)then
            nA = (i+1)*(2*i+1)
        else
            nA = (2*j+1)*(2*i+1)
        end if 

        if (k==l)then
            nB = (k+1)*(2*k+1)
        else
            nB = (2*k+1)*(2*l+1)
        end if 

        if (i==j .and. k==l)then
            n = (i+1)*(2*i+1)*(k+1)*(2*k+1)
        elseif (i.NE.j .and. k==l)then
            n = (2*i+1)*(2*j+1)*(k+1)*(2*k+1)
        elseif (i==j .and. k.NE.l)then
            n = (i+1)*(2*i+1)*(2*k+1)*(2*l+1)   
        else
            n = (2*i+1)*(2*j+1)*(2*k+1)*(2*l+1)  
        end if 

        Allocate(arr(n))
        
        counter=0
        if(nB>nA)then
            do cm=1,nB
                do cl=1,nA
                    counter=counter+1
                    arr(counter) = IndArr(i,j)%arr(cl)//"-"//IndArr(k,l)%arr(cm)
                end do
            end do
        else
            do cl=1,nA
                do cm=1,nB
                    counter=counter+1
                    arr(counter) = IndArr(i,j)%arr(cl)//"-"//IndArr(k,l)%arr(cm)
                end do
            end do
        end if

        

    END SUBROUTINE Gen_Arr_Disp

    SUBROUTINE Get_Disp_Index(i,j,k,l,ci,cj,ck,cl,index)
        
        IMPLICIT NONE
        
        Integer , INTENT(IN) :: i,j,k,l,ci,cj,ck,cl
        integer, INTENT(OUT)::index
        Character(len = 11) ::str1,str2,str3,str4,dispStr
        integer::n,s
       

        str1 =  strNum(ci)//"_"//strNum(cj)//"-"//strNum(ck)//"_"//strNum(cl)
        str2 =  strNum(ci)//"_"//strNum(cj)//"-"//strNum(cl)//"_"//strNum(ck)
        str3 =  strNum(cj)//"_"//strNum(ci)//"-"//strNum(ck)//"_"//strNum(cl)
        str4 =  strNum(cj)//"_"//strNum(ci)//"-"//strNum(cl)//"_"//strNum(ck)



        if (i==j .and. k==l)then
            n = (i+1)*(2*i+1)*(k+1)*(2*k+1)
        elseif (i.NE.j .and. k==l)then
            n = (2*i+1)*(2*j+1)*(k+1)*(2*k+1)
        elseif (i==j .and. k.NE.l)then
            n = (i+1)*(2*i+1)*(2*k+1)*(2*l+1)   
        else
            n = (2*i+1)*(2*j+1)*(2*k+1)*(2*l+1)  
        end if 

       

        do s = 1,n 
            dispStr = DispArr(i,j,k,l)%arr(s)
            if (i==j ) then
                if ( k==l) then
                ! example mm_mm
                    if (     dispStr==str1 &
                        .Or. dispStr==str2 & 
                        .Or. dispStr==str3 &
                        .Or. dispStr==str4 ) Then  

                        index = s 
                        Exit
                    end if 
                else
                    if (     dispStr==str1 &
                        .Or. dispStr==str3 ) Then  ! example mm_mQd

                        index = s 
                        Exit
                    end if 
                end if
            else
                if ( k==l) then
                    if (     dispStr==str1 &
                        .Or. dispStr==str2 ) Then   ! example mQd_mm

                        index = s 
                        Exit
                    end if 
                else
                    if (     dispStr==str1 ) Then   ! example mQd_mQd

                        index = s 
                        Exit
                    end if 
                end if

            end if
        end do 
        

    END SUBROUTINE Get_Disp_Index

 end module Search

 



! SUBROUTINE Gen_Arr_Disp(i,j,k,l,arr)
!     Integer , INTENT(IN) :: i,j,k,l
!     character(len=3),allocatable, INTENT(OUT)::arr(:)
!     integer :: n
!     character()
    

    
    
    
! END SUBROUTINE Gen_Arr_Disp

! ! SUBROUTINE GetIndex_Disp (n,i,k,j,l,ci,ck,cj,cl,index)

! !     Character(len = 2) , dimension(6)::alpha_arr 
! !     Character(len = 20) , dimension(36)::disp_arr 
! !     Integer , INTENT(IN) :: i,j,g,h
! !     Integer, INTENT(INOUT) :: index

! !     Integer :: s,l,m,count
! !     Character(Len=1)::stri,strj,strg,strh
! !     Character(len = 10) ::str1,str2,str3,str4

! !     alpha_arr(1) = "11"
! !     alpha_arr(2) = "12"
! !     alpha_arr(3) = "13"
! !     alpha_arr(4) = "22"
! !     alpha_arr(5) = "23"
! !     alpha_arr(6) = "33"

! !     count=1;

! !     do l=1,6
! !         do m=1,6
! !             disp_arr(count) = alpha_arr(l)//alpha_arr(m)
! !             count=count+1
! !         end do
! !     end do

! !     call GetStr(i,stri)
! !     call GetStr(j,strj)
! !     call GetStr(g,strg)
! !     call GetStr(h,strh)

! !     str1 =  strg//strh//stri//strj
! !     str2 =  strh//strg//stri//strj
! !     str3 =  strg//strh//strj//stri
! !     str4 =  strh//strg//strj//stri

    

! !     do s = 1,36 
        
! !         if (disp_arr(s)==str1 .Or. disp_arr(s)==str2 .Or. disp_arr(s)==str3 .Or. disp_arr(s)==str4) Then

! !             index = s
! !         end if 
! !     end do 



    
! ! END SUBROUTINE GetIndex_mm_mm











! SUBROUTINE GetStr(i,stri)

!     Integer , INTENT(IN) :: i

!     Character(len = 1), INTENT(INOUT) ::stri

!     if (i==0)then
!         stri ="0" 
!         Return
!     elseif (i==1)then
!         stri ="1" 
!         Return
!     elseif (i==2)then
!         stri ="2" 
!         Return
!     elseif (i==3)then
!         stri ="3" 
!         Return
!     elseif (i==4)then
!         stri ="4" 
!         Return                
!     elseif (i==5)then
!         stri ="5" 
!         Return
!     elseif (i==6)then
!         stri ="6" 
!         Return
!     elseif (i==7)then
!         stri ="7" 
!         Return
!     elseif (i==8)then
!         stri ="8" 
!         Return
!     elseif (i==9)then
!         stri ="9" 
!         Return
!     endif

    
! END SUBROUTINE GetStr