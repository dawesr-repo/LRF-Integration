
!  Fortran subroutine for 0.7.1 matlab version Optimization in running time
!  Software version 3.1.1
! last test: 20/08/2023


SUBROUTINE Coordinate_Transformation(coordenates,coord_format,new_coordinates)
    IMPLICIT NONE
    real*8 ,dimension(6), INTENT(IN):: coordenates
    
    Character(len = 20), INTENT(IN) :: coord_format
    real*8 , dimension(6), INTENT(INOUT) :: new_coordinates
    
    new_coordinates=coordenates;

    
    if (coord_format == "Euler_ZYZ") then
        new_coordinates(5) = coordenates(5) - 90
        new_coordinates(6) = coordenates(6) - 90
    end if 
    if (coord_format == "Spherical") then
        new_coordinates(5) = 90 - coordenates(5) 
        new_coordinates(6) = 90 -coordenates(6) 
    end if 
    
End   SUBROUTINE Coordinate_Transformation

SUBROUTINE General_Coordinates_Format(Dim, old_coordinates, general_coodinates)
            IMPLICIT NONE
    
            INTEGER,INTENT(IN) ::  Dim
            real*8 , dimension(6), INTENT(INOUT) ::general_coodinates
            real*8 , dimension (Dim),INTENT(IN) :: old_coordinates
       

            if (Dim==2) then
                    general_coodinates(1) = old_coordinates(1)  !R
                    general_coodinates(2) = old_coordinates(2)  !b1
                    general_coodinates(3) = 0d0  !b2
                    general_coodinates(4) = 0d0  !phi
                    general_coodinates(5) = 0d0  !c1
                    general_coodinates(6) = 0d0  !c2

                elseif (Dim==3) then
                    general_coodinates(1) = old_coordinates(1)  !R
                    general_coodinates(2) = old_coordinates(2)  !b1
                    general_coodinates(3) = 0d0  !b2
                    general_coodinates(4) = 0d0  !phi
                    general_coodinates(5) = old_coordinates(3)  !c1
                    general_coodinates(6) = 0d0  !c2

                elseif (Dim==4) then   
                    general_coodinates(1) = old_coordinates(1)  !R
                    general_coodinates(2) = old_coordinates(2)  !b1
                    general_coodinates(3) = old_coordinates(3)  !b2
                    general_coodinates(4) = old_coordinates(4)  !phi
                    general_coodinates(5) = 0d0  !c1
                    general_coodinates(6) = 0d0  !c2 

                elseif (Dim==5) then   
                        general_coodinates(1) = old_coordinates(1)  !R
                        general_coodinates(2) = old_coordinates(2)  !b1
                        general_coodinates(3) = old_coordinates(3)  !b2
                        general_coodinates(4) = old_coordinates(4)  !phi
                        general_coodinates(5) = old_coordinates(5)  !c1
                        general_coodinates(6) = 0d0  !c2     
                elseif (Dim==6) then   
                        general_coodinates(1) = old_coordinates(1)  !R
                        general_coodinates(2) = old_coordinates(2)  !b1
                        general_coodinates(3) = old_coordinates(3)  !b2
                        general_coodinates(4) = old_coordinates(4)  !phi
                        general_coodinates(5) = old_coordinates(5)  !c1
                        general_coodinates(6) = old_coordinates(6)  !c2       

            end if
  
            RETURN
END SUBROUTINE General_Coordinates_Format

SUBROUTINE GetIndex_mm (i,j,index)

    Character(len = 20) , dimension(6)::alpha_arr 
    Integer , INTENT(IN) :: i , j
    Integer, INTENT(INOUT) :: index

    Integer :: s
    
    Character(len = 5) ::strIJ,strJI
    Character(Len=1)::stri,strj



    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"


    call GetStr(i,stri)
    call GetStr(j,strj)

    strIJ = stri//strj
    strJI = strj//stri


  

    index=1
    do s = 1,6 
        
        if (alpha_arr(s)==strIJ .Or. alpha_arr(s)==strJI) Then
            index = s;
        end if 
    end do 

END SUBROUTINE GetIndex_mm 


SUBROUTINE GetIndex_mQd (i,j,index)

    Character(len = 20) , dimension(15)::alpha_Qd
    Integer , INTENT(IN) :: i , j
    Integer, INTENT(INOUT) :: index
    Integer :: s,u,v,count
    Character(Len=1)::stri,strj
    Character(len = 5) ::strIJ



    alpha_Qd(1) = "11"
    alpha_Qd(2) = "12"
    alpha_Qd(3) = "13"
    alpha_Qd(4) = "14"
    alpha_Qd(5) = "15"
    alpha_Qd(6) = "21"
    alpha_Qd(7) = "22"
    alpha_Qd(8) = "23"
    alpha_Qd(9) = "24"
    alpha_Qd(10) = "25"
    alpha_Qd(11) = "31"
    alpha_Qd(12) = "32"
    alpha_Qd(13) = "33"
    alpha_Qd(14) = "34"
    alpha_Qd(15) = "35"




    call GetStr(i,stri)
    call GetStr(j,strj)

    strIJ = stri//strj
  

    do s = 1,15 
        if (alpha_Qd(s)==strIJ ) Then
            index = s;
        end if 
    end do 

END SUBROUTINE GetIndex_mQd

SUBROUTINE GetIndex_mO (i,j,index)

    Character(len = 20) , dimension(21)::alpha_mO
    Integer , INTENT(IN) :: i , j
    Integer, INTENT(INOUT) :: index

    Integer :: s,u,v,count
    Character(Len=1)::stri,strj
    Character(len = 5) ::strIJ



    alpha_mO(1) = "11"
    alpha_mO(2) = "12"
    alpha_mO(3) = "13"
    alpha_mO(4) = "14"
    alpha_mO(5) = "15"
    alpha_mO(6) = "16"
    alpha_mO(7) = "17"
    alpha_mO(8) = "21"
    alpha_mO(9) = "22"
    alpha_mO(10) = "23"
    alpha_mO(11) = "24"
    alpha_mO(12) = "25"
    alpha_mO(13) = "26"
    alpha_mO(14) = "27"
    alpha_mO(15) = "31"
    alpha_mO(16) = "32"
    alpha_mO(17) = "33"
    alpha_mO(18) = "34"
    alpha_mO(19) = "35"
    alpha_mO(20) = "36"
    alpha_mO(21) = "37"




    call GetStr(i,stri)
    call GetStr(j,strj)

    strIJ = stri//strj
  

    do s = 1,21 
        if (alpha_mO(s)==strIJ ) Then
            index = s;
        end if 
    end do 

END SUBROUTINE GetIndex_mO


SUBROUTINE GetIndex_QdQd (i,j,index)

    Character(len = 20) , dimension(15)::alpha_arr 
    Integer , INTENT(IN) :: i , j
    Integer, INTENT(INOUT) :: index

    Integer :: s
    Character(Len=1)::stri,strj
    Character(len = 5) ::strIJ,strJI



    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "14"
    alpha_arr(5) = "15"
    alpha_arr(6) = "22"
    alpha_arr(7) = "23"
    alpha_arr(8) = "24"
    alpha_arr(9) = "25"
    alpha_arr(10) = "33"
    alpha_arr(11) = "34"
    alpha_arr(12) = "35"
    alpha_arr(13) = "44"
    alpha_arr(14) = "45"
    alpha_arr(15) = "55"



    call GetStr(i,stri)
    call GetStr(j,strj)

    strIJ = stri//strj
    strJI = strj//stri
  


    do s = 1,15 
        ! write(*,*)alpha_arr(s) , strIJ, strJI ,alpha_arr(s)==strIJ , alpha_arr(s)==strJI
        if (alpha_arr(s)==strIJ .Or. alpha_arr(s)==strJI) Then
            index = s;
        end if 
    end do 

END SUBROUTINE GetIndex_QdQd 


SUBROUTINE GetIndex_mm_mm (g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 20) , dimension(36)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index

    Integer :: s,l,m,count
    Character(Len=1)::stri,strj,strg,strh
    Character(len = 10) ::str1,str2,str3,str4

    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"

    count=1;

    do l=1,6
        do m=1,6
            disp_arr(count) = alpha_arr(l)//alpha_arr(m)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str2 =  strh//strg//stri//strj
    str3 =  strg//strh//strj//stri
    str4 =  strh//strg//strj//stri

    

    do s = 1,36 
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str2 .Or. disp_arr(s)==str3 .Or. disp_arr(s)==str4) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mm_mm


SUBROUTINE GetIndex_mQd_mm (g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 2) , dimension(15)::alpha_Qd 
    Character(len = 20) , dimension(90)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index
    Character(Len=1)::stri,strj,strg,strh
    Integer :: s,l,m,count,count_mQ,u,v
    
    Character(len = 10) ::str1,str2,str3,str4

    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"



    
    alpha_Qd(1) = "11"
    alpha_Qd(2) = "12"
    alpha_Qd(3) = "13"
    alpha_Qd(4) = "14"
    alpha_Qd(5) = "15"
    alpha_Qd(6) = "21"
    alpha_Qd(7) = "22"
    alpha_Qd(8) = "23"
    alpha_Qd(9) = "24"
    alpha_Qd(10) = "25"
    alpha_Qd(11) = "31"
    alpha_Qd(12) = "32"
    alpha_Qd(13) = "33"
    alpha_Qd(14) = "34"
    alpha_Qd(15) = "35"
    count=1;

    do l=1,15
        do m=1,6
            disp_arr(count) = alpha_Qd(l)//alpha_arr(m)
            count=count+1
        end do
    end do

    
    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str3 =  strg//strh//strj//stri


    do s = 1,90 
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str3 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mQd_mm


SUBROUTINE GetIndex_mm_mQd (g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 2) , dimension(15)::alpha_Qd 
    Character(len = 20) , dimension(90)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index
    Character(Len=1)::stri,strj,strg,strh
    Integer :: s,l,m,count,count_mQ,u,v
    
    Character(len = 10) ::str1,str2,str3,str4


    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"

    alpha_Qd(1) = "11"
    alpha_Qd(2) = "12"
    alpha_Qd(3) = "13"
    alpha_Qd(4) = "14"
    alpha_Qd(5) = "15"
    alpha_Qd(6) = "21"
    alpha_Qd(7) = "22"
    alpha_Qd(8) = "23"
    alpha_Qd(9) = "24"
    alpha_Qd(10) = "25"
    alpha_Qd(11) = "31"
    alpha_Qd(12) = "32"
    alpha_Qd(13) = "33"
    alpha_Qd(14) = "34"
    alpha_Qd(15) = "35"



    count=1;

    do l=1,15
        do m=1,6
            disp_arr(count) = alpha_arr(m)//alpha_Qd(l)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str2 =  strh//strg//stri//strj



    do s = 1,90 
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str2 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mm_mQd


SUBROUTINE GetIndex_mmm (h,i,j,index)

    Character(len = 20) , dimension(10)::beta_arr 
    Integer , INTENT(IN) :: h,i , j
    Character(Len=1)::stri,strj,strh
    Integer, INTENT(INOUT) :: index

    Integer :: s
    
    Character(len = 3) ::str1,str2,str3,str4,str5,str6



    beta_arr(1) = "111"
    beta_arr(2) = "112"
    beta_arr(3) = "113"
    beta_arr(4) = "122"
    beta_arr(5) = "123"
    beta_arr(6) = "133"
    beta_arr(7) = "222"
    beta_arr(8) = "223"
    beta_arr(9) = "233"
    beta_arr(10) = "333"

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(h,strh)

    str1 =  strh//stri//strj
    str2 =  strh//strj//stri
    str3 =  stri//strh//strj
    str4 =  stri//strj//strh
    str5 =  strj//strh//stri
    str6 =  strj//stri//strh


    index=1
    do s = 1,10
        
        if (beta_arr(s)==str1 .Or. beta_arr(s)==str2 .Or. beta_arr(s)==str3 .Or. beta_arr(s)==str4 &
            .Or. beta_arr(s)==str5 .Or. beta_arr(s)==str6) Then
            index = s;
        end if 
    end do 

END SUBROUTINE GetIndex_mmm 


SUBROUTINE GetIndex_mmQd (h,i,j,index)

    Character(len = 20) , dimension(30)::beta_arr 
    Integer , INTENT(IN) :: h,i , j
    Integer, INTENT(INOUT) :: index

    Integer :: s
    Character(Len=1)::stri,strj,strh
    Character(len = 3) ::str1,str3



    beta_arr(1) = "111"
    beta_arr(2) = "112"
    beta_arr(3) = "113"
    beta_arr(4) = "114"
    beta_arr(5) = "115"

    beta_arr(6) = "121"
    beta_arr(7) = "122"
    beta_arr(8) = "123"
    beta_arr(9) = "124"
    beta_arr(10) = "125"
    
    beta_arr(11) = "131"
    beta_arr(12) = "132"
    beta_arr(13) = "133"
    beta_arr(14) = "134"
    beta_arr(15) = "135"

    beta_arr(16) = "221"
    beta_arr(17) = "222"
    beta_arr(18) = "223"
    beta_arr(19) = "224"
    beta_arr(20) = "225"

    beta_arr(21) = "231"
    beta_arr(22) = "232"
    beta_arr(23) = "233"
    beta_arr(24) = "234"
    beta_arr(25) = "235"

    beta_arr(26) = "331"
    beta_arr(27) = "332"
    beta_arr(28) = "333"
    beta_arr(29) = "334"
    beta_arr(30) = "335"

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(h,strh)

    str1 =  strh//stri//strj
    str3 =  stri//strh//strj


  

    index=1
    do s = 1,30
        
        if (beta_arr(s)==str1 .Or. beta_arr(s)==str3 ) Then
            index = s;
        end if 
    end do 

END SUBROUTINE GetIndex_mmQd


SUBROUTINE GetIndex_mQd_mQd (g,h,i,j,index)

    Character(len = 2) , dimension(15)::alpha_Qd
    Character(len = 20) , dimension(225)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index
    Character(Len=1)::stri,strj,strh,strg
    Integer :: s,l,m,count,count_mQ,u,v
    
    Character(len = 10) ::str1

    alpha_Qd(1) = "11"
    alpha_Qd(2) = "12"
    alpha_Qd(3) = "13"
    alpha_Qd(4) = "14"
    alpha_Qd(5) = "15"
    alpha_Qd(6) = "21"
    alpha_Qd(7) = "22"
    alpha_Qd(8) = "23"
    alpha_Qd(9) = "24"
    alpha_Qd(10) = "25"
    alpha_Qd(11) = "31"
    alpha_Qd(12) = "32"
    alpha_Qd(13) = "33"
    alpha_Qd(14) = "34"
    alpha_Qd(15) = "35"






    alpha_Qd_B = alpha_Qd_A
    count=1;

    do l=1,15
        do m=1,15
            disp_arr(count) = alpha_Qd(l)//alpha_Qd(m)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj


    do s = 1,225 
        
        if (disp_arr(s)==str1 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mQd_mQd


SUBROUTINE GetIndex_mO_mm (g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 2) , dimension(21)::alpha_mO 
    Character(len = 20) , dimension(126)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index
    Character(Len=1)::stri,strj,strh,strg
    Integer :: s,l,m,count,count_mO,u,v
    
    Character(len = 10) ::str1,str2,str3,str4

    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"



    alpha_mO(1) = "11"
    alpha_mO(2) = "12"
    alpha_mO(3) = "13"
    alpha_mO(4) = "14"
    alpha_mO(5) = "15"
    alpha_mO(6) = "16"
    alpha_mO(7) = "17"
    alpha_mO(8) = "21"
    alpha_mO(9) = "22"
    alpha_mO(10) = "23"
    alpha_mO(11) = "24"
    alpha_mO(12) = "25"
    alpha_mO(13) = "26"
    alpha_mO(14) = "27"
    alpha_mO(15) = "31"
    alpha_mO(16) = "32"
    alpha_mO(17) = "33"
    alpha_mO(18) = "34"
    alpha_mO(19) = "35"
    alpha_mO(20) = "36"
    alpha_mO(21) = "37"


    count=1;

    do l=1,21
        do m=1,6
            disp_arr(count) = alpha_mO(l)//alpha_arr(m)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str3 =  strg//strh//strj//stri



    do s = 1,126 
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str3 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mO_mm


SUBROUTINE GetIndex_mm_mO (g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 2) , dimension(21)::alpha_mO 
    Character(len = 20) , dimension(126)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index
    Character(Len=1)::stri,strj,strh,strg
    Integer :: s,l,m,count,count_mO,u,v
    
    Character(len = 10) ::str1,str2


    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"

    alpha_mO(1) = "11"
    alpha_mO(2) = "12"
    alpha_mO(3) = "13"
    alpha_mO(4) = "14"
    alpha_mO(5) = "15"
    alpha_mO(6) = "16"
    alpha_mO(7) = "17"
    alpha_mO(8) = "21"
    alpha_mO(9) = "22"
    alpha_mO(10) = "23"
    alpha_mO(11) = "24"
    alpha_mO(12) = "25"
    alpha_mO(13) = "26"
    alpha_mO(14) = "27"
    alpha_mO(15) = "31"
    alpha_mO(16) = "32"
    alpha_mO(17) = "33"
    alpha_mO(18) = "34"
    alpha_mO(19) = "35"
    alpha_mO(20) = "36"
    alpha_mO(21) = "37"


    count=1;

    do l=1,21
        do m=1,6
            disp_arr(count) = alpha_arr(m)//alpha_mO(l)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str2 =  strh//strg//stri//strj


  

    do s = 1,126
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str2 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mm_mO


SUBROUTINE GetIndex_QQ_mm (g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 2) , dimension(15)::alpha_QQ 
    Character(len = 20) , dimension(90)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index

    Integer :: s,l,m,count
    Character(Len=1)::stri,strj,strh,strg
    Character(len = 10) ::str1,str2,str3,str4

    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"


    alpha_QQ(1) = "11"
    alpha_QQ(2) = "12"
    alpha_QQ(3) = "13"
    alpha_QQ(4) = "14"
    alpha_QQ(5) = "15"
    alpha_QQ(6) = "22"
    alpha_QQ(7) = "23"
    alpha_QQ(8) = "24"
    alpha_QQ(9) = "25"
    alpha_QQ(10) = "33"
    alpha_QQ(11) = "34"
    alpha_QQ(12) = "35"
    alpha_QQ(13) = "44"
    alpha_QQ(14) = "45"
    alpha_QQ(15) = "55"





    count=1;

    do l=1,15
        do m=1,6
            disp_arr(count) = alpha_QQ(l)//alpha_arr(m)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str2 =  strh//strg//stri//strj
    str3 =  strg//strh//strj//stri
    str4 =  strh//strg//strj//stri

    

    do s = 1,90 
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str3 .Or. disp_arr(s)==str2 .Or. disp_arr(s)==str4 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_QQ_mm


SUBROUTINE GetIndex_mm_QQ(g,h,i,j,index)

    Character(len = 2) , dimension(6)::alpha_arr 
    Character(len = 2) , dimension(15)::alpha_QQ 
    Character(len = 20) , dimension(90)::disp_arr 
    Integer , INTENT(IN) :: i,j,g,h
    Integer, INTENT(INOUT) :: index

    Integer :: s,l,m,count
    Character(Len=1)::stri,strj,strh,strg
    Character(len = 10) ::str1,str2,str3,str4

    alpha_arr(1) = "11"
    alpha_arr(2) = "12"
    alpha_arr(3) = "13"
    alpha_arr(4) = "22"
    alpha_arr(5) = "23"
    alpha_arr(6) = "33"


    alpha_QQ(1) = "11"
    alpha_QQ(2) = "12"
    alpha_QQ(3) = "13"
    alpha_QQ(4) = "14"
    alpha_QQ(5) = "15"
    alpha_QQ(6) = "22"
    alpha_QQ(7) = "23"
    alpha_QQ(8) = "24"
    alpha_QQ(9) = "25"
    alpha_QQ(10) = "33"
    alpha_QQ(11) = "34"
    alpha_QQ(12) = "35"
    alpha_QQ(13) = "44"
    alpha_QQ(14) = "45"
    alpha_QQ(15) = "55"





    count=1;

    do l=1,15
        do m=1,6
            disp_arr(count) = alpha_arr(m)//alpha_QQ(l)
            count=count+1
        end do
    end do

    call GetStr(i,stri)
    call GetStr(j,strj)
    call GetStr(g,strg)
    call GetStr(h,strh)

    str1 =  strg//strh//stri//strj
    str2 =  strh//strg//stri//strj
    str3 =  strg//strh//strj//stri
    str4 =  strh//strg//strj//stri
    


    do s = 1,90 
        
        if (disp_arr(s)==str1 .Or. disp_arr(s)==str3 .Or. disp_arr(s)==str2 .Or. disp_arr(s)==str4 ) Then

            index = s
        end if 
    end do 



    
END SUBROUTINE GetIndex_mm_QQ


SUBROUTINE GetStr(i,stri)

    Integer , INTENT(IN) :: i
    Character(len = 1), INTENT(INOUT) ::stri

    if (i==0)then
        stri ="0" 
        Return
    elseif (i==1)then
        stri ="1" 
        Return
    elseif (i==2)then
        stri ="2" 
        Return
    elseif (i==3)then
        stri ="3" 
        Return
    elseif (i==4)then
        stri ="4" 
        Return                
    elseif (i==5)then
        stri ="5" 
        Return
    elseif (i==6)then
        stri ="6" 
        Return
    elseif (i==7)then
        stri ="7" 
        Return
    elseif (i==8)then
        stri ="8" 
        Return
    elseif (i==9)then
        stri ="9" 
        Return
    endif

    
END SUBROUTINE GetStr



!********************************************************
module Tensors_constant  
   implicit none 
   real*8 :: T_Tensor(40,40,40,40),T_Changed(40,40,40,40)

  contains

  subroutine init_Tensors()  
     implicit none 
     integer :: i,j,k,l

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


!********************************************************
SUBROUTINE T_l0(Ar,Br,C,q,M2,IND,lk , result)
        IMPLICIT NONE

        INTEGER ::  i
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind,lk
        real*8 :: q,t_k1,t_k2
        real*8 , dimension(2*lk+1):: M2
        real*8:: eps=EPSILON(result)



        result = 0d0 

        do i = 0,lk
        

            if (ind==0) then

                if (i==0) then
        
                    if (DABS(q)>eps .and. DABS(M2(1))>eps) Then
                        Call T_lk(Ar,Br,C,lk,0,"0",0,0,"0",t_k1)
                        
                        result = result + q*M2(1)*t_k1 
                    end if 
                else

                    if (DABS(q)>eps .and. DABS(M2(2*i))>eps) Then
                        Call T_lk(Ar,Br,C,lk,i,"c",0,0,"0",t_k1)
                        result = result + q*M2(2*i)*t_k1 
                    end if 
                    if ( DABS(q)>eps .and. DABS(M2(2*i+1))>eps) Then
                        Call T_lk(Ar,Br,C,lk,i,"s",0,0,"0",t_k2)
                        result = result + q*M2(2*i+1)*t_k2
                    end if 
                    
                end if
        ELSE
            if (i==0) then

                if (DABS(q)>eps .and. DABS(M2(1))>eps) Then
                    Call T_lk(Ar,Br,C,0,0,"0",lk,0,"0",t_k1)
                    result = result + q*M2(1)*t_k1 
                end if 
            else

                    if (DABS(q)>eps .and. DABS(M2(2*i))>eps) Then
                        Call T_lk(Ar,Br,C,0,0,"0",lk,i,"c",t_k1)
                        result = result + q*M2(2*i)*t_k1 
                    end if 
                    if ( DABS(q)>eps .and. DABS(M2(2*i+1))>eps) Then
                        Call T_lk(Ar,Br,C,0,0,"0",lk,i,"s",t_k2)
                        result = result + q*M2(2*i+1)*t_k2
                    end if 
                    
                end if
                
            end if 
        end do

        RETURN
END SUBROUTINE T_l0



    !********************************************************
SUBROUTINE T_ll(Ar,Br,C,M1,M2,lk1,lk2 , result)
    IMPLICIT NONE

    INTEGER ::  i,j
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
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

                        Call T_lk(Ar,Br,C,lk1,Floor((i*1d0)/2d0),cp1,lk2,Floor((j*1d0)/2d0),cp2,t_k1)
                        result = result + M1(i)*M2(j)*t_k1 
                    end if 
                
                end do
        end do

    RETURN
END SUBROUTINE T_ll


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
    