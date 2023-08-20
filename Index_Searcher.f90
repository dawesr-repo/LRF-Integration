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


    
    ! count =1;
    ! do u=1,3
    !     do v=1,5
    !         write(alpha_Qd(count) , '(i1.1, i1.1)') u,v
    !         count=count+1
    !     end do
    ! end do



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


    ! count =1;
    ! do u=1,3
    !     do v=1,7
    !         write(alpha_mO(count) , '(i1.1, i1.1)') u,v
    !         count=count+1
    !     end do
    ! end do



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




    ! count_mQ =1;
    ! do u=1,3
    !     do v=1,5
    !         write(alpha_Qd(count_mQ) , '(i1.1, i1.1)') u,v
    !         count_mQ=count_mQ+1
    !     end do
    ! end do

    
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
    ! write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i, j
    ! write(str3, '(i1.1, i1.1,i1.1, i1.1)') g,h,j, i

    

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


    ! count_mQ =1;
    
    !     do u=1,3
    !         do v=1,5
        
    !         write(alpha_Qd(count_mQ) , '(i1.1, i1.1)') u,v
    !         count_mQ=count_mQ+1
    !     end do
    ! end do

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

    ! write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i, j
    ! write(str2, '(i1.1, i1.1,i1.1, i1.1)') h,g,i, j

  

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
    ! write(str1, '(i1.1, i1.1, i1.1)') h,i, j
    ! write(str2, '(i1.1, i1.1, i1.1)') h,j, i
    ! write(str3, '(i1.1, i1.1, i1.1)') i, h,j
    ! write(str4, '(i1.1, i1.1, i1.1)') j, h,i
    ! write(str5, '(i1.1, i1.1, i1.1)') i, j,h
    ! write(str6, '(i1.1, i1.1, i1.1)') j, i,h
  

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


    ! count_mQ =1;
    
    !     do u=1,3
    !         do v=1,5
        
    !         write(alpha_Qd_A(count_mQ) , '(i1.1, i1.1)') u,v
    !         count_mQ=count_mQ+1
    !     end do
    ! end do




    alpha_Qd_B = alpha_Qd_A
    count=1;

    do l=1,15
        do m=1,15
            disp_arr(count) = alpha_Qd(l)//alpha_Qd(m)
            count=count+1
        end do
    end do

  !  write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i, j
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

    ! count_mO =1;
    ! do u=1,3
    !     do v=1,7
    !         write(alpha_mO(count_mO) , '(i1.1, i1.1)') u,v
    !         count_mO=count_mO+1
    !     end do
    ! end do

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

    ! write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i, j
    ! write(str3, '(i1.1, i1.1,i1.1, i1.1)') g,h,j, i

    

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

    ! count_mO =1;
    
    !     do u=1,3
    !         do v=1,7
        
    !         write(alpha_mO(count_mO) , '(i1.1, i1.1)') u,v
    !         count_mO=count_mO+1
    !     end do
    ! end do

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

    ! write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i, j
    ! write(str2, '(i1.1, i1.1,i1.1, i1.1)') h,g,i, j

  

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

    ! write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i,j
    ! write(str3, '(i1.1, i1.1,i1.1, i1.1)') g,h,j,i

    ! write(str2, '(i1.1, i1.1,i1.1, i1.1)') h,g,i,j
    ! write(str4, '(i1.1, i1.1,i1.1, i1.1)') h,g,j,i

    

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
    

    ! write(str1, '(i1.1, i1.1,i1.1, i1.1)') g,h,i,j
    ! write(str3, '(i1.1, i1.1,i1.1, i1.1)') g,h,j,i

    ! write(str2, '(i1.1, i1.1,i1.1, i1.1)') h,g,i,j
    ! write(str4, '(i1.1, i1.1,i1.1, i1.1)') h,g,j,i

    

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