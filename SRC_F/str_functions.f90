module str_functions

    implicit none

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function str_cat(string1, string2, string3, string4, string5&
                              , string6, string7, string8, string9&
                              , string10, string11, string12, string13&
                              , string14, string15, verbose) result(stringTot)

        implicit none

        !INPUT
        character (len=*), intent(in) :: string1;
        character (len=*), intent(in), optional  :: string2;
        character (len=*), intent(in), optional  :: string3;
        character (len=*), intent(in), optional  :: string4;
        character (len=*), intent(in), optional  :: string5;
        character (len=*), intent(in), optional  :: string6;
        character (len=*), intent(in), optional  :: string7;
        character (len=*), intent(in), optional  :: string8;
        character (len=*), intent(in), optional  :: string9;
        character (len=*), intent(in), optional  :: string10;
        character (len=*), intent(in), optional  :: string11;
        character (len=*), intent(in), optional  :: string12;
        character (len=*), intent(in), optional  :: string13;
        character (len=*), intent(in), optional  :: string14;
        character (len=*), intent(in), optional  :: string15;
        logical, intent(in), optional  :: verbose;

        !OUTPUT
        character (len=200) :: stringTot;

        !LOCAL
        integer :: i
        logical :: effecVerb

        effecVerb = .false.
        if(present(verbose)) then
            if(verbose) effecVerb = .true.
        end if

        !write(*,*) "WRITE Flag string_join"
        stringTot = ""
        stringTot = trim(adjustL(string1))
        stringTot = adjustL(stringTot)
        !write(*,*) "string1 = ", string1
        !write(*,*) "string2 = ", string2
        !write(*,*) "string3 = ", string3

        if(effecVerb) write(*,*) "string_join_many INSIDE 1 = "

        do i =1, 200
            if(effecVerb) write(*,*) i, "=", stringTot(i:i)
        end do

        if(effecVerb) write(*,*) "len(trim(stringTot)) = ", len(trim(stringTot))

        if(present(string2)) stringTot = str_cat_sub(stringTot, string2, effecVerb)
        if(effecVerb) write(*,*) "string_join_many INSIDE 2 = "

        do i =1, 200
            if(effecVerb) write(*,*) i, "=", stringTot(i:i)
        end do
        if(present(string3)) stringTot = str_cat_sub(stringTot, string3, effecVerb)

        if(effecVerb) write(*,*) "string_join_many INSIDE 3 = "

        do i =1, 200
            if(effecVerb) write(*,*) i, "=", stringTot(i:i)
        end do

        if(present(string4))  stringTot = str_cat_sub(stringTot, string4, effecVerb)
        if(present(string5))  stringTot = str_cat_sub(stringTot, string5, effecVerb)
        if(present(string6))  stringTot = str_cat_sub(stringTot, string6, effecVerb)
        if(present(string7))  stringTot = str_cat_sub(stringTot, string7, effecVerb)
        if(present(string8))  stringTot = str_cat_sub(stringTot, string8, effecVerb)
        if(present(string9))  stringTot = str_cat_sub(stringTot, string9, effecVerb)
        if(present(string10)) stringTot = str_cat_sub(stringTot, string10, effecVerb)
        if(present(string11)) stringTot = str_cat_sub(stringTot, string11, effecVerb)
        if(present(string12)) stringTot = str_cat_sub(stringTot, string12, effecVerb)
        if(present(string13)) stringTot = str_cat_sub(stringTot, string13, effecVerb)
        if(present(string14)) stringTot = str_cat_sub(stringTot, string14, effecVerb)
        if(present(string15)) stringTot = str_cat_sub(stringTot, string15, effecVerb)

        stringTot = adjustL(stringTot)

        if(effecVerb) write(*,*) " stringTot = ", stringTot

    end function str_cat

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function str_cat_sub(string1, string2, verbose) result(stringTot)

        implicit none

        !INPUT
        character (len=*)  , intent(in) :: string1, string2;
        logical, intent(in), optional  :: verbose;

        !OUTPUT
        character (len=200) :: stringTot;

        !LOCAL
        logical :: effecVerb

        effecVerb = .false.
        if(present(verbose)) then
            if(verbose) effecVerb = .true.
        end if

        !write(*,*) "WRITE Flag string_join"

        stringTot = trim(adjustL(string1))//trim(adjustL(string2))
        stringTot = adjustL(stringTot)

        if(effecVerb) write(*,*) "Inside string_join: stringTot = ", stringTot

    end function str_cat_sub

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function strnum_cat(string, number) result(stringTot)

        implicit none

        !INPUT
        character (len=*)  , intent(in) :: string
        integer            , intent(in) :: number

        !OUTPUT
        character (len=100) :: stringTot;

        !LOCAL
        character (len=30)  :: nmbString

        !write(*,*) "WRITE Flag stringNumb_join"

        write(nmbString, fmt='(I8)') number
        stringTot = str_cat(string, nmbString)

        !write(*,*) "WRITE Flag 2 stringNumb_join"

    end function strnum_cat

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function num2str(number, nCharacters) result(stringTot)

        implicit none

        !INPUT
        integer, intent(in) :: number
        integer, intent(in), optional :: nCharacters

        !OUTPUT
        character (len=30) :: stringTot;

        !LOCAL
        character (len=30)  :: nmbString
        integer :: n, i

        write(nmbString, fmt='(I30)') number
        nmbString = adjustL(nmbString)

        if(present(nCharacters)) then
            !write(*,*) "nmbString = ", nmbString
            n = len(trim(nmbString))
            if(n > nCharacters) stop("Inside numb2String the number of characters is to little to represent tis number")

            do i = 1, len(stringTot)
                if(i<=nCharacters) then
                    stringTot(i:i) = "0"
                else
                    stringTot(i:i) = " "
                end if

                !write(*,*) "stringTot = ", stringTot
            end do

            stringTot(nCharacters-n+1:nCharacters) = trim(nmbString)
            stringTot = adjustL(stringTot)

            !write(*,*) "nmbString = ", nmbString
            !write(*,*) "stringTot = ", stringTot
        else
            stringTot = nmbString
        end if

    end function num2str

end module str_functions
