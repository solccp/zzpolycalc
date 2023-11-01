


module getopt_module
    implicit none
    public
    save
    character(len=500) :: optarg
    integer :: optind = 1

    integer, private :: optstr_ind
contains
    character function getopt(optstr)
        character(len=*), intent(in) :: optstr

        integer :: argc
        character(len=500) :: arg
        character :: okey
        integer :: found
        integer :: optstr_len

        optstr_len = len(trim(optstr))

        argc = command_argument_count()
        if( optind > argc ) then
            getopt = '>'
            return
        end if
        call get_command_argument(optind, arg)
        
        if( arg(1:1) == '-') then
            okey = arg(2:2)
            found = 0
            optstr_ind = 1
            do while(optstr_ind <= optstr_len)
                if(optstr(optstr_ind:optstr_ind) == okey) then
                    found = 1
                    if( (optstr_ind + 1) <= optstr_len) then
                        if (optstr(optstr_ind+1:optstr_ind+1) == ':') then
                            optstr_ind = optstr_ind + 1
                            optind = optind + 1
                            call get_command_argument(optind, optarg)
                        end if
                    end if
                    exit
                end if
                optstr_ind = optstr_ind+1
            end do
            if(found > 0) then
                getopt = okey
            else
                getopt = '!'
                optarg = arg
            end if
        else
            getopt = '.'
            optarg = arg
        end if
        optind = optind + 1
        return
    end function

subroutine read_options()
use options_module
implicit none
character :: okey


    do
        okey = getopt('h')
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
        end if


        if (okey == 'h') then
            call print_usage()
            stop
        end if

    end do
end subroutine read_options

subroutine print_usage()
    write(*, '(1x,a)') "Usage: ZZPolyCalc [options] input"
    write(*, '(1x,a)') "Options:"
    write(*, '(1x,10a)') "    ", "-h", "                ",  "Show this message"
end subroutine


end module

