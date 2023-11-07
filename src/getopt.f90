


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

subroutine read_options(input_fname)
  use options_module
  implicit none
  character :: okey
  integer :: argc
  character(len=*), intent(out) :: input_fname

    argc = command_argument_count()

    if (argc < 1) then
        call print_usage()
        stop
    end if


    do
        okey = getopt('ham:pQr:s:uvw:X')
        if(okey == '>') exit
        if(okey == '!') then
            write(*,*) 'unknown option: ', trim(optarg)
            stop
        end if
        if(okey == 'a') then
            is_adjacencyfile = .true.
        end if
        if(okey == 'm') then
            read(optarg, *) maxrecords
        end if

        if(okey == 'p') then
            print_bondlevel = .true.
        end if

        if(okey == 'Q') then
            print_XML = .true.
        end if

        if(okey == 'r') then
            read(optarg, '(a)') read_cache_fname
            has_read_cache_file = .true. 
        end if

        if(okey == 's') then
            read(optarg, *) nbuckets
        end if

        if(okey == 'v') then
            verbose = .true.
        end if

        if(okey == 'u') then
            unsorted_geometry = .true.
        end if

        if(okey == 'w') then
            read(optarg, '(a)') write_cache_fname
            has_write_cache_file = .true. 
        end if

        if(okey == 'X') then
            read_connection_table = .true.
            unsorted_geometry = .true. ! connection table would not work if geometry is sorted
        end if


        if(okey == '.') then
            input_fname = optarg
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
    write(*, '(1x,10a)') "    ", "-a", "                ",  &
                                "Specifies that the input file contains an adjacency matrix instead of XYZ format"
    write(*, '(1x,10a)') "    ", "-m number", "         ",  "Sets the maximum {number} of structures in the cache database"
    write(*, '(1x,10a)') "    ", "-p", "                ",  "Prints intermediate bond-level structures"
    write(*, '(1x,10a)') "    ", "-Q", "                ",  "Print the ZZ polynomial in XML format"
    write(*, '(1x,10a)') "    ", "-r file", "           ",  "Reads cached structures from a {file}"
    write(*, '(1x,10a)') "    ", "-s number", "         ",  "Sets the {number} of buckets in the cache database"
    write(*, '(1x,10a)') "    ", "-u", "                ",  "Uses unmodified input XYZ geometry (sorted by default)"
    write(*, '(1x,10a)') "    ", "-v", "                ",  "Enables verbose printing"
    write(*, '(1x,10a)') "    ", "-w file", "           ",  "Writes cached structures to a {file}"
    write(*, '(1x,10a)') "    ", "-X", "                ",  "Reads connection table from the bottom of the XYZ file"
    write(*, '(1x,10a)') "    ", "-h", "                ",  "Displays this help message"
end subroutine


end module

