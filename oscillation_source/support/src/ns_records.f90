! Provided typed TSV writers, unchanged wire format. No scientific kernels.
module ns_records
    use, intrinsic :: iso_fortran_env, only: real64,int64,output_unit
    implicit none
    private
    public :: emit_header,emit_text,emit_int,emit_real,emit_bool
contains
    subroutine emit_header()
        write(output_unit,'(A)') 'name'//achar(9)//'type'//achar(9)//'value'
    end subroutine emit_header

    function escaped(value) result(out)
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: out
        integer :: i
        out = ''
        do i = 1, len_trim(value)
            select case (iachar(value(i:i)))
            case (9)
                out = out//achar(92)//'t'
            case (10)
                out = out//achar(92)//'n'
            case (13)
                out = out//achar(92)//'r'
            case (92)
                out = out//achar(92)//achar(92)
            case default
                out = out//value(i:i)
            end select
        end do
    end function escaped

    subroutine emit_text(name, value)
        character(len=*), intent(in) :: name, value
        write(output_unit,'(A)') name//achar(9)//'text'//achar(9)//escaped(value)
    end subroutine emit_text

    subroutine emit_int(name, value)
        character(len=*), intent(in) :: name
        integer(int64), intent(in) :: value
        character(len=64) :: buffer
        write(buffer,'(I0)') value
        write(output_unit,'(A)') name//achar(9)//'int'//achar(9)//trim(buffer)
    end subroutine emit_int

    subroutine emit_real(name, value)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: value
        character(len=64) :: buffer
        write(buffer,'(ES26.17E3)') value
        write(output_unit,'(A)') name//achar(9)//'real'//achar(9)//trim(adjustl(buffer))
    end subroutine emit_real

    subroutine emit_bool(name, value)
        character(len=*), intent(in) :: name
        logical, intent(in) :: value
        if (value) then
            write(output_unit,'(A)') name//achar(9)//'bool'//achar(9)//'true'
        else
            write(output_unit,'(A)') name//achar(9)//'bool'//achar(9)//'false'
        end if
    end subroutine emit_bool
end module ns_records
