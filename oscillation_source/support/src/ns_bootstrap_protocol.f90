! Provided T00 transport only. No EOS, integration, mode search or physical constants.
module ns_bootstrap_protocol
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ns_records, only: emit_header,emit_text,emit_int,emit_real,emit_bool
    implicit none
    private
    public :: read_bootstrap_input, emit_header, emit_text, emit_int, emit_real, emit_bool
contains
    subroutine read_bootstrap_input(path, token, real_value, int_value, bool_value, text_value, ok)
        character(len=*), intent(in) :: path
        character(len=*), intent(out) :: token, text_value
        real(real64), intent(out) :: real_value
        integer(int64), intent(out) :: int_value
        logical, intent(out) :: bool_value, ok
        character(len=32) :: request_token
        character(len=256) :: probe_text
        real(real64) :: probe_real
        integer(int64) :: probe_int
        logical :: probe_bool
        integer :: unit, ios
        namelist /request/ request_token
        namelist /parameters/ probe_real, probe_int, probe_bool, probe_text
        token = ''; text_value = ''; real_value = 0.0_real64
        int_value = 0_int64; bool_value = .false.; ok = .false.
        request_token = ''; probe_text = ''; probe_real = 0.0_real64
        probe_int = 0_int64; probe_bool = .false.
        open(newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read(unit, nml=request, iostat=ios)
        if (ios == 0) then
            token = request_token
            read(unit, nml=parameters, iostat=ios)
        end if
        close(unit)
        if (ios /= 0) return
        if (len_trim(request_token) /= 32 .or. .not. ieee_is_finite(probe_real)) return
        real_value = probe_real; int_value = probe_int
        bool_value = probe_bool; text_value = trim(probe_text); ok = .true.
    end subroutine read_bootstrap_input

end module ns_bootstrap_protocol
