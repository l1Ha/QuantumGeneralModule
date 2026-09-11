!> \brief 科学计算 I/O 与诊断可视化工具模块
!> \details 提供标准化多列数据导出、矩阵与时序文件保存、控制台美化横幅及计算进度监测，提升算法库易用性。
!> \author LiHao
module mod_io_utils
    use mod_constants, only: dp
    implicit none
    private

    public :: save_data_table_1d
    public :: save_data_table_2d
    public :: save_matrix_dat
    public :: print_banner
    public :: print_progress_bar

contains

    !> \brief 保存一维双列数据 (x, y) 到标准 ASCII 数据文件
    subroutine save_data_table_1d(filename, x, y, x_name, y_name, header)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: x(:), y(:)
        character(len=*), intent(in), optional :: x_name, y_name, header
        integer :: u, i, n

        n = min(size(x), size(y))
        open(newunit=u, file=trim(adjustl(filename)), status="replace", action="write")

        if (present(header)) write(u, '(A, A)') "# ", trim(header)

        if (present(x_name) .and. present(y_name)) then
            write(u, '(A, A16, A18)') "# ", trim(x_name), trim(y_name)
        else
            write(u, '(A)') "#            x                y"
        end if

        do i = 1, n
            write(u, '(2ES18.8)') x(i), y(i)
        end do
        close(u)
    end subroutine save_data_table_1d

    !> \brief 保存多列时间序列或动力学可观测量数据矩阵
    subroutine save_data_table_2d(filename, x, y_mat, col_names, header)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y_mat(:, :)
        character(len=*), intent(in), optional :: col_names(:), header
        integer :: u, i, j, n_rows, n_cols

        n_rows = min(size(x), size(y_mat, 1))
        n_cols = size(y_mat, 2)

        open(newunit=u, file=trim(adjustl(filename)), status="replace", action="write")
        if (present(header)) write(u, '(A, A)') "# ", trim(header)

        if (present(col_names)) then
            write(u, '(A, A16)', advance='no') "# ", "Coordinate/Time"
            do j = 1, min(n_cols, size(col_names))
                write(u, '(A18)', advance='no') trim(col_names(j))
            end do
            write(u, *)
        end if

        do i = 1, n_rows
            write(u, '(ES18.8)', advance='no') x(i)
            do j = 1, n_cols
                write(u, '(ES18.8)', advance='no') y_mat(i, j)
            end do
            write(u, *)
        end do
        close(u)
    end subroutine save_data_table_2d

    !> \brief 保存实数二维方阵或势能面至数据文件
    subroutine save_matrix_dat(filename, mat, header)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: mat(:, :)
        character(len=*), intent(in), optional :: header
        integer :: u, i, j, nr, nc

        nr = size(mat, 1)
        nc = size(mat, 2)
        open(newunit=u, file=trim(adjustl(filename)), status="replace", action="write")
        if (present(header)) write(u, '(A, A)') "# ", trim(header)

        do i = 1, nr
            do j = 1, nc
                write(u, '(ES16.8)', advance='no') mat(i, j)
            end do
            write(u, *)
        end do
        close(u)
    end subroutine save_matrix_dat

    !> \brief 控制台美化打印标题横幅
    subroutine print_banner(title, width)
        character(len=*), intent(in) :: title
        integer, intent(in), optional :: width
        integer :: w, pad, len_t
        character(len=120) :: bar

        w = 64
        if (present(width)) w = max(20, min(120, width))
        bar = repeat('=', w)
        len_t = len_trim(title)
        pad = max(0, (w - len_t - 2) / 2)

        print '(A)', trim(bar)
        print '(A, A, A)', repeat(' ', pad), trim(title), repeat(' ', max(0, w - pad - len_t - 2))
        print '(A)', trim(bar)
    end subroutine print_banner

    !> \brief 控制台单行更新进度条
    subroutine print_progress_bar(current, total, prefix)
        integer, intent(in) :: current, total
        character(len=*), intent(in), optional :: prefix
        real(dp) :: frac
        integer :: pos, bar_len
        character(len=32) :: pre_str
        character(len=20) :: bar_str

        if (total <= 0) return
        pre_str = "Progress:"
        if (present(prefix)) pre_str = trim(prefix)

        frac = min(1.0_dp, max(0.0_dp, real(current, dp) / real(total, dp)))
        bar_len = 20
        pos = int(frac * real(bar_len, dp))

        bar_str = repeat('=', pos)
        if (pos < bar_len) bar_str(pos+1:pos+1) = '>'
        if (pos + 1 < bar_len) bar_str(pos+2:bar_len) = ' '

        if (current == total .or. mod(current, max(1, total / 20)) == 0) then
            print '(A, A, A, I3, A)', trim(pre_str), " [", trim(bar_str), "] ", int(frac * 100.0_dp), "%"
        end if
    end subroutine print_progress_bar

end module mod_io_utils
