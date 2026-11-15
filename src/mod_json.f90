!> \brief JSON 实验记录写入器模块 (路线图 §15.10「数据格式」项)
!> \details 无外部依赖的最小 JSON (RFC 8259) 写入器, 以可增长字符串缓冲累积文档,
!>          支持对象/数组嵌套, 字符串/整数/实数/布尔/null 标量 (键可选:
!>          无键即数组元素), 以及一维/二维数值数组。字符串按规范转义
!>          (\" \\ 与 \uXXXX 控制字符); 实数以 17 位有效数字写为
!>          1.00000000000000000E+000 形式 (17 位小数 + 3 位指数), 双精度往返无损。
!>
!>          规范化格式约定 (与 Python 镜像 `pygenmod/json_io.py` 一致,
!>          由 tests/fixtures/canonical_record.json 作为跨语言契约固定):
!>            - 每层缩进 2 空格, 对象/数组元素一行一个;
!>            - 空对象 {} 与空数组 [] 单行;
!>            - 元素间逗号 + 换行, 最后一个元素后无逗号;
!>            - 根容器从第 0 列开始。
!> \author LiHao
module mod_json
    use mod_constants, only: dp
    implicit none
    private

    public :: json_writer_t
    public :: json_init
    public :: json_finalize
    public :: json_get_text
    public :: json_begin_object
    public :: json_end_object
    public :: json_begin_array
    public :: json_end_array
    public :: json_add_string
    public :: json_add_int
    public :: json_add_real
    public :: json_add_bool
    public :: json_add_null
    public :: json_add_int_array
    public :: json_add_real_array
    public :: json_add_real_array_2d

    integer, parameter :: INIT_CAP = 4096
    integer, parameter :: MAX_DEPTH = 64

    !> \brief JSON 写入器 (可增长缓冲 + 层级记账)
    type :: json_writer_t
        character(len=:), allocatable :: buf        !< 文档缓冲
        integer :: pos = 0                          !< 当前写入位置
        integer :: depth = 0                        !< 当前打开的容器数
        logical :: after_key = .false.              !< 刚写入键, 值接同行
        logical :: has_elems(0:MAX_DEPTH) = .false.    !< 各层是否已有元素 (深度 0 起)
    end type json_writer_t

contains

    !> \brief 初始化写入器
    subroutine json_init(w)
        type(json_writer_t), intent(out) :: w
        allocate (character(len=INIT_CAP) :: w%buf)
        w%pos = 0
        w%depth = 0
        w%after_key = .false.
        w%has_elems = .false.
    end subroutine json_init

    !> \brief 释放写入器缓冲
    subroutine json_finalize(w)
        type(json_writer_t), intent(inout) :: w
        if (allocated(w%buf)) deallocate (w%buf)
        w%pos = 0
        w%depth = 0
    end subroutine json_finalize

    !> \brief 取回完整文档文本
    pure function json_get_text(w) result(txt)
        type(json_writer_t), intent(in) :: w
        character(len=:), allocatable :: txt
        if (w%pos > 0) then
            txt = w%buf(1:w%pos)
        else
            txt = ""
        end if
    end function json_get_text

    !> \brief 缓冲扩容 (倍增)
    subroutine ensure_capacity(w, need)
        type(json_writer_t), intent(inout) :: w
        character(len=:), allocatable :: newbuf
        integer, intent(in) :: need
        integer :: newcap
        if (w%pos + need <= len(w%buf)) return
        newcap = max(2 * len(w%buf), w%pos + need + 256)
        allocate (character(len=newcap) :: newbuf)
        newbuf(1:w%pos) = w%buf(1:w%pos)
        call move_alloc(newbuf, w%buf)
    end subroutine ensure_capacity

    !> \brief 写入原始片段
    subroutine put(w, txt)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: txt
        call ensure_capacity(w, len(txt))
        w%buf(w%pos + 1:w%pos + len(txt)) = txt
        w%pos = w%pos + len(txt)
    end subroutine put

    !> \brief 换行 + 缩进 (level 层, 每层 2 空格)
    subroutine put_indent(w, level)
        type(json_writer_t), intent(inout) :: w
        integer, intent(in) :: level
        integer :: i
        call put(w, achar(10))
        do i = 1, 2 * level
            call put(w, " ")
        end do
    end subroutine put_indent

    !> \brief 值写入前的逗号/缩进处理
    !> \details 跟在键之后的值与键同行 (after_key 处理);
    !>          文档根位置不缩进; 其余按当前层级 (depth) 处理。
    subroutine pre_value(w)
        type(json_writer_t), intent(inout) :: w
        if (w%after_key) then
            w%after_key = .false.
            return
        end if
        if (w%pos == 0) return
        if (w%has_elems(w%depth)) call put(w, ",")
        call put_indent(w, w%depth)
    end subroutine pre_value

    !> \brief 值写入完成后的层级记账
    subroutine post_value(w)
        type(json_writer_t), intent(inout) :: w
        w%has_elems(w%depth) = .true.
    end subroutine post_value

    !> \brief 写入键 (含引号与转义)
    subroutine put_key(w, key)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        call put(w, '"')
        call put_escaped(w, key)
        call put(w, '": ')
    end subroutine put_key

    !> \brief 写入 JSON 转义字符串内容 (不含外层引号)
    subroutine put_escaped(w, s)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: s
        integer :: i, code
        character(len=1) :: c
        character(len=6) :: esc
        do i = 1, len(s)
            c = s(i:i)
            code = iachar(c)
            select case (c)
            case ('"')
                call put(w, '\"')
            case ('\')
                call put(w, '\\')
            case (achar(8))
                call put(w, '\b')
            case (achar(9))
                call put(w, '\t')
            case (achar(10))
                call put(w, '\n')
            case (achar(12))
                call put(w, '\f')
            case (achar(13))
                call put(w, '\r')
            case default
                if (code < 32) then
                    write (esc, '("\u00", I2.2)') code
                    call put(w, esc)
                else
                    call put(w, c)
                end if
            end select
        end do
    end subroutine put_escaped

    !> \brief 以规范格式写实数 (1.23456789012345670E+000:
    !>        17 位小数 + 恰 3 位指数, 与 Python `f"{v:.17E}"` 规范化一致)
    subroutine put_real(w, val)
        type(json_writer_t), intent(inout) :: w
        real(dp), intent(in) :: val
        character(len=32) :: num
        write (num, '(ES26.17E3)') val
        call put(w, trim(adjustl(num)))
    end subroutine put_real

    !> \brief 开始对象 (键可选; 无键为数组元素)
    subroutine json_begin_object(w, key)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        call put(w, "{")
        w%depth = w%depth + 1
        w%has_elems(w%depth) = .false.
    end subroutine json_begin_object

    !> \brief 结束对象
    subroutine json_end_object(w)
        type(json_writer_t), intent(inout) :: w
        logical :: had
        had = w%has_elems(w%depth)
        w%depth = w%depth - 1
        if (had) call put_indent(w, w%depth)
        call put(w, "}")
        call post_value(w)
    end subroutine json_end_object

    !> \brief 开始数组 (键可选; 无键为数组元素)
    subroutine json_begin_array(w, key)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        call put(w, "[")
        w%depth = w%depth + 1
        w%has_elems(w%depth) = .false.
    end subroutine json_begin_array

    !> \brief 结束数组
    subroutine json_end_array(w)
        type(json_writer_t), intent(inout) :: w
        logical :: had
        had = w%has_elems(w%depth)
        w%depth = w%depth - 1
        if (had) call put_indent(w, w%depth)
        call put(w, "]")
        call post_value(w)
    end subroutine json_end_array

    !> \brief 字符串值 (键可选)
    subroutine json_add_string(w, key, val)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        character(len=*), intent(in) :: val
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        call put(w, '"')
        call put_escaped(w, val)
        call put(w, '"')
        call post_value(w)
    end subroutine json_add_string

    !> \brief 整数值 (键可选)
    subroutine json_add_int(w, key, val)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        integer, intent(in) :: val
        character(len=32) :: num
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        write (num, '(I0)') val
        call put(w, trim(num))
        call post_value(w)
    end subroutine json_add_int

    !> \brief 实数值 (键可选)
    subroutine json_add_real(w, key, val)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        real(dp), intent(in) :: val
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        call put_real(w, val)
        call post_value(w)
    end subroutine json_add_real

    !> \brief 布尔值 (键可选)
    subroutine json_add_bool(w, key, val)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        logical, intent(in) :: val
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        if (val) then
            call put(w, "true")
        else
            call put(w, "false")
        end if
        call post_value(w)
    end subroutine json_add_bool

    !> \brief null 值 (键可选)
    subroutine json_add_null(w, key)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        call put(w, "null")
        call post_value(w)
    end subroutine json_add_null

    !> \brief 整数一维数组 (键可选)
    subroutine json_add_int_array(w, key, vals)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        integer, intent(in) :: vals(:)
        character(len=32) :: num
        integer :: i
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        if (size(vals) == 0) then
            call put(w, "[]")
        else
            call put(w, "[")
            do i = 1, size(vals)
                if (i > 1) call put(w, ",")
                call put_indent(w, w%depth + 1)
                write (num, '(I0)') vals(i)
                call put(w, trim(num))
            end do
            call put_indent(w, w%depth)
            call put(w, "]")
        end if
        call post_value(w)
    end subroutine json_add_int_array

    !> \brief 实数一维数组 (键可选)
    subroutine json_add_real_array(w, key, vals)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        real(dp), intent(in) :: vals(:)
        integer :: i
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        if (size(vals) == 0) then
            call put(w, "[]")
        else
            call put(w, "[")
            do i = 1, size(vals)
                if (i > 1) call put(w, ",")
                call put_indent(w, w%depth + 1)
                call put_real(w, vals(i))
            end do
            call put_indent(w, w%depth)
            call put(w, "]")
        end if
        call post_value(w)
    end subroutine json_add_real_array

    !> \brief 实数二维数组 (键可选; 外层数组每元素为一行)
    subroutine json_add_real_array_2d(w, key, vals)
        type(json_writer_t), intent(inout) :: w
        character(len=*), intent(in), optional :: key
        real(dp), intent(in) :: vals(:, :)
        integer :: i, j, nrow, ncol
        nrow = size(vals, 1)
        ncol = size(vals, 2)
        call pre_value(w)
        if (present(key)) call put_key(w, key)
        if (nrow == 0 .or. ncol == 0) then
            call put(w, "[]")
        else
            call put(w, "[")
            do i = 1, nrow
                if (i > 1) call put(w, ",")
                call put_indent(w, w%depth + 1)
                call put(w, "[")
                do j = 1, ncol
                    if (j > 1) call put(w, ",")
                    call put_indent(w, w%depth + 2)
                    call put_real(w, vals(i, j))
                end do
                call put_indent(w, w%depth + 1)
                call put(w, "]")
            end do
            call put_indent(w, w%depth)
            call put(w, "]")
        end if
        call post_value(w)
    end subroutine json_add_real_array_2d

end module mod_json
