!> \brief TOML 实验记录写入器模块 (路线图 §15.10「数据格式」项, 与 mod_json 同构)
!> \details 无外部依赖的最小 TOML 1.0 写入器: 根表键值对 + [table] /
!>          [table.sub] 节; 字符串 (基本串转义与 JSON 一致), 整数, 实数
!>          (18 位有效数字规范格式, 与 mod_json 相同契约), 布尔, 以及
!>          一/二维数值数组 (多行数组, 每元素一行)。键-值对按写入顺序
!>          输出; 节必须先声明再写键。
!>          规范化约定 (与 Python 镜像 `pygenmod/toml_io.py` 一致,
!>          由 tests/fixtures/canonical_record.toml 固定为跨语言契约):
!>            - "key = value" 单位置格式, 根表键在文件首部;
!>            - 数组: "key = [" + 每元素一行 (2 空格缩进) + "]" ;
!>            - 空数组写 "key = []"; 二维数组为嵌套多行数组;
!>            - 数字: d.dddddddddddddddddE±ddd (与 JSON 契约相同的往返无损格式)。
!> \author LiHao
module mod_toml
    use mod_constants, only: dp
    implicit none
    private

    public :: toml_writer_t
    public :: toml_init
    public :: toml_finalize
    public :: toml_get_text
    public :: toml_begin_table
    public :: toml_begin_array_table
    public :: toml_add_string
    public :: toml_add_int
    public :: toml_add_real
    public :: toml_add_bool
    public :: toml_add_int_array
    public :: toml_add_real_array
    public :: toml_add_real_array_2d

    integer, parameter :: INIT_CAP = 4096

    !> \brief TOML 写入器 (可增长缓冲 + 当前节路径)
    type :: toml_writer_t
        character(len=:), allocatable :: buf
        integer :: pos = 0
        character(len=256) :: cur_table = ""   !< 当前节 (空 = 根表)
        logical :: table_open = .false.        !< 表头后是否已写首个键
    end type toml_writer_t

contains

    subroutine toml_init(w)
        type(toml_writer_t), intent(out) :: w
        allocate (character(len=INIT_CAP) :: w%buf)
        w%pos = 0
        w%cur_table = ""
        w%table_open = .false.
    end subroutine toml_init

    subroutine toml_finalize(w)
        type(toml_writer_t), intent(inout) :: w
        if (allocated(w%buf)) deallocate (w%buf)
        w%pos = 0
    end subroutine toml_finalize

    pure function toml_get_text(w) result(txt)
        type(toml_writer_t), intent(in) :: w
        character(len=:), allocatable :: txt
        if (w%pos > 0) then
            txt = w%buf(1:w%pos)
        else
            txt = ""
        end if
    end function toml_get_text

    subroutine ensure_capacity(w, need)
        type(toml_writer_t), intent(inout) :: w
        character(len=:), allocatable :: newbuf
        integer, intent(in) :: need
        integer :: newcap
        if (w%pos + need <= len(w%buf)) return
        newcap = max(2 * len(w%buf), w%pos + need + 256)
        allocate (character(len=newcap) :: newbuf)
        newbuf(1:w%pos) = w%buf(1:w%pos)
        call move_alloc(newbuf, w%buf)
    end subroutine ensure_capacity

    subroutine put(w, txt)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: txt
        call ensure_capacity(w, len(txt))
        w%buf(w%pos + 1:w%pos + len(txt)) = txt
        w%pos = w%pos + len(txt)
    end subroutine put

    !> \brief 声明节 [name] (name 可为点分嵌套 "a.b"); 空名恢复根表
    subroutine toml_begin_table(w, name)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: name
        if (len_trim(name) == 0) then
            w%cur_table = ""
            w%table_open = .false.
            return
        end if
        if (w%pos > 0) call put(w, achar(10))
        call put(w, "["//trim(name)//"]"//achar(10))
        w%cur_table = trim(name)
        w%table_open = .false.
    end subroutine toml_begin_table

    !> \brief 写键前缀 (自动在值前换行; 根表与节内一致)
    subroutine put_key(w, key)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        ! 换行分隔 (上一字符已是换行时不再重复, 避免表头后空行)
        if (w%pos > 0) then
            if (w%buf(w%pos:w%pos) /= achar(10)) call put(w, achar(10))
        end if
        call put(w, trim(key)//" = ")
    end subroutine put_key

    !> \brief 基本串转义 (与 mod_json 一致)
    subroutine put_escaped(w, s)
        type(toml_writer_t), intent(inout) :: w
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

    !> \brief 实数: 与 mod_json 相同的规范格式 (往返无损)
    subroutine put_real(w, val)
        type(toml_writer_t), intent(inout) :: w
        real(dp), intent(in) :: val
        character(len=32) :: num
        write (num, '(ES26.17E3)') val
        call put(w, trim(adjustl(num)))
    end subroutine put_real

    subroutine toml_add_string(w, key, val)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key, val
        call put_key(w, key)
        call put(w, '"')
        call put_escaped(w, trim(val))
        call put(w, '"')
    end subroutine toml_add_string

    subroutine toml_add_int(w, key, val)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        integer, intent(in) :: val
        character(len=32) :: num
        call put_key(w, key)
        write (num, '(I0)') val
        call put(w, trim(num))
    end subroutine toml_add_int

    subroutine toml_add_real(w, key, val)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: val
        call put_key(w, key)
        call put_real(w, val)
    end subroutine toml_add_real

    subroutine toml_add_bool(w, key, val)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        logical, intent(in) :: val
        call put_key(w, key)
        if (val) then
            call put(w, "true")
        else
            call put(w, "false")
        end if
    end subroutine toml_add_bool

    subroutine toml_add_int_array(w, key, vals)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        integer, intent(in) :: vals(:)
        character(len=32) :: num
        integer :: i
        call put_key(w, key)
        if (size(vals) == 0) then
            call put(w, "[]")
            return
        end if
        call put(w, "[")
        do i = 1, size(vals)
            if (i > 1) call put(w, ",")
            call put(w, achar(10)//"  ")
            write (num, '(I0)') vals(i)
            call put(w, trim(num))
        end do
        call put(w, achar(10)//"]")
    end subroutine toml_add_int_array

    subroutine toml_add_real_array(w, key, vals)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: vals(:)
        integer :: i
        call put_key(w, key)
        if (size(vals) == 0) then
            call put(w, "[]")
            return
        end if
        call put(w, "[")
        do i = 1, size(vals)
            if (i > 1) call put(w, ",")
            call put(w, achar(10)//"  ")
            call put_real(w, vals(i))
        end do
        call put(w, achar(10)//"]")
    end subroutine toml_add_real_array

    subroutine toml_add_real_array_2d(w, key, vals)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: vals(:, :)
        integer :: i, j, nrow, ncol
        nrow = size(vals, 1)
        ncol = size(vals, 2)
        call put_key(w, key)
        if (nrow == 0 .or. ncol == 0) then
            call put(w, "[]")
            return
        end if
        call put(w, "[")
        do i = 1, nrow
            if (i > 1) call put(w, ",")
            call put(w, achar(10)//"  [")
            do j = 1, ncol
                if (j > 1) call put(w, ",")
                call put(w, achar(10)//"    ")
                call put_real(w, vals(i, j))
            end do
            call put(w, achar(10)//"  ]")
        end do
        call put(w, achar(10)//"]")
    end subroutine toml_add_real_array_2d

    !> \brief 声明数组表 [[name]] (数组中的对象元素; 后续键属于该元素)
    subroutine toml_begin_array_table(w, name)
        type(toml_writer_t), intent(inout) :: w
        character(len=*), intent(in) :: name
        if (w%pos > 0) call put(w, achar(10))
        call put(w, "[["//trim(name)//"]]"//achar(10))
        w%cur_table = trim(name)
        w%table_open = .false.
    end subroutine toml_begin_array_table

end module mod_toml
