!> \brief JSON 实验记录写入器验收测试 (路线图 §15.10「数据格式」项)
!> \details 验收组:
!>          组 1: 规范文档与 tests/fixtures/canonical_record.json 逐字节一致
!>                (跨语言契约: 同一 fixture 由 Python stdlib json 独立解析);
!>          组 2: 字符串转义边界 (引号/反斜杠/制表符/控制字符 \u00XX);
!>          组 3: 数值格式往返无损 (ES24.17E2 写出 -> 读回双精度逐位相等);
!>          组 4: 空容器单行 + 嵌套结构 (记录列表/子对象).
program test_json_io
    use general_module
    implicit none

    integer :: n_tests = 0, n_passed = 0

    call print_banner("GeneralModule: JSON Record Writer Test", 60)

    ! ================================================================
    ! 组 1: 规范文档字节一致 (跨语言契约 fixture)
    ! ================================================================
    block
        type(json_writer_t) :: w
        character(len=:), allocatable :: txt, fixture
        real(dp) :: pentad(3), mat(2, 2)
        integer :: unit, ios
        character(len=4096) :: line
        logical :: err_eq

        call json_init(w)
        call json_begin_object(w)
        call json_add_string(w, "module", "mod_json")
        call json_add_string(w, "escape_test", 'quote " backslash \ tab'//achar(9)//'ok')
        call json_add_int(w, "n_grid", 16)
        call json_add_real(w, "saddle_energy_Eh", -0.153705_dp)
        call json_add_real(w, "pi", 3.141592653589793_dp)
        call json_add_bool(w, "converged", .true.)
        call json_add_bool(w, "diverged", .false.)
        call json_add_null(w, "note")
        pentad = [1.0_dp, 2.5_dp, -3.25_dp]
        call json_add_real_array(w, "levels", pentad)
        call json_add_real_array(w, "scaled", [1.2345678901234567e-10_dp])
        call json_add_int_array(w, "shape", [16, 16, 16])
        call json_add_int_array(w, "empty_int", [integer ::])
        call json_add_real_array(w, "empty_real", [real(dp) ::])
        mat(1, 1) = 0.5_dp; mat(2, 1) = -0.25_dp
        mat(1, 2) = 1.0e-3_dp; mat(2, 2) = 42.0_dp
        call json_add_real_array_2d(w, "matrix", mat)
        call json_begin_array(w, "records")
        call json_begin_object(w)
        call json_add_string(w, "tag", "first")
        call json_add_real(w, "weight", 0.75_dp)
        call json_end_object(w)
        call json_begin_object(w)
        call json_add_string(w, "tag", "second")
        call json_add_real(w, "weight", -1.5_dp)
        call json_end_object(w)
        call json_end_array(w)
        call json_begin_object(w, "meta")
        call json_add_string(w, "author", "GeneralModule")
        call json_end_object(w)
        call json_end_object(w)
        txt = json_get_text(w)
        call json_finalize(w)

        ! 读入 fixture: 依次尝试
        !   $GM_SOURCE_DIR/tests/fixtures/ (CMake/套件注入源码根)
        !   fixtures/, tests/fixtures/, ../tests/fixtures/ (相对 CWD)
        fixture = ""
        ios = 1
        block
            character(len=512) :: env_dir
            integer :: elen, est
            call get_environment_variable("GM_SOURCE_DIR", env_dir, elen, est)
            if (est == 0 .and. elen > 0) then
                open (newunit=unit, file=trim(env_dir(:elen))//"/tests/fixtures/canonical_record.json", &
                      status="old", action="read", iostat=ios)
            end if
        end block
        if (ios /= 0) then
            open (newunit=unit, file="fixtures/canonical_record.json", status="old", &
                  action="read", iostat=ios)
        end if
        if (ios /= 0) then
            open (newunit=unit, file="tests/fixtures/canonical_record.json", status="old", &
                  action="read", iostat=ios)
        end if
        if (ios /= 0) then
            open (newunit=unit, file="../tests/fixtures/canonical_record.json", status="old", &
                  action="read", iostat=ios)
        end if
        if (ios == 0) then
            do
                read (unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                if (len_trim(fixture) > 0) fixture = fixture//achar(10)
                fixture = fixture//trim(line)
            end do
            close (unit)
        end if

        err_eq = (fixture == txt)
        n_tests = n_tests + 1
        if (err_eq .and. len(txt) > 500) then
            n_passed = n_passed + 1
            print '(A, I6, A)', " [PASS] canonical document byte-identical to fixture (", len(txt), " chars)"
        else
            print '(A, L3, A, I6, A, I6)', " [FAIL] fixture mismatch: equal=", err_eq, &
                  " doc_len=", len(txt), " fixture_len=", len(fixture)
        end if
    end block

    ! ================================================================
    ! 组 2: 字符串转义边界
    ! ================================================================
    block
        type(json_writer_t) :: w
        character(len=:), allocatable :: txt
        character(len=64) :: tricky
        logical :: ok

        tricky = 'a"b\c'//achar(9)//achar(8)//'d'
        call json_init(w)
        call json_add_string(w, "k", trim(tricky))
        txt = json_get_text(w)
        call json_finalize(w)
        ! 期望: "k": "a\"b\\c\t\bd"
        ok = (index(txt, '"a\"b\\c\t\bd"') > 0)
        n_tests = n_tests + 1
        if (ok) then
            n_passed = n_passed + 1
            print '(A)', " [PASS] string escaping: quote/backslash/tab/backspace"
        else
            print '(A, A)', " [FAIL] escaping: ", txt
        end if
    end block

    ! ================================================================
    ! 组 3: 数值格式往返无损 (写 -> 读回逐位相等)
    ! ================================================================
    block
        type(json_writer_t) :: w
        character(len=:), allocatable :: txt
        real(dp) :: vals(6), back
        integer :: i, unit, ios, n_read
        logical :: mismatch
        character(len=64) :: word
        character(len=512) :: raw

        vals(1) = 0.1_dp
        vals(2) = 3.141592653589793_dp
        vals(3) = -0.153705_dp
        vals(4) = 1.0e-300_dp
        vals(5) = 1.0e308_dp
        vals(6) = 1.2345678901234567e-10_dp

        call json_init(w)
        call json_add_real_array(w, "v", vals)
        txt = json_get_text(w)
        call json_finalize(w)

        open (newunit=unit, file="/tmp/gm_json_roundtrip.txt", status="replace", action="write")
        write (unit, '(A)') txt
        close (unit)

        n_read = 0
        mismatch = .false.
        open (newunit=unit, file="/tmp/gm_json_roundtrip.txt", status="old", action="read")
        do
            read (unit, '(A)', iostat=ios) raw
            if (ios /= 0) exit
            word = adjustl(raw)
            if (len_trim(word) > 0) then
                if (word(1:1) == '-' .or. (word(1:1) >= '0' .and. word(1:1) <= '9')) then
                    ! 去掉行尾逗号并读回
                    i = len_trim(word)
                    if (word(i:i) == ',') word(i:i) = ' '
                    read (word, *, iostat=ios) back
                    if (ios == 0) then
                        n_read = n_read + 1
                        if (n_read <= 6) then
                            if (back /= vals(n_read)) mismatch = .true.
                        end if
                    end if
                end if
            end if
        end do
        close (unit)

        n_tests = n_tests + 1
        if (n_read == 6 .and. .not. mismatch) then
            n_passed = n_passed + 1
            print '(A)', " [PASS] numeric round-trip lossless (6/6 bitwise, incl. 1e-300 / 1e308)"
        else
            print '(A, I6, A, L3)', " [FAIL] numeric round-trip, read =", n_read, " mismatch=", mismatch
        end if
    end block

    ! ================================================================
    ! 组 4: 空容器与嵌套深度
    ! ================================================================
    block
        type(json_writer_t) :: w
        character(len=:), allocatable :: txt
        logical :: ok
        integer :: i

        call json_init(w)
        call json_begin_object(w)
        call json_begin_object(w, "empty_obj")
        call json_end_object(w)
        call json_begin_array(w, "empty_arr")
        call json_end_array(w)
        call json_begin_array(w, "deep")
        do i = 1, 5
            call json_begin_object(w)
            call json_add_int(w, "depth", i)
            if (i == 5) then
                call json_begin_array(w, "tail")
                call json_add_real(w, val=1.5_dp)
                call json_add_null(w)
                call json_end_array(w)
            end if
        end do
        do i = 1, 5
            call json_end_object(w)
        end do
        call json_end_array(w)
        call json_end_object(w)
        txt = json_get_text(w)
        call json_finalize(w)

        ok = (index(txt, '"empty_obj": {}') > 0) .and. (index(txt, '"empty_arr": []') > 0) &
             .and. (index(txt, '"depth": 5') > 0) .and. (index(txt, '"tail": [') > 0)
        ! 5 层嵌套对象各缩进递增
        ok = ok .and. (index(txt, achar(10)//repeat(" ", 14)//'"depth": 5') > 0)
        n_tests = n_tests + 1
        if (ok) then
            n_passed = n_passed + 1
            print '(A)', " [PASS] empty containers inline + 5-level nesting indentation"
        else
            print '(A)', " [FAIL] empty containers / nesting"
        end if
    end block

    print '(A, I2, A, I2, A)', "JSON Writer Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some assertions failed!"
        stop 1
    else
        print '(A)', "SUCCESS: JSON writer test fully PASSED (100%)."
    end if
    print '(A)', "============================================================"

end program test_json_io
