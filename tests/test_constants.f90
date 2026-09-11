!> \brief 物理常数与单位转换单元测试
program test_constants
    use general_module
    implicit none

    integer :: n_tests = 0
    integer :: n_passed = 0
    real(dp) :: tol = 1.0e-10_dp
    real(dp) :: val_in, val_au, val_out, err

    print '(A)', "=================================================="
    print '(A)', "       GeneralModule Unit Tests: Constants        "
    print '(A)', "=================================================="

    ! 1. 基础数学常数校验
    call assert_close("PI definition", PI, 4.0_dp * atan(1.0_dp), tol)
    call assert_close("TWOPI definition", TWOPI, 2.0_dp * PI, tol)
    call assert_close("SQRTPI definition", SQRTPI, sqrt(PI), tol)

    ! 2. 时间单位互转 (fs <-> au)
    val_in = 100.0_dp ! 100 fs
    val_au = to_au(val_in, 'fs')
    val_out = from_au(val_au, 'fs')
    call assert_close("fs -> au -> fs roundtrip", val_in, val_out, tol)

    ! 3. 能量单位互转 (eV <-> au)
    val_in = 13.605693_dp ! 13.6 eV
    val_au = to_au(val_in, 'eV')
    val_out = from_au(val_au, 'eV')
    call assert_close("eV -> au -> eV roundtrip", val_in, val_out, tol)

    ! 4. 波数互转 (cm-1 <-> au)
    val_in = 219474.63_dp
    val_au = to_au(val_in, 'cm-1')
    val_out = from_au(val_au, 'cm-1')
    call assert_close("cm-1 -> au -> cm-1 roundtrip", val_in, val_out, tol)

    ! 5. 长度单位互转 (Angstrom <-> au)
    val_in = 1.54_dp ! 1.54 Angstrom
    val_au = to_au(val_in, 'ang')
    val_out = from_au(val_au, 'ang')
    call assert_close("Angstrom -> au -> Angstrom roundtrip", val_in, val_out, tol)

    ! 6. 偶极矩互转 (Debye <-> au)
    val_in = 1.85_dp ! 水分子偶极矩 ~1.85 Debye
    val_au = to_au(val_in, 'debye')
    val_out = from_au(val_au, 'debye')
    call assert_close("Debye -> au -> Debye roundtrip", val_in, val_out, tol)

    ! 7. 激光光强互转 (W/cm2 <-> au)
    val_in = 1.0e14_dp ! 10^14 W/cm^2
    val_au = to_au(val_in, 'w/cm2')
    val_out = from_au(val_au, 'w/cm2')
    call assert_close("W/cm2 -> au -> W/cm2 roundtrip", val_in, val_out, 1.0e-5_dp)

    ! 8. 质量互转 (amu <-> au)
    val_in = 1.007825_dp ! 氢原子质量 amu
    val_au = to_au(val_in, 'amu')
    val_out = from_au(val_au, 'amu')
    call assert_close("amu -> au -> amu roundtrip", val_in, val_out, tol)

    print '(A)', "--------------------------------------------------"
    print '(A, I2, A, I2, A)', "Constants Tests: ", n_passed, " / ", n_tests, " PASSED."
    if (n_passed /= n_tests) then
        print '(A)', "FAILED: Some tests did not pass!"
        stop 1
    else
        print '(A)', "SUCCESS: All constants tests passed."
    end if

contains

    subroutine assert_close(name, expected, actual, eps)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: expected, actual, eps
        real(dp) :: diff

        n_tests = n_tests + 1
        diff = abs(expected - actual)
        if (diff <= eps) then
            n_passed = n_passed + 1
            print '(A, A35, A)', " [PASS] ", name, ""
        else
            print '(A, A35, A, E12.5, A, E12.5, A, E12.5)', " [FAIL] ", name, &
                " (Exp: ", expected, ", Act: ", actual, ", Diff: ", diff, ")"
        end if
    end subroutine assert_close

end program test_constants
