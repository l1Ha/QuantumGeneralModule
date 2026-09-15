!> \brief 量子力学特殊函数与角动量耦合算法模块
!> \details 包含归一化勒让德多项式、Wigner 3j / Clebsch-Gordan 符号、刚体转动偶极与取向矩阵元。
!> \author LiHao
module mod_special_functions
    use mod_constants, only: dp, PI
    implicit none
    private

    public :: legendre_poly
    public :: assoc_legendre_poly
    public :: wigner_3j
    public :: clebsch_gordan
    public :: wigner_3j_half
    public :: clebsch_gordan_half
    public :: wigner_6j_half
    public :: wigner_9j_half
    public :: rot_matrix_cos_theta
    public :: rot_matrix_cos2_theta

contains

    !> \brief 计算普通勒让德多项式 P_l(x)，x ∈ [-1, 1]
    pure function legendre_poly(l, x) result(p)
        integer, intent(in) :: l
        real(dp), intent(in) :: x
        real(dp) :: p
        integer :: i
        real(dp) :: p0, p1, p2

        if (l < 0) then
            p = 0.0_dp
            return
        else if (l == 0) then
            p = 1.0_dp
            return
        else if (l == 1) then
            p = x
            return
        end if

        p0 = 1.0_dp
        p1 = x
        do i = 2, l
            p2 = ((2.0_dp * real(i, dp) - 1.0_dp) * x * p1 - (real(i, dp) - 1.0_dp) * p0) / real(i, dp)
            p0 = p1
            p1 = p2
        end do
        p = p1
    end function legendre_poly

    !> \brief 计算缔合勒让德函数 P_l^m(x) (带 Condon-Shortley 相位 (-1)^m)
    pure function assoc_legendre_poly(l, m, x) result(plm)
        integer, intent(in) :: l, m
        real(dp), intent(in) :: x
        real(dp) :: plm
        integer :: i, ll, abs_m
        real(dp) :: pmm, pmmp1, pcurrent, somx2, fact

        abs_m = abs(m)
        if (abs_m > l) then
            plm = 0.0_dp
            return
        end if

        ! 计算 P_m^m(x)
        pmm = 1.0_dp
        if (abs_m > 0) then
            somx2 = sqrt(max(0.0_dp, (1.0_dp - x) * (1.0_dp + x)))
            fact = 1.0_dp
            do i = 1, abs_m
                pmm = -pmm * fact * somx2
                fact = fact + 2.0_dp
            end do
        end if

        if (l == abs_m) then
            plm = pmm
            return
        end if

        ! 计算 P_{m+1}^m(x)
        pmmp1 = x * (2.0_dp * real(abs_m, dp) + 1.0_dp) * pmm
        if (l == abs_m + 1) then
            plm = pmmp1
            return
        end if

        ! 递推至 P_l^m(x)
        pcurrent = 0.0_dp
        do ll = abs_m + 2, l
            pcurrent = (x * (2.0_dp * real(ll, dp) - 1.0_dp) * pmmp1 - &
                       (real(ll + abs_m - 1, dp)) * pmm) / real(ll - abs_m, dp)
            pmm = pmmp1
            pmmp1 = pcurrent
        end do
        plm = pcurrent
    end function assoc_legendre_poly

    !> \brief 对数阶乘辅助函数 ln(n!)，避免高量子数阶乘溢出
    pure function log_factorial(n) result(res)
        integer, intent(in) :: n
        real(dp) :: res
        integer :: i

        res = 0.0_dp
        if (n <= 1) return
        do i = 2, n
            res = res + log(real(i, dp))
        end do
    end function log_factorial

    !> \brief 计算 Wigner 3j 符号 (Racah 展开公式)
    pure function wigner_3j(j1, j2, j3, m1, m2, m3) result(w3j)
        integer, intent(in) :: j1, j2, j3, m1, m2, m3
        real(dp) :: w3j
        integer :: t, t_min, t_max
        real(dp) :: delta_coeff, sum_term, sign_t
        real(dp) :: den

        w3j = 0.0_dp

        ! 选择定则校验
        if (m1 + m2 + m3 /= 0) return
        if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3) return
        if (j3 < abs(j1 - j2) .or. j3 > j1 + j2) return

        ! 三角系数 Delta(j1, j2, j3)
        delta_coeff = exp(0.5_dp * ( &
            log_factorial(j1 + j2 - j3) + log_factorial(j1 - j2 + j3) + log_factorial(-j1 + j2 + j3) - &
            log_factorial(j1 + j2 + j3 + 1) + &
            log_factorial(j1 + m1) + log_factorial(j1 - m1) + &
            log_factorial(j2 + m2) + log_factorial(j2 - m2) + &
            log_factorial(j3 + m3) + log_factorial(j3 - m3) ))

        ! 求和范围
        t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
        t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

        sum_term = 0.0_dp
        do t = t_min, t_max
            if (mod(t, 2) == 0) then
                sign_t = 1.0_dp
            else
                sign_t = -1.0_dp
            end if
            den = log_factorial(t) + &
                  log_factorial(j1 + j2 - j3 - t) + &
                  log_factorial(j1 - m1 - t) + &
                  log_factorial(j2 + m2 - t) + &
                  log_factorial(j3 - j2 + m1 + t) + &
                  log_factorial(j3 - j1 - m2 + t)
            sum_term = sum_term + sign_t * exp(-den)
        end do

        if (mod(j1 - j2 - m3, 2) /= 0) then
            w3j = -delta_coeff * sum_term
        else
            w3j = delta_coeff * sum_term
        end if
    end function wigner_3j

    !> \brief 计算 Clebsch-Gordan 系数 <j1 m1 j2 m2 | j3 m3>
    pure function clebsch_gordan(j1, m1, j2, m2, j3, m3) result(cg)
        integer, intent(in) :: j1, m1, j2, m2, j3, m3
        real(dp) :: cg
        real(dp) :: w3j

        if (m1 + m2 /= m3) then
            cg = 0.0_dp
            return
        end if

        w3j = wigner_3j(j1, j2, j3, m1, m2, -m3)
        if (mod(j1 - j2 + m3, 2) /= 0) then
            cg = -sqrt(2.0_dp * real(j3, dp) + 1.0_dp) * w3j
        else
            cg = sqrt(2.0_dp * real(j3, dp) + 1.0_dp) * w3j
        end if
    end function clebsch_gordan

    !> \brief 转动基底偶极矩阵元 <j, m | cos(theta) | j_prime, m>
    pure function rot_matrix_cos_theta(j, j_prime, m) result(mat_elem)
        integer, intent(in) :: j, j_prime, m
        real(dp) :: mat_elem
        real(dp) :: j_max_r, m_r
        integer :: j_max

        mat_elem = 0.0_dp
        if (abs(j - j_prime) /= 1) return
        if (abs(m) > min(j, j_prime)) return

        j_max = max(j, j_prime)
        j_max_r = real(j_max, dp)
        m_r = real(m, dp)

        mat_elem = sqrt((j_max_r**2 - m_r**2) / ((2.0_dp * j_max_r - 1.0_dp) * (2.0_dp * j_max_r + 1.0_dp)))
    end function rot_matrix_cos_theta

    !> \brief 转动基底取向矩阵元 <j, m | cos^2(theta) | j_prime, m>
    pure function rot_matrix_cos2_theta(j, j_prime, m) result(mat_elem)
        integer, intent(in) :: j, j_prime, m
        real(dp) :: mat_elem
        real(dp) :: j_r, m_r
        integer :: j_min

        mat_elem = 0.0_dp
        if (abs(m) > min(j, j_prime)) return
        m_r = real(m, dp)

        if (j == j_prime) then
            j_r = real(j, dp)
            mat_elem = ((j_r + 1.0_dp)**2 - m_r**2) / ((2.0_dp * j_r + 1.0_dp) * (2.0_dp * j_r + 3.0_dp))
            if (j > 0) then
                mat_elem = mat_elem + (j_r**2 - m_r**2) / ((2.0_dp * j_r - 1.0_dp) * (2.0_dp * j_r + 1.0_dp))
            end if
        else if (abs(j - j_prime) == 2) then
            j_min = min(j, j_prime)
            j_r = real(j_min, dp)
            mat_elem = sqrt(((j_r + 1.0_dp)**2 - m_r**2) * ((j_r + 2.0_dp)**2 - m_r**2) / &
                       ((2.0_dp * j_r + 1.0_dp) * (2.0_dp * j_r + 3.0_dp)**2 * (2.0_dp * j_r + 5.0_dp)))
        end if
    end function rot_matrix_cos2_theta

    !> \brief 通用半整数 Wigner 3j 符号 (输入参数均为 2 倍角动量量子数，避免浮点截断)
    pure function wigner_3j_half(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3) result(w3j)
        integer, intent(in) :: two_j1, two_j2, two_j3, two_m1, two_m2, two_m3
        real(dp) :: w3j

        integer :: a, b, c, d, ja_p, ja_m, jb_p, jb_m, jc_p, jc_m
        integer :: t, t_min, t_max, phase
        real(dp) :: log_delta, log_pref, s, den, sign_t, sign_tot

        w3j = 0.0_dp
        if (two_m1 + two_m2 + two_m3 /= 0) return
        if (abs(two_m1) > two_j1 .or. abs(two_m2) > two_j2 .or. abs(two_m3) > two_j3) return
        if (two_j3 < abs(two_j1 - two_j2) .or. two_j3 > two_j1 + two_j2) return
        if (mod(two_j1 + two_j2 + two_j3, 2) /= 0) return
        if (mod(two_j1 - two_m1, 2) /= 0 .or. mod(two_j2 - two_m2, 2) /= 0 .or. mod(two_j3 - two_m3, 2) /= 0) return

        a = (two_j1 + two_j2 - two_j3) / 2
        b = (two_j1 - two_j2 + two_j3) / 2
        c = (-two_j1 + two_j2 + two_j3) / 2
        d = (two_j1 + two_j2 + two_j3) / 2 + 1

        log_delta = 0.5_dp * (log_factorial(a) + log_factorial(b) + log_factorial(c) - log_factorial(d))

        ja_p = (two_j1 + two_m1) / 2
        ja_m = (two_j1 - two_m1) / 2
        jb_p = (two_j2 + two_m2) / 2
        jb_m = (two_j2 - two_m2) / 2
        jc_p = (two_j3 + two_m3) / 2
        jc_m = (two_j3 - two_m3) / 2

        log_pref = log_delta + 0.5_dp * (log_factorial(ja_p) + log_factorial(ja_m) + &
                                         log_factorial(jb_p) + log_factorial(jb_m) + &
                                         log_factorial(jc_p) + log_factorial(jc_m))

        t_min = max(0, (two_j2 - two_j3 - two_m1) / 2, (two_j1 - two_j3 + two_m2) / 2)
        t_max = min(a, ja_m, jb_p)

        s = 0.0_dp
        do t = t_min, t_max
            den = log_factorial(t) + log_factorial(a - t) + &
                  log_factorial(ja_m - t) + log_factorial(jb_p - t) + &
                  log_factorial((two_j3 - two_j2 + two_m1) / 2 + t) + &
                  log_factorial((two_j3 - two_j1 - two_m2) / 2 + t)
            if (mod(t, 2) /= 0) then
                sign_t = -1.0_dp
            else
                sign_t = 1.0_dp
            end if
            s = s + sign_t * exp(-den)
        end do

        phase = (two_j1 - two_j2 - two_m3) / 2
        if (mod(phase, 2) /= 0) then
            sign_tot = -1.0_dp
        else
            sign_tot = 1.0_dp
        end if

        w3j = sign_tot * exp(log_pref) * s
    end function wigner_3j_half

    !> \brief 通用半整数 Clebsch-Gordan 系数 <j1 m1 j2 m2 | j3 m3> (输入均为 2 倍角动量量子数)
    pure function clebsch_gordan_half(two_j1, two_m1, two_j2, two_m2, two_j3, two_m3) result(cg)
        integer, intent(in) :: two_j1, two_m1, two_j2, two_m2, two_j3, two_m3
        real(dp) :: cg
        real(dp) :: w3j, sign_ph
        integer  :: phase

        if (two_m1 + two_m2 /= two_m3) then
            cg = 0.0_dp
            return
        end if

        w3j = wigner_3j_half(two_j1, two_j2, two_j3, two_m1, two_m2, -two_m3)
        phase = (two_j1 - two_j2 + two_m3) / 2
        if (mod(phase, 2) /= 0) then
            sign_ph = -1.0_dp
        else
            sign_ph = 1.0_dp
        end if

        cg = sign_ph * sqrt(real(two_j3 + 1, dp)) * w3j
    end function clebsch_gordan_half

    !> \brief 辅助三角系数 Delta(j1, j2, j3) (输入为 2 倍角动量)
    pure function triangle_half(two_a, two_b, two_c) result(tri)
        integer, intent(in) :: two_a, two_b, two_c
        real(dp) :: tri
        integer  :: a, b, c, d
        real(dp) :: log_d

        tri = 0.0_dp
        if (mod(two_a + two_b + two_c, 2) /= 0) return
        if (two_c < abs(two_a - two_b) .or. two_c > two_a + two_b) return

        a = (two_a + two_b - two_c) / 2
        b = (two_a - two_b + two_c) / 2
        c = (-two_a + two_b + two_c) / 2
        d = (two_a + two_b + two_c) / 2 + 1

        log_d = 0.5_dp * (log_factorial(a) + log_factorial(b) + log_factorial(c) - log_factorial(d))
        tri = exp(log_d)
    end function triangle_half

    !> \brief 通用半整数 Wigner 6j 符号 (输入均为 2 倍角动量量子数)
    pure function wigner_6j_half(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6) result(w6j)
        integer, intent(in) :: two_j1, two_j2, two_j3, two_j4, two_j5, two_j6
        real(dp) :: w6j
        real(dp) :: d1, d2, d3, d4, pref, s, num, den, sign_t
        integer  :: t_min, t_max, t

        w6j = 0.0_dp
        d1 = triangle_half(two_j1, two_j2, two_j3)
        d2 = triangle_half(two_j1, two_j5, two_j6)
        d3 = triangle_half(two_j4, two_j2, two_j6)
        d4 = triangle_half(two_j4, two_j5, two_j3)
        if (d1 < 1.0e-30_dp .or. d2 < 1.0e-30_dp .or. d3 < 1.0e-30_dp .or. d4 < 1.0e-30_dp) return

        pref = d1 * d2 * d3 * d4
        t_min = max((two_j1 + two_j2 + two_j3) / 2, &
                    (two_j1 + two_j5 + two_j6) / 2, &
                    (two_j4 + two_j2 + two_j6) / 2, &
                    (two_j4 + two_j5 + two_j3) / 2)
        t_max = min((two_j1 + two_j2 + two_j4 + two_j5) / 2, &
                    (two_j2 + two_j3 + two_j5 + two_j6) / 2, &
                    (two_j3 + two_j1 + two_j6 + two_j4) / 2)

        s = 0.0_dp
        do t = t_min, t_max
            num = log_factorial(t + 1)
            den = log_factorial(t - (two_j1 + two_j2 + two_j3) / 2) + &
                  log_factorial(t - (two_j1 + two_j5 + two_j6) / 2) + &
                  log_factorial(t - (two_j4 + two_j2 + two_j6) / 2) + &
                  log_factorial(t - (two_j4 + two_j5 + two_j3) / 2) + &
                  log_factorial((two_j1 + two_j2 + two_j4 + two_j5) / 2 - t) + &
                  log_factorial((two_j2 + two_j3 + two_j5 + two_j6) / 2 - t) + &
                  log_factorial((two_j3 + two_j1 + two_j6 + two_j4) / 2 - t)
            if (mod(t, 2) /= 0) then
                sign_t = -1.0_dp
            else
                sign_t = 1.0_dp
            end if
            s = s + sign_t * exp(num - den)
        end do

        w6j = pref * s
    end function wigner_6j_half

    !> \brief 通用半整数 Wigner 9j 符号 (输入均为 2 倍角动量量子数)
    pure function wigner_9j_half(two_j11, two_j12, two_j13, &
                                 two_j21, two_j22, two_j23, &
                                 two_j31, two_j32, two_j33) result(w9j)
        integer, intent(in) :: two_j11, two_j12, two_j13
        integer, intent(in) :: two_j21, two_j22, two_j23
        integer, intent(in) :: two_j31, two_j32, two_j33
        real(dp) :: w9j
        integer  :: two_k, two_k_min, two_k_max
        real(dp) :: w1, w2, w3, sign_k

        w9j = 0.0_dp
        two_k_min = max(abs(two_j11 - two_j33), abs(two_j32 - two_j21), abs(two_j12 - two_j23))
        two_k_max = min(two_j11 + two_j33, two_j32 + two_j21, two_j12 + two_j23)

        do two_k = two_k_min, two_k_max, 2
            w1 = wigner_6j_half(two_j11, two_j21, two_j31, two_j32, two_j33, two_k)
            w2 = wigner_6j_half(two_j12, two_j22, two_j32, two_j21, two_k, two_j23)
            w3 = wigner_6j_half(two_j13, two_j23, two_j33, two_k, two_j11, two_j12)
            if (mod(two_k, 2) /= 0) then
                sign_k = -1.0_dp
            else
                sign_k = 1.0_dp
            end if
            w9j = w9j + real(two_k + 1, dp) * sign_k * w1 * w2 * w3
        end do
    end function wigner_9j_half

end module mod_special_functions
