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
        real(dp) :: delta_coeff, sum_term, term, sign_t
        real(dp) :: num, den

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
        real(dp) :: j_r, m_r

        mat_elem = 0.0_dp
        if (abs(m) > min(j, j_prime)) return

        j_r = real(min(j, j_prime), dp)
        m_r = real(m, dp)

        if (j_prime == j + 1) then
            mat_elem = sqrt(((j_r + 1.0_dp)**2 - m_r**2) / ((2.0_dp * j_r + 1.0_dp) * (2.0_dp * j_r + 3.0_dp)))
        else if (j_prime == j - 1) then
            mat_elem = sqrt((j_r**2 - m_r**2) / ((2.0_dp * j_r - 1.0_dp) * (2.0_dp * j_r + 1.0_dp)))
        end if
    end function rot_matrix_cos_theta

    !> \brief 转动基底取向矩阵元 <j, m | cos^2(theta) | j_prime, m>
    pure function rot_matrix_cos2_theta(j, j_prime, m) result(mat_elem)
        integer, intent(in) :: j, j_prime, m
        real(dp) :: mat_elem
        real(dp) :: j_r, m_r

        mat_elem = 0.0_dp
        if (abs(m) > min(j, j_prime)) return
        j_r = real(j, dp)
        m_r = real(m, dp)

        if (j_prime == j) then
            mat_elem = ((j_r + 1.0_dp)**2 - m_r**2) / ((2.0_dp * j_r + 1.0_dp) * (2.0_dp * j_r + 3.0_dp)) + &
                       (j_r**2 - m_r**2) / ((2.0_dp * j_r - 1.0_dp) * (2.0_dp * j_r + 1.0_dp))
        else if (j_prime == j + 2) then
            mat_elem = sqrt(((j_r + 1.0_dp)**2 - m_r**2) * ((j_r + 2.0_dp)**2 - m_r**2) / &
                       ((2.0_dp * j_r + 1.0_dp) * (2.0_dp * j_r + 3.0_dp)**2 * (2.0_dp * j_r + 5.0_dp)))
        else if (j_prime == j - 2) then
            mat_elem = sqrt((j_r**2 - m_r**2) * ((j_r - 1.0_dp)**2 - m_r**2) / &
                       ((2.0_dp * j_r + 1.0_dp) * (2.0_dp * j_r - 1.0_dp)**2 * (2.0_dp * j_r - 3.0_dp)))
        end if
    end function rot_matrix_cos2_theta

end module mod_special_functions
