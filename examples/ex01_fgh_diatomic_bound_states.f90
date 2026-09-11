!> \brief 示例 1: 使用 Sinc-DVR / FGH 求解双原子分子 Morse 势能面束缚态
!> \details 计算 HF 双原子分子 Morse 振动态能级与本征波函数，并对比 Morse 解析解。
program ex01_fgh_diatomic_bound_states
    use general_module
    implicit none

    type(dvr_1d_t) :: dvr
    integer :: n_pts, stat, v, i, file_unit
    real(dp) :: r_min, r_max, mu_mass
    real(dp) :: d_e, r_e, beta, omega_e, omega_xe
    real(dp) :: e_exact
    real(dp), allocatable :: v_morse(:), eig_vals(:), eig_vecs(:, :)

    print '(A)', "=========================================================="
    print '(A)', "  Example 01: Diatomic Bound States via Sinc-DVR (FGH)   "
    print '(A)', "=========================================================="

    ! 1. 设置双原子分子物理参数 (HF 分子模型参数，原子单位 a.u.)
    d_e = 0.225_dp           ! 势阱深度 D_e (Hartree) ~ 6.12 eV
    r_e = 1.733_dp           ! 平衡核间距 R_e (Bohr) ~ 0.917 Angstrom
    beta = 1.174_dp          ! Morse 势刚度参数 beta (Bohr^-1)
    mu_mass = 1744.5_dp      ! 约化质量 mu (a.u.) ~ 0.957 amu

    ! 解析能级公式常数: E_v = omega_e*(v+1/2) - omega_xe*(v+1/2)^2
    omega_e = beta * sqrt(2.0_dp * d_e / mu_mass)
    omega_xe = (omega_e**2) / (4.0_dp * d_e)

    print '(A, F10.6, A)', " Morse Well Depth De: ", d_e, " Hartree"
    print '(A, F10.2, A)', " Harmonic Frequency:  ", omega_e * AU2CM, " cm^-1"
    print '(A, F10.2, A)', " Anharmonicity we*xe: ", omega_xe * AU2CM, " cm^-1"
    print '(A)', "----------------------------------------------------------"

    ! 2. 初始化 Sinc-DVR 空间网格
    r_min = 0.8_dp
    r_max = 6.0_dp
    n_pts = 256
    call dvr_sinc_init(r_min, r_max, n_pts, mu_mass, dvr)

    ! 3. 构造离散 Morse 势能曲线 V(R) = D_e * [1 - exp(-beta*(R - R_e))]^2
    allocate(v_morse(n_pts))
    allocate(eig_vals(n_pts))
    allocate(eig_vecs(n_pts, n_pts))

    do i = 1, n_pts
        v_morse(i) = d_e * (1.0_dp - exp(-beta * (dvr%x(i) - r_e)))**2
    end do

    ! 4. 调用 FGH 求解器求本征值与本征函数
    call fgh_solve_bound_states(dvr, v_morse, eig_vals, eig_vecs, stat)
    if (stat /= 0) then
        print '(A, I4)', "Error: FGH diagonalization failed with status: ", stat
        stop 1
    end if

    ! 5. 打印前 6 个束缚态能级与解析值对比
    print '(A5, A16, A16, A16)', "v", "FGH Energy(au)", "Exact Morse(au)", "Diff(cm^-1)"
    print '(A)', "----------------------------------------------------------"
    do v = 0, 5
        e_exact = omega_e * (real(v, dp) + 0.5_dp) - omega_xe * (real(v, dp) + 0.5_dp)**2
        print '(I4, 2F16.8, F16.4)', v, eig_vals(v + 1), e_exact, (eig_vals(v + 1) - e_exact) * AU2CM
    end do
    print '(A)', "----------------------------------------------------------"

    ! 6. 保存基态与前 3 激发态波函数至数据文件
    open(newunit=file_unit, file="diatomic_bound_states.dat", status="replace", action="write")
    write(file_unit, '(A)') "# R(Bohr)  V_pot(au)  Psi_v0  Psi_v1  Psi_v2  Psi_v3"
    do i = 1, n_pts
        write(file_unit, '(6ES16.8)') dvr%x(i), v_morse(i), &
            eig_vecs(i, 1) / sqrt(dvr%dx), &
            eig_vecs(i, 2) / sqrt(dvr%dx), &
            eig_vecs(i, 3) / sqrt(dvr%dx), &
            eig_vecs(i, 4) / sqrt(dvr%dx)
    end do
    close(file_unit)
    print '(A)', "Saved wavefunctions to: diatomic_bound_states.dat"

end program ex01_fgh_diatomic_bound_states
