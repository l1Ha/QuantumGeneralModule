# GeneralModule: 现代量子动力学通用算法库 (Fortran 2008 / Python)

[![CI](https://github.com/l1Ha/QuantumGeneralModule/actions/workflows/ci.yml/badge.svg)](https://github.com/l1Ha/QuantumGeneralModule/actions/workflows/ci.yml)
[![Fortran 2008](https://img.shields.io/badge/Fortran-2008-734f96.svg)](https://fortran-lang.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Tests: 115/115 Pass](https://img.shields.io/badge/Tests-115%2F115%20Pass%20(100%25)-brightgreen.svg)](tests/)

`GeneralModule` 是一个面向超快强场物理、分子光物理与量子动力学模拟的现代化通用科学计算算法库。该库遵循严格的 **Fortran 2008 规范**，具备高数值精度、零外部动态库强依赖、模块化架构与出色的 AI Agent 友好性。

---

## 🌟 核心特性

1. **零外部库依赖 (Zero External Dependencies)**
   - 内部集成高精度 Householder QL 实对称矩阵本征求解器、Gauss-Jordan 全主元实/复方阵求逆与 Cooley-Tukey 1D/2D 快速傅里叶变换（FFT）。
   - 纯 Fortran 自包含样条插值、Lindblad 主方程、Krotov 最优控制、通用多通道定态密耦（Johnson Log-Derivative）与含时连续态散射矩阵求解器，无需强制链接外部 LAPACK/BLAS 或 FFTW，开箱即用。
2. **现代 Fortran 2008 标准设计**
   - 统一强类型参数定义（`real(dp) => real64`）。
   - 纯函数（`pure function`）与显式 `intent(in/out/inout)` 契约，杜绝隐式全局变量副作用。
3. **AI 友好型结构化接口 (AI-Friendly)**
   - 算法模块支持统一顶层聚合入口：`use general_module`。
   - 参数配置采用清晰的派生类型（Derived Types，如 `pulse_config_t`, `absorbing_boundary_t`, `dvr_1d_t`, `spline_1d_t`, `scattering_state_t`, `multichannel_result_t`, `td_scattering_channel_t`），自解释、低耦合、便于大语言模型精确构造与调用。
4. **全链路双语生态支持**
   - 附带标准 Python 伴侣分析包 `pygenmod`，无缝衔接参数预计算、波包与散射长度可视化、发表级绘图（含一键动力学出图流水线 `plot_rovibrational_dynamics.py`）。
5. **全自动 CI/CD 持续集成**
   - 内置 GitHub Actions 跨平台持续集成（Ubuntu / macOS），全自动化执行 12 大测试套件与 Python 验证。

---

## 📂 目录结构与架构

```text
GeneralModule/
├── .github/workflows/ci.yml       # 🚀 GitHub Actions 跨平台 CI 持续集成工作流
├── CONFIG_GUIDE.md                # 📖 详尽配置与部署指南（全编译平台、参数全典、AI提示词模板）
├── fpm.toml                       # Fortran Package Manager 配置文件
├── CMakeLists.txt                 # CMake 跨平台构建系统
├── README.md                      # 本文档
├── .gitignore                     # Git 忽略规则
├── src/                           # 核心 Fortran 源代码 (20 核心模块 + 1 聚合入口)
│   ├── mod_constants.f90          # 1. 物理常数与各单位 a.u. 双向转换
│   ├── mod_special_functions.f90  # 2. 勒让德、Wigner 3j、CG、转动偶极/取向矩阵元
│   ├── mod_linear_algebra.f90     # 3. 对称矩阵本征求解 (EISPACK TRED2/TQL2)、实/复方阵求逆、1D/2D FFT
│   ├── mod_dvr_grid.f90           # 4. Sinc-DVR、Legendre-DVR、FGH 束缚态求解、格点期望值
│   ├── mod_laser_pulse.f90        # 5. 超快强场脉冲时域合成、矢量势、椭偏场与便捷构造器
│   ├── mod_absorbing_boundary.f90 # 6. 复吸收势边界（CAP）、概率流与存活范数
│   ├── mod_thermal_ensemble.f90   # 7. 玻尔兹曼转振分布、热系综统计平均
│   ├── mod_wavepacket_propagator.f90 # 8. 分裂算符 (1D/2D Split-Operator)、RK4、Bloch、ABM4
│   ├── mod_coulomb_atomic.f90     # 9. 强场原子模型势、Keldysh参数、有质动力能与ADK电离率
│   ├── mod_hhg_spectra.f90        # 10. Ehrenfest偶极加速度、HHG谐波谱、Gabor时频与SFA模型
│   ├── mod_chebyshev_propagator.f90 # 11. 切比雪夫多项式大步长推进器与能谱窗算子
│   ├── mod_multistate_coupling.f90 # 12. 多势能面非绝热耦合动力学与Landau-Zener跃迁
│   ├── mod_rovibrational.f90      # 13. 分子转振耦合、Franck-Condon因子、转动常数、偶极矩阵与STIRAP脉冲
│   ├── mod_io_utils.f90           # 14. 科学数据多列保存、矩阵导出、控制台横幅与进度条监测
│   ├── mod_interpolation.f90      # 15. 高精度自然/固定导数三次样条插值、解析一阶/二阶导数与势能面渐近外推
│   ├── mod_photofragment_flux.f90 # 16. 自相关函数与吸收截面谱、渐近散射振幅与光解离碎片动能释放谱(KER)
│   ├── mod_open_quantum.f90       # 17. 开放量子系统 Lindblad 耗散主方程、自发跃迁/退相位弛豫、纯度与冯·诺依曼熵
│   ├── mod_optimal_control.f90    # 18. 量子最优控制理论 Krotov 算法、目标保真度与激光场原位迭代优化
│   ├── mod_ti_scattering.f90      # 19. 非含时散射理论、零能 Numerov/Log-Derivative 散射长度、相移、K/S/T 矩阵与通用多通道密耦
│   ├── mod_td_scattering.f90      # 20. 含时波包散射理论、连续态能量通量透射率 T(E)、含时 S 矩阵提取与 Möller 动量投影
│   └── general_module.f90         # 顶层聚合入口模块 (use general_module)
├── tests/                         # 自动化单元测试套件 (12 个套件，100% 全部通过)
│   ├── test_constants.f90
│   ├── test_special_functions.f90
│   ├── test_dvr_grid.f90
│   ├── test_laser_pulse.f90
│   ├── test_propagators.f90
│   ├── test_atomic_hhg.f90
│   ├── test_laser_rovibrational_control.f90 # 激光调控分子转振态布居转移综合测试
│   ├── test_interpolation.f90     # 三次样条插值与渐近外推测试
│   ├── test_photofragment_flux.f90# 自相关吸收谱与碎片 KER 分支比测试
│   ├── test_open_quantum_opt.f90  # Lindblad 耗散退相干与 Krotov 最优控制测试
│   ├── test_ti_scattering.f90     # 非含时散射长度、相移、S矩阵与多通道密耦测试
│   ├── test_td_scattering.f90     # 含时波包散射透射谱、S矩阵元与 Möller 投影测试
│   └── run_all_tests.sh           # 自动化测试运行脚本 (100% Pass, 115/115 断言)
├── examples/                      # 典型物理应用算例 (6 大完整前沿算例)
│   ├── ex01_fgh_diatomic_bound_states.f90 # 双原子 Morse 势能级与波函数求解
│   ├── ex02_pulse_synthesis.f90           # 啁啾、双色、太赫兹脉冲时频生成
│   ├── ex03_split_operator_1d.f90         # 1D 波包动力学演化与 CAP 吸收边界
│   ├── ex04_field_free_orientation.f90    # 刚体转子无场定向与玻尔兹曼热平均
│   ├── ex05_hhg_lewenstein_spectrum.f90   # 强场阿秒高次谐波发射与半经典截止能
│   ├── ex06_two_state_nonadiabatic.f90    # 双态避差穿越非绝热动力学与分支比
│   └── build_examples.sh          # 算例编译运行脚本
└── python/                        # Python 辅助分析与可视化套件 (pygenmod)
    ├── pyproject.toml
    ├── test_pygenmod.py           # Python 单元测试 (100% Pass, 10/10 测试)
    ├── plot_rovibrational_dynamics.py # 出版级分子转振受控动力学一键绘图管道
    └── pygenmod/
        ├── __init__.py
        ├── constants.py
        ├── pulse.py
        ├── dvr.py
        ├── coulomb.py
        ├── hhg.py
        ├── multistate.py
        ├── rovibrational.py       # 转振态索引映射、FC因子、转动常数与跃迁偶极
        ├── scattering.py          # 散射长度(Numerov/方势阱/范德华)与波函数渐近线绘图
        └── visualizer.py          # 发表级科学绘图工具
```

---

## 🧩 核心模块 API 说明

### 1. 物理常数与单位转换 (`mod_constants`)
提供 CODATA 推荐的最新高精度基础物理常数，以及原子单位（a.u.）与常用实验单位的双向纯函数转换：
- **常数**：`PI`, `TWOPI`, `HALFPI`, `SQRTPI`, `EYE`, `C_LIGHT`, `HBAR`, `M_E`, `CHARGE_E`, `KB`, `AMU2AU` 等。
- **纯函数接口**：
  ```fortran
  ! 将实际物理单位数值转换为原子单位 a.u.
  val_au = to_au(val, "fs")       ! 支持: fs, ps, s, eV, cm-1, J, K, Angstrom, nm, m, amu, MV/cm, Debye, W/cm2
  ! 将原子单位 a.u. 数值转换为实际物理单位
  val_si = from_au(val_au, "eV")
  ```

### 2. 量子力学特殊函数 (`mod_special_functions`)
- `legendre_poly(l, x)`: 计算普通勒让德多项式 $P_l(x)$。
- `assoc_legendre_poly(l, m, x)`: 计算缔合勒让德多项式 $P_l^m(x)$（含 Condon-Shortley 相位）。
- `wigner_3j(j1, j2, j3, m1, m2, m3)`: 计算 Wigner 3j 符号，自动校验三角定则与磁量子数守恒。
- `clebsch_gordan(j1, m1, j2, m2, j3, m3)`: 计算 Clebsch-Gordan 耦合系数 $\langle j_1 m_1 j_2 m_2 | j_3 m_3 \rangle$。
- `rot_matrix_cos_theta(j, j_prime, m)`: 刚体转子偶极跃迁矩阵元 $\langle j, m | \cos\theta | j', m \rangle$。
- `rot_matrix_cos2_theta(j, j_prime, m)`: 刚体转子极化取向矩阵元 $\langle j, m | \cos^2\theta | j', m \rangle$。

### 3. 线性代数与 FFT (`mod_linear_algebra`)
- `diag_symmetric_matrix(n, a_in, d, z, stat)`: 基于 Householder 三对角化与隐式位移 QL 迭代求解 $n \times n$ 实对称矩阵的全部本征值与本征向量，自动按升序排列。
- `fft_1d(data_vec, isign)`: 原位一维复数快速傅里叶变换（Cooley-Tukey 算法，`isign = -1` 为时域到频域，`+1` 为逆变换）。
- `fft_2d(data_mat, isign)`: 二维复数快速傅里叶变换。

### 4. 离散变量表象网格 (`mod_dvr_grid`)
- `dvr_sinc_init(x_min, x_max, n_pts, mass, dvr)`: 初始化 Colbert-Miller Sinc-DVR 网格与解析动能算符矩阵 $T_{ij}$。
- `dvr_legendre_init(n_pts, dvr)`: 使用 Newton-Raphson 法初始化 Gauss-Legendre DVR 零点、高斯积分权重与转动动能矩阵。
- `fgh_solve_bound_states(dvr, v_pot, eig_vals, eig_vecs, stat)`: Fourier Grid Hamiltonian (FGH) 求解任意分子一维势能曲线的束缚态能级与波函数。

### 5. 激光脉冲合成 (`mod_laser_pulse`)
- 派生类型：`type(pulse_config_t)`
  - 支持形态枚举：`PULSE_GAUSSIAN`, `PULSE_SIN2`, `PULSE_FLATTOP`, `PULSE_CHIRP`, `PULSE_TWOCOLOR`, `PULSE_THZ_TRAIN`。
- `pulse_envelope(t, cfg)`: 瞬时无量纲包络函数 $f(t) \in [0, 1]$。
- `pulse_electric_field(t, cfg)`: 瞬时激光电场强度 $E(t)$。
- `pulse_stark_shift(t, cfg, alpha_parallel, alpha_perp, theta)`: 瞬时动力学 AC Stark 能级位移。
- `pulse_generate_timeseries(t_arr, cfg, e_arr, env_arr)`: 矢量化批量生成时域电场与包络。

### 6. 吸收边界与通量监测 (`mod_absorbing_boundary`)
- 派生类型：`type(absorbing_boundary_t)`（支持 `CAP_SIN2` 与 `CAP_POLYNOMIAL`）。
- `cap_init(r_start, r_end, strength, cap_type, cap_obj)`: 初始化复吸收势边界。
- `cap_evaluate(r, cap_obj)`: 计算特定网格位置的虚部吸收势 $-i W(r)$。
- `cap_apply_mask(r_grid, dt, cap_obj, psi)`: 施加平滑吸收掩膜衰减外部波包 $\psi(r) \leftarrow \psi(r) \cdot \exp(-W(r) \Delta t)$。
- `calculate_probability_flux(r_grid, mass, psi, idx_detect)`: 采用二阶中心差分计算渐近监测面处的瞬时量子概率流密度。
- `calculate_norm_inside(r_grid, psi, r_boundary)`: 计算反应区（$R \le R_0$）内部波包存活总几率。

### 7. 热统计力学与系综平均 (`mod_thermal_ensemble`)
- `boltzmann_rotational_weights(temp_k, b_rot, j_max, weights, z_rot)`: 计算有限温度 $T$ 下刚体转动各初态 $J$ 的玻尔兹曼权重与配分函数 $Z_{rot}$。
- `boltzmann_vibrational_weights(temp_k, omega_e, v_max, weights)`: 计算简谐振动态的玻尔兹曼权重。
- `thermal_average_1d(observables, weights)` / `thermal_average_2d(obs_matrix, weights, avg_series)`: 针对静态或动力学含时矩阵进行统计加权热平均。
- `bose_einstein_factor(omega_au, temp_k)`: 开放系统热浴关联的 Bose-Einstein 分布因子。

### 8. 波包推进器与积分器 (`mod_wavepacket_propagator`)
- `propagate_split_operator_1d(psi, v_pot, dx, mass, dt)`: 基于 FFT 的一维二阶辛对称分裂算符（Split-Operator）单步波包推进器。
- `propagate_split_operator_2d(psi, v_pot, dx, dy, mass_x, mass_y, dt)`: 二维自由度 Split-Operator FFT 波包推进器。
- `rk4_step(y, t, dt, f_deriv)`: 通用 4 阶 Runge-Kutta 状态向量积分器。
- `solve_bloch_two_level(r_bloch, rabi_freq, detuning, dt)`: 光学 Bloch 方程数值求解器，保模长推进 Bloch 矢量 $[u, v, w]$。
- `abm4_step(y, fn, fn_m1, fn_m2, fn_m3, dt, t_next, f_eval)`: 4 阶 Adams-Bashforth-Moulton 预估-校正多步法推进器。

### 9. 强场原子模型与电离率 (`mod_coulomb_atomic`)
- `get_atom_config(name, cfg)`: 调取 H, He, Ne, Ar, Kr, Xe 等目标原子的单活性电子（SAE）电离能与软核参数。
- `soft_core_coulomb_potential(x, soft_a, z_eff)`: 一维软核库仑势 $V(x) = -Z / \sqrt{x^2 + a^2}$。
- `soft_core_coulomb_derivative(x, soft_a, z_eff)`: 软核库仑势空间一阶导数（用于计算偶极受力）。
- `keldysh_parameter(omega, e_peak, ip_au)`: Keldysh 强场电离区段判据 $\gamma$。
- `ponderomotive_energy(e_peak, omega)`: 电子在激光场中的有质动力能 $U_p = E_0^2 / (4\omega^2)$。
- `hhg_cutoff_energy(ip_au, e_peak, omega)`: 半经典高次谐波截断能量定律 $E_{cutoff} = I_p + 3.17 U_p$。
- `adk_ionization_rate(e_field, ip_au, z_eff, l, m)`: 静态/准静态场 ADK 隧穿电离率。

### 10. 高次谐波与偶极时频分析 (`mod_hhg_spectra`)
- `calculate_dipole_length(x_grid, dx, psi)`: 长度表象瞬时偶极矩 $d(t) = \langle\psi| x |\psi\rangle$。
- `calculate_dipole_acceleration(dv_dx, dx, psi, e_field)`: 基于 Ehrenfest 定理的瞬时偶极加速度 $a(t) = -\langle\psi| \partial V/\partial x + E(t) |\psi\rangle$。
- `hhg_power_spectrum(a_t, dt, omega_arr, spectrum)`: 高次谐波发射功率谱 $S(\omega) = |\text{FFT}[a(t) W(t)]|^2$（内建平滑加窗抑制泄露）。
- `gabor_transform_point(t_arr, a_t, dt, t_0, omega, sigma)`: Gabor 时频小波变换，解析各阶阿秒谐波发射时间轮廓。
- `lewenstein_sfa_dipole(t, t_arr, a_field, e_field, ip_au, dt)`: 基于 Lewenstein 强场近似（SFA）鞍点行动量积分的半经典含时偶极矩。

### 11. 切比雪夫推进器与能谱滤波 (`mod_chebyshev_propagator`)
- `chebyshev_propagate_step(psi, dt, e_min, e_max, h_mult, order)`: 切比雪夫多项式展开单步推进器，支持超大时间步长与机器精度级范数守恒。
- `window_operator_pes(psi, dx, e_k, gamma, eig_vals, eig_vecs)`: Schafer-Kulander 能量窗算子，精确提取能量分辨的光电子动能谱（PES）。

### 12. 多势能面非绝热动力学 (`mod_multistate_coupling`)
- `propagate_split_operator_2channel(psi1, psi2, v11, v22, v12, dx, mass, dt)`: 双通道非绝热耦合核波包二阶辛对称推进器（解析 $2 \times 2$ 幺正矩阵指数）。
- `landau_zener_probability(v12, velocity, delta_slope)`: 经典 Landau-Zener 避差穿越跃迁几率计算。
- `calculate_channel_populations(psi1, psi2, dx, pop1, pop2, ratio)`: 各电子通道波包总几率与非绝热转移分支比。

### 13. 分子转振耦合与激光调控 (`mod_rovibrational`)
- `calc_franck_condon_factors(chi_a, chi_b, dx, fc_mat)`: 计算双原子分子振动态间 Franck-Condon 重叠因子矩阵 $FC(v, v') = |\langle\chi_v | \chi_{v'}\rangle|^2$。
- `calc_vibrational_dipole_matrix(chi, dipole_grid, dx, dip_mat)`: 计算核间距依赖偶极矩在振动态基底下的跃迁积分 $M(v, v') = \langle\chi_v | \mu(R) | \chi_{v'}\rangle$。
- `calc_rotational_constants_bv(chi, r_grid, dx, mass, b_v)`: 严格数值积分各振动态有效转动常数 $B_v = \langle\chi_v | \frac{\hbar^2}{2\mu R^2} | \chi_v\rangle$。
- `build_rovibrational_hamiltonian(v_max, j_max, e_vib, b_v, h_diag)`: 构造 $|v, J\rangle$ 转振空间本征能级对角哈密顿量 $E(v, J) = E_v + B_v J(J+1)$。
- `build_rovibrational_dipole_matrix(v_max, j_max, dip_vib, dip_mat)`: 构造严格满足偶极选择定则 $\Delta J = \pm 1, \Delta M = 0$ 的全转振跃迁矩阵。
- `build_rovibrational_polarizability_matrix(v_max, j_max, alpha_vib, polar_mat)`: 构造满足极化选择定则 $\Delta J = 0, \pm 2$ 的激光取向耦合矩阵。
- `create_stirap_pulses(peak_p, peak_s, dur_p, dur_s, delay, w_p, w_s, cfg_p, cfg_s)`: 便捷生成受激拉曼绝热通道（STIRAP）Stokes 超前 Pump 脉冲对。
- `rovibrational_state_index(v, j, j_max)` / `rovibrational_state_unindex(idx, j_max, v, j)`: 二维量子数 $(v, J)$ 与一维基底线性索引快速双向互转。

### 14. 科学计算 I/O 与诊断工具 (`mod_io_utils`)
- `save_data_table_1d(filename, x, y, ...)`: 导出双列 ASCII 科学数据表，自动生成格式化注释头。
- `save_data_table_2d(filename, x, y_mat, col_names, header)`: 导出多自由度动力学时序演化数据表（多列二维矩阵）。
- `save_matrix_dat(filename, mat, header)`: 导出二维实对称或势能面方阵数据。
- `print_banner(title, width)`: 控制台美化打印计算任务标题横幅。
- `print_progress_bar(current, total, prefix)`: 动力学时域演化单行就地刷新进度条。

### 15. 高精度样条插值与势能面外推 (`mod_interpolation`)
- 派生类型：`type(spline_1d_t)`（支持 `BC_NATURAL` 与 `BC_CLAMPED`）。
- `spline_1d_init(x, y, spline, bc_type, yp_0, yp_n)`: 内部采用三对角追赶法（Thomas 算法）构建自然三次样条（Natural BC, $y''=0$）或固定导数三次样条（Clamped BC）。
- `spline_1d_eval(spline, x_eval)`: 样条内插求值，自动二分查找区间，无外部矩阵求逆依赖。
- `spline_1d_deriv(spline, x_eval)` / `spline_1d_deriv2(spline, x_eval)`: 解析一阶与二阶导数连续计算，误差达到机器浮点极限，直接用于势能面受力与曲率分析。
- `potential_extrapolate_1d(x_eval, spline, r_min, r_max, a_rep, b_rep, c6_disp, v_inf)`: 科学势能面专用接合外推器，短程指数排斥 $A e^{-B R}$，长程范德华 $V_\infty - C_6/R^6$ 渐近平滑接合。

### 16. 光碎片动力学与通量谱分析 (`mod_photofragment_flux`)
- `calc_autocorrelation(psi_0, psi_t, dx)`: 波包初态与含时态重叠自相关函数 $C(t) = \langle\psi(0)|\psi(t)\rangle$。
- `heller_absorption_spectrum(t_arr, c_t, dt, gamma_damp, omega_arr, spectrum)`: 基于 Heller 理论的连续吸收截面光波谱 $\sigma(\omega) \propto \omega \text{Re}\int_0^\infty C(t) e^{i(\omega+E_0)t - \gamma t} dt$。
- `photofragment_energy_amplitude(t_arr, psi_at_r_det, dt, e_arr, amp_e)`: 渐近监测面 $R_{det}$ 处含时散射波的时间-能量傅里叶散射振幅 $A(E) = \frac{1}{\sqrt{2\pi}}\int \psi(R_{det}, t) e^{i E t} dt$。
- `fragment_kinetic_energy_release(e_photon, v_asymptote, e_bound, ker_spectrum, n_pts)`: 光解离碎片动能释放谱（KER）与能量守恒分析。
- `photofragment_branching_ratio(flux_channels, n_channels, ratios)`: 多通道渐近概率流积分与光化学反应分支比。

### 17. 开放量子系统与 Lindblad 耗散主方程 (`mod_open_quantum`)
- 派生类型：`type(lindblad_system_t)`。
- `lindblad_init(n_levels, n_channels, sys)`: 初始化开放量子系统密度矩阵维度与耗散通道。
- `lindblad_add_decay_channel(sys, i_from, j_to, rate_gamma)`: 添加自发跃迁弛豫算符 $L = \sqrt{\gamma} |j\rangle\langle i|$。
- `lindblad_add_dephasing_channel(sys, level_idx, rate_gamma_d)`: 添加纯退相位跃迁算符 $L = \sqrt{\gamma_d} |i\rangle\langle i|$。
- `lindblad_rhs(rho, h_eff, sys, drho_dt)`: 计算密度矩阵主方程导数 $\frac{d\rho}{dt} = -i[H, \rho] + \sum_k \gamma_k (L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\})$。
- `propagate_lindblad_rk4(rho, h_eff, sys, dt)`: 4 阶保迹保埃尔米特 Runge-Kutta 密度矩阵时域推进器。
- `density_matrix_purity(rho)`: 量子态纯度 $\text{Tr}(\rho^2)$（纯态为 1，最大混合态为 $1/N$）。
- `von_neumann_entropy(rho)`: 冯·诺依曼量子信息熵 $S = -\text{Tr}(\rho \ln \rho)$。
- `quantum_coherence_l1(rho)`: 全局 $l_1$-范数量子相干度 $C_{l_1}(\rho) = \sum_{i \ne j} |\rho_{ij}|$。

### 18. 量子最优控制理论 Krotov 算法 (`mod_optimal_control`)
- 派生类型：`type(oct_config_t)`。
- `oct_config_init(n_steps, dt, alpha_penalty, max_iter, tol_fidelity, cfg)`: 配置控制时域、场强惩罚权重 $\alpha_0$、收敛阈值。
- `state_transfer_fidelity(psi_final, psi_target)`: 目标态转移保真度 $F = |\langle\psi_{target}|\psi(T)\rangle|^2$。
- `oct_shape_function(t, t_total, shape_type)`: 开关整形约束函数 $S(t) = \sin^2(\pi t / T)$，确保激光场在脉冲起止点平滑归零。
- `oct_krotov_step(h0, mu, psi_forward, chi_backward, field_old, field_new, cfg, delta_j)`: 单步执行 Krotov 正向-反向共轭态交替推进并原位更新控制激光场 $\Delta E(t) = -\frac{S(t)}{\alpha_0} \text{Im}\langle\chi(t)|\mu|\psi(t)\rangle$。
- `oct_optimize_pulse(h0, mu, psi_init, psi_target, cfg, field_opt, final_fidelity, stat)`: 高层封装端到端激光脉冲自动优化循环。

### 19. 非含时散射理论与超冷碰撞 (`mod_ti_scattering`)
- 派生类型：`type(scattering_state_t)`, `type(ere_result_t)`, `type(resonance_info_t)`。
- `riccati_bessel_neumann(l, x, jl, nl, d_jl, d_nl)`: 任意轨道角动量 $l$ 的 Riccati-Bessel 函数 $\hat{j}_l(x) = x j_l(x)$ 与 Riccati-Neumann 函数 $\hat{n}_l(x) = x n_l(x)$ 及其解析导数计算，机器精度满足 Wronskian 恒等式 $W[\hat{j}, \hat{n}] = 1.0$。
- `calc_scattering_length_numerov(r_grid, v_pot, mass, a_s, u_zero, stat)`: 零能 Numerov 算法，自包含外推 $u(r) \to C(r - a_s)$ 提取 s-波散射长度 $a_s = r_N - u(r_N)/u'(r_N)$。
- `calc_scattering_length_logder(r_grid, v_pot, mass, a_s, stat)`: Johnson 对数导数比值法，彻底免疫深吸引阱区波函数在经典禁区指数上溢，极低温散射长度黄金标准算法。
- `calc_phase_shift_single_l(r_grid, v_pot, mass, energy, l, delta, k_mat, s_mat, t_mat, stat)`: 有限正能量定态薛定谔方程积分与渐近匹配，提取分波相移 $\delta_l$、反应矩阵 $K_l = \tan\delta_l$、幺正散射矩阵 $S_l = e^{2i\delta_l}$ 与跃迁矩阵 $T_l = S_l - 1$。
- `calc_scattering_wavefunction_ti(r_grid, v_pot, mass, energy, l, norm_type, u_wf, phase_shift, stat)`: **非含时连续谱散射能量本征波函数 $u_{l, E}(r)$ 求解器**。支持三种物理归一化规范：`NORM_ENERGY`（$\delta(E-E')$ 能量归一化，渐近振幅 $\sqrt{\frac{2\mu}{\pi \hbar^2 k}}$）、`NORM_MOMENTUM`（$\delta(k-k')$ 动量归一化，渐近振幅 $\sqrt{2/\pi}$）以及 `NORM_UNIT_AMPLITUDE`（驻波单位振幅 1.0），精确匹配外边界 Riccati 函数提取相移并确保波函数全局相位严格对齐。
- `calc_partial_wave_cross_sections(r_grid, v_pot, mass, energy, l_max, delta_arr, sigma_part, sigma_tot, stat)`: 分波弹性散射截面 $\sigma_l = \frac{4\pi}{k^2}(2l+1)\sin^2\delta_l$ 与总截面 $\sigma_{tot}$。
- `optical_theorem_cross_section(k_wave, delta_arr, l_max)`: 光学定理自洽校验 $\sigma_{optical} = \frac{4\pi}{k} \text{Im}[f(0)]$。
- `calc_differential_cross_section(energy, mass, delta_arr, l_max, theta_grid, dsigma_domega)`: 角度分辨微分散射截面 $\frac{d\sigma}{d\Omega}(\theta) = |f(\theta)|^2$（Legendre 级数展开）。
- `calc_differential_cross_section_identical(energy, mass, delta_arr, l_max, theta_grid, particle_stat, dsigma_domega)`: 考虑全同粒子量子统计干涉的微分散射截面（`STAT_DISTINGUISHABLE`, `STAT_IDENTICAL_BOSON` 在 $\pi/2$ 处干涉增强 4 倍, `STAT_IDENTICAL_FERMION` 在 $\pi/2$ 处奇宇称严格相消归零, `STAT_FERMION_UNPOLARIZED` 自旋统计混合）。
- `calc_transport_cross_sections(energy, mass, delta_arr, l_max, sigma_m, sigma_v)`: 输运截面计算，包含动量传输截面 $\sigma_m = \int (1-\cos\theta) d\sigma$ 与粘滞截面 $\sigma_v = \int \sin^2\theta d\sigma$。
- `calc_differential_legendre_expansion(theta_grid, dsigma_domega, n_theta, k_max, a_k, a_fb)`: 微分散射截面各向异性勒让德展开多极矩 $A_K$ 与前后散射非对称度参数 $A_{FB} = (\sigma_F - \sigma_B)/(\sigma_F + \sigma_B)$。
- `calc_cross_section_spectrum(r_grid, v_pot, mass, energy_grid, n_e, l_max, sigma_tot_spectrum, stat)`: 宽能量范围散射截面与分波相移连续能谱扫描。
- `calc_generalized_cross_sections(k_wave, s_matrix_diag, l_max, sigma_el, sigma_inel, sigma_tot)`: 广义吸收复势弹性截面、非弹性吸收截面与光学总截面。
- `fit_effective_range_expansion(r_grid, v_pot, mass, k_list, n_k, a_s, r_0, stat)`: 超低能区有效力程展开 $k\cot\delta_0(k) = -1/a_s + \frac{1}{2} r_0 k^2$ 最小二乘拟合。
- `van_der_waals_mean_length(mass, c6_au)` / `gribakin_flambaum_length(mass, c6_au, phase_phi)`: 范德华长程色散平均散射长度 $\bar{a} \approx 0.4779888 (2\mu C_6)^{1/4}$ 与 Gribakin-Flambaum 半经典解析散射长度。
- `analyze_shape_resonance(energy_grid, delta_grid, n_pts, hbar, res_info, stat)`: 形状共振 Wigner 散射时延 $\tau(E) = 2\hbar \frac{d\delta}{dE}$ 峰值追踪与 Breit-Wigner 参数（共振能量 $E_R$、线宽 $\Gamma$、准束缚态寿命）提取。
- `calc_coupled_channel_smatrix_2x2(r_grid, v11, v22, v12, mass, total_energy, delta_e, s_matrix, inelastic_prob, stat)`: 双通道非绝热耦合密耦定态散射矩阵求解器，基于 Cayley 变换构建 $2 \times 2$ 严格幺正 $\mathbf{S}$ 矩阵并计算非弹性转移几率 $P_{1\to 2} = |S_{12}|^2$。
- `calc_multichannel_close_coupling_logder(r_grid, v_mat, mass, total_energy, thresholds, l_channels, res, stat)`: **通用任意 $N$ 通道定态密耦求解器 (Johnson Matrix Log-Derivative Method)**。支持开通道（$E > E_i$）与闭通道（$E \le E_i$）任意混合系统，基于 Feshbach 投影 / Schur 补消除闭通道数值发散，构造通量归一化实对称反应矩阵 $\mathbf{K}_{oo}$ 与绝对幺正散射矩阵 $\mathbf{S}_{oo}$（机器精度满足 $\mathbf{S}^\dagger \mathbf{S} = \mathbf{I}$），直接输出各态-态跃迁几率矩阵 $P_{i\to j} = |S_{ij}|^2$、部分弹性与非弹性截面 $\sigma_{i\to j}$ 及特征相移和 $\delta_{sum}$。
- `calc_feshbach_resonance_scan(r_grid, v_mat, mass, energy_grid, n_energies, thresholds, l_channels, s_wave_length, eigenphase_sums, stat)`: 跨 Feshbach 共振能量扫描，追踪闭通道准束缚态引起的开通道散射长度极点发散 $a_s(E) \to \pm \infty$ 与相移特征跳变。

### 20. 含时波包散射理论与 S-矩阵 (`mod_td_scattering`)
- 派生类型：`type(td_scattering_channel_t)`。
- `gaussian_wavepacket_1d(x_grid, x0, sigma_x, k0, psi_0)`: 构造空间与动量严格归一化的入射高斯散射波包。
- `gaussian_momentum_amplitude(k, x0, sigma_x, k0)`: 解析动量表象振幅 $g(k)$ 与入射通量权重。
- `accumulate_flux_amplitude(t, psi_at_det, dt, energy_grid, n_energies, hbar, amp_accum)`: 渐近监测面时间-能量傅里叶振幅 $A(E) = \frac{1}{\sqrt{2\pi}} \int_0^T \psi(x_{det}, t) e^{i E t/\hbar} dt$。
- `calculate_td_transmission(energy_grid, n_energies, amp_trans, mass, hbar, x0, sigma_x, k0, t_prob)`: 单次含时波包模拟直接提取连续能域透射几率谱 $T(E) = \frac{\hbar k_E}{\mu} \frac{|A_{trans}(E)|^2}{|g(k_E)|^2}$。
- `calculate_td_smatrix_element(amp_scatter, amp_free, n_energies, s_matrix_e, phase_shift_e)`: 比对有势碰撞与自由对照演化，提取全能量散射矩阵元 $S(E) = A_{scatter}(E)/A_{free}(E)$ 与散射相移 $\delta(E) = \frac{1}{2}\arg(S(E))$。
- `project_wavepacket_to_smatrix(x_grid, dx, psi_final, mass, hbar, k0, sigma_x, x0, energy_grid, n_energies, t_prob, r_prob)`: Möller 动量表象渐近投影法，末态波包动量空间分解提取连续态透射与反射概率。
- `multichannel_td_smatrix_elements(...)`: 多通道含时概率通量提取非绝热碰撞非弹性 S-矩阵元 $|S_{ij}(E)|^2$。
- `wavepacket_centroid_position(x_grid, dx, psi)` / `wavepacket_wigner_delay(...)`: 波包质心轨迹追踪与含时 Wigner 散射时延 $\tau_W(E) = 2\hbar \frac{d\delta}{dE}$。
- `calculate_td_differential_cross_section_2d(x_grid, y_grid, psi_2d, mass, hbar, theta_grid, dsigma_dtheta)`: 二维连续态含时波包散射角分布微分散射截面 $\frac{d\sigma}{d\theta}(\theta)$。
- `accumulate_wavefunction_spectral_projection(psi_t, t, dt, energy, hbar, psi_energy_accum)`: **含时动力学全空间谱投影原位累积器**。在波包推进主循环中无缝累积时间-能量半傅里叶变换 $\int_0^T \Psi(x, t) e^{iEt/\hbar} dt$。
- `extract_td_scattering_wavefunction(x_grid, psi_energy_accum, energy, mass, hbar, x0, sigma_x, k0, psi_energy_norm, stat)`: **含时谱投影连续谱能量本征波函数提取器**。严格消除入射波包动量权重 $g(k_E)$ 与态密度变换因子，直接恢复满足严格 $\delta(E-E')$ 能量正交归一化的定态连续能量本征函数 $\psi_E(x)$。

---

## 📖 详细配置手册

本算法库配备了详尽的配置与环境搭建指南：
👉 **[CONFIG_GUIDE.md](file:///Users/lihao/Library/CloudStorage/SynologyDrive-aecho/Codes/Fortran/GeneralModule/CONFIG_GUIDE.md)**
- **全平台编译器配置**：GCC/gfortran (9~15)、Intel oneAPI (ifx/ifort)、macOS Xcode 许可绕过说明。
- **构建系统深度指南**：fpm, CMake, 通用 Makefile, 纯命令行打包实战。
- **算法参数配置全典**：`pulse_config_t`, `dvr_1d_t`, `absorbing_boundary_t`, `atom_config_t` 取值范围与物理单位换算。
- **AI 编程与提示词工程模板**：提供直接粘贴给 AI 的 System Prompts 与常用物理场景任务 Prompts。

---

## 🚀 快速上手与集成指南

### 1. 在其他 Fortran 程序中调用

仅需一行代码引入顶层聚合模块即可访问全部功能：

```fortran
program my_simulation
    use general_module
    implicit none

    type(dvr_1d_t) :: dvr
    real(dp) :: dt, t_val
    type(pulse_config_t) :: laser
    type(atom_config_t) :: atom

    ! 1. 简便的单位转换
    dt = to_au(0.1_dp, "fs")

    ! 2. 激光脉冲与靶原子配置
    call get_atom_config("Ar", atom)
    laser%shape_type = PULSE_GAUSSIAN
    laser%field_peak = 0.05_dp
    laser%duration = to_au(30.0_dp, "fs")

    ! 3. 网格初始化与计算
    call dvr_sinc_init(-10.0_dp, 10.0_dp, 256, 1.0_dp, dvr)

    print *, "Current e-field at t=0:", pulse_electric_field(0.0_dp, laser)
    print *, "HHG Cutoff Energy (eV):", hhg_cutoff_energy(atom%ip_au, laser%field_peak, 0.057_dp) * AU2EV
end program my_simulation
```

### 2. 编译与链接方式

#### 方式 A: 直接使用 gfortran 编译
```bash
# 1. 编译 GeneralModule 模块为静态库
cd GeneralModule/src
DEVELOPER_DIR=/Library/Developer/CommandLineTools gfortran -O3 -fPIC -c *.f90
DEVELOPER_DIR=/Library/Developer/CommandLineTools ar rcs libgeneral_module.a *.o

# 2. 编译并链接您的程序
DEVELOPER_DIR=/Library/Developer/CommandLineTools gfortran -O3 -I/path/to/GeneralModule/src my_code.f90 /path/to/GeneralModule/src/libgeneral_module.a -o my_code
```

#### 方式 B: 使用 Fortran Package Manager (`fpm`)
在您的项目 `fpm.toml` 中声明依赖：
```toml
[dependencies]
general_module = { path = "/path/to/GeneralModule" }
```
然后直接执行：
```bash
fpm build
```

#### 方式 C: 使用 CMake
在您的 `CMakeLists.txt` 中引入：
```cmake
add_subdirectory(/path/to/GeneralModule GeneralModule_build)
target_link_libraries(my_executable PRIVATE GeneralModule_static)
```

---

## 🧪 自动化测试套件

算法库内建完备的单元测试，覆盖物理常数往返转换、角动量耦合、FGH 谐振子能级、激光脉冲包络、Bloch / Split-Operator 模长守恒、强场原子模型与多态非绝热耦合、转振态激光调控、三次样条插值与外推、自相关吸收谱与碎片 KER 分支比、Lindblad 耗散主方程与 Krotov 最优控制、非含时散射长度与 S-矩阵、含时波包散射连续态透射谱：

```bash
cd GeneralModule/tests
chmod +x run_all_tests.sh
./run_all_tests.sh
```

**测试结果清单**：
```text
================================================================
          Running GeneralModule Test Suite Suite                
================================================================
Constants Tests:             12 / 12 PASSED
Special Function Tests:      10 / 10 PASSED
DVR Grid Tests:               7 /  7 PASSED
Laser Pulse Tests:            9 /  9 PASSED
Propagator Tests:             8 /  8 PASSED
Extended Atomic Tests:       14 / 14 PASSED
Rovibrational Control Tests:  8 /  8 PASSED
Interpolation Tests:          7 /  7 PASSED
Photofragment & Flux Tests:   6 /  6 PASSED
Open Quantum & OCT Tests:    10 / 10 PASSED
TI Scattering Tests:         16 / 16 PASSED
TD Scattering Tests:          9 /  9 PASSED
----------------------------------------------------------------
ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass, 117/117 断言通过)
================================================================
```

---

## 📊 典型物理算例 (Examples)

位于 `GeneralModule/examples/`，一键编译运行全部 6 大物理算例：
```bash
cd GeneralModule/examples
chmod +x build_examples.sh
./build_examples.sh
```

1. **`ex01_fgh_diatomic_bound_states.f90`**
   - 求解双原子分子（HF）Morse 势能面束缚态，与解析解对比精度优于 $10^{-8}$ a.u.（差值小于 $0.0001\text{ cm}^{-1}$）。
2. **`ex02_pulse_synthesis.f90`**
   - 演示高斯啁啾脉冲、双色合成光场与 THz 脉冲串的时域电场与 FFT 频谱转换。
3. **`ex03_split_operator_1d.f90`**
   - 演示一维高斯波包穿过势垒发生透射与反射，进入边界复吸收势（CAP）平滑衰减，无人工边界反射伪影。
4. **`ex04_field_free_orientation.f90`**
   - 模拟极性分子 CO 刚体转动波包的非绝热无场取向与排列，演示 $30\text{ K}$ 玻尔兹曼转动系综热平均。
5. **`ex05_hhg_lewenstein_spectrum.f90`**
   - 强场激光作用下氩原子的阿秒高次谐波发射仿真，通过 Lewenstein 强场近似生成偶极加速度并精确展现谐波平台与截止能。
6. **`ex06_two_state_nonadiabatic.f90`**
   - 模拟核波包穿过经典 Tully 避差交叉势能面（Avoided Crossing），计算无绝热跃迁与各通道末态布居比。

---

## 🐍 Python 辅助分析套件 (`pygenmod`)

提供轻量 Python 库，用于数据交互、前处理计算与出版级可视化：

```bash
cd GeneralModule/python
python3 test_pygenmod.py   # 运行 9 大单元测试 (100% Pass)
```

### 1. 超冷散射长度与零能波函数渐近线可视化 (`scattering.py`)
计算任意相互作用势的零能波函数 $u(r)$，解析并可视化渐近线 $C(r - a_s)$ 与 $r$ 轴截距所确定的散射长度 $a_s$：
```python
from pygenmod import calc_scattering_length_numerov, plot_scattering_length_wavefunction
import numpy as np

r = np.linspace(0.01, 12.0, 600)
v_pot = np.where(r <= 2.0, -1.0, 0.0)  # 吸引方势阱
a_s, u_wf = calc_scattering_length_numerov(r, v_pot, mass=1.0)
plot_scattering_length_wavefunction(r, v_pot, u_wf, a_s, filename="scattering_length.png")
```

### 2. 微分散射截面与全同粒子干涉绘图 (`scattering.py`)
计算可区分粒子、全同玻色子、极化费米子角分布，极坐标与直角坐标双面板出版级成图：
```python
from pygenmod import calc_differential_cross_section_identical, plot_differential_cross_sections
import numpy as np

theta = np.linspace(0.0, np.pi, 181)
delta = np.array([0.45, 0.15, 0.05])  # l=0, 1, 2 分波相移
ds_dist = calc_differential_cross_section_identical(0.05, 1.0, delta, theta, 'distinguishable')
ds_boson = calc_differential_cross_section_identical(0.05, 1.0, delta, theta, 'boson')
ds_fermion = calc_differential_cross_section_identical(0.05, 1.0, delta, theta, 'fermion')

plot_differential_cross_sections(theta, {
    "Distinguishable": ds_dist,
    "Identical Bosons": ds_boson,
    "Polarized Fermions": ds_fermion
}, filename="result_differential_cross_section.png")
```

### 3. 转振态激光调控动力学一键绘图 (`plot_rovibrational_dynamics.py`)
结合 Fortran 计算生成的 `test_rovibrational_dynamics.dat`，一键绘制 4 面板 300 DPI 矢量级高清科研图（激光电场、各振动态布居演化、转振精细结构布居、全系综总几率守恒验证）：
```bash
python3 python/plot_rovibrational_dynamics.py
# 生成高分辨率图表: python/result_rovibrational_dynamics.png
```

### 4. 多通道密耦 S-矩阵热图与 Feshbach 共振绘图 (`scattering.py`)
计算多通道定态跃迁几率热图矩阵与 Feshbach 共振极点扫描：
```python
from pygenmod import calc_multichannel_close_coupling, plot_multichannel_smatrix, plot_feshbach_resonance
import numpy as np

# 1. 求解多通道耦合与状态转移几率矩阵
res = calc_multichannel_close_coupling(r_grid, v_mat, mass=1.0, total_energy=0.5, thresholds=np.array([0.0, 0.2]))
plot_multichannel_smatrix(res['prob_matrix'], ["Ch 1", "Ch 2"], filename="multichannel_smatrix.png")

# 2. 绘制 Feshbach 共振散射长度发散与特征相移跃升
plot_feshbach_resonance(energy_grid, scattering_lengths, eigenphase_sums, filename="feshbach_resonance.png")
```

### 5. Sinc-DVR 束缚态求解与波函数绘图示例：
```python
from pygenmod import dvr_sinc_init, fgh_solve_bound_states, plot_wavefunctions

# 1. 初始化 Sinc-DVR
dvr = dvr_sinc_init(x_min=-8.0, x_max=8.0, n_points=256, mass=1.0)

# 2. 定义势能曲线 (如谐振子势)
v_pot = 0.5 * (dvr.x ** 2)

# 3. 求解本征态
eig_vals, wavefuncs = fgh_solve_bound_states(dvr, v_pot)

# 4. 生成发表级对比图
plot_wavefunctions(dvr.x, v_pot, eig_vals, wavefuncs, n_states=4, filename="harmonic_oscillator.png")
```

---

## 🤖 AI 编程与 Agent 使用指南

若您在日常开发中配合 AI 编程助手（如 Antigravity, Claude, Copilot）工作，可以直接指示 AI：

> *"请使用 `general_module` 编写一个 Fortran 程序，调用 `dvr_sinc_init` 和 `fgh_solve_bound_states` 求解一维双阱势 $V(x) = x^4 - 2x^2$ 的前 4 个本征态与能级劈裂，并用 `to_au` 统一输入物理单位。"*

AI 将直接利用自包含的类型定义和确定的函数签名生成无错代码，无需额外的繁琐胶水代码或外部环境配置。
