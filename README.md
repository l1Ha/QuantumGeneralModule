# GeneralModule: 现代量子动力学通用算法库 (Fortran 2008 / Python)

`GeneralModule` 是一个面向超快强场物理、分子光物理与量子动力学模拟的现代化通用科学计算算法库。该库遵循严格的 **Fortran 2008 规范**，具备高数值精度、零外部动态库强依赖、模块化架构与出色的 AI Agent 友好性。

---

## 🌟 核心特性

1. **零外部库依赖 (Zero External Dependencies)**
   - 内部集成高精度 Householder QL 实对称矩阵本征求解器与 Cooley-Tukey 1D/2D 快速傅里叶变换（FFT）。
   - 无需强制链接外部 LAPACK/BLAS 或 FFTW，开箱即用，支持跨平台一键编译。
2. **现代 Fortran 2008 标准设计**
   - 统一强类型参数定义（`real(dp) => real64`）。
   - 纯函数（`pure function`）与显式 `intent(in/out/inout)` 契约，杜绝隐式全局变量副作用。
3. **AI 友好型结构化接口 (AI-Friendly)**
   - 算法模块支持统一顶层聚合入口：`use general_module`。
   - 参数配置采用清晰的派生类型（Derived Types，如 `pulse_config_t`, `absorbing_boundary_t`, `dvr_1d_t`），自解释、低耦合、便于大语言模型精确构造与调用。
4. **全链路双语生态支持**
   - 附带标准 Python 包装分析包 `pygenmod`，无缝衔接参数预计算、波包可视化与发表级绘图。

---

## 📂 目录结构与架构

```text
GeneralModule/
├── CONFIG_GUIDE.md                # 📖 详尽配置与部署指南（全编译平台、参数全典、AI提示词模板）
├── fpm.toml                       # Fortran Package Manager 配置文件
├── CMakeLists.txt                 # CMake 跨平台构建系统
├── README.md                      # 本文档
├── .gitignore                     # Git 忽略规则
├── src/                           # 核心 Fortran 源代码 (14 核心模块 + 1 聚合入口)
│   ├── mod_constants.f90          # 1. 物理常数与各单位 a.u. 双向转换
│   ├── mod_special_functions.f90  # 2. 勒让德、Wigner 3j、CG、转动偶极/取向矩阵元
│   ├── mod_linear_algebra.f90     # 3. 对称矩阵本征求解 (EISPACK TRED2/TQL2)、1D/2D FFT
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
│   └── general_module.f90         # 顶层聚合入口模块 (use general_module)
├── tests/                         # 自动化单元测试套件 (7 个套件，100% 全部通过)
│   ├── test_constants.f90
│   ├── test_special_functions.f90
│   ├── test_dvr_grid.f90
│   ├── test_laser_pulse.f90
│   ├── test_propagators.f90
│   ├── test_atomic_hhg.f90
│   ├── test_laser_rovibrational_control.f90 # 激光调控分子转振态布居转移综合测试
│   └── run_all_tests.sh           # 自动化测试运行脚本 (100% Pass, 69/69 断言)
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
    ├── test_pygenmod.py           # Python 单元测试 (100% Pass)
    └── pygenmod/
        ├── __init__.py
        ├── constants.py
        ├── pulse.py
        ├── dvr.py
        ├── coulomb.py
        ├── hhg.py
        ├── multistate.py
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

算法库内建完备的单元测试，覆盖物理常数往返转换、角动量耦合、FGH 谐振子能级、激光脉冲包络、Bloch / Split-Operator 模长守恒、强场原子模型与多态非绝热耦合：

```bash
cd GeneralModule/tests
chmod +x run_all_tests.sh
./run_all_tests.sh
```

**测试结果示例**：
```text
Constants Tests:             10 / 10 PASSED
Special Function Tests:      12 / 12 PASSED
DVR Grid Tests:               8 /  8 PASSED
Laser Pulse Tests:            9 /  9 PASSED
Propagator Tests:             8 /  8 PASSED
Extended Atomic Tests:       14 / 14 PASSED
Rovibrational Control Tests:  8 /  8 PASSED
---------------------------------------------
ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass, 69/69 tests)
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

提供轻量 Python 库，用于数据交互与前后期处理：

```bash
cd GeneralModule/python
python3 test_pygenmod.py
```

### Python 调用示例：
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
