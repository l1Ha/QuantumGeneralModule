# GeneralModule: 现代量子动力学通用算法库 (Fortran 2008 / Python)

[![CI](https://github.com/l1Ha/QuantumGeneralModule/actions/workflows/ci.yml/badge.svg)](https://github.com/l1Ha/QuantumGeneralModule/actions/workflows/ci.yml)
[![Fortran 2008](https://img.shields.io/badge/Fortran-2008-734f96.svg)](https://fortran-lang.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Tests: 335/335 Pass](https://img.shields.io/badge/Tests-335%2F335%20Pass%20(100%25)-brightgreen.svg)](tests/)
[![Literature: 35 Topics](https://img.shields.io/badge/Literature-35%20Topics%20(PRL%2FPRA%2FRMP%2FScience%2FNature)-blue.svg)](LITERATURE.md)

`GeneralModule` 是一个面向超快强场物理、分子光物理、前沿量子散射与表面动力学模拟的现代化通用科学计算算法库。该库遵循严格的 **Fortran 2008 规范**，具备高数值精度、零外部动态库强依赖、模块化架构与出色的 AI Agent 友好性。

---

## 🌟 核心特性

1. **零外部库依赖 (Zero External Dependencies)**
   - 内部集成高精度 Householder QL 实对称矩阵本征求解器、Gauss-Jordan 全主元实/复方阵求逆与 Cooley-Tukey 1D/2D 快速傅里叶变换（FFT）。
   - 纯 Fortran 自包含样条插值、Lindblad 主方程、Krotov 最优控制、通用多通道定态密耦（Johnson Log-Derivative）、外场多基组散射、各向异性偶极耦合、超冷光缔合速率、少体 Efimov 物理、低维光晶格 CIR、自电离 Fano/CCR、交叉电磁场、三原子反应 PES、旋量 BEC 自旋动力学、三原子超球面反应动力学、偶极量子液滴 LHY、强场非顺序双电离 (NSDI)、磁/光 Feshbach 束缚态、阿秒瞬态吸收光谱 (ATAS)、双色反向圆偏振 PECD、超冷极性分子偶极遮蔽、里德堡原子阻塞、2D 表面量子散射与 SAR、气-固催化 Eley-Rideal 反应、金属表面非绝热电子摩擦 (GLE)、掠入射快原子衍射 (GIFAD)、冷离子-中性原子杂化极化散射、Tully 最少开关表面跳跃 (FSSH)、强场分子定向与超转子动力学、光晶格 Bose-Hubbard 映射、反应路径哈密顿量 (RPH/CVT)、相对论径向狄拉克方程以及共振非弹性 X 射线散射 (RIXS) 求解器，无需强制链接外部 LAPACK/BLAS 或 FFTW，开箱即用。
2. **现代 Fortran 2008 标准设计**
   - 统一强类型参数定义（`real(dp) => real64`）。
   - 纯函数（`pure function`）与显式 `intent(in/out/inout)` 契约，杜绝隐式全局变量副作用。
3. **AI 友好型结构化接口 (AI-Friendly)**
   - 算法模块支持统一顶层聚合入口：`use general_module`。
   - 参数配置采用清晰的派生类型（Derived Types，如 `pulse_config_t`, `surface_lattice_t`, `ion_atom_system_t`, `tully_model_t`, `rotor_molecule_t`, `optical_lattice_t`, `rph_path_t`, `dirac_state_t`, `rixs_system_t` 等），自解释、低耦合、便于大语言模型精确构造与调用。
4. **全链路双语生态支持**
   - 附带标准 Python 伴侣分析包 `pygenmod`，无缝衔接参数预计算、波包与散射长度可视化、Breit-Rabi 能级图、发表级绘图（含一键动力学出图流水线 `plot_rovibrational_dynamics.py`）。
5. **全自动 CI/CD 持续集成**
   - 内置 GitHub Actions 跨平台持续集成（Ubuntu / macOS），全自动化执行 40 大测试套件（335 个单元断言 100% 通过）与 35 大物理应用工程算例。

---

## 📂 目录结构与架构

```text
GeneralModule/
├── .github/workflows/ci.yml       # 🚀 GitHub Actions 跨平台 CI 持续集成工作流
├── CONFIG_GUIDE.md                # 📖 详尽配置与部署指南（全编译平台、参数全典、AI提示词模板）
├── LITERATURE.md                  # 📚 科学文献典藏与理论映射全典 (PRL/PRA/RMP/JCP 权威论文与代码映射)
├── fpm.toml                       # Fortran Package Manager 配置文件
├── CMakeLists.txt                 # CMake 跨平台构建系统
├── README.md                      # 本文档
├── .gitignore                     # Git 忽略规则
├── src/                           # 核心 Fortran 源代码 (48 核心模块 + 1 聚合入口)
│   ├── mod_constants.f90          # 1. 物理常数与各单位 a.u. 双向转换
│   ├── mod_special_functions.f90  # 2. 勒让德、Wigner 3j/6j/9j、CG半整数代数、转动偶极/取向矩阵元
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
│   ├── mod_ti_scattering.f90      # 19. 非含时散射理论、零能 Numerov/Log-Derivative 散射长度、相移、分段网格与多通道密耦
│   ├── mod_td_scattering.f90      # 20. 含时波包散射理论、连续态能量通量透射率 T(E)、含时 S 矩阵提取与 Möller 动量投影
│   ├── mod_field_scattering.f90   # 21. 外加电磁场超冷散射、四大经典基组严格幺正变换、Breit-Rabi本征态与磁Feshbach共振扫描
│   ├── mod_dipolar_scattering.f90 # 22. 各向异性磁偶极/电偶极散射、自旋弛豫截面/热速率、极性分子Stark诱导偶极与耦合势
│   ├── mod_photoassociation.f90   # 23. 超冷光缔合谱学、自由-束缚态Franck-Condon重叠积分、受激线宽与热平均光缔合速率
│   ├── mod_three_body_recombination.f90 # 24. 超冷三体复合碰撞与 Efimov 少体物理、普适 K3 速率与离散标度不变性
│   ├── mod_confined_scattering.f90 # 25. 低维光晶格受限量子散射与约束诱导共振 CIR、有效 1D 相互作用与 Tonks 气体
│   ├── mod_autoionization_fano.f90 # 26. 自电离体系、Fano 组态相互作用线型与复坐标旋转法 CCR 共振寿命
│   ├── mod_crossed_field_scattering.f90 # 27. 交叉静电磁场 E x B 转振-自旋动力学、宇称破缺与倾角扫描
│   ├── mod_triatomic_geometry.f90 # 28. 三原子反应散射 Jacobi 坐标、LEPS 反应势能面与锥形交叉 Berry 几何相位
│   ├── mod_spinor_bec.f90         # 29. 超冷旋量玻色爱因斯坦凝聚 F=1 多体自旋动力学、铁磁/极性相与相干自旋混合
│   ├── mod_hyperspherical_reactive.f90 # 30. 三原子超球面反应动力学、Delves 质量标度、反应偏角、Eckart 隧穿与反应速率 k(T)
│   ├── mod_dipolar_droplets_lhy.f90    # 31. 超冷偶极量子液滴、Pelster-Lima 涨落积分 Q5、自束缚平衡密度 n0 与 eGPE 能量密度
│   ├── mod_strong_field_nsdi.f90       # 32. 强场非顺序双电离、Corkum 三步模型、3.17 Up 截断、2D 动量关联谱与双电离膝盖结构
│   ├── mod_feshbach_bound_states.f90   # 33. 磁与光 Feshbach 共振、Coupled-channel 弱束缚分子态 Eb(B)、闭通道权重与光致损耗
│   ├── mod_attosecond_transient_absorption.f90 # 34. 阿秒瞬态吸收光谱 (ATAS)、相位微扰模型 (PPM)、动态 Fano 参数与光诱导态
│   ├── mod_bicircular_pecd.f90         # 35. 双色反向旋转圆偏振场、C3 动力学对称性、手性四面体势与光电子圆二色性 (PECD)
│   ├── mod_ultracold_reaction_shielding.f90    # 36. 超冷极性分子反应动力学、微波/静电偶极遮蔽势垒、WKB 隧穿抑制与 gamma > 100
│   ├── mod_rydberg_blockade.f90        # 37. 里德堡原子阻塞、C6 标度律、双原子动力学与 1D 阵列 PXP 量子多体疤痕
│   ├── mod_surface_scattering.f90      # 38. 表面量子散射、2D 晶格相干衍射、硬波纹表面 (HCS) 与选择性吸附共振 (SAR)
│   ├── mod_surface_reaction_er.f90     # 39. 气-固表面催化与 Eley-Rideal 提取机理、超热放热能量分配与振动反转
│   ├── mod_surface_electronic_friction.f90 # 40. 金属表面非绝热动力学、局域密度摩擦 (LDFA)、广义朗之万 (GLE) 与电子-空穴对耗散
│   ├── mod_grazing_fast_atom_diffraction.f90 # 41. 掠入射快原子表面量子衍射 (GIFAD)、轴向沟道快慢解耦与亚皮米波纹反演
│   ├── mod_ion_atom_scattering.f90     # 42. 冷离子-中性原子杂化散射、1/r^4 极化势、Langevin 截面与 Paul 阱微运动致热
│   ├── mod_surface_hopping_fssh.f90    # 43. Tully 最少开关表面跳跃 (FSSH)、Velocity Verlet 核推进、能量守恒动量重标度与 NACV
│   ├── mod_molecular_alignment.f90     # 44. 强场非绝热分子定向/取向、转动复苏序参量与光学离心机超转子动力学
│   ├── mod_optical_lattice_hubbard.f90  # 45. 超冷光晶格、Mathieu 能带结构、Wannier 轨道、Bose-Hubbard 映射与 Bloch 振荡
│   ├── mod_reaction_path_hamiltonian.f90 # 46. 多原子反应路径哈密顿量 (RPH)、内禀反应坐标 (IRC)、变分过渡态理论 (CVT) 与 Eckart 隧穿
│   ├── mod_relativistic_atomic.f90     # 47. 相对论原子结构、径向狄拉克方程、精细结构分裂 (Cs D1/D2) 与核心极化模型势
│   ├── mod_resonant_xray_scattering.f90 # 48. 共振非弹性 X 射线散射 (RIXS)、Kramers-Heisenberg 二阶截面、XAS/XES 与 Huang-Rhys 振动级数
│   └── general_module.f90          # 顶层聚合入口模块 (use general_module)
├── tests/                          # 自动化单元测试套件 (40 个套件，100% 全部通过，335/335 断言)
│   ├── test_constants.f90
│   ├── test_special_functions.f90
│   ├── test_dvr_grid.f90
│   ├── test_laser_pulse.f90
│   ├── test_propagators.f90
│   ├── test_atomic_hhg.f90
│   ├── test_laser_rovibrational_control.f90 # 激光调控分子转振态布居转移综合测试
│   ├── test_interpolation.f90      # 三次样条插值与渐近外推测试
│   ├── test_photofragment_flux.f90 # 自相关吸收谱与碎片 KER 分支比测试
│   ├── test_open_quantum_opt.f90   # Lindblad 耗散退相干与 Krotov 最优控制测试
│   ├── test_ti_scattering.f90      # 非含时散射长度、相移、S矩阵与分段网格多通道测试
│   ├── test_td_scattering.f90      # 含时波包散射透射谱、S矩阵元与 Möller 投影测试
│   ├── test_field_scattering.f90   # 外场四大基组幺正变换、Breit-Rabi解析与数值比对、磁Feshbach共振拟合测试
│   ├── test_dipolar_scattering.f90 # 各向异性偶极张量、自旋弛豫截面、极性分子Stark感应与耦合势测试
│   ├── test_photoassociation.f90   # 自由-束缚重叠积分、受激线宽、热光缔合速率与Raman双光子测试
│   ├── test_three_body_recombination.f90 # Efimov 超越方程根、普适 K3 速率、干涉极小值与有限温度幂律
│   ├── test_confined_scattering.f90       # 准 1D 波导 CIR 极点、结合能、Tonks-Girardeau 与各向异性分裂
│   ├── test_autoionization_fano.f90       # Fano 抗共振零点、自电离寿命与复坐标旋转 CCR 极点
│   ├── test_crossed_field_scattering.f90  # 交叉静电磁场非共线态混合、Stark 定向度与避免交叉能谱
│   ├── test_triatomic_geometry.f90        # Jacobi 坐标可逆映射、LEPS 反应势能面与锥形交叉 Berry 几何相位
│   ├── test_spinor_bec.f90                # F=1 旋量凝聚体铁磁/极性相、保全几率与保磁化强度 RK4 自旋混合
│   ├── test_hyperspherical_reactive.f90   # 三原子反应偏角、Eckart 隧穿传递几率、累积反应几率 N(E) 与热速率 k(T)
│   ├── test_dipolar_droplets_lhy.f90      # 162Dy 偶极长度、Pelster-Lima Q5 积分、自束缚平衡密度与 eGPE 负化学势
│   ├── test_strong_field_nsdi.f90         # 强场 3.17 Up 回碰截止、(e,2e) Lotz 截面、2D 平行动量关联与双电离膝盖结构
│   ├── test_feshbach_bound_states.f90     # 6Li/87Rb 磁 Feshbach 弱束缚态能谱 Eb(B)、闭通道权重 Z(B) 与 OFR 双体损耗率 K2
│   ├── test_attosecond_transient_absorption.f90 # 氦原子 2s2p 自电离、动态 Fano q 与 2D ATAS 时延谱不对称性测试
│   ├── test_bicircular_pecd.f90           # 双色反向圆偏振场 C3 对称、手性四面体不变量 chi 与 PECD 前后发射不对称测试
│   ├── test_ultracold_reaction_shielding.f90     # KRb 微波遮蔽排斥势垒、WKB 隧穿几率与蒸发冷却比值 gamma > 100 测试
│   ├── test_rydberg_blockade.f90          # 87Rb 70S 阻塞半径 Rb、双原子强阻塞抑制与 1D 阵列量子多体疤痕测试
│   ├── test_surface_scattering.f90        # 2D 晶格 HCS 程函衍射幺正性、SAR Fano 线型与 Debye-Waller 声子衰减测试
│   ├── test_surface_reaction_er.f90       # H+H/Cu(111) ER 反应放热量超热分配、振动布居反转与反应截面/速率常数测试
│   ├── test_surface_electronic_friction.f90  # Au(111) 电子摩擦系数空间衰减、GLE 辛步进、电子-空穴对能损与振动寿命测试
│   ├── test_grazing_fast_atom_diffraction.f90 # 快轴向沟道解耦、经典彩虹偏转角、1D 横向 Bragg 衍射谱与亚皮米波纹度反演测试
│   ├── test_ion_atom_scattering.f90       # 极化尺度 R*、Langevin 标度律、MERE 相移与微运动致热率测试
│   ├── test_surface_hopping_fssh.f90      # Tully SAC 能隙、核动能守恒、系综分支比与 Ehrenfest 模范数测试
│   ├── test_molecular_alignment.f90       # 转动复苏周期、瞬态定向峰值、光学离心机加速与离心破键测试
│   ├── test_optical_lattice_hubbard.f90   # Bloch 能带极限、超流-Mott 判据、带隙演化与引力 Bloch 振荡测试
│   ├── test_reaction_path_hamiltonian.f90 # RPH 鞍点势垒、CVT 变分界、Eckart 隧穿与高温极限测试
│   ├── test_relativistic_atomic.f90       # 狄拉克氢 1s 基态、2p 精细结构劈裂、Cs 核心极化与相对论 E1 振子强度测试
│   ├── test_resonant_xray_scattering.f90  # Cu L3 XAS 吸收峰、Kramers-Heisenberg 共振放大、2D RIXS 图谱与 Laguerre 声子级数测试
│   └── run_all_tests.sh            # 自动化测试运行脚本 (100% Pass, 335/335 断言, 40 测试套件)
├── examples/                       # 典型物理应用算例 (35 大完整前沿算例)
│   ├── ex01_fgh_diatomic_bound_states.f90 # 双原子 Morse 势能级与波函数求解
│   ├── ex02_pulse_synthesis.f90            # 啁啾、双色、太赫兹脉冲时频生成
│   ├── ex03_split_operator_1d.f90          # 1D 波包动力学演化与 CAP 吸收边界
│   ├── ex04_field_free_orientation.f90     # 刚体转子无场定向与玻尔兹曼热平均
│   ├── ex05_hhg_lewenstein_spectrum.f90    # 强场阿秒高次谐波发射与半经典截止能
│   ├── ex06_two_state_nonadiabatic.f90     # 双态避差穿越非绝热动力学与分支比
│   ├── ex07_scattering_wavefunctions_ti_td.f90 # 连续态能量本征波函数非含时与含时双向求解对比
│   ├── ex08_ultracold_feshbach_segmented.f90   # 多扇区分段网格超冷磁 Feshbach 共振与散射长度扫描
│   ├── ex09_dipolar_relaxation_scattering.f90  # 各向异性磁偶极自旋弛豫与极性分子外电场 Stark 诱导偶极
│   ├── ex10_photoassociation_spectroscopy.f90  # 超冷光缔合跃迁谱学、受激饱和展宽与双光子 Raman 缔合
│   ├── ex11_three_body_efimov_recombination.f90 # 超冷三体 Efimov 复合速率与干涉极小值/共振峰扫描
│   ├── ex12_confined_cir_scattering.f90         # 准一维光晶格波导中 87Rb 约束诱导共振 CIR 与分子结合能
│   ├── ex13_autoionization_fano_resonance.f90   # 自电离 Fano 不对称吸收谱与复坐标旋转 CCR 寿命提取
│   ├── ex14_spinor_bec_dynamics.f90             # 87Rb 与 23Na 凝聚体宏观自旋混合动力学与铁磁相图
│   ├── ex15_crossed_field_stark_zeeman.f90      # 交叉电磁场中极性顺磁分子 Stark-Zeeman 态混合与空间定向
│   ├── ex16_triatomic_reaction_berry_phase.f90  # 三原子反应路径 LEPS 势能面与锥形交叉 Berry 几何相位
│   ├── ex17_hyperspherical_reaction_rates.f90   # 三原子超球面反应动力学、Eckart 隧穿累积反应几率与正则热速率常数
│   ├── ex18_dipolar_quantum_droplets.f90        # 162Dy 偶极量子液滴自束缚平衡密度、负化学势与气-液滴相变
│   ├── ex19_strong_field_nsdi_recollision.f90   # 800nm 强场电子重碰撞动能 3.17 Up 截断、2D 动量关联谱与双电离膝盖结构
│   ├── ex20_feshbach_molecular_bound_states.f90 # 6Li 磁 Feshbach 晕轮二聚体结合能、闭通道权重与 87Rb 光 Feshbach 损耗
│   ├── ex21_attosecond_transient_absorption.f90 # 氦原子 2s2p 自电离、动态 Fano q、LIS 态与 2D ATAS 谱图
│   ├── ex22_bicircular_pecd_chiral.f90          # 双色反向圆偏振场、手性分子对映体能量分辨 PECD 谱与 2D PAD
│   ├── ex23_ultracold_molecule_shielding.f90    # KRb 超冷极性分子微波遮蔽势垒、隧穿抑制与蒸发冷却速率比 gamma
│   ├── ex24_rydberg_blockade_dynamics.f90       # 87Rb 70S 阻塞半径、二原子双激发抑制与 1D 阵列量子多体疤痕动力学
│   ├── ex25_surface_corrugated_diffraction.f90  # He/LiF(001) 2D 晶格衍射谱、选择性吸附共振 SAR 与 Debye-Waller 扫描
│   ├── ex26_eley_rideal_surface_reaction.f90    # H+H/Cu(111) 气-固 ER 催化反应、超热放热能量分配与高振动激发态反转分布
│   ├── ex27_surface_electronic_friction_gle.f90 # NO/Au(111) 表面非绝热散射、广义朗之万 GLE 碰撞轨迹与电子-空穴对能损
│   ├── ex28_grazing_fast_atom_diffraction.f90   # keV He 掠入射快原子表面量子衍射 (GIFAD)、经典彩虹角与亚皮米波纹度反演
│   ├── ex29_cold_ion_atom_scattering.f90        # Yb+/Li 与 Ba+/Rb 冷离子-原子杂化碰撞截面与 Paul 阱微运动致热率
│   ├── ex30_tully_surface_hopping.f90           # Tully SAC 与 DAC 双避免交叉斯托克斯干涉系综非绝热动力学
│   ├── ex31_molecular_alignment_revival.f90     # N2 与 CO2 飞秒无场转动复苏与光学离心机超转子加速
│   ├── ex32_optical_lattice_bose_hubbard.f90    # 87Rb 光晶格深度扫描、Bose-Hubbard U/J 相变与引力 Bloch 振荡
│   ├── ex33_rph_variational_transition_state.f90 # 多原子反应路径变分 CVT 速率与 Eckart 量子隧穿增强因子
│   ├── ex34_relativistic_dirac_cesium.f90       # 铯原子 6s/6p/5d 相对论狄拉克能级、精细结构分裂与 D1/D2 振子强度
│   ├── ex35_rixs_core_level_spectroscopy.f90    # 铜氧化物 Cu L3 共振非弹性 X 射线散射 2D 能损图谱与声子伴峰
│   └── build_examples.sh           # 算例编译运行脚本 (全 35 算例编译运行通过)
└── python/                        # Python 辅助分析与可视化套件 (pygenmod)
    ├── pyproject.toml
    ├── test_pygenmod.py           # Python 单元测试 (100% Pass, 12/12 测试)
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
        ├── field_scattering.py    # 四大基组幺正变换、Breit-Rabi图谱与磁Feshbach色散拟合
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
- `calc_feshbach_resonance_scan(r_grid, v_mat, mass, energy_grid, n_energies, thresholds, l_channels, s_wave_length, eigenphase_sums, stat)`: 跨 Feshbach 共振能量扫描，追踪闭通道准束缚态引起的开通道散射长度极点发散 $a_s(E) \to \pm \infty$ 与相移特征跳变。
- `segmented_grid_t`: 多扇区自适应变步长分段径向网格派生类型，支持短程深阱区密格点、长程弱渐近区稀疏格点的高效离散。
- `create_segmented_grid(r_start, r_bounds, dr_steps, grid, stat)`: 多扇区分段径向网格构造器，自动分段对齐与全局单调展平。
- `calc_scattering_length_segmented_numerov(grid, v_pot, mass, a_s, u_zero, stat)`: 分段网格零能 Numerov 散射长度求解器，各扇区交接面采用 4 阶 Taylor 导数光滑桥接。
- `calc_scattering_wavefunction_segmented_ti(grid, v_pot, mass, energy, l, norm_type, u_wf, phase_shift, stat)`: 分段网格连续散射态能量本征波函数求解器。
- `calc_phase_shift_segmented(grid, v_pot, mass, energy, l, delta, stat)`: 分段网格分波相移与反应矩阵求解。
- `calc_multichannel_close_coupling_segmented_logder(grid, v_mat, mass, total_energy, thresholds, l_channels, res, stat)`: **分段网格多通道定态密耦 Johnson 矩阵对数导数求解器**，扇区间局域对数导数矩阵精确传递，大幅节约深势阱多通道计算量。

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

### 21. 外场电磁场超冷散射与四大基组变换 (`mod_field_scattering`)
- 派生类型：`type(cold_atom_t)`, `type(field_channel_t)`, `type(field_feshbach_result_t)`。
- 基组常量：`BASIS_UNCOUPLED`（非耦合基 $|s_1 m_{s1} i_1 m_{i1} s_2 m_{s2} i_2 m_{i2} l m_l\rangle$）、`BASIS_F_COUPLED`（单体自旋耦合基 $|f_1 m_{f1} f_2 m_{f2} l m_l\rangle$）、`BASIS_TOTAL_SPIN`（两体总自旋基 $|(s_1 s_2)S (i_1 i_2)I F M_F l m_l\rangle$）、`BASIS_FIELD_DRESSED`（外场渐近本征态缀饰基）。
- `get_cold_atom_preset(name, atom, stat)`: 获取预置碱金属同位素原子参数（支持 $^6\text{Li}, ^7\text{Li}, ^{23}\text{Na}, ^{39}\text{K}, ^{40}\text{K}, ^{87}\text{Rb}, ^{133}\text{Cs}$ 等）。
- `calc_breit_rabi_energies(atom, b_field_gauss, energies, stat)`: 任意磁场下单原子 Zeeman-超精细 Breit-Rabi 能级全数值求解。
- `build_field_collision_channels(atom1, atom2, basis_type, two_Mtot, l_max, channels, stat)`: 根据总磁量子数守恒与分波截断构建碰撞通道空间。
- `calc_basis_transform_matrix(atom1, atom2, channels_src, channels_dst, basis_src, basis_dst, u_mat, stat)`: 精确构造四大物理基组之间的幺正变换矩阵 $\mathbf{U}$。
- `calc_zeeman_hyperfine_hamiltonian(atom1, atom2, b_field_gauss, channels, h_zeeman, stat)`: 组装外磁场下单体与两体塞曼-超精细哈密顿矩阵。
- `calc_magnetic_feshbach_resonance_scan(atom1, atom2, b_grid, n_b, e_col, v_singlet, v_triplet, res, stat)`: 扫描磁场并拟合磁 Feshbach 共振中心 $B_0$、宽度 $\Delta B$、背景散射长度 $a_{\text{bg}}$ 与零交叉点 $B_{\text{zero}}$。

### 22. 各向异性偶极超冷散射与自旋弛豫 (`mod_dipolar_scattering`)
- 派生类型：`type(polar_molecule_t)`, `type(dipolar_channel_t)`。
- `c2q_spherical_harmonic_tensor(theta, phi, q)`: 秩-2 球谐空间张量算符 $C_{2, q}(\theta, \phi) = \sqrt{\frac{4\pi}{5}} Y_{2, q}(\theta, \phi)$。
- `c2q_orbital_matrix_element(l1, m1, l2, m2, q)`: 轨道角动量分波球谐张量矩阵元 $\langle l_1, m_1 | C_{2, q} | l_2, m_2 \rangle$（含 Wigner 3j 符号与奇偶性校验）。
- `spin_tensor_coupled_matrix_element(s1, s2, s_tot, ms_tot, s_prime, ms_prime, q)`: 秩-2 两体自旋张量算符 $[\mathbf{s}_1 \otimes \mathbf{s}_2]^{(2)}_q$ 耦合矩阵元（严格自旋阶梯递推与相因子）。
- `calc_mddi_coupling_strength(mu1_bohr, mu2_bohr)`: 磁偶极-偶极相互作用 (MDDI) 强度 $C_{\text{dd}} = \frac{\mu_0 \mu_1 \mu_2}{4\pi}$ (a.u.)。
- `calc_dipolar_relaxation_cross_section(energy, mass, mu_mag, delta_m, sigma_rel)`: 磁阱超冷碰撞两体非弹性自旋弛豫截面 $\sigma_{\text{rel}}(E)$。
- `calc_dipolar_relaxation_thermal_rate(temp_kelvin, mass, mu_mag, delta_m, k_rel)`: 玻尔兹曼系综平均超冷磁偶极自旋弛豫速率系数 $K_{\text{rel}}(T)$。
- `calc_stark_induced_dipole(mol, e_field_dc, d_ind)`: 刚体极性分子外加直流电场 Stark 诱导电偶极矩 $d_{\text{ind}}(\mathcal{E})$。
- `calc_electric_dipolar_length(mass, dipole_debye, a_d)`: 电偶极相互作用特征偶极长度 $a_d = \frac{m d^2}{2 \hbar^2}$。
- `build_dipolar_channel_basis(s_tot, l_max, channels, stat)`: 构建包含轨道多分波 $(l, m_l)$ 与两体自旋 $(S, M_S)$ 的全通道基矢。
- `calc_dipolar_potential_matrix(r, channels, mu1_bohr, mu2_bohr, v_mat, stat)`: 构造各向异性偶极多通道耦合势矩阵，驱动 s-波与 d-波间的本征角动量转移。

### 23. 超冷光缔合谱学与分子生成 (`mod_photoassociation`)
- 派生类型：`type(pa_transition_t)`, `type(pa_rate_result_t)`。
- `calc_free_bound_fc_overlap(r_grid, psi_free, psi_bound, overlap)`: 能量归一化自由散射态 $\psi_E(r)$ 与激发振动态 $\psi_v(r)$ 自由-束缚态 Franck-Condon 空间重叠积分 $I_{\text{FB}} = \int_0^\infty \psi_E(r) \psi_v(r) dr$。
- `calc_free_bound_fc_density(energy, overlap, f_fb)`: 自由-束缚跃迁态密度 $f_{\text{FB}}(E) = |I_{\text{FB}}(E)|^2$。
- `calc_pa_stimulated_linewidth(laser_intensity_w_cm2, d_trans_debye, overlap, gamma_stim_au)`: 激光辐射诱导受激光缔合跃迁线宽 $\hbar \Gamma_{\text{stim}} = 2\pi \left(\frac{I}{2\varepsilon_0 c}\right) |d_{\text{el}} I_{\text{FB}}|^2$。
- `calc_pa_cross_section(trans, energy, detuning, cross_sec_au)`: 单能量入射超冷原子光缔合吸收散射截面 $\sigma_{\text{PA}}(E, \Delta)$（幺正共振洛伦兹线型）。
- `calc_pa_thermal_rate_coefficient(trans, temp_kelvin, detuning, k_pa)`: 麦克斯韦-玻尔兹曼热平衡系综平均光缔合速率系数 $K_{\text{PA}}(T, \Delta)$。
- `calc_pa_detuning_scan(trans, temp_kelvin, detunings, n_pts, rates, peak_detuning, stat)`: 激光失谐频率连续扫描谱线生成与半高全宽（FWHM）/共振峰位提取。
- `calc_twophoton_raman_coupling(omega1, omega2, delta1, omega_eff)`: 双光子 Raman / STIRAP 绝热受激跃迁生成振转基态分子的有效拉比频率 $\Omega_{\text{eff}} = \frac{\Omega_1 \Omega_2}{2\Delta_1}$。

### 24. 超冷三体复合与 Efimov 少体物理 (`mod_three_body_recombination`)
- `solve_efimov_s0_identical_bosons()`: 精确超越方程求解全同玻色子 Efimov 标度指数 $s_0 \approx 1.00624$ 与标度常数 $e^{\pi/s_0} \approx 22.7$。
- `calc_three_body_recombination_a_positive(param, a_scat)`: $a>0$ 侧普适 $a^4$ 复合损失速率 $K_3$ 与干涉极小值。
- `calc_three_body_recombination_a_negative(param, a_scat)`: $a<0$ 侧三体复合共振峰 $a_-^{(n)}$。
- `calc_unitary_three_body_loss_temperature(temp_kelvin, mass)`: 酉极限幺正饱和温度依赖复合速率 $K_3 \propto T^{-2}$。

### 25. 低维光晶格受限散射与约束诱导共振 (`mod_confined_scattering`)
- `init_waveguide_1d(omega_perp_hz, mass_amu, wg)`: 初始化准一维横向简谐光阱波导与振荡长度 $a_\perp$。
- `calc_olshanii_cir_parameters(wg, a_3d_bohr, g_1d, a_1d)`: Olshanii 约束诱导共振极点 $a_{\text{CIR}} \approx 1.0326 a_\perp$ 与有效 1D 耦合强度 $g_{\text{1D}}$。
- `calc_confined_dimer_binding_energy(wg, a_3d_bohr)`: 准一维横向陷阱修饰分子结合能。
- `calc_lieb_liniger_parameter(n_1d, g_1d, mass, gamma_ll)`: Lieb-Liniger 关联参数 $\gamma_{\text{LL}}$，判别 Tonks-Girardeau 强关联费米化气体。

### 26. 自电离体系、Fano 共振与复坐标旋转 (`mod_autoionization_fano`)
- `calc_fano_profile(energy, e0, gamma, q_param)`: 计算 Fano 不对称吸收线型 $\sigma(\epsilon) = \sigma_0 \frac{(q+\epsilon)^2}{1+\epsilon^2}$ 与反共振零点。
- `calc_autoionization_lifetime(gamma_au)`: 自电离共振寿命 $\tau = \hbar / \Gamma$（飞秒 fs）。
- `solve_ccr_resonance_model(e_discrete, v_coupl, theta_rot)`: 复坐标旋转法（CCR）非厄米哈密顿量本征求解，提取复本征能量 $E - i\Gamma/2$。

### 27. 交叉静电磁场转振-自旋动力学 (`mod_crossed_field_scattering`)
- `init_crossed_field_config(e_kv_cm, b_gauss, tilt_angle_deg, cfg)`: 配置任意倾角 $\beta$ 的交叉静电磁场 $\mathbf{E} \times \mathbf{B}$。
- `solve_crossed_field_eigenstates(cfg, j_max, s_spin, energies, states)`: 严格对角化包含宇称破缺 Stark 效应与 Zeeman 效应的分子转动自旋哈密顿量。
- `calc_crossed_field_observables(state, dipole_debye, cos_theta, sz)`: 计算实验室系诱导电取向度 $\langle\cos\theta\rangle$ 与自旋投影。

### 28. 三原子反应散射、Jacobi 坐标与几何相位 (`mod_triatomic_geometry`)
- `jacobi_to_internuclear(jacobi, masses, dist)` / `internuclear_to_jacobi`: 任意三原子质量体系 Jacobi 坐标与核间距解析双向变换。
- `calc_leps_potential(dist, leps)`: London-Eyring-Polanyi-Sato (LEPS) 全势能面评估与对称鞍点反应势垒计算。
- `calc_conical_intersection_adiabats(x, y, ci, e_lower, e_upper)`: 线性锥形交叉（CI）双绝热势能面与能隙分裂。
- `calc_berry_phase_around_ci(ci, radius, n_steps)`: 闭合回路数值积分提取精确拓扑几何相位 $\Phi_B = \pi$。

### 29. 旋量玻色爱因斯坦凝聚自旋动力学 (`mod_spinor_bec`)
- `init_spinor_preset(preset_name, b_field_gauss, density_cm3, param)`: 装载 $^{87}\text{Rb}$（铁磁相 $c_2<0$）与 $^{23}\text{Na}$（反铁磁相 $c_2>0$）实验参数。
- `calc_quadratic_zeeman_shift(b_gauss, atom_name)`: Breit-Rabi 二阶塞曼能量位移 $q_Z(B) \propto B^2$。
- `propagate_spinor_sma_rk4(param, dt_au, state)`: 单模近似（SMA）下高精度 RK4 保全几率与保磁化强度相干自旋演化。

### 30. 三原子超球面反应动力学与热速率常数 (`mod_hyperspherical_reactive`)
- `init_reaction_mass(ma, mb, mc, rmass)`: Delves 质量标度超球面坐标变换因子 $d$ 与反应偏转角 $\beta_{skew}$（$\text{H}+\text{H}_2 \to 60^\circ$）。
- `calc_eckart_transmission(energy, v_b, omega_im)`: 鞍点 Eckart 势垒精确量子隧穿传递几率 $P(E)$。
- `calc_cumulative_reaction_probability(ts, energy)`: 跨势垒多振动态通道求和累积反应几率 $N(E)$。
- `calc_canonical_rate_constant(ts, rmass, temp_k, n_steps)`: 正则全量子热反应速率常数 $k(T)$ 严格玻尔兹曼积分。
- `calc_tst_wigner_rate(v_b, omega_im, temp_k, prefactor)`: 经典过渡态理论（TST）与 Wigner 势垒隧穿修正速率。

### 31. 超冷偶极量子液滴与李-黄-杨量子涨落 (`mod_dipolar_droplets_lhy`)
- `init_dipolar_droplet_param(atom_name, a_scat_bohr, param)`: 磁偶极特征长度 $a_{dd}$ 与相对偶极强度 $\epsilon_{dd} = a_{dd} / a_s$。
- `calc_pelster_lima_q5(eps_dd)`: Pelster-Lima 超越平均场 LHY 零点量子涨落修正积分 $Q_5(\epsilon_{dd})$。
- `calc_equilibrium_droplet_density(param)`: 自由空间零压平衡平顶自束缚量子液滴核心密度 $n_0$。
- `calc_droplet_chemical_potential(param, density)`: 自束缚负化学势 $\mu(n_0) < 0$ 稳定性判据。
- `calc_critical_atom_number(param)`: 3D 自束缚量子液滴相变与蒸发临界原子数 $N_{crit}$。

### 32. 强场非顺序双电离与电子重碰撞动量谱 (`mod_strong_field_nsdi`)
- `init_nsdi_laser(wavelength_nm, intensity_w_cm2, laser)`: 激光峰值电场 $F_0$、角频率 $\omega$ 与有质动力势 $U_p = F_0^2 / (4\omega^2)$。
- `calc_recollision_trajectory(phi_0, up_au, phi_r, e_rec, ok)`: 求解经典电子轨道回碰根 $x(\phi_r) = 0$ 与 $3.173 U_p$ 动能截断。
- `calc_lotz_cross_section(e_rec, ip2)`: 电子碰撞电离 $(e, 2e)$ Lotz 经验截面。
- `calc_nsdi_2d_momentum_dist(laser, target, n_pts, p_max, grid, dist_2d, corr)`: 生成全周期纵向平行动量关联谱 $P(p_{z1}, p_{z2})$ 与 COLTRIMS 正向关联特征系数 $C_{corr} > 0$。
- `calc_double_ion_yield_curve(wavelength, target, n_int, imin, imax, ints, y_nsdi, y_sdi)`: 扫描光强并重现双电离非顺序“膝盖平台结构”（Knee Structure）。

### 33. 磁与光 Feshbach 共振与分子弱束缚态 (`mod_feshbach_bound_states`)
- `init_mfr_preset(name, mfr)`: 装载 $^{6}\text{Li}$（宽共振 $s_{res} \approx 51$）与 $^{87}\text{Rb}$（窄共振 $s_{res} \approx 0.13$）物理参数。
- `calc_mfr_scattering_length(mfr, b_gauss)`: 外磁场依赖 s 波有效散射长度 $a(B) = a_{bg}(1 - \Delta B / (B - B_0))$。
- `calc_mfr_bound_energy_coupled(mfr, b_gauss)`: 耦合通道有限相互作用程 $R^*$ 精确分子结合能 $E_b(B)$。
- `calc_mfr_closed_channel_fraction(mfr, b_gauss)`: Hellmann-Feynman 定理分子态闭通道成分占比 $Z(B) = 1 - 1/\sqrt{1 + 2R^*/a(B)}$。
- `calc_ofr_complex_scattering_length(ofr, delta_hz, a_re, a_im)`: 光 Feshbach 共振（OFR）色散复散射长度 $\tilde{a}(\Delta_L)$。
- `calc_ofr_inelastic_loss_rate(ofr, delta_hz)`: 实验可测光致双体非弹性损失速率常数 $K_2(\Delta_L)$（$\text{cm}^3/\text{s}$）。

### 34. 阿秒瞬态吸收光谱与光诱导态自电离干涉 (`mod_attosecond_transient_absorption`)
- `init_atas_helium_benchmark(bright_state)`: 氦原子 $2s2p (^1P)$ 经典双激发态基准 ($E_0 = 60.15\text{ eV}, \Gamma = 0.037\text{ eV}, q_0 = -2.80$)。
- `calc_laser_dressed_fano_q(q0, phase_shift)`: 相位微扰模型 (PPM) 下激光修饰含时动态 Fano 不对称参数 $q(\tau)$。
- `calc_light_induced_state_energy(e_bright, e_dark, omega_nir, rabi)`: 强激光缀饰诱导的光诱导态 (LIS) 特征能级位置。
- `calc_quantum_beat_period_fs(delta_e_ev)`: 明暗态相干干涉拍频特征周期 $T_{\text{beat}} = h / \Delta E$。
- `calc_atas_spectrum(state, nir_i, nir_lam, ne, emin, emax, ntau, taumin, taumax, e_grid, tau_grid, spec_2d)`: 全量计算二维能量-时延瞬态吸收差分光密度矩阵 $\Delta\text{OD}(\omega, \tau)$。

### 35. 双色反向旋转圆偏振场与分子光电子圆二色性 (`mod_bicircular_pecd`)
- `init_bicircular_field(field, omega1, r_freq, i1, i2, h1, h2, phi1, phi2, fwhm, env)`: 双色椭圆/圆偏振场合成与时频参数配置。
- `calc_dynamical_symmetry_fold(h1, h2, freq_ratio)`: 计算离散动力学旋转对称度（反向旋转 $\omega+2\omega$ 输出 $C_3$ 三叶草对称）。
- `init_chiral_tetrahedral_molecule(mol, enantiomer)`: 构建四中心手性对映体分子（$R$-型与 $S$-型）。
- `calc_chirality_measure(mol)`: 计算并验证伪标量手性不变量 $\chi(R) = -\chi(S)$。
- `calc_chiral_beta1_model(mol, e_ev, photon_ev)`: 手性四面体势奇宇称不对称参数 $\beta_1(E)$ 模型。
- `calc_forward_backward_asymmetry(beta1)`: Ritchie 光电子前后发射不对称度百分比 $G_{\text{PECD}} = \beta_1 / 2$。
- `calc_pecd_pad_spectrum(b1, b2, gamma33, n_th, n_ph, th_grid, ph_grid, pad)`: 合成三维手性螺旋浆光电子角分布 (PAD)。

### 36. 超冷极性分子反应动力学与微波/静电偶极遮蔽 (`mod_ultracold_reaction_shielding`)
- `init_ultracold_molecule_preset(mol, name)`: 初始化 $^{40}\text{K}^{87}\text{Rb}$, $^{23}\text{Na}^{87}\text{Rb}$, $^{23}\text{Na}^{40}\text{K}$ 实验质量与电偶极矩。
- `calc_effective_shielding_potential(mol, cfg, r, l, v_eff)`: 计算微波蓝失谐免交叉排斥偶极屏蔽相互作用势 $V_{\text{eff}}(R)$。
- `calc_shielding_barrier_height(mol, cfg, r_bar, v_bar)`: 计算长程工程化排斥势垒位置 $R_{\text{bar}} \sim 350\ a_0$ 与势垒高度 $V_{\text{bar}}$。
- `calc_wkb_tunneling_probability(mol, cfg, e_coll, t_tunnel)`: WKB 量子隧穿衰减因数，验证短程反应区量子反射。
- `calc_shielded_scattering_rates(mol, cfg, temp, k2_el, k2_inel, gamma)`: 提取弹性与反应损耗速率，满足蒸发冷却判据 $\gamma = K_2^{(\text{el})} / K_2^{(\text{inel})} > 100$。

### 37. 里德堡原子阻塞、PXP 约束模型与量子多体疤痕 (`mod_rydberg_blockade`)
- `init_rydberg_atom(atom, species, n, l)`: $^{87}\text{Rb}$ 高主量子数 $nS$ 态量子缺陷校正与 $C_6 \propto n^{11}$ 标度相互作用。
- `calc_rydberg_blockade_radius(atom, rabi)`: 解析计算里德堡阻塞半径 $R_b = (|C_6| / \hbar\Omega)^{1/6}$。
- `calc_two_atom_dynamics(atom, spacing, rabi, delta, tmax, n_steps, t_arr, pg, ps, pd)`: 细化自适应步长 RK4 积分，验证强阻塞区双激发 $|rr\rangle$ 深度抑制与 $\sqrt{2}\Omega$ 集体纠缠态振荡。
- `calc_z2_order_parameter(n_atoms, occ)`: 计算一维反铁磁 Néel 空间交错电荷密度波序参量 $\mathcal{O}_{Z_2}$。
- `calc_rydberg_scar_dynamics(cfg, atom, tmax, n_steps, t_arr, z2_arr)`: 模拟 PXP 拓扑约束链量子多体疤痕长寿命宏观相干复苏。

### 38. 表面量子散射与选择性吸附共振 (`mod_surface_scattering`)
- `init_surface_lattice(lat, ax, ay, zetax, zetay, m_sub, theta_d)`: 2D 晶格常数、倒格矢 $\mathbf{G}$ 与表面波纹度配置。
- `init_surface_potential_morse(pot, well_depth, alpha, mass, n_bound)`: 表面吸引阱深与 Morse 束缚态能级求解。
- `calc_surface_diffraction_channels(lat, mass, energy, theta, phi, max_m, n_ch, channels)`: 2D 衍射通道开/闭判据与出射角。
- `calc_hcs_diffraction_probabilities(lat, mass, energy, theta, phi, max_m, n_ch, channels)`: 硬波纹表面程函近似贝塞尔衍射强度与幺正归一化。
- `calc_selective_adsorption_resonance(lat, pot, mass, energy, theta, phi, m, n, v, is_near, delta_e, fano_ratio)`: 闭通道束缚态共振与 Fano 线型调制。
- `calc_surface_debye_waller(lat, mass, kiz, kzg, temp)`: 表面声子非弹性热激发 Debye-Waller 衰减因子。

### 39. 气-固表面催化与 Eley-Rideal 反应动力学 (`mod_surface_reaction_er`)
- `init_er_reaction_system(sys, name)`: 预置直接 Eley-Rideal 反应体系参数与巨大放热量 $\Delta E_{\text{exo}}$。
- `calc_er_potential_2d(sys, r, z_cm, v_pot)`: 反应路径 2D 势能面 $V(r, Z_{\text{cm}})$。
- `calc_er_energy_partitioning(sys, e_inc, partition)`: 超热放热能量通道分配（振动 $\sim 50\%$、平动 $\sim 35\%$、基底 $\sim 10\%$）。
- `calc_er_vibrational_populations(sys, e_inc, max_v, v_dist)`: 产物分态振动布居反转分布 $P(v)$。
- `calc_er_reaction_cross_section(sys, e_inc)`: 入射动能依赖反应截面 $\sigma_{\text{ER}}(E_i)$。
- `calc_er_thermal_rate_constant(sys, temp)`: 温度依赖准无势垒热催化速率常数 $k_{\text{ER}}(T)$。

### 40. 金属表面非绝热动力学与电子摩擦耗散 (`mod_surface_electronic_friction`)
- `init_metal_surface(surf, name, temp)`: 金属基底 (Au, Cu, Pt) 费米能与峰值摩擦系数 $\eta_0$。
- `calc_electronic_friction_coeff(surf, z)`: 局域密度摩擦近似 (LDFA) 空间指数衰减分布 $\eta(z)$。
- `calc_surface_morse_force(z, well, alpha, ze, v_pot, force)`: 绝热表面保守 Morse 力与排斥壁。
- `integrate_gle_scattering_trajectory(surf, mass, e_inc, dt, n_steps, t_arr, z_arr, v_arr, loss_res)`: 广义朗之万方程 (GLE) 动力学积分与电子-空穴对能损 $\Delta E_{\text{loss}}$。
- `calc_vibrational_relaxation_rate(surf, z, mass)`: 表面吸附分子高频化学键振动弛豫速率与寿命 $\tau_{\text{vib}}$。

### 41. 掠入射快原子表面量子衍射与彩虹散射 (`mod_grazing_fast_atom_diffraction`)
- `init_gifad_experiment(cfg, projectile, mass, e_kev, theta_deg, ax, zeta)`: keV 准直束轴向沟道快慢自由度解耦与有效横向动能 $E_\perp \sim \text{eV}$。
- `calc_gifad_transverse_kinematics(cfg, e_perp, lambda_perp)`: 横向垂直能量与德布罗意量子波长 $\lambda_\perp$。
- `calc_gifad_rainbow_angle(cfg)`: 表面经典彩虹散射角 $\theta_R$。
- `calc_gifad_diffraction_spectrum(cfg, max_order, spec)`: 1D 横向量子 Bragg 衍射谱与彩虹包络调制。
- `calc_surface_corrugation_from_rainbow(ax, theta_r)`: 亚皮米级表面波纹幅度 $\zeta$ 逆向反演重构。

### 42. 冷离子-中性原子杂化散射与极化阱动力学 (`mod_ion_atom_scattering`)
- `init_ion_atom_system(sys, m_ion, m_atom, q_ion, alpha_atom, stat)`: 初始化杂化冷碰撞系统特征尺度 $R^*$ 与能量尺度 $E^*$。
- `calc_langevin_cross_section(sys, e_coll_ev)` / `calc_langevin_rate_coefficient(sys)`: 经典 Langevin 螺旋俘获临界碰撞参数 $b_c$、截面 $\sigma_L(E)$ 与能量无关速率 $K_L$。
- `calc_ion_atom_phase_shift(sys, l_orb, e_coll_ev, r_match)`: $1/r^4$ 极化势分波相移与微分散射求解。
- `calc_mere_phase_shift_s_wave(sys, a_s, r_eff, e_coll_ev)`: 修正有效力程展开 (MERE) 极化奇异散射修正。
- `calc_rf_micromotion_heating(sys, omega_rf, q_param, temp_ion, temp_atom, heat_rate, t_limit)`: Paul 射频阱中微运动导致的离-原碰撞致热率与平衡极限温度。

### 43. 最少开关表面跳跃与非绝热混合量子-经典动力学 (`mod_surface_hopping_fssh`)
- `init_tully_model(model, model_type, stat)`: Tully 经典非绝热三大基准模型 (SAC 单避免交叉, DAC 双避免交叉与斯托克斯干涉, ECR 扩展耦合反射)。
- `calc_adiabatic_surface_and_nacv(model, r, v_adia, d_nacv)`: 绝热势能面与非绝热导数耦合矢量 (NACV) 解析投影。
- `init_fssh_trajectory(traj, r_init, p_init, init_state)`: 初始化 Velocity Verlet 核推进与电子态波函数相干矢量。
- `propagate_fssh_step(model, traj, dt, jumped)`: Tully 最少开关几率判断、能量守恒核动量沿 NACV 重标度与禁阻跳跃反射修正。
- `run_fssh_ensemble(...)` / `propagate_ehrenfest_step(...)`: 蒙特卡洛系综分支比统计与 Ehrenfest 平均场动力学比对。

### 44. 强场分子定向、取向与超转子动力学 (`mod_molecular_alignment`)
- `init_rotor_molecule(mol, name, b_rot, delta_alpha, dipole, d_e, stat)`: 刚体线性分子转动结构、各向异性极化率与电偶极矩。
- `calc_cos2_matrix_elements(...)` / `calc_cos_matrix_elements(...)`: 球谐角动量基底偶极矩阵元与外场极化耦合矩阵。
- `simulate_laser_induced_alignment(mol, dt, t_max, i_laser, duration, t_grid, cos2_trace, stat)`: 飞秒非绝热激光脉冲驱动转动波包形成、无场复苏序参量 $\langle\cos^2\theta\rangle(t)$。
- `calc_optical_centrifuge_kick(mol, alpha_chirp, duration, j_super)`: 光学离心机恒定角加速度将分子加速至极端高角动量超转子态 ($J \gg 1$)。
- `calc_superrotor_dissociation(mol, j_rot, is_dissociated, e_rot)`: 极端离心旋转势垒高度与分子共价键机械破键判定。

### 45. 超冷光晶格与玻色-哈伯德微观映射 (`mod_optical_lattice_hubbard`)
- `init_optical_lattice(latt, mass, lambda_nm, s_depth, stat)`: 光晶格激光波矢、反冲能量 $E_R$ 与周期驻波深度 $V_0/E_R$。
- `calc_bloch_band_energies(latt, q_quasi, n_bands, energies)`: Mathieu 方程动量空间展开、Bloch 能带结构与带隙。
- `calc_bose_hubbard_parameters(latt, a_s_nm, bh, stat)`: Wannier 局域轨道重叠积分、严格紧束缚参数（跃迁能 $J$ 与在位能 $U$）以及超流-Mott 相变判据 $(U/J)_c$。
- `calc_bloch_oscillation_dynamics(latt, force, t_bloch, nu_bloch, p_lz)`: 恒定外力下布洛赫振荡周期 $T_B$ 与第一激发带 Landau-Zener 带间隧穿几率 $P_{\text{LZ}}$。

### 46. 多原子反应路径哈密顿量与变分过渡态理论 (`mod_reaction_path_hamiltonian`)
- `init_rph_benchmark_reaction(path, reaction_type, stat)`: 内禀反应坐标 (IRC) 路径曲率、能量梯度与垂直正交振动模频率演化。
- `calc_generalized_tst_rate(path, s_idx, temp_k, rate_gtst)`: 沿反应坐标任意分界面的广义过渡态理论 (GTST) 正则速率。
- `calc_cvt_rate_constant(path, temp_k, s_opt, rate_cvt, stat)`: 正则变分过渡态理论 (CVT) 自由能最大瓶颈寻优 $k^{\text{CVT}}(T) = \min_s k^{\text{GTST}}(T, s)$。
- `calc_eckart_tunneling_factor(v_barrier, imag_freq, temp_k, kappa_tunnel)`: 解析不对称 Eckart 势垒半经典量子穿透系数 $\kappa(T)$。

### 47. 相对论原子结构与径向狄拉克方程 (`mod_relativistic_atomic`)
- `calc_dirac_model_potential(z_nuclear, z_ion, alpha_core, r_cut, r)`: Norcross-Klapisch 相对论核心屏蔽库仑势与极化模型势。
- `solve_radial_dirac_eigenvalue(z_nuclear, z_ion, alpha_core, r_cut, n_princ, kappa, state, stat)`: 求解径向狄拉克方程，获得 Sommerfeld 相对论本征能量与有效量子亏损 $\mu$。
- `calc_dirac_fine_structure_splitting(state_lower, state_upper, delta_fs_ev, delta_fs_cm1)`: 纯相对论天然自旋-轨道耦合劈裂能（如类氢 $2p_{3/2} - 2p_{1/2}$ 与重碱金属双线）。
- `calc_dirac_e1_matrix_element(state_i, state_f, r_overlap, osc_strength)`: 大小分量径向重叠积分与相对论电偶极 (E1) 跃迁吸收振子强度 $f_{if}$。

### 48. 共振非弹性 X 射线散射与内壳层光谱 (`mod_resonant_xray_scattering`)
- `init_rixs_system(sys, e_init, e_inter, gamma_core, d_in, e_fin, gamma_fin, d_out, stat)`: 配置初态、核心激发中间态与低能终态流形。
- `calc_xas_cross_section(sys, omega_in_ev)`: 基于光学定理计算一阶 X 射线吸收截面 (XAS)。
- `calc_kramers_heisenberg_cross_section(sys, omega_in_ev, omega_loss_ev)`: Kramers-Heisenberg 二阶微扰公式相干干涉与能损截面。
- `calc_rixs_2d_map(sys, n_in, w_in_grid, n_loss, w_loss_grid, rixs_map)`: 生成高精度 2D RIXS 入射能量-能量损失散射强度图谱。
- `calc_huang_rhys_vibrational_rixs(omega_0, s_factor, gamma_core, detuning, n_max, loss_intensity)`: 基于相联拉盖尔多项式严格解析求解电-声耦合 Huang-Rhys 振动 Franck-Condon 伴峰级数。

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

算法库内建完备的单元测试，覆盖物理常数往返转换、半整数角动量耦合、FGH 谐振子能级、激光脉冲包络、Bloch / Split-Operator 模长守恒、强场原子模型与多态非绝热耦合、转振态激光调控、三次样条插值与外推、自相关吸收谱与碎片 KER 分支比、Lindblad 耗散主方程与 Krotov 最优控制、非含时散射长度与 S-矩阵、含时波包散射连续态透射谱、外加电磁场四大基组幺正变换与磁 Feshbach 共振拟合、各向异性偶极张量与自旋弛豫、超冷光缔合受激线宽与分子生成：

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
 1. Constants Tests:                    10 / 10 PASSED
 2. Special Function Tests:             12 / 12 PASSED
 3. DVR Grid Tests:                      8 /  8 PASSED
 4. Laser Pulse Tests:                   9 /  9 PASSED
 5. Propagator Tests:                    8 /  8 PASSED
 6. Extended Atomic Tests:              14 / 14 PASSED
 7. Rovibrational Control Tests:         8 /  8 PASSED
 8. Interpolation Tests:                 7 /  7 PASSED
 9. Photofragment & Flux Tests:          6 /  6 PASSED
10. Open Quantum & OCT Tests:           10 / 10 PASSED
11. TI Scattering Tests:                19 / 19 PASSED
12. TD Scattering Tests:                 9 /  9 PASSED
13. Field Scattering Tests:             42 / 42 PASSED
14. Dipolar Scattering Tests:           21 / 21 PASSED
15. Photoassociation Tests:             14 / 14 PASSED
16. Three-Body Recombination Tests:      5 /  5 PASSED
17. Confined Scattering CIR Tests:       5 /  5 PASSED
18. Autoionization Fano/CCR Tests:       5 /  5 PASSED
19. Crossed-Field Scattering Tests:      5 /  5 PASSED
20. Triatomic Geometry & Berry Tests:    5 /  5 PASSED
21. Spinor BEC Dynamics Tests:           5 /  5 PASSED
22. Hyperspherical Reactive Tests:       5 /  5 PASSED
23. Dipolar Droplets LHY Tests:          5 /  5 PASSED
24. Strong-Field NSDI Tests:             5 /  5 PASSED
25. Feshbach Bound States Tests:         5 /  5 PASSED
26. ATAS Transient Absorption Tests:     5 /  5 PASSED
27. Bicircular PECD Tests:               5 /  5 PASSED
28. Ultracold Reaction Shielding Tests:  5 /  5 PASSED
29. Rydberg Blockade Tests:              5 /  5 PASSED
30. Surface Scattering & SAR Tests:      5 /  5 PASSED
31. Surface Eley-Rideal Reaction Tests:  5 /  5 PASSED
32. Surface Electronic Friction Tests:   5 /  5 PASSED
33. GIFAD Fast Atom Diffraction Tests:   5 /  5 PASSED
34. Cold Ion-Atom Scattering Tests:      5 /  5 PASSED
35. Surface Hopping FSSH Tests:          5 /  5 PASSED
36. Molecular Alignment Tests:           5 /  5 PASSED
37. Optical Lattice & Hubbard Tests:     5 /  5 PASSED
38. RPH & Variational TST Tests:         5 /  5 PASSED
39. Relativistic Atomic & Dirac Tests:   5 /  5 PASSED
40. Resonant X-ray RIXS Tests:           5 /  5 PASSED
----------------------------------------------------------------
ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass, 335/335 断言通过)
================================================================
```

---

## 📊 典型物理算例 (Examples)

位于 `GeneralModule/examples/`，一键编译运行全部 35 大物理前沿算例：
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
7. **`ex07_scattering_wavefunctions_ti_td.f90`**
   - **连续谱散射能量本征波函数定量求解与交叉验证**：在同一一维势垒上对比非含时逆向 Numerov 匹配法与含时高斯波包 Split-Operator 谱投影法，直接给出高精度的连续态空间本征波函数 $\psi_E(x)$、隧穿透射几率 $T(E)$ 与反射几率 $R(E)$，两套独立算法相对偏差 $< 1.9\%$。
8. **`ex08_ultracold_feshbach_segmented.f90`**
   - **分段网格多通道超冷磁 Feshbach 共振色散扫描与拟合**：以 $^{87}\text{Rb} + ^{87}\text{Rb}$ 双通道碰撞系统为例，采用 3 扇区分段网格（短程深阱区 $dr=0.005\,a_0$、中间区 $dr=0.02\,a_0$、长程渐近区 $dr=0.1\,a_0$），在磁场范围 $B \in [70, 90]\,\text{Gauss}$ 内高精度扫描色散散射长度 $a_s(B)$，并通过非线性拟合准确提取磁共振中心位置 $B_0$、共振宽度 $\Delta B$、背景散射长度 $a_{\text{bg}}$ 与零散射点 $B_{\text{zero}}$。
9. **`ex09_dipolar_relaxation_scattering.f90`**
   - **超冷磁偶极自旋弛豫与极性分子外电场 Stark 诱导偶极扫描**：模拟弱磁阱中 $^{87}\text{Rb}$ 原子受各向异性磁偶极相互作用驱动的二阶超精细两体自旋弛豫过程，计算碰撞能量扫描与微开尔文热平均弛豫速率 $K_{\text{rel}}(T)$；同步模拟刚体超冷极性双原子分子 $^{40}\text{K}^{87}\text{Rb}$ 在直流外电场 $\mathcal{E}$ 下的本征转动态 Stark 混叠，追踪诱导电偶极矩 $d_{\text{ind}}(\mathcal{E})$ 趋于永久偶极矩饱和，并计算电偶极相互作用特征长度 $a_d$。
10. **`ex10_photoassociation_spectroscopy.f90`**
    - **超冷原子光缔合光谱学与双光子 Raman 缔合分子态**：模拟超冷 $^{87}\text{Rb}$ 碰撞对由自由连续态受激跃迁至激发束缚态分子（$0_u^+ / 1_g$ 振动态）的光缔合吸收过程，计算自由-束缚态 Franck-Condon 重叠积分 $I_{\text{FB}}$、激光强度依赖的受激线宽 $\Gamma_{\text{stim}}$、单能量吸收截面与不同温度下的光缔合速率系数 $K_{\text{PA}}(T, \Delta)$ 洛伦兹-不对称展宽能谱，并计算利用双光子 STIRAP 绝热转移至超冷振转基态分子所需的有效双光子拉比耦合强度 $\Omega_{\text{eff}}$。
11. **`ex11_three_body_efimov_recombination.f90`**
    - **超冷三体 Efimov 复合速率与干涉极小值/共振峰扫描**：精确解算 Efimov 超越代数方程超越根 $s_0 = 1.00624$，在 $a < 0$ 区域扫描 Efimov 三体束缚态引起的共振复合损耗峰 $a_-^{(n)}$ 与在 $a > 0$ 区域由不同 Stueckelberg 几何路径干涉导致的复合速率极小值（Efimov 窗口 $a_+^{(n)}$），并给出强幺正极限下有限温度普适 $T^{-2}$ 速率标度。
12. **`ex12_confined_cir_scattering.f90`**
    - **准一维光晶格波导中 87Rb 约束诱导共振 CIR 与分子结合能**：解算横向简谐受限势阱下 Olshanii 约束诱导共振（CIR）发散极点，模拟外加磁场调控下 3D 自由空间散射长度越过临界阈值 $a_s = a_\perp / C$ 时有效 1D 相互作用强度 $g_{1D}$ 的极点发散与束缚分子二聚体结合能 $E_b$，展示低维强关联 Tonks-Girardeau 极限。
13. **`ex13_autoionization_fano_resonance.f90`**
    - **自电离 Fano 不对称吸收谱与复坐标旋转 CCR 寿命提取**：模拟连续态与离散态组态干涉产生的 Fano 非对称吸收轮廓（抗共振零点与不对称参数 $q$），结合复坐标旋转法（CCR）实现无反射人工边界，在复能量平面上准确对角化并分离非物理连续态与物理自电离极点共振寿命 $\Gamma$。
14. **`ex14_spinor_bec_dynamics.f90`**
    - **87Rb 与 23Na 凝聚体宏观自旋混合动力学与铁磁相图**：针对 $F=1$ 旋量玻色爱因斯坦凝聚体，计算单模近似（SMA）下相干自旋振荡动力学，精确保持总粒子数归一化与纵向磁化强度守恒，重现反铁磁（极性）与铁磁基态相图。
15. **`ex15_crossed_field_stark_zeeman.f90`**
    - **交叉电磁场中极性顺磁分子 Stark-Zeeman 态混合与空间定向**：模拟同时处在外加静电场与静磁场（夹角 $\beta$ 可调）中的双原子自由基分子，精确对角化非共线转动-超精细哈密顿量，揭示空间宇称破坏、交错免交叉能谱与空间取向度演化。
16. **`ex16_triatomic_reaction_berry_phase.f90`**
    - **三原子反应路径 LEPS 势能面与锥形交叉 Berry 几何相位**：基于质心 Jacobi 坐标严格映射三原子核构型，构建 London-Eyring-Polanyi-Sato (LEPS) 反应过渡态势垒，并在势能面避免交叉附近沿闭合回路数值积分非绝热几何相，精确提取 $\pi$ 模 Berry 几何相位。
17. **`ex17_hyperspherical_reaction_rates.f90`**
    - **三原子超球面反应动力学、Eckart 隧穿累积反应几率与正则热速率常数**：利用超球面坐标与质量标度坐标处理三原子重排反应，计算 Delves 反应偏角、Eckart 势垒传递几率、累积反应几率 $N(E)$ 以及宽温区正则反应热速率常数 $k(T)$。
18. **`ex18_dipolar_quantum_droplets.f90`**
    - **162Dy 偶极量子液滴自束缚平衡密度、负化学势与气-液滴相变**：求解含 Lee-Huang-Yang (LHY) 量子涨落修正与 Pelster-Lima $Q_5$ 各向异性极化积分的扩展 Gross-Pitaevskii 方程 (eGPE)，模拟超冷偶极原子自束缚量子液滴相变与平衡密度。
19. **`ex19_strong_field_nsdi_recollision.f90`**
    - **800nm 强场电子重碰撞动能 3.17 Up 截断、2D 动量关联谱与双电离膝盖结构**：基于 Corkum 三步半经典模型与电离相位抽样，模拟强激光场下隧穿电子回碰母离子过程，复现 $3.17\,U_p$ 经典重碰撞动能截止、(e,2e) 碰撞电离截面与双电离电离率“膝盖 (Knee)”台阶结构。
20. **`ex20_feshbach_molecular_bound_states.f90`**
    - **6Li 磁 Feshbach 晕轮二聚体结合能、闭通道权重与 87Rb 光 Feshbach 损耗**：求解双通道密耦哈密顿量，给出外加磁场下开-闭通道耦合对弱束缚分子结合能 $E_b(B)$ 与闭通道占据权重 $Z(B)$ 的调制，并计算光 Feshbach 共振受激吸收与非弹性光致双体损耗速率常数 $K_2$。
21. **`ex21_attosecond_transient_absorption.f90`**
    - **氦原子 2s2p 自电离、动态 Fano q、LIS 态与 2D ATAS 谱图**：结合阿秒极紫外 (XUV) 脉冲与强近红外 (NIR) 激光场，基于相位微扰模型 (PPM) 计算双激发态自电离吸收谱，演示动态 Fano 参数演化、光诱导态 (LIS) 杂化与 2D 延迟吸收图谱。
22. **`ex22_bicircular_pecd_chiral.f90`**
    - **双色反向圆偏振场、手性分子对映体能量分辨 PECD 谱与 2D PAD**：构建角动量对称性具备 $C_3$ 动力学对称的双色反向旋转椭圆偏振激光场，结合具有手性势能面的分子模型，计算光电子前后发射不对称度与能量分辨光电子圆二色性 (PECD)。
23. **`ex23_ultracold_molecule_shielding.f90`**
    - **KRb 超冷极性分子微波遮蔽势垒、隧穿抑制与蒸发冷却速率比 gamma**：针对微波场驱动的分子旋转态偶极相互作用工程，构建具有短程排斥势垒的长程修饰势，利用 WKB 积分评估短程反应非弹性损失抑制，计算弹性-非弹性碰撞截面比 $\gamma > 100$ 的蒸发冷却可行域。
24. **`ex24_rydberg_blockade_dynamics.f90`**
    - **87Rb 70S 阻塞半径、二原子双激发抑制与 1D 阵列量子多体疤痕动力学**：模拟主量子数 $n=70$ 里德堡原子 $C_6/R^6$ 范德华相互作用导致的激光激发阻塞，演示两原子能级阻塞抑制，并在 1D Rydberg 链上求解 PXP 受限哈密顿量，重现高纠缠下的非热化量子多体疤痕相干振荡。
25. **`ex25_surface_corrugated_diffraction.f90`**
    - **He/LiF(001) 2D 晶格衍射谱、选择性吸附共振 SAR 与 Debye-Waller 扫描**：模拟热中性原子掠射至硬波纹周期表面晶格 (HCS)，计算程函近似 2D Bragg 衍射分支几率、表面束缚能级选择性吸附共振 (SAR) Fano 线型与晶格温度依赖的 Debye-Waller 声子非弹性衰减。
26. **`ex26_eley_rideal_surface_reaction.f90`**
    - **H+H/Cu(111) 气-固 ER 催化反应、超热放热能量分配与高振动激发态反转分布**：构建二维反应表面，模拟入射气相原子与表面吸附质直接碰撞反应的 Eley-Rideal 机理，输出超热反应焓在产物振动、转动与平动能中的非统计分配，并展现强烈的振动态布居反转 (Population Inversion)。
27. **`ex27_surface_electronic_friction_gle.f90`**
    - **NO/Au(111) 表面非绝热散射、广义朗之万 GLE 碰撞轨迹与电子-空穴对能损**：采用局域密度摩擦近似 (LDFA) 计算金属表面自由电子气摩擦张量，运行辛步进广义朗之万方程 (GLE) 轨迹，定量计算分子碰撞过程中电子-空穴对激发所造成的超快能量耗散与振动弛豫寿命。
28. **`ex28_grazing_fast_atom_diffraction.f90`**
    - **keV He 掠入射快原子表面量子衍射 (GIFAD)、经典彩虹角与亚皮米波纹度反演**：模拟高能快原子掠入射解耦，分离轴向高速沟道平动与亚电子伏横向量子干涉运动，计算经典表面彩虹偏转角，并通过彩虹极大值解析重构晶体表面亚皮米量级的微观波纹振幅。
29. **`ex29_cold_ion_atom_scattering.f90`**
    - **Yb+/Li 与 Ba+/Rb 冷离子-原子杂化碰撞截面与 Paul 阱微运动致热率**：模拟 $1/r^4$ 极化势引导的超冷离子-原子混合系统，解算有效微运动致热速率、临界 Langevin 碰撞参数与微开尔文极限温度。
30. **`ex30_tully_surface_hopping.f90`**
    - **Tully SAC 与 DAC 双避免交叉斯托克斯干涉系综非绝热动力学**：运行最少开关表面跳跃 (FSSH) 蒙特卡洛系综，计算核动量重标度、非绝热跃迁概率与斯托克斯量子相干干涉条纹。
31. **`ex31_molecular_alignment_revival.f90`**
    - **N2 与 CO2 飞秒无场转动复苏与光学离心机超转子加速**：模拟强激光场诱导分子非绝热取向对齐，追踪全复苏与分式复苏波包演化，并演示光学离心机角加速度驱动至超转子态引起离心破键。
32. **`ex32_optical_lattice_bose_hubbard.f90`**
    - **87Rb 光晶格深度扫描、Bose-Hubbard U/J 相变与引力 Bloch 振荡**：数值求解 Mathieu 方程计算晶格能带与 Wannier 局域化基底，映射 Bose-Hubbard 跃迁参数 $J$ 与在位相互作用 $U$，模拟引力场驱动的布洛赫振荡。
33. **`ex33_rph_variational_transition_state.f90`**
    - **多原子反应路径变分 CVT 速率与 Eckart 量子隧穿增强因子**：沿着内禀反应坐标 (IRC) 寻优自由能瓶颈，计算正则变分过渡态理论 (CVT) 反应速率与 Eckart 势垒半经典隧穿修正因子 $\kappa(T)$。
34. **`ex34_relativistic_dirac_cesium.f90`**
    - **铯原子 6s/6p/5d 相对论狄拉克能级、精细结构分裂与 D1/D2 振子强度**：数值求解全相对论径向狄拉克方程，精确计算自旋-轨道耦合天然精细结构裂分与电偶极吸收振子强度。
35. **`ex35_rixs_core_level_spectroscopy.f90`**
    - **铜氧化物 Cu L3 共振非弹性 X 射线散射 2D 能损图谱与声子伴峰**：利用 Kramers-Heisenberg 二阶截面公式模拟同步辐射共振非弹性 X 射线散射，生成入射能量-能量损失 2D 光谱，并计算 Huang-Rhys 声子级数。

---

## 🐍 Python 辅助分析套件 (`pygenmod`)

提供轻量 Python 库，用于数据交互、前处理计算与出版级可视化：

```bash
cd GeneralModule/python
python3 test_pygenmod.py   # 运行 12 大单元测试 (100% Pass)
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

### 6. 外场超冷原子碰撞与 Breit-Rabi / 磁 Feshbach 绘图 (`field_scattering.py`)
计算碱金属单原子在任意磁场下的 Zeeman-超精细 Breit-Rabi 能级劈裂，或构建四大基组（非耦合、f-耦合、总自旋、场缀饰）幺正变换矩阵并绘制磁 Feshbach 共振色散：
```python
from pygenmod import (
    get_cold_atom_preset, calc_breit_rabi_energies,
    build_field_collision_channels, calc_basis_transform_matrix,
    BASIS_UNCOUPLED, BASIS_TOTAL_SPIN,
    plot_breit_rabi_diagram, plot_magnetic_feshbach_resonance
)
import numpy as np

# 1. 绘制 87Rb 单原子 Breit-Rabi 能级图 (0 ~ 200 Gauss)
rb87 = get_cold_atom_preset("87Rb")
plot_breit_rabi_diagram(rb87, b_max_gauss=200.0, save_path="breit_rabi_rb87.png")

# 2. 四大基组幺正变换矩阵计算 (例如: 非耦合基组 -> 总自旋耦合基组)
ch_unc = build_field_collision_channels(rb87, rb87, BASIS_UNCOUPLED, two_Mtot=2, l_max=0)
ch_spin = build_field_collision_channels(rb87, rb87, BASIS_TOTAL_SPIN, two_Mtot=2, l_max=0)
U_spin_unc = calc_basis_transform_matrix(rb87, rb87, ch_unc, ch_spin, BASIS_UNCOUPLED, BASIS_TOTAL_SPIN)

# 3. 绘制磁 Feshbach 共振色散曲线 a_s(B) 与解析拟合
b_grid = np.linspace(50.0, 110.0, 100)
a_s = 100.0 * (1.0 - 5.0 / (b_grid - 80.0))  # 示例共振峰
plot_magnetic_feshbach_resonance(b_grid, a_s, save_path="feshbach_resonance_fit.png")
```

---

## 📚 科学文献典藏与权威理论全典

本算法库的所有物理模型、数值微分/积分格式、渐近边界匹配与对角化算法均严格对齐国际主流顶级物理期刊（PRL, PRA, PR, RMP, JCP, CPC 等）经典文献。

项目根目录下建立了完整的独立文献全典：👉 **[LITERATURE.md](LITERATURE.md)**

### 核心涵盖的物理领域与经典学术奠基：
1. **非含时散射理论、分波相移与有效力程展开**
   - E. P. Wigner, *Phys. Rev.* **73**, 1002 (1948) [分波相移与低能极限]
   - H. A. Bethe, *Phys. Rev.* **76**, 38 (1949) [有效力程展开 ERE]
   - E. P. Wigner, *Phys. Rev.* **98**, 145 (1955); F. T. Smith, *Phys. Rev.* **118**, 349 (1960) [Wigner-Smith 时延矩阵]
2. **多通道密耦对数导数法与多扇区分段网格**
   - B. R. Johnson, *J. Comput. Phys.* **13**, 445 (1973); *J. Chem. Phys.* **67**, 4086 (1977) [矩阵比值递推]
   - D. E. Manolopoulos, *J. Chem. Phys.* **85**, 6425 (1986); *J. Comput. Phys.* **105**, 169 (1993) [分段扇区局域传播]
   - J. M. Hutson & C. R. Le Sueur, *Comput. Phys. Commun.* **241**, 9 (2019) [MOLSCAT 分子密耦系统]
3. **长程范德华色散与半经典解析散射长度**
   - G. F. Gribakin & V. V. Flambaum, *Phys. Rev. A* **48**, 546 (1993) [平均散射长度 $\bar{a}$]
   - Bo Gao, *Phys. Rev. A* **58**, 4222 (1998); *Phys. Rev. A* **72**, 042719 (2005) [纯范德华解析解与 MQDT]
4. **超冷原子自旋相互作用与外场磁 Feshbach 共振**
   - H. Feshbach, *Ann. Phys.* **5**, 357 (1958); U. Fano, *Phys. Rev.* **124**, 1866 (1961) [共振与组态相互作用]
   - H. T. C. Stoof et al., *Phys. Rev. B* **38**, 4688 (1988); E. Tiesinga et al., *Phys. Rev. A* **47**, 4114 (1993) [自旋交换]
   - C. Chin, R. Grimm, P. S. Julienne, and E. Tiesinga, *Rev. Mod. Phys.* **82**, 1225 (2010) [超冷原子 Feshbach 共振全景综述]
5. **塞曼-超精细 Breit-Rabi 能谱与四大经典基组变换**
   - G. Breit & I. I. Rabi, *Phys. Rev.* **38**, 2082 (1931) [单电子 Breit-Rabi 解析公式]
   - D. A. Varshalovich et al., *Quantum Theory of Angular Momentum*, World Scientific (1988) [Racah 角动量耦合代数]
   - J. P. Burke, Jr., Ph.D. thesis, Univ. of Colorado (1999) [四大基组表象与幺正投影]
6. **强场超快物理、高次谐波发射与隧穿电离**
   - L. V. Keldysh, *Sov. Phys. JETP* **20**, 1307 (1965); M. V. Ammosov et al., *Sov. Phys. JETP* **64**, 1191 (1986) [ADK]
   - P. B. Corkum, *Phys. Rev. Lett.* **71**, 1994 (1993); M. Lewenstein et al., *Phys. Rev. A* **49**, 2117 (1994) [SFA 强场近似]
7. **离散变量表象 (DVR)、虚时间与分裂算符波包动力学**
   - D. T. Colbert & W. H. Miller, *J. Chem. Phys.* **96**, 1982 (1992) [Sinc-DVR]
   - C. C. Marston & G. G. Balint-Kurti, *J. Chem. Phys.* **91**, 3571 (1989) [Fourier Grid Hamiltonian, FGH]
   - M. D. Feit, J. A. Fleck, Jr., & A. Steiger, *J. Comput. Phys.* **47**, 412 (1982) [Split-Operator 算法]
8. **开放量子系统 Lindblad 耗散与 Krotov 最优控制**
   - G. Lindblad, *Commun. Math. Phys.* **48**, 119 (1976); V. Gorini et al., *J. Math. Phys.* **17**, 821 (1976) [Lindblad 动力学]
   - V. F. Krotov, *Global Methods in Optimal Control Theory* (1996); R. Somlói et al., *Chem. Phys.* **172**, 85 (1993) [Krotov 算法]
9. **各向异性磁偶极/电偶极超冷散射与自旋弛豫**
   - S. Hensler et al., *Appl. Phys. B* **77**, 765 (2003); J. Stuhler et al., *Phys. Rev. Lett.* **95**, 150406 (2005) [磁偶极散射与超冷自旋弛豫]
   - M. Marinescu & L. You, *Phys. Rev. Lett.* **81**, 4596 (1998) [超冷极性分子电偶极-偶极相互作用 EDDI 与 Stark 能级]
   - K.-K. Ni et al., *Science* **322**, 231 (2008); S. Ospelkaus et al., *Science* **327**, 853 (2010) [超冷极性双原子分子碰撞与化学反应控制]
10. **超冷原子光缔合谱学与自由-束缚态分子生成**
    - H. R. Thorsheim, J. Weiner, & P. S. Julienne, *Phys. Rev. Lett.* **58**, 2420 (1987) [超冷原子光缔合 PA 奠基之作]
    - K. M. Jones, E. Tiesinga, P. D. Lett, & P. S. Julienne, *Rev. Mod. Phys.* **78**, 483 (2006) [超冷原子光缔合光谱学权威综述]
    - J. G. Danzl et al., *Science* **321**, 1062 (2008); *Nat. Phys.* **6**, 265 (2010) [双光子 STIRAP 绝热受激跃迁制备超冷振转基态分子]

---

## 🤖 AI 编程与 Agent 使用指南

若您在日常开发中配合 AI 编程助手（如 Antigravity, Claude, Copilot）工作，可以直接指示 AI：

> *"请使用 `general_module` 编写一个 Fortran 程序，调用 `dvr_sinc_init` 和 `fgh_solve_bound_states` 求解一维双阱势 $V(x) = x^4 - 2x^2$ 的前 4 个本征态与能级劈裂，并用 `to_au` 统一输入物理单位。"*

AI 将直接利用自包含的类型定义和确定的函数签名生成无错代码，无需额外的繁琐胶水代码或外部环境配置。
