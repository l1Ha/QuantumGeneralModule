# GeneralModule: 现代量子动力学通用算法库 (Fortran 2008 / Python)

[![CI](https://github.com/l1Ha/QuantumGeneralModule/actions/workflows/ci.yml/badge.svg)](https://github.com/l1Ha/QuantumGeneralModule/actions/workflows/ci.yml)
[![Fortran 2008](https://img.shields.io/badge/Fortran-2008-734f96.svg)](https://fortran-lang.org/)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Tests: 340/340 Pass](https://img.shields.io/badge/Tests-340%2F340%20Pass%20(100%25)-brightgreen.svg)](tests/)
[![Literature: 42 Topics](https://img.shields.io/badge/Literature-42%20Topics%20(PRL%2FPRA%2FRMP%2FScience%2FNature)-blue.svg)](LITERATURE.md)

`GeneralModule` 是一个面向超快强场物理、分子光物理、前沿量子散射与表面动力学模拟的现代化通用科学计算算法库。该库遵循严格的 **Fortran 2008 规范**，具备高数值精度、零外部动态库强依赖、模块化架构与出色的 AI Agent 友好性。

---

## 🌟 核心特性

1. **零外部库依赖 (Zero External Dependencies)**
   - 内部集成高精度 Householder QL 实对称矩阵本征求解器、Gauss-Jordan 全主元实/复方阵求逆与 Cooley-Tukey 1D/2D 快速傅里叶变换（FFT）。
   - 纯 Fortran 自包含样条插值、Lindblad 主方程、Krotov 最优控制、通用多通道定态密耦（Johnson Log-Derivative）、外场多基组散射、各向异性偶极耦合、超冷光缔合速率、少体 Efimov 物理、低维光晶格 CIR、自电离 Fano/CCR、交叉电磁场、三原子反应 PES、旋量 BEC 自旋动力学、三原子超球面反应动力学、偶极量子液滴 LHY、强场非顺序双电离 (NSDI)、磁/光 Feshbach 束缚态、阿秒瞬态吸收光谱 (ATAS)、双色反向圆偏振 PECD、超冷极性分子偶极遮蔽、里德堡原子阻塞、2D 表面量子散射与 SAR、气-固催化 Eley-Rideal 反应、金属表面非绝热电子摩擦 (GLE)、掠入射快原子衍射 (GIFAD)、冷离子-中性原子杂化极化散射、Tully 最少开关表面跳跃 (FSSH)、强场分子定向与超转子动力学、光晶格 Bose-Hubbard 映射、反应路径哈密顿量 (RPH/CVT)、相对论径向狄拉克方程、共振非弹性 X 射线散射 (RIXS) 以及亚稳态原子碰撞潘宁电离与缔合电离 (Penning & Associative Ionization) 求解器，无需强制链接外部 LAPACK/BLAS 或 FFTW，开箱即用。
2. **现代 Fortran 2008 标准设计**
   - 统一强类型参数定义（`real(dp) => real64`）。
   - 纯函数（`pure function`）与显式 `intent(in/out/inout)` 契约，杜绝隐式全局变量副作用。
3. **AI 友好型结构化接口 (AI-Friendly)**
   - 算法模块支持统一顶层聚合入口：`use general_module`。
   - 参数配置采用清晰的派生类型（Derived Types，如 `pulse_config_t`, `surface_lattice_t`, `ion_atom_system_t`, `tully_model_t`, `rotor_molecule_t`, `optical_lattice_t`, `rph_path_t`, `dirac_state_t`, `rixs_system_t`, `penning_system_t` 等），自解释、低耦合、便于大语言模型精确构造与调用。
4. **全链路双语生态支持**
   - 附带标准 Python 伴侣分析包 `pygenmod`，无缝衔接参数预计算、波包与散射长度可视化、Breit-Rabi 能级图、发表级绘图（含一键动力学出图流水线 `plot_rovibrational_dynamics.py`）。
5. **全自动 CI/CD 持续集成**
   - 内置 GitHub Actions 跨平台持续集成（Ubuntu / macOS），全自动化执行 41 大测试套件（340 个单元断言 100% 通过）与 36 大物理应用工程算例。

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
│   ├── mod_penning_associative_ionization.f90 # 49. 亚稳态原子潘宁电离与缔合化学电离动力学、光学势自电离宽度、PIES 电子发射能谱与自旋抑制
│   └── general_module.f90          # 顶层聚合入口模块 (use general_module, 49 大核心物理模块)
├── tests/                          # 自动化单元测试套件 (41 个套件，100% 全部通过，340/340 断言)
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
│   ├── test_penning_associative_ionization.f90 # 潘宁电离与缔合电离势能面、自电离宽度、低能 AI 俘获主导至高能 PI 转变与 PIES 能谱测试
│   └── run_all_tests.sh            # 自动化测试运行脚本 (100% Pass, 340/340 断言, 41 测试套件)
├── examples/                       # 典型物理应用算例 (36 大完整前沿算例)
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
│   ├── ex36_penning_associative_ionization.f90  # 亚稳态 He*(2^3S)+Ar 潘宁电离与缔合电离空间分流比、PIES 电子能谱与自旋极化抑制
│   └── build_examples.sh           # 算例编译运行脚本 (全 36 算例编译运行通过)
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
- **理论基础与物理机制**：采用 CODATA 2018/2022 基础物理常数推荐值。微观原子分子与超快强场物理计算以原子单位制（Hartree a.u.）为基准（约定 $\hbar = m_e = e = 4\pi\epsilon_0 = 1$）。
- **详细数学表达式**：
  - 玻尔半径与长度单位：$1\text{ a.u.} = a_0 = \frac{4\pi\epsilon_0 \hbar^2}{m_e e^2} \approx 0.5291772109\text{ \AA} = 5.291772109 \times 10^{-11}\text{ m}$
  - 哈特里能量与能量转换：$1\text{ a.u.} = E_h = \frac{\hbar^2}{m_e a_0^2} \approx 27.211386246\text{ eV} \approx 4.359744722 \times 10^{-18}\text{ J} \approx 219474.63\text{ cm}^{-1} \approx 3.15775 \times 10^5\text{ K}$
  - 原子时间尺度：$1\text{ a.u.} = \tau_0 = \frac{\hbar}{E_h} \approx 0.024188843265\text{ fs} = 2.4188843265 \times 10^{-17}\text{ s}$
  - 原子电场强度：$1\text{ a.u.} = \mathcal{E}_0 = \frac{E_h}{e a_0} \approx 5.1422067 \times 10^{11}\text{ V/m} \approx 51.422\text{ MV/cm}$
  - 强场激光峰值光强：$I_{\text{a.u.}} = \frac{1}{2}\epsilon_0 c \mathcal{E}_0^2 \approx 3.509445 \times 10^{16}\text{ W/cm}^2$（即 $E_0 = \sqrt{I / (3.5094 \times 10^{16}\text{ W/cm}^2)}$）
  - 磁感应强度：$1\text{ a.u.} = B_0 = \frac{\hbar}{e a_0^2} \approx 2.350517568 \times 10^5\text{ Tesla} \approx 2.3505 \times 10^9\text{ Gauss}$
- **核心 API 映射**：
  - `to_au(val, unit)`: 将实验单位数值（fs, ps, s, eV, cm-1, J, K, Angstrom, nm, m, amu, MV/cm, Debye, W/cm2, Gauss, Tesla）严格转换为原子单位。
  - `from_au(val_au, unit)`: 将原子单位数值转换回国际标准或实验光谱常用单位。
  - 基础常数：`PI`, `TWOPI`, `HALFPI`, `SQRTPI`, `EYE`, `C_LIGHT`, `HBAR`, `M_E`, `CHARGE_E`, `KB`, `AMU2AU`, `AU2EV`, `EV2AU`, `AU2ANG`, `ANG2AU`, `AU2FS`, `FS2AU`。

### 2. 量子力学特殊函数与角动量代数 (`mod_special_functions`)
- **理论基础与物理机制**：量子体系的空间旋转对称性、球对称中心力场展开以及多角动量耦合严格遵循 $SO(3)$ 与 $SU(2)$ 李代数。
- **详细数学表达式**：
  - 勒让德多项式与缔合勒让德函数（含 Condon-Shortley 相位 $(-1)^m$）：
    $$P_l(x) = \frac{1}{2^l l!} \frac{d^l}{dx^l}(x^2 - 1)^l, \quad P_l^m(x) = (-1)^m (1 - x^2)^{m/2} \frac{d^m}{dx^m} P_l(x)$$
  - 正则球面调和函数：$Y_{lm}(\theta, \phi) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}} P_l^m(\cos\theta) e^{i m \phi}$
  - Wigner 3j 符号与 Clebsch-Gordan 耦合系数恒等式：
    $$\langle j_1 m_1 j_2 m_2 | j_3 m_3 \rangle = (-1)^{j_1 - j_2 + m_3} \sqrt{2j_3 + 1} \begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}$$
    严格服从三角定则 $|j_1 - j_2| \le j_3 \le j_1 + j_2$ 与磁量子数加和守恒 $m_1 + m_2 = m_3$。
  - 刚体转子取向矩阵元（基于 Wigner-Eckart 定理求得的精确解析闭式解）：
    $$\langle j, m | \cos\theta | j', m \rangle = \sqrt{\frac{2j'+1}{2j+1}} \langle j' m 1 0 | j m \rangle \langle j' 0 1 0 | j 0 \rangle = \begin{cases} \sqrt{\frac{j^2 - m^2}{(2j-1)(2j+1)}}, & j' = j - 1 \\ \sqrt{\frac{(j+1)^2 - m^2}{(2j+1)(2j+3)}}, & j' = j + 1 \\ 0, & \text{otherwise} \end{cases}$$
    $$\langle j, m | \cos^2\theta | j', m \rangle = \frac{1}{3}\delta_{j j'} + \frac{2}{3}\sqrt{\frac{2j'+1}{2j+1}} \langle j' m 2 0 | j m \rangle \langle j' 0 2 0 | j 0 \rangle$$
- **核心 API 映射**：`legendre_poly`, `assoc_legendre_poly`, `spherical_harmonic`, `wigner_3j`, `clebsch_gordan`, `wigner_6j`, `wigner_9j`, `rot_matrix_cos_theta`, `rot_matrix_cos2_theta`。

### 3. 线性代数与快速傅里叶变换 (`mod_linear_algebra`)
- **理论基础与物理机制**：量子离散哈密顿算符的本征能级求解对应实对称稠密矩阵谱分解，采用 Householder 正交相似变换降低带宽，结合隐式位移 QL 迭代收敛；动量与坐标空间表象转换基于分治 Cooley-Tukey 快速傅里叶变换。
- **详细数学表达式**：
  - Householder 镜像反射正交矩阵：
    $$\mathbf{P}_k = \mathbf{I} - 2 \frac{\mathbf{u}_k \mathbf{u}_k^T}{\mathbf{u}_k^T \mathbf{u}_k}, \quad \mathbf{T} = \mathbf{P}_{n-2} \cdots \mathbf{P}_1 \mathbf{A} \mathbf{P}_1 \cdots \mathbf{P}_{n-2}$$
    将 $n \times n$ 实对称矩阵 $\mathbf{A}$ 严格正交相似变换为三对角矩阵 $\mathbf{T}$。
  - 隐式位移 QL 迭代算法：
    $$\mathbf{T}_k - s_k \mathbf{I} = \mathbf{Q}_k \mathbf{L}_k \implies \mathbf{T}_{k+1} = \mathbf{L}_k \mathbf{Q}_k + s_k \mathbf{I} \xrightarrow{k \to \infty} \mathbf{\Lambda} = \text{diag}(\lambda_1, \lambda_2, \dots, \lambda_n)$$
  - Cooley-Tukey 1D 离散傅里叶变换（蝶形运算）：
    $$X_k = \sum_{n=0}^{N-1} x_n \exp\left( -i \frac{2\pi k n}{N} \right) = E_k + e^{-i \frac{2\pi k}{N}} O_k, \quad k = 0, \dots, N/2 - 1$$
- **核心 API 映射**：`diag_symmetric_matrix`（升序输出全部本征值与本征向量）, `inv_real_matrix`, `inv_complex_matrix`, `fft_1d`, `fft_2d`。

### 4. 离散变量表象与谱方法网格 (`mod_dvr_grid`)
- **理论基础与物理机制**：离散变量表象 (Discrete Variable Representation, DVR) 将无限维量子算符在局域正交格点空间投影，使得任意局域势能算符严格对角化 $V_{ij} = V(x_i)\delta_{ij}$，而动能算符具有全局解析形式。
- **详细数学表达式**：
  - Colbert-Miller Sinc-DVR 基函数定义：
    $$\theta_i(x) = \frac{1}{\sqrt{\Delta x}} \text{sinc}\left(\frac{\pi(x - x_i)}{\Delta x}\right) = \frac{\sin[\pi(x - x_i)/\Delta x]}{\pi (x - x_i)/\sqrt{\Delta x}}$$
  - 动能算符解析矩阵元（二阶微分算子）：
    $$T_{ij} = -\frac{\hbar^2}{2m} \int_{-\infty}^\infty \theta_i(x) \frac{d^2}{dx^2} \theta_j(x) dx = \frac{\hbar^2}{2m \Delta x^2} \begin{cases} \frac{\pi^2}{3}, & i = j \\ \frac{2(-1)^{i-j}}{(i - j)^2}, & i \ne j \end{cases}$$
  - Gauss-Legendre DVR 节点 $x_i$ 与积分权重：$P_N(x_i) = 0, \quad w_i = \frac{2}{(1 - x_i^2)[P_N'(x_i)]^2}$
  - Fourier Grid Hamiltonian (FGH) 束缚态本征方程：
    $$\sum_{j=1}^N \left( T_{ij} + V(x_i)\delta_{ij} \right) \psi_j^{(n)} = E_n \psi_i^{(n)}$$
- **核心 API 映射**：`dvr_sinc_init`, `dvr_legendre_init`, `fgh_solve_bound_states`。

### 5. 激光脉冲合成与动力学 Stark 效应 (`mod_laser_pulse`)
- **理论基础与物理机制**：超短强激光脉冲与物质相互作用，包含时变包络、载波包络相位（CEP）、线性啁啾、双色反向旋转场以及非共振 AC Stark 诱导能级移动。
- **详细数学表达式**：
  - 瞬时电场与包络函数：
    $$E(t) = E_0 f(t) \cos(\omega(t)(t - t_0) + \phi_{\text{CEP}})$$
    - 高斯包络：$f(t) = \exp\left(-2\ln 2 \frac{(t - t_0)^2}{\tau^2}\right)$（$\tau$ 为半高全宽 FWHM）
    - $\sin^2$ 包络：$f(t) = \sin^2\left(\frac{\pi t}{T}\right) \Theta(t)\Theta(T - t)$
  - 线性频率啁啾：$\omega(t) = \omega_0 + \beta (t - t_0)$，矢势 $A(t) = -\int_{-\infty}^t E(t') dt'$
  - 双色合成相干光场：$E(t) = E_1 f_1(t) \cos(\omega t) + E_2 f_2(t) \cos(2\omega t + \phi_{12})$
  - 动力学极化张量 AC Stark 位移：
    $$\Delta E_{\text{Stark}}(t, \theta) = -\frac{1}{2} E(t)^2 \left[ \Delta\alpha \cos^2\theta + \alpha_\perp \right], \quad \Delta\alpha = \alpha_\parallel - \alpha_\perp$$
- **核心 API 映射**：`pulse_envelope`, `pulse_electric_field`, `pulse_vector_potential`, `pulse_stark_shift`, `pulse_generate_timeseries`。

### 6. 复吸收势边界与量子概率流密度 (`mod_absorbing_boundary`)
- **理论基础与物理机制**：在有限坐标网格内演化波包时，向外扩散或电离的波包会在边界产生非物理反射；引入复吸收势 (Complex Absorbing Potential, CAP) 破坏时间反演对称性并平滑吸收外行波。
- **详细数学表达式**：
  - 有效非厄米哈密顿量：$\hat{H}_{\text{eff}} = \hat{H}_0 - i W(r)$
  - 多项式型与 $\sin^2$ 型 CAP：
    $$W_{\text{poly}}(r) = \eta \left(\frac{r - r_{\text{start}}}{r_{\text{end}} - r_{\text{start}}}\right)^n \Theta(r - r_{\text{start}}), \quad W_{\sin^2}(r) = \eta \sin^2\left(\frac{\pi(r - r_{\text{start}})}{2(r_{\text{end}} - r_{\text{start}})}\right)$$
  - 平滑吸收掩膜作用算符：$\psi(r, t + \Delta t) \leftarrow \psi(r, t) \cdot \exp\left( - \frac{W(r)\Delta t}{\hbar} \right)$
  - 量子概率流密度与连续性方程衰减律：
    $$\mathbf{j}(r, t) = \frac{\hbar}{\mu} \text{Im}\left[ \psi^*(r, t) \nabla \psi(r, t) \right], \quad \frac{\partial |\psi|^2}{\partial t} + \nabla \cdot \mathbf{j} = -\frac{2}{\hbar} W(r) |\psi|^2$$
- **核心 API 映射**：`cap_init`, `cap_evaluate`, `cap_apply_mask`, `calculate_probability_flux`, `calculate_norm_inside`。

### 7. 热统计力学与系综平均 (`mod_thermal_ensemble`)
- **理论基础与物理机制**：处于热平衡温度 $T$ 下的气相分子或粒子处于正则统计系综，各量子态按玻尔兹曼因子分配权重，宏观可观测量由量子力学期待值的系综统计加权给出。
- **详细数学表达式**：
  - 刚体转动配分函数与初态权重（考虑 $(2J+1)$ 空间简并度）：
    $$Z_{\text{rot}}(T) = \sum_{J=0}^{J_{\max}} (2J+1) \exp\left( - \frac{B J(J+1)}{k_B T} \right), \quad w_J(T) = \frac{(2J+1)}{Z_{\text{rot}}(T)} \exp\left( - \frac{B J(J+1)}{k_B T} \right)$$
  - 简谐振动配分函数：$Z_{\text{vib}}(T) = \sum_{v=0}^{v_{\max}} \exp\left( - \frac{\hbar\omega_e (v + 1/2)}{k_B T} \right) = \frac{e^{-\hbar\omega_e/(2k_B T)}}{1 - e^{-\hbar\omega_e/(k_B T)}}$
  - 热系综物理可观测量加权平均：
    $$\langle \hat{O} \rangle(T, t) = \sum_{J=0}^{J_{\max}} w_J(T) \langle \psi_J(t) | \hat{O} | \psi_J(t) \rangle$$
  - Bose-Einstein 玻色子环境统计因子：$n_{\text{BE}}(\omega, T) = \left[ \exp\left(\frac{\hbar\omega}{k_B T}\right) - 1 \right]^{-1}$
- **核心 API 映射**：`boltzmann_rotational_weights`, `boltzmann_vibrational_weights`, `thermal_average_1d`, `thermal_average_2d`, `bose_einstein_factor`。

### 8. 波包推进器与光学 Bloch 方程 (`mod_wavepacket_propagator`)
- **理论基础与物理机制**：含时薛定谔方程严格酉算符推进采用二阶辛对称 Strang 分裂算符算法；二能级系统相干驱动与耗散采用保模长光学 Bloch 矢量推进。
- **详细数学表达式**：
  - 2 阶辛对称 Strang 分裂算符（Trotter 分解）：
    $$\hat{U}(\Delta t) = \exp\left( -\frac{i\hat{H}\Delta t}{\hbar} \right) = \exp\left( -\frac{i\hat{V}\Delta t}{2\hbar} \right) \exp\left( -\frac{i\hat{T}\Delta t}{\hbar} \right) \exp\left( -\frac{i\hat{V}\Delta t}{2\hbar} \right) + \mathcal{O}(\Delta t^3)$$
    其中动能演化在动量空间通过 FFT 对角作用：$\tilde{\psi}(p) = \mathcal{F}[\psi(x)], \quad \tilde{\psi}'(p) = \tilde{\psi}(p) e^{-i \frac{p^2}{2m}\frac{\Delta t}{\hbar}}$。
  - 光学 Bloch 矢量方程（$\mathbf{R} = [u, v, w]^T$）：
    $$\frac{du}{dt} = -\Delta v - \frac{u}{T_2}, \quad \frac{dv}{dt} = \Delta u + \Omega_R(t) w - \frac{v}{T_2}, \quad \frac{dw}{dt} = -\Omega_R(t) v - \frac{w - w_0}{T_1}$$
    其中 $u = 2\text{Re}(\rho_{12})$ 为同相色散分量，$v = 2\text{Im}(\rho_{21})$ 为正交吸收分量，$w = \rho_{22} - \rho_{11}$ 为反转粒子数差。
- **核心 API 映射**：`propagate_split_operator_1d`, `propagate_split_operator_2d`, `rk4_step`, `solve_bloch_two_level`, `abm4_step`。

### 9. 强场原子模型与隧穿电离 (`mod_coulomb_atomic`)
- **理论基础与物理机制**：在超快强激光场（$I \ge 10^{14}\text{ W/cm}^2$）下，外电场强度可与原子核库仑场比拟，束缚电子发生非微扰隧穿电离；采用单活性电子（SAE）模型势与经典三步模型标度。
- **详细数学表达式**：
  - 一维软核库仑模型势与受力：
    $$V_{\text{soft}}(x) = -\frac{Z_{\text{eff}}}{\sqrt{x^2 + a^2}}, \quad F_{\text{soft}}(x) = -\frac{dV_{\text{soft}}}{dx} = -\frac{Z_{\text{eff}} x}{(x^2 + a^2)^{3/2}}$$
  - Keldysh 绝热参数（划分多光子电离 $\gamma \gg 1$ 与准静态隧穿电离 $\gamma \ll 1$）：
    $$\gamma = \frac{\omega \sqrt{2 I_p}}{F_0} = \sqrt{\frac{I_p}{2 U_p}}, \quad U_p = \frac{F_0^2}{4\omega^2} = \frac{e^2 \mathcal{E}_0^2}{4 m_e \omega^2}$$
  - 高次谐波与高能重碰撞电子截止能量定律：$E_{\text{cutoff}} = I_p + 3.17 U_p$
  - 准静态 ADK (Ammosov-Delone-Krainov) 隧穿电离率公式：
    $$W_{\text{ADK}}(F) = C_{n^* l}^2 f(l, m) I_p \left( \frac{2(2I_p)^{3/2}}{F} \right)^{2n^* - |m| - 1} \exp\left( - \frac{2(2I_p)^{3/2}}{3 F} \right)$$
    其中有效主量子数 $n^* = Z_{\text{eff}} / \sqrt{2 I_p}$，系数 $C_{n^* l}^2 = \frac{2^{2n^*}}{n^* \Gamma(n^* + l^* + 1) \Gamma(n^* - l^*)}$。
- **核心 API 映射**：`get_atom_config`, `soft_core_coulomb_potential`, `soft_core_coulomb_derivative`, `keldysh_parameter`, `ponderomotive_energy`, `hhg_cutoff_energy`, `adk_ionization_rate`。

### 10. 高次谐波与偶极时频分析 (`mod_hhg_spectra`)
- **理论基础与物理机制**：强场高次谐波发射（HHG）来源于电离电子在激光场中的三步物理过程：隧穿电离、激光加速与回碰复合；瞬时偶极辐射通过 Ehrenfest 定理或 Lewenstein 强场近似（SFA）求解。
- **详细数学表达式**：
  - 长度表象与加速度表象瞬时偶极响应（Ehrenfest 定理）：
    $$d(t) = \langle \psi(t) | \hat{x} | \psi(t) \rangle, \quad a(t) = \frac{d^2 d}{dt^2} = -\langle \psi(t) | \left( \frac{\partial V}{\partial x} + E(t) \right) | \psi(t) \rangle$$
  - Lewenstein SFA 强场近似鞍点重碰撞偶极矩闭合公式：
    $$d(t) = i \int_0^\infty d\tau \left( \frac{\pi}{\epsilon + i\tau/2} \right)^{3/2} d_x^*(p_{\text{st}}(t, \tau) + A(t)) e^{-i S(p_{\text{st}}, t, \tau)} E(t - \tau) d_x(p_{\text{st}}(t, \tau) + A(t - \tau)) + \text{c.c.}$$
    其中准经典准自由连续态平稳动量为 $p_{\text{st}}(t, \tau) = -\frac{1}{\tau}\int_{t-\tau}^t A(t') dt'$，准经典作用量为 $S(p_{\text{st}}, t, \tau) = \int_{t-\tau}^t \left[ \frac{(p_{\text{st}} + A(t'))^2}{2} + I_p \right] dt'$。
  - 阿秒脉冲高次谐波辐射发射功率谱：
    $$S(\omega) = \left| \frac{1}{\sqrt{2\pi}} \int_0^T a(t) W_{\text{Hann}}(t) e^{-i \omega t} dt \right|^2$$
  - Gabor 小波时频变换（解析阿秒量子轨道与正反向啁啾路径）：
    $$G(t_0, \omega) = \int a(t) \exp\left( - \frac{(t - t_0)^2}{2\sigma^2} \right) e^{-i \omega t} dt$$
- **核心 API 映射**：`calculate_dipole_length`, `calculate_dipole_acceleration`, `hhg_power_spectrum`, `gabor_transform_point`, `lewenstein_sfa_dipole`。

### 11. 切比雪夫推进器与能谱滤波 (`mod_chebyshev_propagator`)
- **理论基础与物理机制**：切比雪夫多项式展开是全域全局近似演化算符 $\hat{U}(\Delta t) = e^{-i\hat{H}\Delta t/\hbar}$ 的最高精度方案，时间步长可跨越数百飞秒而保持机器精度的酉性与范数守恒；结合 Schafer-Kulander 能量窗算子可从单次含时波包演化中高精度滤波提取定态连续光电子动能谱 (PES)。
- **详细数学表达式**：
  - 谱重标度哈密顿量（将谱域映射至 $[-1, 1]$）：
    $$\hat{H}_{\text{norm}} = \frac{\hat{H} - \bar{E}}{\Delta E}, \quad \bar{E} = \frac{E_{\max} + E_{\min}}{2}, \quad \Delta E = \frac{E_{\max} - E_{\min}}{2}$$
  - 第一类切比雪夫多项式递推展开含时演化算符：
    $$e^{-i \hat{H} \Delta t / \hbar} |\psi(t)\rangle = e^{-i \bar{E} \Delta t / \hbar} \sum_{n=0}^M c_n\left(\frac{\Delta E \Delta t}{\hbar}\right) T_n(-i \hat{H}_{\text{norm}}) |\psi(t)\rangle$$
    其中递推基底为 $T_0(x) = 1, T_1(x) = x, T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)$，展开系数为第一类 Bessel 函数 $c_n(\alpha) = (2 - \delta_{n0}) (-i)^n J_n(\alpha)$。
  - Schafer-Kulander 能量窗投影算子（动能谱 PES 提取）：
    $$\hat{\mathcal{P}}(E_k, \gamma) = \frac{\gamma^{2^m}}{(\hat{H} - E_k)^{2^m} + \gamma^{2^m}} \implies P(E_k) = \langle \psi | \hat{\mathcal{P}}(E_k, \gamma) | \psi \rangle$$
- **核心 API 映射**：`chebyshev_propagate_step`, `window_operator_pes`。

### 12. 多势能面非绝热动力学 (`mod_multistate_coupling`)
- **理论基础与物理机制**：多电子态分子体系在避免交叉（Avoided Crossing）或锥形交叉区域，核运动与电子运动的 Born-Oppenheimer 绝热近似失效；核波包在多势能面间发生相干分束与无辐射非绝热跃迁。
- **详细数学表达式**：
  - 透热哈密顿量与绝热势能面变换：
    $$\mathbf{H}_{\text{dia}}(R) = \begin{pmatrix} V_{11}(R) & V_{12}(R) \\ V_{12}(R) & V_{22}(R) \end{pmatrix} \xrightarrow{\mathbf{U}(R)} \mathbf{V}_{\text{adia}}(R) = \begin{pmatrix} E_-(R) & 0 \\ 0 & E_+(R) \end{pmatrix}$$
    绝热本征能级为：$E_\pm(R) = \frac{V_{11} + V_{22}}{2} \pm \sqrt{\left(\frac{V_{11} - V_{22}}{2}\right)^2 + V_{12}^2}$，最小避免交叉能隙为 $\Delta E_{\min} = 2 |V_{12}(R_c)|$。
  - 经典 Landau-Zener 非绝热跃迁概率公式：
    $$P_{\text{LZ}} = \exp\left( - \frac{2\pi |V_{12}(R_c)|^2}{\hbar v |\Delta F|} \right), \quad \Delta F = \left.\left| \frac{dV_{11}}{dR} - \frac{dV_{22}}{dR} \right|\right|_{R = R_c}$$
  - 双通道非绝热核波包推进格式（利用解析矩阵指数）：
    $$\begin{pmatrix} \psi_1(t+\Delta t) \\ \psi_2(t+\Delta t) \end{pmatrix} = \exp\left( -\frac{i \mathbf{V}(R)\Delta t}{2\hbar} \right) \exp\left( -\frac{i \hat{T}\Delta t}{\hbar} \right) \exp\left( -\frac{i \mathbf{V}(R)\Delta t}{2\hbar} \right) \begin{pmatrix} \psi_1(t) \\ \psi_2(t) \end{pmatrix}$$
- **核心 API 映射**：`propagate_split_operator_2channel`, `landau_zener_probability`, `calculate_channel_populations`。

### 13. 分子转振耦合、STIRAP 与 Franck-Condon 谱学 (`mod_rovibrational`)
- **理论基础与物理机制**：双原子分子转动与振动自由度耦合，激光场诱导电偶极跃迁遵循角动量选择定则 $\Delta J = \pm 1$；利用受激拉曼绝热通道（STIRAP）通过相干暗态实现高保真度无激发态损耗的基态制备。
- **详细数学表达式**：
  - 双原子分子 Morse 势能与离心势能：
    $$V(R) = D_e \left[ 1 - e^{-a(R - R_e)} \right]^2, \quad V_{\text{eff}}(R) = V(R) + \frac{\hbar^2 J(J+1)}{2\mu R^2}$$
  - 转振态对角能级与态依有效转动常数：
    $$E(v, J) = E_v + B_v J(J+1), \quad B_v = \langle \chi_v | \frac{\hbar^2}{2\mu R^2} | \chi_v \rangle$$
  - Franck-Condon 因子与振动跃迁偶极矩：
    $$FC(v, v') = |\langle \chi_v | \chi_{v'} \rangle|^2, \quad M(v, v') = \langle \chi_v | \mu(R) | \chi_{v'} \rangle$$
  - STIRAP 三能级相干暗态（Dark State）无辐射传输：
    $$|D(t)\rangle = \cos\Theta(t) |1\rangle - \sin\Theta(t) |3\rangle, \quad \tan\Theta(t) = \frac{\Omega_P(t)}{\Omega_S(t)}$$
    在逆直觉时序（Stokes 脉冲 $\Omega_S(t)$ 先于 Pump 脉冲 $\Omega_P(t)$ 入射）下，体系严格沿着暗态演化，中间损耗激发态 $|2\rangle$ 始终零布居。
- **核心 API 映射**：`morse_potential`, `calc_franck_condon_factors`, `calc_vibrational_dipole_matrix`, `calc_rotational_constants_bv`, `build_rovibrational_hamiltonian`, `build_rovibrational_dipole_matrix`, `build_rovibrational_polarizability_matrix`, `create_stirap_pulses`, `rovibrational_state_index`。

### 14. 科学计算 I/O 与诊断工具 (`mod_io_utils`)
- **理论基础与物理机制**：科学数据流的无损持久化、多维时空演化矩阵导出与工业级运行诊断工具。
- **详细数学表达式**：
  - 动力学时序数据输出规范：$D_{ij} = y_j(t_i)$，自动校验格点单调性与浮点有效数字精度（双精度 `ES24.16E3`）。
  - 矩阵范数与保真度诊断：$\| \mathbf{A} \|_F = \sqrt{\sum_{i,j} |A_{ij}|^2}$。
- **核心 API 映射**：`save_data_table_1d`, `save_data_table_2d`, `save_matrix_dat`, `print_banner`, `print_progress_bar`。

### 15. 高精度样条插值与势能面渐近外推 (`mod_interpolation`)
- **理论基础与物理机制**：从离散从头算量子化学电子结构点重构光滑连续的势能面与力场，采用自然边界或导数钳位三次样条，并在物理极限区平滑连接短程斥力核与长程多极色散。
- **详细数学表达式**：
  - 三次样条段式插值多项式（$x \in [x_i, x_{i+1}]$）：
    $$S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3$$
    满足一阶导数连续 $S_i'(x_{i+1}) = S_{i+1}'(x_{i+1})$ 与二阶曲率连续 $S_i''(x_{i+1}) = S_{i+1}''(x_{i+1})$。
  - 三对角矩阵方程（Thomas 算法 $O(N)$ 极速求解）：
    $$h_{i-1} c_{i-1} + 2(h_{i-1} + h_i) c_i + h_i c_{i+1} = 3\left( \frac{y_{i+1} - y_i}{h_i} - \frac{y_i - y_{i-1}}{h_{i-1}} \right)$$
  - 渐近平滑接合势外推函数：
    $$V_{\text{extrap}}(R) = \begin{cases} A_{\text{rep}} e^{-B_{\text{rep}} R}, & R < R_{\min} \\ S(R), & R_{\min} \le R \le R_{\max} \\ V_\infty - \frac{C_6}{R^6} - \frac{C_8}{R^8}, & R > R_{\max} \end{cases}$$
- **核心 API 映射**：`spline_1d_init`, `spline_1d_eval`, `spline_1d_deriv`, `spline_1d_deriv2`, `potential_extrapolate_1d`。

### 16. 分子光解离动力学与碎片动能释放谱 (`mod_photofragment_flux`)
- **理论基础与物理机制**：光子激发中性分子至排斥态后发生单分子碎裂；根据波包自相关函数傅里叶变换解析光吸收谱，并在渐近反应通道监测量子概率流提取碎片动能释放谱 (Kinetic Energy Release, KER)。
- **详细数学表达式**：
  - 波包自相关函数与 Heller 连续吸收截面（光跃迁振子强度）：
    $$C(t) = \langle \psi(0) | \psi(t) \rangle, \quad \sigma_{\text{abs}}(\omega) \propto \omega \int_{-\infty}^\infty C(t) e^{i(E_0 + \hbar\omega)t/\hbar} e^{-\gamma |t|} dt$$
  - 碎片动能释放（KER）能量守恒律：
    $$E_{\text{KER}} = \hbar\omega - D_0 - E_{\text{int}}(A) - E_{\text{int}}(B)$$
  - 渐近边界监测面 $R_{\text{det}}$ 处时间-能量傅里叶散射振幅：
    $$A(E) = \frac{1}{\sqrt{2\pi\hbar}} \int_0^\infty \psi(R_{\text{det}}, t) e^{i E t/\hbar} dt, \quad \frac{dP}{dE_{\text{KER}}} \propto \frac{\hbar k_E}{\mu} |A(E_{\text{KER}})|^2$$
- **核心 API 映射**：`calc_autocorrelation`, `heller_absorption_spectrum`, `photofragment_energy_amplitude`, `fragment_kinetic_energy_release`, `photofragment_branching_ratio`。

### 17. 开放量子系统与 Lindblad 耗散主方程 (`mod_open_quantum`)
- **理论基础与物理机制**：真实量子系统不可避免地受到环境库（真空电磁场、热声子浴）的耗散与退相干影响；密度矩阵在弱耦合与玻恩-马尔可夫近似下遵循完全正定保迹 (CPTP) Lindblad 主方程。
- **详细数学表达式**：
  - Lindblad 超算符主方程：
    $$\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}, \hat{\rho}] + \sum_k \gamma_k \left( \hat{L}_k \hat{\rho} \hat{L}_k^\dagger - \frac{1}{2} \{ \hat{L}_k^\dagger \hat{L}_k, \hat{\rho} \} \right)$$
    其中弛豫算符 $\hat{L}_{i \to j} = |j\rangle\langle i|$ 对应自发辐射跃迁，退相位算符 $\hat{L}_{\text{deph}} = |i\rangle\langle i|$ 对应纯退相干。
  - 量子信息度量函数：
    - 纯度（Purity）：$\mathcal{P} = \text{Tr}(\hat{\rho}^2) \in [1/N, 1]$
    - 冯·诺依曼熵（von Neumann Entropy）：$S_{\text{vN}} = -\text{Tr}(\hat{\rho} \ln\hat{\rho})$
    - 全局 $l_1$-范数量子相干度：$\mathcal{C}_{l_1}(\hat{\rho}) = \sum_{i \ne j} |\rho_{ij}|$
- **核心 API 映射**：`lindblad_init`, `lindblad_add_decay_channel`, `lindblad_add_dephasing_channel`, `lindblad_rhs`, `propagate_lindblad_rk4`, `density_matrix_purity`, `von_neumann_entropy`, `quantum_coherence_l1`。

### 18. 量子最优控制理论 Krotov 算法 (`mod_optimal_control`)
- **理论基础与物理机制**：量子态工程要求设计形状受限的超快激光脉冲，使得量子体系从初态以最高保真度转移至预定目标态；Krotov 算法通过引入伴随协态保证每步迭代代价泛函严格单调无振荡提升。
- **详细数学表达式**：
  - 目标代价泛函（终态投影保真度与场强能耗惩罚）：
    $$J[\psi, \epsilon] = |\langle \psi(T) | \phi_{\text{target}} \rangle|^2 - \int_0^T \frac{\alpha_0}{S(t)} [\epsilon(t) - \epsilon_{\text{ref}}(t)]^2 dt$$
    其中 $S(t) = \sin^2(\pi t / T)$ 为脉冲端点包络约束函数。
  - 伴随协态反向传播方程与终态边界条件：
    $$i\hbar \frac{\partial |\chi(t)\rangle}{\partial t} = \hat{H}^\dagger |\chi(t)\rangle, \quad |\chi(T)\rangle = \langle \phi_{\text{target}} | \psi(T) \rangle |\phi_{\text{target}}\rangle$$
  - Krotov 激光电场原位更新公式（严格单调收敛 $\Delta J \ge 0$）：
    $$\epsilon^{(k+1)}(t) = \epsilon^{(k)}(t) + \frac{S(t)}{\alpha_0} \text{Im}\left[ \langle \chi^{(k)}(t) | \hat{\mu} | \psi^{(k+1)}(t) \rangle \right]$$
- **核心 API 映射**：`oct_config_init`, `state_transfer_fidelity`, `oct_shape_function`, `oct_krotov_step`, `oct_optimize_pulse`。

### 19. 非含时散射理论与超冷碰撞 (`mod_ti_scattering`)
- **理论基础与物理机制**：两体量子碰撞在渐近区遵从分波展开与光学定理；深阱区波函数通过 Johnson 矩阵对数导数法和 Manolopoulos 变步长分段扇区传播递推，彻底免疫经典禁区闭通道指数发散；在极低能区通过有效力程展开（ERE）与长程范德华解析色散提取散射长度。
- **详细数学表达式**：
  - 径向定态薛定谔方程与 Riccati 渐近边界条件：
    $$u_l''(r) + \left[ k^2 - \frac{l(l+1)}{r^2} - \frac{2\mu}{\hbar^2} V(r) \right] u_l(r) = 0, \quad u_l(r) \xrightarrow{r \to \infty} A_l \left[ \hat{j}_l(kr) \cos\delta_l - \hat{n}_l(kr) \sin\delta_l \right]$$
    其中弹性散射反应矩阵 $K_l = \tan\delta_l$，幺正散射矩阵元 $S_l = e^{2i\delta_l}$，跃迁矩阵元 $T_l = S_l - 1$。
  - 分波弹性截面、总截面与光学定理：
    $$\sigma_l = \frac{4\pi}{k^2}(2l+1)\sin^2\delta_l, \quad \sigma_{\text{tot}} = \sum_{l=0}^\infty \sigma_l = \frac{4\pi}{k} \text{Im}[f(0)]$$
  - 超低动能区 $s$-波有效力程展开 (ERE)：
    $$k \cot\delta_0(k) = -\frac{1}{a_s} + \frac{1}{2} r_0 k^2 - P_r r_0^3 k^4 + \mathcal{O}(k^6)$$
  - Gribakin-Flambaum 范德华平均散射长度与半经典散射长度：
    $$\bar{a} = \frac{2\pi}{\Gamma(1/4)^2}\left( \frac{2\mu C_6}{\hbar^2} \right)^{1/4} \approx 0.4779888 \cdot \left( \frac{2\mu C_6}{\hbar^2} \right)^{1/4}, \quad a_s = \bar{a}\left[ 1 - \tan\left( \Phi - \frac{\pi}{8} \right) \right]$$
  - Johnson 矩阵比值对数导数递推格式：
    $$\mathbf{R}_{i+1} = \mathbf{M}_i - \mathbf{R}_i^{-1}, \quad \mathbf{M}_i = 12 \mathbf{Q}_i^{-1} - 10 \mathbf{I}, \quad \mathbf{Q}_i = \mathbf{I} - \frac{h^2}{12}\mathbf{W}(r_i)$$
  - Wigner-Smith 碰撞时延与 Breit-Wigner 形状共振线型：
    $$\tau(E) = 2\hbar \frac{d\delta_l(E)}{dE} = -i\hbar S_l^\dagger(E) \frac{dS_l(E)}{dE} \approx \frac{2\hbar \Gamma}{(E - E_R)^2 + (\Gamma/2)^2}$$
- **核心 API 映射**：`riccati_bessel_neumann`, `calc_scattering_length_numerov`, `calc_scattering_length_logder`, `calc_phase_shift_single_l`, `calc_scattering_wavefunction_ti`, `calc_partial_wave_cross_sections`, `optical_theorem_cross_section`, `calc_differential_cross_section`, `fit_effective_range_expansion`, `van_der_waals_mean_length`, `gribakin_flambaum_length`, `analyze_shape_resonance`, `calc_coupled_channel_smatrix_2x2`, `calc_feshbach_resonance_scan`, `create_segmented_grid`, `calc_multichannel_close_coupling_segmented_logder`。

### 20. 含时波包散射理论与 S-矩阵 (`mod_td_scattering`)
- **理论基础与物理机制**：含时波包动力学通过单次推进高斯波包直接求解全连续能域散射信息；利用时间-能量傅里叶半变换在渐近区精确提取透射几率谱 $T(E)$，并利用 Möller 动量投影算符提取非弹性 $S$ 矩阵元与非绝热通道分支比。
- **详细数学表达式**：
  - 入射最小不确定度高斯散射波包与动量谱分布：
    $$\psi(x, 0) = (2\pi\sigma_x^2)^{-1/4} \exp\left( -\frac{(x - x_0)^2}{4\sigma_x^2} + i k_0 x \right), \quad g(k) = (2\sigma_x^2/\pi)^{1/4} \exp\left( -\sigma_x^2 (k - k_0)^2 - i k x_0 \right)$$
  - 渐近边界通量监测面 $x_{\text{det}}$ 处时间-能量傅里叶散射振幅：
    $$A(E) = \frac{1}{\sqrt{2\pi\hbar}} \int_0^\infty \psi(x_{\text{det}}, t) e^{i E t/\hbar} dt$$
  - 能量分辨连续谱透射几率与反射几率：
    $$T(E) = \frac{\hbar k_E}{\mu |g(k_E)|^2} |A_{\text{trans}}(E)|^2, \quad R(E) = \frac{\hbar k_E}{\mu |g(k_E)|^2} |A_{\text{refl}}(E)|^2, \quad T(E) + R(E) = 1.0$$
  - 动力学散射矩阵元与含时 Wigner 散射时延：
    $$S(E) = \frac{A_{\text{scatter}}(E)}{A_{\text{free}}(E)} = e^{2i\delta(E)}, \quad \tau_W(E) = 2\hbar \frac{d\delta(E)}{dE}$$
  - 全空间连续能量本征函数原位半傅里叶谱投影提取：
    $$\psi_E(x) = \frac{\hbar k_E}{\mu \sqrt{2\pi} g(k_E)} \int_0^\infty \Psi(x, t) e^{i E t/\hbar} dt \implies \hat{H}\psi_E(x) = E \psi_E(x)$$
- **核心 API 映射**：`gaussian_wavepacket_1d`, `gaussian_momentum_amplitude`, `accumulate_flux_amplitude`, `calculate_td_transmission`, `calculate_td_smatrix_element`, `project_wavepacket_to_smatrix`, `multichannel_td_smatrix_elements`, `wavepacket_centroid_position`, `extract_td_scattering_wavefunction`。

### 21. 外场电磁场超冷散射与四大基组变换 (`mod_field_scattering`)
- **理论基础与物理机制**：超冷碱金属碰撞中，外加磁场打破单原子超精细简并，驱动单重态 $V_0(R)$（自旋 $S=0$）与三重态 $V_1(R)$（自旋 $S=1$）势能面间的塞曼自旋交换；四大经典物理基组（非耦合基、单体自旋基、总自旋基与场缀饰基）之间的严格酉变换是多通道磁 Feshbach 共振计算的基石。
- **详细数学表达式**：
  - 单原子 Zeeman-超精细 Breit-Rabi 解析哈密顿量：
    $$\hat{H}_{\text{atom}} = A_{\text{hfs}} \mathbf{I} \cdot \mathbf{S} + (g_J \mu_B S_z - g_I \mu_N I_z) B$$
    $$E(F = I \pm 1/2, M) = -\frac{\Delta E_{\text{hfs}}}{2(2I+1)} - g_I \mu_N B M \pm \frac{\Delta E_{\text{hfs}}}{2} \sqrt{1 + \frac{4M x}{2I+1} + x^2}, \quad x = \frac{(g_J \mu_B + g_I \mu_N) B}{\Delta E_{\text{hfs}}}$$
  - 两体自旋交换势算符分解（投影单重态与三重态）：
    $$\hat{V}_{\text{spin}}(R) = V_0(R) \hat{\mathcal{P}}_0 + V_1(R) \hat{\mathcal{P}}_1 = \bar{V}(R) + \Delta V(R) \mathbf{S}_1 \cdot \mathbf{S}_2$$
    其中 $\bar{V}(R) = \frac{V_0(R) + 3 V_1(R)}{4}, \Delta V(R) = V_1(R) - V_0(R)$。
  - 四大基组幺正变换算符 $\mathbf{U} = \langle \text{basis}_A | \text{basis}_B \rangle$（如非耦合基到总自旋耦合基）：
    $$\langle s_1 m_{s1} i_1 m_{i1} s_2 m_{s2} i_2 m_{i2} | (s_1 s_2)S M_S (i_1 i_2)I M_I \rangle = \langle s_1 m_{s1} s_2 m_{s2} | S M_S \rangle \langle i_1 m_{i1} i_2 m_{i2} | I M_I \rangle$$
  - 磁 Feshbach 共振色散拟合公式：
    $$a_s(B) = a_{\text{bg}} \left( 1 - \frac{\Delta B}{B - B_0} \right)$$
- **核心 API 映射**：`get_cold_atom_preset`, `calc_breit_rabi_energies`, `build_field_collision_channels`, `calc_basis_transform_matrix`, `calc_zeeman_hyperfine_hamiltonian`, `calc_magnetic_feshbach_resonance_scan`。

### 22. 各向异性偶极超冷散射与自旋弛豫 (`mod_dipolar_scattering`)
- **理论基础与物理机制**：磁性原子（如 Cr, Dy, Er）或极性分子（如 KRb, NaK）具有强各向异性偶极-偶极相互作用，破坏单轨道角动量守恒，驱动 $s$ 分波与 $d$ 分波间的强偶极混合，并在磁阱中产生非弹性两体自旋弛豫加热损耗。
- **详细数学表达式**：
  - 电子磁偶极-偶极相互作用 (MDDI) 秩-2 球谐张量展开：
    $$\hat{V}_{\text{dd}}(\mathbf{r}) = \frac{\mu_0 g^2 \mu_B^2}{4\pi r^3} \left[ \mathbf{S}_1 \cdot \mathbf{S}_2 - 3(\mathbf{S}_1 \cdot \hat{r})(\mathbf{S}_2 \cdot \hat{r}) \right] = -\frac{\mu_0 g^2 \mu_B^2}{4\pi r^3} \sqrt{\frac{24\pi}{5}} \sum_{q=-2}^2 (-1)^q Y_{2,-q}(\hat{r}) [\mathbf{S}_1 \otimes \mathbf{S}_2]^{(2)}_q$$
  - 轨道球谐角动量多极矩阵元：
    $$\langle l m_l | C_{2, q} | l' m_l' \rangle = (-1)^{m_l} \sqrt{(2l+1)(2l'+1)} \begin{pmatrix} l & 2 & l' \\ 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} l & 2 & l' \\ -m_l & q & m_l' \end{pmatrix}$$
  - 极性分子外加直流电场 Stark 诱导偶极矩与特征电偶极长度：
    $$d_{\text{ind}}(\mathcal{E}) = d_0 \langle \cos\theta \rangle_{\mathcal{E}}, \quad a_d = \frac{m d_{\text{ind}}^2}{2 \hbar^2}$$
  - 超冷磁阱中两体偶极自旋弛豫截面与玻尔兹曼热平均速率系数：
    $$\sigma_{\text{rel}}(E) = \frac{8\pi}{15 k^2} \left( \frac{\mu C_{\text{dd}}}{\hbar^2} \right)^2 \frac{k_f}{k_i}, \quad K_{\text{rel}}(T) = \sqrt{\frac{8 k_B T}{\pi \mu}} \int_0^\infty \left( \frac{E}{k_B T} \right) \sigma_{\text{rel}}(E) e^{-E/(k_B T)} \frac{dE}{k_B T}$$
- **核心 API 映射**：`c2q_spherical_harmonic_tensor`, `c2q_orbital_matrix_element`, `spin_tensor_coupled_matrix_element`, `calc_mddi_coupling_strength`, `calc_dipolar_relaxation_cross_section`, `calc_dipolar_relaxation_thermal_rate`, `calc_stark_induced_dipole`, `calc_electric_dipolar_length`, `calc_dipolar_potential_matrix`。

### 23. 超冷光缔合谱学与分子生成 (`mod_photoassociation`)
- **理论基础与物理机制**：超冷原子碰撞过程中吸收红失谐激光光子，跃迁至激发态二聚体分子的长程弱束缚振动态；利用双光子 STIRAP 或反转拉曼跃迁可高效率合成绝对振转基态极性分子。
- **详细数学表达式**：
  - 能量归一化自由态 $\psi_E(R)$ 与束缚分子态 $\psi_v(R)$ 空间 Franck-Condon 重叠积分：
    $$I_{\text{FB}}(E) = \int_0^\infty \psi_{\text{free}}(E, R) \mu(R) \psi_{\text{bound}}(R) dR, \quad f_{\text{FB}}(E) = |I_{\text{FB}}(E)|^2$$
  - 激光强度驱动受激展宽线宽与单能量光缔合截面：
    $$\hbar \Gamma_{\text{stim}}(E) = 2\pi \left( \frac{I}{2\varepsilon_0 c} \right) |d_{\text{el}} I_{\text{FB}}(E)|^2, \quad \sigma_{\text{PA}}(E, \Delta) = \frac{\pi}{k^2} \frac{\hbar \Gamma_{\text{stim}}(E) \gamma_{\text{sp}}}{(E - \hbar\Delta)^2 + [(\gamma_{\text{sp}} + \hbar\Gamma_{\text{stim}}(E))/2]^2}$$
  - Bohn-Julienne 麦克斯韦-玻尔兹曼热平衡光缔合速率常数：
    $$K_{\text{PA}}(T, \Delta) = \left( \frac{2\pi\hbar^2}{\mu k_B T} \right)^{3/2} \frac{1}{h} \int_0^\infty e^{-E/(k_B T)} \sigma_{\text{PA}}(E, \Delta) \frac{2E}{\hbar} dE$$
  - 双光子 STIRAP 绝热受激跃迁有效拉比耦合频率：$\Omega_{\text{eff}} = \frac{\Omega_1 \Omega_2}{2\Delta_1}$。
- **核心 API 映射**：`calc_free_bound_fc_overlap`, `calc_free_bound_fc_density`, `calc_pa_stimulated_linewidth`, `calc_pa_cross_section`, `calc_pa_thermal_rate_coefficient`, `calc_pa_detuning_scan`, `calc_twophoton_raman_coupling`。

### 24. 超冷三体复合与 Efimov 少体物理 (`mod_three_body_recombination`)
- **理论基础与物理机制**：在大散射长度 $|a| \gg r_{\text{vdW}}$ 强相互作用极限下，三体系统在超径向展现离散标度不变性（Efimov 物理效应）；三体碰撞形成深束缚分子并释放动能，造成超冷原子捕获阱的特征三体复合原子损耗。
- **详细数学表达式**：
  - Efimov 超径向薛定谔超越代数方程（全同玻色子）：
    $$\frac{8}{\sqrt{3}} \frac{\sin(s_0 \pi / 6)}{s_0 \cos(s_0 \pi / 2)} = 1 \implies s_0 \approx 1.00624, \quad \lambda = e^{\pi / s_0} \approx 22.694$$
  - Braaten-Hammer 普适三体复合损失速率公式：
    $$K_3(a > 0) = \frac{128\pi^2 (4\pi - 3\sqrt{3})\hbar}{m} a^4 \left[ \sin^2\left( s_0 \ln\frac{a}{a_+} \right) + \sinh^2\eta_+ \right]$$
    $$K_3(a < 0) = \frac{4590 \sinh(2\eta_-)}{\sin^2\left( s_0 \ln\frac{|a|}{a_-} \right) + \sinh^2\eta_-} \frac{\hbar |a|^4}{m}$$
    在 $a > 0$ 呈现不同通道量子干涉极小值窗口 $a_+^{(n)}$，在 $a < 0$ 侧呈现 Efimov 三聚体束缚态引起的巨大共振损耗峰 $a_-^{(n)}$。
  - 强相互作用幺正饱和极限有限温度幂律：$K_3^{\text{unitary}}(T) \approx \frac{36\sqrt{3}\pi^2 \hbar^5}{m^3 (k_B T)^2}$。
- **核心 API 映射**：`solve_efimov_s0_identical_bosons`, `calc_three_body_recombination_a_positive`, `calc_three_body_recombination_a_negative`, `calc_unitary_three_body_loss_temperature`。

### 25. 低维光晶格受限散射与约束诱导共振 (`mod_confined_scattering`)
- **理论基础与物理机制**：强二维横向紧聚焦光晶格波导中，质心横向零点振荡长度 $a_\perp = \sqrt{\hbar/(\mu\omega_\perp)}$ 限制了分子的空间发散；当 3D 自由空间散射长度与横向束缚尺度匹配时，准一维有效散射耦合发生几何共振发散（Olshanii CIR）。
- **详细数学表达式**：
  - 横向简谐束缚势与零点振荡尺度：$V_\perp(\rho) = \frac{1}{2}\mu \omega_\perp^2 \rho^2, \quad a_\perp = \sqrt{\frac{\hbar}{\mu\omega_\perp}}$
  - Olshanii 约束诱导共振 (CIR) 临界发散关系式：
    $$g_{\text{1D}} = \frac{2\hbar^2 a_s}{\mu a_\perp^2} \frac{1}{1 - C \frac{a_s}{a_\perp}}, \quad C = -\frac{\zeta(1/2)}{\sqrt{2}} \approx 1.0326$$
    共振极点出现在 $a_s = a_{\text{CIR}} = a_\perp / C$；1D 有效散射长度满足 $a_{\text{1D}} = -\frac{a_\perp^2}{2 a_s}\left( 1 - C\frac{a_s}{a_\perp} \right)$。
  - 受限波导分子二聚体结合能：$E_b^{\text{1D}} = \frac{\hbar^2}{2\mu (a_{\text{1D}})^2}$
  - Lieb-Liniger 强关联费米化参数（Tonks-Girardeau 极限）：$\gamma_{\text{LL}} = \frac{m g_{\text{1D}}}{\hbar^2 n_{\text{1D}}} \gg 1$。
- **核心 API 映射**：`init_waveguide_1d`, `calc_olshanii_cir_parameters`, `calc_confined_dimer_binding_energy`, `calc_lieb_liniger_parameter`。

### 26. 自电离体系、Fano 共振与复坐标旋转 (`mod_autoionization_fano`)
- **理论基础与物理机制**：当离散准束缚态（如原子的双激发态）的能量落在连续电离谱之内时，电子组态相互作用引起两条干涉跃迁路径：离散跃迁路径与连续跃迁路径相干叠加，产生特征性非对称 Fano 吸收轮廓与抗共振零点。
- **详细数学表达式**：
  - 组态相互作用 Fano 线型公式与反共振极小：
    $$\sigma(\epsilon) = \sigma_0 \frac{(q + \epsilon)^2}{1 + \epsilon^2}, \quad \epsilon = \frac{E - E_r}{\Gamma/2}$$
    其中不对称因子为 $q = \frac{\langle \Phi | \hat{T} | i \rangle}{\pi V_E^* \langle \psi_E | \hat{T} | i \rangle}$，自电离宽度为 $\Gamma = 2\pi |V_E|^2 = 2\pi |\langle \psi_E | \hat{H} | \Phi \rangle|^2$。
  - 复坐标旋转法 (Complex Coordinate Rotation, CCR) 非厄米谱分解：
    $$r \to r e^{i\theta}, \quad \hat{H}(\theta) = e^{-2i\theta} \hat{T} + \hat{V}(r e^{i\theta}) \implies E_{\text{res}} = E_r - i \frac{\Gamma}{2}$$
    连续谱沿负虚轴旋转 $2\theta$，共振准束缚态极点暴露于复能量下半平面，衰变寿命为 $\tau = \hbar / \Gamma$。
- **核心 API 映射**：`calc_fano_profile`, `calc_autoionization_lifetime`, `solve_ccr_resonance_model`。

### 27. 交叉静电磁场转振-自旋动力学 (`mod_crossed_field_scattering`)
- **理论基础与物理机制**：同时处在外加静电场 $\mathbf{E}$ 与外加静磁场 $\mathbf{B}$ 中的极性开壳层分子，外场相互作用与内部转动-超精细耦合相互竞争；任意非共线倾角 $\beta$ 彻底破坏分子空间宇称与投影对称性，产生复杂避免交叉能谱与空间三维定向取向。
- **详细数学表达式**：
  - 交叉电磁场分子有效哈密顿量：
    $$\hat{H} = B_e \hat{\mathbf{J}}^2 + \gamma_{\text{sr}} \hat{\mathbf{J}} \cdot \hat{\mathbf{S}} - \boldsymbol{\mu}_e \cdot \mathbf{E} - \boldsymbol{\mu}_m \cdot \mathbf{B}$$
    其中静电场取沿 $z$ 轴 $\mathbf{E} = E \hat{z}$，静磁场位于 $xz$ 平面 $\mathbf{B} = B(\sin\beta \hat{x} + \cos\beta \hat{z})$。
  - 实验室系分子空间电取向度与自旋极化分量：
    $$\langle \cos\theta \rangle_n = \langle \psi_n | \cos\theta | \psi_n \rangle, \quad \langle S_z \rangle_n = \langle \psi_n | \hat{S}_z | \psi_n \rangle, \quad \langle S_x \rangle_n = \langle \psi_n | \hat{S}_x | \psi_n \rangle$$
- **核心 API 映射**：`init_crossed_field_config`, `solve_crossed_field_eigenstates`, `calc_crossed_field_observables`, `scan_tilt_angle_spectrum`。

### 28. 三原子反应散射、Jacobi 坐标与几何相位 (`mod_triatomic_geometry`)
- **理论基础与物理机制**：气相三原子反应碰撞 $A + BC \to AB + C$ 基于质心分离质心 Jacobi 反应坐标体系展开；势能面（PES）采用经典 LEPS 形式构建过渡态活化势垒；沿势能面避差交叉闭合回路环绕将诱导非平凡的 Longuet-Higgins / Berry 几何相位。
- **详细数学表达式**：
  - 质心 Jacobi 反应坐标向三原子核间距的可逆保模变换：
    $$\mathbf{r} = \mathbf{r}_B - \mathbf{r}_A, \quad \mathbf{R} = \mathbf{r}_C - \frac{m_A \mathbf{r}_A + m_B \mathbf{r}_B}{m_A + m_B}$$
    $$R_{AB} = |\mathbf{r}|, \quad R_{BC} = \left| \mathbf{R} - \frac{m_A}{m_A + m_B}\mathbf{r} \right|, \quad R_{AC} = \left| \mathbf{R} + \frac{m_B}{m_A + m_B}\mathbf{r} \right|$$
  - Sato 修正 London-Eyring-Polanyi-Sato (LEPS) 势能面解析表达：
    $$V(r_1, r_2, r_3) = \sum_{i=1}^3 \frac{Q_i}{1 + S_i} - \sqrt{\frac{1}{2} \left[ \left(\frac{J_1}{1+S_1} - \frac{J_2}{1+S_2}\right)^2 + \left(\frac{J_2}{1+S_2} - \frac{J_3}{1+S_3}\right)^2 + \left(\frac{J_3}{1+S_3} - \frac{J_1}{1+S_1}\right)^2 \right]}$$
  - 锥形交叉 (Conical Intersection, CI) 与拓扑 Berry 几何相位闭路积分：
    $$\Phi_B = \oint_C \mathbf{A}(\mathbf{R}) \cdot d\mathbf{R} = \oint_C \langle \psi_{\text{adia}}(\mathbf{R}) | \nabla_{\mathbf{R}} | \psi_{\text{adia}}(\mathbf{R}) \rangle \cdot d\mathbf{R} = \pi$$
- **核心 API 映射**：`jacobi_to_internuclear`, `internuclear_to_jacobi`, `calc_leps_potential`, `calc_conical_intersection_adiabats`, `calc_berry_phase_around_ci`。

### 29. 旋量玻色爱因斯坦凝聚自旋动力学 (`mod_spinor_bec`)
- **理论基础与物理机制**：$F=1$ 旋量玻色-爱因斯坦凝聚体（Spinor BEC）由三组分超冷原子波函数构型；具有接触自旋无关常数 $c_0$ 与自旋交换常数 $c_2$；在单模近似（SMA）下凝聚体空间波函数锁定，自旋动力学展现为保全几率与保纵向磁化强度的非线性约瑟夫森相干自旋振荡。
- **详细数学表达式**：
  - 多组分含时 Gross-Pitaevskii 方程组：
    $$i\hbar \frac{\partial \psi_m}{\partial t} = \left(-\frac{\hbar^2\nabla^2}{2M} + V_{\text{trap}}(\mathbf{r}) + q m^2 - p m\right)\psi_m + c_0 n(\mathbf{r}) \psi_m + c_2 n(\mathbf{r}) \sum_{\alpha=x,y,z} (\mathbf{F}_\alpha)_{mm'} \psi_{m'} \cdot \mathbf{F}(\mathbf{r})$$
  - 相互作用参数与低能 s 波散射长度关系：
    $$c_0 = \frac{4\pi\hbar^2}{M}\frac{a_0 + 2a_2}{3}, \quad c_2 = \frac{4\pi\hbar^2}{M}\frac{a_2 - a_0}{3}$$
  - 单模近似 (Single-Mode Approximation, SMA) 下自旋振荡方程组：
    $$i\hbar \frac{d\psi_{\pm 1}}{dt} = \left[ c_0 n + c_2 n (|\psi_{\pm 1}|^2 + |\psi_0|^2 - |\psi_{\mp 1}|^2) \pm p + q \right] \psi_{\pm 1} + c_2 n \psi_0^2 \psi_{\mp 1}^*$$
    $$i\hbar \frac{d\psi_0}{dt} = \left[ c_0 n + c_2 n (|\psi_1|^2 + |\psi_{-1}|^2) \right] \psi_0 + 2 c_2 n \psi_1 \psi_{-1} \psi_0^*$$
  - Breit-Rabi 二阶塞曼位移：
    $$q(B) = \frac{(g_I - g_J)^2 \mu_B^2 B^2}{16 \Delta E_{\text{hfs}}}$$
- **核心 API 映射**：`init_spinor_preset`, `calc_quadratic_zeeman_shift`, `propagate_spinor_sma_rk4`。

### 30. 三原子超球面反应动力学与热速率常数 (`mod_hyperspherical_reactive`)
- **理论基础与物理机制**：气相三体反应体系在 Delves 质量标度超球面坐标下将三体散射解耦为超半径演化与角向超角运动；基于过渡态鞍点解析 Eckart 势垒精确刻画量子隧穿与反射效应；通过累积反应几率 (CRP) 玻尔兹曼热积分输出微观可逆热反应速率常数。
- **详细数学表达式**：
  - Delves 质量标度因子与超角偏转角：
    $$d = \left( \frac{m_A m_C}{m_{AB} m_{ABC}} \right)^{1/4}, \quad \beta_{\text{skew}} = \arctan\left( \sqrt{\frac{m_B(m_A+m_B+m_C)}{m_A m_C}} \right)$$
  - 不对称 Eckart 势垒透射几率严格解析解：
    $$V_{\text{Eckart}}(x) = \frac{A y}{1-y} + \frac{B y}{(1-y)^2}, \quad y = -e^{\alpha x}$$
    $$P(E) = \frac{\cosh[2\pi(k_1 + k_2)] - \cosh[2\pi(k_1 - k_2)]}{\cosh[2\pi(k_1 + k_2)] + \cosh[2\pi d]}$$
    $$k_1 = \frac{\sqrt{2\mu E}}{\hbar}, \quad k_2 = \frac{\sqrt{2\mu(E - V_0 + V_1)}}{\hbar}, \quad d = \frac{1}{2}\sqrt{\frac{8\mu V_0}{\alpha^2\hbar^2} - 1}$$
  - 全量子累积反应几率 (CRP) 与正则热速率常数玻尔兹曼积分：
    $$N(E) = \sum_{v, J} P_{v, J}(E), \quad k(T) = \frac{1}{2\pi\hbar Q_R(T)} \int_0^\infty N(E) e^{-E / (k_B T)} dE$$
  - Wigner 势垒量子穿透修正：
    $$\kappa_{\text{Wigner}}(T) = 1 + \frac{1}{24}\left( \frac{\hbar \omega^{\ddagger}}{k_B T} \right)^2$$
- **核心 API 映射**：`init_reaction_mass`, `calc_eckart_transmission`, `calc_cumulative_reaction_probability`, `calc_canonical_rate_constant`, `calc_tst_wigner_rate`。

### 31. 超冷偶极量子液滴与李-黄-杨量子涨落 (`mod_dipolar_droplets_lhy`)
- **理论基础与物理机制**：各向异性长程偶极-偶极相互作用导致玻色凝聚体在平均场平均引力下坍塌；通过引入 Lee-Huang-Yang (LHY) 零点量子涨落超越平均场排斥项，在自由空间形成零压平衡、内部密度平顶的自束缚超冷量子液滴。
- **详细数学表达式**：
  - 扩展 Gross-Pitaevskii 方程 (eGPE)：
    $$i\hbar \frac{\partial \psi}{\partial t} = \left[ -\frac{\hbar^2\nabla^2}{2M} + V_{\text{ext}}(\mathbf{r}) + g |\psi|^2 + \Phi_{\text{dd}}(\mathbf{r}) + \gamma_{\text{LHY}} |\psi|^3 \right] \psi$$
  - 磁偶极特征长度、相对偶极强度与 Pelster-Lima $Q_5$ 涨落积分：
    $$a_{\text{dd}} = \frac{\mu_0 \mu_{\text{mag}}^2 M}{12\pi\hbar^2}, \quad \epsilon_{\text{dd}} = \frac{a_{\text{dd}}}{a_s}, \quad \gamma_{\text{LHY}} = \frac{128\sqrt{\pi}\hbar^2 a_s^{5/2}}{3M} Q_5(\epsilon_{\text{dd}})$$
    $$Q_5(\epsilon) = \frac{1}{2}\int_0^1 dx \, (1 - \epsilon + 3\epsilon x^2)^{5/2}$$
  - 自由空间零压平顶平衡核心密度与自束缚负化学势判据：
    $$n_0 = \frac{25\pi}{16384} \frac{(1 - \epsilon_{\text{dd}})^2}{a_s^5 Q_5(\epsilon_{\text{dd}})^2}, \quad \mu(n_0) = g n_0 \left(1 - \frac{4}{3}\epsilon_{\text{dd}}\right) + \frac{5}{2}\gamma_{\text{LHY}} n_0^{3/2} < 0$$
- **核心 API 映射**：`init_dipolar_droplet_param`, `calc_pelster_lima_q5`, `calc_equilibrium_droplet_density`, `calc_droplet_chemical_potential`, `calc_critical_atom_number`。

### 32. 强场非顺序双电离与电子重碰撞动量谱 (`mod_strong_field_nsdi`)
- **理论基础与物理机制**：强激光场原子电离电子在时变交变场中发生反向加速运动并回碰母离子，通过经典重碰撞散射激发或碰撞电离第二电子；再电离电子动量在 COLTRIMS 符合测量谱上展现出特征性的同向平行关联分布与非顺序电离产率“膝盖平台结构”。
- **详细数学表达式**：
  - 经典电子动力学轨道与二次回碰相位根：
    $$v(t) = \frac{e F_0}{m \omega}(\sin\omega t - \sin\phi_0), \quad x(t) = \frac{e F_0}{m \omega^2}[\cos\phi_0 - \cos\omega t - (\omega t - \phi_0)\sin\phi_0]$$
    $$x(\phi_r) = 0 \implies E_{\text{rec}}(\phi_r) = \frac{1}{2} m v^2(\phi_r) \le 3.173 U_p, \quad U_p = \frac{e^2 F_0^2}{4 m \omega^2}$$
  - 二次激发/电离 Lotz 碰撞截面：
    $$\sigma_{\text{Lotz}}(E) = \sum_i a_i q_i \frac{\ln(E/I_i)}{E \cdot I_i} \left[ 1 - b_i \exp\left( -c_i \left(\frac{E}{I_i} - 1\right) \right) \right]$$
  - 双电子纵向动量关联函数与皮尔逊关联系数：
    $$P(p_{z1}, p_{z2}) = \iint W_{\text{ADK}}(\phi_0) \sigma_{\text{rec}}(E_{\text{rec}}) \delta(p_{z1} + p_{z2} - P_z) d\phi_0, \quad C_{\text{corr}} = \frac{\langle p_{z1} p_{z2} \rangle}{\sqrt{\langle p_{z1}^2 \rangle \langle p_{z2}^2 \rangle}} > 0$$
- **核心 API 映射**：`init_nsdi_laser`, `calc_recollision_trajectory`, `calc_lotz_cross_section`, `calc_nsdi_2d_momentum_dist`, `calc_double_ion_yield_curve`。

### 33. 磁与光 Feshbach 共振与分子弱束缚态 (`mod_feshbach_bound_states`)
- **理论基础与物理机制**：超冷两体碰撞散射中，外磁场或外光场驱动闭通道束缚分子能级扫描至与开通道散射渐近能量简并；开-闭通道多重耦合实现微观散射长度 $a$ 从 $-\infty$ 到 $+\infty$ 的任意调谐，并在共振点近旁生成弱束缚 Feshbach 缔合二聚体分子。
- **详细数学表达式**：
  - 磁 Feshbach 共振有效散射长度色散公式：
    $$a(B) = a_{\text{bg}} \left( 1 - \frac{\Delta B}{B - B_0} \right)$$
  - 包含有限相互作用程 $R^*$ 的双通道弱束缚态结合能与闭通道成分比率：
    $$\sqrt{\frac{2\mu |E_b|}{\hbar^2}} = \frac{-1 + \sqrt{1 + 4 R^* / a(B)}}{2 R^*}, \quad R^* = \frac{\hbar^2}{2\mu a_{\text{bg}} \delta\mu \Delta B}$$
    $$Z(B) = 1 - \frac{1}{\sqrt{1 + 2 R^* / a(B)}} = \frac{1}{\delta\mu} \frac{\partial E_b}{\partial B}$$
  - 光 Feshbach 共振 (OFR) 复散射长度与光致非弹性两体损失速率：
    $$\tilde{a}(\Delta_L) = a_{\text{bg}} + \frac{l_{\text{opt}} \Gamma_{\text{mol}} / 2}{\Delta_L + i \Gamma_{\text{mol}} / 2}, \quad K_2(\Delta_L) = \frac{4\pi\hbar}{\mu} \text{Im}[\tilde{a}(\Delta_L)] = \frac{2\pi\hbar}{\mu} \frac{l_{\text{opt}} \Gamma_{\text{mol}}^2}{\Delta_L^2 + (\Gamma_{\text{mol}}/2)^2}$$
- **核心 API 映射**：`init_mfr_preset`, `calc_mfr_scattering_length`, `calc_mfr_bound_energy_coupled`, `calc_mfr_closed_channel_fraction`, `calc_ofr_complex_scattering_length`, `calc_ofr_inelastic_loss_rate`。

### 34. 阿秒瞬态吸收光谱与光诱导态自电离干涉 (`mod_attosecond_transient_absorption`)
- **理论基础与物理机制**：超快孤立极紫外 (XUV) 阿秒脉冲激发原子内壳层双激发自电离态，强红外 (NIR) 激光控制场施加动态 AC Stark 调制与光诱导态 (LIS) 耦合；在时间-能量二维瞬态吸收谱上展现出动态 Fano 线型演化与超快量子拍频现象。
- **详细数学表达式**：
  - 相位微扰模型 (PPM) 下含时偶极相位调制与动态 Fano 不对称参数：
    $$d(t) \propto e^{-i E_0 t/\hbar - \Gamma t / (2\hbar)} e^{i \Delta\phi(t, \tau)}, \quad \Delta\phi(t, \tau) = -\frac{1}{\hbar} \int_\tau^t \Delta E_{\text{AC}}(t') dt'$$
    $$q(\tau) = \frac{q_0 + \tan[\Delta\phi(\tau)]}{1 - q_0 \tan[\Delta\phi(\tau)]}$$
  - 光诱导态 (LIS) 能量准能级与明暗态量子拍频周期：
    $$E_{\text{LIS}} = E_{\text{dark}} \pm \hbar\omega_{\text{NIR}} + \alpha_{\text{Stark}} I_{\text{NIR}}, \quad T_{\text{beat}} = \frac{h}{|E_{\text{bright}} - E_{\text{LIS}}|}$$
  - 阿秒瞬态吸收差分光密度二维矩阵：
    $$\Delta\text{OD}(\omega, \tau) = -\log_{10}\left( \frac{I_{\text{trans}}(\omega, \tau)}{I_0(\omega)} \right) \propto -\text{Im}\left[ \frac{\tilde{d}(\omega, \tau)}{\tilde{E}_{\text{XUV}}(\omega)} \right]$$
- **核心 API 映射**：`init_atas_helium_benchmark`, `calc_laser_dressed_fano_q`, `calc_light_induced_state_energy`, `calc_quantum_beat_period_fs`, `calc_atas_spectrum`。

### 35. 双色反向旋转圆偏振场与分子光电子圆二色性 (`mod_bicircular_pecd`)
- **理论基础与物理机制**：由频率成有理比的反向旋转圆偏振光场复合而成的双色场具有精确的离散空间-时间对称性；当手性四面体分子在手性光场或圆偏振光下电离时，电离光电子角分布 (PAD) 展现出特征的激光传播方向前后不对称性 (PECD)。
- **详细数学表达式**：
  - 双色椭圆/圆偏振场电场矢量合成：
    $$\mathbf{E}(t) = \frac{F_1}{\sqrt{2}} \left[ \cos(\omega_1 t) \hat{\mathbf{x}} + \sigma_1 \sin(\omega_1 t) \hat{\mathbf{y}} \right] + \frac{F_2}{\sqrt{2}} \left[ \cos(\omega_2 t + \phi) \hat{\mathbf{x}} + \sigma_2 \sin(\omega_2 t + \phi) \hat{\mathbf{y}} \right]$$
  - 动力学旋转对称度折叠数（如反向旋转 $\omega+2\omega$ 呈 $C_3$ 对称）：
    $$C_N: \quad N = p + q \quad (\text{当 } \omega_1 : \omega_2 = p : q, \; \sigma_1 \sigma_2 = -1)$$
  - 四面体几何手性不变量与 Ritchie 光电子前后不对称参数：
    $$\chi_{\text{mol}} = (\mathbf{r}_1 - \mathbf{r}_4) \cdot [(\mathbf{r}_2 - \mathbf{r}_4) \times (\mathbf{r}_3 - \mathbf{r}_4)] \prod_{i < j} (Z_i - Z_j)$$
    $$I(\theta, \phi) = \frac{\sigma_{\text{tot}}}{4\pi} \left[ 1 + \beta_1 P_1(\cos\theta) + \beta_2 P_2(\cos\theta) + \cdots \right], \quad G_{\text{PECD}} = \frac{\beta_1}{2}$$
- **核心 API 映射**：`init_bicircular_field`, `calc_dynamical_symmetry_fold`, `init_chiral_tetrahedral_molecule`, `calc_chirality_measure`, `calc_chiral_beta1_model`, `calc_forward_backward_asymmetry`, `calc_pecd_pad_spectrum`。

### 36. 超冷极性分子反应动力学与微波/静电偶极遮蔽 (`mod_ultracold_reaction_shielding`)
- **理论基础与物理机制**：超冷双原子极性分子由于长程各向异性电偶极吸引易发生非弹性碰撞猝灭与化学反应损失；通过施加微波蓝失谐缀饰场诱导分子间长程免交叉有效排斥势垒，将碰撞分子有效屏蔽在短程反应区外以保护蒸发冷却。
- **详细数学表达式**：
  - 微波蓝失谐缀饰态有效排斥屏蔽势：
    $$V_{\text{eff}}(R) = \frac{\hbar\Delta}{2} + \sqrt{\left(\frac{\hbar\Delta}{2}\right)^2 + \left(\frac{d_{\text{mol}}^2}{4\pi\epsilon_0 R^3}\right)^2} - \frac{C_6}{R^6}$$
  - 短程反应区 WKB 半经典量子隧穿几率：
    $$P_{\text{WKB}}(E) = \exp\left( -2 \int_{R_{\text{in}}}^{R_{\text{out}}} \sqrt{\frac{2\mu}{\hbar^2}\max(0, V_{\text{eff}}(R) - E)} \, dR \right)$$
  - 超冷两体弹性与非弹性损失速率比判据：
    $$\sigma_{\text{el}}(E) = \frac{4\pi}{k^2} \sin^2 \delta_0(E), \quad K_2^{(\text{inel})} = \frac{2h}{\mu} \langle P_{\text{WKB}}(E) \rangle_T, \quad \gamma = \frac{K_2^{(\text{el})}}{K_2^{(\text{inel})}} > 100$$
- **核心 API 映射**：`init_ultracold_molecule_preset`, `calc_effective_shielding_potential`, `calc_shielding_barrier_height`, `calc_wkb_tunneling_probability`, `calc_shielded_scattering_rates`。

### 37. 里德堡原子阻塞、PXP 约束模型与量子多体疤痕 (`mod_rydberg_blockade`)
- **理论基础与物理机制**：高度激发的里德堡原子拥有巨大的长程范德华相互作用，使得阻塞半径内的双激发被强能级失谐阻断；该物理机制将希尔伯特空间投影到无临近激发的 PXP 约束子空间中，初态 Néel 态展现出长寿命相干复苏的量子多体疤痕现象。
- **详细数学表达式**：
  - 范德华相互作用与里德堡阻塞半径：
    $$V_{\text{vdW}}(R) = \frac{C_6}{R^6}, \quad C_6 \propto n^{11}, \quad R_b = \left( \frac{|C_6|}{\hbar \Omega_{\text{Rabi}}} \right)^{1/6}$$
  - 拓扑约束 PXP 哈密顿量：
    $$\hat{H}_{\text{PXP}} = \frac{\hbar\Omega}{2} \sum_{i=1}^L \hat{P}_{i-1} \hat{\sigma}_x^{(i)} \hat{P}_{i+1} - \hbar\Delta \sum_{i=1}^L \hat{n}_i, \quad \hat{P}_i = |g_i\rangle\langle g_i| = 1 - \hat{n}_i$$
  - 反铁磁 Néel 序参量与集体相干拉比振荡态演化：
    $$\mathcal{O}_{\mathbb{Z}_2}(t) = \frac{2}{L}\sum_{i=1}^L (-1)^i \langle \psi(t) | \hat{n}_i | \psi(t) \rangle, \quad |\psi_{2\text{-atom}}(t)\rangle = \cos\left(\frac{\sqrt{2}\Omega t}{2}\right)|gg\rangle - i \sin\left(\frac{\sqrt{2}\Omega t}{2}\right)\frac{|gr\rangle+|rg\rangle}{\sqrt{2}}$$
- **核心 API 映射**：`init_rydberg_atom`, `calc_rydberg_blockade_radius`, `calc_two_atom_dynamics`, `calc_z2_order_parameter`, `calc_rydberg_scar_dynamics`。

### 38. 表面量子散射与选择性吸附共振 (`mod_surface_scattering`)
- **理论基础与物理机制**：热能原子/分子束在晶体表面散射时受周期性点阵势作用发生量子布拉格衍射；在硬波纹表面 (HCS) 模型下衍射强度由第一类贝塞尔函数支配；当入射粒子动能与表面 Morse 束缚态能级发生微观共振耦合时诱发选择性吸附共振 (SAR) 并产生特征 Fano 线型调制。
- **详细数学表达式**：
  - 二维硬波纹表面几何与 Bragg 动量守恒：
    $$\zeta(\mathbf{R}) = \zeta_x \cos\left( \frac{2\pi x}{a_x} \right) + \zeta_y \cos\left( \frac{2\pi y}{a_y} \right), \quad \mathbf{k}_{\parallel, \mathbf{G}} = \mathbf{k}_{\parallel} + \mathbf{G} = \mathbf{k}_{\parallel} + m\mathbf{b}_x + n\mathbf{b}_y$$
  - 程函近似衍射散射 $S$ 矩阵元：
    $$S_{mn} = \frac{1}{a_x a_y} \int_0^{a_x} dx \int_0^{a_y} dy \, \exp\left[ -i \mathbf{G}\cdot\mathbf{R} - i (k_{z, \mathbf{G}} + k_{iz})\zeta(\mathbf{R}) \right] = (-i)^{|m|+|n|} J_m(c_x) J_n(c_y)$$
  - 选择性吸附共振 (SAR) 束缚能量匹配条件与声子热衰减 Debye-Waller 因子：
    $$k_{z, \mathbf{G}}^2 = \frac{2M}{\hbar^2} (E_{\text{inc}} - V_0) - |\mathbf{k}_{\parallel} + \mathbf{G}|^2 = \frac{2M}{\hbar^2} E_v^{\text{Morse}} < 0$$
    $$I_{\mathbf{G}}(T) = I_{\mathbf{G}}(0) \exp\left[ -2 W_{\mathbf{G}}(T) \right] = I_{\mathbf{G}}(0) \exp\left[ -\frac{3\hbar^2 (k_{iz} + k_{z, \mathbf{G}})^2 T}{M k_B \Theta_D^2} \right]$$
- **核心 API 映射**：`init_surface_lattice`, `init_surface_potential_morse`, `calc_surface_diffraction_channels`, `calc_hcs_diffraction_probabilities`, `calc_selective_adsorption_resonance`, `calc_surface_debye_waller`。

### 39. 气-固表面催化与 Eley-Rideal 反应动力学 (`mod_surface_reaction_er`)
- **理论基础与物理机制**：气相入射原子直接与吸附在固体表面上的化学吸附原子发生瞬态单次碰撞并结合脱附生成气相分子（直接 Eley-Rideal 反应通道）；巨大放热量 $\Delta E_{\text{exo}}$ 在飞秒至皮秒尺度内非统计分配到产物各自由度，驱动新生分子展现出极端的振动态布居反转与超热平动动能分布。
- **详细数学表达式**：
  - 二维反应势能面与总可用能量守恒：
    $$V(r, Z_{\text{cm}}) = V_{\text{gas}}(r) + V_{\text{chem}}(Z_{\text{ads}}) + V_{\text{int}}(r, Z_{\text{cm}})$$
    $$E_{\text{avail}} = E_{\text{inc}} + E_{\text{bind}} + \Delta E_{\text{exo}} = \langle E_{\text{vib}} \rangle + \langle E_{\text{rot}} \rangle + \langle E_{\text{trans}} \rangle + \Delta E_{\text{bath}}$$
  - 产物分态振动布居反转高斯分布模型：
    $$P(v) = \frac{1}{\sqrt{2\pi \sigma_v^2}} \exp\left[ -\frac{(v - v_{\text{peak}})^2}{2\sigma_v^2} \right], \quad v_{\text{peak}} \approx \frac{\alpha_{\text{vib}} E_{\text{avail}}}{\hbar \omega_e}$$
  - 入射动能依赖反应截面与微观热速率常数：
    $$\sigma_{\text{ER}}(E_i) = \sigma_0 \left(1 - \frac{V_{\text{act}}}{E_i}\right) \Theta(E_i - V_{\text{act}}), \quad k_{\text{ER}}(T) = \sqrt{\frac{8 k_B T}{\pi \mu}} \int_{V_{\text{act}}}^\infty \sigma_{\text{ER}}(E) \frac{E}{(k_B T)^2} e^{-E / (k_B T)} dE$$
- **核心 API 映射**：`init_er_reaction_system`, `calc_er_potential_2d`, `calc_er_energy_partitioning`, `calc_er_vibrational_populations`, `calc_er_reaction_cross_section`, `calc_er_thermal_rate_constant`。

### 40. 金属表面非绝热动力学与电子摩擦耗散 (`mod_surface_electronic_friction`)
- **理论基础与物理机制**：分子在金属表面散射或化学吸附过程中，核运动诱发费米能级附近的低能电子-空穴对激发（e-h pairs），破坏 Born-Oppenheimer 绝热假定；在局域密度摩擦近似（LDFA）与广义朗之万方程（GLE）下，体系展现为黏滞电子摩擦阻尼力与高斯白噪声热涨落，驱动核动能耗散与吸附键特征振动寿命衰减。
- **详细数学表达式**：
  - 广义朗之万动力学方程 (Generalized Langevin Equation, GLE)：
    $$M \ddot{z}(t) = -\frac{\partial V(z)}{\partial z} - M \int_0^t \gamma(t - t') \dot{z}(t') dt' + \xi(t), \quad \langle \xi(t) \xi(t') \rangle = 2 M \eta(z) k_B T \delta(t - t')$$
  - 局域密度摩擦近似 (LDFA) 空间依赖摩擦系数：
    $$\eta(z) = \eta_0 \exp\left[ -\beta (z - z_{\text{surf}}) \right] = \frac{4\pi}{3} k_F n_0 \sum_l (2l + 1) \sin^2(\delta_l - \delta_{l+1})$$
  - 非绝热散射电子-空穴对累积能量损失与高频振动弛豫寿命：
    $$\Delta E_{\text{loss}} = \int_0^{t_{\text{final}}} M \eta(z(t)) \dot{z}^2(t) dt, \quad \tau_{\text{vib}} = \frac{1}{\eta(z_{\text{eq}})}$$
- **核心 API 映射**：`init_metal_surface`, `calc_electronic_friction_coeff`, `calc_surface_morse_force`, `integrate_gle_scattering_trajectory`, `calc_vibrational_relaxation_rate`。

### 41. 掠入射快原子表面量子衍射与彩虹散射 (`mod_grazing_fast_atom_diffraction`)
- **理论基础与物理机制**：能量达数 keV 的快轻原子（He, Ne, H）以极小掠角 $\theta \ll 1^\circ$ 入射至平整单晶表面；沿晶轴方向的快运动与表面沟道势相互作用解耦为经典运动，而垂直晶轴的横向慢运动（$E_\perp \sim \text{meV}\sim\text{eV}$）发生高相干量子 Bragg 衍射，并在边缘展现出经典彩虹折射极大，可实现亚皮米级表面波纹幅度的高精度反演。
- **详细数学表达式**：
  - 快慢运动自由度解耦与有效垂直德布罗意波长：
    $$E_\perp = E_{\text{beam}} \sin^2\theta, \quad \lambda_\perp = \frac{h}{\sqrt{2 M E_\perp}} = \frac{h}{\sqrt{2 M E_{\text{beam}}} \sin\theta}$$
  - 表面经典彩虹散射角与表面几何极值斜率对应：
    $$\theta_R = 2 \arctan\left( \max_{x} \left| \frac{\partial \zeta(x)}{\partial x} \right| \right) \approx \frac{4\pi \zeta}{a_x}$$
  - 一维横向 Bragg 衍射峰位与亚皮米波纹幅度逆向反演：
    $$\sin\theta_m - \sin\theta_{\text{in}} = m \frac{\lambda_\perp}{a_x}, \quad \zeta = \frac{a_x \theta_R}{4\pi}$$
- **核心 API 映射**：`init_gifad_experiment`, `calc_gifad_transverse_kinematics`, `calc_gifad_rainbow_angle`, `calc_gifad_diffraction_spectrum`, `calc_surface_corrugation_from_rainbow`。

### 42. 冷离子-中性原子杂化散射与极化阱动力学 (`mod_ion_atom_scattering`)
- **理论基础与物理机制**：单离子与超冷中性原子在杂化阱碰撞中由长程诱导偶极极化势 $V(r) = -C_4 / (2r^4)$ 决定相互作用；高能区呈现经典无势垒螺旋俘获的 Langevin 动力学，超冷能区呈现修正有效力程展开（MERE）的量子多波散射；射频 Paul 阱微运动碰撞引发非平衡致热效应。
- **详细数学表达式**：
  - 极化相互作用长程特征尺度与特征能量：
    $$V(r) = -\frac{C_4}{2 r^4} = -\frac{q^2 \alpha_{\text{pol}}}{8\pi\epsilon_0 r^4}, \quad R^* = \sqrt{\frac{2\mu C_4}{\hbar^2}}, \quad E^* = \frac{\hbar^2}{2\mu (R^*)^2}$$
  - 经典 Langevin 螺旋俘获临界碰撞参数与速率系数：
    $$b_c(E) = \left( \frac{2 C_4}{E} \right)^{1/4}, \quad \sigma_L(E) = \pi b_c^2 = \pi \sqrt{\frac{2 C_4}{E}}, \quad K_L = v \sigma_L(E) = 2\pi \sqrt{\frac{C_4}{\mu}}$$
  - 极化势修正有效力程展开 (Modified Effective Range Expansion, MERE)：
    $$k \cot \delta_0 = -\frac{1}{a_s} + \frac{\pi}{3 R^*} k + \frac{4}{3} \frac{k^2}{R^*} \ln\left( \frac{k R^*}{4} \right) + \frac{1}{2} r_{\text{eff}} k^2 + \mathcal{O}(k^3)$$
  - 射频微运动诱导碰撞致热率与平衡极限温度：
    $$\frac{d\langle E_{\text{ion}} \rangle}{dt} = \kappa_{\text{rf}} \cdot q_{\text{Mathieu}}^2 \cdot K_L n_{\text{atom}} (E_{\text{ion}} - E_{\text{atom}}), \quad T_{\text{limit}} \propto T_{\text{atom}} \left(\frac{m_{\text{ion}}}{m_{\text{atom}}}\right)^\nu$$
- **核心 API 映射**：`init_ion_atom_system`, `calc_langevin_cross_section`, `calc_langevin_rate_coefficient`, `calc_ion_atom_phase_shift`, `calc_mere_phase_shift_s_wave`, `calc_rf_micromotion_heating`。

### 43. 最少开关表面跳跃与非绝热混合量子-经典动力学 (`mod_surface_hopping_fssh`)
- **理论基础与物理机制**：Tully 最少开关表面跳跃（FSSH）是描述非绝热多势能面分子动力学的经典-量子混合框架；原子核自由度沿单一绝热势能面作牛顿力学运动，电子态相干波函数沿含时薛定谔方程推进；在非绝热导数耦合矢量 (NACV) 显著区域，核轨迹以概率瞬时发生随机跳跃，并通过沿耦合矢量重标度核动量以保证总能量守恒。
- **详细数学表达式**：
  - 电子含时密度矩阵运动方程与非绝热导数耦合矢量 (NACV)：
    $$i\hbar \dot{\rho}_{jk} = (V_j - V_k)\rho_{jk} - i\hbar \sum_l \left( \dot{\mathbf{R}} \cdot \mathbf{d}_{jl} \rho_{lk} - \rho_{jl} \dot{\mathbf{R}} \cdot \mathbf{d}_{lk} \right)$$
    $$\mathbf{d}_{jk}(\mathbf{R}) = \frac{\langle \psi_j | \nabla_{\mathbf{R}} \hat{H}_{\text{el}} | \psi_k \rangle}{V_k(\mathbf{R}) - V_j(\mathbf{R})}$$
  - Tully 最少开关跳跃转移几率 (Fewest Switches Probability)：
    $$g_{k \to j} = \max\left( 0, \; \frac{2 \Delta t \, \text{Re}\left( \rho_{jk}^* \dot{\mathbf{R}} \cdot \mathbf{d}_{jk} \right)}{\rho_{kk}} \right)$$
  - 沿 NACV 矢量瞬时动量重标度与禁阻跳跃能量守恒修正：
    $$\mathbf{P}_{\text{new}} = \mathbf{P}_{\text{old}} - \gamma \mathbf{d}_{jk}, \quad \frac{(\mathbf{P}_{\text{new}})^2}{2M} + V_j = \frac{(\mathbf{P}_{\text{old}})^2}{2M} + V_k$$
    $$\gamma = \frac{\mathbf{P} \cdot \mathbf{d}_{jk}}{M} - \text{sgn}(\mathbf{P} \cdot \mathbf{d}_{jk}) \sqrt{\left(\frac{\mathbf{P} \cdot \mathbf{d}_{jk}}{M}\right)^2 - \frac{2 (V_j - V_k)}{M |\mathbf{d}_{jk}|^2}}$$
- **核心 API 映射**：`init_tully_model`, `calc_adiabatic_surface_and_nacv`, `init_fssh_trajectory`, `propagate_fssh_step`, `run_fssh_ensemble`, `propagate_ehrenfest_step`。

### 44. 强场分子定向、取向与超转子动力学 (`mod_molecular_alignment`)
- **理论基础与物理机制**：非共振强飞秒激光脉冲通过诱导极化率各向异性施加角向拉曼受激扭矩，激发宽带转动波包，在激光脉冲熄灭后在真空演化中展现出周期性无场宏观空间定向与取向复苏；利用光学离心机恒定角加速度光场可将分子连续加速至极端高角动量量子数超转子态（$J \gg 1$），最终诱发离心解离破键。
- **详细数学表达式**：
  - 极化各向异性激光相互作用势与定向序参量：
    $$V_{\text{laser}}(\theta, t) = -\frac{1}{4} \mathcal{E}^2(t) \left[ (\alpha_\parallel - \alpha_\perp) \cos^2\theta + \alpha_\perp \right] - \mu_0 \mathcal{E}(t) \cos\theta$$
    $$\langle \cos^2\theta \rangle(t) = \sum_{J, M} \rho_J \left| \sum_{J'} c_{J'}^{(J)}(t) \langle Y_{J' M} | \cos^2\theta | Y_{J M} \rangle \right|^2$$
  - 刚体分子无场相干复苏周期：
    $$T_{\text{rev}} = \frac{1}{2 B_{\text{rot}} c} = \frac{\pi \hbar}{B_{\text{rot}}}$$
  - 光学离心机恒定角加速度强迫激发与超转子离心有效势：
    $$\omega_{\text{rot}}(t) = 2 \beta t, \quad J_{\text{super}} \approx \frac{\beta \tau_{\text{pulse}}}{B_{\text{rot}}}, \quad V_{\text{eff}}(R, J) = V_{\text{Morse}}(R) + \frac{\hbar^2 J(J+1)}{2 \mu R^2}$$
- **核心 API 映射**：`init_rotor_molecule`, `calc_cos2_matrix_elements`, `calc_cos_matrix_elements`, `simulate_laser_induced_alignment`, `calc_optical_centrifuge_kick`, `calc_superrotor_dissociation`。

### 45. 超冷光晶格与玻色-哈伯德微观映射 (`mod_optical_lattice_hubbard`)
- **理论基础与物理机制**：反向对射相干激光形成周期性光学点阵驻波场；超冷原子在周期势中形成能带结构，在深阱紧束缚极限下投影至正交 Wannier 轨道，微观映射为玻色-哈伯德（Bose-Hubbard）模型；在外加恒定力作用下展示动量空间的布洛赫振荡，并在高能带边界发生 Landau-Zener 带间跃迁。
- **详细数学表达式**：
  - 光晶格驻波势与单粒子能带 Mathieu 方程：
    $$V(x) = V_0 \sin^2(k_L x), \quad \left[ -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V_0 \sin^2(k_L x) \right] \phi_{n, q}(x) = E_n(q) \phi_{n, q}(x)$$
  - 玻色-哈伯德紧束缚跃迁常数 $J$ 与在位排斥能 $U$：
    $$J = -\int dx \, w^*(x - x_i) \left[ -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x) \right] w(x - x_{i+1}) \approx \frac{4}{\sqrt{\pi}} E_R \left(\frac{V_0}{E_R}\right)^{3/4} \exp\left( -2\sqrt{\frac{V_0}{E_R}} \right)$$
    $$U = \frac{4\pi\hbar^2 a_s}{m} \int dx \, |w(x)|^4 \approx \sqrt{\frac{8}{\pi}} k_L a_s E_R \left(\frac{V_0}{E_R}\right)^{3/4}$$
  - 布洛赫振荡周期与第一激发带 Landau-Zener 隧穿几率：
    $$T_B = \frac{2\hbar k_L}{F_{\text{ext}}}, \quad P_{\text{LZ}} = \exp\left( -\frac{\pi \Delta_{\text{gap}}^2}{4 \hbar v_F F_{\text{ext}}} \right)$$
- **核心 API 映射**：`init_optical_lattice`, `calc_bloch_band_energies`, `calc_bose_hubbard_parameters`, `calc_bloch_oscillation_dynamics`。

### 46. 多原子反应路径哈密顿量与变分过渡态理论 (`mod_reaction_path_hamiltonian`)
- **理论基础与物理机制**：沿质量加权 Fukui 内禀反应坐标 (IRC) 将多自由度反应体系严格投影为一维大振幅最小能量路径 (MEP) 与 $3N-7$ 个正交振动简正模；正则变分过渡态理论 (CVT) 通过极小化沿路径各分界面的广义吉布斯自由能瓶颈，并耦合解析不对称 Eckart 势垒穿透因子，实现反应速率的高精度第一性原理预测。
- **详细数学表达式**：
  - Miller-Handy-Adams 反应路径哈密顿量 (RPH)：
    $$H_{\text{RPH}}(s, p_s, \{\mathbf{Q}, \mathbf{P}\}) = \frac{\left( p_s - \sum_{k, l} Q_k P_l B_{k, l}(s) \right)^2}{2 \left[ 1 + \sum_k Q_k B_{k, s}(s) \right]^2} + V_0(s) + \sum_{k=1}^{3N-7} \left( \frac{1}{2} P_k^2 + \frac{1}{2} \omega_k^2(s) Q_k^2 \right)$$
  - 路径曲率耦合张量与切线导数：
    $$B_{k, s}(s) = -\mathbf{L}_k^T(s) \frac{d\mathbf{t}(s)}{ds}$$
  - 正则变分过渡态 (CVT) 自由能瓶颈优化与 Eckart 隧穿校正速率：
    $$k^{\text{CVT}}(T) = \min_{s} k^{\text{GTST}}(T, s) = \min_s \left\{ \frac{k_B T}{h} \frac{Q^{\ddagger}(T, s)}{Q^R(T)} \exp\left[ -\frac{V_0(s)}{k_B T} \right] \right\}$$
    $$k^{\text{CVT/Eckart}}(T) = \kappa(T) \cdot k^{\text{CVT}}(T), \quad \kappa(T) = \frac{1}{k_B T} \int_0^\infty P_{\text{Eckart}}(E) e^{-E / (k_B T)} dE$$
- **核心 API 映射**：`init_rph_benchmark_reaction`, `calc_generalized_tst_rate`, `calc_cvt_rate_constant`, `calc_eckart_tunneling_factor`。

### 47. 相对论原子结构与径向狄拉克方程 (`mod_relativistic_atomic`)
- **理论基础与物理机制**：重元素体系中相对论效应（质量-速度修正、Darwin 项、自旋-轨道耦合）不可忽视；采用双分量径向狄拉克方程与 Norcross-Klapisch 极化模型势，求解 Sommerfeld 相对论单电子本征能级、微观精细结构劈裂常数以及相对论电偶极 (E1) 振子强度。
- **详细数学表达式**：
  - 径向狄拉克大小分量一阶耦合常微分方程组：
    $$\frac{d}{dr} \begin{pmatrix} P(r) \\ Q(r) \end{pmatrix} = \begin{pmatrix} -\frac{\kappa}{r} & \frac{1}{c} \left( 2 c^2 + E - V(r) \right) \\ -\frac{1}{c} \left( E - V(r) \right) & \frac{\kappa}{r} \end{pmatrix} \begin{pmatrix} P(r) \\ Q(r) \end{pmatrix}$$
    $$\kappa = -(j + 1/2) \cdot \text{sgn}(j - l)$$
  - Sommerfeld 精细结构相对论能量本征值：
    $$E_{n j} = m_e c^2 \left[ \left( 1 + \left( \frac{Z \alpha}{n - (j + 1/2) + \sqrt{(j + 1/2)^2 - (Z\alpha)^2}} \right)^2 \right)^{-1/2} - 1 \right]$$
  - 相对论电偶极 (E1) 径向矩阵元与吸收振子强度：
    $$R_{i \to f} = \int_0^\infty \left[ P_i(r) P_f(r) + Q_i(r) Q_f(r) \right] r \, dr, \quad f_{if} = \frac{2 m_e}{3 \hbar^2} (E_f - E_i) \frac{\max(j_i, j_f)}{2 j_i + 1} |R_{i \to f}|^2$$
- **核心 API 映射**：`calc_dirac_model_potential`, `solve_radial_dirac_eigenvalue`, `calc_dirac_fine_structure_splitting`, `calc_dirac_e1_matrix_element`。

### 48. 共振非弹性 X 射线散射与内壳层光谱 (`mod_resonant_xray_scattering`)
- **理论基础与物理机制**：共振非弹性 X 射线散射 (RIXS) 为光子入-光子出的二阶共振光谱过程；初态芯电子被 X 射线光子跃迁激发至中间导带或自电离态，随后高能价电子退激跃迁填补芯孔并辐射发射出次级光子；借助 Kramers-Heisenberg 二阶极化微扰公式求解电子-声子耦合 Huang-Rhys 振动伴线展开与低能电子元激发。
- **详细数学表达式**：
  - Kramers-Heisenberg 二阶极化散射微扰截面：
    $$\frac{d^2 \sigma}{d\Omega d\omega_2} = \frac{\omega_2}{\omega_1} \sum_f \left| \sum_m \frac{\langle f | \hat{\mathbf{e}}_2^* \cdot \hat{\mathbf{D}} | m \rangle \langle m | \hat{\mathbf{e}}_1 \cdot \hat{\mathbf{D}} | i \rangle}{E_i - E_m + \hbar\omega_1 + i \Gamma_m / 2} \right|^2 \delta(E_i - E_f + \hbar\omega_1 - \hbar\omega_2)$$
  - 光学定理共振 X 射线吸收截面 (XAS)：
    $$\sigma_{\text{XAS}}(\omega_1) = 4\pi^2 \alpha \hbar\omega_1 \sum_m |\langle m | \hat{\mathbf{e}}_1 \cdot \hat{\mathbf{D}} | i \rangle|^2 \frac{\Gamma_m / (2\pi)}{(E_m - E_i - \hbar\omega_1)^2 + (\Gamma_m / 2)^2}$$
  - 电-声耦合相联拉盖尔多项式 Franck-Condon 伴线强度展开：
    $$I(n, \omega_{\text{loss}}) \propto \exp(-S) \frac{S^n}{n!} \left| \sum_{m=0}^\infty \frac{e^{-S} (-1)^m \sqrt{m! n!} \sum_{l=0}^{\min(m, n)} \frac{S^{l} (-1)^l}{l! (m-l)! (n-l)!}}{\Delta_m + i \Gamma_m / 2} \right|^2$$
- **核心 API 映射**：`init_rixs_system`, `calc_xas_cross_section`, `calc_kramers_heisenberg_cross_section`, `calc_rixs_2d_map`, `calc_huang_rhys_vibrational_rixs`。

### 49. 亚稳态原子碰撞潘宁电离与缔合电离 (`mod_penning_associative_ionization`)
- **理论基础与物理机制**：亚稳态原子 $A^*$（如 $\text{He}^*(2^3S, 2^1S)$）的巨大电子激发能超过靶原子/分子 $B$ 的电离能 $I_p$，在微观碰撞过程中发生自发无辐射自电离；反应通道解离分支为潘宁电离（PI：$A^* + B \to A + B^+ + e^-$）与缔合电离（AI：$A^* + B \to AB^+ + e^-$）；利用复光学势 $V_{\text{opt}}(R) = V_*(R) - \frac{i}{2}\Gamma(R)$ 求解半经典存活几率与电离电子能谱 (PIES)；超冷能区中，自旋极化可使非弹性化学电离损失速率被压制数个数量级。
- **详细数学表达式**：
  - 复光学势与半经典初态碰撞存活几率：
    $$V_{\text{opt}}(R) = V_*(R) - \frac{i}{2}\Gamma(R), \quad P_{\text{surv}}(b, E) = \exp\left( -2 \int_{R_{\text{turn}}}^\infty \frac{\Gamma(R)}{\hbar v_r(R)} dR \right)$$
    $$v_r(R) = \sqrt{\frac{2}{\mu} \left( E - V_*(R) - \frac{E b^2}{R^2} \right)}$$
  - 潘宁电离 (PI)、缔合电离 (AI) 与总化学电离截面：
    $$\sigma_{\text{tot}}(E) = 2\pi \int_0^\infty b [1 - P_{\text{surv}}(b, E)] db = \sigma_{\text{PI}}(E) + \sigma_{\text{AI}}(E)$$
    $$\sigma_{\text{AI}}(E) = 2\pi \int_0^\infty b \, db \int_{R_{\text{turn}}}^{R_c(b)} \frac{\Gamma(R)}{\hbar v_r(R)} \exp\left( -2 \int_{R_{\text{turn}}}^R \frac{\Gamma(R')}{\hbar v_r(R')} dR' \right) dR, \quad \left(V_+(R_c) + \frac{E b^2}{R_c^2} = E\right)$$
  - 潘宁电离电子能谱 (PIES) 局域静止相条件：
    $$E_{\text{elec}}(R) = V_*(R) - V_+(R), \quad \frac{d\sigma}{dE_e} \propto \sum_{R_*} \frac{R_*^2 \Gamma(R_*)}{|\frac{d}{dR}[V_*(R) - V_+(R)]|_{R_*}} \sqrt{1 - \frac{V_*(R_*)}{E}}$$
  - 超冷复散射长度与自旋极化自电离抑制比：
    $$a = \alpha - i\beta, \quad K_{\text{loss}} = \frac{4\pi \hbar}{\mu} \beta, \quad \rho_{\text{suppress}} = \frac{K_{\text{loss}}(\text{unpolarized})}{K_{\text{loss}}(\text{spin-polarized})} \sim 10^3 \sim 10^5$$
- **核心 API 映射**：`init_penning_system`, `calc_penning_classical_turning_point`, `calc_penning_cross_sections`, `calc_pies_spectrum`, `calc_penning_thermal_rate`, `calc_ultracold_penning_rates`。

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
41. Penning & Associative Ioniz Tests:  5 /  5 PASSED
----------------------------------------------------------------
ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass, 340/340 断言通过)
================================================================
```

---

## 📊 典型物理算例 (Examples)

位于 `GeneralModule/examples/`，一键编译运行全部 36 大物理前沿算例：
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
36. **`ex36_penning_associative_ionization.f90`**
    - **亚稳态 He*(2^3S)+Ar 潘宁电离与缔合电离空间分流比、PIES 电子能谱与自旋极化抑制**：构建高激发亚稳态中性原子碰撞光学势模型，演示极低碰撞动能下缔合电离 $AB^+$ 形成占据绝对主导（分流比 $>85\%$），随动能增大平滑渡越至解离潘宁电离 $A+B^++e^-$ 主导（高能区 AI 占比降至 0%）；精确输出潘宁电离电子发射能谱 (PIES) 奇异峰；计算微开尔文超冷区复散射长度损耗，复现自旋极化禁阻自旋反平行电离通道高达 $10^4$ 倍的超冷量子气体寿命抑制效应。

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
