# GeneralModule 算法库详细配置与集成部署指南

本文档为 `GeneralModule` 现代量子动力学算法库的权威配置、编译构建、参数调优及多语言混合编程指南。无论您是高校/科研院所的研究人员、高性能计算（HPC）工程师，还是配合 AI Coding Agent（如 Antigravity、Cursor、Copilot）进行算法开发，本手册都提供了完备的配置说明。

> [!NOTE]
> 理论公式推导、物理机理与国际权威期刊参考文献（PRL, PRA, RMP, JCP 等）详见独立全典：
> 👉 **[LITERATURE.md](LITERATURE.md)**

---

## 目录
1. [系统环境与编译器配置](#1-系统环境与编译器配置)
   - [1.1 编译器版本要求](#11-编译器版本要求)
   - [1.2 macOS 环境特殊配置（Xcode 许可绕过）](#12-macos-环境特殊配置xcode-许可绕过)
   - [1.3 Linux (Ubuntu / CentOS / Arch) 配置](#13-linux-ubuntu--centos--arch-配置)
   - [1.4 Windows (WSL2 / MSYS2) 配置](#14-windows-wsl2--msys2-配置)
   - [1.5 Intel oneAPI (ifx/ifort) 配置](#15-intel-oneapi-ifxifort-配置)
2. [多构建系统配置指南](#2-多构建系统配置指南)
   - [2.1 Fortran Package Manager (fpm) 配置](#21-fortran-package-manager-fpm-配置)
   - [2.2 CMake 现代目标化构建配置](#22-cmake-现代目标化构建配置)
   - [2.3 零外部依赖通用 Makefile 模板](#23-零外部依赖通用-makefile-模板)
   - [2.4 纯命令行快速编译与打包](#24-纯命令行快速编译与打包)
3. [算法库参数配置全典 (Configuration Reference)](#3-算法库参数配置全典-configuration-reference)
   - [3.1 激光脉冲配置 (`pulse_config_t`)](#31-激光脉冲配置-pulse_config_t)
   - [3.2 离散变量网格配置 (`dvr_1d_t`, `dvr_legendre_t`)](#32-离散变量网格配置-dvr_1d_t-dvr_legendre_t)
   - [3.3 复吸收边界配置 (`absorbing_boundary_t`)](#33-复吸收边界配置-absorbing_boundary_t)
   - [3.4 强场原子模型配置 (`atom_config_t`)](#34-强场原子模型配置-atom_config_t)
   - [3.5 分子转振态基底与激光动力学耦合配置 (`mod_rovibrational`)](#35-分子转振态基底与激光动力学耦合配置-mod_rovibrational)
   - [3.6 科学计算 I/O 与诊断可视化工具 (`mod_io_utils`)](#36-科学计算-io-与诊断可视化工具-mod_io_utils)
   - [3.7 时间推进步长稳定性判据 (CFL 条件)](#37-时间推进步长稳定性判据-cfl-条件)
   - [3.8 高精度三次样条插值与势能面外推配置 (`mod_interpolation`)](#38-高精度三次样条插值与势能面外推配置-mod_interpolation)
   - [3.9 自相关函数、光吸收谱与光碎片通量配置 (`mod_photofragment_flux`)](#39-自相关函数光吸收谱与光碎片通量配置-mod_photofragment_flux)
   - [3.10 开放量子系统与 Lindblad 耗散主方程配置 (`mod_open_quantum`)](#310-开放量子系统与-lindblad-耗散主方程配置-mod_open_quantum)
   - [3.11 量子最优控制理论 Krotov 算法配置 (`mod_optimal_control`)](#311-量子最优控制理论-krotov-算法配置-mod_optimal_control)
   - [3.12 非含时散射理论与超冷散射长度配置 (`mod_ti_scattering`)](#312-非含时散射理论与超冷散射长度配置-mod_ti_scattering)
   - [3.13 含时波包散射理论与 S-矩阵元提取配置 (`mod_td_scattering`)](#313-含时波包散射理论与-s-矩阵元提取配置-mod_td_scattering)
   - [3.14 外加电磁场超冷散射与多基组密耦配置 (`mod_field_scattering`)](#314-外加电磁场超冷散射与多基组密耦配置-mod_field_scattering)
   - [3.15 各向异性偶极超冷散射与极性分子配置 (`mod_dipolar_scattering`)](#315-各向异性偶极超冷散射与极性分子配置-mod_dipolar_scattering)
   - [3.16 超冷光缔合谱学与分子生成配置 (`mod_photoassociation`)](#316-超冷光缔合谱学与分子生成配置-mod_photoassociation)
   - [3.17 超冷少体物理与 Efimov 三体复合配置 (`mod_three_body_recombination`)](#317-超冷少体物理与-efimov-三体复合配置-mod_three_body_recombination)
   - [3.18 低维光晶格受限量子散射与 CIR 配置 (`mod_confined_scattering`)](#318-低维光晶格受限量子散射与-cir-配置-mod_confined_scattering)
   - [3.19 自电离体系与 Fano 共振 / 复坐标旋转法配置 (`mod_autoionization_fano`)](#319-自电离体系与-fano-共振--复坐标旋转法配置-mod_autoionization_fano)
   - [3.20 交叉静电磁场分子量子动力学配置 (`mod_crossed_field_scattering`)](#320-交叉静电磁场分子量子动力学配置-mod_crossed_field_scattering)
   - [3.21 三原子反应散射 Jacobi 几何与 LEPS 势能面配置 (`mod_triatomic_geometry`)](#321-三原子反应散射-jacobi-几何与-leps-势能面配置-mod_triatomic_geometry)
   - [3.22 超冷旋量玻色爱因斯坦凝聚自旋动力学配置 (`mod_spinor_bec`)](#322-超冷旋量玻色爱因斯坦凝聚自旋动力学配置-mod_spinor_bec)
   - [3.23 三原子超球面反应动力学与热速率常数配置 (`mod_hyperspherical_reactive`)](#323-三原子超球面反应动力学与热速率常数配置-mod_hyperspherical_reactive)
   - [3.24 超冷偶极量子液滴与李-黄-杨量子涨落配置 (`mod_dipolar_droplets_lhy`)](#324-超冷偶极量子液滴与李-黄-杨量子涨落配置-mod_dipolar_droplets_lhy)
   - [3.25 强场非顺序双电离与电子重碰撞相关动量谱配置 (`mod_strong_field_nsdi`)](#325-强场非顺序双电离与电子重碰撞相关动量谱配置-mod_strong_field_nsdi)
   - [3.26 磁与光 Feshbach 共振与弱束缚分子态配置 (`mod_feshbach_bound_states`)](#326-磁与光-feshbach-共振与弱束缚分子态配置-mod_feshbach_bound_states)
   - [3.27 阿秒瞬态吸收光谱与光诱导态自电离干涉 (`mod_attosecond_transient_absorption`)](#327-阿秒瞬态吸收光谱与光诱导态自电离干涉-mod_attosecond_transient_absorption)
   - [3.28 双色反向圆偏振场与分子光电子圆二色性 PECD (`mod_bicircular_pecd`)](#328-双色反向圆偏振场与分子光电子圆二色性-pecd-mod_bicircular_pecd)
   - [3.29 超冷极性分子反应动力学与微波/静电偶极遮蔽 (`mod_ultracold_reaction_shielding`)](#329-超冷极性分子反应动力学与微波静电偶极遮蔽-mod_ultracold_reaction_shielding)
   - [3.30 里德堡原子阻塞、PXP 约束模型与量子多体疤痕 (`mod_rydberg_blockade`)](#330-里德堡原子阻塞pxp-约束模型与量子多体疤痕-mod_rydberg_blockade)
   - [3.31 表面量子散射与选择性吸附共振 (`mod_surface_scattering`)](#331-表面量子散射与选择性吸附共振-mod_surface_scattering)
   - [3.32 气-固界面催化反应与 Eley-Rideal 提取机理 (`mod_surface_reaction_er`)](#332-气-固界面催化反应与-eley-rideal-提取机理-mod_surface_reaction_er)
   - [3.33 金属表面非绝热动力学与电子摩擦耗散 (`mod_surface_electronic_friction`)](#333-金属表面非绝热动力学与电子摩擦耗散-mod_surface_electronic_friction)
   - [3.34 掠入射快原子表面量子衍射与彩虹散射 (`mod_grazing_fast_atom_diffraction`)](#334-掠入射快原子表面量子衍射与彩虹散射-mod_grazing_fast_atom_diffraction)
   - [3.35 冷离子-中性原子杂化散射与极化阱动力学 (`mod_ion_atom_scattering`)](#335-冷离子-中性原子杂化散射与极化阱动力学-mod_ion_atom_scattering)
   - [3.36 最少开关表面跳跃与非绝热混合量子-经典动力学 (`mod_surface_hopping_fssh`)](#336-最少开关表面跳跃与非绝热混合量子-经典动力学-mod_surface_hopping_fssh)
   - [3.37 强场分子定向、取向与超转子动力学 (`mod_molecular_alignment`)](#337-强场分子定向取向与超转子动力学-mod_molecular_alignment)
   - [3.38 超冷光晶格与玻色-哈伯德微观映射 (`mod_optical_lattice_hubbard`)](#338-超冷光晶格与玻色-哈伯德微观映射-mod_optical_lattice_hubbard)
   - [3.39 多原子反应路径哈密顿量与变分过渡态理论 (`mod_reaction_path_hamiltonian`)](#339-多原子反应路径哈密顿量与变分过渡态理论-mod_reaction_path_hamiltonian)
   - [3.40 相对论原子结构与径向狄拉克方程 (`mod_relativistic_atomic`)](#340-相对论原子结构与径基狄拉克方程-mod_relativistic_atomic)
   - [3.41 共振非弹性 X 射线散射与内壳层光谱 (`mod_resonant_xray_scattering`)](#341-共振非弹性-x-射线散射与内壳层光谱-mod_resonant_xray_scattering)
   - [3.42 亚稳态原子碰撞潘宁电离与缔合电离 (`mod_penning_associative_ionization`)](#342-亚稳态原子碰撞潘宁电离与缔合电离-mod_penning_associative_ionization)
4. [Python 伴侣库 `pygenmod` 配置与混合编程](#4-python-伴侣库-pygenmod-配置与混合编程)
   - [4.1 本地可编辑模式安装](#41-本地可编辑模式安装)
   - [4.2 数据交互规范（.dat 与无损二进制）](#42-数据交互规范-dat-与无损二进制)
   - [4.3 自动化发表级绘图管道](#43-自动化发表级绘图管道)
5. [AI Agent 自动化编程与提示词工程模板](#5-ai-agent-自动化编程与提示词工程模板)
   - [5.1 给 AI 的上下文 System Prompt 模板](#51-给-ai-的上下文-system-prompt-模板)
   - [5.2 常见物理场景任务提示词 (User Prompts)](#52-常见物理场景任务提示词-user-prompts)
   - [5.3 AI 编写 Fortran 常见陷阱自查表](#53-ai-编写-fortran-常见陷阱自查表)
6. [常见错误排查自查表 (Troubleshooting & FAQ)](#6-常见错误排查自查表-troubleshooting--faq)

---

## 1. 系统环境与编译器配置

### 1.1 编译器版本要求
`GeneralModule` 采用现代 **Fortran 2008** 标准编写，无须任何第三方 C/C++ 运行时。
- **GNU Fortran (`gfortran`)**: 推荐 GCC 9.0 或更高版本（已在 GCC 12, 13, 14, 15 上完成全面验证）。
- **Intel Fortran (`ifx` / `ifort`)**: 推荐 Intel oneAPI HPC Toolkit 2023.0 或更新版本。
- **LLVM Flang**: Flang 16+ 支持。

### 1.2 macOS 环境特殊配置（Xcode 许可绕过）
在 macOS 平台（尤其使用 Homebrew 安装的 GCC/gfortran 时），系统汇编器 `as` 或归档器 `ar` 可能会因系统未接受 Xcode 许可而报错：
```text
You have not agreed to the Xcode license agreements. Please run 'sudo xcodebuild -license'...
```
**标准解决方案**：在 shell 环境中指定 CommandLineTools 路径（无需 root 权限或同意完整 Xcode 协议）：
```bash
# 临时生效
export DEVELOPER_DIR="/Library/Developer/CommandLineTools"

# 永久生效（写入 ~/.zshrc 或 ~/.bash_profile）
echo 'export DEVELOPER_DIR="/Library/Developer/CommandLineTools"' >> ~/.zshrc
source ~/.zshrc
```
在此环境变量下，所有 `gfortran`、`ar`、`ranlib` 命令均可直接以极速本地执行。

### 1.3 Linux (Ubuntu / CentOS / Arch) 配置
- **Ubuntu / Debian**:
  ```bash
  sudo apt-get update
  sudo apt-get install -y gfortran cmake build-essential
  ```
- **CentOS / RHEL / Rocky Linux**:
  ```bash
  sudo dnf install -y gcc-gfortran cmake make
  ```
- **Arch Linux**:
  ```bash
  sudo pacman -S gcc-fortran cmake
  ```

### 1.4 Windows (WSL2 / MSYS2) 配置
推荐在 Windows 10/11 上通过 **WSL2 (Ubuntu 22.04/24.04 LTS)** 使用，体验与原生 Linux 一致。若使用原生 Windows 环境，推荐安装 **MSYS2 UCRT64**：
```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-cmake make
```

### 1.5 Intel oneAPI (ifx/ifort) 配置
若需追求极限向量化性能，可使用 Intel 现代编译器 `ifx`：
```bash
source /opt/intel/oneapi/setvars.sh
ifx -O3 -xHost -fPIC -c src/*.f90
```

---

## 2. 多构建系统配置指南

### 2.1 Fortran Package Manager (fpm) 配置
`fpm` 是 Fortran 官方现代包管理器，类似 Rust 的 `cargo` 或 Python 的 `poetry`。

#### 1) 在 `GeneralModule` 根目录执行
```bash
# 构建整个库
fpm build

# 运行全套自动化单元测试
fpm test

# 运行特定测试（例如测试原子强场物理模块）
fpm test test_atomic_hhg

# 运行特定物理算例（例如运行高次谐波模拟）
fpm run --example ex05_hhg_lewenstein_spectrum
```

#### 2) 在外部项目中使用 `fpm` 引用 `GeneralModule`
在您自己工程的 `fpm.toml` 中添加本库为本地或 Git 依赖：
```toml
[dependencies]
general_module = { path = "/绝对路径或相对路径/to/GeneralModule" }
# 或通过 Git 仓库依赖:
# general_module = { git = "https://github.com/your-username/GeneralModule.git", tag = "v1.1.0" }
```

---

### 2.2 CMake 现代目标化构建配置
`GeneralModule` 根目录提供了标准 `CMakeLists.txt`。

#### 1) 编译静态库与动态库
```bash
cd GeneralModule
mkdir build && cd build

# macOS / Linux 配置
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j4

# 运行测试
ctest --output-on-failure
```
构建完成后将在 `build/` 目录下生成：
- 静态库：`libgeneral_module.a`
- 动态共享库：`libgeneral_module.dylib` (macOS) 或 `libgeneral_module.so` (Linux)
- 模块定义文件：`general_module.mod` 及各子模块 `.mod` 文件

#### 2) 在下游 CMake 项目中直接引用
在您的项目的 `CMakeLists.txt` 中添加：
```cmake
cmake_minimum_required(VERSION 3.14)
project(MyQuantumSimulation LANGUAGES Fortran)

# 添加 GeneralModule 子目录
add_subdirectory(/path/to/GeneralModule GeneralModule_build)

# 声明您的可执行目标
add_executable(run_dynamics src/main.f90)

# 直接链接 GeneralModule 目标（自动继承 include 路径与编译参数）
target_link_libraries(run_dynamics PRIVATE GeneralModule_static)
```

---

### 2.3 零外部依赖通用 Makefile 模板
对于习惯使用经典 `Makefile` 的项目，可直接使用以下模板：

```makefile
# ==============================================================================
# 通用量子动力学工程 Makefile
# ==============================================================================
FC = gfortran
FFLAGS = -O3 -fPIC -Wall

# GeneralModule 路径配置
GENMOD_DIR = /path/to/GeneralModule
GENMOD_SRC = $(GENMOD_DIR)/src
GENMOD_LIB = $(GENMOD_SRC)/libgeneral_module.a

# 目标可执行程序
TARGET = my_simulation
SRCS = main.f90

all: $(TARGET)

# 自动检查并构建 GeneralModule 静态库
$(GENMOD_LIB):
	@echo ">> Compiling GeneralModule library..."
	cd $(GENMOD_SRC) && $(FC) $(FFLAGS) -c *.f90
	cd $(GENMOD_SRC) && ar rcs libgeneral_module.a *.o

# 链接用户工程
$(TARGET): $(SRCS) $(GENMOD_LIB)
	$(FC) $(FFLAGS) -I$(GENMOD_SRC) $(SRCS) $(GENMOD_LIB) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o *.mod

.PHONY: all clean
```

---

### 2.4 纯命令行快速编译与打包
如果您无需构建系统，仅想在几秒内构建并测试：
```bash
cd GeneralModule/src
# 1. 编译全部模块生成目标文件与 .mod
DEVELOPER_DIR=/Library/Developer/CommandLineTools gfortran -O3 -fPIC -c *.f90

# 2. 打包为静态库
DEVELOPER_DIR=/Library/Developer/CommandLineTools ar rcs libgeneral_module.a *.o

# 3. 编译任意调用程序
DEVELOPER_DIR=/Library/Developer/CommandLineTools gfortran -O3 -I/path/to/GeneralModule/src my_prog.f90 /path/to/GeneralModule/src/libgeneral_module.a -o my_prog
```

---

## 3. 算法库参数配置全典 (Configuration Reference)

### 3.1 激光脉冲配置 (`pulse_config_t`)
用于 `mod_laser_pulse` 模块中各种形状超快脉冲的精确时域合成：

```fortran
use general_module
type(pulse_config_t) :: laser
```

| 字段名称 | 类型 | 物理意义与默认值 | 推荐配置取值与单位 |
| :--- | :--- | :--- | :--- |
| `shape_type` | `integer` | 脉冲包络形态，默认 `PULSE_GAUSSIAN` | `PULSE_GAUSSIAN` (1), `PULSE_SIN2` (2), `PULSE_FLATTOP` (3), `PULSE_CHIRP` (4), `PULSE_TWOCOLOR` (5), `PULSE_THZ_TRAIN` (6) |
| `field_peak` | `real(dp)` | 激光电场峰值 $E_0$，默认 `0.0_dp` | $0.01 \sim 0.15\text{ a.u.}$ ($10^{13} \sim 10^{15}\text{ W/cm}^2$)。转换公式：$E_0 = \sqrt{I / 3.51 \times 10^{16}}$ |
| `freq_central` | `real(dp)` | 载波中心角频率 $\omega_0$，默认 `0.0_dp` | $800\text{ nm} \rightarrow 0.057\text{ a.u.}$, $400\text{ nm} \rightarrow 0.114\text{ a.u.}$, $1300\text{ nm} \rightarrow 0.035\text{ a.u.}$ |
| `duration` | `real(dp)` | 半高全宽（FWHM）持续时间 $\tau$ | 使用 `to_au(30.0_dp, "fs")`，常用 $5 \sim 100\text{ fs}$ |
| `t_center` | `real(dp)` | 脉冲峰值中心时刻 $t_0$ | 常用 `0.0_dp` |
| `chirp_rate` | `real(dp)` | 瞬时频率线性啁啾率 $\beta$ | 正啁啾 $>0$，负啁啾 $<0$（如 `1.0e-5_dp a.u.`） |
| `cep_phase` | `real(dp)` | 载波包络相位 $\phi_{CEP}$ | $0.0 \sim 2\pi\text{ rad}$（调控亚周期波形极性） |
| `ellipticity` | `real(dp)` | 椭偏率 $\epsilon \in [-1, 1]$ | `0.0`（纯线偏振）, `1.0`（左旋圆偏振）, `-1.0`（右旋圆偏振） |
| `two_color_ratio` | `real(dp)` | 双色场二次谐波场强比 $E_{2\omega} / E_\omega$ | 常用 $0.1 \sim 0.4$（破坏空间反演对称性） |
| `two_color_phase` | `real(dp)` | 双色相对相位 $\Delta\phi$ | 常用 `0.0` 或 `HALFPI`（控制电子不对称定向发射） |
| `train_count` | `integer` | 太赫兹单周期脉冲个数 | 常用 $1 \sim 5$ |
| `train_delay` | `real(dp)` | 脉冲序列脉冲间延迟 $\Delta t_{delay}$ | 通常设为转动周期的分数或整倍数（如 $100\text{ fs} \sim 2\text{ ps}$） |

---

### 3.2 离散变量网格配置 (`dvr_1d_t`, `dvr_legendre_t`)

#### 1) 一维 Sinc-DVR 网格 (`dvr_1d_t`)
基于 Colbert-Miller 正弦基矢表象，用于分子核间距 $R$ 或电子坐标 $x$：
```fortran
call dvr_sinc_init(x_min, x_max, n_pts, mass, dvr)
```
- `x_min`, `x_max`: 空间物理边界（原子单位 Bohr）。束缚态模拟时应超出势能经典禁区足够远（例如振动态取 $0.8 \sim 8.0\text{ Bohr}$，强场光解离取 $-50 \sim 100\text{ Bohr}$）。
- `n_pts`: 网格点数。推荐取 $128, 256, 512$ 等 $2^N$ 形式（以便后续无缝衔接 FFT）。
- `mass`: 体系约化质量（原子单位电子质量 $m_e$）。$1\text{ amu} \approx 1822.8885\text{ a.u.}$。调用 `to_au(mass_amu, "amu")` 传入。

#### 2) 角向 Legendre-DVR 网格 (`dvr_legendre_t`)
用于分子极角 $\theta$ 取向与转动动力学：
```fortran
call dvr_legendre_init(n_pts, dvr_leg)
```
- `n_pts`: 角向 Gauss-Legendre 节点数（通常 $32 \sim 64$ 个节点足以精准展开高达 $J=30$ 的转动波包）。

---

### 3.3 复吸收边界配置 (`absorbing_boundary_t`)
用于在有限空间网格边界消除非物理波包反射，模拟光解离、光电离及连续态吸收：
```fortran
call cap_init(r_start, r_end, strength, cap_type, cap_obj)
```

| 参数 | 物理意义 | 典型推荐设置值 |
| :--- | :--- | :--- |
| `r_start` | 吸收起始位置 $R_0$ | 必须设置在感兴趣的内部反应区外侧（如 $R_{max} - 15\text{ a.u.}$） |
| `r_end` | 吸收终止边界 $R_{max}$ | 空间网格的最大外边界 |
| `strength` | 吸收强度系数 $A$ | 通常取 $0.1 \sim 0.5\text{ a.u.}$。**注意**：过大可能引起边界反射势垒，过小则吸收不充分 |
| `cap_type` | 吸收函数形态 | `CAP_SIN2` (首选，平滑无突变) 或 `CAP_POLYNOMIAL` |

---

### 3.4 强场原子模型配置 (`atom_config_t`)
内置高精度单活性电子（SAE）模型参数：
```fortran
call get_atom_config("Ar", atom)
```
内置预设包含：`"H"` ($I_p=13.606\text{ eV}$), `"He"` ($24.587\text{ eV}$), `"Ne"` ($21.565\text{ eV}$), `"Ar"` ($15.760\text{ eV}$), `"Kr"` ($13.999\text{ eV}$), `"Xe"` ($12.130\text{ eV}$)。
- 自定义势阱参数：直接修改 `atom%ip_au` 与软核参数 `atom%soft_core_a`。

---

### 3.5 分子转振态基底与激光动力学耦合配置 (`mod_rovibrational`)
针对极性双原子分子超快红外与强场激光相干调控，算法库提供了完整的 $|v, J\rangle$ 转振空间基底与跃迁矩阵构建：

1. **转振基底索引映射**：
   ```fortran
   ! 将 (v, J) 映射至一维线性索引 k ∈ [1, (v_max+1)*(j_max+1)]
   k = rovibrational_state_index(v, j, j_max)
   ! 逆向反解
   call rovibrational_state_unindex(k, j_max, v, j)
   ```
2. **转振能级与偶极跃迁矩阵构建**：
   ```fortran
   ! 构造无场本征能级: E(v, J) = E_vib(v) + B_v * J * (J + 1)
   call build_rovibrational_hamiltonian(v_max, j_max, e_vib, b_v, h_diag)
   ! 构造偶极矩阵元 (严格满足 Delta J = +/-1, Delta M = 0)
   call build_rovibrational_dipole_matrix(v_max, j_max, dip_vib, dip_mat)
   ```
3. **STIRAP 绝热通道脉冲对配置**：
   ```fortran
   ! 自动生成 Stokes 先于 Pump 的反直觉时序脉冲对
   call create_stirap_pulses(peak_p, peak_s, dur_p_fs, dur_s_fs, delay_fs, &
                             freq_p_au, freq_s_au, cfg_pump, cfg_stokes)
   ```

---

### 3.6 科学计算 I/O 与诊断可视化工具 (`mod_io_utils`)
消除低效的手工格式化文件读写与控制台杂乱输出，内置发表级数据交换与监测例程：
- `save_data_table_1d(filename, x, y, x_name, y_name, header)`：标准化一维双列 ASCII 保存。
- `save_data_table_2d(filename, time, pop_mat, col_names, header)`：多自由度或转振态布居含时序列矩阵导出。
- `save_matrix_dat(filename, mat, header)`：矩阵元或势能面网格导出。
- `print_banner(title, width)`：控制台美化横幅输出。
- `print_progress_bar(step, total, prefix)`：含时动力学推进循环中单行就地刷新进度百分比与进度条。

---

### 3.7 时间推进步长稳定性判据 (CFL 条件)
使用 `propagate_split_operator_1d` 或 `propagate_split_operator_2d` 进行时间演化时，为保证波包演化的辛对称幺正性与相位精度，时间步长 $\Delta t$ 必须满足空间动能离散的 Courant-Friedrichs-Lewy (CFL) 上限：
$$\Delta t \le \frac{2 m \Delta x^2}{\pi \hbar}$$
- **电子动力学（$m = 1\text{ a.u.}$）**：若 $\Delta x = 0.2\text{ Bohr}$，则 $\Delta t \le \frac{2 \times 0.04}{3.14} \approx 0.025\text{ a.u.} \approx 0.6\text{ 自动单位 (约 } 0.0006\text{ fs)}$。
- **核动力学（$m \sim 2000\text{ a.u.}$）**：步长可大幅放宽至 $\Delta t \sim 1.0 - 5.0\text{ a.u.} (0.02 - 0.1\text{ fs})$。
- **切比雪夫推进器（`chebyshev_propagate_step`）**：时间步长不受上述高频震荡严格限制，单步可达几个飞秒且保持机器精度。

---

### 3.8 高精度三次样条插值与势能面外推配置 (`mod_interpolation`)
针对量子化学第一性原理计算离散单点能（Ab Initio PES）或从头算偶极矩曲面的高精度平滑重构：
1. **构造自然/固定导数样条**：
   ```fortran
   type(spline_1d_t) :: pes_spline
   ! 自然边界条件 (BC_NATURAL: y''(x_0) = y''(x_n) = 0)
   call spline_1d_init(r_grid, v_ab_initio, pes_spline, bc_type=BC_NATURAL)
   ! 固定一阶导数边界条件 (BC_CLAMPED)
   call spline_1d_init(r_grid, v_ab_initio, pes_spline, bc_type=BC_CLAMPED, yp_0=0.0_dp, yp_n=0.0_dp)
   ```
2. **内插与连续导数求值**：
   ```fortran
   v_val = spline_1d_eval(pes_spline, r_curr)       ! 高精度函数值
   force = -spline_1d_deriv(pes_spline, r_curr)     ! 解析一阶导数（受力）
   curv  = spline_1d_deriv2(pes_spline, r_curr)     ! 解析二阶导数（曲率）
   ```
3. **势能面全域物理外推（短程排斥 + 长程范德华）**：
   ```fortran
   ! 当 r < r_min 时采用 A*exp(-B*r) 指数排斥，当 r > r_max 时采用 V_inf - C6/r^6 平滑过渡
   v_extrap = potential_extrapolate_1d(r_curr, pes_spline, r_min=1.0_dp, r_max=12.0_dp, &
                                       a_rep=100.0_dp, b_rep=2.5_dp, c6_disp=25.0_dp, v_inf=0.0_dp)
   ```

---

### 3.9 自相关函数、光吸收谱与光碎片通量配置 (`mod_photofragment_flux`)
面向分子光解离（Photodissociation）、光缔合与超快激发态动力学终态测量分析：
1. **波包自相关函数与 Heller 吸收截面谱**：
   ```fortran
   ! 计算每步含时重叠 C(t) = <psi(0) | psi(t)>
   c_t(it) = calc_autocorrelation(psi_init, psi_curr, dx)

   ! 快速傅里叶时频积分提取吸收截面 sigma(omega)
   call heller_absorption_spectrum(t_arr, c_t, dt, gamma_damp=0.002_dp, &
                                   omega_arr=omega_arr, spectrum=sigma_abs, e_zero=e_init)
   ```
2. **渐近面动能释放谱（KER）与多通道光解离分支比**：
   ```fortran
   ! 在探测边界 R_det 处提取连续态能量谱振幅 A(E)
   call photofragment_energy_amplitude(t_arr, psi_at_det, dt, e_grid, amp_e)

   ! 积分各出射通道渐近概率通量计算分支比 (Branching Ratio)
   call photofragment_branching_ratio(flux_channels, n_channels=2, ratios=branch_ratios)
   ```

---

### 3.10 开放量子系统与 Lindblad 耗散主方程配置 (`mod_open_quantum`)
模拟受热浴环境耗散、自发辐射衰减（$T_1$ 弛豫）与介质碰撞纯退相位（$T_2^*$ 退相干）影响的非幺正密度矩阵演化：
1. **初始化耗散系统与跃迁通道**：
   ```fortran
   type(lindblad_system_t) :: sys
   call lindblad_init(n_levels=3, n_channels=2, sys=sys)

   ! 添加能级 2 到能级 1 的自发跃迁弛豫通道 (rate = gamma_1)
   call lindblad_add_decay_channel(sys, i_from=2, j_to=1, rate_gamma=1.0e-4_dp)

   ! 添加能级 2 的纯退相位通道 (rate = gamma_dephasing)
   call lindblad_add_dephasing_channel(sys, level_idx=2, rate_gamma_d=5.0e-5_dp)
   ```
2. **Runge-Kutta 4 阶密度矩阵积分推进**：
   ```fortran
   ! 单步演化: d(rho)/dt = -i[H, rho] + D[rho]
   call propagate_lindblad_rk4(rho, h_eff, sys, dt)
   ```
3. **量子态纯度、信息熵与相干度诊断**：
   ```fortran
   purity  = density_matrix_purity(rho)      ! Tr(rho^2) ∈ [1/N, 1]
   entropy = von_neumann_entropy(rho)        ! -Tr(rho * ln(rho))
   coher   = quantum_coherence_l1(rho)       ! sum_{i != j} |rho_ij|
   ```

---

### 3.11 量子最优控制理论 Krotov 算法配置 (`mod_optimal_control`)
用于高保真度目标态布居转移（State Preparation）或量子逻辑门激光脉冲波形自适应设计：
1. **控制参数配置**：
   ```fortran
   type(oct_config_t) :: oct_cfg
   ! n_steps: 离散时步, dt: 推进步长
   ! alpha_penalty: 激光能量惩罚因子 (建议原子单位取 20.0 ~ 100.0)
   ! max_iter: 最大迭代轮数, tol_fidelity: 目标保真度阈值 (如 0.999)
   call oct_config_init(n_steps=1000, dt=0.5_dp, alpha_penalty=50.0_dp, &
                        max_iter=50, tol_fidelity=0.999_dp, cfg=oct_cfg)
   ```
2. **端到端激光脉冲自动优化**：
   ```fortran
   call oct_optimize_pulse(h0_mat, mu_mat, psi_init, psi_target, &
                           oct_cfg, field_opt, final_fidelity, stat)
   print *, "Optimization converged with final fidelity:", final_fidelity
   ```

---

### 3.12 非含时散射理论与超冷散射长度配置 (`mod_ti_scattering`)
面向超冷原子/分子碰撞、磁 Feshbach 共振、散射截面及分波相移计算：
1. **零能 Numerov 与 Johnson 对数导数法计算散射长度 $a_s$**：
   ```fortran
   ! 零能 Numerov 算法 (同时导出零能径向波函数 u_zero)
   call calc_scattering_length_numerov(r_grid, v_pot, mass, a_s, u_zero, stat)

   ! Johnson 比值法 (抗数值上溢，专用于深吸引阱体系)
   call calc_scattering_length_logder(r_grid, v_pot, mass, a_s, stat)
   ```
2. **正能量分波相移与全套散射矩阵元 ($K_l, S_l, T_l$)**：
   ```fortran
   ! 求解单分波相移 delta_l、K 矩阵元、幺正 S 矩阵元与 T 矩阵元
   call calc_phase_shift_single_l(r_grid, v_pot, mass, energy=0.05_dp, l=0, &
                                  delta=delta_0, k_mat=k_0, s_mat=s_0, t_mat=t_0)

   ! 多分波截面计算与光学定理自洽校验
   call calc_partial_wave_cross_sections(r_grid, v_pot, mass, energy=0.05_dp, l_max=4, &
                                         delta_arr=delta_arr, sigma_partial=sigma_part, &
                                         sigma_total=sigma_tot)
   sigma_opt = optical_theorem_cross_section(k_wave, delta_arr, l_max=4)
   ```
3. **微分散射截面、全同粒子干涉与输运截面**：
   ```fortran
   ! 微分散射截面 (区分粒子 |f(theta)|^2)
   call calc_differential_cross_section(energy, mass, delta_arr, l_max, theta_grid, ds_dist)

   ! 全同玻色子 (偶数分波增强 4 倍，theta=pi/2 处干涉相长)
   call calc_differential_cross_section_identical(energy, mass, delta_arr, l_max, theta_grid, &
                                                 STAT_IDENTICAL_BOSON, ds_boson)

   ! 极化全同费米子 (奇数分波，theta=pi/2 处由于偶极/宇称相消严格归零)
   call calc_differential_cross_section_identical(energy, mass, delta_arr, l_max, theta_grid, &
                                                 STAT_IDENTICAL_FERMION, ds_fermion)

   ! 动量传输截面 sigma_m 与粘滞截面 sigma_v
   call calc_transport_cross_sections(energy, mass, delta_arr, l_max, sigma_m, sigma_v)

   ! Legendre 多极各向异性展开 A_K 与前后非对称度 A_FB
   call calc_differential_legendre_expansion(theta_grid, ds_dist, n_theta, k_max=4, &
                                            a_k=a_k, a_fb=a_fb)
   ```
4. **超低能区有效力程展开 (ERE: $a_s, r_0$)**：
   ```fortran
   ! 自动在多个动量点提取相移并进行 k*cot(delta_0) = -1/a_s + 0.5*r_0*k^2 拟合
   call fit_effective_range_expansion(r_grid, v_pot, mass, k_list, n_k=4, &
                                      a_s=as_fit, r_0=r0_fit)
   ```
5. **形状共振 (Shape Resonance) Wigner 时延分析**：
   ```fortran
   type(resonance_info_t) :: res
   ! 分析相移跃升峰值提取 Wigner 散射时延与共振线宽 Gamma
   call analyze_shape_resonance(e_grid, delta_grid, n_pts, hbar=1.0_dp, res_info=res)
   print *, "Resonance Energy:", res%e_res, "Width Gamma:", res%gamma_width
   ```
6. **双通道非绝热耦合密耦定态 S-矩阵**：
   ```fortran
   complex(dp) :: s_2x2(2, 2)
   real(dp) :: p_inelastic
   call calc_coupled_channel_smatrix_2x2(r_grid, v11, v22, v12, mass, &
                                         total_energy=0.5_dp, delta_e=0.1_dp, &
                                         s_matrix=s_2x2, inelastic_prob=p_inelastic)
   ```
7. **通用任意 $N$ 通道定态密耦求解器 (Johnson Log-Derivative) 与 Feshbach 共振**：
   ```fortran
   type(multichannel_result_t) :: mc_res
   real(dp) :: v_mat(3, 3, 500), thresholds(3)
   integer  :: l_channels(3)

   thresholds = [0.0_dp, 0.05_dp, 0.40_dp]  ! 渐近通道阈值
   l_channels = [0, 0, 0]                   ! 通道轨道角动量

   ! 求解全通道定态密耦 (自动处理开通道与闭通道 Schur 补变换)
   call calc_multichannel_close_coupling_logder(r_grid, v_mat, mass=1.0_dp, &
                                                total_energy=0.15_dp, &
                                                thresholds=thresholds, &
                                                l_channels=l_channels, &
                                                res=mc_res)

   print *, "开通道数:", mc_res%n_open, "闭通道数:", mc_res%n_closed
   print *, "态-态跃迁几率 P(1->2):", mc_res%prob_matrix(2, 1)
   print *, "S-矩阵幺正性 |S11|^2 + |S12|^2:", mc_res%prob_matrix(1, 1) + mc_res%prob_matrix(2, 1)

   ! 跨 Feshbach 共振能区连续扫描
   call calc_feshbach_resonance_scan(r_grid, v_mat, mass=1.0_dp, &
                                     energy_grid=e_scan, n_energies=100, &
                                     thresholds=thresholds, l_channels=l_channels, &
                                     s_wave_length=as_scan, eigenphase_sums=delta_scan)
   ```
8. **定态连续谱散射能量本征波函数求解与归一化**：
   ```fortran
   real(dp), allocatable :: u_wf(:)
   real(dp) :: delta_phase
   allocate(u_wf(n_pts))

   ! norm_type 可选: NORM_ENERGY (delta(E-E') 能量归一化),
   !                 NORM_MOMENTUM (delta(k-k') 动量归一化),
   !                 NORM_UNIT_AMPLITUDE (渐近振幅为 1.0)
   call calc_scattering_wavefunction_ti(r_grid, v_pot, mass=1.0_dp, energy=0.20_dp, &
                                        l=0, norm_type=NORM_ENERGY, &
                                        u_wf=u_wf, phase_shift=delta_phase)
   ```
9. **多扇区分段网格 (Segmented Grid) 与自适应步长密耦推进**：
   ```fortran
   type(segmented_grid_t) :: grid
   real(dp) :: r_bounds(3), dr_steps(3)
   type(multichannel_result_t) :: mc_res

   ! 1. 构造 3 扇区分段网格：近核深势阱区步长极密，长程色散区稀疏放缩
   r_bounds = [3.0_dp, 15.0_dp, 100.0_dp]
   dr_steps = [0.005_dp, 0.02_dp, 0.10_dp]
   call create_segmented_grid(r_start=0.05_dp, r_bounds=r_bounds, dr_steps=dr_steps, grid=grid)

   ! 2. 分段网格零能散射长度 (Numerov + 4阶 Taylor 导数跨扇区光滑桥接)
   call calc_scattering_length_segmented_numerov(grid, v_pot, mass=1.0_dp, a_s=as_val)

   ! 3. 分段网格多通道定态密耦 (Johnson 矩阵对数导数局域无缝传递)
   call calc_multichannel_close_coupling_segmented_logder( &
       grid, v_mat, mass=1.0_dp, total_energy=0.15_dp, &
       thresholds=thresholds, l_channels=l_channels, res=mc_res)
   ```

---

### 3.13 含时波包散射理论与 S-矩阵元提取配置 (`mod_td_scattering`)
面向单次含时波包动力学推进提取全连续能量谱散射观测量：
1. **构造入射高斯散射波包与动量谱**：
   ```fortran
   ! 坐标表象构造入射包
   call gaussian_wavepacket_1d(x_grid, x0=-12.0_dp, sigma_x=1.5_dp, k0=1.2_dp, psi_0=psi)
   ! 动量表象解析权重
   gk = gaussian_momentum_amplitude(k_val, x0=-12.0_dp, sigma_x=1.5_dp, k0=1.2_dp)
   ```
2. **渐近透射面通量时间-能量傅里叶振幅与透射谱 $T(E)$**：
   ```fortran
   ! 在演化循环中累积透射面振幅 A(E)
   call accumulate_flux_amplitude(t_curr, psi(idx_det), dt, energy_grid, n_energies, hbar, amp_trans)
   ! 碰撞结束后直接计算连续能域透射几率 T(E)
   call calculate_td_transmission(energy_grid, n_energies, amp_trans, mass, hbar, &
                                  x0=-12.0_dp, sigma_x=1.5_dp, k0=1.2_dp, t_prob=t_prob)
   ```
3. **比对自由对照波包提取全能量散射矩阵元 $S(E)$ 与散射相移 $\delta(E)$**：
   ```fortran
   call calculate_td_smatrix_element(amp_scatter, amp_free, n_energies, s_mat_e, phase_shift_e)
   ! 含时 Wigner 散射时延
   call wavepacket_wigner_delay(energy_grid, phase_shift_e, n_energies, hbar, delay_w)
   ```
4. **Möller 渐近动量投影法**：
   ```fortran
   ! 末态波包动量空间分解，一步提取连续态透射与反射概率
   call project_wavepacket_to_smatrix(x_grid, dx, psi_final, mass, hbar, &
                                      k0, sigma_x, x0, energy_grid, n_energies, &
                                      t_prob, r_prob)
   ```
5. **二维含时波包角分布微分散射截面 $\frac{d\sigma}{d\theta}(\theta)$**：
   ```fortran
   ! 提取 2D 波包在检测半径 r_det 处的角度积分微分散射截面
   call calculate_td_differential_cross_section_2d(x_grid=x_grid, y_grid=y_grid, psi_2d=psi_2d, &
                                                  mass=1.0_dp, hbar=1.0_dp, theta_grid=theta_grid, &
                                                  dsigma_dtheta=ds_theta)
   ```
6. **含时波包动力学全空间谱投影提取连续本征波函数 $\psi_E(x)$**：
   ```fortran
   complex(dp), allocatable :: psi_accum(:), psi_energy_norm(:)
   allocate(psi_accum(nx), psi_energy_norm(nx))
   psi_accum = (0.0_dp, 0.0_dp)

   ! 1. 在含时演化推进循环中原位积分累积半傅里叶谱投影
   do step = 1, total_steps
       call accumulate_wavefunction_spectral_projection(psi, t_curr, dt, target_energy, hbar, psi_accum)
       call propagate_split_operator_1d(psi, v_pot, dx, mass, dt)
   end do

   ! 2. 演化结束后归一化提取严格 delta(E-E') 能量归一化连续能量本征函数
   call extract_td_scattering_wavefunction(x_grid, psi_accum, target_energy, mass, hbar, &
                                           x0, sigma_x, k0, psi_energy_norm)
   ```

---

### 3.14 外加电磁场超冷散射与多基组密耦配置 (`mod_field_scattering`)
面向外加磁场（Zeeman 效应）与直流电场（Stark 效应）下的碱金属原子与极性分子超冷量子碰撞散射：
1. **四大表象基组选择与应用场景**：
   - `BASIS_UNCOUPLED`：非耦合基组 $|m_{s1}, m_{i1}, m_{s2}, m_{i2}, L, M_L\rangle$。在外磁场下 Zeeman 相互作用严格对角，适于强场/Paschen-Back 机制与全势矩阵直接组装。
   - `BASIS_F_COUPLED`：单原子超精细耦合基组 $|(s_1 i_1) f_1 m_{f1}, (s_2 i_2) f_2 m_{f2}, L, M_L\rangle$。在零磁场下超精细能级对角，适于弱场与超精细能级态制备。
   - `BASIS_TOTAL_SPIN`：总自旋耦合基组 $|(s_1 s_2) S, (i_1 i_2) I, F M_F, L, M_L\rangle$。电子自旋单重态 $V_0(r)$ 与三重态 $V_1(r)$ 相互作用势直接对角，自旋交换矩阵元 $\hat{\mathbf{s}}_1 \cdot \hat{\mathbf{s}}_2$ 严格对角。
   - `BASIS_FIELD_DRESSED`：渐近本征通道基组 $|\alpha_1(B), \alpha_2(B), L, M_L\rangle$。在无限远渐近区严格对角化单原子外场哈密顿量，消除长程人工非绝热耦合。
2. **冷原子同位素参数库与 Breit-Rabi 能级求解**：
   ```fortran
   type(cold_atom_t) :: rb87
   real(dp), allocatable :: e_levels(:), states(:, :)
   integer :: n_states

   ! 提取内置同位素 (支持 6Li, 7Li, 23Na, 40K, 87Rb, 133Cs, 40K87Rb)
   call get_cold_atom_preset("87Rb", rb87)

   ! 求解 B = 50.0 Gauss 下的单原子 Zeeman-超精细 Breit-Rabi 本征态
   call calc_breit_rabi_energies(rb87, 50.0_dp, e_levels, states, n_states)
   ```
3. **通道自动枚举与四大基组间机器精度（$<10^{-14}$）幺正变换**：
   ```fortran
   type(field_channel_t), allocatable :: ch_unc(:), ch_spin(:)
   real(dp), allocatable :: U_spin_unc(:, :)
   integer :: n_ch

   ! 给定守恒总磁量子数 2*M_tot 与分波 l_max (如 s-波 l=0)
   call build_field_collision_channels(rb87, rb87, BASIS_UNCOUPLED, two_Mtot=2, l_max=0, &
                                       channels=ch_unc, n_channels=n_ch)
   call build_field_collision_channels(rb87, rb87, BASIS_TOTAL_SPIN, two_Mtot=2, l_max=0, &
                                       channels=ch_spin, n_channels=n_ch)

   ! 构建 U_{spin <- unc} 变换矩阵 (满足 U * U^T = I)
   allocate(U_spin_unc(n_ch, n_ch))
   call calc_basis_transform_matrix(rb87, rb87, ch_unc, ch_spin, n_ch, &
                                    BASIS_UNCOUPLED, BASIS_TOTAL_SPIN, 0.0_dp, U_spin_unc)
   ```
4. **外场多通道相互作用势矩阵与磁 Feshbach 共振色散扫描**：
   ```fortran
   type(field_feshbach_result_t) :: fb_res

   ! 自动组装 V_{ij}(r; B) 并通过 Johnson 矩阵对数导数求解多通道密耦，扫描磁 Feshbach 共振
   call calc_magnetic_feshbach_resonance_scan( &
       r_grid, v_singlet, v_triplet, rb87, rb87, &
       basis_type=BASIS_UNCOUPLED, two_Mtot=2, l_max=0, &
       b_min=50.0_dp, b_max=110.0_dp, n_b=61, incident_energy=1.0e-9_dp, &
       incident_channel=1, res=fb_res)

   print *, "检测到磁 Feshbach 共振极点 B_0 (Gauss):", fb_res%b_res_pole
   print *, "共振宽度 Delta_B (Gauss):", fb_res%delta_b
   print *, "背景散射长度 a_bg (a0):", fb_res%a_bg
   ```

> [!TIP]
> 完整的多扇区分段网格外场超冷散射与磁 Feshbach 共振拟合端到端工程示例，详见：
> 👉 [`examples/ex08_ultracold_feshbach_segmented.f90`](examples/ex08_ultracold_feshbach_segmented.f90)
> 对应理论参考文献与公式映射，详见：
> 👉 [`LITERATURE.md`](LITERATURE.md) 第 4 节与第 5 节。

---

### 3.15 各向异性偶极散射与超冷自旋弛豫配置 (`mod_dipolar_scattering`)
面向超冷偶极气体（如磁性原子 $^{52}\text{Cr}, ^{164}\text{Dy}, ^{168}\text{Er}$ 或偶极玻色子 $^{87}\text{Rb}$）与极性双原子分子（如 $^{40}\text{K}^{87}\text{Rb}, ^{23}\text{Na}^{40}\text{K}$）的长程各向异性相互作用：
1. **各向异性磁偶极 (MDDI) 算符矩阵元与两体自旋张量**：
   ```fortran
   real(dp) :: c2q_val, spin_elem, c_dd
   ! 1. 空间秩-2 球谐张量矩阵元 <l1, m1 | C_{2, q} | l2, m2>
   c2q_val = c2q_orbital_matrix_element(l1=0, m1=0, l2=2, m2=0, q=0)
   ! 2. 秩-2 两体自旋张量矩阵元 <S, Ms | [s1 x s2]^{(2)}_q | S', Ms'>
   spin_elem = spin_tensor_coupled_matrix_element(s1=1.0_dp, s2=1.0_dp, &
                                                   s_tot=1, ms_tot=0, s_prime=1, ms_prime=0, q=0)
   ! 3. 磁偶极特征耦合系数 C_dd (a.u.)
   c_dd = calc_mddi_coupling_strength(mu1_bohr=1.0_dp, mu2_bohr=1.0_dp)
   ```
2. **超冷原子磁阱自旋弛豫截面 $\sigma_{\text{rel}}$ 与热速率 $K_{\text{rel}}(T)$**：
   ```fortran
   real(dp) :: sigma_rel, k_rel
   ! 计算能量 E 下两体自旋弛豫散射截面
   call calc_dipolar_relaxation_cross_section(energy=1.0e-9_dp, mass=87.0_dp*AMU2AU, &
                                              mu_mag=1.0_dp, delta_m=1, sigma_rel=sigma_rel)
   ! 计算微开尔文 (T = 1.0 uK) 玻尔兹曼热平衡系综平均弛豫速率 (cm^3/s)
   call calc_dipolar_relaxation_thermal_rate(temp_kelvin=1.0e-6_dp, mass=87.0_dp*AMU2AU, &
                                             mu_mag=1.0_dp, delta_m=1, k_rel=k_rel)
   ```
3. **极性分子直流 Stark 诱导电偶极矩与特征偶极长度**：
   ```fortran
   type(polar_molecule_t) :: krb
   real(dp) :: d_ind, a_d
   ! 初始化 KRb 极性分子 (转动常数 B_e = 1.114 GHz, 永久偶极矩 d_0 = 0.574 Debye)
   krb%rot_constant_ghz = 1.114_dp
   krb%dipole_moment_debye = 0.574_dp
   krb%mass_amu = 127.0_dp

   ! 施加 DC 外电场 (10.0 kV/cm) 求解诱导电偶极矩
   call calc_stark_induced_dipole(krb, e_field_dc_kv_cm=10.0_dp, d_ind_debye=d_ind)
   ! 计算电偶极特征散射长度 a_d (a.u.)
   call calc_electric_dipolar_length(krb%mass_amu * AMU2AU, d_ind, a_d)
   ```
4. **偶极多通道耦合势矩阵组装**：
   ```fortran
   type(dipolar_channel_t), allocatable :: channels(:)
   real(dp), allocatable :: v_dd(:, :)
   ! 构建自旋 S=1、分波截断 l_max=2 的耦合通道基
   call build_dipolar_channel_basis(s_tot=1, l_max=2, channels=channels, n_channels=n_ch)
   ! 组装距离 r 处的各向异性耦合势能矩阵
   call calc_dipolar_potential_matrix(r=50.0_dp, channels=channels, n_channels=n_ch, &
                                      mu1_bohr=1.0_dp, mu2_bohr=1.0_dp, v_mat=v_dd)
   ```

---

### 3.16 超冷光缔合谱学与分子生成配置 (`mod_photoassociation`)
面向超冷激光受激自由-束缚态光缔合 (Photoassociation, PA) 与 STIRAP 超冷分子制备：
1. **自由-束缚 Franck-Condon 重叠积分与跃迁态密度**：
   ```fortran
   real(dp) :: overlap, f_fb
   ! 散射能量本征波函数与激发态分子束缚态波函数积分
   call calc_free_bound_fc_overlap(r_grid, psi_free_norm, psi_bound, overlap)
   ! 计算态密度 f_FB(E) = |<psi_E | psi_v>|^2
   call calc_free_bound_fc_density(energy=1.0e-9_dp, overlap=overlap, f_fb=f_fb)
   ```
2. **激光功率依赖的受激线宽与光缔合散射截面**：
   ```fortran
   type(pa_transition_t) :: trans
   real(dp) :: gamma_stim, sigma_pa
   trans%laser_intensity_w_cm2 = 100.0_dp      ! 激光强度 100 W/cm^2
   trans%trans_dipole_debye = 2.5_dp           ! 跃迁电子偶极矩 2.5 Debye
   trans%gamma_nat_mhz = 6.0_dp                ! 分子自发辐射天然线宽 6 MHz
   trans%fc_overlap = overlap

   ! 受激跃迁线宽 hbar * Gamma_stim
   call calc_pa_stimulated_linewidth(trans%laser_intensity_w_cm2, trans%trans_dipole_debye, &
                                     trans%fc_overlap, gamma_stim)
   ! 求解单能量入射散射吸收截面 (失谐 Delta = 0)
   call calc_pa_cross_section(trans, energy=1.0e-9_dp, detuning_au=0.0_dp, cross_sec_au=sigma_pa)
   ```
3. **热平衡系综平均光缔合速率与光谱扫描**：
   ```fortran
   real(dp) :: k_pa, detuning_grid(200), rates(200), peak_det
   ! 计算 T = 50 uK 下的热平均光缔合速率系数 K_PA(T, Delta)
   call calc_pa_thermal_rate_coefficient(trans, temp_kelvin=50.0e-6_dp, &
                                         detuning_au=0.0_dp, k_pa=k_pa)
   ! 激光频率失谐连续扫描生成完整光缔合吸收能谱
   call calc_pa_detuning_scan(trans, temp_kelvin=50.0e-6_dp, detunings=detuning_grid, &
                              n_pts=200, rates=rates, peak_detuning=peak_det)
   ```
4. **双光子 Raman / STIRAP 缔合基态分子有效耦合**：
   ```fortran
   real(dp) :: omega_eff
   ! 给定泵浦光与斯托克斯光拉比频率 Omega_1, Omega_2 及激发态中间失谐 Delta_1
   call calc_twophoton_raman_coupling(omega1=10.0_dp, omega2=15.0_dp, delta1=100.0_dp, &
                                      omega_eff=omega_eff)
   ```

> [!TIP]
> 完整的各向异性偶极自旋弛豫与光缔合分子生成前沿工程算例，详见：
> 👉 [`examples/ex09_dipolar_relaxation_scattering.f90`](examples/ex09_dipolar_relaxation_scattering.f90)
> 👉 [`examples/ex10_photoassociation_spectroscopy.f90`](examples/ex10_photoassociation_spectroscopy.f90)
> 对应理论推导与学术文献全典，详见：
> 👉 [`LITERATURE.md`](LITERATURE.md) 第 9 节与第 10 节。

---

### 3.17 超冷少体物理与 Efimov 三体复合配置 (`mod_three_body_recombination`)
面向超冷原子少体物理、Efimov 普适三聚体能级与三体复合损耗速率：
1. **Efimov 超径向超越方程根与离散标度因子**：
   ```fortran
   use mod_three_body_recombination
   real(dp) :: s0, lambda_scale
   ! 求解全同玻色子 Efimov 超越方程普遍根 (s_0 ~ 1.00624)
   s0 = solve_efimov_s0_identical_bosons()
   ! 离散标度因子 lambda = exp(pi / s_0) ~ 22.694
   lambda_scale = calc_efimov_scale_factor(s0)
   ```
2. **Braaten-Hammer 普适三体复合速率 $K_3(a)$**：
   ```fortran
   type(efimov_param_t) :: param
   type(three_body_loss_t) :: loss
   integer :: stat
   param%s0 = s0
   param%scale_factor = lambda_scale
   param%a_star = 200.0_dp   ! a_* 三体参数
   param%a_minus = -100.0_dp ! a_- 三聚体共振极点
   param%eta_star = 0.06_dp  ! 非弹性耗散因子

   ! 计算正散射长度 a > 0 下的复合损耗（呈现 Efimov 干涉相消极小值）
   call calc_three_body_recombination_a_positive(a_scat=500.0_dp, mass_atom=87.0_dp*AMU2AU, &
                                                 param=param, res=loss, stat=stat)
   print *, "K_3(a>0) [cm^6/s] = ", loss%k3_si
   ```
3. **幺正极限有限温度饱和幂律 $K_3 \propto T^{-2}$**：
   ```fortran
   real(dp) :: k3_t_au, k3_t_si
   ! 计算 T = 1.0 uK 下幺正极限三体复合损耗速率
   call calc_unitary_three_body_loss_temperature(temp_kelvin=1.0e-6_dp, mass_atom=87.0_dp*AMU2AU, &
                                                 eta_star=0.06_dp, k3_au=k3_t_au, k3_si=k3_t_si, stat=stat)
   ```

---

### 3.18 低维光晶格受限量子散射与 CIR 配置 (`mod_confined_scattering`)
面向光晶格一维谐振波导与二维平面囚禁中的低维超冷量子散射与约束诱导共振：
1. **一维光波导初始化与 Olshanii CIR 极点**：
   ```fortran
   use mod_confined_scattering
   type(waveguide_1d_t) :: wg
   type(cir_result_t)   :: cir
   integer :: stat
   ! 约束频率 omega_perp = 2*pi * 20 kHz，折合质量 mu = 43.5 amu
   call init_waveguide_1d(omega_trans_au=7.61e-12_dp, reduced_mass_au=43.5_dp*AMU2AU, &
                          wg=wg, stat=stat)
   print *, "谐振子长度 a_perp (a0) = ", wg%a_perp_au
   print *, "Olshanii CIR 极点 a_CIR = ", wg%a_cir_au
   ```
2. **重整化 1D 相互作用强度 $g_{\text{1D}}$ 与结合能 $E_b$**：
   ```fortran
   ! 计算 3D 散射长度 a_s 处的 1D 有效耦合常数与 Tonks 态判据
   call calc_olshanii_cir_parameters(wg, a_scat_3d_au=500.0_dp, linear_density_au=5.0e-5_dp, &
                                     res=cir, stat=stat)
   ! 求解受限诱导双原子分子结合能 E_b
   e_b = calc_confined_dimer_binding_energy(wg, a_scat_3d_au=500.0_dp)
   ```

---

### 3.19 自电离体系与 Fano 共振 / 复坐标旋转法配置 (`mod_autoionization_fano`)
面向原子分子双激发态自电离、组态相互作用与复能级共振寿命：
1. **Fano 不对称吸收线型与抗共振零点**：
   ```fortran
   use mod_autoionization_fano
   type(fano_profile_t) :: fano
   real(dp) :: sigma, tau_au, tau_fs
   integer :: stat
   fano%e_resonance_au = 2.22_dp
   fano%gamma_width_au = 0.00137_dp
   fano%q_parameter    = -2.80_dp
   fano%sigma_0_au     = 1.0_dp
   call calc_autoionization_lifetime(fano%gamma_width_au, tau_au, tau_fs, stat)
   sigma = calc_fano_profile(fano, energy_au=2.22_dp)
   ```
2. **复坐标旋转法 (CCR) 提取复本征能量 $E_R - i\Gamma/2$**：
   ```fortran
   type(ccr_resonance_t) :: ccr
   ! 坐标复旋转角度 theta = 0.30 rad 求解准束缚共振态
   call solve_ccr_resonance_model(e_bound_0=2.22_dp, e_cont_0=2.20_dp, v_coupl=0.015_dp, &
                                  theta_rad=0.30_dp, res=ccr, stat=stat)
   print *, "共振能量 E_R = ", ccr%e_r_au, " 衰变宽度 Gamma = ", ccr%gamma_au
   ```

---

### 3.20 交叉静电磁场分子量子动力学配置 (`mod_crossed_field_scattering`)
面向开壳层极性顺磁分子在非共线 $\mathbf{E} \times \mathbf{B}$ 交叉外场中的态混合与取向调控：
1. **交叉场参数初始化与哈密顿量对角化**：
   ```fortran
   use mod_crossed_field_scattering
   type(crossed_field_config_t) :: cfg
   type(crossed_field_state_t)  :: state
   integer :: stat
   ! 电场 12 kV/cm，磁场 1000 G，夹角 theta_EB = 45 deg，转动常数 1.114 GHz，偶极 0.574 D
   call init_crossed_field_config(e_field_kv_cm=12.0_dp, b_field_gauss=1000.0_dp, &
                                  theta_eb_deg=45.0_dp, rot_ghz=1.114_dp, &
                                  dipole_d=0.574_dp, j_max=2, cfg=cfg, stat=stat)
   call solve_crossed_field_eigenstates(cfg, state, stat)
   call calc_crossed_field_observables(cfg, state)
   print *, "基态 Stark 定向度 <cos theta> = ", state%orientation(1)
   ```

---

### 3.21 三原子反应散射 Jacobi 几何与 LEPS 势能面配置 (`mod_triatomic_geometry`)
面向三原子反应动力学散射网格、LEPS 反应势能面与锥形交叉 Berry 几何相位：
1. **Jacobi 反应坐标与核间距双向转换**：
   ```fortran
   use mod_triatomic_geometry
   type(jacobi_coord_t)      :: jac_in, jac_out
   type(internuclear_dist_t) :: dist
   ! Jacobi (r, R, gamma) -> 核间距 (r12, r23, r31)
   call jacobi_to_internuclear(jac_in, dist)
   ! 核间距 -> Jacobi 坐标
   call internuclear_to_jacobi(dist, m_a, m_b, m_c, jac_out)
   ```
2. **基准 H3 LEPS 反应势能面与鞍点活化能垒**：
   ```fortran
   type(leps_param_t) :: leps_h3
   call init_default_h3_leps(leps_h3)
   v_pot = calc_leps_potential(dist, leps_h3)
   ```
3. **锥形交叉 (CI) 绝热分裂与拓扑 Berry 几何相位**：
   ```fortran
   type(conical_intersection_t) :: ci
   ci%x_ci = 0.0_dp; ci%y_ci = 0.0_dp; ci%kappa_tuning = 0.5_dp; ci%lambda_coupl = 0.5_dp
   ! 环绕锥形交叉闭合回路积分获得拓扑 Berry 相位 Phi_B = pi
   phi_b = calc_berry_phase_around_ci(ci, radius=0.20_dp, n_steps=1000)
   ```

---

### 3.22 超冷旋量玻色爱因斯坦凝聚自旋动力学配置 (`mod_spinor_bec`)
面向 $F=1$ 旋量凝聚体（$^{87}\text{Rb}$ 铁磁相 / $^{23}\text{Na}$ 极性反铁磁相）相干自旋混合与相变：
1. **旋量凝聚体参数预设与分类**：
   ```fortran
   use mod_spinor_bec
   type(spinor_param_t) :: param
   integer :: stat
   ! 装载 87Rb 凝聚体参数 (磁场 B = 0.25 G, 密度 n = 1e14 cm^-3)
   call init_spinor_preset("87Rb", b_field_gauss=0.25_dp, density_cm3=1.0e14_dp, param=param, stat=stat)
   print *, "自旋交换常数 c_2 = ", param%c2_au, " 是否铁磁: ", param%is_ferromagnet
   ```
2. **RK4 保全几率与保磁化强度相干自旋振荡演化**：
   ```fortran
   type(spinor_state_t) :: state
   ! 单步步长 dt_au 推进单模近似 (SMA) 旋量动力学方程
   call propagate_spinor_sma_rk4(param, dt_au=100.0_dp, state=state)
   print *, "m=0 分量布居: ", state%pop_0, " 磁化强度: ", state%magnetization_mz
   ```

---

### 3.23 三原子超球面反应动力学与热速率常数配置 (`mod_hyperspherical_reactive`)
面向多原子化学反应动力学 $A + BC \to AB + C$ 的超球面坐标几何与量子隧穿速率常数：
1. **反应质量标度与反应偏角 $\beta_{skew}$**：
   ```fortran
   use mod_hyperspherical_reactive
   type(reaction_mass_t) :: rmass
   ! 初始化 H + H2 质量运动学参数
   call init_reaction_mass(mass_a_amu=1.0078_dp, mass_b_amu=1.0078_dp, mass_c_amu=1.0078_dp, rmass=rmass)
   print *, "Delves 标度因子 d = ", rmass%scale_factor_d, " 偏角 beta = ", rmass%skew_angle_deg, " deg"
   ```
2. **Eckart 鞍点量子隧穿与累积反应几率 $N(E)$**：
   ```fortran
   type(transition_state_t) :: ts
   real(dp) :: prob_t, n_e
   ts%v_barrier_au = 0.425_dp * EV2AU; ts%omega_im_au = 1500.0_dp * CM2AU
   ts%omega_bend_au = 900.0_dp * CM2AU; ts%omega_symm_au = 2050.0_dp * CM2AU; ts%n_trans_states = 5
   ! 计算总能量 E 下的 Eckart 势垒隧穿传递几率与多通道累积反应几率
   prob_t = calc_eckart_transmission(energy_au=0.45_dp*EV2AU, v_barrier_au=ts%v_barrier_au, omega_im_au=ts%omega_im_au)
   n_e    = calc_cumulative_reaction_probability(ts, energy_au=0.45_dp*EV2AU)
   ```
3. **正则热反应速率常数 $k(T)$ 与 Wigner 隧穿修正**：
   ```fortran
   real(dp) :: k_exact, k_tst
   ! 室温 300 K 下基于 N(E) 严格玻尔兹曼积分的全量子速率与 Wigner 修正 TST 速率 (cm^3/s)
   k_exact = calc_canonical_rate_constant(ts, rmass, temp_kelvin=300.0_dp, n_e_steps=500)
   k_tst   = calc_tst_wigner_rate(ts%v_barrier_au, ts%omega_im_au, temp_kelvin=300.0_dp, prefactor=1.0e-10_dp)
   ```

---

### 3.24 超冷偶极量子液滴与李-黄-杨量子涨落配置 (`mod_dipolar_droplets_lhy`)
面向磁性稀薄玻色气体（$^{162}\text{Dy}, ^{166}\text{Er}$）中自束缚偶极量子液滴的宏观平顶密度与相平衡：
1. **偶极长度与 Pelster-Lima 涨落积分 $Q_5(\epsilon_{dd})$**：
   ```fortran
   use mod_dipolar_droplets_lhy
   type(dipolar_droplet_param_t) :: param
   integer :: stat
   ! 装载 162Dy (磁矩 10 mu_B) 在散射长度 a_s = 70 a0 下的特征参数
   call init_dipolar_droplet_param("162Dy", a_scat_bohr=70.0_dp, param=param, stat=stat)
   print *, "偶极长度 a_dd = ", param%a_dd_au, " a0, 偶极强度 eps_dd = ", param%epsilon_dd
   print *, "量子涨落增强因子 Q_5 = ", param%q5_factor
   ```
2. **自由空间自束缚平衡密度 $n_0$ 与负化学势 $\mu(n_0) < 0$**：
   ```fortran
   real(dp) :: n0_au, n0_cm3, mu_eq, n_crit
   n0_au  = calc_equilibrium_droplet_density(param)
   n0_cm3 = n0_au / ((5.29177210903e-9_dp)**3)
   mu_eq  = calc_droplet_chemical_potential(param, n0_au)
   n_crit = calc_critical_atom_number(param)
   print *, "平顶平衡密度 n_0 = ", n0_cm3, " cm^-3, 化学势 mu = ", mu_eq, " a.u., 临界原子数 = ", n_crit
   ```

---

### 3.25 强场非顺序双电离与电子重碰撞相关动量谱配置 (`mod_strong_field_nsdi`)
面向强红外激光脉冲（如 800 nm, $10^{14}-10^{15} \text{ W/cm}^2$）驱动稀有气体原子（$\text{He}, \text{Ar}, \text{Ne}$）强场双电离：
1. **激光场与靶原子初始化**：
   ```fortran
   use mod_strong_field_nsdi
   type(nsdi_laser_t)  :: laser
   type(nsdi_target_t) :: target
   integer :: stat
   call init_nsdi_laser(wavelength_nm=800.0_dp, intensity_w_cm2=2.5e14_dp, laser=laser, stat=stat)
   call init_nsdi_target("Ar", target=target, stat=stat)
   print *, "有质动力能 Up = ", laser%up_au * AU2EV, " eV, 3.17 Up 截断能 = ", 3.173_dp * laser%up_au * AU2EV, " eV"
   ```
2. **经典回碰轨道与 2D 双电子平行动量分布 $P(p_{z1}, p_{z2})$**：
   ```fortran
   real(dp) :: grid_p(31), dist_2d(31, 31), corr_coeff
   ! 数值积分生成 COLTRIMS 实验可测的双电子纵向动量谱并提取关联系数
   call calc_nsdi_2d_momentum_dist(laser, target, n_pts=31, p_max=2.5_dp, &
                                   grid_p=grid_p, dist_2d=dist_2d, corr_coeff=corr_coeff)
   print *, "一三象限相关系数 C_corr = ", corr_coeff, " (> 0 表明同向关联出射)"
   ```
3. **双电离产率光强依赖曲线与非顺序“膝盖结构”**：
   ```fortran
   real(dp) :: intensities(10), y_nsdi(10), y_sdi(10)
   call calc_double_ion_yield_curve(800.0_dp, target, 10, 1.5e14_dp, 8.0e14_dp, intensities, y_nsdi, y_sdi)
   ```

---

### 3.26 磁与光 Feshbach 共振与弱束缚分子态配置 (`mod_feshbach_bound_states`)
面向超冷原子磁场 Feshbach 共振（MFR）分子态能谱与激光驱动光 Feshbach 共振（OFR）：
1. **磁 Feshbach 共振分类与两通道分子结合能 $E_b(B)$**：
   ```fortran
   use mod_feshbach_bound_states
   type(mfr_param_t) :: li6_mfr
   real(dp) :: a_b, eb_coup, z_closed
   integer :: stat
   ! 装载 6Li 832 G 极宽共振参数
   call init_mfr_preset("6Li", li6_mfr, stat)
   print *, "共振极点 B0 = ", li6_mfr%b0_gauss, " G, 强度参数 s_res = ", li6_mfr%s_res, " (>> 1 为宽共振)"
   ! 计算 B = 800 G 下的散射长度、耦合通道分子结合能与闭通道占比 Z(B)
   a_b      = calc_mfr_scattering_length(li6_mfr, 800.0_dp)
   eb_coup  = calc_mfr_bound_energy_coupled(li6_mfr, 800.0_dp)
   z_closed = calc_mfr_closed_channel_fraction(li6_mfr, 800.0_dp)
   print *, "散射长度 a(B) = ", a_b, " a0, 结合能 = ", eb_coup * AU2EV, " eV, 闭通道权重 Z = ", z_closed
   ```
2. **光 Feshbach 共振 (OFR) 复散射长度与光致损耗率 $K_2$**：
   ```fortran
   type(ofr_param_t) :: ofr
   real(dp) :: a_re, a_im, k2_loss
   call init_ofr_param("87Rb", mass_amu=86.91_dp, a_bg_bohr=100.0_dp, gamma_hz=1.0e7_dp, &
                       l_opt_bohr=50.0_dp, ofr=ofr, stat=stat)
    ! 计算激光失谐 Delta_L = +5.0 MHz 下的复散射长度与双体非弹性损失常数 (cm^3/s)
    call calc_ofr_complex_scattering_length(ofr, delta_hz=5.0e6_dp, a_real_au=a_re, a_imag_au=a_im)
    k2_loss = calc_ofr_inelastic_loss_rate(ofr, delta_hz=5.0e6_dp)
    print *, "调谐后散射长度 Re(a) = ", a_re, " a0, 双体损失率 K_2 = ", k2_loss, " cm^3/s"
    ```

---

### 3.27 阿秒瞬态吸收光谱与光诱导态自电离干涉 (`mod_attosecond_transient_absorption`)
面向 XUV 阿秒单脉冲与强 NIR 飞秒探针激光场操纵自电离态吸收动力学：
1. **初始化氦原子自电离态与光诱导态 (LIS)**：
   ```fortran
   use mod_attosecond_transient_absorption
   type(atas_state_t) :: he_state
   real(dp) :: e_lis, t_beat
   integer :: stat
   call init_atas_helium_benchmark(he_state, stat)
   e_lis  = calc_light_induced_state_energy(he_state%energy_ev, e_dark_ev=58.60_dp, &
                                            omega_nir_ev=1.55_dp, rabi_ev=0.15_dp)
   t_beat = calc_quantum_beat_period_fs(abs(he_state%energy_ev - 58.60_dp))
   print *, "LIS 态能量 = ", e_lis, " eV, 量子拍频周期 = ", t_beat, " fs"
   ```
2. **全量计算二维时延瞬态吸收谱 $\Delta\text{OD}(\omega, \tau)$**：
   ```fortran
   real(dp) :: e_grid(31), tau_grid(25), spec_2d(31, 25)
   call calc_atas_spectrum(he_state, nir_intensity_w_cm2=2.0e12_dp, nir_wavelength_nm=800.0_dp, &
                           n_energy=31, e_min_ev=59.5_dp, e_max_ev=60.8_dp, &
                           n_delay=25, tau_min_fs=-30.0_dp, tau_max_fs=30.0_dp, &
                           e_grid_ev=e_grid, tau_grid_fs=tau_grid, spec_2d=spec_2d)
   ```

---

### 3.28 双色反向圆偏振场与分子光电子圆二色性 PECD (`mod_bicircular_pecd`)
面向 $\omega + 2\omega$ 双色旋转场与手性四面体势单/多光子光电子角分布不对称性：
1. **双色圆偏振场合成与离散动力学对称性**：
   ```fortran
   use mod_bicircular_pecd
   type(bicircular_field_t) :: field
   integer :: n_fold
   call init_bicircular_field(field, omega1_au=0.057_dp, r_freq=2.0_dp, &
                              i1_wcm2=1.0e14_dp, i2_wcm2=5.0e13_dp, &
                              h1=1, h2=-1, phi1=0.0_dp, phi2=0.0_dp, &
                              fwhm_fs=25.0_dp, envelope_type=1)
   n_fold = calc_dynamical_symmetry_fold(field%h1, field%h2, freq_ratio=2)
   print *, "合成场离散旋转对称度: C_", n_fold  ! 反向旋转输出 C_3 (三叶草形)
   ```
2. **手性分子对映体初始化与 PECD 前后发射不对称度**：
   ```fortran
   type(chiral_tetrahedral_molecule_t) :: mol_r, mol_s
   real(dp) :: chi_r, chi_s, b1_r, g_pecd_r
   call init_chiral_tetrahedral_molecule(mol_r, "R")
   call init_chiral_tetrahedral_molecule(mol_s, "S")
   chi_r = calc_chirality_measure(mol_r)
   chi_s = calc_chirality_measure(mol_s)  ! 满足严格反号 chi_s = -chi_r
   b1_r  = calc_chiral_beta1_model(mol_r, energy_ev=5.0_dp, photon_energy_ev=10.0_dp)
   g_pecd_r = calc_forward_backward_asymmetry(b1_r)
   print *, "R 对映体手性不对称度 G_PECD = ", g_pecd_r * 100.0_dp, " %"
   ```

---

### 3.29 超冷极性分子反应动力学与微波/静电偶极遮蔽 (`mod_ultracold_reaction_shielding`)
面向超冷极性分子（KRb, NaRb）微波蓝失谐免交叉排斥势垒与非弹性损失抑制：
1. **分子预设与屏蔽场配置**：
   ```fortran
   use mod_ultracold_reaction_shielding
   type(ultracold_molecule_t) :: krb
   type(shielding_config_t)   :: cfg
   real(dp) :: r_bar, v_bar_k
   call init_ultracold_molecule_preset(krb, "KRb")
   call init_shielding_config(cfg, method=1, detuning_mhz=15.0_dp, rabi_mhz=5.0_dp, &
                             e_field_kv_cm=0.0_dp, y_loss=1.0_dp)
   call calc_shielding_barrier_height(krb, cfg, r_bar, v_bar_k)
   print *, "遮蔽势垒半径 R_bar = ", r_bar, " a0, 势垒高度 = ", v_bar_k * 1.0e6_dp, " uK"
   ```
2. **WKB 隧穿抑制与蒸发冷却优良因子 $\gamma$**：
   ```fortran
   real(dp) :: k2_el, k2_inel, gamma_ratio
   call calc_shielded_scattering_rates(krb, cfg, temp_uk=0.5_dp, &
                                      k2_el_cm3s=k2_el, k2_inel_cm3s=k2_inel, &
                                      gamma_ratio=gamma_ratio)
   print *, "弹性散射率 K2_el = ", k2_el, " cm^3/s, 非弹性损失率 K2_inel = ", k2_inel, " cm^3/s"
   print *, "冷却因子 gamma = ", gamma_ratio, " (>> 100 满足玻色/费米简并蒸发冷却准则)"
   ```

---

### 3.30 里德堡原子阻塞、PXP 约束模型与量子多体疤痕 (`mod_rydberg_blockade`)
面向里德堡原子量子模拟器二原子阻塞与 1D 链多体疤痕相干振荡：
1. **里德堡态 $C_6 \propto n^{11}$ 标度律与阻塞半径 $R_b$**：
   ```fortran
   use mod_rydberg_blockade
   type(rydberg_atom_t) :: rb70
   real(dp) :: r_b
   call init_rydberg_atom(rb70, "87Rb", n_principal=70, l_orbital=0)
   r_b = calc_rydberg_blockade_radius(rb70, rabi_mhz=2.0_dp)
   print *, "87Rb 70S C6/h = ", rb70%c6_mhz_um6, " MHz*um^6, 阻塞半径 R_b = ", r_b, " um"
   ```
2. **双原子阻塞动力学与 1D 阵列量子多体疤痕**：
   ```fortran
    type(rydberg_array_config_t) :: chain
    real(dp) :: t_arr(100), z2_arr(100)
    ! 模拟 10 原子阵列从交错 Néel 态出发的 PXP 疤痕周期复苏
    call init_rydberg_array(chain, n_atoms=10, spacing_um=5.0_dp, &
                            rabi_mhz=2.0_dp, detuning_mhz=0.0_dp, boundary_cond=2)
    call calc_rydberg_scar_dynamics(chain, rb70, t_max_us=2.5_dp, n_steps=100, &
                                    t_arr=t_arr, z2_order_arr=z2_arr)
    print *, "Z2 交错序参量初始 = ", z2_arr(1), " 周期复苏幅值 = ", maxval(z2_arr(20:))
    ```

---

### 3.31 表面量子散射与选择性吸附共振 (`mod_surface_scattering`)
面向低能轻原子（如热 He 束）在晶体表面的相干弹性/非弹性衍射与吸附束缚态共振：
1. **表面晶格与 Morse 吸引阱束缚态初始化**：
   ```fortran
   use mod_surface_scattering
   type(surface_lattice_t) :: lif
   type(surface_potential_t) :: he_pot
   integer :: n_bound
   ! LiF(001): ax=ay=2.84 A, corrugation=0.06 A, M_sub=25.94 amu, Theta_D=730 K
   call init_surface_lattice(lif, ax_ang=2.84_dp, ay_ang=2.84_dp, &
                             zeta_x_ang=0.06_dp, zeta_y_ang=0.06_dp, &
                             m_sub_amu=25.94_dp, debye_temp_k=730.0_dp)
   call init_surface_potential_morse(he_pot, well_depth_mev=7.5_dp, &
                                     range_inv_ang=1.1_dp, mass_amu=4.0026_dp, &
                                     n_bound=n_bound)
   ```
2. **2D 衍射通道与硬波纹表面 (HCS) 程函几率**：
   ```fortran
   type(diffraction_beam_t) :: channels(25)
   integer :: n_ch
   call calc_hcs_diffraction_probabilities(lif, mass_amu=4.0026_dp, energy_ev=0.020_dp, &
                                           theta_i_deg=40.0_dp, phi_i_deg=0.0_dp, &
                                           max_order=1, n_channels=n_ch, channels=channels)
   ```
3. **选择性吸附共振 (SAR) Fano 线型与声子 Debye-Waller 衰减**：
   ```fortran
   real(dp) :: de_mev, fano_ratio, dw_factor
   logical :: is_near
   call calc_selective_adsorption_resonance(lif, he_pot, mass_amu=4.0026_dp, &
                                            energy_ev=0.020_dp, theta_i_deg=55.0_dp, &
                                            phi_i_deg=0.0_dp, m_res=-1, n_res=0, &
                                            v_bound=0, is_near_res=is_near, &
                                            delta_e_mev=de_mev, fano_specular_ratio=fano_ratio)
   dw_factor = calc_surface_debye_waller(lif, mass_amu=4.0026_dp, k_iz_au=1.8_dp, &
                                         kz_g_au=1.8_dp, temp_k=300.0_dp)
   ```

---

### 3.32 气-固界面催化反应与 Eley-Rideal 提取机理 (`mod_surface_reaction_er`)
面向气相超热原子与表面化学吸附原子的直接碰撞提取反应动力学：
1. **预置反应体系初始化与放热量计算**：
   ```fortran
   use mod_surface_reaction_er
   type(er_reaction_system_t) :: er_sys
   call init_er_reaction_system(er_sys, "H+H/Cu(111)")
   print *, "反应放热量 Delta E = ", er_sys%delta_e_exo_ev, " eV"
   ```
2. **放热能量分配与产物振动布居反转**：
   ```fortran
   type(er_energy_partition_t) :: ep
   real(dp) :: p_vib(0:5)
   call calc_er_energy_partitioning(er_sys, e_incident_ev=0.15_dp, partition=ep)
   call calc_er_vibrational_populations(er_sys, e_incident_ev=0.15_dp, max_v=5, v_dist=p_vib)
   print *, "振动激发能量 E_vib = ", ep%e_vib_ev, " eV, P(v=2) = ", p_vib(2)
   ```
3. **反应截面与温度依赖催化速率常数**：
   ```fortran
   real(dp) :: sigma_er, k_rate
   sigma_er = calc_er_reaction_cross_section(er_sys, e_incident_ev=0.15_dp)
   k_rate   = calc_er_thermal_rate_constant(er_sys, temp_k=300.0_dp)
   print *, "截面 sigma = ", sigma_er, " A^2, 速率常数 k = ", k_rate, " cm^3/s"
   ```

---

### 3.33 金属表面非绝热动力学与电子摩擦耗散 (`mod_surface_electronic_friction`)
面向分子撞击金属表面时电子-空穴对激发与广义朗之万方程 (GLE) 动力学：
1. **金属基底初始化与指数空间摩擦系数分布**：
   ```fortran
   use mod_surface_electronic_friction
   type(metal_surface_t) :: au111
   real(dp) :: eta_val
   call init_metal_surface(au111, "Au(111)", temp_k=300.0_dp)
   eta_val = calc_electronic_friction_coeff(au111, z_bohr=2.0_dp)
   ```
2. **广义朗之万方程 (GLE) 轨迹步进与非绝热能损**：
   ```fortran
   type(scattering_loss_result_t) :: res
   real(dp) :: t_fs(200), z_ang(200), v_ms(200)
   call integrate_gle_scattering_trajectory(au111, mass_amu=30.006_dp, e_incident_ev=0.50_dp, &
                                           dt_fs=0.5_dp, n_steps=200, t_arr_fs=t_fs, &
                                           z_arr_ang=z_ang, v_arr_ms=v_ms, loss_res=res)
   print *, "单次碰撞电子-空穴对能损 = ", res%e_lost_ev, " eV, 最近距离 = ", res%z_turnaround_ang, " A"
   ```
3. **吸附分子振动弛豫寿命计算**：
   ```fortran
   real(dp) :: gamma_vib, tau_ps
   gamma_vib = calc_vibrational_relaxation_rate(au111, z_ang=1.058_dp, mass_amu=30.006_dp)
   tau_ps = 1.0_dp / gamma_vib
   print *, "表面振动耗散率 = ", gamma_vib, " ps^-1, 寿命 = ", tau_ps, " ps"
   ```

---

### 3.34 掠入射快原子表面量子衍射与彩虹散射 (`mod_grazing_fast_atom_diffraction`)
面向 keV 准直束快原子在低指数轴向沟道中的量子衍射与彩虹调制：
1. **GIFAD 实验参数与快慢自由度解耦**：
   ```fortran
   use mod_grazing_fast_atom_diffraction
   type(gifad_experiment_t) :: exp_cfg
   real(dp) :: e_perp, lambda_perp
   ! 1.0 keV He 束以 1.0 度掠射角入射 LiF(001) <110> 沟道 (ax=2.84 A, zeta=0.05 A)
   call init_gifad_experiment(exp_cfg, projectile="He", mass_amu=4.0026_dp, &
                             e_kev=1.0_dp, theta_deg=1.0_dp, &
                             ax_ang=2.84_dp, corrugation_ang=0.05_dp)
   call calc_gifad_transverse_kinematics(exp_cfg, e_perp, lambda_perp)
   print *, "横向垂直能量 E_perp = ", e_perp, " eV, 横向德布罗意波长 = ", lambda_perp, " A"
   ```
2. **经典表面彩虹偏转角与 1D 横向量子衍射谱**：
   ```fortran
   type(gifad_spectrum_t) :: spec
   real(dp) :: th_rainbow
   th_rainbow = calc_gifad_rainbow_angle(exp_cfg)
   call calc_gifad_diffraction_spectrum(exp_cfg, max_order=15, spec=spec)
   print *, "经典表面彩虹角 theta_R = ", th_rainbow, " deg, 开通道数 = ", spec%n_open_orders
   ```
3. **表面亚皮米波纹度逆向反演重构**：
   ```fortran
   real(dp) :: zeta_recon
   zeta_recon = calc_surface_corrugation_from_rainbow(exp_cfg%ax_channel_ang, th_rainbow)
   print *, "逆向重构波纹幅度 zeta = ", zeta_recon, " A (精度达亚皮米量级)"
   ```

---

### 3.35 冷离子-中性原子杂化散射与极化阱动力学 (`mod_ion_atom_scattering`)
面向 Paul 射频阱与光偶极阱杂化系统中的超冷带电离子与中性原子碰撞：
1. **杂化系统初始化与极化物理特征尺度**：
   ```fortran
   use mod_ion_atom_scattering
   type(ion_atom_system_t) :: sys
   integer :: stat
   ! Yb+ 离子 (174 amu, +1e) 与 6Li 原子 (6 amu, alpha = 164 a.u.)
   call init_ion_atom_system(sys, m_ion_amu=174.0_dp, m_atom_amu=6.0_dp, &
                             charge_ion=1.0_dp, alpha_atom_au=164.0_dp, stat=stat)
   print *, "特征极化长度 R* (a0) = ", sqrt(2.0_dp * sys%mu_au * sys%c4_au)
   ```
2. **Langevin 经典反应俘获截面与速率常数**：
   ```fortran
   real(dp) :: b_crit, sigma_langevin, k_langevin
   ! 计算碰撞能量 1.0 meV 下的临界碰撞参数与截面
   call calc_langevin_critical_impact_parameter(sys, e_coll_ev=1.0e-3_dp, b_crit_au=b_crit)
   sigma_langevin = calc_langevin_cross_section(sys, e_coll_ev=1.0e-3_dp)
   k_langevin = calc_langevin_rate_coefficient(sys)
   ```
3. **Paul 阱射频微运动非弹性致热与平衡极限温度**：
   ```fortran
   real(dp) :: dE_dt, t_limit
   ! 射频频率 2*pi * 2.0 MHz, Mathieu 稳定性参数 q = 0.25
   call calc_rf_micromotion_heating(sys, omega_rf_hz=1.256e7_dp, q_param=0.25_dp, &
                                    temp_ion_k=1.0e-3_dp, temp_atom_k=1.0e-6_dp, &
                                    heating_rate_k_per_s=dE_dt, temp_limit_k=t_limit)
   ```

---

### 3.36 最少开关表面跳跃与非绝热混合量子-经典动力学 (`mod_surface_hopping_fssh`)
面向光化学反应通道分支与非绝热势能面跃迁动力学：
1. **Tully 经典非绝热模型初始化与 NACV 导数耦合**：
   ```fortran
   use mod_surface_hopping_fssh
   type(tully_model_t) :: model
   real(dp) :: v_adia(2), d_nacv(2, 2)
   integer :: stat
   call init_tully_model(model, TULLY_SAC, stat)
   call calc_adiabatic_surface_and_nacv(model, r=0.0_dp, v_adia=v_adia, d_nacv=d_nacv)
   ```
2. **单条轨迹 Velocity Verlet 与 Tully 表面跳跃推进**：
   ```fortran
   type(fssh_trajectory_t) :: traj
   logical :: hopped
   call init_fssh_trajectory(traj, r_init=-5.0_dp, p_init=20.0_dp, init_state=1)
   call propagate_fssh_step(model, traj, dt=0.5_dp, jumped=hopped)
   ```
3. **蒙特卡洛系综与透射/反射分支比统计**：
   ```fortran
   real(dp) :: t1, t2, r1, r2
   call run_fssh_ensemble(model, n_trajectories=1000, r_init=-5.0_dp, p_init=20.0_dp, &
                          init_state=1, dt=0.5_dp, t_final=1000.0_dp, &
                          trans_1=t1, trans_2=t2, refl_1=r1, refl_2=r2)
   ```

---

### 3.37 强场分子定向、取向与超转子动力学 (`mod_molecular_alignment`)
面向短脉冲激光驱动的非绝热分子转动对齐与光学离心机超转子加速：
1. **转子分子初始化与转动复苏周期**：
   ```fortran
   use mod_molecular_alignment
   type(rotor_molecule_t) :: mol
   integer :: stat
   call init_rotor_molecule(mol, "N2", b_rot_cm1=1.9982_dp, delta_alpha_au=6.70_dp, &
                            dipole_debye=0.0_dp, d_dissoc_ev=9.76_dp, stat=stat)
   print *, "转动完全复苏周期 T_rev (ps) = ", mol%t_rev_ps
   ```
2. **飞秒激光诱导对齐与无场序参量 $\langle\cos^2\theta\rangle(t)$**：
   ```fortran
   real(dp) :: t_grid(100), cos2_trace(100)
   call simulate_laser_induced_alignment(mol, dt_fs=10.0_dp, t_max_ps=10.0_dp, &
                                        i_laser_wcm2=3.0e13_dp, duration_fs=100.0_dp, &
                                        t_grid_ps=t_grid, cos2_trace=cos2_trace, stat=stat)
   ```
3. **光学离心机角加速度驱动至超转子态与离心势垒破键**：
   ```fortran
   integer :: j_super
   logical :: is_broken
   real(dp) :: e_rot
   ! 线性啁啾加速度 alpha = 0.5 THz/ps, 脉宽 100 ps
   call calc_optical_centrifuge_kick(mol, alpha_chirp_thz_ps=0.5_dp, duration_ps=100.0_dp, &
                                     j_superrotor=j_super)
   call calc_superrotor_dissociation(mol, j_rot=j_super, is_dissociated=is_broken, e_rot_ev=e_rot)
   ```

---

### 3.38 超冷光晶格与玻色-哈伯德微观映射 (`mod_optical_lattice_hubbard`)
面向驻波光晶格中的周期量子输运与强关联玻色-哈伯德超流-Mott 绝缘体相变：
1. **光晶格系统初始化与 Mathieu 能带展开**：
   ```fortran
   use mod_optical_lattice_hubbard
   type(optical_lattice_t) :: latt
   real(dp) :: band_energies(4)
   integer :: stat
   ! 87Rb 原子，激光波长 1064 nm，晶格深度 s = 10 E_R
   call init_optical_lattice(latt, mass_amu=86.9_dp, lambda_nm=1064.0_dp, s_depth=10.0_dp, stat=stat)
   call calc_bloch_band_energies(latt, q_quasi=0.0_dp, n_bands=4, band_energies=band_energies)
   ```
2. **Bose-Hubbard 微观参数 $J$ 与 $U$ 映射**：
   ```fortran
   type(bose_hubbard_param_t) :: bh
   call calc_bose_hubbard_parameters(latt, a_s_nm=5.28_dp, bh=bh, stat=stat)
   print *, "跃迁矩阵元 J/h (Hz) = ", bh%tunneling_j_hz
   print *, "在位相互作用 U/h (Hz) = ", bh%onsite_u_hz
   print *, "超流-Mott 判定比值 U/J = ", bh%u_over_j
   ```
3. **引力布洛赫振荡周期与 Landau-Zener 带间隧穿几率**：
   ```fortran
   real(dp) :: t_bloch, nu_bloch, p_lz
   call calc_bloch_oscillation_dynamics(latt, force_grav_au=1.0e-25_dp, &
                                        t_bloch_s=t_bloch, nu_bloch_hz=nu_bloch, p_lz=p_lz)
   ```

---

### 3.39 多原子反应路径哈密顿量与变分过渡态理论 (`mod_reaction_path_hamiltonian`)
面向气相与凝聚相化学反应速率常数变分优化及量子隧穿穿透：
1. **反应路径参数化与垂直简正模频率**：
   ```fortran
   use mod_reaction_path_hamiltonian
   type(rph_path_t) :: path
   integer :: stat
   call init_rph_benchmark_reaction(path, reaction_type=1, stat=stat)
   ```
2. **变分过渡态理论 (CVT) 正则速率常数寻优**：
   ```fortran
   real(dp) :: s_bottleneck, rate_cvt
   call calc_cvt_rate_constant(path, temp_k=300.0_dp, s_opt=s_bottleneck, &
                               rate_cvt=rate_cvt, stat=stat)
   ```
3. **Eckart 势垒半经典量子隧穿修正系数**：
   ```fortran
   real(dp) :: kappa_tunnel
   call calc_eckart_tunneling_factor(path%barrier_height_ev, path%imag_freq_ts_cm1, &
                                     temp_k=300.0_dp, kappa_tunnel=kappa_tunnel)
   print *, "量子隧穿增强总速率 k(T) = ", kappa_tunnel * rate_cvt
   ```

---

### 3.40 相对论原子结构与径基狄拉克方程 (`mod_relativistic_atomic`)
面向重原子/高离化态相对论效应、自旋-轨道耦合与多极辐射矩阵元：
1. **Sommerfeld 狄拉克本征态与有效量子亏损**：
   ```fortran
   use mod_relativistic_atomic
   type(dirac_state_t) :: state_6s, state_6p12, state_6p32
   integer :: stat
   ! 铯 Cs (Z=55, z_ion=1.0, alpha_core=19.0 a.u., r_cut=2.0 a.u.)
   call solve_radial_dirac_eigenvalue(z_nuclear=55.0_dp, z_ion=1.0_dp, &
                                      alpha_core=19.0_dp, r_cut=2.0_dp, &
                                      n_princ=6, kappa=-1, state=state_6s, stat=stat)
   call solve_radial_dirac_eigenvalue(55.0_dp, 1.0_dp, 19.0_dp, 2.0_dp, 6, 1, state_6p12)
   call solve_radial_dirac_eigenvalue(55.0_dp, 1.0_dp, 19.0_dp, 2.0_dp, 6, -2, state_6p32)
   ```
2. **天然自旋-轨道耦合精细结构分裂 $\Delta E_{\text{FS}}$**：
   ```fortran
   real(dp) :: delta_ev, delta_cm1
   call calc_dirac_fine_structure_splitting(state_6p12, state_6p32, delta_ev, delta_cm1)
   print *, "6p 双重态精细分裂 = ", delta_cm1, " cm^-1"
   ```
3. **相对论电偶极 (E1) 振子强度 $f_{if}$**：
   ```fortran
   real(dp) :: f_d1, f_d2
   call calc_dirac_e1_matrix_element(state_6s, state_6p12, r_overlap_au=3.8_dp, osc_strength=f_d1)
   call calc_dirac_e1_matrix_element(state_6s, state_6p32, r_overlap_au=3.8_dp, osc_strength=f_d2)
   ```

---

### 3.41 共振非弹性 X 射线散射与内壳层光谱 (`mod_resonant_xray_scattering`)
面向同步辐射与 X 射线自由电子激光 (XFEL) 探测关联电子材料、低能集体激发与声子边带：
1. **RIXS 多能级系统配置**：
   ```fortran
   use mod_resonant_xray_scattering
   type(rixs_system_t) :: sys
   integer :: stat
   ! 配置 Cu L3 边 (931.5 eV, Gamma=0.35 eV) 与 dd 轨道激发终态 (1.80 eV)
   call init_rixs_system(sys, e_init=0.0_dp, e_inter=[931.5_dp], gamma_core=[0.35_dp], &
                         d_in=[1.0_dp], e_fin=[0.0_dp, 1.80_dp], gamma_fin=[0.05_dp, 0.08_dp], &
                         d_out=reshape([0.85_dp, 0.70_dp], [2, 1]), stat=stat)
   ```
2. **Kramers-Heisenberg 二阶散射截面与 2D RIXS Map**：
   ```fortran
   real(dp) :: sigma_rixs, rixs_2d(5, 5)
   sigma_rixs = calc_kramers_heisenberg_cross_section(sys, omega_in_ev=931.5_dp, omega_loss_ev=1.80_dp)
   call calc_rixs_2d_map(sys, 5, w_in_grid, 5, w_loss_grid, rixs_2d)
   ```
3. **电-声耦合 Huang-Rhys 振动 Franck-Condon 伴峰级数**：
   ```fortran
   real(dp) :: vib_intensity(0:4)
   ! 振动特征能量 70 meV, Huang-Rhys 常数 S = 0.40
   call calc_huang_rhys_vibrational_rixs(omega_0_ev=0.070_dp, s_factor=0.40_dp, &
                                        gamma_core_ev=0.35_dp, detuning_ev=0.0_dp, &
                                        n_max_loss=4, loss_intensity=vib_intensity)
   ```

---

### 3.42 亚稳态原子碰撞潘宁电离与缔合电离 (`mod_penning_associative_ionization`)
面向高激发态亚稳态稀有气体原子（如 $\text{He}^*(2^3S), \text{Ne}^*(^3P_{0,2})$）与靶原子/分子的化学电离、反应动力学与超冷量子气体寿命控制：
1. **潘宁与缔合电离光学势体系配置**：
   ```fortran
   use mod_penning_associative_ionization
   type(penning_system_t) :: sys
   integer :: stat
   ! 配置 He*(2^3S) + Ar 体系: Morse 势入口阱深 De=0.005 eV, 离子阱深 De=4.0 eV, 自电离宽度 A=5.0 eV, alpha=2.0 /A
   call init_penning_system(sys, &
                            v_star=morse_param_t(de_ev=0.005_dp, re_ang=4.5_dp, a_ang=1.2_dp), &
                            v_plus=morse_param_t(de_ev=4.0_dp, re_ang=2.4_dp, a_ang=1.8_dp), &
                            gamma_w=gamma_width_t(a_ev=5.0_dp, alpha_ang=2.0_dp), &
                            red_mass_amu=3.636_dp, r_max_au=30.0_dp, stat=stat)
   ```
2. **半经典光学势存活几率与分流截面**：
   ```fortran
   real(dp) :: sig_pi, sig_ai, sig_tot
   ! 求解 E_coll = 0.05 eV 下的潘宁电离截面与缔合离子生成截面
   call calc_penning_cross_sections(sys, e_coll_ev=0.05_dp, b_max_ang=6.0_dp, nb=1000, &
                                    sigma_pi=sig_pi, sigma_ai=sig_ai, sigma_tot=sig_tot)
   ```
3. **潘宁电离电子发射能谱 (PIES) 与超冷复散射长度**：
   ```fortran
   real(dp) :: pies(200), e_elec(200), k_el, k_loss
   call calc_pies_spectrum(sys, e_coll_ev=0.05_dp, n_pts=200, e_elec_grid=e_elec, pies_spec=pies)
   ! 超冷温度 1 uK 下基于复散射长度 a = alpha - i*beta 评估自旋极化寿命与损失速率
   call calc_ultracold_penning_rates(alpha_scat_m=1.0e-9_dp, beta_loss_m=1.0e-13_dp, &
                                     mass_kg=6.64e-27_dp, temp_uk=1.0_dp, &
                                     k_elastic=k_el, k_loss=k_loss)
   ```

> [!TIP]
> 包含全部前沿复杂体系的应用工程算例（11 至 36），详见：
> 👉 [`examples/ex11_three_body_efimov_recombination.f90`](examples/ex11_three_body_efimov_recombination.f90)
> 👉 [`examples/ex12_confined_cir_scattering.f90`](examples/ex12_confined_cir_scattering.f90)
> 👉 [`examples/ex13_autoionization_fano_resonance.f90`](examples/ex13_autoionization_fano_resonance.f90)
> 👉 [`examples/ex14_spinor_bec_dynamics.f90`](examples/ex14_spinor_bec_dynamics.f90)
> 👉 [`examples/ex15_crossed_field_stark_zeeman.f90`](examples/ex15_crossed_field_stark_zeeman.f90)
> 👉 [`examples/ex16_triatomic_reaction_berry_phase.f90`](examples/ex16_triatomic_reaction_berry_phase.f90)
> 👉 [`examples/ex17_hyperspherical_reaction_rates.f90`](examples/ex17_hyperspherical_reaction_rates.f90)
> 👉 [`examples/ex18_dipolar_quantum_droplets.f90`](examples/ex18_dipolar_quantum_droplets.f90)
> 👉 [`examples/ex19_strong_field_nsdi_recollision.f90`](examples/ex19_strong_field_nsdi_recollision.f90)
> 👉 [`examples/ex20_feshbach_molecular_bound_states.f90`](examples/ex20_feshbach_molecular_bound_states.f90)
> 👉 [`examples/ex21_attosecond_transient_absorption.f90`](examples/ex21_attosecond_transient_absorption.f90)
> 👉 [`examples/ex22_bicircular_pecd_chiral.f90`](examples/ex22_bicircular_pecd_chiral.f90)
> 👉 [`examples/ex23_ultracold_molecule_shielding.f90`](examples/ex23_ultracold_molecule_shielding.f90)
> 👉 [`examples/ex24_rydberg_blockade_dynamics.f90`](examples/ex24_rydberg_blockade_dynamics.f90)
> 👉 [`examples/ex25_surface_corrugated_diffraction.f90`](examples/ex25_surface_corrugated_diffraction.f90)
> 👉 [`examples/ex26_eley_rideal_surface_reaction.f90`](examples/ex26_eley_rideal_surface_reaction.f90)
> 👉 [`examples/ex27_surface_electronic_friction_gle.f90`](examples/ex27_surface_electronic_friction_gle.f90)
> 👉 [`examples/ex28_grazing_fast_atom_diffraction.f90`](examples/ex28_grazing_fast_atom_diffraction.f90)
> 👉 [`examples/ex29_cold_ion_atom_scattering.f90`](examples/ex29_cold_ion_atom_scattering.f90)
> 👉 [`examples/ex30_tully_surface_hopping.f90`](examples/ex30_tully_surface_hopping.f90)
> 👉 [`examples/ex31_molecular_alignment_revival.f90`](examples/ex31_molecular_alignment_revival.f90)
> 👉 [`examples/ex32_optical_lattice_bose_hubbard.f90`](examples/ex32_optical_lattice_bose_hubbard.f90)
> 👉 [`examples/ex33_rph_variational_transition_state.f90`](examples/ex33_rph_variational_transition_state.f90)
> 👉 [`examples/ex34_relativistic_dirac_cesium.f90`](examples/ex34_relativistic_dirac_cesium.f90)
> 👉 [`examples/ex35_rixs_core_level_spectroscopy.f90`](examples/ex35_rixs_core_level_spectroscopy.f90)
> 👉 [`examples/ex36_penning_associative_ionization.f90`](examples/ex36_penning_associative_ionization.f90)
> 对应理论文献详见：👉 [`LITERATURE.md`](LITERATURE.md) 第 11 节至第 36 节。

---

## 4. Python 伴侣库 `pygenmod` 配置与混合编程

`GeneralModule/python` 目录提供了一个符合 PEP 517/518 标准的 Python 纯粹伴侣分析库 `pygenmod`，用于快速完成参数预计算、波包与谱线生成。

### 4.1 本地可编辑模式安装
在您的 Python 环境（推荐 Conda 或 Python venv）中：
```bash
cd GeneralModule/python
pip install -e .
```
安装后即可在任意 Python 脚本中直接 `import pygenmod`。

### 4.2 数据交互规范（.dat 与无损二进制）
Fortran 程序推荐将计算所得的高维张量、含时序列或能量谱输出为标准空格/Tab 分隔的 `.dat` 文本文件：
```fortran
open(newunit=u, file="wavefunction.dat", status="replace", action="write")
write(u, '(A)') "# x(au)  Re(Psi)  Im(Psi)  ProbDensity"
do i = 1, n
    write(u, '(4ES16.8)') x(i), real(psi(i)), aimag(psi(i)), abs(psi(i))**2
end do
close(u)
```
Python 中使用 NumPy 一行代码加载：
```python
import numpy as np
data = np.loadtxt("wavefunction.dat")
x = data[:, 0]
prob = data[:, 3]
```

### 4.3 自动化发表级绘图管道
内置的 `visualizer.py` 提供自动美化的图表输出：
```python
from pygenmod import dvr_sinc_init, fgh_solve_bound_states, plot_wavefunctions

dvr = dvr_sinc_init(-6.0, 6.0, 256, mass=1.0)
v_pot = 0.5 * dvr.x**2
eig_vals, wavefuncs = fgh_solve_bound_states(dvr, v_pot)

# 生成 300 DPI 出版质量图片
plot_wavefunctions(dvr.x, v_pot, eig_vals, wavefuncs, n_states=4, filename="harmonic_oscillator.png")
```

---

## 5. AI Agent 自动化编程与提示词工程模板

在配合大模型（LLM）进行科学计算时，给 AI 提供高质量的上下文约束与输入规范，能保证 AI **一次性生成 100% 语法无错且物理准确的代码**。

### 5.1 给 AI 的上下文 System Prompt 模板
当您启动新的 AI 对话或配置 Cursor / Antigravity 规则时，可直接将以下提示词粘贴到提示框或 `.cursorrules` 中：

```markdown
你是一个量子动力学与超快强场物理的 Fortran/Python 专家。
在本项目中，必须遵循以下规范：
1. 始终使用 `use general_module` 引入算法接口，不要从头重复编写单位转换、勒让德多项式、三对角对角化或分裂算符 FFT。
2. 变量定义一律采用 `real(dp)`，其中 `dp` 来自 `general_module`。
3. 物理单位换算必须显式调用 `to_au(val, 'unit')` 或 `from_au(val, 'unit')`，严禁在代码中手写不透明的常数硬编码。
4. 激光脉冲使用 `type(pulse_config_t)` 配置，网格使用 `type(dvr_1d_t)` 初始化，双态耦合使用 `propagate_split_operator_2channel`。
5. 过程子程序严格遵循 Fortran 2008 标准，所有参数明确标注 `intent(in/out/inout)`。
```

### 5.2 常见物理场景任务提示词 (User Prompts)

#### 算例 A: 任意分子势能面束缚态求解
> *"请使用 `general_module` 编写一个 Fortran 程序，利用 Sinc-DVR（`dvr_sinc_init`）和 FGH 求解器（`fgh_solve_bound_states`），计算一维双阱势 $V(x) = \frac{1}{2} x^4 - 2 x^2$ 的前 6 个能级本征值与对称/反对称波函数，并将能级劈裂输出到控制台。"*

#### 算例 B: 强场高次谐波发射与截止能计算
> *"请使用 `general_module` 构建一个 800nm、30fs 高斯脉冲作用于氖原子（Neon）的模型，调用 `get_atom_config('Ne', atom)`，利用 Lewenstein 强场近似（`lewenstein_sfa_dipole`）生成时域偶极加速度，并用 `hhg_power_spectrum` 输出高次谐波功率谱，对比半经典截断能 $I_p + 3.17 U_p$。"*

#### 算例 C: 极性分子太赫兹场无场定向与玻尔兹曼系综平均
> *"请使用 `general_module` 编写一个刚体转子分子在室温（$T=300\text{ K}$）下的定向动力学程序。调用 `boltzmann_rotational_weights` 计算各 $J$ 态的初态布居权重，并在外加脉冲后模拟自由演化量子拍，调用 `thermal_average_2d` 输出系综平均定向度 $\langle\cos\theta\rangle(t)$。"*

#### 算例 D: 激光脉冲调控分子转振态相干布居转移与梯级跃迁
> *"请使用 `general_module` 编写一个超快红外激光脉冲调控双原子分子转振态布居的完整程序。首先用 Sinc-DVR 求解 Morse 势束缚态波函数并积分各态转动常数 $B_v$ 与偶极矩阵，调用 `build_rovibrational_hamiltonian` 与 `build_rovibrational_dipole_matrix` 构建全空间 $|v, J\rangle$ 哈密顿量与跃迁偶极，配置共振 $\sin^2$ 激光脉冲驱动体系从基态 $|v=0, J=0\rangle$ 发生相干布居转移，RK4 演化 TDSE 并使用 `save_data_table_2d` 保存全转振态时序数据。"*

#### 算例 E: 超冷磁偶极自旋弛豫与极性分子电偶极 Stark 效应
> *"请使用 `general_module` 编写一个超冷偶极体系动力学分析程序。首先调用 `calc_dipolar_relaxation_cross_section` 和 `calc_dipolar_relaxation_thermal_rate` 计算弱磁阱中 $^{87}\text{Rb}$ 原子在 $T=1.0\,\mu\text{K}$ 下的各向异性磁偶极两体自旋弛豫截面与热平均速率常数；接着配置极性双原子分子 $^{40}\text{K}^{87}\text{Rb}$，在外加直流电场 $\mathcal{E} \in [0, 20]\,\text{kV/cm}$ 下调用 `calc_stark_induced_dipole` 扫描诱导电偶极矩 $d_{\text{ind}}(\mathcal{E})$ 并计算对应的特征电偶极长度 $a_d$。"*

#### 算例 F: 超冷原子光缔合谱学与双光子 Raman 缔合
> *"请使用 `general_module` 编写一个超冷光缔合（Photoassociation）谱学模拟程序。配置基态碰撞波函数与激发态束缚振动态，调用 `calc_free_bound_fc_overlap` 计算自由-束缚态 Franck-Condon 因子，并利用 `calc_pa_detuning_scan` 在温度 $T=50\,\mu\text{K}$ 下模拟激光失谐范围 $[-1.0, 1.0]\,\text{GHz}$ 的光缔合吸收谱，提取共振峰值与热展宽线宽，最后调用 `calc_twophoton_raman_coupling` 评估制备振转基态分子的双光子 STIRAP 有效耦合拉比频率。"*

### 5.3 AI 编写 Fortran 常见陷阱自查表
| 常见 AI 陷阱 | 产生原因 | 正确做法 |
| :--- | :--- | :--- |
| **数组从 0 开始索引** | 习惯 Python/C 导致越界 | Fortran 默认下标从 1 开始（转动态除外，转动量子数可显式声明 `(0:j_max)`） |
| **整数相除精度丢失** | 如写成 `1 / 2 * mass` 得 `0` | 必须书写浮点字面量：`0.5_dp * mass` 或 `1.0_dp / 2.0_dp` |
| **Pure 函数内修改变量** | 违背 `pure` 语义导致编译报错 | `pure` 函数内不能包含任何 `intent(out/inout)`、I/O 打印语句或全局状态 |
| **忘记分配或重复分配** | `allocate()` 逻辑不严密 | 分配前使用 `if (.not. allocated(arr)) allocate(arr(...))` 进行防重保护 |

---

## 6. 常见错误排查自查表 (Troubleshooting & FAQ)

### Q1: 编译报 `Fatal Error: Cannot open module file 'general_module.mod'`
- **原因**：编译器未在包含路径中找到已编译的 `.mod` 模块接口文件。
- **解决办法**：
  1. 确保先在 `GeneralModule/src` 目录下运行过编译命令。
  2. 编译主程序时加上 `-I/path/to/GeneralModule/src` 包含路径。

### Q2: 链接报 `Undefined symbols for architecture ...`
- **原因**：仅包含了 `.mod` 文件，但链接阶段未传入对应的目标文件 `.o` 或静态库 `libgeneral_module.a`。
- **解决办法**：
  在编译命令行末尾加上静态库文件：`gfortran -I.../src main.f90 .../src/libgeneral_module.a -o main`。

### Q3: 波包演化一段时间后在边界产生剧烈的锯齿震荡或反射
- **原因**：未开启复吸收边界（CAP），或吸收区起点设置过靠近网格边界，导致波包撞墙反弹。
- **解决办法**：
  1. 增大物理网格空间范围。
  2. 确保在演化循环中每一步调用 `cap_apply_mask(x_grid, dt, cap, psi)`。
  3. 调整 `cap_init` 中的强度至 $0.1 \sim 0.3$。

### Q4: 谐波谱（HHG）在低频区出现大面积白噪声或截断不明显
- **原因**：时域偶极加速度在脉冲开关时存在微小的数值直流偏移（DC offset），傅里叶变换产生强烈的频谱泄漏。
- **解决办法**：
  直接使用内置的 `hhg_power_spectrum`，该例程已自动内置 Hanning 平滑时间窗消除端点泄漏。
