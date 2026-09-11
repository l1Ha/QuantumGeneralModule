# GeneralModule 算法库详细配置与集成部署指南

本文档为 `GeneralModule` 现代量子动力学算法库的权威配置、编译构建、参数调优及多语言混合编程指南。无论您是高校/科研院所的研究人员、高性能计算（HPC）工程师，还是配合 AI Coding Agent（如 Antigravity、Cursor、Copilot）进行算法开发，本手册都提供了完备的配置说明。

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
   - [3.5 时间推进步长稳定性判据 (CFL 条件)](#35-时间推进步长稳定性判据-cfl-条件)
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

### 3.5 时间推进步长稳定性判据 (CFL 条件)
使用 `propagate_split_operator_1d` 或 `propagate_split_operator_2d` 进行时间演化时，为保证波包演化的辛对称幺正性与相位精度，时间步长 $\Delta t$ 必须满足空间动能离散的 Courant-Friedrichs-Lewy (CFL) 上限：
$$\Delta t \le \frac{2 m \Delta x^2}{\pi \hbar}$$
- **电子动力学（$m = 1\text{ a.u.}$）**：若 $\Delta x = 0.2\text{ Bohr}$，则 $\Delta t \le \frac{2 \times 0.04}{3.14} \approx 0.025\text{ a.u.} \approx 0.6\text{ 自动单位 (约 } 0.0006\text{ fs)}$。
- **核动力学（$m \sim 2000\text{ a.u.}$）**：步长可大幅放宽至 $\Delta t \sim 1.0 - 5.0\text{ a.u.} (0.02 - 0.1\text{ fs})$。
- **切比雪夫推进器（`chebyshev_propagate_step`）**：时间步长不受上述高频震荡严格限制，单步可达几个飞秒且保持机器精度。

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
