---
title: 量子动力学数值计算
subtitle: 基于 GeneralModule 的 Fortran 2008 理论与实践
author: 李皓
email: LIH_ao@outlook.com
edition: 第三版
date: 2026-09-25
---

# 关于本书

本书围绕 GeneralModule 开源代码库，系统讨论量子动力学问题的数学表述、离散化方法、数值算法和工程实现。本书的写作目标不是罗列 API，而是建立一条从物理模型到可验证代码的完整链条。

本书强调三个原则：

1. 理论先行。每一个数值方法都从算符、表象和边界条件出发。
2. 算法明确。每一个实现都说明其数学假设、复杂度和适用范围。
3. 结果可验证。每一个物理量都必须有解析极限、守恒量或独立实现作为对照。

本书对应的代码仓库为：

- 仓库主页：<https://github.com/l1Ha/QuantumGeneralModule>
- 源码二维码：

![](qr/repo_root.png){width=2.0cm}

本书使用 GeneralModule 的 MIT 许可版本。所有公式和算法均给出文献来源或项目内源码定位。若未来引入第三方代码，必须在仓库许可证和变更记录中明确说明其来源、许可证和修改范围。

# 符号与约定

本书使用以下符号：

| 符号 | 含义 |
|---|---|
| $\hbar$ | 约化普朗克常数 |
| $\hat H$ | 哈密顿算符 |
| $\hat T$ | 动能算符 |
| $\hat V$ | 势能算符 |
| $\psi$ | 波函数 |
| $\rho$ | 密度矩阵 |
| $a_s$ | 散射长度 |
| $S_{ij}$ | 散射矩阵元 |
| $J$ | 总角动量量子数 |
| $L$ | 轨道角动量量子数 |
| $DVR$ | 离散变量表象 |
| $FGH$ | Fourier Grid Hamiltonian |

除特别说明外，本书在原子单位下工作：

$$
\hbar=1,\qquad m_e=1,\qquad e=1,\qquad 4\pi\varepsilon_0=1.
$$

能量单位为 Hartree，长度单位为 Bohr。若使用实验单位，应先通过 `mod_constants` 中的接口转换为原子单位。

# 第一部分　量子力学与坐标表象

# 第 1 章　态、算符与演化

## 1.1　希尔伯特空间与内积

量子态属于复希尔伯特空间 $\mathcal{H}$。设 $|\psi\rangle,|\phi\rangle\in\mathcal{H}$，内积为

$$
\langle \psi|\phi\rangle
=
\int_{\Omega}\psi^{*}(\mathbf{x})\phi(\mathbf{x})\,d\mathbf{x}.
$$

在离散网格上，该内积变为

$$
\langle \psi|\phi\rangle_N
=
\sum_{i=1}^{N}w_i\psi_i^{*}\phi_i,
$$

其中 $w_i$ 是数值积分权重。若忽略权重，即使本征值看起来收敛，期望值和归一化仍可能错误。

## 1.2　算符与厄米性

可观测量由自伴算符表示。若

$$
\hat A=\hat A^{\dagger},
$$

则其本征值为实数，且不同本征值对应的本征态正交。对任意态 $|\psi\rangle$，

$$
\langle A\rangle
=
\langle\psi|\hat A|\psi\rangle.
$$

哈密顿算符的自伴性保证封闭系统演化保持范数：

$$
\frac{d}{dt}\langle\psi|\psi\rangle=0.
$$

若离散哈密顿量不满足 $H=H^{\dagger}$，则演化算符不再幺正，必须检查离散化过程是否破坏对称性。

## 1.3　演化算符

对不含时哈密顿量，含时薛定谔方程

$$
i\hbar\frac{\partial}{\partial t}|\psi(t)\rangle
=
\hat H|\psi(t)\rangle
$$

的形式解为

$$
|\psi(t)\rangle
=
\hat U(t,t_0)|\psi(t_0)\rangle,
$$

$$
\hat U(t,t_0)
=
\exp\left[-\frac{i}{\hbar}\hat H(t-t_0)\right].
$$

若哈密顿量随时间变化，则

$$
\hat U(t,t_0)
=
\mathcal{T}
\exp\left[
-\frac{i}{\hbar}
\int_{t_0}^{t}\hat H(t')\,dt'
\right],
$$

其中 $\mathcal{T}$ 为时间排序算符。数值传播器的作用是近似该时间排序指数。

## 1.4　守恒量与数值诊断

若 $[\hat H,\hat A]=0$，则

$$
\frac{d}{dt}\langle A\rangle=0.
$$

在数值计算中，至少应检查以下守恒量或近似守恒量：

- 波函数范数 $\langle\psi|\psi\rangle$；
- 能量 $\langle H\rangle$；
- 总角动量 $\hat J^2$；
- 总磁量子数 $M$；
- 密度矩阵迹 $\operatorname{Tr}\rho$。

若某守恒量出现漂移，应先区分物理近似（例如吸收边界）和数值误差（例如时间步长过大）。

## 1.5　代码定位

实对称本征求解和矩阵求逆位于：

- [src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90)
- 二维码：

![](qr/src__mod_linear_algebra.f90.png){width=1.8cm}

# 第 2 章　单位制与数值精度

## 2.1　原子单位

原子单位定义为

$$
\hbar=1,\qquad m_e=1,\qquad e=1,\qquad 4\pi\varepsilon_0=1.
$$

由此得到 Hartree 能量

$$
E_h=\frac{\hbar^2}{m_e a_0^2},
$$

和原子单位时间

$$
t_0=\frac{\hbar}{E_h}.
$$

若将哈密顿量写为

$$
\hat H=E_h\,\tilde{\hat H},
$$

则薛定谔方程变为

$$
i\frac{\partial\tilde\psi}{\partial\tilde t}
=
\tilde{\hat H}\tilde\psi.
$$

## 2.2　单位转换接口

GeneralModule 的 `mod_constants` 提供 `to_au` 和 `from_au` 接口。核心源码为：

- [src/mod_constants.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90)
- 二维码：

![](qr/src__mod_constants.f90.png){width=1.8cm}

该模块集中定义了：

- 长度：Bohr、Å、nm、m；
- 能量：Hartree、eV、$\mathrm{cm^{-1}}$、K、J；
- 时间：原子单位、fs、ps、s；
- 电场：V/m、MV/cm；
- 磁场：Tesla、Gauss。

推荐做法是：输入端立即转换为原子单位，计算全程保持原子单位，输出端再转换为目标单位。这样可以将单位错误限制在输入和输出层。

## 2.3　误差传播

设观测量为 $f(\mathbf{x})$，输入扰动为 $\delta x_i$，则一阶误差传播为

$$
\delta f
=
\sum_i
\frac{\partial f}{\partial x_i}\delta x_i.
$$

对矩阵问题，若哈密顿量存在扰动 $\delta H$，则本征值扰动满足

$$
|\delta\lambda_i|
\le
\|\delta H\|_2.
$$

本征向量扰动由能隙控制：

$$
\sin\theta_i
\lesssim
\frac{\|\delta H\|_2}{\operatorname{gap}_i},
\qquad
\operatorname{gap}_i=\min_{j\ne i}|\lambda_i-\lambda_j|.
$$

因此，近简并态的本征向量对数值扰动敏感，但其子空间仍然稳定。

# 第 3 章　角动量理论

## 3.1　角动量算符

角动量算符满足

$$
[\hat J_i,\hat J_j]
=
i\hbar\sum_k\varepsilon_{ijk}\hat J_k.
$$

共同本征态为

$$
\hat J^2|jm\rangle
=
\hbar^2j(j+1)|jm\rangle,
$$

$$
\hat J_z|jm\rangle
=
\hbar m|jm\rangle.
$$

升降算符为

$$
\hat J_\pm=\hat J_x\pm i\hat J_y,
$$

且

$$
\hat J_\pm|jm\rangle
=
\hbar\sqrt{j(j+1)-m(m\pm1)}\,|j,m\pm1\rangle.
$$

## 3.2　Clebsch–Gordan 系数

两个角动量耦合为

$$
|j_1m_1\rangle|j_2m_2\rangle
=
\sum_{JM}
C^{JM}_{j_1m_1j_2m_2}|JM\rangle.
$$

非零系数要求

$$
M=m_1+m_2,
$$

$$
|j_1-j_2|\le J\le j_1+j_2.
$$

Wigner $3j$ 符号与 Clebsch–Gordan 系数的关系为

$$
C^{JM}_{j_1m_1j_2m_2}
=
(-1)^{j_1-j_2+M}\sqrt{2J+1}
\begin{pmatrix}
j_1 & j_2 & J\\
m_1 & m_2 & -M
\end{pmatrix}.
$$

## 3.3　重耦合与 6j、9j 符号

三个角动量的重耦合系数可用 $6j$ 符号表示：

$$
\langle (j_1j_2)J_{12},j_3;JM|j_1,(j_2j_3)J_{23};JM\rangle
=
(-1)^{j_1+j_2+j_3+J}
\sqrt{(2J_{12}+1)(2J_{23}+1)}
\begin{Bmatrix}
j_1 & j_2 & J_{12}\\
j_3 & J & J_{23}
\end{Bmatrix}.
$$

四个角动量重耦合则由 $9j$ 符号描述。GeneralModule 的 `mod_special_functions` 实现了 Legendre 多项式、缔合 Legendre 函数、Wigner $3j$、Clebsch–Gordan、$6j$ 和 $9j$ 符号。

- 源码：

[src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)

- 二维码：

![](qr/src__mod_special_functions.f90.png){width=1.8cm}

# 第 4 章　空间固定与体固定坐标

## 4.1　空间固定坐标

空间固定坐标也称实验室坐标。其原点通常取分子质心，坐标轴由外部实验装置定义。设分子轴与空间固定 $Z$ 轴夹角为 $\theta$，则外场耦合常写为

$$
\hat H_{\mathrm{ext}}
=
-d\mathcal{E}\cos\theta
+
g\mu_B B S_Z.
$$

空间固定表象适合描述外场方向、散射方向、光子偏振和实验观测角分布。

## 4.2　体固定坐标

体固定坐标随分子一起转动。空间固定坐标 $\mathbf{r}^{(S)}$ 与体固定坐标 $\mathbf{r}^{(B)}$ 通过 Euler 角 $\alpha,\beta,\gamma$ 联系：

$$
\mathbf{r}^{(S)}
=
\mathsf{R}(\alpha,\beta,\gamma)\mathbf{r}^{(B)}.
$$

球张量算符在两套坐标之间的变换为

$$
T^{(k)}_q(S)
=
\sum_{q'=-k}^{k}
D^{(k)*}_{q'q}(\alpha,\beta,\gamma)
T^{(k)}_{q'}(B).
$$

体固定表象适合描述分子内部振动、转动常数、Coriolis 耦合、离心畸变和分子对称性。

## 4.3　Watson 转振哈密顿量

对非线性分子，体固定坐标中的转振哈密顿量可写为

$$
\hat H_{\mathrm{vr}}
=
\hat T_{\mathrm{vib}}
+
\hat T_{\mathrm{rot}}
+
\hat T_{\mathrm{Cor}}
+
V(\mathbf q)
+
\hat V_{\mathrm{cd}}.
$$

振动动能为

$$
\hat T_{\mathrm{vib}}
=
-\frac{\hbar^2}{2}
\sum_{i,j}^{3N-6}
\frac{\partial}{\partial q_i}
G_{ij}(\mathbf q)
\frac{\partial}{\partial q_j}.
$$

转动动能为

$$
\hat T_{\mathrm{rot}}
=
\frac12
(\hat{\mathbf J}-\hat{\boldsymbol{\pi}})^T
\mathbf A(\mathbf q)
(\hat{\mathbf J}-\hat{\boldsymbol{\pi}}),
$$

其中 $\hat{\boldsymbol{\pi}}$ 是振动角动量，$\mathbf A(\mathbf q)$ 是逆惯量张量。展开得到

$$
\hat T_{\mathrm{rot}}
+
\hat T_{\mathrm{Cor}}
=
\frac12
\hat{\mathbf J}^T\mathbf A\hat{\mathbf J}
-
\hat{\mathbf J}^T\mathbf A\hat{\boldsymbol{\pi}}
+
\frac12
\hat{\boldsymbol{\pi}}^T\mathbf A\hat{\boldsymbol{\pi}}.
$$

中间项为 Coriolis 耦合。刚性转子极限下，

$$
\hat H_{\mathrm{rot}}=B\hat{\mathbf J}^2.
$$

对称陀螺分子满足

$$
\hat H_{\mathrm{rot}}
=
B\hat{\mathbf J}^2
+
(A-B)\hat J_z^2.
$$

离心畸变的最低阶修正为

$$
\hat V_{\mathrm{cd}}
=
-D_J\hat{\mathbf J}^4
-
D_{JK}\hat{\mathbf J}^2\hat J_z^2
-
D_K\hat J_z^4.
$$

## 4.4　Jacobi 坐标与超球坐标

三原子反应 $A+BC$ 常用 Jacobi 坐标

$$
(r,R,\gamma),
$$

其中 $r$ 是 $BC$ 键长，$R$ 是原子 $A$ 到 $BC$ 质心的距离，$\gamma$ 是两者夹角。三原子核间距 $(r_{AB},r_{BC},r_{CA})$ 与 Jacobi 坐标存在解析映射。

在超球坐标中，总尺度 $\rho$ 和超角 $\Omega$ 满足

$$
d\mathbf{R}
=
\rho^2\sin\chi\,d\rho\,d\chi\,d\Omega.
$$

这使多个反应通道可以组织为单一超径向坐标与一组角坐标。

- Jacobi 坐标源码：

[src/mod_triatomic_geometry.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_triatomic_geometry.f90)

- 二维码：

![](qr/src__mod_triatomic_geometry.f90.png){width=1.8cm}

- 超球反应源码：

[src/mod_hyperspherical_reactive.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_hyperspherical_reactive.f90)

- 二维码：

![](qr/src__mod_hyperspherical_reactive.f90.png){width=1.8cm}

# 第二部分　离散表示与谱方法

# 第 5 章　基组展开与 DVR

## 5.1　基组展开

设

$$
|\psi\rangle=\sum_{n=1}^{N}c_n|\chi_n\rangle.
$$

将薛定谔方程投影到基组得广义本征值问题

$$
\sum_j H_{ij}c_j
=
E\sum_j S_{ij}c_j,
$$

其中

$$
H_{ij}=\langle\chi_i|\hat H|\chi_j\rangle,
\qquad
S_{ij}=\langle\chi_i|\chi_j\rangle.
$$

若基组正交，则 $S_{ij}=\delta_{ij}$，问题退化为

$$
H\mathbf c=E\mathbf c.
$$

## 5.2　离散变量表象

DVR 的核心性质是势能矩阵近似对角：

$$
V_{ij}
\approx
V(x_i)\delta_{ij}.
$$

于是哈密顿矩阵为

$$
H_{ij}
=
T_{ij}
+
V(x_i)\delta_{ij}.
$$

这与传统基函数展开不同：在谐振子基或 Morse 基中，势能矩阵通常是稠密的；在 DVR 中，势能矩阵是对角的，而动能矩阵的结构由所选 DVR 决定。

# 第 6 章　Sinc-DVR、Fourier 网格与 Legendre-DVR

## 6.1　Sinc-DVR

对均匀网格

$$
x_i=x_{\min}+i\Delta x,
\qquad
i=1,\ldots,N,
$$

Colbert–Miller Sinc-DVR 的动能矩阵元为

$$
T_{ij}
=
\frac{\hbar^2}{2m\Delta x^2}
\begin{cases}
\pi^2/3, & i=j,\\[4pt]
2(-1)^{i-j}/(i-j)^2, & i\ne j.
\end{cases}
$$

在原子单位下，

$$
T_{ij}
=
\frac{1}{2m\Delta x^2}
\begin{cases}
\pi^2/3, & i=j,\\[4pt]
2(-1)^{i-j}/(i-j)^2, & i\ne j.
\end{cases}
$$

总哈密顿量为

$$
H_{ij}
=
T_{ij}
+
V(x_i)\delta_{ij}.
$$

若离散内积权重为 $w_i=\Delta x$，则连续归一化波函数为

$$
\psi(x_i)=\frac{c_i}{\sqrt{\Delta x}}.
$$

## 6.2　Fourier 网格与 FGH

长度为 $L$ 的周期域上，Fourier 基为

$$
\phi_k(x)
=
\frac{1}{\sqrt L}e^{ikx},
\qquad
k=\frac{2\pi n}{L}.
$$

动量算符在该基下对角：

$$
\hat p\,\phi_k(x)=\hbar k\,\phi_k(x),
$$

因此

$$
T_k=\frac{\hbar^2k^2}{2m}.
$$

FGH 方法使用坐标网格表示势能、动量网格表示动能，从而得到

$$
H_{ij}
=
T_{ij}
+
V(x_i)\delta_{ij}.
$$

FGH 与 Split-Operator 方法共享同一傅里叶离散结构：本征态计算使用矩阵对角化，时间演化使用 FFT。

## 6.3　Legendre-DVR

分子取向问题常用 $x=\cos\theta$ 作为坐标。Legendre 多项式满足

$$
(1-x^2)P_n''(x)-2xP_n'(x)+n(n+1)P_n(x)=0,
$$

且

$$
\int_{-1}^{1}P_n(x)P_m(x)\,dx
=
\frac{2}{2n+1}\delta_{nm}.
$$

Gauss-Legendre 求积为

$$
\int_{-1}^{1}f(x)\,dx
\approx
\sum_{i=1}^{N}w_if(x_i).
$$

在角度坐标中，

$$
\int_0^\pi f(\theta)\sin\theta\,d\theta
=
\int_{-1}^{1}f(\arccos x)\,dx
\approx
\sum_iw_if(\arccos x_i).
$$

转子动能算符为

$$
\hat J^2
=
-\hbar^2
\left[
\frac{1}{\sin\theta}
\frac{\partial}{\partial\theta}
\left(
\sin\theta
\frac{\partial}{\partial\theta}
\right)
\right].
$$

在球谐基 $|lm\rangle$ 中，

$$
\hat J^2|lm\rangle
=
\hbar^2l(l+1)|lm\rangle.
$$

在 DVR 网格上，它变为矩阵 $J^2_{ij}$。

## 6.4　DVR 不止一种

DVR 是一族方法，不是单一算法。常见形式包括：

| DVR 类型 | 坐标 | 动能矩阵 | 典型用途 |
|---|---|---|---|
| Sinc-DVR | 等距 $x$ | 解析长程衰减矩阵 | 一维局域势能面 |
| Fourier/FGH | 周期 $x$ | 动量空间对角 | 周期或大范围波包 |
| Gauss-Legendre DVR | $x=\cos\theta$ | 角动量矩阵 | 分子取向与转子 |
| Gauss-Hermite DVR | 无限域谐振子权重 | 谐振子结构 | 振动模式 |
| Laguerre-DVR | 半无限径向域 | 径向结构 | 径向库仑或长程势 |
| Distributed Gaussian DVR | 可调局域高斯 | 动能重叠矩阵 | 多维非均匀势能面 |

因此，将 DVR 等同于 Sinc-DVR 是不完整的。GeneralModule 目前提供 Sinc-DVR 与 Gauss-Legendre DVR，并可通过同一框架扩展到其他正交多项式网格。

- DVR 源码：

[src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90)

- 二维码：

![](qr/src__mod_dvr_grid.f90.png){width=1.8cm}

# 第 7 章　张量积网格与收敛性

## 7.1　多维张量积 DVR

若体系有 $d$ 个自由度，则张量积网格为

$$
(q_{1,i_1},q_{2,i_2},\ldots,q_{d,i_d}).
$$

总网格数为

$$
N_{\mathrm{tot}}=\prod_{\alpha=1}^{d}N_\alpha.
$$

若动能可分离，

$$
\hat H=\sum_{\alpha}\hat T_\alpha+\hat V(q_1,\ldots,q_d),
$$

则

$$
\hat T_1\rightarrow T_1\otimes I_2,
\qquad
\hat T_2\rightarrow I_1\otimes T_2.
$$

可分离势能项为

$$
V_1(q_1)
\rightarrow
\operatorname{diag}[V_1(q_{1,i})]\otimes I_2.
$$

耦合势能项在网格点集上仍是对角的，但以全维度索引展开。

## 7.2　复杂度

若每个维度有 $N$ 个点，$d$ 维张量积网格大小为 $N^d$。稠密对角化复杂度为

$$
O(N^{3d}),
$$

存储复杂度为

$$
O(N^{2d}).
$$

因此，直接对角化只适用于低维或截断后的小问题。高维问题需要稀疏矩阵、Krylov 方法、Chebyshev 传播、张量分解或缩减基方法。

## 7.3　收敛性

设精确本征值为 $E$，数值本征值为 $E_N$。对光滑势函数，谱方法通常表现出快速收敛：

$$
E_N-E
\sim
C_Ne^{-\alpha N}.
$$

若势函数不可微或有奇点，则收敛可能退化为代数阶：

$$
E_N-E
\sim
CN^{-p}.
$$

若使用 Richardson 外推，可由三种网格估计阶数：

$$
p
=
\frac{
\ln\left[
\frac{E_{N_1}-E_{N_2}}
{E_{N_2}-E_{N_3}}
\right]
}
{\ln(N_2/N_1)}.
$$

若 $p$ 不稳定，则说明尚未进入渐近区，或存在边界、奇点、混叠等系统误差。

## 7.4　误差来源

| 误差来源 | 数学表现 | 诊断方法 |
|---|---|---|
| 空间截断 | 边界反射 | 改变区间 |
| 网格欠采样 | 混叠、高频振荡 | 改变 $\Delta x$ |
| 势能奇点 | 局部误差放大 | 局部加密或变换坐标 |
| 非厄米离散 | 范数漂移 | 检查 $H=H^T$ |
| 本征求解截断 | 只得到部分谱 | 与解析极限比较 |

- 测试源码：

[tests/test_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_dvr_grid.f90)

- 二维码：

![](qr/tests__test_dvr_grid.f90.png){width=1.8cm}

# 第三部分　含时演化、开放系统与量子控制

# 第 8 章　含时薛定谔方程

## 8.1　演化算符

含时薛定谔方程为

$$
i\hbar\frac{\partial}{\partial t}|\psi(t)\rangle
=
\hat H(t)|\psi(t)\rangle.
$$

若 $\hat H$ 不含时，则

$$
|\psi(t)\rangle
=
\hat U(t,t_0)|\psi(t_0)\rangle,
$$

$$
\hat U(t,t_0)
=
\exp\left[-\frac{i}{\hbar}\hat H(t-t_0)\right].
$$

若 $\hat H$ 含时，则

$$
\hat U(t,t_0)
=
\mathcal{T}
\exp\left[
-\frac{i}{\hbar}
\int_{t_0}^{t}\hat H(t')\,dt'
\right].
$$

## 8.2　可观测量

坐标期望值为

$$
\langle x\rangle(t)
=
\int x|\psi(x,t)|^2\,dx.
$$

动量期望值为

$$
\langle p\rangle(t)
=
\int \psi^{*}(x,t)
\left(-i\hbar\frac{\partial}{\partial x}\right)
\psi(x,t)\,dx.
$$

在动量空间中，

$$
\langle p\rangle(t)
=
\int \hbar k\,|\tilde\psi(k,t)|^2\,dk.
$$

# 第 9 章　Split-Operator 方法

## 9.1　分裂公式

设

$$
\hat H=\hat T+\hat V.
$$

由于 $\hat T$ 与 $\hat V$ 不对易，

$$
e^{-i(\hat T+\hat V)\Delta t/\hbar}
\ne
e^{-i\hat T\Delta t/\hbar}
e^{-i\hat V\Delta t/\hbar}.
$$

一阶分裂为

$$
e^{-i(\hat T+\hat V)\Delta t/\hbar}
=
e^{-i\hat T\Delta t/\hbar}
e^{-i\hat V\Delta t/\hbar}
+
O(\Delta t^2).
$$

二阶对称 Strang 分裂为

$$
\hat U(\Delta t)
=
e^{-i\hat V\Delta t/(2\hbar)}
e^{-i\hat T\Delta t/\hbar}
e^{-i\hat V\Delta t/(2\hbar)}
+
O(\Delta t^3).
$$

Feit、Fleck 和 Steiger 将该谱方法系统用于薛定谔方程求解【Feit, Fleck & Steiger 1982】。Kosloff 对含时量子动力学方法作了系统综述【Kosloff 1988】。

## 9.2　坐标与动量表象

在坐标表象中，

$$
[\hat V\psi](x)=V(x)\psi(x).
$$

在动量表象中，

$$
[\hat T\tilde\psi](k)=\frac{\hbar^2k^2}{2m}\tilde\psi(k).
$$

因此一步传播为

$$
\psi(x,t+\Delta t)
=
e^{-iV(x)\Delta t/(2\hbar)}
\mathcal{F}^{-1}
\left[
e^{-i\hbar k^2\Delta t/(2m)}
\mathcal{F}
\left[
e^{-iV(x)\Delta t/(2\hbar)}
\psi(x,t)
\right]
\right].
$$

在原子单位下，

$$
\psi(x,t+\Delta t)
=
e^{-iV(x)\Delta t/2}
\mathcal{F}^{-1}
\left[
e^{-ik^2\Delta t/(2m)}
\mathcal{F}
\left[
e^{-iV(x)\Delta t/2}
\psi(x,t)
\right]
\right].
$$

## 9.3　数值步骤

1. 计算势能半步相位 $e^{-iV(x_i)\Delta t/2}$；
2. 将波函数乘以该相位；
3. FFT 变换到动量空间；
4. 乘以动能相位 $e^{-ik^2\Delta t/(2m)}$；
5. 逆 FFT 回到坐标空间；
6. 再次乘以势能半步相位。

每步需要两次 FFT，复杂度为 $O(N\log N)$，显著优于稠密矩阵指数的 $O(N^3)$。

## 9.4　范数与稳定性

每个指数相位因子都是幺正算符。因此在无穷精度下，Split-Operator 严格保持范数：

$$
\|\psi(t+\Delta t)\|_2=\|\psi(t)\|_2.
$$

有限精度误差来源包括：

- FFT 归一化不一致；
- 网格混叠；
- 复势导致非幺正演化；
- 吸收边界；
- 长时间舍入误差；
- 不当归一化掩盖真实误差。

因此应监测

$$
\epsilon_{\mathrm{norm}}(t)
=
\left|\|\psi(t)\|_2-1\right|.
$$

- 波包传播源码：

[src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)

- 二维码：

![](qr/src__mod_wavepacket_propagator.f90.png){width=1.8cm}

# 第 10 章　FFT、混叠与吸收边界

## 10.1　动量网格

均匀网格 $x_i=x_0+i\Delta x$ 对应动量网格

$$
k_j
=
\frac{2\pi}{L}
\begin{cases}
j, & 0\le j<N/2,\\
j-N, & N/2\le j<N.
\end{cases}
$$

其中 $L=N\Delta x$。最大可表示动量为

$$
k_{\max}=\frac{\pi}{\Delta x}.
$$

若波函数在动量空间的有效宽度超过 $|k|<\pi/\Delta x$，则发生混叠：

$$
k_{\mathrm{alias}}=k-\frac{2m\pi}{\Delta x},
\qquad
m\in\mathbb{Z}.
$$

## 10.2　复吸收边界

当波包离开有限计算域时，硬边界会产生非物理反射。复吸收势在边界区域加入虚部势：

$$
V_{\mathrm{eff}}(x)=V(x)-iW(x),
\qquad
W(x)\ge0.
$$

常用 $\sin^2$ 形式为

$$
W(x)
=
\eta
\sin^2
\left[
\frac{\pi}{2}
\frac{x-x_a}{L_a}
\right],
\qquad
x_a\le x\le x_a+L_a.
$$

加入 CAP 后，哈密顿量变为非厄米：

$$
\hat H_{\mathrm{eff}}=\hat T+\hat V-i\hat W.
$$

内部区域波函数的模平方不再守恒，其损失应与流出概率一致：

$$
\frac{d}{dt}
\int_{\Omega_{\mathrm{in}}}|\psi|^2\,dx
=
-
\int_{\partial\Omega_{\mathrm{in}}}J_n\,dS
-
\int_{\Omega_{\mathrm{CAP}}}2W(x)|\psi|^2\,dx.
$$

概率流为

$$
J(x,t)
=
\frac{\hbar}{m}
\operatorname{Im}
\left[
\psi^{*}(x,t)
\frac{\partial\psi(x,t)}{\partial x}
\right].
$$

- 吸收边界源码：

[src/mod_absorbing_boundary.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_absorbing_boundary.f90)

- 二维码：

![](qr/src__mod_absorbing_boundary.f90.png){width=1.8cm}

# 第 11 章　Runge–Kutta、多步方法与 Bloch 方程

## 11.1　RK4

将波函数或密度矩阵写成常微分方程

$$
\frac{d\mathbf y}{dt}=f(t,\mathbf y).
$$

经典四阶 Runge–Kutta 方法为

$$
\begin{aligned}
k_1&=f(t_n,y_n),\\
k_2&=f\left(t_n+\frac{\Delta t}{2},y_n+\frac{\Delta t}{2}k_1\right),\\
k_3&=f\left(t_n+\frac{\Delta t}{2},y_n+\frac{\Delta t}{2}k_2\right),\\
k_4&=f(t_n+\Delta t,y_n+\Delta t k_3),\\
y_{n+1}
&=
y_n+
\frac{\Delta t}{6}
(k_1+2k_2+2k_3+k_4).
\end{aligned}
$$

其局部截断误差为 $O(\Delta t^5)$，全局误差为 $O(\Delta t^4)$。

## 11.2　ABM4

Adams–Bashforth 四阶预估式为

$$
y_{n+1}^{(p)}
=
y_n
+
\frac{\Delta t}{24}
(55f_n-59f_{n-1}+37f_{n-2}-9f_{n-3}).
$$

在预估点计算

$$
f_{n+1}^{(p)}=f(t_{n+1},y_{n+1}^{(p)}),
$$

然后使用 Adams–Moulton 四阶校正式

$$
y_{n+1}
=
y_n
+
\frac{\Delta t}{24}
(9f_{n+1}^{(p)}+19f_n-5f_{n-1}+f_{n-2}).
$$

多步方法依赖历史点的一致性和启动策略。若历史导数来自不同时间步长或不同坐标系，阶数会退化。

## 11.3　光学 Bloch 方程

二能级系统可写为

$$
|\psi(t)\rangle=c_1(t)|1\rangle+c_2(t)|2\rangle.
$$

定义 Bloch 矢量

$$
\mathbf R=(u,v,w)^T.
$$

在无耗散旋转波近似下，

$$
\frac{d\mathbf R}{dt}
=
\boldsymbol{\Omega}\times\mathbf R,
\qquad
\boldsymbol{\Omega}=(\Omega_R,0,\Delta)^T.
$$

展开为

$$
\dot u=-\Delta v,
$$

$$
\dot v=\Delta u-\Omega_R w,
$$

$$
\dot w=\Omega_R v.
$$

脉冲面积为

$$
\Theta=\int\Omega_R(t)\,dt.
$$

当 $\Theta=\pi$ 时发生布居反转。

- 测试源码：

[tests/test_propagators.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_propagators.f90)

- 二维码：

![](qr/tests__test_propagators.f90.png){width=1.8cm}

# 第 12 章　Chebyshev 传播

对不含时哈密顿量，可利用 Chebyshev 多项式展开演化算符：

$$
e^{-i\hat Ht/\hbar}
=
\sum_{n=0}^{\infty}a_n(t)T_n(\tilde{\hat H}).
$$

若哈密顿量谱满足

$$
E_{\min}\le\hat H\le E_{\max},
$$

则令

$$
\tilde{\hat H}
=
\frac{2\hat H-(E_{\max}+E_{\min})I}{E_{\max}-E_{\min}}.
$$

Chebyshev 展开系数为

$$
a_n(t)
=
(2-\delta_{n0})(-i)^n
J_n\left(\frac{\Delta E\,t}{2\hbar}\right),
$$

其中 $J_n$ 是第一类 Bessel 函数，$\Delta E=E_{\max}-E_{\min}$。所需项数约为

$$
N_{\mathrm{Cheb}}
\sim
\frac{\Delta E\,t}{2\hbar}
+
c\ln\epsilon^{-1}.
$$

Chebyshev 方法适合大时间步和长时传播，但要求谱界估计准确。若谱外态存在，可能出现非物理放大。

- 源码：

[src/mod_chebyshev_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_chebyshev_propagator.f90)

- 二维码：

![](qr/src__mod_chebyshev_propagator.f90.png){width=1.8cm}

# 第 13 章　Lindblad 开放量子系统

## 13.1　密度矩阵

封闭系统的纯态对应密度矩阵

$$
\hat\rho=|\psi\rangle\langle\psi|.
$$

一般混合态为

$$
\hat\rho
=
\sum_\alpha p_\alpha|\psi_\alpha\rangle\langle\psi_\alpha|,
\qquad
p_\alpha\ge0,
\qquad
\sum_\alpha p_\alpha=1.
$$

密度矩阵满足

$$
\hat\rho=\hat\rho^{\dagger},
\qquad
\operatorname{Tr}\hat\rho=1,
\qquad
\hat\rho\ge0.
$$

## 13.2　Lindblad 方程

封闭系统的刘维尔方程为

$$
\frac{d\hat\rho}{dt}
=
-\frac{i}{\hbar}[\hat H,\hat\rho].
$$

马尔可夫近似下的 Lindblad 主方程为

$$
\frac{d\hat\rho}{dt}
=
-\frac{i}{\hbar}[\hat H,\hat\rho]
+
\sum_k\gamma_k
\left(
\hat L_k\hat\rho\hat L_k^{\dagger}
-
\frac12
\{\hat L_k^{\dagger}\hat L_k,\hat\rho\}
\right).
$$

该形式由 Gorini–Kossakowski–Sudarshan 和 Lindblad 给出【Gorini, Kossakowski & Sudarshan 1976；Lindblad 1976】。

## 13.3　耗散通道

自发辐射跃迁算符为

$$
\hat L_{\mathrm{rel}}=\sqrt{\gamma}|g\rangle\langle e|.
$$

纯退相位算符为

$$
\hat L_{\phi}=\sqrt{\gamma_\phi}|e\rangle\langle e|.
$$

二能级弛豫满足

$$
\frac{d\rho_{ee}}{dt}=-\gamma\rho_{ee},
$$

$$
\frac{d\rho_{eg}}{dt}
=
-\left(i\omega_0+\frac{\gamma}{2}+\gamma_\phi\right)\rho_{eg}.
$$

因此布居弛豫时间为 $1/\gamma$，相干衰减率为

$$
T_2^{-1}=\frac{1}{2T_1}+\gamma_\phi.
$$

## 13.4　纯度、相干度与熵

纯度为

$$
P=\operatorname{Tr}(\hat\rho^2).
$$

冯·诺依曼熵为

$$
S=-\operatorname{Tr}(\hat\rho\ln\hat\rho).
$$

$l_1$ 范数相干度为

$$
C_{l1}(\hat\rho)
=
\sum_{i\ne j}|\rho_{ij}|.
$$

数值实现必须检查：

1. 迹守恒 $|\operatorname{Tr}\hat\rho-1|<\epsilon$；
2. 厄米性 $\|\hat\rho-\hat\rho^{\dagger}\|<\epsilon$；
3. 正定性；
4. 物理耗散趋势。

- 源码：

[src/mod_open_quantum.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_open_quantum.f90)

- 二维码：

![](qr/src__mod_open_quantum.f90.png){width=1.8cm}

# 第 14 章　Krotov 量子最优控制

量子控制的目标是设计外场 $\epsilon(t)$，使演化态在终端时刻逼近目标态。保真度为

$$
F=|\langle\phi_{\mathrm{target}}|\psi(T)\rangle|^2.
$$

最优控制泛函可写为

$$
J
=
F
-
\int_0^T\lambda(t)\epsilon^2(t)\,dt.
$$

Krotov 方法引入伴随态 $|\chi(t)\rangle$，满足反向方程

$$
i\hbar\frac{\partial}{\partial t}|\chi(t)\rangle
=
\hat H(t)|\chi(t)\rangle,
$$

终端条件为

$$
|\chi(T)\rangle=|\phi_{\mathrm{target}}\rangle.
$$

场更新可写为

$$
\epsilon^{(k+1)}(t)
=
\epsilon^{(k)}(t)
+
\frac{S(t)}{\alpha}
\operatorname{Im}
\left[
\langle\chi^{(k)}(t)|\hat\mu|\psi^{(k+1)}(t)\rangle
\right],
$$

其中 $S(t)$ 是脉冲开关函数，$\alpha$ 是场强惩罚参数。Somlói、Kazakov 和 Tannor 给出了分子运动的广义弛豫控制方法【Somlói, Kazakov & Tannor 1993】；Reich、Ndong 和 Koch 讨论了单调收敛条件【Reich, Ndong & Koch 2012】。

- 源码：

[src/mod_optimal_control.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_optimal_control.f90)

- 二维码：

![](qr/src__mod_optimal_control.f90.png){width=1.8cm}

# 第四部分　散射、分子与非绝热动力学

# 第 15 章　定态散射与 Numerov 方法

## 15.1　径向薛定谔方程

对中心势 $V(r)$，径向波函数 $u_l(r)$ 满足

$$
-\frac{\hbar^2}{2\mu}
\frac{d^2u_l}{dr^2}
+
\left[
V(r)
+
\frac{\hbar^2l(l+1)}{2\mu r^2}
\right]u_l(r)
=
Eu_l(r).
$$

在原子单位下，

$$
-\frac{1}{2\mu}u_l''(r)
+
\left[
V(r)
+
\frac{l(l+1)}{2\mu r^2}
\right]u_l(r)
=
Eu_l(r).
$$

## 15.2　Numerov 方法

Numerov 方法用于二阶常微分方程

$$
u''(r)=f(r)u(r).
$$

其三点递推为

$$
\left(1-\frac{\Delta r^2}{12}f_{i+1}\right)u_{i+1}
=
2
\left(1+\frac{5\Delta r^2}{12}f_i\right)u_i
-
\left(1-\frac{\Delta r^2}{12}f_{i-1}\right)u_{i-1}.
$$

对零能散射问题，

$$
-\frac{\hbar^2}{2\mu}u''(r)+V(r)u(r)=0.
$$

若势能在远区趋于零，则渐近形式为

$$
u(r)\rightarrow C(r-a_s).
$$

因此

$$
a_s=r-\frac{u(r)}{u'(r)}.
$$

- 散射源码：

[src/mod_ti_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_ti_scattering.f90)

- 二维码：

![](qr/src__mod_ti_scattering.f90.png){width=1.8cm}

## 15.3　Johnson 对数导数方法

直接传播 $u(r)$ 在深势阱中容易指数溢出。Johnson 方法传播相邻点比值

$$
R_i=\frac{u_i}{u_{i-1}}.
$$

Numerov 递推可改写为

$$
R_{i+1}
=
\frac{
2(1-5w_i)-(1+w_{i-1})/R_i
}
{1+w_{i+1}},
$$

其中

$$
w_i=\frac{\Delta r^2}{12}f_i.
$$

渐近截距可由

$$
a_s
=
\frac{R_Nr_{N-1}-r_N}{R_N-1}
$$

得到。

# 第 16 章　相移、S 矩阵与光学定理

## 16.1　分波展开

散射波函数可展开为

$$
\psi(\mathbf r)
=
\sum_{lm}
\frac{u_l(r)}{r}
Y_{lm}(\hat{\mathbf r}).
$$

在远区，

$$
u_l(r)
\rightarrow
\frac{1}{k}
\sin\left(kr-\frac{l\pi}{2}+\delta_l\right).
$$

散射矩阵元为

$$
S_l=e^{2i\delta_l}.
$$

反应矩阵为

$$
K_l=\tan\delta_l.
$$

跃迁矩阵为

$$
T_l=S_l-I.
$$

分波弹性截面为

$$
\sigma_l
=
\frac{4\pi}{k^2}(2l+1)\sin^2\delta_l.
$$

总弹性截面为

$$
\sigma_{\mathrm{el}}
=
\sum_l\sigma_l.
$$

## 16.2　光学定理

总截面与向前散射振幅满足

$$
\sigma_{\mathrm{tot}}
=
\frac{4\pi}{k}
\operatorname{Im}f(0).
$$

该关系可用于检验散射振幅、相移和截面计算的一致性。

## 16.3　有效力程展开

低能散射的有效力程展开为

$$
k\cot\delta_0(k)
=
-\frac{1}{a_s}
+
\frac12r_0k^2
+
O(k^4).
$$

零能截面为

$$
\sigma(0)=4\pi a_s^2.
$$

对范德华长程势，Gribakin–Flambaum 半经典平均散射长度为

$$
\bar a
=
\frac{2}{\sqrt\pi}
\frac{
\Gamma(5/6)
}{
\Gamma(2/3)
}
\left(
\frac{2\mu C_6}{\hbar^2}
\right)^{1/4}.
$$

# 第 17 章　多通道散射与 Feshbach 共振

## 17.1　多通道耦合方程

设通道波函数为 $\mathbf F(r)$，则耦合径向方程为

$$
-\frac{\hbar^2}{2\mu}
\mathbf F''(r)
+
\left[
\mathbf V(r)
+
\mathbf E_{\mathrm{th}}
+
\frac{\hbar^2\mathbf L^2}{2\mu r^2}
\right]
\mathbf F(r)
=
E\mathbf F(r).
$$

在原子单位下，

$$
-\frac{1}{2\mu}
\mathbf F''(r)
+
\left[
\mathbf V(r)
+
\mathbf E_{\mathrm{th}}
+
\frac{\mathbf L^2}{2\mu r^2}
\right]
\mathbf F(r)
=
E\mathbf F(r).
$$

## 17.2　Johnson 矩阵对数导数方法

多通道 Log-Derivative 方法传播矩阵

$$
\mathbf Y(r)
=
\mathbf F'(r)\mathbf F^{-1}(r).
$$

在匹配半径处，由 $\mathbf Y$ 构造 $\mathbf K$ 矩阵，再由 Cayley 变换得到幺正散射矩阵：

$$
\mathbf S
=
(\mathbf I+i\mathbf K)(\mathbf I-i\mathbf K)^{-1}.
$$

跃迁矩阵为

$$
\mathbf T=\mathbf S-\mathbf I.
$$

态到态跃迁概率为

$$
P_{i\rightarrow j}=|S_{ij}|^2.
$$

## 17.3　Feshbach 共振

闭通道束缚态与开通道连续谱耦合会产生 Feshbach 共振。散射长度随磁场变化可写为

$$
a_s(B)
=
a_{\mathrm{bg}}
\left[
1-\frac{\Delta B}{B-B_0}
\right].
$$

共振宽度为 $\Delta B$，共振位置为 $B_0$。对超冷原子，Feshbach 共振是调控相互作用强度的核心手段。

- 外场散射源码：

[src/mod_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_field_scattering.f90)

- 二维码：

![](qr/src__mod_field_scattering.f90.png){width=1.8cm}

# 第 18 章　含时波包散射

## 18.1　高斯波包

初始高斯波包为

$$
\psi(x,0)
=
\left(\frac{1}{2\pi\sigma^2}\right)^{1/4}
\exp\left[
-\frac{(x-x_0)^2}{4\sigma^2}
+ik_0x
\right].
$$

其动量空间振幅为

$$
\tilde\psi(k,0)
=
\left(\frac{2\sigma^2}{\pi}\right)^{1/4}
\exp\left[
-\sigma^2(k-k_0)^2
-ikx_0
\right].
$$

## 18.2　透射与反射

定义左、右通量积分：

$$
T(t)
=
\int_{x>x_m}|\psi(x,t)|^2\,dx,
$$

$$
R(t)
=
\int_{x<x_l}|\psi(x,t)|^2\,dx.
$$

当 $t\rightarrow\infty$ 时，$T+R$ 应接近 1，除非存在吸收或非弹性通道。

## 18.3　定态与含时方法的交叉验证

定态散射给出能量分辨的 $S(E)$；含时波包通过傅里叶变换也可得到能量分辨量。两者应满足：

- 透射峰位置一致；
- 共振宽度一致；
- 相移或延迟一致；
- 总通量守恒。

GeneralModule 的 `ex07` 示例专门用于定态和含时散射波函数交叉验证。

- 含时散射源码：

[src/mod_td_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_td_scattering.f90)

- 二维码：

![](qr/src__mod_td_scattering.f90.png){width=1.8cm}

# 第 19 章　分子转振光谱与光碎片

## 19.1　转振哈密顿量

双分子振转哈密顿量可近似为

$$
\hat H
=
-\frac{\hbar^2}{2\mu}
\frac{d^2}{dR^2}
+
\frac{\hat{\mathbf J}^2}{2\mu R^2}
+
V(R).
$$

在原子单位下，

$$
\hat H
=
-\frac{1}{2\mu}
\frac{d^2}{dR^2}
+
\frac{\hat{\mathbf J}^2}{2\mu R^2}
+
V(R).
$$

## 19.2　Morse 势

Morse 势为

$$
V(R)
=
D_e
\left[
1-e^{-\beta(R-R_e)}
\right]^2.
$$

解析能级为

$$
E_v
=
\hbar\omega_e\left(v+\frac12\right)
-
\hbar\omega_ex_e\left(v+\frac12\right)^2,
$$

其中

$$
\omega_e=\beta\sqrt{\frac{2D_e}{\mu}},
$$

$$
\omega_ex_e=\frac{\omega_e^2}{4D_e}.
$$

## 19.3　跃迁矩阵元

电偶极跃迁矩阵元为

$$
M_{v'v''J'J''}
=
\int
\chi_{v'J'}^{*}(R)
\mu(R)
\chi_{v''J''}(R)\,dR.
$$

Franck–Condon 因子为

$$
F_{v'v''}
=
\left|
\int
\chi_{v'}^{*}(R)\chi_{v''}(R)\,dR
\right|^2.
$$

## 19.4　光碎片动能释放谱

自相关函数为

$$
C(t)=\langle\psi(0)|\psi(t)\rangle.
$$

吸收谱为

$$
\sigma(\omega)
\propto
\omega
\operatorname{Re}
\int_0^\infty
C(t)e^{i\omega t}W(t)\,dt,
$$

其中 $W(t)$ 是时间窗函数。动能释放谱为

$$
P(E_k)
=
\sqrt{\frac{2E_k}{m}}
|A(E_k)|^2.
$$

- 转振源码：

[src/mod_rovibrational.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rovibrational.f90)

- 二维码：

![](qr/src__mod_rovibrational.f90.png){width=1.8cm}

- 光碎片源码：

[src/mod_photofragment_flux.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_photofragment_flux.f90)

- 二维码：

![](qr/src__mod_photofragment_flux.f90.png){width=1.8cm}

# 第 20 章　非绝热动力学与 FSSH

## 20.1　电子—核耦合

设电子态为 $|\phi_1\rangle,|\phi_2\rangle$，核坐标为 $R$，非绝热耦合为

$$
d_{12}(R)
=
\langle\phi_1|\partial_R|\phi_2\rangle.
$$

核速度为 $v$ 时，电子振幅方程为

$$
i\hbar\dot c_1
=
E_1c_1
-
i\hbar v d_{12}c_2,
$$

$$
i\hbar\dot c_2
=
E_2c_2
+
i\hbar v d_{12}c_1.
$$

## 20.2　Tully 最少开关表面跳跃

FSSH 使用经典核轨迹和量子电子振幅。在活性态 $j$ 上，向态 $k$ 的跳跃概率为

$$
g_{j\rightarrow k}
=
\max
\left[
0,
-\frac{2\Delta t}{|c_j|^2}
\operatorname{Re}
\left(
c_j^{*}c_k\,\mathbf v\cdot\mathbf d_{jk}
\right)
\right].
$$

若随机数小于 $g_{j\rightarrow k}$，则尝试跳跃。跳跃后需要沿非绝热耦合方向重标动量以保持总能量：

$$
\mathbf p'\cdot\hat{\mathbf e}_{jk}
=
\operatorname{sign}
\left(
\mathbf p\cdot\hat{\mathbf e}_{jk}
\right)
\sqrt{
(\mathbf p\cdot\hat{\mathbf e}_{jk})^2
-
2M\Delta E
}.
$$

若动量不足，则为 frustrated hop。

## 20.3　Ehrenfest 平均场动力学

Ehrenfest 方法使用平均场力：

$$
\mathbf F
=
-\left\langle\Psi\left|\frac{\partial\hat H}{\partial\mathbf R}\right|\Psi\right\rangle.
$$

核仍沿经典轨迹运动，电子态保持相干演化。与 FSSH 相比，Ehrenfest 不能自然描述分支化过程。

- FSSH 源码：

[src/mod_surface_hopping_fssh.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_hopping_fssh.f90)

- 二维码：

![](qr/src__mod_surface_hopping_fssh.f90.png){width=1.8cm}

# 第 21 章　高级专题地图

GeneralModule 还包含多个研究级专题模块。它们不是孤立 API，而是同一套数学结构的延伸。

| 专题 | 代表模块 | 核心问题 |
|---|---|---|
| Efimov 物理 | `mod_three_body_recombination` | 三体共振、离散标度不变性 |
| CIR | `mod_confined_scattering` | 低维受限诱导共振 |
| Fano/CCR | `mod_autoionization_fano` | 自电离线型与复标度 |
| 强场 NSDI | `mod_strong_field_nsdi` | 隧穿、重碰撞、双电离 |
| ATAS | `mod_attosecond_transient_absorption` | 阿秒瞬态吸收 |
| PECD | `mod_bicircular_pecd` | 手性光电子圆二色性 |
| Rydberg blockade | `mod_rydberg_blockade` | 里德堡阻塞与多体 |
| 表面散射 | `mod_surface_scattering` | 表面衍射与吸附 |
| Eley–Rideal | `mod_surface_reaction_er` | 气—固反应 |
| 电子摩擦 | `mod_surface_electronic_friction` | 非绝热耗散 |
| GIFAD | `mod_grazing_fast_atom_diffraction` | 掠入射衍射 |
| Dirac 原子 | `mod_relativistic_atomic` | 相对论精细结构 |
| RIXS | `mod_resonant_xray_scattering` | 共振非弹性 X 射线散射 |
| Penning 电离 | `mod_penning_associative_ionization` | 亚稳态碰撞电离 |

这些模块的共同验证策略是：先在解析极限或低维模型中验证，再逐步加入耦合、外场、耗散和多通道效应。

# 第五部分　工程实现、验证与复现

# 第 22 章　Fortran 2008 工程结构

GeneralModule 采用四层结构：

1. 物理常数与单位；
2. 数学基础与线性代数；
3. 通用动力学与离散方法；
4. 物理专题求解器。

所有模块默认 `private`，只导出明确列出的公共实体。大多数过程使用 `intent(in)`、`intent(out)` 和 `intent(inout)` 说明副作用。大量纯数学过程标记为 `pure`。

顶层门面模块为：

- [src/general_module.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/general_module.f90)

它通过 `use` 重新导出各子模块，使用户可以写

```fortran
use general_module
implicit none
```

但内部模块不应反向依赖聚合模块，以避免循环依赖。

# 第 23 章　构建、测试与持续集成

## 23.1　三套构建入口

GeneralModule 提供三种构建方式：

| 入口 | 用途 | 特点 |
|---|---|---|
| `Makefile` | 快速构建、测试、示例 | 显式拓扑顺序 |
| `fpm.toml` | Fortran 包管理 | 显式列出测试与示例 |
| `CMakeLists.txt` | 跨平台与 IDE | 对象库、静态库、共享库、CTest |

- 构建配置：

[fpm.toml](https://github.com/l1Ha/QuantumGeneralModule/blob/main/fpm.toml)

- 二维码：

![](qr/fpm.toml.png){width=1.8cm}

- CMake：

[CMakeLists.txt](https://github.com/l1Ha/QuantumGeneralModule/blob/main/CMakeLists.txt)

- 二维码：

![](qr/CMakeLists.txt.png){width=1.8cm}

## 23.2　测试体系

项目当前包含 41 个 Fortran 测试程序和 36 个物理示例。测试框架由项目自实现，每个测试程序维护计数器和断言逻辑。测试应覆盖：

| 层级 | 检查内容 |
|---|---|
| 单位与常数 | 往返转换误差 |
| 特殊函数 | 选择定则、归一化、对称性 |
| 本征问题 | 正交性、解析谱、条件数 |
| DVR | 网格收敛、边界效应 |
| 传播 | 范数、能量、动量守恒 |
| 散射 | 光学定理、幺正性、解析极限 |
| 开放系统 | 迹、厄米性、正定性 |
| 非绝热 | 分支比、能量守恒、统计收敛 |

- 测试脚本：

[tests/run_all_tests.sh](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/run_all_tests.sh)

- 二维码：

![](qr/tests__run_all_tests.sh.png){width=1.8cm}

## 23.3　持续集成

GitHub Actions 配置位于：

[.github/workflows/ci.yml](https://github.com/l1Ha/QuantumGeneralModule/blob/main/.github/workflows/ci.yml)

其作用是自动化编译、测试和示例运行，防止平台相关错误进入主干。

# 第 24 章　Python 交叉验证与可视化

Python 包 `pygenmod` 使用 NumPy 和 Matplotlib 对部分算法进行独立实现、交叉验证和可视化。它不是 Fortran 库的替代品，而是验证与绘图工具。

- NumPy 源码：

[python/pygenmod/visualizer.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py)

- 二维码：

![](qr/python__pygenmod__visualizer.py.png){width=1.8cm}

Python 测试位于：

[python/test_pygenmod.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/test_pygenmod.py)

- 二维码：

![](qr/python__test_pygenmod.py.png){width=1.8cm}

推荐流程是：

1. Fortran 输出原始数据；
2. Python 独立计算参考值；
3. 绘制对比图；
4. 保存输入参数、随机种子、版本和命令；
5. 将图和数据放入同一实验目录。

# 第 25 章　引用、开源与可复现性

## 25.1　引用原则

本书引用的文献用于说明理论来源和算法背景。GeneralModule 的实现是独立编写的；引用不表示复制第三方源码。若未来引入第三方代码，必须：

1. 保留其许可证；
2. 在源码头注明来源；
3. 在 `CHANGELOG` 中说明修改范围；
4. 在文档中给出原始仓库或论文链接。

## 25.2　可复现实验记录

数值实验应至少记录：

| 字段 | 内容 |
|---|---|
| 实验名称 | 物理问题和目标 |
| 源码版本 | commit、文件时间或 tag |
| 环境 | 编译器、系统、CPU、Python 版本 |
| 参数 | 单位、网格、初态、时间步 |
| 随机性 | 种子或确定性声明 |
| 理论基准 | 解析解、文献值、守恒量 |
| 输出 | 能量、截面、范数、分支比 |
| 收敛判据 | 相对误差和阈值 |
| 复现命令 | 完整脚本 |

# 第 26 章　源码索引与二维码

下表给出关键源码、测试、示例和文档的入口。二维码图像存放在 `docs/qr/`，命名规则是将源码路径中的 `/` 替换为 `__`。

| 内容 | 源码路径 | 二维码 |
|---|---|---|
| 仓库主页 | <https://github.com/l1Ha/QuantumGeneralModule> | `docs/qr/repo_root.png` |
| README | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/README.md> | `docs/qr/README.md.png` |
| 配置指南 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/CONFIG_GUIDE.md> | `docs/qr/CONFIG_GUIDE.md.png` |
| 文献映射 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/LITERATURE.md> | `docs/qr/LITERATURE.md.png` |
| 常数与单位 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90> | `docs/qr/src__mod_constants.f90.png` |
| 特殊函数 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90> | `docs/qr/src__mod_special_functions.f90.png` |
| 线性代数 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90> | `docs/qr/src__mod_linear_algebra.f90.png` |
| DVR | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90> | `docs/qr/src__mod_dvr_grid.f90.png` |
| 波包传播 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90> | `docs/qr/src__mod_wavepacket_propagator.f90.png` |
| 吸收边界 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_absorbing_boundary.f90> | `docs/qr/src__mod_absorbing_boundary.f90.png` |
| Chebyshev | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_chebyshev_propagator.f90> | `docs/qr/src__mod_chebyshev_propagator.f90.png` |
| Lindblad | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_open_quantum.f90> | `docs/qr/src__mod_open_quantum.f90.png` |
| Krotov | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_optimal_control.f90> | `docs/qr/src__mod_optimal_control.f90.png` |
| 定态散射 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_ti_scattering.f90> | `docs/qr/src__mod_ti_scattering.f90.png` |
| 含时散射 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_td_scattering.f90> | `docs/qr/src__mod_td_scattering.f90.png` |
| 外场散射 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_field_scattering.f90> | `docs/qr/src__mod_field_scattering.f90.png` |
| 转振光谱 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rovibrational.f90> | `docs/qr/src__mod_rovibrational.f90.png` |
| 光碎片 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_photofragment_flux.f90> | `docs/qr/src__mod_photofragment_flux.f90.png` |
| FSSH | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_hopping_fssh.f90> | `docs/qr/src__mod_surface_hopping_fssh.f90.png` |
| Jacobi 坐标 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_triatomic_geometry.f90> | `docs/qr/src__mod_triatomic_geometry.f90.png` |
| 超球反应 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_hyperspherical_reactive.f90> | `docs/qr/src__mod_hyperspherical_reactive.f90.png` |
| FGH 示例 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex01_fgh_diatomic_bound_states.f90> | `docs/qr/examples__ex01_fgh_diatomic_bound_states.f90.png` |
| Split 示例 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex03_split_operator_1d.f90> | `docs/qr/examples__ex03_split_operator_1d.f90.png` |
| TI/TD 示例 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex07_scattering_wavefunctions_ti_td.f90> | `docs/qr/examples__ex07_scattering_wavefunctions_ti_td.f90.png` |
| Feshbach 示例 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex08_ultracold_feshbach_segmented.f90> | `docs/qr/examples__ex08_ultracold_feshbach_segmented.f90.png` |
| FSSH 示例 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex30_tully_surface_hopping.f90> | `docs/qr/examples__ex30_tully_surface_hopping.f90.png` |
| Python 可视化 | <https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py> | `docs/qr/python__pygenmod__visualizer.py.png` |

# 参考文献

1. J. J. Sakurai and J. Napolitano, *Modern Quantum Mechanics*, 3rd ed., Cambridge University Press, Cambridge, 2020.
2. C. Cohen-Tannoudji, B. Diu, and F. Laloë, *Quantum Mechanics*, Wiley, New York, 1977.
3. D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific, Singapore, 1988. DOI: 10.1142/0270.
4. E. B. Wilson, Jr., J. C. Decius, and P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra*, McGraw-Hill, New York, 1955.
5. P. R. Bunker and P. Jensen, *Molecular Symmetry and Spectroscopy*, 2nd ed., NRC Research Press, Ottawa, 1998.
6. G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, Baltimore, 2013.
7. J. H. Wilkinson, *The Algebraic Eigenvalue Problem*, Clarendon Press, Oxford, 1965.
8. D. T. Colbert and W. H. Miller, “A novel discrete variable representation for quantum mechanical reactive scattering via the S-matrix Kohn method”, *J. Chem. Phys.* **96**, 1982 (1992). DOI: 10.1063/1.462125.
9. C. C. Marston and G. G. Balint-Kurti, “The Fourier grid Hamiltonian method for bound state eigenvalues and eigenfunctions”, *J. Chem. Phys.* **91**, 3571 (1989). DOI: 10.1063/1.456888.
10. D. Baye and P.-H. Heenen, “Generalised meshes for quantum mechanical problems”, *J. Phys. B: At. Mol. Opt. Phys.* **19**, 1991 (1986). DOI: 10.1088/0022-3700/19/14/006.
11. M. D. Feit, J. A. Fleck, Jr., and A. Steiger, “Solution of the Schrödinger equation by a spectral method”, *J. Comput. Phys.* **47**, 412 (1982). DOI: 10.1016/0021-9991(82)90091-2.
12. R. Kosloff, “Time-dependent quantum-mechanical methods for molecular dynamics”, *J. Phys. Chem.* **92**, 2087 (1988). DOI: 10.1021/j100319a003.
13. V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, “Completely positive dynamical semigroups of N-level systems”, *J. Math. Phys.* **17**, 821 (1976). DOI: 10.1063/1.522979.
14. G. Lindblad, “On the generators of quantum dynamical semigroups”, *Commun. Math. Phys.* **48**, 119 (1976). DOI: 10.1007/BF01608499.
15. H.-P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems*, Oxford University Press, Oxford, 2002.
16. R. Somlói, J. Kazakov, and D. J. Tannor, “A generalized relaxation method for optimal control of molecular motion”, *Chem. Phys.* **172**, 85 (1993). DOI: 10.1016/0301-0104(93)80108-L.
17. D. M. Reich, M. Ndong, and C. P. Koch, “Monotonically convergent optimal control theory of quantum systems”, *J. Chem. Phys.* **136**, 104103 (2012). DOI: 10.1063/1.3691827.
18. B. R. Johnson, “The multichannel log-derivative method for treating reactive collisions”, *J. Comput. Phys.* **13**, 445 (1973). DOI: 10.1016/0021-9991(73)90049-1.
19. D. E. Manolopoulos, “An improved log-derivative method for solving the radial Schrödinger equation”, *J. Chem. Phys.* **85**, 6425 (1986). DOI: 10.1063/1.451472.
20. H. Feshbach, “Unified theory of nuclear reactions”, *Ann. Phys. (N.Y.)* **5**, 357 (1958). DOI: 10.1016/0003-4916(58)90007-1.
21. C. Chin, R. Grimm, P. S. Julienne, and E. Tiesinga, “Feshbach resonances in ultracold gases”, *Rev. Mod. Phys.* **82**, 1225 (2010). DOI: 10.1103/RevModPhys.82.1225.
22. J. C. Tully, “Molecular dynamics with electronic transitions”, *J. Chem. Phys.* **93**, 1061 (1990). DOI: 10.1063/1.459170.
23. C. R. Harris et al., “Array programming with NumPy”, *Nature* **585**, 357 (2020). DOI: 10.1038/s41586-020-2649-2.
24. J. D. Hunter, “Matplotlib: A 2D graphics environment”, *Comput. Sci. Eng.* **9**, 90 (2007). DOI: 10.1109/MCSE.2007.55.
25. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>
