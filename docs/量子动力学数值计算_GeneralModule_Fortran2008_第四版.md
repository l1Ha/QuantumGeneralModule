---
title: 量子动力学数值计算
subtitle: 基于 GeneralModule 的 Fortran 2008 理论与实践
author: 李皓
email: LIH_ao@outlook.com
edition: 第四版
date: 2026-09-25
---

# 关于本书

本书围绕 GeneralModule 开源代码库，系统讨论量子动力学问题的数学表述、离散化方法、数值算法和工程实现。本书强调理论、算法与代码验证的完整链条。

本书对应的代码仓库为：

- <https://github.com/l1Ha/QuantumGeneralModule>

- 仓库二维码：

![](qr/repo_root.png){width=2.0cm}

# 符号与约定

除特别说明外，本书使用原子单位：

$$
\hbar=1,\qquad m_e=1,\qquad e=1,\qquad 4\pi\varepsilon_0=1.
$$

能量单位为 Hartree，长度单位为 Bohr。

# GeneralModule 第四版 · 扩展章节（第一部分）

## 量子力学数学结构、角动量代数、分子坐标理论与数值方法基础

本章为 GeneralModule 第四版新增的理论扩展章节，面向计算量子动力学、分子光谱学与超冷碰撞领域的研究生读者及算法开发者。全章沿四条主线递进展开：第一章论述希尔伯特空间与算符谱理论，这是有限基组与离散变量表象方法的数学根基；第二章系统建立角动量代数、Clebsch–Gordan 系数与 Wigner 3j/6j/9j 符号的重耦合理论；第三章讨论空间固定坐标系与体固定坐标系、Euler 角、Jacobi 坐标、超球坐标，以及 Watson 转振哈密顿量、Coriolis 耦合、离心畸变与振动角动量；第四章阐述单位制、实对称本征问题的 Householder–QL 算法、条件数、误差传播与收敛阶。

每一节均按统一体例组织：**定义**、**公式与推导**、**物理含义**、**数值陷阱**、**GeneralModule 实现映射**（源码路径、函数、GitHub 直链与二维码）、**小型数值实验**与**练习**。所有数值实验结果均由本库源码直接编译运行获得（GNU Fortran 11.4，`real64` 双精度，随机数种子取默认值），可复现、可作为读者自建验收测试的基线。行文约定：公式采用 Pandoc 兼容的 LaTeX 记号，行内公式以 $...$ 界定，独立公式以 $$...$$ 界定。

---

## 第一章 希尔伯特空间与算符谱理论

量子力学的全部数学结构建立在内积空间及其完备化之上。计算量子动力学的每一个数值算法——无论是有限基组展开、离散变量表象（DVR）还是波包传播——都是在希尔伯特空间的某个有限维截断子空间内工作。因此，严格理解无限维理论与有限维近似之间的对应与偏差，是判别数值结果可靠性的前提。本节内容由浅入深：先给出内积与完备性的基本定义，继而讨论基组展开与表象变换，然后陈述自伴算符的谱定理，最后落实到幺正演化的数值实现。

### 1.1 内积、范数与希尔伯特空间

**定义。** 设 $\mathcal{H}$ 为复数域 $\mathbb{C}$ 上的线性空间。若映射 $\langle\,\cdot\,|\,\cdot\,\rangle:\mathcal{H}\times\mathcal{H}\to\mathbb{C}$ 满足对第一变元共轭线性、对第二变元线性、共轭对称 $\langle\psi|\phi\rangle=\langle\phi|\psi\rangle^{*}$ 以及正定性 $\langle\psi|\psi\rangle\ge 0$（等号当且仅当 $\psi=0$），则称其为 $\mathcal{H}$ 上的内积。由内积诱导范数 $\|\psi\|=\sqrt{\langle\psi|\psi\rangle}$，并满足 Cauchy–Schwarz 不等式 $|\langle\psi|\phi\rangle|\le\|\psi\|\,\|\phi\|$ 与三角不等式。若 $\mathcal{H}$ 在该范数下完备，即任意 Cauchy 序列 $\{\psi_n\}$（满足 $\|\psi_n-\psi_m\|\to 0$）都收敛于 $\mathcal{H}$ 中元素，则称 $\mathcal{H}$ 为希尔伯特空间。量子力学标准模型进一步要求 $\mathcal{H}$ 可分，即存在可数的稠密子集。

**公式与推导。** 微观体系中最常见的两类实现为：有限维空间 $\mathbb{C}^{n}$，内积 $\langle u|v\rangle=\sum_i u_i^{*}v_i$；平方可积函数空间 $L^{2}(\mathbb{R}^{3})$，内积 $\langle\psi|\phi\rangle=\int_{\mathbb{R}^{3}}\psi^{*}(\mathbf{r})\phi(\mathbf{r})\,d^{3}\mathbf{r}$，其元素满足 $\int|\psi|^{2}\,d^{3}\mathbf{r}<\infty$。Cauchy–Schwarz 不等式的证明取 $\lambda\in\mathbb{C}$ 并利用 $0\le\|\psi+\lambda\phi\|^{2}=\|\psi\|^{2}+\lambda\langle\psi|\phi\rangle+\lambda^{*}\langle\phi|\psi\rangle+|\lambda|^{2}\|\phi\|^{2}$：当 $\phi\neq 0$ 时取 $\lambda=-\langle\phi|\psi\rangle/\|\phi\|^{2}$，判别式条件立即给出 $|\langle\psi|\phi\rangle|^{2}\le\|\psi\|^{2}\|\phi\|^{2}$。对散射问题，严格而言连续谱态（平面波、散射波）不属于 $L^{2}$，须借助装备希尔伯特空间 $\Phi\subset\mathcal{H}\subset\Phi^{\times}$（rigged Hilbert space）处理；数值计算中则以大盒子边界条件将连续谱离散化，使其纳入可分希尔伯特空间框架。

**物理含义。** 量子态是 $\mathcal{H}$ 中的矢量，测量概率由 Born 规则 $p(a)=\langle\psi|P_{a}|\psi\rangle$ 给出；两态之间的交叠 $\langle\psi|\phi\rangle$ 度量其相干叠加的强度，其模方在光谱学中对应跃迁线强度中的重叠因子。范数平方 $\|\psi\|^{2}=1$ 表述总概率守恒，是任何数值传播格式必须保持的第一守恒量。

**数值陷阱。** 其一，离散化后内积的度量约定不一：将 $L^{2}$ 内积 $\int\psi^{*}\phi\,dx$ 以梯形或矩形公式离散，得到 $\sum_i\psi_i^{*}\phi_i\,\Delta x$；而本征求解器（如本库 `fgh_solve_bound_states`）返回的本征矢量在"格点 Kronecker 度量" $\sum_i z_i^{*}z_i'=1$ 下归一化，物理波函数须按 $\psi(x_i)=z_i/\sqrt{\Delta x}$ 还原。混淆两种度量是转振计算中出现量级为 $\Delta x$ 的系统性错误的头号来源，3.4 节的数值实验将实际演示该错误（转动常数被缩小约 90 倍）。其二，长时间传播中归一化以 $\mathcal{O}(\varepsilon)$ 每步的速度漂移，须定期监测 $\langle\psi|\psi\rangle$。

**GeneralModule 实现映射。** 源码 `src/mod_dvr_grid.f90`：派生类型 `dvr_1d_t`（格点坐标 `x`、动能矩阵 `t_mat`）、`dvr_expectation_value`（以 $\sum_i|\psi_i|^{2}\,O(x_i)\,\Delta x$ 计算期望值）、`dvr_matrix_element`（两态跃迁矩阵元）。GitHub 直链：[src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90)。二维码：`qr/src__mod_dvr_grid.f90.png`。

![DVR 网格模块二维码](qr/src__mod_dvr_grid.f90.png)

**小型数值实验。** 取 Gauss 波包 $\psi(x)=\pi^{-1/4}\exp(-x^{2}/2)$，在区间 $[-10,10]$ 上以 Sinc-DVR 网格（$n=256$）离散并计算 $\sum_i|\psi(x_i)|^{2}\Delta x$，所得值为 $0.999999999997$，偏差 $3\times 10^{-12}$，与 $L^{2}$ 内积的求积误差一致；若省略权重 $\Delta x$，结果变为 $0.0323$，即偏差达两个数量级以上，直观演示度量约定的重要性。

**练习。** (1) 证明平行四边形定律 $\|\psi+\phi\|^{2}+\|\psi-\phi\|^{2}=2\|\psi\|^{2}+2\|\phi\|^{2}$，并说明内积空间的范数必然满足它。(2) 对三维谐振子基态波函数 $\psi(\mathbf{r})=(\alpha/\pi)^{3/4}e^{-\alpha r^{2}/2}$ 计算归一化常数并验证 $\langle\psi|\psi\rangle=1$。(3) 说明散射态为何不属于 $L^{2}$，并解释盒子归一化 $\psi\propto\sin(kr)$ 如何在 $L\to\infty$ 时恢复连续谱。

### 1.2 完备性、基组展开与表象变换

**定义。** 设 $\{\varphi_n\}_{n=1}^{\infty}\subset\mathcal{H}$ 为正交归一序列，即 $\langle\varphi_m|\varphi_n\rangle=\delta_{mn}$。若对任意 $\psi\in\mathcal{H}$ 均有 $\psi=\sum_n c_n\varphi_n$（按范数收敛），$c_n=\langle\varphi_n|\psi\rangle$，则称 $\{\varphi_n\}$ 为完备正交归一基。完备性等价于单位分解（resolution of identity）

$$\sum_{n=1}^{\infty}|\varphi_n\rangle\langle\varphi_n|=\mathbb{1},$$

也等价于 Parseval（Plancherel）等式 $\|\psi\|^{2}=\sum_n|c_n|^{2}$。给定基组即选定一个表象；基组间的幺正矩阵 $U_{mn}=\langle\varphi_m|\chi_n\rangle$ 实现表象变换，算符矩阵按 $A^{(\chi)}=U^{\dagger}A^{(\varphi)}U$ 变换。

**公式与推导。** 由单位分解立即得到 Parseval 等式：$\|\psi\|^{2}=\langle\psi|\mathbb{1}|\psi\rangle=\sum_n|\langle\varphi_n|\psi\rangle|^{2}$。变分原理给出截断近似的理论保障：对任意归一化 $\psi$，$\langle\psi|\hat{H}|\psi\rangle\ge E_0$（Ritz 变分原理），因此在 $N$ 维截断基中对角化 $\hat{H}$ 得到的最低本征值 $E_0^{(N)}$ 单调不增地自上方逼近精确基态能量；Hylleraas–Undheim–MacDonald 定理进一步保证第 $k$ 条能级 $E_k^{(N)}$ 随基组增大单调不增且与精确谱交错。对于 Sinc-DVR，Shannon 采样定理保证带限函数的展开是精确的，动能矩阵元具有解析形式

$$T_{ij}=\frac{\hbar^{2}}{2\mu\,\Delta x^{2}}
\begin{cases}
\dfrac{\pi^{2}}{3}, & i=j,\\[2mm]
\dfrac{2(-1)^{i-j}}{(i-j)^{2}}, & i\neq j,
\end{cases}$$

这正是本库 `dvr_sinc_init` 所实现的公式（Colbert–Miller, 1992）。

**物理含义。** 选取表象对应实验上的测量基选择：坐标表象波函数描述空间概率幅，能量表象展开系数描述定态布居，角动量表象描述取向分布。Fourier 网格哈密顿量（FGH）与 DVR 方法的物理内涵是：在正交格点上以"离散坐标本征态"近似连续坐标本征态，势能算符严格对角，动能算符解析非对角，二者相加后对角化即得束缚能级。

**数值陷阱。** 其一，截断基组的完备性只在变分子空间内成立：若势能面在网格外仍有可观振幅，能级不随网格增大收敛（基组套设错误，box error）；近解离能级尤其敏感，3.4 节实验中 Morse 势在盒子 $[0.4,6.0]\,a_0$ 内只得 15 条束缚能级，而半经典估计为 17 条，即为例证。其二，Fourier 网格存在 Nyquist 混叠：波数 $|k|>\pi/\Delta x$ 的分量被折叠为低频伪影，故网格间距必须满足 $\Delta x\le\pi/k_{\max}$。其三，若动能矩阵元以数值微分构造而非解析式，破坏了厄米性（$T_{ij}\neq T_{ji}^{*}$），变分上界性质随之失效。

**GeneralModule 实现映射。** 源码 `src/mod_dvr_grid.f90`：`dvr_sinc_init`（Colbert–Miller Sinc-DVR 动能矩阵）、`fgh_solve_bound_states`（FGH 束缚态求解，组装 $H_{ij}=T_{ij}+V(x_i)\delta_{ij}$ 并调用对称对角化）、`dvr_legendre_init`（Gauss–Legendre 角向 DVR）。文献依据见 Colbert 与 Miller（J. Chem. Phys. 96, 1982 (1992), DOI: 10.1063/1.462125）以及 Marston 与 Balint-Kurti（J. Chem. Phys. 91, 3571 (1989), DOI: 10.1063/1.456888）。GitHub 直链：[src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90)。二维码：`qr/src__mod_dvr_grid.f90.png`。

**小型数值实验。** 以 1.1 节的同一线型分子势（见 3.4 节参数）考察网格收敛性：将网格点数从 $n=250$ 增至 $n=500$，基态能量 $E_0$ 的变化小于 $10^{-9}\ E_h$，表明低能级已收敛；而第 15 条能级随盒子从 $6.0\ a_0$ 扩至 $8.0\ a_0$ 下移约 $4\ \mathrm{cm^{-1}}$，显示近解离能级的盒子误差。此实验说明"逐条能级、逐个参数地检验收敛"是变分计算不可省略的工序。

**练习。** (1) 由单位分解证明 Parseval 等式。(2) 证明 Sinc 函数族 $\mathrm{sinc}[\pi(x-x_i)/\Delta x]$ 在带限函数空间内正交归一。(3) 对宽度为 $L$ 的无限深势阱，估计用 Sinc-DVR 重现前 10 条能级所需的最小格距。

### 1.3 自伴算符与谱定理

**定义。** 稠密定义域 $\mathcal{D}(A)$ 上的线性算符 $A$ 称为对称的（厄米的），若 $\langle\psi|A\phi\rangle=\langle A\psi|\phi\rangle$ 对所有 $\psi,\phi\in\mathcal{D}(A)$ 成立；若进一步满足 $\mathcal{D}(A)=\mathcal{D}(A^{\dagger})$，则称 $A$ 为自伴的（self-adjoint）。在有限维空间中二者无区别，但在无限维空间中定义域之差是本质的。自伴算符的谱由离散谱（孤立有限重数本征值）、连续谱与残余谱组成；谱定理断言：自伴算符 $A$ 幺正等价于其谱测度，即

$$A=\int_{\sigma(A)}\lambda\,dE_{\lambda},\qquad \mathbb{1}=\int_{\sigma(A)}dE_{\lambda},$$

其中 $\{E_{\lambda}\}$ 为投影值测度；对纯离散谱退化为 $A=\sum_k\lambda_k|k\rangle\langle k|$。

**公式与推导。** 有限维谱定理可用 Rayleigh 商的极值刻画（Courant–Fischer 极小极大原理）：$E_0=\min_{\psi\neq0}\langle\psi|A|\psi\rangle/\langle\psi|\psi\rangle$，取极小的 $\psi$ 满足变分方程 $\delta\{\langle\psi|A|\psi\rangle-\lambda\langle\psi|\psi\rangle\}=0$，即 $A\psi=\lambda\psi$；对子空间逐级重复即得全部本征对。自伴性给出三条直接推论：本征值全为实数（$\lambda\langle\psi|\psi\rangle=\langle\psi|A|\psi\rangle\in\mathbb{R}$）；属于不同本征值的本征矢正交；任意对称实矩阵可由正交矩阵对角化 $A=Z\Lambda Z^{T}$。反例说明定义域的重要性：半直线 $[0,\infty)$ 上的动量算符 $-i\hbar\,d/dx$ 在满足边界条件的定义域上对称但不自伴，其本征值可任意复化，不存在物理可观测对应的谱分解。

**物理含义。** 可观测量与自伴算符一一对应：能量、角动量分量、偶极矩均为自伴算符，故其测量值实、本征态可构成正交完备基、含时演化保概率。数值上，DVR 或基组方法组装的哈密顿矩阵正是自伴算符在有限维子空间的投影，其对角化结果（能级与波函数）构成一切后续光谱与动力学分析的原材料。复吸收势（CAP）则刻意将哈密顿变为非厄米以吸收出射流，此时谱不再全实，须以非厄米谱理论解释（参见本库 `mod_absorbing_boundary.f90`）。

**数值陷阱。** 其一，矩阵组装代码只填上三角或只填下三角而忘记镜像，得到形式上"对称"实则不对称的输入，`tred2`/`tql2` 假定对称性后给出的结果不可预测；验收办法是显式检验 $\max_{ij}|A_{ij}-A_{ji}|=0$。其二，简并本征值对应的本征矢在对角化后是简并子空间内任意正交基，依赖于舍入噪声，跨程序比对本征矢无意义，比对投影算符 $\sum_{k\in\text{简并}}|k\rangle\langle k|$ 才有意义。其三，投影矩阵 $A=PAP^{\dagger}$ 仍是自伴的，但若投影后基组不正交（如重叠矩阵未同时处理），广义本征问题 $Hc=ESc$ 被误当作标准问题对角化，能级系统性偏离。

**GeneralModule 实现映射。** 源码 `src/mod_linear_algebra.f90`：`diag_symmetric_matrix`（三对角化 `tred2` 加 QL 迭代 `tql2` 的统一入口，输出升序本征值 $d$ 与正交本征矢矩阵 $z$）、`inv_real_matrix`/`inv_complex_matrix`（全主元 Gauss–Jordan 求逆，`stat=-1` 指示奇异矩阵）。GitHub 直链：[src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90)。二维码：`qr/src__mod_linear_algebra.f90.png`。

![线性代数模块二维码](qr/src__mod_linear_algebra.f90.png)

**小型数值实验。** 取 $n=300$ 的随机实对称矩阵（元素均匀分布于 $[-1,1]$），调用 `diag_symmetric_matrix` 后逐列计算残余 $\|A z_k-\lambda_k z_k\|_{\infty}$ 与正交性偏差 $\max_{ij}|(Z^{T}Z-\mathbb{1})_{ij}|$。实测结果：最大残余 $4.4\times10^{-14}$，最大正交性偏差 $2.2\times10^{-14}$。二者均为 $\mathcal{O}(n\varepsilon\|A\|)$ 量级（$\varepsilon=2.2\times10^{-16}$，$n\varepsilon\approx6.6\times10^{-14}$），符合舍入误差的预期标度，说明算法达到后向稳定。

**练习。** (1) 证明有限维空间中对称算符必然自伴。(2) 验证半直线上动量算符的对称性并找出使其非自伴的边界条件缺口。(3) 对 $n=10$ 的随机对称矩阵计算投影算符残差 $\|\sum_k|k\rangle\langle k|-\mathbb{1}\|_{F}$ 并解释其量级。

### 1.4 幺正演化与波包传播

**定义。** 算符 $U$ 称为幺正的，若 $U^{\dagger}U=U U^{\dagger}=\mathbb{1}$。含时薛定谔方程

$$i\hbar\frac{\partial}{\partial t}|\psi(t)\rangle=\hat{H}(t)|\psi(t)\rangle$$

在 $\hat{H}$ 不含时时的形式解为 $|\psi(t)\rangle=U(t)|\psi(0)\rangle$，$U(t)=\exp(-i\hat{H}t/\hbar)$；Stone 定理断言：强连续单参数幺正群 $U(t)$ 的生成元必为自伴算符，反之亦然。幺正性使 $\|\psi(t)\|$ 守恒，对应概率守恒。

**公式与推导。** 数值传播的核心是对 $U(\Delta t)$ 的高精度近似。Baker–Campbell–Hausdorff（BCH）展开给出对称分裂（Trotter–Strang 分裂）

$$e^{-i(\hat{T}+\hat{V})\Delta t/\hbar}
=e^{-i\hat{V}\Delta t/2\hbar}\,e^{-i\hat{T}\Delta t/\hbar}\,e^{-i\hat{V}\Delta t/2\hbar}
+\mathcal{O}(\Delta t^{3}),$$

误差首项为 $-\frac{i\Delta t^{3}}{24\hbar^{3}}[\hat{V},[\hat{V},\hat{T}]]+\frac{i\Delta t^{3}}{24\hbar^{3}}[\hat{T},[\hat{T},\hat{V}]]$ 型双对易子。由于三个因子均为幺正，对称分裂格式对任意步长都严格幺正——这不是近似而是精确性质，因此范数在机器精度内守恒；其三阶局部误差对束缚体系表现为相位误差而非振幅误差，长时间不产生久期增长。在坐标网格上，$\hat{V}$ 对角、$\hat{T}$ 在动量表象对角，故每次推进只需两次快速傅里叶变换；对大步长与长时间演化，Chebyshev 多项式展开 $e^{-i\hat{H}t/\hbar}\approx\sum_{n}a_n T_n(\hat{H}_{\mathrm{norm}})$ 提供谱收敛的替代方案（Kosloff, 1988）。

**物理含义。** 波包传播直接给出光解通量、散射矩阵、含时对齐度等可观测量；幺正性对应封闭体系的概率守恒，任何非幺正格式（如未加修正的显式欧拉法）都会在数千步内累积出物理上无意义的增益或衰减。分裂算符格式的辛对称性还保证了相空间体积守恒的量子对应——能级布居在长时间平均下不漂移。

**数值陷阱。** 其一，FFT 隐含周期边界：波包到达网格边缘后从另一侧回卷（wrap-around），须以足够大的盒子或复吸收势（本库 `mod_absorbing_boundary.f90`）消除。其二，时间步长须同时满足分裂误差与Nyquist 条件：$\Delta t$ 过大时高动量分量的相位 $e^{-iT(p)\Delta t/\hbar}$ 采样不足，产生高频混叠。其三，尽管格式严格幺正，浮点舍入仍以每步 $\mathcal{O}(\varepsilon)$ 引入微小的范数漂移，万步量级后可观测；宜在传播中周期性检查 $\langle\psi|\psi\rangle$ 并记录。

**GeneralModule 实现映射。** 源码 `src/mod_wavepacket_propagator.f90`：`propagate_split_operator_1d`（坐标半步势能相位、FFT 至动量空间、动能整步相位、逆变换、势能后半步，内部调用 `mod_linear_algebra.f90` 的 `fft_1d`）、`propagate_split_operator_2d`、`rk4_step`、`abm4_step`；`src/mod_chebyshev_propagator.f90` 提供切比雪夫大步长推进器。文献依据：Feit, Fleck 与 Steiger（J. Comput. Phys. 47, 412 (1982), DOI: 10.1016/0021-9991(82)90091-2）、Kosloff（J. Phys. Chem. 92, 2087 (1988), DOI: 10.1021/j100319a003）。GitHub 直链：[src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)。二维码：`qr/src__mod_wavepacket_propagator.f90.png`。

![波包传播模块二维码](qr/src__mod_wavepacket_propagator.f90.png)

**小型数值实验。** 一维谐振子（$m=1$，$\omega=0.02\ E_h/\hbar$），网格 $n=1024$、$\Delta x=0.05\ a_0$、$\Delta t=0.5\ a.u.$，以高斯波包为初态连续推进 $20000$ 步（相当于 $10000\ a.u.\approx242\ \mathrm{fs}$）。实测终态范数偏差 $|\ \|\psi\|^{2}-1\ |=8.8\times10^{-12}$，确认了格式的严格幺正性与舍入漂移的 $\mathcal{O}(N_{\mathrm{step}}\varepsilon)$ 标度。

**练习。** (1) 用 BCH 展开推导对称分裂的 $\mathcal{O}(\Delta t^{3})$ 误差项并证明其为厄米算符的虚倍数。(2) 证明 $e^{-i\hat{V}\Delta t/2\hbar}e^{-i\hat{T}\Delta t/\hbar}e^{-i\hat{V}\Delta t/2\hbar}$ 对任意 $\Delta t$ 幺正。(3) 对谐振子精确解比较分裂格式的相位误差随 $\Delta t$ 的标度，并估计达到 $10^{-6}$ 相位精度所需步长。
# 第二部分　谱方法、离散表示与收敛理论

## 2.1　从基组展开到矩阵问题

设哈密顿量为

$$
\hat H=\hat T+\hat V.
$$

将态展开为

$$
|\psi\rangle=\sum_{n=1}^{N}c_n|\chi_n\rangle.
$$

代入薛定谔方程并在基组上投影，得到广义本征值问题

$$
\sum_j H_{ij}c_j
=
E\sum_j S_{ij}c_j,
$$

其中

$$
H_{ij}=\langle \chi_i|\hat H|\chi_j\rangle,
\qquad
S_{ij}=\langle \chi_i|\chi_j\rangle.
$$

若基组正交，则 $S_{ij}=\delta_{ij}$，问题退化为标准本征值问题

$$
H\mathbf c=E\mathbf c.
$$

若基组非正交但线性无关，可以对重叠矩阵作 Cholesky 分解

$$
S=LL^T,
$$

并定义

$$
\tilde H=L^{-1}H L^{-T}.
$$

广义问题变为对称标准问题

$$
\tilde H\mathbf y=E\mathbf y,
\qquad
\mathbf c=L^{-T}\mathbf y.
$$

数值上要避免直接求 $S^{-1}$，因为当基组近线性相关时，$S$ 的最小本征值接近零，$\kappa(S)$ 会急剧增大。

## 2.2　DVR 的基本结构

离散变量表象不是一种具体算法，而是一类满足以下条件的表示：

1. 网格点 $\{x_i\}_{i=1}^{N}$ 对应正交函数或正交多项式；
2. 离散内积由权重 $\{w_i\}$ 定义；
3. 局域算符在网格上近似对角；
4. 动能矩阵由解析或半解析公式给出。

离散内积为

$$
\langle f|g\rangle_N
=
\sum_{i=1}^{N}w_i f(x_i)^{*}g(x_i).
$$

对局域势能算符，

$$
\langle f|\hat V|g\rangle
=
\int f^{*}(x)V(x)g(x)\,dx
\approx
\sum_i w_i f(x_i)^{*}V(x_i)g(x_i),
$$

因此

$$
V_{ij}=V(x_i)\delta_{ij}.
$$

这是 DVR 的核心优势。它把高维问题中的势能矩阵由稠密矩阵变成对角矩阵，使矩阵存储从 $O(N^2)$ 降至 $O(N)$。

## 2.3　Fourier 网格与 FGH

均匀周期域 $x\in[0,L)$ 上的 Fourier 基为

$$
\phi_k(x)=\frac{1}{\sqrt L}e^{ikx},
\qquad
k=\frac{2\pi n}{L},
$$

$$
n=-\frac{N}{2},\ldots,\frac{N}{2}-1.
$$

动量算符满足

$$
\hat p\phi_k(x)=\hbar k\phi_k(x).
$$

因此动能矩阵在动量空间对角：

$$
T_k=\frac{\hbar^2k^2}{2m}.
$$

势能在坐标网格上近似为

$$
V_{ij}=V(x_i)\delta_{ij}.
$$

这构成 Fourier Grid Hamiltonian 方法。Marston 和 Balint-Kurti 给出了 FGH 在分子束缚态计算中的系统表述【Marston & Balint-Kurti 1989】。

动量分辨率和最大可表示动量为

$$
\Delta k=\frac{2\pi}{L},
\qquad
k_{\max}=\frac{\pi}{\Delta x}.
$$

若波函数动量宽度为 $\sigma_k$，为避免混叠通常要求

$$
k_{\max}\gtrsim 6\sigma_k.
$$

空间分辨率应满足

$$
\Delta x
\ll
\frac{1}{k_{\max}^{\mathrm{local}}}.
$$

在数值实验中，最直接的检验是把 $N$ 增加一倍，观察本征值和期望值是否稳定。

## 2.4　Sinc-DVR 的推导与误差结构

对均匀网格

$$
x_i=x_{\min}+i\Delta x,
\qquad
i=1,\ldots,N,
$$

Colbert 和 Miller 给出的 Sinc-DVR 动能矩阵元为

$$
T_{ij}
=
\frac{\hbar^2}{2m\Delta x^2}
\begin{cases}
\pi^2/3, & i=j,\\[4pt]
2(-1)^{i-j}/(i-j)^2, & i\ne j.
\end{cases}
$$

原子单位下 $\hbar=1$，因此

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
H_{ij}=T_{ij}+V(x_i)\delta_{ij}.
$$

若离散内积权重为 $w_i=\Delta x$，则连续归一化波函数为

$$
\psi(x_i)=\frac{c_i}{\sqrt{\Delta x}}.
$$

### 误差来源

Sinc-DVR 的误差主要来自四类：

1. 空间截断：计算域太窄导致边界反射；
2. 网格欠采样：波函数高频分量超过 $k_{\max}$；
3. 势能奇点：局部梯度过大；
4. 非厄米离散：动能和势能使用了不一致的内积。

对光滑势函数，谱方法通常呈指数收敛：

$$
E_N-E\sim C_Ne^{-\alpha N}.
$$

若势能不可微或存在奇点，收敛可能退化为代数阶：

$$
E_N-E\sim CN^{-p}.
$$

因此，不能仅凭两次计算结果接近判断收敛，应至少使用三个网格密度做 Richardson 外推。

## 2.5　Legendre-DVR 与角向表示

分子取向问题常采用

$$
x=\cos\theta.
$$

Legendre 多项式满足

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

在角坐标中，

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

在 DVR 网格上，它变为矩阵 $J^2_{ij}$。GeneralModule 的 `dvr_legendre_init` 计算 Gauss-Legendre 零点、权重和 $J^2$ 矩阵。

### 对应源码

- [src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90)
- 二维码：`qr/src__mod_dvr_grid.f90.png`

![](qr/src__mod_dvr_grid.f90.png){width=2.0cm}

## 2.6　DVR 类型的适用范围

DVR 不是单一算法。不同 DVR 的差异在于权重、动能矩阵和边界条件。

| DVR 类型 | 坐标 | 动能矩阵 | 适用问题 | 主要风险 |
|---|---|---|---|---|
| Sinc-DVR | 等距 $x$ | 解析长程衰减矩阵 | 一维局域势能面 | 边界反射 |
| Fourier/FGH | 周期 $x$ | 动量空间对角 | 周期或大范围波包 | 周期性假设 |
| Gauss-Legendre DVR | $x=\cos\theta$ | 角动量矩阵 | 分子取向与转子 | 端点处理 |
| Gauss-Hermite DVR | 无限域 | 谐振子结构 | 振动模式 | 势能范围超出权重衰减 |
| Laguerre-DVR | 半无限径向域 | 径向结构 | 长程势 | 奇点处理 |
| Distributed Gaussian DVR | 可调局域高斯 | 动能重叠矩阵 | 多维非均匀势能面 | 线性相关 |

因此，将 DVR 等同于 Sinc-DVR 是不完整的。GeneralModule 提供 Sinc-DVR 与 Gauss-Legendre DVR，并可通过同一框架扩展到其他正交多项式网格。

## 2.7　多维张量积 DVR

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

### 复杂度

若每个维度有 $N$ 个点，$d$ 维张量积网格大小为 $N^d$。稠密对角化复杂度为

$$
O(N^{3d}),
$$

存储复杂度为

$$
O(N^{2d}).
$$

因此，直接对角化只适用于低维或截断后的小问题。高维问题需要稀疏矩阵、Krylov 方法、Chebyshev 传播、张量分解或缩减基方法。

## 2.8　收敛性与 Richardson 外推

设精确本征值为 $E$，数值本征值为 $E_N$。若误差服从

$$
E_N-E
=
CN^{-p}+O(N^{-(p+1)}),
$$

则三种网格 $N_1,N_2,N_3$ 可估计阶数 $p$：

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

### 数值实验

1. 对 Morse 势分别取 $N=128,256,512,1024$，绘制前六个能级误差。
2. 固定 $\Delta x$，改变区间长度，观察高激发态边界反射。
3. 将 Morse 势替换为谐振子势，比较 Sinc-DVR 和解析谱。
4. 用 Legendre-DVR 计算刚性转子 $J=0,1,2$ 能级，并与 $BJ(J+1)$ 比较。

## 2.9　从矩阵结构到算法选择

在 DVR 中，哈密顿量矩阵的结构决定算法：

| 矩阵结构 | 算法 | 复杂度 |
|---|---|---|
| 中小型稠密对称矩阵 | Householder + QL | $O(N^3)$ |
| 大型稀疏矩阵 | Lanczos / Arnoldi | 接近 $O(kN)$ |
| FFT 可分离动能 | Split-Operator | $O(N\log N)$ 每步 |
| 大步长谱传播 | Chebyshev 展开 | $O(kN)$ 每步 |
| 高维张量结构 | Tensor train / MCTDH | 问题相关 |

GeneralModule 的 `mod_linear_algebra` 提供自包含稠密对称本征求解和 FFT；`mod_dvr_grid` 提供网格构造；`mod_wavepacket_propagator` 和 `mod_chebyshev_propagator` 提供时间演化。三者构成一个完整的谱方法教学与验证链路。

## 2.10　Python 交叉验证

`pygenmod` 使用 NumPy 和 Matplotlib 对 DVR、散射和谱方法进行独立实现与可视化。这有两个作用：

1. 用另一种语言检验 Fortran 数值实现；
2. 将 Fortran 输出转化为发表级图形。

NumPy 的数组操作和 FFT 实现可作为参考算法，Matplotlib 用于绘制能级、波函数和误差曲线【Harris et al. 2020；Hunter 2007】。此类交叉验证不是替代 Fortran 主计算，而是数值验证工具。

### 对应源码

- [python/pygenmod/visualizer.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py)
- 二维码：`qr/python__pygenmod__visualizer.py.png`

![](qr/python__pygenmod__visualizer.py.png){width=2.0cm}

## 2.11　参考文献

1. D. T. Colbert and W. H. Miller, “A novel discrete variable representation for quantum mechanical reactive scattering via the S-matrix Kohn method”, *J. Chem. Phys.* **96**, 1982 (1992). DOI: 10.1063/1.462125.
2. C. C. Marston and G. G. Balint-Kurti, “The Fourier grid Hamiltonian method for bound state eigenvalues and eigenfunctions”, *J. Chem. Phys.* **91**, 3571 (1989). DOI: 10.1063/1.456888.
3. D. Baye and P.-H. Heenen, “Generalised meshes for quantum mechanical problems”, *J. Phys. B: At. Mol. Opt. Phys.* **19**, 1991 (1986). DOI: 10.1088/0022-3700/19/14/006.
4. J. C. Light and T. Carrington, Jr., “Discrete-variable representations and their utilization in quantum dynamics calculations”, *Adv. Chem. Phys.* **114**, 263 (2000). DOI: 10.1002/9780470141731.ch5.
5. C. R. Harris et al., “Array programming with NumPy”, *Nature* **585**, 357 (2020). DOI: 10.1038/s41586-020-2649-2.
6. J. D. Hunter, “Matplotlib: A 2D graphics environment”, *Comput. Sci. Eng.* **9**, 90 (2007). DOI: 10.1109/MCSE.2007.55.
7. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>
# 第三部分　含时演化、开放系统与量子控制

## 3.1　时间演化与守恒结构

含时薛定谔方程为

$$
i\hbar\frac{\partial}{\partial t}|\psi(t)\rangle
=
\hat H(t)|\psi(t)\rangle.
$$

若哈密顿量不含时，则形式解为

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

若哈密顿量含时，则

$$
\hat U(t,t_0)
=
\mathcal{T}
\exp\left[
-\frac{i}{\hbar}
\int_{t_0}^{t}\hat H(t')\,dt'
\right],
$$

其中 $\mathcal{T}$ 是时间排序算符。数值方法的任务是用有限步近似该时间排序指数。

对不含时厄米哈密顿量，

$$
\hat U^{\dagger}\hat U=I,
$$

即

$$
\langle\psi(t)|\psi(t)\rangle
=
\langle\psi(t_0)|\psi(t_0)\rangle.
$$

数值传播器必须保持或至少诊断该性质。

### 可观测量

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

若使用动量空间波函数 $\tilde\psi(k,t)$，则

$$
\langle p\rangle(t)
=
\int \hbar k\,|\tilde\psi(k,t)|^2\,dk.
$$

## 3.2　Split-Operator 方法

对于不含时哈密顿量

$$
\hat H=\hat T+\hat V,
$$

若 $\hat T$ 与 $\hat V$ 不对易，则

$$
e^{-i(\hat T+\hat V)\Delta t/\hbar}
\ne
e^{-i\hat T\Delta t/\hbar}
e^{-i\hat V\Delta t/\hbar}.
$$

利用 Baker–Campbell–Hausdorff 展开可得一阶分裂

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

Feit、Fleck 和 Steiger 将该谱方法系统用于薛定谔方程求解【Feit, Fleck & Steiger 1982】。Kosloff 对含时量子动力学中的谱传播方法作了系统综述【Kosloff 1988】。

在坐标表象中，$\hat V$ 是局域算符：

$$
[\hat V\psi](x)=V(x)\psi(x).
$$

在动量表象中，$\hat T=p^2/(2m)$ 是对角算符：

$$
[\hat T\tilde\psi](k)=\frac{\hbar^2 k^2}{2m}\tilde\psi(k).
$$

因此，一步传播为

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

在原子单位下，$\hbar=1$，公式简化为

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

### 数值步骤

1. 在坐标网格上计算势能半步相位 $e^{-iV(x_i)\Delta t/2}$；
2. 将波函数乘以该相位；
3. 使用 FFT 变换到动量空间；
4. 乘以动能相位 $e^{-ik^2\Delta t/(2m)}$；
5. 使用逆 FFT 回到坐标空间；
6. 再次乘以势能半步相位。

每步需要两次 FFT，复杂度为 $O(N\log N)$，显著优于稠密矩阵指数的 $O(N^3)$。

### 范数与稳定性

每个指数相位因子都是幺正算符。因此，在无穷精度下 Split-Operator 严格保持范数：

$$
\|\psi(t+\Delta t)\|_2=\|\psi(t)\|_2.
$$

有限精度误差来源包括：

- FFT 归一化不一致；
- 网格混叠；
- 复势导致非幺正演化；
- 吸收边界；
- 长时间积累的舍入误差；
- 不当归一化掩盖真实误差。

因此，应监测

$$
\epsilon_{\mathrm{norm}}(t)
=
\left|
\|\psi(t)\|_2-1
\right|,
$$

而不是在每步强制归一化。

### 对应源码

- [src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)
- 二维码：`qr/src__mod_wavepacket_propagator.f90.png`

![](qr/src__mod_wavepacket_propagator.f90.png){width=2.0cm}

## 3.3　FFT、动量网格与混叠

均匀网格 $x_i=x_0+i\Delta x$，$i=0,\ldots,N-1$ 对应动量网格

$$
k_j
=
\frac{2\pi}{L}
\begin{cases}
j, & 0\le j<N/2,\\
j-N, & N/2\le j<N,
\end{cases}
$$

其中 $L=N\Delta x$。这样，$k$ 覆盖

$$
-\frac{\pi}{\Delta x}
\le
k
<
\frac{\pi}{\Delta x}.
$$

离散傅里叶变换为

$$
\tilde\psi(k_j)
=
\sum_{n=0}^{N-1}
\psi(x_n)e^{-ik_jx_n}\Delta x,
$$

逆变换为

$$
\psi(x_n)
=
\frac{1}{L}
\sum_{j=0}^{N-1}
\tilde\psi(k_j)e^{ik_jx_n}.
$$

若波函数在动量空间的有效宽度超过 $|k|<\pi/\Delta x$，则发生混叠：

$$
k_{\mathrm{alias}}
=
k-2m\pi/\Delta x,
\qquad
m\in\mathbb{Z}.
$$

避免混叠的条件是

$$
\Delta x
<
\frac{\pi}{k_{\max}}.
$$

在数值实验中应通过逐步减小 $\Delta x$ 检查可观测量是否收敛。

## 3.4　Runge–Kutta 与多步方法

将波函数或密度矩阵写成常微分方程

$$
\frac{d\mathbf{y}}{dt}=f(t,\mathbf{y}),
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

其局部截断误差为 $O(\Delta t^5)$，全局误差为 $O(\Delta t^4)$。RK4 对非幺正演化或非线性常微分方程非常通用，但对长时幺正演化可能产生范数漂移。

Adams–Bashforth–Moulton 四阶预估—校正方法使用历史导数。Adams–Bashforth 预估式为

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
f_{n+1}^{(p)}
=
f(t_{n+1},y_{n+1}^{(p)}),
$$

然后使用 Adams–Moulton 校正式

$$
y_{n+1}
=
y_n
+
\frac{\Delta t}{24}
(9f_{n+1}^{(p)}+19f_n-5f_{n-1}+f_{n-2}).
$$

多步方法依赖历史点的一致性和启动策略。若历史导数来自不同时间步长或不同坐标系，阶数会退化。

### 对应源码

- [src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)
- 二维码：`qr/src__mod_wavepacket_propagator.f90.png`

![](qr/src__mod_wavepacket_propagator.f90.png){width=2.0cm}

### 练习

1. 将谐振子波包分别用 Split-Operator 和 RK4 传播，比较范数误差随 $\Delta t$ 的标度。
2. 构造非光滑右端项，比较 RK4 和 ABM4 的稳定性。
3. 改变 FFT 网格密度，观察动量空间尾部截断对 $\langle x\rangle$、$\langle p\rangle$ 的影响。

## 3.5　二能级系统与 Bloch 方程

二能级系统的态可写为

$$
|\psi(t)\rangle
=
c_1(t)|1\rangle+c_2(t)|2\rangle.
$$

在旋转波近似下，拉比频率为 $\Omega_R(t)$，失谐为 $\Delta(t)$，Bloch 矢量为

$$
\mathbf{R}
=
(u,v,w)^T.
$$

光学 Bloch 方程的无耗散部分为

$$
\frac{d\mathbf{R}}{dt}
=
\boldsymbol{\Omega}\times\mathbf{R},
$$

其中

$$
\boldsymbol{\Omega}
=
(\Omega_R,0,\Delta)^T.
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

$w=\rho_{22}-\rho_{11}$ 表示布居反转。若 $\Delta=0$ 且脉冲面积为

$$
\Theta
=
\int \Omega_R(t)\,dt,
$$

则完全反转为

$$
\Theta=\pi.
$$

该方程是检验 RK4、保范数传播和激光脉冲合成的基础模型。

### 对应源码

- [src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)
- 二维码：`qr/src__mod_wavepacket_propagator.f90.png`

![](qr/src__mod_wavepacket_propagator.f90.png){width=2.0cm}

## 3.6　复吸收边界

当波包离开有限计算域时，硬边界会产生非物理反射。复吸收势（Complex Absorbing Potential, CAP）在边界区域加入虚部势：

$$
V_{\mathrm{eff}}(x)
=
V(x)-iW(x),
$$

$$
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
\hat H_{\mathrm{eff}}
=
\hat T+\hat V-i\hat W.
$$

内部区域波函数的模平方不再守恒，其损失应与流出概率一致：

$$
\frac{d}{dt}
\int_{\Omega_{\mathrm{in}}}
|\psi|^2\,dx
=
-
\int_{\partial\Omega_{\mathrm{in}}}
J_n\,dS
-
\int_{\Omega_{\mathrm{CAP}}}
2W(x)|\psi|^2\,dx.
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

好的 CAP 参数应满足：

1. 内部物理区域几乎不受 CAP 影响；
2. 反射振幅足够小；
3. 吸收区长度和强度不过度增大计算域；
4. 吸收量与通量诊断一致。

### 对应源码

- [src/mod_absorbing_boundary.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_absorbing_boundary.f90)
- 二维码：`qr/src__mod_absorbing_boundary.f90.png`

![](qr/src__mod_absorbing_boundary.f90.png){width=2.0cm}

## 3.7　Chebyshev 传播与谱方法大步长

对不含时哈密顿量，可利用 Chebyshev 多项式展开演化算符：

$$
e^{-i\hat H t/\hbar}
=
\sum_{n=0}^{\infty}
a_n(t)T_n(\tilde{\hat H}),
$$

其中 $\tilde{\hat H}$ 是缩放到 $[-1,1]$ 的哈密顿量。若哈密顿量谱满足

$$
E_{\min}\le \hat H\le E_{\max},
$$

则令

$$
\tilde{\hat H}
=
\frac{2\hat H-(E_{\max}+E_{\min})I}
{E_{\max}-E_{\min}}.
$$

Chebyshev 展开系数为

$$
a_n(t)
=
(2-\delta_{n0})(-i)^n
J_n(\Delta E\,t/2\hbar),
$$

其中 $J_n$ 是第一类 Bessel 函数，$\Delta E=E_{\max}-E_{\min}$。所需项数大约为

$$
N_{\mathrm{Cheb}}
\sim
\frac{\Delta E\,t}{2\hbar}
+
c\ln\epsilon^{-1}.
$$

Chebyshev 方法适合大时间步和长时传播，但要求谱界估计准确。若谱外态存在，可能出现非物理放大。

### 对应源码

- [src/mod_chebyshev_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_chebyshev_propagator.f90)
- 二维码：`qr/src__mod_chebyshev_propagator.f90.png`

![](qr/src__mod_chebyshev_propagator.f90.png){width=2.0cm}

## 3.8　密度矩阵与开放量子系统

封闭系统的态可用纯态 $|\psi\rangle$ 描述，对应密度矩阵

$$
\hat\rho=|\psi\rangle\langle\psi|.
$$

一般混合态为

$$
\hat\rho
=
\sum_\alpha p_\alpha
|\psi_\alpha\rangle\langle\psi_\alpha|,
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

封闭系统的刘维尔方程为

$$
\frac{d\hat\rho}{dt}
=
-\frac{i}{\hbar}[\hat H,\hat\rho].
$$

与环境耦合后，马尔可夫近似下的 Lindblad 主方程为

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

该形式由 Gorini–Kossakowski–Sudarshan 和 Lindblad 给出【Gorini, Kossakowski & Sudarshan 1976；Lindblad 1976】。Breuer 和 Petruccione 对开放量子系统理论作了系统论述【Breuer & Petruccione 2002】。

### 完全正定性与迹守恒

Lindblad 形式保证

$$
\frac{d}{dt}\operatorname{Tr}\hat\rho=0,
$$

且演化映射完全正定。数值实现必须检查：

1. 迹守恒：
   $|\operatorname{Tr}\hat\rho-1|<\epsilon$；
2. 厄米性：
   $\|\hat\rho-\hat\rho^{\dagger}\|<\epsilon$；
3. 正定性：
   最小本征值不显著小于 $-\epsilon$；
4. 纯度单调性或物理耗散趋势。

### 常见耗散通道

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
\frac{d\rho_{ee}}{dt}
=
-\gamma\rho_{ee},
$$

$$
\frac{d\rho_{eg}}{dt}
=
-\left(i\omega_0+\frac{\gamma}{2}+\gamma_\phi\right)\rho_{eg}.
$$

因此，布居弛豫时间为 $1/\gamma$，相干衰减时间为

$$
T_2^{-1}
=
\frac{1}{2T_1}+\gamma_\phi.
$$

### 纯度、相干度与熵

纯度为

$$
P=\operatorname{Tr}(\hat\rho^2).
$$

纯态满足 $P=1$；最大混合态满足 $P=1/d$。冯·诺依曼熵为

$$
S=-\operatorname{Tr}(\hat\rho\ln\hat\rho).
$$

$l_1$ 范数相干度为

$$
C_{l1}(\hat\rho)
=
\sum_{i\ne j}|\rho_{ij}|.
$$

这三个量分别描述混合程度、非对角相干和信息含量。

### RK4 数值推进

将 Lindblad 方程写成矩阵值常微分方程

$$
\dot{\hat\rho}=F(\hat\rho),
$$

则 RK4 步为

$$
\hat\rho_{n+1}
=
\hat\rho_n
+
\frac{\Delta t}{6}
(K_1+2K_2+2K_3+K_4).
$$

每一步之后可检查厄米性和迹，但不应在参数研究开始前就强制归一化，否则会掩盖真实时间步长误差。

### 对应源码

- [src/mod_open_quantum.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_open_quantum.f90)
- 二维码：`qr/src__mod_open_quantum.f90.png`

![](qr/src__mod_open_quantum.f90.png){width=2.0cm}

## 3.9　量子最优控制与 Krotov 方法

量子控制的目标是设计外场 $\epsilon(t)$，使演化态在终端时刻逼近目标态。保真度为

$$
F
=
|\langle\phi_{\mathrm{target}}|\psi(T)\rangle|^2.
$$

最优控制问题可写成泛函

$$
J
=
F
-
\int_0^T
\lambda(t)\epsilon^2(t)\,dt.
$$

Krotov 方法引入伴随态 $|\chi(t)\rangle$，满足反向方程

$$
i\hbar\frac{\partial}{\partial t}|\chi(t)\rangle
=
\hat H(t)|\chi(t)\rangle,
$$

终端条件为

$$
|\chi(T)\rangle
=
|\phi_{\mathrm{target}}\rangle.
$$

场的更新形式可写为

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

### Krotov 单次迭代流程

1. 以当前场 $\epsilon^{(k)}(t)$ 正向传播 $|\psi(t)\rangle$；
2. 设置伴随态终端条件 $|\chi(T)\rangle$；
3. 反向传播 $|\chi(t)\rangle$；
4. 在每一时刻计算耦合项
   $\operatorname{Im}\langle\chi|\hat\mu|\psi\rangle$；
5. 更新 $\epsilon(t)$；
6. 用新场开始下一轮迭代。

收敛判据应同时包括保真度增量、场变化和物理约束。

### 对应源码

- [src/mod_optimal_control.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_optimal_control.f90)
- 二维码：`qr/src__mod_optimal_control.f90.png`

![](qr/src__mod_optimal_control.f90.png){width=2.0cm}

## 3.10　验证与收敛实验

含时量子动力学的验证应包含以下层次：

| 验证层次 | 检查量 | 目标 |
|---|---|---|
| 解析极限 | 自由粒子、谐振子、$\pi$ 脉冲 | 算法和单位正确 |
| 范数守恒 | $\|\psi(t)\|$ | 幺正性与 FFT 一致 |
| 可观测守恒 | $\langle H\rangle$、$\langle J\rangle$ | 对称性未被破坏 |
| 收敛性 | $E$、$T(E)$、$F$ 随 $N,\Delta t$ 变化 | 误差可控 |
| CAP 一致性 | 内部范数损失与流出通量 | 吸收边界不污染物理区 |
| 开放系统 | 迹、厄米性、正定性 | Lindblad 结构保持 |
| 控制收敛 | $F$、$\epsilon(t)$ 变化 | 最优控制收敛 |

一个基本数值实验是自由粒子波包。设初始高斯波包为

$$
\psi(x,0)
=
\left(\frac{1}{2\pi\sigma^2}\right)^{1/4}
\exp\left[
-\frac{(x-x_0)^2}{4\sigma^2}
+ik_0x
\right].
$$

无外场时，

$$
\langle x\rangle(t)
=
x_0+\frac{\hbar k_0}{m}t.
$$

数值传播应再现该线性运动和动量守恒。

### 对应测试

- [tests/test_propagators.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_propagators.f90)
- [tests/test_open_quantum_opt.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_open_quantum_opt.f90)
- 二维码：`qr/tests__test_propagators.f90.png`

![](qr/tests__test_propagators.f90.png){width=2.0cm}

## 3.11　参考文献

1. M. D. Feit, J. A. Fleck, Jr., and A. Steiger, “Solution of the Schrödinger equation by a spectral method”, *J. Comput. Phys.* **47**, 412 (1982). DOI: 10.1016/0021-9991(82)90091-2.
2. R. Kosloff, “Time-dependent quantum-mechanical methods for molecular dynamics”, *J. Phys. Chem.* **92**, 2087 (1988). DOI: 10.1021/j100319a003.
3. V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, “Completely positive dynamical semigroups of N-level systems”, *J. Math. Phys.* **17**, 821 (1976). DOI: 10.1063/1.522979.
4. G. Lindblad, “On the generators of quantum dynamical semigroups”, *Commun. Math. Phys.* **48**, 119 (1976). DOI: 10.1007/BF01608499.
5. H.-P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems*, Oxford University Press, Oxford, 2002.
6. R. Somlói, J. Kazakov, and D. J. Tannor, “A generalized relaxation method for optimal control of molecular motion”, *Chem. Phys.* **172**, 85 (1993). DOI: 10.1016/0301-0104(93)80108-L.
7. D. M. Reich, M. Ndong, and C. P. Koch, “Monotonically convergent optimal control theory of quantum systems”, *J. Chem. Phys.* **136**, 104103 (2012). DOI: 10.1063/1.3691827.
8. D. J. Tannor, *Introduction to Quantum Mechanics: A Time-Dependent Perspective*, University Science Books, Sausalito, 2007.
9. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>
# 第四部分　散射、分子动力学与复杂体系

## 4.1　定态散射的基本问题

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

在远区 $V(r)\rightarrow0$ 时，

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

## 4.2　Numerov 方法与散射长度

Numerov 方法适用于二阶常微分方程

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

GeneralModule 的 `calc_scattering_length_numerov` 实现该流程，并在远区用三点公式计算对数导数：

$$
u'(r_N)
\approx
\frac{3u_N-4u_{N-1}+u_{N-2}}{2\Delta r}.
$$

### 对应源码

- [src/mod_ti_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_ti_scattering.f90)
- 二维码：`qr/src__mod_ti_scattering.f90.png`

![](qr/src__mod_ti_scattering.f90.png){width=2.0cm}

## 4.3　Johnson 对数导数方法

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

得到。Johnson 方法在深势阱和长程势中数值稳定性更好。

## 4.4　相移与光学定理

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

总截面与向前散射振幅满足光学定理

$$
\sigma_{\mathrm{tot}}
=
\frac{4\pi}{k}
\operatorname{Im}f(0).
$$

该关系可用于检验散射振幅、相移和截面计算的一致性。

## 4.5　有效力程展开

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

GeneralModule 的 `fit_effective_range_expansion` 和 `gribakin_flambaum_length` 分别实现这两类低能分析。

## 4.6　多通道耦合与 Johnson 矩阵对数导数

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

## 4.7　分段网格与局部扇区

在长程势或强局域势问题中，单一均匀网格往往效率低。分段网格将径向区间划分为多个局部扇区：

$$
[r_{\min},r_{\max}]
=
\bigcup_{s=1}^{S}
[r_s,r_{s+1}].
$$

每个扇区可使用不同步长 $\Delta r_s$。这种方法在近区加密、远区稀疏，可以显著减少总网格数，同时保持误差可控。

GeneralModule 的 `segmented_grid_t` 存储每个扇区的起点、终点、步长、点数和全局索引。相关函数包括：

- `create_segmented_grid`
- `calc_scattering_length_segmented_numerov`
- `calc_scattering_wavefunction_segmented_ti`
- `calc_phase_shift_segmented`
- `calc_multichannel_close_coupling_segmented_logder`

## 4.8　Feshbach 共振

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

GeneralModule 的 `calc_feshbach_resonance_scan` 在磁场网格上计算 $a_s(B)$，并拟合共振参数。

### 对应源码

- [src/mod_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_field_scattering.f90)
- 二维码：`qr/src__mod_field_scattering.f90.png`

![](qr/src__mod_field_scattering.f90.png){width=2.0cm}

## 4.9　含时波包散射

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

定态散射给出能量分辨的 $S(E)$；含时波包通过傅里叶变换也可得到能量分辨量。两者应满足：

- 透射峰位置一致；
- 共振宽度一致；
- 相移或延迟一致；
- 总通量守恒。

GeneralModule 的 `ex07` 示例专门用于定态和含时散射波函数交叉验证。

### 对应源码

- [src/mod_td_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_td_scattering.f90)
- 二维码：`qr/src__mod_td_scattering.f90.png`

![](qr/src__mod_td_scattering.f90.png){width=2.0cm}

## 4.10　分子转振哈密顿量与光谱

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

### 对应源码

- [src/mod_rovibrational.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rovibrational.f90)
- 二维码：`qr/src__mod_rovibrational.f90.png`

![](qr/src__mod_rovibrational.f90.png){width=2.0cm}

## 4.11　光碎片动能释放谱

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

### 对应源码

- [src/mod_photofragment_flux.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_photofragment_flux.f90)
- 二维码：`qr/src__mod_photofragment_flux.f90.png`

![](qr/src__mod_photofragment_flux.f90.png){width=2.0cm}

## 4.12　非绝热动力学与 FSSH

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

若动量不足，则为 frustrated hop。Ehrenfest 方法使用平均场力：

$$
\mathbf F
=
-\left\langle\Psi\left|\frac{\partial\hat H}{\partial\mathbf R}\right|\Psi\right\rangle.
$$

核仍沿经典轨迹运动，电子态保持相干演化。与 FSSH 相比，Ehrenfest 不能自然描述分支化过程。

- FSSH 源码：

[src/mod_surface_hopping_fssh.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_hopping_fssh.f90)

- 二维码：

![](qr/src__mod_surface_hopping_fssh.f90.png){width=2.0cm}

## 4.13　高级专题的统一研究框架

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

## 4.14　参考文献

1. B. R. Johnson, “The multichannel log-derivative method for treating reactive collisions”, *J. Comput. Phys.* **13**, 445 (1973). DOI: 10.1016/0021-9991(73)90049-1.
2. D. E. Manolopoulos, “An improved log-derivative method for solving the radial Schrödinger equation”, *J. Chem. Phys.* **85**, 6425 (1986). DOI: 10.1063/1.451472.
3. H. Feshbach, “Unified theory of nuclear reactions”, *Ann. Phys. (N.Y.)* **5**, 357 (1958). DOI: 10.1016/0003-4916(58)90007-1.
4. C. Chin, R. Grimm, P. S. Julienne, and E. Tiesinga, “Feshbach resonances in ultracold gases”, *Rev. Mod. Phys.* **82**, 1225 (2010). DOI: 10.1103/RevModPhys.82.1225.
5. J. C. Tully, “Molecular dynamics with electronic transitions”, *J. Chem. Phys.* **93**, 1061 (1990). DOI: 10.1063/1.459170.
6. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>

# 第五部分　工程实现、验证与复现

## 5.1　代码架构的分层原则

GeneralModule 的工程结构不是按物理主题随机堆叠，而是按依赖方向分层：

1. 物理常数与单位：`mod_constants`；
2. 数学基础：`mod_special_functions`、`mod_linear_algebra`；
3. 表示与离散：`mod_dvr_grid`、`mod_interpolation`、`mod_io_utils`；
4. 通用动力学：波包传播、开放系统、最优控制；
5. 物理专题：散射、表面、相对论、RIXS、Penning 等；
6. 聚合门面：`general_module`。

低层模块不应反向依赖高层模块，也不应依赖聚合模块 `general_module`。这保证单元测试可以独立编译，也避免循环依赖导致编译顺序不可维护。

模块接口通常遵循以下规范：

- 所有公共实体显式 `public`；
- 数值类型统一使用 `dp => real64`；
- 输入输出通过 `intent(in)`、`intent(out)`、`intent(inout)` 说明；
- 大多数无副作用计算标记为 `pure`；
- 动态数组使用 `allocatable`；
- 错误通过状态码 `stat` 返回，而不是直接终止进程。

这种接口风格便于编译器检查，也便于 AI 辅助生成代码时减少隐式错误。

### 对应源码

- [src/mod_constants.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90)
- 二维码：

![](qr/src__mod_constants.f90.png){width=2.0cm}

- [src/general_module.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/general_module.f90)

## 5.2　三套构建体系的分工

GeneralModule 同时提供 `Makefile`、`fpm.toml` 和 `CMakeLists.txt`。三者不是重复实现，而是服务于不同工作流。

| 构建入口 | 主要用途 | 优点 | 局限 |
|---|---|---|---|
| Makefile | 快速本地构建 | 明确拓扑顺序 | 平台相关性较强 |
| fpm.toml | Fortran 包管理 | 标准化、依赖清晰 | 自动发现关闭 |
| CMakeLists.txt | 跨平台与 IDE | 对象库、静态库、共享库、CTest | 配置复杂 |

`Makefile` 使用显式模块顺序编译，保证 `.mod` 文件依赖正确。`fpm.toml` 显式列出 41 个测试和 36 个示例，避免隐式自动发现。`CMakeLists.txt` 使用对象库组织源码，并生成静态库、共享库和可执行文件。

- Makefile：

[Makefile](https://github.com/l1Ha/QuantumGeneralModule/blob/main/Makefile)

- 二维码：

![](qr/Makefile.png){width=2.0cm}

- fpm 配置：

[fpm.toml](https://github.com/l1Ha/QuantumGeneralModule/blob/main/fpm.toml)

- 二维码：

![](qr/fpm.toml.png){width=2.0cm}

- CMake：

[CMakeLists.txt](https://github.com/l1Ha/QuantumGeneralModule/blob/main/CMakeLists.txt)

- 二维码：

![](qr/CMakeLists.txt.png){width=2.0cm}

## 5.3　测试体系的验证层次

GeneralModule 当前包含 41 个 Fortran 测试程序和 36 个物理示例。测试框架由项目自实现，每个测试程序维护计数器和断言逻辑。

测试应覆盖以下层次：

| 层次 | 检查内容 | 典型失败模式 |
|---|---|---|
| 单位与常数 | 往返转换误差 | 错误量纲、数量级错 |
| 特殊函数 | 选择定则、归一化、对称性 | 相位约定错误 |
| 本征问题 | 正交性、解析谱、条件数 | 非厄米离散 |
| DVR | 网格收敛、边界效应 | 网格欠采样 |
| 传播 | 范数、能量、动量守恒 | 时间步过大 |
| 散射 | 光学定理、幺正性、解析极限 | 匹配半径不当 |
| 开放系统 | 迹、厄米性、正定性 | 强制归一化掩盖误差 |
| 非绝热 | 分支比、能量守恒、统计收敛 | 轨迹数不足 |
| 工程接口 | 编译、链接、可执行性 | 模块顺序错误 |

### 单元测试示例

常数和单位转换测试应包含

$$
\left|
\operatorname{from}_{au}(\operatorname{to}_{au}(x,u),u)-x
\right|
\le
\epsilon|x|.
$$

DVR 测试应检查

$$
\sum_i w_i\psi_i^{(k)}\psi_i^{(l)}
=
\delta_{kl}.
$$

散射测试应检查

$$
\sigma_{\mathrm{tot}}
=
\frac{4\pi}{k}
\operatorname{Im}f(0).
$$

开放系统测试应检查

$$
\left|\operatorname{Tr}\rho(t)-1\right|<\epsilon,
\qquad
\|\rho(t)-\rho^{\dagger}(t)\|<\epsilon.
$$

### 对应源码

- [tests/run_all_tests.sh](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/run_all_tests.sh)
- 二维码：

![](qr/tests__run_all_tests.sh.png){width=2.0cm}

- [tests/test_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_dvr_grid.f90)

- 二维码：

![](qr/tests__test_dvr_grid.f90.png){width=2.0cm}

- [tests/test_ti_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_ti_scattering.f90)

- 二维码：

![](qr/tests__test_ti_scattering.f90.png){width=2.0cm}

## 5.4　持续集成与平台验证

GitHub Actions 配置位于：

[.github/workflows/ci.yml](https://github.com/l1Ha/QuantumGeneralModule/blob/main/.github/workflows/ci.yml)

其作用是自动化编译、测试和示例运行，防止平台相关错误进入主干。CI 应至少检查：

1. 源码能否编译；
2. 静态库和共享库能否生成；
3. 单元测试是否通过；
4. 物理示例是否运行；
5. Python 测试是否通过；
6. 输出文件是否存在且非空。

CI 不能替代物理收敛研究，但可以快速暴露接口和平台相关问题。

## 5.5　Python 交叉验证与可视化

Python 包 `pygenmod` 使用 NumPy 和 Matplotlib 对部分算法进行独立实现、交叉验证和可视化。它不是 Fortran 库的替代品，而是验证与绘图工具。

- NumPy 源码：

[python/pygenmod/visualizer.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py)

- 二维码：

![](qr/python__pygenmod__visualizer.py.png){width=2.0cm}

Python 测试位于：

[python/test_pygenmod.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/test_pygenmod.py)

- 二维码：

![](qr/python__test_pygenmod.py.png){width=2.0cm}

推荐流程是：

1. Fortran 输出原始数据；
2. Python 独立计算参考值；
3. 绘制对比图；
4. 保存输入参数、随机种子、版本和命令；
5. 将图和数据放入同一实验目录。

## 5.6　可复现实验记录

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

没有上述记录，即使结果看起来正确，也难以复现或定位误差。

## 5.7　引用与开源规范

本书引用的文献用于说明理论来源和算法背景。GeneralModule 的实现是独立编写的；引用不表示复制第三方源码。若未来引入第三方代码，必须：

1. 保留其许可证；
2. 在源码头注明来源；
3. 在 `CHANGELOG` 中说明修改范围；
4. 在文档中给出原始仓库或论文链接。

# 第六部分　高级专题与代码映射

## 6.1　三体复合与 Efimov 物理

在共振两体相互作用极限下（$|a|\gg r_0$），全同玻色子三体超径向波函数方程导出特征超越方程

$$
s\cosh\left(\frac{\pi s}{2}\right)
-
\frac{8}{\sqrt{3}}\sinh\left(\frac{\pi s}{6}\right)
=
0.
$$

该方程在实轴上具有唯一正根

$$
s_0\approx1.0062378.
$$

它决定了有效吸引超径向势

$$
V_{\mathrm{eff}}(R)
=
-\frac{s_0^2+1/4}{2\mu R^2}.
$$

三体体系破缺连续标度对称性，展现离散标度不变性：

$$
\lambda=e^{\pi/s_0}\approx22.69438.
$$

Braaten 和 Hammer 对少体普适性作了系统综述【Braaten & Hammer 2006】。

### 对应源码

- [src/mod_three_body_recombination.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_three_body_recombination.f90)

## 6.2　低维受限散射与 CIR

在横向谐振子紧束缚光波导中，

$$
V_\perp(\rho)=\frac12\mu\omega_\perp^2\rho^2,
$$

振子特征长度为

$$
a_\perp=\sqrt{\frac{\hbar}{\mu\omega_\perp}}.
$$

准一维有效接触势耦合常数由横向模态多体格林函数重整化：

$$
g_{\mathrm{1D}}(a_s)
=
\frac{2\hbar^2a_s}{\mu a_\perp^2}
\frac{1}{1-C a_s/a_\perp},
\qquad
C=\frac{|\zeta(1/2)|}{\sqrt2}\approx1.0326.
$$

当三维散射长度逼近

$$
a_{\mathrm{CIR}}=\frac{a_\perp}{C}
\approx0.96843a_\perp
$$

时，$g_{\mathrm{1D}}\rightarrow\pm\infty$，发生约束诱导共振。

Olshanii 给出了准一维约束诱导共振理论【Olshanii 1998】。

### 对应源码

- [src/mod_confined_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_confined_scattering.f90)

## 6.3　Fano 共振与复坐标旋转

当离散准束缚态 $|\phi\rangle$ 与连续态 $|\psi_E\rangle$ 耦合时，吸收截面表现为非对称 Fano 线型：

$$
\sigma(\epsilon)
=
\sigma_0\frac{(q+\epsilon)^2}{1+\epsilon^2},
\qquad
\epsilon=\frac{E-E_0}{\Gamma/2}.
$$

自电离寿命为

$$
\tau=\frac{\hbar}{\Gamma}.
$$

复坐标旋转法通过非厄米变换

$$
U(\theta)=\exp(i\theta r\partial_r)
$$

将坐标旋转为

$$
r\rightarrow re^{i\theta}.
$$

连续谱旋转后，共振极点

$$
E_{\mathrm{res}}=E_R-i\frac{\Gamma}{2}
$$

可从第二黎曼叶暴露出来。

Fano 给出了组态相互作用理论【Fano 1961】；Reinhardt 和 Moiseyev 讨论了复坐标旋转方法【Reinhardt 1982；Moiseyev 1998】。

### 对应源码

- [src/mod_autoionization_fano.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_autoionization_fano.f90)

## 6.4　交叉电磁场与分子取向

极性开壳层分子在非共线静电场和磁场中的哈密顿量可写为

$$
\hat H
=
B_{\mathrm{rot}}\hat{\mathbf J}^2
-
d\mathcal E\cos\theta_R
+
g_S\mu_B
\left(
B_z\hat S_z+B_x\hat S_x
\right).
$$

横向磁场分量 $B_x=B\sin\theta_{EB}$ 产生非对角算符

$$
\hat S_x=\frac{\hat S_++\hat S_-}{2},
$$

驱动 $\Delta M_S=\pm1$ 耦合，破坏沿电场轴的柱对称性，形成可由倾角连续调控的 Stark–Zeeman 反交叉谱。

### 对应源码

- [src/mod_crossed_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_crossed_field_scattering.f90)

## 6.5　三原子反应与 Berry 相位

在 $A+BC$ 反应中，Jacobi 坐标 $(r,R,\gamma)$ 与三原子核间距 $(r_{12},r_{23},r_{31})$ 存在解析映射。London–Eyring–Polanyi–Sato 势能面为

$$
V(r_{12},r_{23},r_{31})
=
Q_1+Q_2+Q_3
-
\sqrt{
\frac12
\left[
(J_1-J_2)^2
+
(J_2-J_3)^2
+
(J_3-J_1)^2
\right]
}.
$$

在锥形交叉附近，绝热波函数环绕退化点一周获得 Berry 相位

$$
\Phi_B=\oint\mathbf A\cdot d\mathbf R=\pi.
$$

Longuet-Higgins 等讨论了 Jahn–Teller 与几何相位【Longuet-Higgins et al. 1958】；Berry 给出了量子几何相位的一般表述【Berry 1984】。

### 对应源码

- [src/mod_triatomic_geometry.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_triatomic_geometry.f90)

## 6.6　旋量 BEC 与自旋混合

自旋 $F=1$ 玻色气体的低能碰撞由总自旋 $F_{\mathrm{tot}}=0,2$ 两个通道的散射长度决定：

$$
c_0=\frac{4\pi\hbar^2(a_0+2a_2)}{3m},
$$

$$
c_2=\frac{4\pi\hbar^2(a_2-a_0)}{3m}.
$$

在单模近似下，旋量分量满足非线性方程

$$
i\hbar\frac{d\zeta_{\pm1}}{dt}
=
\left[
q_Z+c_2n(\rho_0+\rho_{\pm1}-\rho_{\mp1})
\right]\zeta_{\pm1}
+
c_2n\zeta_0^2\zeta_{\mp1}^{*},
$$

$$
i\hbar\frac{d\zeta_0}{dt}
=
c_2n(\rho_{+1}+\rho_{-1})\zeta_0
+
2c_2n\zeta_{+1}\zeta_{-1}\zeta_0^{*}.
$$

系统守恒总几率与磁化强度

$$
m_z=|\zeta_{+1}|^2-|\zeta_{-1}|^2.
$$

Ho 和 Ohmi–Machida 给出了旋量 BEC 理论【Ho 1998；Ohmi & Machida 1998】。

### 对应源码

- [src/mod_spinor_bec.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_spinor_bec.f90)

## 6.7　表面散射与电子摩擦

表面散射需要同时处理表面周期势、入射动量和吸附通道。Eley–Rideal 机理描述气相原子或分子直接撞击表面吸附物并发生反应。电子摩擦则将金属表面电子—空穴对激发等效为阻尼力：

$$
\mathbf F_{\mathrm{fric}}=-\eta(\mathbf R)\dot{\mathbf R}.
$$

这类问题通常需要结合势能面、非绝热耦合和经典轨迹。

### 对应源码

- [src/mod_surface_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_scattering.f90)
- [src/mod_surface_reaction_er.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_reaction_er.f90)
- [src/mod_surface_electronic_friction.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_electronic_friction.f90)

## 6.8　相对论原子与 RIXS

相对论径向 Dirac 方程为

$$
\left[
c\boldsymbol{\alpha}\cdot\mathbf p
+
\beta mc^2
+
V(r)
\right]
\psi
=
E\psi.
$$

它自然给出自旋轨道分裂和精细结构。RIXS 使用 Kramers–Heisenberg 二阶散射截面：

$$
I_{fi}(\omega)
\propto
\left|
\sum_n
\frac{
\langle f|\hat D^{\dagger}|n\rangle
\langle n|\hat D|i\rangle
}
{E_i+\hbar\omega-E_n+i\Gamma_n/2}
\right|^2.
$$

这类模块主要用于内壳层激发、共振非弹性散射和精细结构分析。

### 对应源码

- [src/mod_relativistic_atomic.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_relativistic_atomic.f90)
- [src/mod_resonant_xray_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_resonant_xray_scattering.f90)

## 6.9　Penning 与缔合电离

Penning 电离发生在亚稳态原子与中性原子碰撞时：

$$
A^{*}+B\rightarrow A+B^{+}+e^{-}.
$$

缔合电离则为

$$
A^{*}+B\rightarrow AB^{+}+e^{-}.
$$

GeneralModule 的 `mod_penning_associative_ionization` 使用入口中性态 Morse 势、离子态 Morse 势和距离依赖自电离宽度 $\Gamma(R)$，并计算经典存活概率

$$
S(b)
=
\exp\left[
-2\int_{R_0}^{\infty}
\frac{\Gamma(R)}{\hbar v_R(R)}\,dR
\right].
$$

### 对应源码

- [src/mod_penning_associative_ionization.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_penning_associative_ionization.f90)

## 6.10　强场物理与 HHG

强场物理中常用 Keldysh 参数

$$
\gamma=\frac{\omega\sqrt{2I_p}}{E_0},
$$

和有质动力能

$$
U_p=\frac{E_0^2}{4\omega^2}.
$$

HHG 截止能量近似为

$$
E_{\mathrm{cut}}
=
I_p+3.17U_p.
$$

Corkum 提出三步模型【Corkum 1993】；Lewenstein 等给出强场近似下的偶极响应理论【Lewenstein et al. 1994】。

### 对应源码

- [src/mod_coulomb_atomic.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_coulomb_atomic.f90)
- [src/mod_hhg_spectra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_hhg_spectra.f90)
- [src/mod_strong_field_nsdi.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_strong_field_nsdi.f90)

## 6.11　里德堡阻塞与量子多体

里德堡阻塞源于强偶极相互作用引起的能级移位：

$$
V(R)\sim\frac{C_6}{R^6}
\quad\text{或}\quad
V(R)\sim\frac{C_3}{R^3}.
$$

当相互作用能大于激发激光线宽时，一个里德堡激发会阻止邻近原子再被激发。该机制用于量子门、多体动力学和量子多体疤痕。

### 对应源码

- [src/mod_rydberg_blockade.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rydberg_blockade.f90)

## 6.12　参考文献

1. E. Braaten and H.-W. Hammer, “Universality in few-body systems with large scattering length”, *Phys. Rep.* **428**, 259 (2006). DOI: 10.1016/j.physrep.2006.03.001.
2. M. Olshanii, “Atomic scattering in the presence of an external confinement and a gas of hard-core bosons”, *Phys. Rev. Lett.* **81**, 938 (1998). DOI: 10.1103/PhysRevLett.81.938.
3. U. Fano, “Effects of configuration interaction on intensities and phase shifts”, *Phys. Rev.* **124**, 1866 (1961). DOI: 10.1103/PhysRev.124.1866.
4. W. P. Reinhardt, “Complex coordinates in the theory of atomic and molecular structure and dynamics”, *Annu. Rev. Phys. Chem.* **33**, 223 (1982). DOI: 10.1146/annurev.pc.33.100182.001255.
5. N. Moiseyev, “Quantum theory of resonances: calculating energies, widths and cross-sections by complex scaling”, *Phys. Rep.* **302**, 212 (1998). DOI: 10.1016/S0370-1573(98)00002-7.
6. H. C. Longuet-Higgins et al., “Studies of the Jahn–Teller effect. II. The dynamical problem”, *Proc. R. Soc. Lond. A* **244**, 1 (1958). DOI: 10.1098/rspa.1958.0022.
7. M. V. Berry, “Quantal phase factors accompanying adiabatic changes”, *Proc. R. Soc. Lond. A* **392**, 45 (1984). DOI: 10.1098/rspa.1984.0023.
8. T.-L. Ho, “Spinor Bose condensates in optical traps”, *Phys. Rev. Lett.* **81**, 742 (1998). DOI: 10.1103/PhysRevLett.81.742.
9. T. Ohmi and K. Machida, “Bose–Einstein condensation with internal degrees of freedom in alkali atom gases”, *J. Phys. Soc. Jpn.* **67**, 1822 (1998). DOI: 10.1143/JPSJ.67.1822.
10. P. B. Corkum, “Plasma perspective on strong field multiphoton ionization”, *Phys. Rev. Lett.* **71**, 1994 (1993). DOI: 10.1103/PhysRevLett.71.1994.
11. M. Lewenstein, P. Balcou, M. Y. Ivanov, A. L’Huillier, and P. B. Corkum, “Theory of high-harmonic generation by low-frequency laser fields”, *Phys. Rev. A* **49**, 2117 (1994). DOI: 10.1103/PhysRevA.49.2117.
12. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>
# 第七部分　典型数值实验与代码映射

## 7.1　Morse 势能面的 FGH 本征态

### 物理模型

双原子分子 Morse 势为

$$
V(R)=D_e\left[1-e^{-\beta(R-R_e)}\right]^2.
$$

一维径向哈密顿量为

$$
\hat H
=
-\frac{\hbar^2}{2\mu}
\frac{d^2}{dR^2}
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
V(R).
$$

### 数值实验

1. 选择计算域 $[R_{\min},R_{\max}]$，例如 $[0.8,6.0]$ Bohr。
2. 选择网格数 $N$，例如 $128,256,512,1024$。
3. 用 Sinc-DVR 构造动能矩阵。
4. 在网格点上计算 $V(R_i)$。
5. 组装 $H_{ij}=T_{ij}+V(R_i)\delta_{ij}$。
6. 调用实对称本征求解器。
7. 将本征向量除以 $\sqrt{\Delta x}$ 得到连续归一化波函数。
8. 将前若干能级与 Morse 解析式比较。

### 解析能级

Morse 势的解析能级为

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
\qquad
\omega_ex_e=\frac{\omega_e^2}{4D_e}.
$$

### 验证指标

- 能级误差随 $N$ 收敛；
- 波函数在边界处趋于零；
- 离散内积满足
  $\sum_i\Delta x\,|\psi_i|^2=1$；
- 期望值 $\langle R\rangle$ 和 $\langle R^2\rangle$ 稳定。

### 对应源码

- [examples/ex01_fgh_diatomic_bound_states.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex01_fgh_diatomic_bound_states.f90)
- 二维码：

![](qr/examples__ex01_fgh_diatomic_bound_states.f90.png){width=2.0cm}

- [src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90)

- 二维码：

![](qr/src__mod_dvr_grid.f90.png){width=2.0cm}

## 7.2　一维波包的 Split-Operator 传播

### 物理模型

考虑自由粒子或一维势阱中的高斯波包：

$$
\psi(x,0)
=
\left(\frac{1}{2\pi\sigma^2}\right)^{1/4}
\exp\left[
-\frac{(x-x_0)^2}{4\sigma^2}
+ik_0x
\right].
$$

无外场时，

$$
\langle x\rangle(t)=x_0+\frac{\hbar k_0}{m}t.
$$

### 数值实验

1. 在均匀网格上初始化高斯波包。
2. 设置势能数组 $V(x)$。
3. 使用 Split-Operator 传播多个时间步。
4. 每步后计算 $\langle x\rangle$、$\langle p\rangle$、$\langle H\rangle$ 和 $\|\psi\|$。
5. 改变 $\Delta t$ 和 $N$，绘制误差曲线。

### 验证指标

- 自由粒子时 $\langle x\rangle(t)$ 线性；
- 动量 $\langle p\rangle$ 守恒；
- 能量 $\langle H\rangle$ 守恒；
- 范数 $\|\psi\|$ 守恒；
- 边界处波函数无明显反射。

### 对应源码

- [examples/ex03_split_operator_1d.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex03_split_operator_1d.f90)
- 二维码：

![](qr/examples__ex03_split_operator_1d.f90.png){width=2.0cm}

- [src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)

- 二维码：

![](qr/src__mod_wavepacket_propagator.f90.png){width=2.0cm}

## 7.3　定态与含时散射的交叉验证

### 物理模型

一维散射问题的定态形式为

$$
-\frac{\hbar^2}{2m}\psi''(x)+V(x)\psi(x)=E\psi(x).
$$

含时形式为

$$
i\hbar\frac{\partial\psi}{\partial t}
=
-\frac{\hbar^2}{2m}\psi''(x)+V(x)\psi(x).
$$

### 数值实验

1. 在定态方法中扫描能量 $E$，得到透射率 $T(E)$。
2. 在含时方法中用高斯波包入射，通过通量积分得到 $T(E)$。
3. 比较共振峰位置、宽度和幅度。
4. 改变波包宽度 $\sigma$，检查能量分辨率。
5. 改变 CAP 长度和强度，检查吸收误差。

### 验证指标

- 峰位一致；
- 峰宽一致；
- 透射与反射之和接近 1；
- 不同 $\sigma$ 下结果稳定。

### 对应源码

- [examples/ex07_scattering_wavefunctions_ti_td.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex07_scattering_wavefunctions_ti_td.f90)
- 二维码：

![](qr/examples__ex07_scattering_wavefunctions_ti_td.f90.png){width=2.0cm}

- [src/mod_td_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_td_scattering.f90)

- 二维码：

![](qr/src__mod_td_scattering.f90.png){width=2.0cm}

## 7.4　多通道 Feshbach 扫描

### 物理模型

多通道散射矩阵由开通道 K 矩阵构造：

$$
\mathbf S=(\mathbf I+i\mathbf K)(\mathbf I-i\mathbf K)^{-1}.
$$

Feshbach 共振的散射长度随磁场变化：

$$
a_s(B)
=
a_{\mathrm{bg}}
\left[
1-\frac{\Delta B}{B-B_0}
\right].
$$

### 数值实验

1. 选择磁场范围 $[B_{\min},B_{\max}]$。
2. 对每个 $B$ 构造耦合势矩阵。
3. 调用多通道 Log-Derivative 求解器。
4. 提取入射道散射长度 $a_s(B)$。
5. 拟合 $a_{\mathrm{bg}}$、$\Delta B$ 和 $B_0$。
6. 改变磁场步长和径向网格，检查拟合参数稳定。

### 验证指标

- $S$ 矩阵幺正；
- 开通道与闭通道划分正确；
- 共振位置对网格不敏感；
- 拟合残差可控。

### 对应源码

- [examples/ex08_ultracold_feshbach_segmented.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex08_ultracold_feshbach_segmented.f90)
- 二维码：

![](qr/examples__ex08_ultracold_feshbach_segmented.f90.png){width=2.0cm}

- [src/mod_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_field_scattering.f90)

- 二维码：

![](qr/src__mod_field_scattering.f90.png){width=2.0cm}

## 7.5　Lindblad 弛豫与退相位

### 物理模型

二能级 Lindblad 主方程为

$$
\dot\rho
=
-\frac{i}{\hbar}[H,\rho]
+
\gamma\left(
\sigma_-\rho\sigma_+
-
\frac12\{\sigma_+\sigma_-,\rho\}
\right)
+
\gamma_\phi
\left(
\sigma_z\rho\sigma_z
-
\rho
\right).
$$

### 数值实验

1. 初始化激发态 $\rho_{ee}(0)=1$。
2. 设置弛豫率 $\gamma$ 和退相位率 $\gamma_\phi$。
3. 用 RK4 传播若干时间步。
4. 记录 $\rho_{ee}(t)$、$\rho_{eg}(t)$、$\operatorname{Tr}\rho$、纯度 $P$。
5. 改变 $\Delta t$，检查收敛。

### 验证指标

- $\rho_{ee}(t)$ 指数衰减；
- $|\rho_{eg}(t)|$ 按 $e^{-(\gamma/2+\gamma_\phi)t}$ 衰减；
- 迹守恒；
- 纯度单调下降。

### 对应源码

- [tests/test_open_quantum_opt.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/test_open_quantum_opt.f90)
- 二维码：

![](qr/tests__test_open_quantum_opt.f90.png){width=2.0cm}

- [src/mod_open_quantum.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_open_quantum.f90)

- 二维码：

![](qr/src__mod_open_quantum.f90.png){width=2.0cm}

## 7.6　Tully FSSH 分支比

### 物理模型

FSSH 在活性态 $j$ 上向态 $k$ 的跳跃概率为

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

### 数值实验

1. 初始化核坐标和动量。
2. 选择初始电子态。
3. 对每条轨迹使用 Velocity-Verlet 推进核。
4. 使用 FSSH 规则判断跳跃。
5. 统计透射和反射分支比。
6. 改变轨迹数 $N_{\mathrm{traj}}$，估计统计误差。

### 验证指标

- 分支比随 $N_{\mathrm{traj}}$ 收敛；
- 能量守恒；
- frustrated hop 比例合理；
- 不同随机种子结果稳定。

### 对应源码

- [examples/ex30_tully_surface_hopping.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex30_tully_surface_hopping.f90)
- 二维码：

![](qr/examples__ex30_tully_surface_hopping.f90.png){width=2.0cm}

- [src/mod_surface_hopping_fssh.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_hopping_fssh.f90)

- 二维码：

![](qr/src__mod_surface_hopping_fssh.f90.png){width=2.0cm}

## 7.7　Python 可视化与交叉验证

Fortran 输出通常为文本或二进制数据。Python 可用于：

1. 读取数据；
2. 独立复算关键量；
3. 绘制能级、波函数、截面和误差曲线；
4. 生成发表级图像。

推荐使用 NumPy 处理数组，Matplotlib 绘图【Harris et al. 2020；Hunter 2007】。

### 对应源码

- [python/pygenmod/visualizer.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py)
- 二维码：

![](qr/python__pygenmod__visualizer.py.png){width=2.0cm}

## 7.8　实验记录模板

| 字段 | 建议内容 |
|---|---|
| 实验名称 | 物理问题和目标 |
| 源码版本 | commit、tag 或文件时间 |
| 编译器 | gfortran/ifx/flang 版本 |
| 系统 | OS、CPU、内存 |
| 参数 | 单位、网格、初态、时间步 |
| 随机性 | 种子或确定性声明 |
| 理论基准 | 解析解、文献值、守恒量 |
| 输出 | 能量、截面、范数、分支比 |
| 收敛判据 | 相对误差和阈值 |
| 复现命令 | 完整脚本 |

# 第八部分　源码阅读、模块扩展与工程实践

## 8.1　阅读 GeneralModule 源码的顺序

阅读大型科学计算库时，建议按以下顺序：

1. 先看 `README.md`，确认项目定位、模块清单和示例结构。
2. 再看 `LITERATURE.md`，了解每个模块对应的理论来源。
3. 进入 `src/mod_constants.f90`，掌握类型别名和单位转换。
4. 阅读 `src/mod_linear_algebra.f90`，理解本征、求逆和 FFT。
5. 阅读 `src/mod_dvr_grid.f90`，掌握 DVR 与 FGH。
6. 阅读 `src/mod_wavepacket_propagator.f90`，理解 Split-Operator、RK4、ABM4 和 Bloch 方程。
7. 进入散射和专题模块，结合测试与示例阅读。
8. 最后阅读 `general_module.f90`，了解公共接口如何聚合。

不要从物理专题模块直接入手。若缺少底层数学和表示方法的基础，很容易把实现细节误解为物理结构。

- README：

[README.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/README.md)

- 二维码：

![](qr/README.md.png){width=2.0cm}

- LITERATURE：

[LITERATURE.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/LITERATURE.md)

- 二维码：

![](qr/LITERATURE.md.png){width=2.0cm}

## 8.2　模块模板与接口设计

一个新模块应遵循以下结构：

```fortran
module mod_example
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: example_config_t
    public :: example_result_t
    public :: solve_example

    type :: example_config_t
        real(dp) :: mass = 1.0_dp
        real(dp) :: energy = 0.0_dp
        integer  :: n_grid = 256
    end type example_config_t

    type :: example_result_t
        real(dp), allocatable :: energies(:)
        real(dp), allocatable :: wavefunctions(:, :)
        integer :: stat = 0
    end type example_result_t

contains

    subroutine solve_example(cfg, res)
        type(example_config_t), intent(in)  :: cfg
        type(example_result_t), intent(out) :: res

        res%stat = 0
        if (cfg%n_grid < 8) then
            res%stat = -1
            return
        end if

        allocate(res%energies(cfg%n_grid))
        allocate(res%wavefunctions(cfg%n_grid, cfg%n_grid))
    end subroutine solve_example

end module mod_example
```

关键要求：

- 不要在库内部直接写 `stop`；
- 不要把调试输出写死在核心算法中；
- 不要让模块依赖 `general_module`；
- 不要用全局变量传递状态；
- 不要在接口中隐式转换单位；
- 所有可能失败的分配或求解都应返回状态码。

## 8.3　数值算法的验证模板

每个新算法至少应准备三类测试：

1. 解析极限；
2. 守恒量；
3. 收敛阶。

例如对一个新 DVR，应验证：

$$
\sum_i w_i\psi_i^{(k)}\psi_i^{(l)}=\delta_{kl},
$$

$$
E_N-E
=
CN^{-p}+O(N^{-(p+1)}),
$$

以及对称矩阵的厄米性

$$
H=H^T.
$$

对传播器，应验证：

$$
\|\psi(t+\Delta t)\|_2=\|\psi(t)\|_2,
$$

$$
\langle H\rangle(t)=\langle H\rangle(0).
$$

对散射，应验证：

$$
\sigma_{\mathrm{tot}}=\frac{4\pi}{k}\operatorname{Im}f(0).
$$

## 8.4　测试模板

一个 Fortran 测试程序应包含：

- 测试计数器；
- 断言函数；
- 明确的失败输出；
- 返回码；
- 可重复运行的参数。

```fortran
program test_example
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer :: n_pass, n_total
    real(dp) :: x, y

    n_pass = 0
    n_total = 0

    x = 1.0_dp
    y = x*x
    n_total = n_total + 1
    if (abs(y - 1.0_dp) < 1.0e-12_dp) then
        n_pass = n_pass + 1
    else
        print *, 'FAIL: square'
    end if

    print '(A,I0,A,I0,A)', 'Passed ', n_pass, '/', n_total, ' tests.'
    if (n_pass /= n_total) stop 1

end program test_example
```

测试应避免只检查“程序能运行”。应检查具体数值、守恒量和误差阶数。

## 8.5　CI 中的检查层次

持续集成应至少包含：

| 层次 | 内容 |
|---|---|
| 编译 | 静态库、共享库、可执行文件 |
| 单元测试 | 数学函数、本征、DVR、传播 |
| 物理示例 | 输出存在且非空 |
| Python 测试 | NumPy 交叉验证 |
| 文档 | README、API、公式可读 |
| 复现 | 关键参数和随机种子 |

GitHub Actions 配置位于：

[.github/workflows/ci.yml](https://github.com/l1Ha/QuantumGeneralModule/blob/main/.github/workflows/ci.yml)

## 8.6　如何扩展新物理模块

推荐步骤：

1. 明确物理模型和哈密顿量；
2. 明确坐标、表象和边界条件；
3. 写出核心公式和近似；
4. 搜索文献并确认理论来源；
5. 设计配置类型和结果类型；
6. 实现最小可用算法；
7. 编写解析极限测试；
8. 加入收敛实验；
9. 加入 Python 交叉验证；
10. 更新 README、LITERATURE 和示例索引。

## 8.7　文档与引用规范

每个模块应说明：

- 物理问题；
- 关键假设；
- 输入输出；
- 单位制；
- 算法名称；
- 文献来源；
- 测试位置；
- 示例位置；
- 已知限制。

引用外部文献时，应给出作者、题目、期刊、卷、页码和 DOI。若引用开源软件，应给出项目地址和许可证。若未来引入第三方代码，必须在源码和变更记录中明确说明。

## 8.8　版本与可复现性

推荐使用语义化版本：

- `MAJOR`：接口不兼容变更；
- `MINOR`：新增功能或模块；
- `PATCH`：修复错误或改进文档。

每次重要数值变更都应记录：

1. 修改了哪些算法；
2. 对哪些测试有影响；
3. 是否改变默认精度；
4. 是否改变输出格式；
5. 是否需要重新标定物理参数。

# 第九部分　路线图

## 9.1　算法方向

- 更多正交多项式 DVR；
- 稀疏 Krylov 传播；
- 多维张量分解；
- 高阶分裂算符；
- 自适应时间步长；
- 并行散射计算。

## 9.2　物理方向

- 更多表面反应模型；
- 非绝热耦合场；
- 开放系统谱方法；
- 高精度相对论结构；
- 多通道 Feshbach 系统；
- 少体与多体扩展。

## 9.3　工程方向

- 统一测试断言库；
- 更多 CI 平台；
- 自动生成 API 文档；
- 自动化收敛报告；
- 数据格式标准化；
- 示例输出归档。

## 9.4　文档方向

- 每个模块增加理论说明；
- 每个示例增加收敛数据；
- 将二维码索引扩展到全部源码；
- 增加中英双语 API 说明；
- 增加可复现实验清单。

# 第十部分　源码索引与二维码

以下二维码均指向 GeneralModule 仓库中的对应文件。二维码图像保存在 `docs/qr/`，文件名规则是把源码路径中的 `/` 替换为 `__`。

| 内容 | 路径 | 二维码 |
|---|---|---|
| 仓库主页 | <https://github.com/l1Ha/QuantumGeneralModule> | `docs/qr/repo_root.png` |
| README | [README.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/README.md) | `docs/qr/README.md.png` |
| 配置指南 | [CONFIG_GUIDE.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/CONFIG_GUIDE.md) | `docs/qr/CONFIG_GUIDE.md.png` |
| 文献映射 | [LITERATURE.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/LITERATURE.md) | `docs/qr/LITERATURE.md.png` |
| 常数与单位 | [src/mod_constants.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90) | `docs/qr/src__mod_constants.f90.png` |
| 特殊函数 | [src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90) | `docs/qr/src__mod_special_functions.f90.png` |
| 线性代数 | [src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90) | `docs/qr/src__mod_linear_algebra.f90.png` |
| DVR | [src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90) | `docs/qr/src__mod_dvr_grid.f90.png` |
| 波包传播 | [src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90) | `docs/qr/src__mod_wavepacket_propagator.f90.png` |
| 吸收边界 | [src/mod_absorbing_boundary.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_absorbing_boundary.f90) | `docs/qr/src__mod_absorbing_boundary.f90.png` |
| Chebyshev | [src/mod_chebyshev_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_chebyshev_propagator.f90) | `docs/qr/src__mod_chebyshev_propagator.f90.png` |
| Lindblad | [src/mod_open_quantum.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_open_quantum.f90) | `docs/qr/src__mod_open_quantum.f90.png` |
| Krotov | [src/mod_optimal_control.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_optimal_control.f90) | `docs/qr/src__mod_optimal_control.f90.png` |
| 定态散射 | [src/mod_ti_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_ti_scattering.f90) | `docs/qr/src__mod_ti_scattering.f90.png` |
| 含时散射 | [src/mod_td_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_td_scattering.f90) | `docs/qr/src__mod_td_scattering.f90.png` |
| 外场散射 | [src/mod_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_field_scattering.f90) | `docs/qr/src__mod_field_scattering.f90.png` |
| 转振光谱 | [src/mod_rovibrational.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rovibrational.f90) | `docs/qr/src__mod_rovibrational.f90.png` |
| 光碎片 | [src/mod_photofragment_flux.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_photofragment_flux.f90) | `docs/qr/src__mod_photofragment_flux.f90.png` |
| FSSH | [src/mod_surface_hopping_fssh.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_hopping_fssh.f90) | `docs/qr/src__mod_surface_hopping_fssh.f90.png` |
| Jacobi 坐标 | [src/mod_triatomic_geometry.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_triatomic_geometry.f90) | `docs/qr/src__mod_triatomic_geometry.f90.png` |
| 超球反应 | [src/mod_hyperspherical_reactive.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_hyperspherical_reactive.f90) | `docs/qr/src__mod_hyperspherical_reactive.f90.png` |
| 三体复合 | [src/mod_three_body_recombination.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_three_body_recombination.f90) | `docs/qr/src__mod_three_body_recombination.f90.png` |
| CIR | [src/mod_confined_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_confined_scattering.f90) | `docs/qr/src__mod_confined_scattering.f90.png` |
| Fano/CCR | [src/mod_autoionization_fano.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_autoionization_fano.f90) | `docs/qr/src__mod_autoionization_fano.f90.png` |
| 交叉场 | [src/mod_crossed_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_crossed_field_scattering.f90) | `docs/qr/src__mod_crossed_field_scattering.f90.png` |
| 旋量 BEC | [src/mod_spinor_bec.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_spinor_bec.f90) | `docs/qr/src__mod_spinor_bec.f90.png` |
| 表面散射 | [src/mod_surface_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_scattering.f90) | `docs/qr/src__mod_surface_scattering.f90.png` |
| Eley–Rideal | [src/mod_surface_reaction_er.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_reaction_er.f90) | `docs/qr/src__mod_surface_reaction_er.f90.png` |
| 电子摩擦 | [src/mod_surface_electronic_friction.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_electronic_friction.f90) | `docs/qr/src__mod_surface_electronic_friction.f90.png` |
| 相对论原子 | [src/mod_relativistic_atomic.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_relativistic_atomic.f90) | `docs/qr/src__mod_relativistic_atomic.f90.png` |
| RIXS | [src/mod_resonant_xray_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_resonant_xray_scattering.f90) | `docs/qr/src__mod_resonant_xray_scattering.f90.png` |
| Penning | [src/mod_penning_associative_ionization.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_penning_associative_ionization.f90) | `docs/qr/src__mod_penning_associative_ionization.f90.png` |
| 强场 NSDI | [src/mod_strong_field_nsdi.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_strong_field_nsdi.f90) | `docs/qr/src__mod_strong_field_nsdi.f90.png` |
| 里德堡阻塞 | [src/mod_rydberg_blockade.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rydberg_blockade.f90) | `docs/qr/src__mod_rydberg_blockade.f90.png` |
| FGH 示例 | [examples/ex01_fgh_diatomic_bound_states.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex01_fgh_diatomic_bound_states.f90) | `docs/qr/examples__ex01_fgh_diatomic_bound_states.f90.png` |
| Split 示例 | [examples/ex03_split_operator_1d.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex03_split_operator_1d.f90) | `docs/qr/examples__ex03_split_operator_1d.f90.png` |
| TI/TD 示例 | [examples/ex07_scattering_wavefunctions_ti_td.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex07_scattering_wavefunctions_ti_td.f90) | `docs/qr/examples__ex07_scattering_wavefunctions_ti_td.f90.png` |
| Feshbach 示例 | [examples/ex08_ultracold_feshbach_segmented.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex08_ultracold_feshbach_segmented.f90) | `docs/qr/examples__ex08_ultracold_feshbach_segmented.f90.png` |
| FSSH 示例 | [examples/ex30_tully_surface_hopping.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex30_tully_surface_hopping.f90) | `docs/qr/examples__ex30_tully_surface_hopping.f90.png` |
| Python 可视化 | [python/pygenmod/visualizer.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py) | `docs/qr/python__pygenmod__visualizer.py.png` |

# 结语

本书的目标是把 GeneralModule 从“可编译的代码库”推进到“可验证、可复现、可教学、可扩展的科研资料”。后续版本应继续沿三个方向推进：

1. 理论上补齐更多严格推导和文献映射；
2. 算法上扩展高阶方法和多维表示；
3. 工程上完善测试、CI、可复现性和模块模板。

若发现 Bug、公式错误、代码不一致或复现步骤缺失，请通过仓库 Issue 或以下联系方式反馈：

- 维护者：李皓
- 电子邮箱：LIH_ao@outlook.com
- 仓库：<https://github.com/l1Ha/QuantumGeneralModule>

# 第十一部分　高阶方法与实现细节

## 11.1　Krylov 子空间方法

当哈密顿矩阵过大而无法稠密对角化时，可以使用 Krylov 子空间方法。给定初态 $|\psi_0\rangle$，Krylov 子空间为

$$
\mathcal K_m(\hat H,|\psi_0\rangle)
=
\operatorname{span}
\{
|\psi_0\rangle,
\hat H|\psi_0\rangle,
\ldots,
\hat H^{m-1}|\psi_0\rangle
\}.
$$

Lanczos 算法构造正交基

$$
|\psi_0\rangle,|\psi_1\rangle,\ldots,|\psi_{m-1}\rangle,
$$

使得哈密顿量在该子空间内表示为三对角矩阵

$$
T_m
=
\begin{pmatrix}
\alpha_1 & \beta_1 &        &        \\
\beta_1 & \alpha_2 & \beta_2 &        \\
        & \beta_2 & \ddots & \ddots \\
        &        & \ddots & \alpha_m
\end{pmatrix}.
$$

于是

$$
|\psi(t)\rangle
\approx
V_m
\exp(-i T_m t/\hbar)
e_1,
$$

其中 $V_m$ 的列向量是 Krylov 基。

Lanczos 方法适合求极端本征值和时间演化，但对重复本征值和近简并子空间需要块方法或重启策略。

## 11.2　高阶分裂算符

二阶 Strang 分裂为

$$
S_2(\Delta t)
=
e^{-i\hat V\Delta t/(2\hbar)}
e^{-i\hat T\Delta t/\hbar}
e^{-i\hat V\Delta t/(2\hbar)}.
$$

四阶方法可以使用 Suzuki 分数步长：

$$
S_4(\Delta t)
=
S_2(p_1\Delta t)
S_2(p_2\Delta t)
S_2(p_1\Delta t),
$$

其中

$$
p_1=\frac{1}{2-2^{1/3}},
\qquad
p_2=1-2p_1.
$$

更高阶方法虽然每步代价增加，但对长时间积分和严格要求幺正性的问题更有效。

## 11.3　稀疏矩阵与存储策略

在 DVR 中，若只保留非零矩阵元，可显著降低存储需求。设平均每行非零元数为 $z$，则稀疏存储需求为

$$
O(zN),
$$

而稠密存储需求为

$$
O(N^2).
$$

当 $z\ll N$ 时，稀疏矩阵向量乘法复杂度为

$$
O(zN),
$$

远低于稠密矩阵向量乘法的 $O(N^2)$。

适合稀疏存储的结构包括：

- 局域势能的对角矩阵；
- 近邻耦合动能矩阵；
- 分块对角角动量矩阵；
- 单电子耦合算符。

## 11.4　张量分解与高维问题

高维问题的网格数为

$$
N_{\mathrm{tot}}=\prod_{\alpha=1}^{d}N_\alpha.
$$

直接稠密存储代价为

$$
O(N^{2d}).
$$

张量分解将高阶张量表示为低秩因子乘积，例如 Tucker 分解：

$$
A_{i_1\ldots i_d}
=
\sum_{r_1,\ldots,r_d}
G_{r_1\ldots r_d}
U^{(1)}_{i_1r_1}
\cdots
U^{(d)}_{i_dr_d}.
$$

Tensor-train 表示则为

$$
A_{i_1\ldots i_d}
=
G_1(i_1)
G_2(i_2)
\cdots
G_d(i_d).
$$

若张量秩较低，存储和计算代价可从指数级降低到多项式级。

## 11.5　多通道散射的匹配与渐近展开

在匹配半径 $R_m$ 处，径向波函数可写成入射波和出射波叠加：

$$
\mathbf F(R)
=
\mathbf J(R)\mathbf A
+
\mathbf N(R)\mathbf B,
$$

其中 $\mathbf J$ 和 $\mathbf N$ 分别是规则和不规则径向函数矩阵。定义对数导数

$$
\mathbf Y(R_m)
=
\mathbf F'(R_m)\mathbf F^{-1}(R_m).
$$

由 $\mathbf Y$ 可以构造 $\mathbf K$ 矩阵：

$$
\mathbf K
=
\left[
\mathbf J'(R_m)-\mathbf Y(R_m)\mathbf J(R_m)
\right]^{-1}
\left[
\mathbf Y(R_m)\mathbf N(R_m)-\mathbf N'(R_m)
\right].
$$

再由 Cayley 变换得到散射矩阵：

$$
\mathbf S
=
(\mathbf I+i\mathbf K)(\mathbf I-i\mathbf K)^{-1}.
$$

## 11.6　分段网格的误差控制

分段网格将径向域划分为若干扇区：

$$
[r_{\min},r_{\max}]
=
\bigcup_{s=1}^{S}[r_s,r_{s+1}].
$$

每个扇区可使用不同步长 $\Delta r_s$。若局部误差为

$$
\epsilon_s
\approx
C_s\Delta r_s^{p},
$$

则总误差可近似为

$$
\epsilon_{\mathrm{tot}}
\approx
\sum_s\epsilon_s.
$$

为了等误差分布，可令

$$
C_s\Delta r_s^p
=
\mathrm{const}.
$$

这样在强变化区域使用细网格，在平滑区域使用粗网格。

## 11.7　FSSH 的统计误差

FSSH 是随机轨迹方法。若透射概率为 $p$，轨迹数为 $N$，则标准误差为

$$
\sigma_p
=
\sqrt{\frac{p(1-p)}{N}}.
$$

更稳健的置信区间可使用 Wilson 区间估计：

$$
\frac{
\hat p+\frac{z^2}{2N}
\pm
z\sqrt{
\frac{\hat p(1-\hat p)}{N}
+
\frac{z^2}{4N^2}
}
}
{
1+\frac{z^2}{N}
},
$$

其中 $z$ 是正态分位数。对 95% 置信区间，$z\approx1.96$。

增加轨迹数时，统计误差按 $N^{-1/2}$ 下降。

## 11.8　开放系统的谱方法

Lindblad 主方程为

$$
\dot\rho
=
\mathcal L\rho.
$$

若 $\mathcal L$ 不含时，则形式解为

$$
\rho(t)=e^{\mathcal L t}\rho(0).
$$

可以把超算符 $\mathcal L$ 表示为矩阵并做谱分解：

$$
\mathcal L
=
\sum_\lambda
\lambda |R_\lambda\rangle\langle L_\lambda|,
$$

其中 $|R_\lambda\rangle$ 和 $\langle L_\lambda|$ 分别是右本征算符和左本征算符。于是

$$
\rho(t)
=
\sum_\lambda
e^{\lambda t}
\langle L_\lambda|\rho(0)\rangle
|R_\lambda\rangle.
$$

实部为负的本征值对应衰减模式，纯虚本征值对应相干振荡。

## 11.9　量子控制中的约束

最优控制不仅要最大化保真度，还要满足物理约束：

- 峰值场强 $|\epsilon(t)|\le\epsilon_{\max}$；
- 频谱带宽约束；
- 脉冲面积约束；
- 系统是否能实现目标态；
- 避免激发不可控态；
- 避免超过损伤阈值。

常见惩罚泛函为

$$
J_\epsilon
=
\int_0^T
\lambda(t)\epsilon^2(t)\,dt.
$$

若使用硬约束，则可采用投影或罚函数方法。

## 11.10　数值实验建议

| 主题 | 建议实验 |
|---|---|
| Krylov | 比较不同子空间维数下的传播误差 |
| 高阶分裂 | 比较 $S_2$ 与 $S_4$ 的 $\Delta t$ 标度 |
| 稀疏矩阵 | 对比稠密与稀疏存储的内存和时间 |
| 张量分解 | 在二维势能面上测试秩与误差关系 |
| 多通道散射 | 改变匹配半径，检查 $S$ 矩阵稳定性 |
| 分段网格 | 改变局部步长，检查总误差 |
| FSSH | 增加轨迹数，绘制置信区间 |
| Lindblad | 改变耗散率，检查谱分解 |
| 最优控制 | 改变惩罚参数，比较保真度和场强 |

## 11.11　参考文献

1. G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, Baltimore, 2013.
2. Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, Philadelphia, 2003.
3. M. H. Kalos and P. A. Whitlock, *Monte Carlo Methods*, 2nd ed., Wiley-VCH, Weinheim, 2008.
4. D. J. Tannor, *Introduction to Quantum Mechanics: A Time-Dependent Perspective*, University Science Books, Sausalito, 2007.
5. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>