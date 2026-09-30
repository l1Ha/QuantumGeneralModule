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

每一节均按统一体例组织：**定义**、**公式与推导**、**物理含义**、**数值陷阱**、**GeneralModule 实现映射**（源码路径、函数、GitHub 直链与二维码）、**小型数值实验**与**练习**。所有数值实验结果均由本库源码直接编译运行获得（GNU Fortran 11.4，`real64` 双精度，随机数种子取默认值），可复现、可作为读者自建验收测试的基线。
值得强调的是验收测试的分层思想：第一层是代数恒等式（正交性、对称性、选择定则），应达到机器精度；第二层是解析特例（谐振子、Morse 势、刚转子、Laplacian 精确谱），应达到截断误差与舍入误差的预期标度；第三层才是与实验数据的比对（如 $\mathrm{H}_2$ 的 $B_v$、$D_J$），其偏差属于物理模型的系统误差而非数值误差。三层界限分明，才能在结果异常时迅速定位问题属于代码、算法还是模型。行文约定：公式采用 Pandoc 兼容的 LaTeX 记号，行内公式以 $...$ 界定，独立公式以 $$...$$ 界定。
四章之间存在紧密的逻辑依赖，建议按序阅读。第一章的谱定理为第三章一切变分计算提供合法性依据；第二章的角动量代数在第三章以 Wigner D 函数、Legendre 展开与重耦合网络的形式反复出现；第四章的谱与误差理论则是前三章全部数值实验的度量衡。反过来，第四章的验收测试方法——以已知解析谱校验数值实现——在每一章的小型数值实验中都有实例，读者可将全部实验脚本合并为一个回归测试套件纳入自身的持续集成体系，这也正是本库 340 项单元断言的构建思路。

---

## 第一章 希尔伯特空间与算符谱理论

量子力学的全部数学结构建立在内积空间及其完备化之上。计算量子动力学的每一个数值算法——无论是有限基组展开、离散变量表象（DVR）还是波包传播——都是在希尔伯特空间的某个有限维截断子空间内工作。因此，严格理解无限维理论与有限维近似之间的对应与偏差，是判别数值结果可靠性的前提。本节内容由浅入深：先给出内积与完备性的基本定义，继而讨论基组展开与表象变换，然后陈述自伴算符的谱定理，最后落实到幺正演化的数值实现。

### 1.1 内积、范数与希尔伯特空间

**定义。** 设 $\mathcal{H}$ 为复数域 $\mathbb{C}$ 上的线性空间。若映射 $\langle\,\cdot\,|\,\cdot\,\rangle:\mathcal{H}\times\mathcal{H}\to\mathbb{C}$ 满足对第一变元共轭线性、对第二变元线性、共轭对称 $\langle\psi|\phi\rangle=\langle\phi|\psi\rangle^{*}$ 以及正定性 $\langle\psi|\psi\rangle\ge 0$（等号当且仅当 $\psi=0$），则称其为 $\mathcal{H}$ 上的内积。由内积诱导范数 $\|\psi\|=\sqrt{\langle\psi|\psi\rangle}$，并满足 Cauchy–Schwarz 不等式 $|\langle\psi|\phi\rangle|\le\|\psi\|\,\|\phi\|$ 与三角不等式。若 $\mathcal{H}$ 在该范数下完备，即任意 Cauchy 序列 $\{\psi_n\}$（满足 $\|\psi_n-\psi_m\|\to 0$）都收敛于 $\mathcal{H}$ 中元素，则称 $\mathcal{H}$ 为希尔伯特空间。量子力学标准模型进一步要求 $\mathcal{H}$ 可分，即存在可数的稠密子集。

**公式与推导。** 微观体系中最常见的两类实现为：有限维空间 $\mathbb{C}^{n}$，内积 $\langle u|v\rangle=\sum_i u_i^{*}v_i$；平方可积函数空间 $L^{2}(\mathbb{R}^{3})$，内积 $\langle\psi|\phi\rangle=\int_{\mathbb{R}^{3}}\psi^{*}(\mathbf{r})\phi(\mathbf{r})\,d^{3}\mathbf{r}$，其元素满足 $\int|\psi|^{2}\,d^{3}\mathbf{r}<\infty$。Cauchy–Schwarz 不等式的证明取 $\lambda\in\mathbb{C}$ 并利用 $0\le\|\psi+\lambda\phi\|^{2}=\|\psi\|^{2}+\lambda\langle\psi|\phi\rangle+\lambda^{*}\langle\phi|\psi\rangle+|\lambda|^{2}\|\phi\|^{2}$：当 $\phi\neq 0$ 时取 $\lambda=-\langle\phi|\psi\rangle/\|\phi\|^{2}$，判别式条件立即给出 $|\langle\psi|\phi\rangle|^{2}\le\|\psi\|^{2}\|\phi\|^{2}$。对散射问题，严格而言连续谱态（平面波、散射波）不属于 $L^{2}$，须借助装备希尔伯特空间 $\Phi\subset\mathcal{H}\subset\Phi^{\times}$（rigged Hilbert space）处理；数值计算中则以大盒子边界条件将连续谱离散化，使其纳入可分希尔伯特空间框架。
希尔伯特空间的抽象框架还包括两条对计算至关重要的结构性事实。其一，任何可分无穷维希尔伯特空间在选定正交归一基后与平方可和序列空间 $\ell^{2}$ 等距同构：$L^{2}(\mathbb{R})$ 中的函数问题与 $\ell^{2}$ 中的序列问题可以互相翻译，数值离散化正是这一同构的有限截断。其二，多自由度体系的态空间按张量积组织：$\mathcal{H}=\mathcal{H}_{\mathrm{vib}}\otimes\mathcal{H}_{\mathrm{rot}}\otimes\mathcal{H}_{\mathrm{spin}}$，基组按乘积构造，算符按 $\hat{A}\otimes\hat{B}$ 组合。转振哈密顿量在该乘积基上的矩阵呈现特征性的分块结构——振动指标内近对角、转动指标内带状（$\Delta J=0,\pm1,\pm2$ 型耦合）——这正是第三章各模块矩阵组装代码的骨架，也是稀疏存储与分块对角化得以实施的根源。

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
变分上界性质值得完整推导：由谱定理 $\hat{H}=\sum_k E_k|k\rangle\langle k|$，任意归一化 $\psi=\sum_k c_k|k\rangle$ 满足 $\langle\psi|\hat{H}|\psi\rangle=\sum_k|c_k|^{2}E_k\ge E_0\sum_k|c_k|^{2}=E_0$，等号当且仅当 $\psi$ 为基态。截断子空间内的对角化等价于在该子空间内取 Rayleigh 商极小，故 $E_0^{(N)}\ge E_0$ 且随 $N$ 增大单调不增；Hylleraas–Undheim–MacDonald 交错定理进一步保证 $(N+1)$ 维子空间的第 $k$ 条本征值夹于 $N$ 维结果与精确值之间。收敛速率方面，由于能量是波函数误差的二次泛函，谱方法的能级误差以基组误差的平方量级下降，这解释了 DVR/FGH 计算中能级总是先于波函数收敛的普遍经验，也提示验收时应以波函数形状（节点数、渐近行为）而非仅以能级为收敛判据。Gauss 型求积的另一要点是其代数精度：$n$ 点 Gauss–Legendre 求积对次数不超过 $2n-1$ 的多项式严格成立，由此保证 Legendre DVR 中求积权重与正交归一的严格相容，2.1 节的 $\hat{J}^{2}$ 谱检验正建立在这一性质之上。

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
谱定理的价值在泛函演算中完整体现：对任意合理定义的函数 $f$，$f(\hat{A})=\sum_k f(\lambda_k)|k\rangle\langle k|$，于是传播子 $\exp(-i\hat{H}t/\hbar)$、虚时间算符 $\exp(-\beta\hat{H})$ 与能量投影算符全部化为逐本征值的标量运算——1.4 节的切比雪夫传播器正是该演算的多项式逼近。谱隙 $\mathrm{gap}=|\lambda_{k+1}-\lambda_k|$ 同样具有双重身份：它既控制本征矢的扰动敏感度（4.4 节的 Davis–Kahan 界），又控制含时动力学的绝热性——能隙越大，避免非绝热跃迁所需的脉冲缓变条件越宽松。这把线性代数中的条件数概念与强场物理中的 Landau–Zener 行为直接联系起来。

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
切比雪夫方案的关键构造如下：将哈密顿量按谱半径归一 $\hat{H}_{\mathrm{norm}}=(\hat{H}-\bar{E}\,\mathbb{1})/\Delta E\in[-1,1]$（$\bar{E}$ 与 $\Delta E$ 分别为谱中心与半宽度，可由数次 Lanczos 迭代估计），则传播子的切比雪夫展开系数由第一类 Bessel 函数解析给出：

$$\psi(t)=e^{-i\bar{E}\tau/\hbar}\sum_{n=0}^{N}\left(2-\delta_{n0}\right)(-i)^{n}J_{n}(\tau)\,T_{n}(\hat{H}_{\mathrm{norm}})\,\psi(0),\qquad \tau=\frac{\Delta E\,t}{\hbar}.$$

截断误差由 Bessel 函数尾部 $J_{N+1}(\tau)$ 控制，所需阶数约为 $\tau$ 本身的量级，单步即可覆盖任意长度的演化时间，这是分裂算符格式（每步误差 $\mathcal{O}(\Delta t^{2})$）无法比拟的。其代价有二：每阶需要一次哈密顿–矢量乘法；截断使格式不严格幺正，范数偏差为 $\mathcal{O}(\tau^{N+1}/(N+1)!)$，须以阶数控制。对显式含时哈密顿 $\hat{H}(t)$，上述两类格式均需配合中间点采样，其系统误差由 Magnus 展开的首个对易子 $[\hat{H}(t_1),\hat{H}(t_2)]$ 决定；强场问题中该对易子来自不相互对易的 Stark 项与动能项，是步长收敛检验的必查项。

**物理含义。** 波包传播直接给出光解通量、散射矩阵、含时对齐度等可观测量；幺正性对应封闭体系的概率守恒，任何非幺正格式（如未加修正的显式欧拉法）都会在数千步内累积出物理上无意义的增益或衰减。分裂算符格式的辛对称性还保证了相空间体积守恒的量子对应——能级布居在长时间平均下不漂移。
从可观测量的角度，波包传播的输出直接对应实验时间序列：自相关函数的 Fourier 变换给出吸收谱（本库 `mod_photofragment_flux.f90` 的光解谱模块即按此构建），瞬态布居对应泵浦–探测信号，渐近通量对应散射矩阵元。因此传播格式的每一项保真指标——范数守恒、能量漂移、边界吸收强度——都会以确定的方式映射到光谱线形上：范数泄漏表现为基线漂移，相位误差表现为频移，吸收不完全表现为共振展宽。把数值参数与谱学不确定度定量联系起来，是从"算得出"走向"算得可信"的关键一步。

**数值陷阱。** 其一，FFT 隐含周期边界：波包到达网格边缘后从另一侧回卷（wrap-around），须以足够大的盒子或复吸收势（本库 `mod_absorbing_boundary.f90`）消除。其二，时间步长须同时满足分裂误差与Nyquist 条件：$\Delta t$ 过大时高动量分量的相位 $e^{-iT(p)\Delta t/\hbar}$ 采样不足，产生高频混叠。其三，尽管格式严格幺正，浮点舍入仍以每步 $\mathcal{O}(\varepsilon)$ 引入微小的范数漂移，万步量级后可观测；宜在传播中周期性检查 $\langle\psi|\psi\rangle$ 并记录。

**GeneralModule 实现映射。** 源码 `src/mod_wavepacket_propagator.f90`：`propagate_split_operator_1d`（坐标半步势能相位、FFT 至动量空间、动能整步相位、逆变换、势能后半步，内部调用 `mod_linear_algebra.f90` 的 `fft_1d`）、`propagate_split_operator_2d`、`rk4_step`、`abm4_step`；`src/mod_chebyshev_propagator.f90` 提供切比雪夫大步长推进器。文献依据：Feit, Fleck 与 Steiger（J. Comput. Phys. 47, 412 (1982), DOI: 10.1016/0021-9991(82)90091-2）、Kosloff（J. Phys. Chem. 92, 2087 (1988), DOI: 10.1021/j100319a003）。GitHub 直链：[src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90)。二维码：`qr/src__mod_wavepacket_propagator.f90.png`。

![波包传播模块二维码](qr/src__mod_wavepacket_propagator.f90.png)

**小型数值实验。** 一维谐振子（$m=1$，$\omega=0.02\ E_h/\hbar$），网格 $n=1024$、$\Delta x=0.05\ a_0$、$\Delta t=0.5\ a.u.$，以高斯波包为初态连续推进 $20000$ 步（相当于 $10000\ a.u.\approx242\ \mathrm{fs}$）。实测终态范数偏差 $|\ \|\psi\|^{2}-1\ |=8.8\times10^{-12}$，确认了格式的严格幺正性与舍入漂移的 $\mathcal{O}(N_{\mathrm{step}}\varepsilon)$ 标度。

**练习。** (1) 用 BCH 展开推导对称分裂的 $\mathcal{O}(\Delta t^{3})$ 误差项并证明其为厄米算符的虚倍数。(2) 证明 $e^{-i\hat{V}\Delta t/2\hbar}e^{-i\hat{T}\Delta t/\hbar}e^{-i\hat{V}\Delta t/2\hbar}$ 对任意 $\Delta t$ 幺正。(3) 对谐振子精确解比较分裂格式的相位误差随 $\Delta t$ 的标度，并估计达到 $10^{-6}$ 相位精度所需步长。

---

## 第二章 角动量代数与重耦合理论

角动量是量子力学中唯一具有普适代数结构的守恒量族：无论体系是原子、分子还是超冷碰撞通道，其耦合与变换均由同一套 Racah 代数支配。本章从转动群与角动量本征基出发，依次建立 Clebsch–Gordan 系数与 Wigner 3j 符号、重耦合理论的 6j/9j 符号，最后以 Wigner–Eckart 定理统一处理张量算符矩阵元与光谱选择定则。GeneralModule 在 `mod_special_functions.f90` 中实现了完整的半整数推广符号计算器，是全库外场散射、超精细耦合与转振跃迁模块的代数引擎。
学习本章的推荐路径是：先掌握 2.1 节的代数推导——它只用到线性代数；再以 2.2 节的 CG 系数为基本计算单元，理解 3j 符号作为其对称化包装在求和与对称性质上的优势；随后进入 2.3 节的重耦合结构，体会 6j/9j 符号如何把多体角动量的表象变换压缩为单个标量函数；最后在 2.4 节看到全部抽象如何凝结为光谱学的选择定则。贯穿全章的方法论提示是：角动量代数的每一个公式都具有双重身份——代数身份（对易关系与幺正性）与几何身份（旋转下的变换性质），二者互相校验，任何一方的偏离都应视为实现错误的信号。

### 2.1 转动群与角动量本征基

**定义。** 角动量算符 $\hat{\mathbf{J}}=(\hat{J}_x,\hat{J}_y,\hat{J}_z)$ 由对易关系

$$[\hat{J}_i,\hat{J}_j]=i\hbar\,\varepsilon_{ijk}\hat{J}_k$$

代数地定义，与具体表示（轨道、自旋、分子转动）无关。$\hat{J}^{2}$ 与 $\hat{J}_z$ 对易，故可取共同本征基 $|jm\rangle$：$\hat{J}^{2}|jm\rangle=\hbar^{2}j(j+1)|jm\rangle$，$\hat{J}_z|jm\rangle=\hbar m|jm\rangle$。升降算符 $\hat{J}_{\pm}=\hat{J}_x\pm i\hat{J}_y$ 满足 $[\hat{J}_z,\hat{J}_{\pm}]=\pm\hbar\hat{J}_{\pm}$。

**公式与推导。** 由 $\hat{J}^{2}-\hat{J}_z^{2}=\hat{J}_x^{2}+\hat{J}_y^{2}\ge0$ 得 $m^{2}\le j(j+1)$，故 $m$ 有界；升降算符改变 $m$ 而不改变 $j$，故 $m$ 的取值范围为 $-j,\ldots,j$ 的等差序列，从而 $2j\in\{0,1,2,\ldots\}$，$j$ 可取整数或半奇数。升降系数由模方恒等式 $\hat{J}_{\pm}\hat{J}_{\mp}=\hat{J}^{2}-\hat{J}_z^{2}\pm\hbar\hat{J}_z$ 给出：

$$\hat{J}_{\pm}|jm\rangle=\hbar\sqrt{j(j+1)-m(m\pm1)}\;|j,m\pm1\rangle.$$

轨道角动量的坐标表示为 $\hat{J}^{2}Y_{lm}=\hbar^{2}l(l+1)Y_{lm}$，球谐函数

$$Y_{lm}(\theta,\varphi)=(-1)^{m}\sqrt{\frac{2l+1}{4\pi}\,\frac{(l-m)!}{(l+m)!}}\;P_{l}^{m}(\cos\theta)\,e^{im\varphi}$$

采用 Condon–Shortley 相位约定，其中缔合勒让德函数 $P_l^m$ 由三项递推自 $P_m^m$ 逐级生成。$SO(3)$ 仅容许整数表示，而 $SU(2)$ 双覆盖容许半奇数表示，此即自旋的数学来源。
从群论视角看，$\hat{J}^{2}$ 的不变性源于转动群 Casimir 算符的地位：不可约表示由单一标记 $j$ 刻画，$2j+1$ 维表示空间内的任何矢量在整体转动下仅在该子空间内旋转，这解释了 $m$ 简并为何不受球对称势影响。半奇数表示与整数表示的区别还体现在转动 $2\pi$ 的效应上：整数表示严格复原，半奇数表示获得 $(-1)^{2j}$ 的变号，自旋统计与同核双分子的核自旋统计权重（ortho/para 结构）均由此而来。此外，Wigner D 函数是 $SU(2)$ 表示矩阵在 Euler 角参数化下的显式形式，$d^{J}_{MK}(\beta)$ 可表达为缔合勒让德型函数的组合，因此本节的 $P_l^m$ 递推同时是 3.1 节转动矩阵的数值构件，两节共享同一套稳定递推代码。

**物理含义。** 转动不变性保证 $\hat{J}^{2}$ 与 $\hat{J}_z$ 是好量子数，能级具有 $2j+1$ 重 $m$ 简并；外场沿 $z$ 轴时简并被部分解除（Zeeman 效应）， $m$ 成为光谱标识。分子转动能级 $E_J=B\,J(J+1)$ 正是 $\hat{J}^{2}$ 本征值的直接后果。

**数值陷阱。** 其一，半奇数量子数若以浮点数存储，会因 $0.5$ 的二进制精确表示而侥幸可用，但 $j\pm1/2$ 的链式运算在多模块间传递时仍应以两倍整数（$2j$）贯穿全库——本库 `wigner_3j_half` 系列正是强制采用两倍整数接口以杜绝截断错误。其二，缔合勒让德函数的向上递推对 $|m|\le l$ 稳定，但当 $m$ 接近 $l$ 且 $l$ 很大时 $P_m^m\propto(1-x^{2})^{m/2}$ 急剧减小，与归一化因子的大阶乘相乘后易发生下溢，应在对数域完成阶乘组合（本库 `log_factorial`）。其三，相位约定必须全库统一：Condon–Shortley 相位若在球谐函数与 Wigner 符号两处不一致，交叉矩阵元会出现 $(-1)$ 级别的符号错误，且该错误在强度谱中不可见，仅在对称性校验中暴露。

**GeneralModule 实现映射。** 源码 `src/mod_special_functions.f90`：`legendre_poly`（三项递推计算 $P_l(x)$）、`assoc_legendre_poly`（含 Condon–Shortley 相位的 $P_l^m(x)$）；`src/mod_dvr_grid.f90`：`dvr_legendre_init` 以 Newton 迭代求 Gauss–Legendre 节点并构造 $\hat{J}^{2}$ 角向矩阵。GitHub 直链：[src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)。二维码：`qr/src__mod_special_functions.f90.png`。

![特殊函数模块二维码](qr/src__mod_special_functions.f90.png)

**小型数值实验。** 以 $n=40$ 的 Gauss–Legendre 节点构造 $\hat{J}^{2}$ 的谱求和表示（见 2.2 节），对角化后与本征值序列 $j(j+1)$，$j=0,\ldots,39$ 比对：最大偏差 $1.8\times10^{-12}$。该实验同时验证了 Legendre DVR 与角动量代数的相容性，是角向网格的必备验收测试。

**练习。** (1) 完成升降算符系数的推导。(2) 证明 $Y_{ll}(\theta,\varphi)\propto\sin^{l}\theta\,e^{il\varphi}$ 并归一化。(3) 验证对易关系 $[\hat{J}_x,\hat{J}_y]=i\hbar\hat{J}_z$ 在 $j=1/2$ 的 Pauli 矩阵表示下成立。

### 2.2 Clebsch–Gordan 系数与 Wigner 3j 符号

**定义。** 两个角动量 $\hat{\mathbf{J}}_1,\hat{\mathbf{J}}_2$ 耦合为总角动量 $\hat{\mathbf{J}}=\hat{\mathbf{J}}_1+\hat{\mathbf{J}}_2$ 时，非耦合基 $|j_1m_1\rangle|j_2m_2\rangle$ 与耦合基 $|(j_1j_2)jm\rangle$ 之间的幺正变换系数即 Clebsch–Gordan（CG）系数 $\langle j_1m_1,j_2m_2|jm\rangle$。Wigner 3j 符号是其对称化写法：

$$\langle j_1m_1,j_2m_2|j_3m_3\rangle
=(-1)^{j_1-j_2+m_3}\sqrt{2j_3+1}
\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}.$$

**公式与推导。** 非零条件（选择定则）有三：磁量子数守恒 $m_1+m_2=m_3$；三角条件 $|j_1-j_2|\le j_3\le j_1+j_2$；以及 $j_i\ge|m_i|$。Racah 闭式给出

$$\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}
=(-1)^{j_1-j_2-m_3}\,\Delta(j_1j_2j_3)
\sqrt{\mathcal{N}}\sum_{t}(-1)^{t}\,\mathcal{D}_t^{-1},$$

其中三角系数 $\Delta(j_1j_2j_3)=\sqrt{\dfrac{(j_1+j_2-j_3)!\,(j_1-j_2+j_3)!\,(-j_1+j_2+j_3)!}{(j_1+j_2+j_3+1)!}}$，归一化因子 $\mathcal{N}=\prod_{i=1}^{3}(j_i+m_i)!\,(j_i-m_i)!$，求和项分母为六个阶乘之积

$$\mathcal{D}_t=t!\,(j_1+j_2-j_3-t)!\,(j_1-m_1-t)!\,(j_2+m_2-t)!\,(j_3-j_2+m_1+t)!\,(j_3-j_1-m_2+t)!,$$

$t$ 取使所有阶乘自变量非负的整数区间 $[t_{\min},t_{\max}]$。对称性质：交换两列偶次不变、奇次乘以 $(-1)^{j_1+j_2+j_3}$；全体 $m\to-m$ 反射不变。正交归一关系

$$\sum_{m_1m_2}\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3\end{pmatrix}
\begin{pmatrix} j_1 & j_2 & j_3' \\ m_1 & m_2 & m_3'\end{pmatrix}
=\frac{\delta_{j_3j_3'}\,\delta_{m_3m_3'}}{2j_3+1}$$

与完备性关系 $\sum_{j_3m_3}(2j_3+1)(\cdots)(\cdots)=\delta_{m_1m_1'}\delta_{m_2m_2'}$ 直接源于耦合变换的幺正性，是任何 3j 实现的首选验收测试。
CG 系数可通过最高权态构造法从升降算符完整导出：耦合空间的最高权态 $|j_3,j_3\rangle$ 由要求 $\hat{J}_{+}$ 作用为零唯一确定（允许相差一个整体相位），随后反复施加 $\hat{J}_{-}$ 即生成整个多重态的展开系数；Condon–Shortley 相位约定正是通过规定最高权系数为正实数而锁定。这一构造解释了 Racah 公式中交错求和的来源，也说明 3j 符号的对称性质（列置换、$m$ 反射、Regge 对称）并非独立公理，而是幺正性与相位约定的推论。实用中还应熟记两个特例族：$j_3=0$ 时 3j 符号退化为 $(-1)^{j_1-m_1}\delta_{j_1j_2}\delta_{m_1,-m_2}/\sqrt{2j_1+1}$，常用于标量积矩阵元与约化矩阵元的化简；三个磁量子数全为零时 3j 符号仅在 $j_1+j_2+j_3$ 为偶数时非零，此性质广泛用于宇称选择定则。

**物理含义。** CG 系数支配一切两体角动量耦合：自旋–轨道耦合、原子超精细结构 $(\mathbf{s}\,\mathbf{i})\mathbf{f}$、双原子碰撞通道耦合 $(s_1i_1)f_1,(s_2i_2)f_2$，以及转振跃迁的方向余弦矩阵元。单重态与三重态 $(\tfrac12\,\tfrac12)S=0,1$ 的展开系数 $1/\sqrt2$ 即最简单的 CG 系数，是磁 Feshbach 共振多通道计算的入口（本库 LITERATURE 全典第 5.2 节的四大基组变换即建立在 CG 与 9j 符号之上）。
分子光谱中最直观的例子是碱金属二聚体的单重–三重结构：电子自旋交换作用 $\hat{\mathbf{s}}_1\cdot\hat{\mathbf{s}}_2$ 在总自旋基下对角，本征值由 $\tfrac12[S(S+1)-s_1(s_1+1)-s_2(s_2+1)]$ 给出，单重态与三重态的势能曲线劈裂由此完全确定；而从实验可测的塞曼子能级出发重构耦合方案时，CG 系数正是连接两套表象的词典。另一个例子是转振跃迁的 Hönl–London 因子：线强度中转动部分的全部信息包含于 $\langle J'|\cos\theta|J\rangle$ 型矩阵元的模方，其 $J$ 依赖完全由 3j 符号的几何因子决定，与振动态无关。

**数值陷阱。** 其一，阶乘溢出：$j\gtrsim 20$ 时 $(2j)!$ 超出双精度范围，必须以对数阶乘 $\ln n!$ 累加后再指数化——本库以 `log_factorial` 贯穿全部符号计算。其二，Racah 求和为交错级数，$j$ 较大时出现严重相消，相对误差可达 $10^{-8}$ 以上；高 $j$ 场合应改用递推关系或扩充足位精度。其三，相位约定：CG 与 3j 之间的相位因子 $(-1)^{j_1-j_2+m_3}$ 若与库内其他模块不一致，将使干涉项符号翻转。

**GeneralModule 实现映射。** 源码 `src/mod_special_functions.f90`：`wigner_3j` 与 `clebsch_gordan`（整数角动量）、`wigner_3j_half` 与 `clebsch_gordan_half`（两倍整数接口，覆盖半奇数自旋）、`log_factorial`（对数阶乘辅助）。文献依据：D. A. Varshalovich, A. N. Moskalev 与 V. K. Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific, Singapore (1988), DOI: 10.1142/0270。GitHub 直链：[src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)。二维码：`qr/src__mod_special_functions.f90.png`。

**小型数值实验。** 四项独立校验：(a) 值检验：$\begin{pmatrix}1&1&2\\0&0&0\end{pmatrix}=0.3651483717$，与精确值 $\sqrt{2/15}=0.3651483717$ 一致；(b) 对称性检验：$\begin{pmatrix}3&2&4\\1&1&-2\end{pmatrix}=+0.06299408$ 而 $\begin{pmatrix}2&3&4\\1&1&-2\end{pmatrix}=-0.06299408$，比值恰为 $(-1)^{j_1+j_2+j_3}=(-1)^{9}=-1$；(c) 正交性检验：$\sum_{m_1}\bigl[\begin{pmatrix}3&2&4\\m_1&-m_1&0\end{pmatrix}\bigr]^{2}=0.1111111$，与 $1/(2\cdot4+1)=1/9$ 一致；(d) CG 检验：$\langle 1\,1,1\,0|2\,1\rangle=0.70710678=1/\sqrt2$。全部达到机器精度。

**练习。** (1) 由 Racah 公式计算 $\begin{pmatrix}1&1&0\\m&-m&0\end{pmatrix}$ 并验证其值为 $(-1)^{1-m}/\sqrt3$。(2) 从耦合变换的幺正性出发推导正交关系。(3) 编写程序验证 $j_1=3,j_2=2$ 时对全部 $(j_3,m_3)$ 的完备性关系，统计求和项数以估计工作量。

### 2.3 Wigner 6j 与 9j 符号：重耦合理论

**定义。** 三个角动量 $j_1,j_2,j_3$ 耦合到总角动量 $J$ 存在多条路径：先耦合 $(j_1j_2)\to j_{12}$ 再与 $j_3$ 耦合，或先耦合 $(j_2j_3)\to j_{23}$ 再与 $j_1$ 耦合。两套耦合基之间的幺正变换系数由 Wigner 6j 符号（Racah W 系数）给出：

$$\bigl\langle (j_1j_2)j_{12},j_3;JM \big| j_1,(j_2j_3)j_{23};JM \bigr\rangle
=(-1)^{j_1+j_2+j_3+J}\sqrt{(2j_{12}+1)(2j_{23}+1)}
\begin{Bmatrix} j_1 & j_2 & j_{12} \\ j_3 & J & j_{23} \end{Bmatrix}.$$

四个角动量的重耦合（如原子结构中的 $LS$ 与 $jj$ 耦合方案互换）则由 Wigner 9j 符号承担：

$$\bigl\langle (j_{11}j_{12})j_{13},(j_{21}j_{22})j_{23};J \big| (j_{11}j_{21})j_{31},(j_{12}j_{22})j_{32};J \bigr\rangle
=\sqrt{(2j_{13}+1)(2j_{23}+1)(2j_{31}+1)(2j_{32}+1)}
\begin{Bmatrix} j_{11} & j_{12} & j_{13} \\ j_{21} & j_{22} & j_{23} \\ j_{31} & j_{32} & J \end{Bmatrix}.$$

**公式与推导。** 6j 符号满足四重三角条件（四面体四面的三角关系）与正交关系

$$\sum_{x}(2x+1)\begin{Bmatrix} a & b & x \\ c & d & p \end{Bmatrix}\begin{Bmatrix} a & b & x \\ c & d & q \end{Bmatrix}
=\frac{\delta_{pq}}{2p+1},$$

该关系直接源于重耦合变换的幺正性。9j 符号可展开为 6j 符号的三重求和（Varshalovich, 第 10 章）：

$$\begin{Bmatrix} j_{11} & j_{12} & j_{13} \\ j_{21} & j_{22} & j_{23} \\ j_{31} & j_{32} & j_{33} \end{Bmatrix}
=\sum_{x}(-1)^{2x}(2x+1)
\begin{Bmatrix} j_{11} & j_{21} & j_{31} \\ j_{32} & j_{33} & x \end{Bmatrix}
\begin{Bmatrix} j_{12} & j_{22} & j_{32} \\ j_{21} & x & j_{23} \end{Bmatrix}
\begin{Bmatrix} j_{13} & j_{23} & j_{33} \\ x & j_{11} & j_{12} \end{Bmatrix},$$

求和变量 $x$ 取遍同时满足三个三角条件的（整数或半奇数）值。本库实现严格遵循该式，其中相位因子在两倍整数接口下表现为 $(2x)$ 奇偶性的判别。
重耦合理论还有两条常用的进阶关系。其一，Biedenharn–Elliott 求和恒等式约束四个 6j 符号乘积的求和，是推导高阶递推关系与渐近分析的工具。其二，Ponzano–Regge 半经典渐近式把 6j 符号与欧氏四面体的几何联系起来：当六个角动量远大于 1 时，6j 符号由四面体体积决定的振荡相位乘以 $V^{-1/2}$ 型振幅给出；这为密耦计算中的通道截断提供了物理直觉——大角动量通道的重耦合系数快速振荡，对总截面的贡献相消。9j 符号另具 72 个元素的对称群（行列置换与转置），高效实现可利用该对称群将重复求值减少约一个量级。

**物理含义。** 重耦合理论是现代少体与碰撞物理的枢纽。在超冷原子碰撞中，双原子通道基 $|(s_1i_1)f_1,(s_2i_2)f_2\rangle$ 与总自旋基 $|(s_1s_2)S,(i_1i_2)I\rangle$ 之间的幺正变换由 9j 符号给出：

$$\langle (f_1f_2)F | (SI)F \rangle
=\sqrt{(2f_1+1)(2f_2+1)(2S+1)(2I+1)}
\begin{Bmatrix} s_1 & i_1 & f_1 \\ s_2 & i_2 & f_2 \\ S & I & F \end{Bmatrix},$$

本库 `mod_field_scattering.f90` 的 `calc_basis_transform_matrix` 以机器精度实现该变换。在超球坐标反应散射中，广义角动量本征函数（超球谐函数）的耦合与排列对称化同样由 3j/6j 网络承担； Efimov 三体物理中的超径向方程分离变量亦依赖 6j 重耦合。

**数值陷阱。** 其一，三角条件判断必须与求和上下界严格一致：本库 `triangle_half` 在三角条件不满足时返回零，`wigner_6j_half` 依此决定 $t$ 求和区间；若调用者以浮点传入量子数，区间端点的整数运算会静默丢失半奇数信息。其二，6j 求和同样为交错级数，$j\gtrsim 30$ 时相消显著。其三，9j 的计算量为 $\mathcal{O}(j)$ 个 6j、每个 6j 为 $\mathcal{O}(j)$ 项求和，总计 $\mathcal{O}(j^{2})$；在大通道数密耦计算中应缓存重复调用的符号值。

**GeneralModule 实现映射。** 源码 `src/mod_special_functions.f90`：`triangle_half`（辅助三角系数）、`wigner_6j_half`（Racah 单求和公式）、`wigner_9j_half`（三重 6j 展开）；配套 `src/mod_field_scattering.f90` 的 `calc_basis_transform_matrix` 与 `build_field_collision_channels`。GitHub 直链：[src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)。二维码：`qr/src__mod_special_functions.f90.png`。

**小型数值实验。** 三项校验：(a) 正交性：$\sum_x(2x+1)\begin{Bmatrix}1&1&x\\2&2&2\end{Bmatrix}\begin{Bmatrix}1&1&x\\2&2&4\end{Bmatrix}$ 实测为 $0$（$<10^{-30}$），符合 $\delta_{24}/5=0$；(b) 值检验：$\begin{Bmatrix}\tfrac12&\tfrac12&1\\ \tfrac12&\tfrac12&1\end{Bmatrix}=0.1666666667=1/6$；(c) 独立重构检验：9j 符号 $\begin{Bmatrix}\tfrac12&\tfrac12&1\\ \tfrac12&\tfrac12&0\\ 1&0&1\end{Bmatrix}$ 实测 $0.1666666667$；与此同时，不经 6j 而以 CG 系数四级联直接在非耦合基（四个 $1/2$ 自旋的 16 维空间）中重构重耦合系数 $\langle(f_1f_2)F|(SI)F\rangle$，得 $0.5000000000$；按定义式 $\sqrt{(2f_1+1)(2f_2+1)(2S+1)(2I+1)}\times\text{9j}=3\times(1/6)=0.5$，两条完全独立的计算路径在机器精度内吻合。

**练习。** (1) 证明 6j 正交关系由重耦合幺正性直接导出。(2) 用 CG 级联数值重构 $\langle(f_1f_2)F|(SI)F\rangle$ 在 $f_1=1,f_2=2,S=1,I=2,F=3$ 的值，并与 9j 公式比对。(3) 证明 9j 符号任一行或列满足三角条件。

### 2.4 不可约张量算符、Wigner–Eckart 定理与选择定则

**定义。** 相对于角动量 $\hat{\mathbf{J}}$ 的 $k$ 阶不可约张量算符是一组 $2k+1$ 个分量 $\{T^{(k)}_q\}$，$q=-k,\ldots,k$，满足

$$[\hat{J}_z,T^{(k)}_q]=\hbar q\,T^{(k)}_q,\qquad
[\hat{J}_{\pm},T^{(k)}_q]=\hbar\sqrt{k(k+1)-q(q\pm1)}\,T^{(k)}_{q\pm1}.$$

Wigner–Eckart 定理断言其矩阵元分解为几何因子与约化矩阵元之积：

$$\langle j'm'|T^{(k)}_q|jm\rangle
=(-1)^{j'-m'}\begin{pmatrix} j' & k & j \\ -m' & q & m \end{pmatrix}\langle j'\|T^{(k)}\|j\rangle.$$

**公式与推导。** 该定理的证明思路是：算符 $T^{(k)}_q|jm\rangle$ 按 $\hat{\mathbf{J}}$ 的表示变换的方式与 $|kq\rangle|jm\rangle$ 的耦合乘积相同，故其在磁量子数子空间内的系数必为 CG 系数；与 $m$ 无关的约化矩阵元吸收全部动力学信息。
高阶张量可由低阶张量经耦合乘积递归构造：

$$\bigl[\mathbf{A}^{(k_1)}\otimes\mathbf{B}^{(k_2)}\bigr]^{(k)}_{q}
=\sum_{q_1q_2}\langle k_1q_1,k_2q_2|kq\rangle\,A^{(k_1)}_{q_1}B^{(k_2)}_{q_2},$$

其约化矩阵元由 6j 符号给出（张量解耦公式）。该构造在多体问题中递归进行：三体相互作用的矩阵元最终化为 6j/9j 符号网络，与 2.3 节的重耦合理论合流。这体现了角动量代数的结构性优势：以少数几个基本符号及其代数关系，覆盖任意多体体系的全部角动量矩阵元，且每一步都可机器验证。直接推论为选择定则：$|j-k|\le j'\le j+k$ 且 $m'=m+q$。线性分子沿场轴的偶极作用 $\hat{H}_{d}=-dE\cos\theta$ 中，$\cos\theta=C^{(1)}_0(\theta,\varphi)$ 是秩 1 张量的 $q=0$ 分量，故平行跃迁满足 $\Delta J=\pm1$、$\Delta m=0$；非共振激光极化作用 $\frac14\Delta\alpha E^{2}\cos^{2}\theta$ 分解为零秩与二秩两部分，$\cos^{2}\theta$ 的矩阵元仅在 $\Delta J=0,\pm2$ 时非零（拉曼与分子对齐选择定则）。闭式矩阵元为

$$\langle J m|\cos\theta|J\pm1,m\rangle
=\sqrt{\frac{(J\pm1)^{2}-m^{2}}{(2J\pm1)(2J\pm3)}},\qquad
\langle J m|\cos^{2}\theta|J m\rangle
=\frac{(J+1)^{2}-m^{2}}{(2J+1)(2J+3)}+\frac{J^{2}-m^{2}}{(2J-1)(2J+1)}.$$

**物理含义。** 选择定则是光谱学的语言：红外平行带 $\Delta J=\pm1$、拉曼散射 $\Delta J=0,\pm2$、多极相互作用对分波的约束（如磁偶极–偶极作用按二秩张量展开驱动 $L\leftrightarrow L\pm2$ 耦合但严格守恒总 $M$），全部由 Wigner–Eckart 定理统一给出。约化矩阵元与几何因子的分离还意味着：一旦计算了与取向无关的动力学量，任意取向依赖的矩阵元都是零成本的代数运算。
给出一个完整的计算链条示例：双原子长程磁偶极–偶极作用的二秩张量展开 $\hat{V}_{dd}\propto-\frac{\sqrt{6}}{r^{3}}\sum_{q}(-1)^{q}C_{2,-q}(\hat{\mathbf{r}})\,[\mathbf{s}_1\otimes\mathbf{s}_2]^{(2)}_{q}$ 中，空间因子 $\langle L'M_L'|C_{2q}|L M_L\rangle$ 化为 3j 符号组合 $\sqrt{(2L+1)(2L'+1)}\begin{pmatrix} L' & 2 & L \\ 0 & 0 & 0\end{pmatrix}\begin{pmatrix} L' & 2 & L \\ -M_L' & q & M_L\end{pmatrix}$，其中 $(L'\,2\,L;0\,0\,0)$ 因子强制 $L+L'+2$ 为偶数，故 $s$ 波（$L=0$）仅与 $d$ 波（$L'=2$）耦合；自旋张量因子的矩阵元则由 6j 与 CG 网络给出。全部矩阵元在机器精度内解析可得，这正是各向异性偶极自旋弛豫截面计算得以高度自动化的原因，也说明 Wigner–Eckart 定理不仅给出选择定则，还直接给出矩阵元的完整算法。

**数值陷阱。** 其一，约化矩阵元约定不一：文献中存在相差 $\sqrt{2j+1}$ 因子或相位的多种定义，跨文献移植公式时必须逐条核对；本库以显式矩阵元函数规避隐式约定。其二，$\cos^{2}\theta$ 矩阵元的两条闭式在 $J=0$ 边界涉及 $(-1)!$ 型奇点，实现中须以分支判断屏蔽（本库 `rot_matrix_cos2_theta` 对 $J=0$ 的对角项单独处理）。其三，$m$ 截断条件 $|m|\le\min(j,j')$ 忘记判断时会产生虚假的非零矩阵元，破坏 $m$ 分块对角结构。

**GeneralModule 实现映射。** 源码 `src/mod_special_functions.f90`：`rot_matrix_cos_theta`（$\langle j m|\cos\theta|j' m\rangle$，强制 $\Delta j=\pm1$）、`rot_matrix_cos2_theta`（$\Delta j=0,\pm2$）；`src/mod_rovibrational.f90`：`build_rovibrational_dipole_matrix`（按 $\Delta J=\pm1$ 组装全转振偶极矩阵）、`build_rovibrational_polarizability_matrix`（按 $\Delta J=0,\pm2$ 组装极化率矩阵）。GitHub 直链：[src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)。二维码：`qr/src__mod_special_functions.f90.png`。

**小型数值实验。** 调用 `rot_matrix_cos_theta(1,2,0)` 得 $0.51639778$，与精确值 $2/\sqrt{15}=0.5163977795$ 吻合；调用 `rot_matrix_cos_theta(1,1,0)` 得 $0$，确认 $\Delta J=0$ 偶极禁阻；调用 `rot_matrix_cos2_theta(2,2,0)` 得 $0.52380952$，与解析式 $9/35+4/15=0.5238095238$ 一致。三例共同验证了张量矩阵元实现的选择定则与归一化。

**练习。** (1) 由 Wigner–Eckart 定理导出 $\langle J m|\cos\theta|J+1,m\rangle$ 的闭式。(2) 证明 $\cos^{2}\theta$ 的矩阵元仅在 $\Delta J=0,\pm2$ 时非零，并给出 $\Delta J=\pm2$ 项的显式表达式。(3) 对电四极跃迁（秩 2 张量）写出选择定则并讨论与拉曼跃迁的异同。

---

## 第三章 分子坐标体系与转振哈密顿量

分子的转振动力学在数学上等价于在合适的广义坐标下分离并处理动能算符。坐标选择的优劣直接决定哈密顿量矩阵的稀疏性、耦合的强弱与基组收敛速度。本章由浅入深地讨论四套坐标体系：以 Euler 角连接的空间固定系与体固定系（3.1 节）、描述反应碰撞的 Jacobi 坐标（3.2 节）、以质量标度实现通道对称的超球坐标（3.3 节），以及统一转振自由度的 Watson 哈密顿量及其微扰展开——Coriolis 耦合、离心畸变与振动角动量（3.4 节）。
坐标选择的原则可以概括为两条：其一，动能算符的结构应尽可能简单，即对角主导、耦合项低阶；其二，哈密顿量的分块应与守恒量或近守恒量对齐，以实现最大程度的解耦。SF/BF 分解牺牲前者换取后者——动能中出现 Coriolis 交叉项，但势能面获得体固定系下的不变表达；Jacobi 坐标与超球坐标则针对反应散射的通道结构优化分块，使三套排列通道在同一组变量下对称出现。理解每种坐标体系以何种代价换取何种可解释性，是构造高效转振与反应散射算法的先决条件。

### 3.1 空间固定系、体固定系、Euler 角与 Wigner D 函数

**定义。** 空间固定系（SF，space-fixed frame）$\{X,Y,Z\}$ 以实验室为参照，外场方向通常取 $Z$ 轴；体固定系（BF，body-fixed frame）$\{x,y,z\}$ 固连于分子，通常取惯性主轴。两系之间的转动以三个 Euler 角 $(\alpha,\beta,\gamma)$（$z$-$y$-$z$ 约定）参数化，对应的完备正交函数族是 Wigner D 函数

$$D^{J}_{MK}(\alpha,\beta,\gamma)=e^{-iM\alpha}\,d^{J}_{MK}(\beta)\,e^{-iK\gamma},$$

满足正交归一 $\frac{2J+1}{8\pi^{2}}\int D^{J'*}_{M'K'}D^{J}_{MK}\,d\Omega=\delta_{JJ'}\delta_{MM'}\delta_{KK'}$，其中 $d\Omega=d\alpha\,\sin\beta\,d\beta\,d\gamma$。$M$ 与 $K$ 分别为总角动量在 SF 轴与 BF 轴上的投影量子数。

**公式与推导。** 刚性转子的经典动能为 $T=\frac12\sum_{\alpha}I_{\alpha}\omega_{\alpha}^{2}$（$\alpha=a,b,c$ 为主轴），量子化后得不对称陀螺哈密顿量

$$\hat{H}_{\mathrm{rot}}=A\hat{J}_a^{2}+B\hat{J}_b^{2}+C\hat{J}_c^{2},\qquad A=\frac{\hbar^{2}}{2I_a},\ B=\frac{\hbar^{2}}{2I_b},\ C=\frac{\hbar^{2}}{2I_c}.$$

在 $|JK\rangle$ 基下利用 $\hat{J}_{\pm}$ 的梯性质，$\hat{J}_a^{2},\hat{J}_c^{2}$ 产生 $\Delta K=\pm2$ 的非对角元，矩阵在固定 $J$ 的 $2J+1$ 维 $K$ 子空间内对角化即得不对称陀螺能级；对称陀螺（$B=C$）情形解析可解：

$$E_{JK}=B\,J(J+1)+(A-B)K^{2}.$$
不对称陀螺在 $|JK\rangle$ 基下的非对角矩阵元为

$$\langle J,K|\hat{J}_a^{2}|J,K\pm2\rangle=\frac{\hbar^{2}}{4}\sqrt{(J\mp K)(J\mp K-1)(J\pm K+1)(J\pm K+2)},$$

于是固定 $J$ 的哈密顿量在 $K$ 空间内是带宽为二的带状矩阵，可按 $K$ 的奇偶分块，再以宇称组合 $|K\rangle\pm|-K\rangle$ 将规模减半。核自旋统计进一步把允许的转动态按置换对称性分族——例如水分子 $K_a+K_c$ 的奇偶对应 ortho 与 para 两族——各族能级的布居比由核自旋简并度冻结，这正是转振光谱拟合中统计权重因子的来源，也是同核分子（如 $\mathrm{H}_2$ 的奇偶 $J$ 族）红外谱缺失的深层原因。

SF 与 BF 角动量分量之间以方向余弦（即 D 函数）相连：$\hat{J}_Z=\sum_K D^{J}_{MK}(\alpha\beta\gamma)\,\hat{J}_{z'}$ 型关系是全部体固定相互作用（势能面、偶极矩、极化率张量）与实验室可观测量之间换算的枢纽。线型分子的取向分布退化为 $K=0$ 的特殊情形，此时 $D^{J}_{M0}(\alpha\beta\gamma)\propto Y_{JM}(\beta,\alpha)$，方向余弦矩阵元 $\langle J M|\cos\theta|J' M\rangle$ 即 2.4 节的闭式。

**物理含义。** 势能面与偶极面是分子内禀属性，天然表达在体固定系中；而光谱跃迁强度、外场对齐与散射边界条件表达在空间固定系中。两套表象的取舍是转振计算的核心策略：以 $|J M K\rangle$（对称陀螺）或宇称组合 $|J M K\rangle\pm|J M,-K\rangle$ 为基，可将哈密顿量分块至 $(J,M)$ 或 $(J,\text{宇称})$ 子空间，使矩阵维数下降一至两个量级。
实验上可调的物理量——外场方向、偏振、波长——定义在空间固定系，而势能面与偶极面是分子内禀属性；两套语言之间转译的精度决定了模拟与实验可比对的深度。以激光对齐为例，非共振场作用 $\tfrac14\Delta\alpha E^{2}\cos^{2}\theta$ 在空间固定系中是沿实验室 $Z$ 轴固定的张量，转换到体固定系后化为对分子轴的标量作用，正是这一转译使得对齐动力学可以按 $(J,M)$ 分块求解；对齐度 $\langle\cos^{2}\theta\rangle$ 的时间演化则由少数几个分块的本征展开叠加而成。同理，Stark 移位与 Zeeman 移位的方向依赖性也全部由 D 函数承载。

**数值陷阱。** 其一，Euler 角在 $\beta=0,\pi$ 存在坐标奇异性，经典轨迹法在极区必须切换参数化，而基函数法则（D 函数展开）自动规避；两套方法衔接时须显式检验。其二，转动约定（主动/被动、$z$-$y$-$z$ 与 $x$-$y$-$z$ 约定、D 函数与旋转矩阵的复共轭关系）在不同教材之间不一致，跨文献移植公式时极易引入 $\beta\to\pi-\beta$ 或复共轭级别的错误；验收手段是对已知解析矩阵元（如 $\langle 10|\cos\theta|20\rangle$）逐一比对。其三，$K$ 截断须与势能面的各向异性匹配：$K_{\max}$ 过小会系统性高估转动激发能，且收敛不一定从下方单调。

**GeneralModule 实现映射。** 源码 `src/mod_special_functions.f90`：`assoc_legendre_poly`（$K=0$ 线型转子的 $d^{J}_{M0}\propto P^{M}_{J}$ 构件）、`rot_matrix_cos_theta` 与 `rot_matrix_cos2_theta`（方向余弦矩阵元）；`src/mod_dvr_grid.f90`：`dvr_legendre_init`（$\theta$ 方向 Gauss–Legendre DVR，可组装 $\hat{J}^{2}$ 与 $\cos\theta$ 矩阵）；`src/mod_rovibrational.f90`：`build_rovibrational_hamiltonian`（$E(v,J)=E_{\mathrm{vib}}(v)+B_vJ(J+1)$）。文献依据：Bunker 与 Jensen, *Molecular Symmetry and Spectroscopy*, 2nd ed., NRC Research Press, Ottawa (1998)；Varshalovich 等 (1988)。GitHub 直链：[src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)。二维码：`qr/src__mod_special_functions.f90.png`。

**小型数值实验。** 以 2.2 节验证过的 $\hat{J}^{2}$ 谱求和矩阵表示线型刚性转子 $\hat{H}=B\hat{J}^{2}$，取 $\mathrm{H}_2$ 的振动平均转动常数 $B_0=60.18\ \mathrm{cm^{-1}}$（见 3.4 节实验）：对角化得 $E_J/B_0=j(j+1)$，$E_1=120.36\ \mathrm{cm^{-1}}$、$E_2=360.9\ \mathrm{cm^{-1}}$，与实验转动能级（$118.5$、$354.0\ \mathrm{cm^{-1}}$，计入离心畸变后下移）的偏差在一阶刚转子近似预期之内。该实验演示了"角动量代数 + DVR + 谱定理"三位一体的最小转动能级计算。

**练习。** (1) 证明 Wigner D 函数的正交归一关系。(2) 对不对称陀螺写出 $\langle JK|\hat{J}_a^{2}|J,K\pm2\rangle$ 的显式表达式并讨论 $K$ 结构。(3) 由 $D^{J}_{MK}$ 的完备性推导立体角积分化为 D 函数耦合系数的公式。

### 3.2 Jacobi 坐标与三原子反应几何

**定义。** 对 $A+BC$ 型三原子体系，质心分离后取原子–双原子 Jacobi 坐标 $(r,R,\gamma)$：$r$ 为双原子核间距，$R$ 为原子 $A$ 至双原子质心的距离，$\gamma$ 为 $\mathbf{r}$ 与 $\mathbf{R}$ 的夹角。另一组内禀坐标是三条核间距 $(r_{12},r_{23},r_{31})$，势能面天然对称地表达于此。二者存在封闭的解析双向变换。

**公式与推导。** 记双原子质心到 $B$、$C$ 的距离为 $d_b=\dfrac{m_c}{m_b+m_c}r$ 与 $d_c=\dfrac{m_b}{m_b+m_c}r$，由余弦定理

$$r_{AB}^{2}=R^{2}+d_b^{2}-2Rd_b\cos\gamma,\qquad
r_{AC}^{2}=R^{2}+d_c^{2}+2Rd_c\cos\gamma.$$

逆变换中 $R^{2}$ 由 Stewart 定理给出。推导如下：将 $m_b d_b=m_c d_c$ 代入加权平均

$$\frac{m_b r_{AB}^{2}+m_c r_{AC}^{2}}{m_b+m_c}
=R^{2}+\frac{m_b d_b^{2}+m_c d_c^{2}}{m_b+m_c}
+2R\cos\gamma\,\frac{-m_b d_b+m_c d_c}{m_b+m_c},$$

交叉项因 $m_b d_b=m_c d_c$ 而严格消失，且 $m_b d_b^{2}+m_c d_c^{2}=d_bd_c(m_b+m_c)$，故

$$R^{2}=\frac{m_b r_{AB}^{2}+m_c r_{AC}^{2}}{m_b+m_c}-d_b d_c,$$

$\gamma$ 则由 $\cos\gamma=(R^{2}+d_b^{2}-r_{AB}^{2})/(2Rd_b)$ 恢复。转振动能算符在 Jacobi 坐标下为

$$\hat{T}=-\frac{\hbar^{2}}{2\mu_r}\frac{1}{r}\frac{\partial^{2}}{\partial r^{2}}r
-\frac{\hbar^{2}}{2\mu_R}\frac{1}{R}\frac{\partial^{2}}{\partial R^{2}}R
+\frac{(\hat{\mathbf{J}}-\hat{\mathbf{j}})^{2}}{2\mu_R R^{2}}
+\frac{\hat{\mathbf{j}}^{2}}{2\mu_r r^{2}},$$

其中 $\hat{\mathbf{j}}$ 为双原子转动角动量，$\hat{\mathbf{J}}$ 为总角动量，耦合项 $(\hat{\mathbf{J}}-\hat{\mathbf{j}})^{2}$ 展开后产生轨道–转动耦合，是反应散射密耦方程的中心结构。
该耦合的显式结构值得展开：$(\hat{\mathbf{J}}-\hat{\mathbf{j}})^{2}=\hat{J}^{2}+\hat{j}^{2}-2\,\hat{\mathbf{J}}\cdot\hat{\mathbf{j}}$，标量积以 $\hat{J}_z\hat{j}_z+\tfrac12(\hat{J}_{+}\hat{j}_{-}+\hat{J}_{-}\hat{j}_{+})$ 作用于体固定基，产生 $\Delta\Omega=0,\pm1$ 的螺旋性阶梯耦合，其矩阵元为 $\sqrt{J(J+1)-\Omega(\Omega\pm1)}$ 型的纯代数因子。势能面亦按 Legendre 级数 $V(r,R,\gamma)=\sum_{\lambda}V_{\lambda}(r,R)P_{\lambda}(\cos\gamma)$ 展开，展开系数由 Gauss–Legendre 求积获得，与 1.2 节的 DVR 求积完全同构；耦合矩阵元 $\langle j\Omega|P_{\lambda}|j'\Omega'\rangle$ 再次化为 3j 符号网络。于是整条反应散射计算链——坐标变换、势能展开、角动量耦合、通道推进——的每一个环节都落在前两章建立的工具之上。

**物理含义。** 每个反应通道（$A+BC$、$AB+C$、$AC+B$）有各自的 Jacobi 坐标集；反应路径本质上是三套坐标之间经由内禀核间距 $(r_{12},r_{23},r_{31})$ 的连续变形。London–Eyring–Polanyi–Sato（LEPS）解析势能面以 Morse 单重态与反 Morse 三重态曲线构造库仑积分 $Q_i$ 与交换积分 $J_i$：

$$V=Q_1+Q_2+Q_3-\sqrt{\tfrac12\bigl[(J_1-J_2)^{2}+(J_2-J_3)^{2}+(J_3-J_1)^{2}\bigr]},$$

在 $\mathrm{H}+\mathrm{H}_2$ 基准体系上复现解离渐近与共线鞍点。此外，环绕锥形交叉的核置换回路的 Berry 几何相位（$\Phi_B=\pi$）亦表达在 $(r_{12},r_{23},r_{31})$ 空间的拓扑结构中（本库 `calc_berry_phase_around_ci`）。

**数值陷阱。** 其一，反三角函数定义域：数值误差可使 $\cos\gamma$ 越出 $[-1,1]$，本库以截断处理，但该截断会掩盖上游误差的累积，规范做法是同时监测截断频率。其二，共线构型 $\gamma=0,\pi$ 处 $\mathbf{r}\times\mathbf{R}=0$，角动量耦合的体固定展开出现奇异性，弯曲基组须以 Legendre 或超球谐函数表达。其三，核置换对称性：同核体系（如 $\mathrm{H}_3$）要求势能面与波函数对置换群不变，坐标变换若未与置换算符一致排序，对称化误差会以 $10^{-3}$ 量级污染能级。

**GeneralModule 实现映射。** 源码 `src/mod_triatomic_geometry.f90`：`jacobi_to_internuclear` 与 `internuclear_to_jacobi`（双向解析变换，含 $\cos\gamma$ 域截断）、`calc_leps_potential` 与 `init_default_h3_leps`（标准 $\mathrm{H}_3$ LEPS 面参数：$D_e=0.1744\ E_h$，$r_e=1.401\ a_0$，$\beta=1.044\ a_0^{-1}$，Sato 参数 $\Delta=0.10$）、`calc_conical_intersection_adiabats` 与 `calc_berry_phase_around_ci`。文献依据：Truhlar 与 Horowitz（J. Chem. Phys. 68, 2466 (1978), DOI: 10.1063/1.436019）、Sato（J. Chem. Phys. 23, 592 (1955), DOI: 10.1063/1.1742050）、Berry（Proc. R. Soc. Lond. A 392, 45 (1984), DOI: 10.1098/rspa.1984.0023）。GitHub 直链：[src/mod_triatomic_geometry.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_triatomic_geometry.f90)。二维码：`qr/src__mod_triatomic_geometry.f90.png`。

![三原子几何模块二维码](qr/src__mod_triatomic_geometry.f90.png)

**小型数值实验。** 三项校验：(a) 往返一致性：随机取 $(r,R,\gamma)=(1.4,3.2,1.1)\ (a_0,\ a_0,\ \mathrm{rad})$，经 `jacobi_to_internuclear` 再经 `internuclear_to_jacobi` 返回，三个坐标的恢复误差分别为 $0$、$4.4\times10^{-16}$、$0$，达到机器精度；(b) LEPS 双原子渐近：沿 $r_{23}=r_{31}=8\ a_0$ 扫描 $r_{12}$，极小值位于 $r_{12}=1.4011\ a_0$，深度 $V_{\min}=-0.174359\ E_h$，与输入 Morse 参数 $D_e,r_e$ 精确一致；(c) 共线鞍点：在共线平面 $r_{31}=r_{12}+r_{23}$ 上以二维 Newton 法求解 $\nabla V=0$，得鞍点 $r_{12}=r_{23}=1.775\ a_0$（$r_{31}=3.55\ a_0$），$V^{\ddagger}=-0.15371\ E_h$，相对 $\mathrm{H}_2$ 极小的势垒高度为 $12.96\ \mathrm{kcal/mol}$；文献中 London 面（未加 Sato 修正）的参考值约为 $9.8\ \mathrm{kcal/mol}$，差异来自 Sato 参数 $\Delta$ 的选择，属模型系统的固有性质而非数值误差。

**练习。** (1) 完成 Stewart 定理中交叉项消失与 $m_b d_b^{2}+m_c d_c^{2}=d_bd_c(m_b+m_c)$ 的详细代数。(2) 推导 Jacobi 坐标下动能算符的显式形式并指出 $\gamma$ 依赖项。(3) 以本库 LEPS 面扫描 Sato 参数 $\Delta\in[0.05,0.20]$，讨论势垒高度对 $\Delta$ 的单调性。

### 3.3 Delves 质量标度超球坐标

**定义。** 将 Jacobi 坐标作质量标度变换

$$S=d\,R,\qquad s=\frac{r}{d},\qquad d=\left(\frac{\mu_{A,BC}}{\mu_{BC}}\right)^{1/4},$$

定义超半径 $\rho=\sqrt{S^{2}+s^{2}}$ 与超角 $\alpha=\arctan(s/S)$。三维体系共六个内部加转动自由度，完整超球坐标为 $(\rho,\alpha,\theta,\phi)$ 加整体转动，其中广义角动量算符 $\hat{\Lambda}^{2}$（grand angular momentum）生成五维超球面上的 Laplace–Beltrami 算符。
$\hat{\Lambda}^{2}$ 的本征函数是超球谐函数，按超角动量量子数 $K$ 组织，简并度随 $K$ 以四次多项式增长，故基组截断必须与物理通道数匹配以避免维度爆炸。该基组的物理优越性在于：势能各向异性随 $\rho$ 演化——渐近区退化为二体分波结构，强相互作用区各通道自然混合——而超球基组恰好按这一物理演化组织耦合矩阵，使非对角耦合集中于物理上真正相关的通道之间。动能取各向同性的紧凑形式

$$\hat{T}=-\frac{\hbar^{2}}{2\mu}\left(\frac{\partial^{2}}{\partial\rho^{2}}+\frac{5}{\rho}\frac{\partial}{\partial\rho}-\frac{\hat{\Lambda}^{2}}{\rho^{2}}\right),
\qquad \mu=\sqrt{\mu_R\mu_r}=\sqrt{\frac{m_A m_B m_C}{M}}.$$

**公式与推导。** 标度因子 $d$ 的选取条件是动能的各向同性：以 $S=dR$、$s=r/d$ 变换动量 $P_S=P_R/d$、$P_s=P_r\,d$，动能化为 $T=P_S^{2}d^{2}/(2\mu_R)+P_s^{2}/(2\mu_r d^{2})$；令 $d^{4}=\mu_R/\mu_r$ 则两项系数相等并等于 $1/(2\mu)$，其中 $\mu=\mu_R/d^{2}=\sqrt{\mu_R\mu_r}$。恒等式 $\mu_R\mu_r=\dfrac{m_A(m_B+m_C)}{M}\cdot\dfrac{m_Bm_C}{m_B+m_C}=\dfrac{m_A m_B m_C}{M}$ 给出三体约化质量的紧凑表达，与库内 `init_reaction_mass` 的实现一致。不同排列通道（$A+BC$、$AB+C$）的 Jacobi 坐标系在 $(S,s)$ 平面内相差一个由质量决定的刚性旋转，通道轴之间的夹角即反应偏角（reaction skew angle）：

$$\tan\beta_{\mathrm{skew}}=\sqrt{\frac{m_B M}{m_A m_C}},\qquad
\sin\beta_{\mathrm{skew}}=\sqrt{\frac{m_B M}{(m_A+m_B)(m_B+m_C)}},\qquad
\cos\beta_{\mathrm{skew}}=\sqrt{\frac{m_A m_C}{(m_A+m_B)(m_B+m_C)}}.$$

三个表达式自洽：$\sin^{2}\beta+\cos^{2}\beta=\dfrac{m_B(m_A+m_B+m_C)+m_A m_C}{(m_A+m_B)(m_B+m_C)}=\dfrac{(m_A+m_B)(m_B+m_C)}{(m_A+m_B)(m_B+m_C)}=1$，代数上严格成立。

**物理含义。** 超球坐标把三原子体系的全部构型纳入一个统一的 $(\rho,\alpha)$ 平面：$\rho$ 度量整体尺寸，$\alpha$ 在三个通道扇区之间连续插值，反应即沿 $\rho$ 的推进与 $\alpha$ 的偏转。同核反应 $\mathrm{H}+\mathrm{H}_2\to\mathrm{H}_2+\mathrm{H}$ 的偏角恰为 $60^{\circ}$，三通道在超球面上完全对称，这是超球方法处理反应散射与 Efimov 物理的根本优势；本库 `mod_three_body_recombination.f90` 的 Efimov 模块与 `mod_hyperspherical_reactive.f90` 的过渡态速率模块共享同一套质量标度几何。
超球方法的核心数值对象是固定 $\rho$ 处的表面绝热通道：将 $\hat{\Lambda}^{2}/(2\mu\rho^{2})+V(\rho,\Omega)$ 在超角基上对角化，得到绝热势曲线 $U_{\nu}(\rho)$ 与绝热通道函数 $\Phi_{\nu}(\rho;\Omega)$；沿 $\rho$ 推进时通道间由非对角导数耦合 $\langle\Phi_{\nu'}|\partial_\rho\Phi_{\nu}\rangle$ 连接，Johnson 的对角化修正可将其中的纯几何项解析扣除，使剩余耦合反映真实的非绝热物理。广义角动量基的截断 $K_{\max}$ 与通道数 $N_{\mathrm{ch}}$ 的收敛检验同样遵循 1.2 节的变分逻辑：通道数增加时累积反应几率 $N(E)$ 单调收敛，为速率常数的误差棒提供了方向性依据。

**数值陷阱。** 其一，$\rho=0$（三体聚心）是动能的真正奇点，边界条件须按超球谐展开的正则性处理；普通网格方法应使网格内边界远离该点并检验能量对内边界位置的不敏感性。其二，质量极不对称体系（如 $\mathrm{Mu}+\mathrm{H}_2$）的偏角趋近极端值，通道扇区高度压缩，$K$（超角动量）基组收敛显著变慢，宜采用可变超角基或绝热通道对角化。其三，标度因子 $d$ 与通道约化质量在不同文献中存在倒数约定之别，调用跨库代码时第一项校验应是往返变换与偏角数值。

**GeneralModule 实现映射。** 源码 `src/mod_hyperspherical_reactive.f90`：`init_reaction_mass`（$d$、$\mu$、$\beta_{\mathrm{skew}}$ 的自动生成）、`jacobi_to_hyperspherical` 与 `hyperspherical_to_jacobi`（双向变换）、`calc_eckart_transmission`、`calc_cumulative_reaction_probability`、`calc_canonical_rate_constant`、`calc_tst_wigner_rate`。文献依据：B. R. Johnson（J. Chem. Phys. 73, 5051 (1980), DOI: 10.1063/1.440058）、R. T. Pack 与 G. A. Parker（J. Chem. Phys. 87, 3888 (1987), DOI: 10.1063/1.452944）。GitHub 直链：[src/mod_hyperspherical_reactive.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_hyperspherical_reactive.f90)。二维码：`qr/src__mod_hyperspherical_reactive.f90.png`。

![超球反应模块二维码](qr/src__mod_hyperspherical_reactive.f90.png)

**小型数值实验。** (a) 对 $\mathrm{H}+\mathrm{H}_2$（三质量均为 $1.0078\ \mathrm{amu}$）调用 `init_reaction_mass`，得偏角 $60.000000^{\circ}$，与解析值一致；(b) 对 $\mathrm{D}+\mathrm{H}_2$（$m_A=2.0141$，$m_B=m_C=1.0078\ \mathrm{amu}$）得 $54.740638^{\circ}$；若以理想整数质量 $2:1:1$ 计算，则解析极限为 $\arccos(1/\sqrt3)=54.735610^{\circ}$，微差来自同位素质量的非严格整数比，说明偏角对质量比连续敏感；(c) 以返回的 $d$ 作超球往返变换 $(r,R)=(1.4,3.2)\ a_0$，恢复误差 $0$ 与 $4.4\times10^{-16}$。

**练习。** (1) 完成动能各向同性条件 $d^{4}=\mu_R/\mu_r$ 的推导。(2) 证明偏角表达式满足 $\sin^{2}\beta+\cos^{2}\beta=1$。(3) 对 $\mathrm{H}_3$ 计算等边三角形构型的超半径与超角，并验证三个通道扇区的等价性。

### 3.4 Watson 转振哈密顿量、Coriolis 耦合、离心畸变与振动角动量

**定义。** 在满足 Eckart 条件的体固定系中，精确的振转哈密顿量由 Watson 给出：

$$\hat{H}_{\mathrm{vr}}=\frac{1}{2}\sum_{\alpha\beta}\mu_{\alpha\beta}\bigl(\hat{J}_{\alpha}-\hat{\pi}_{\alpha}\bigr)\bigl(\hat{J}_{\beta}-\hat{\pi}_{\beta}\bigr)
+\frac{1}{2}\sum_{k}\hat{p}_{k}^{2}+V(\mathbf{Q})
-\frac{\hbar^{2}}{8}\sum_{\alpha}\mu_{\alpha\alpha},$$

其中 $\mu_{\alpha\beta}$ 为有效惯量逆张量（简正坐标的函数），$\hat{J}_{\alpha}$ 为体固定角动量分量，$\hat{\pi}_{\alpha}=\sum_{kl}\zeta^{\alpha}_{kl}Q_k\hat{p}_l$ 为振动角动量，$\zeta^{\alpha}_{kl}$ 为 Coriolis ζ 常数（对 $k,l$ 反对称），末项为 Watson 伪势。Eckart 条件

$$\sum_i m_i\mathbf{a}_i=\mathbf{0},\qquad \sum_i m_i\,\mathbf{r}_{i,e}\times\mathbf{a}_i=\mathbf{0}$$

（$\mathbf{a}_i$ 为位移，$\mathbf{r}_{i,e}$ 为平衡位形）保证振动与转动的动量耦合在平衡位形处线性消失。

**公式与推导。** 经典动能 $T=\frac12\sum_{\alpha\beta}\mu_{\alpha\beta}(J_{\alpha}-\pi_{\alpha})(J_{\beta}-\pi_{\beta})+\frac12\sum_k p_k^{2}$ 的展开给出三类耦合：其一，$\mu_{\alpha\beta}$ 对 $Q$ 的线性展开产生离心畸变项，对不对称陀螺按约化处理（A 或 S reduction）写入四次幂算符

$$\hat{H}_{\mathrm{cd}}=-D_J\hat{J}^{4}-D_{JK}\hat{J}^{2}\hat{J}_z^{2}-D_K\hat{J}_z^{4}+d_1\hat{J}^{2}(\hat{J}_{+}^{2}+\hat{J}_{-}^{2})+d_2(\hat{J}_{+}^{4}+\hat{J}_{-}^{4});$$

对双原子分子，$D_v=4B_v^{3}/\omega_v^{2}$ 可由转动–振动二级微扰严格导出：$\hat{H}_{\mathrm{rot}}=B(\hat{R})\hat{J}^{2}$ 中 $B(\hat{R})$ 的非对角部分通过 $\omega_v$ 量级的振动能隙与转动能级虚耦合，二级位移整理为 $-D_vJ^{2}(J+1)^{2}$。其二，$-\hat{J}_{\alpha}\hat{\pi}_{\beta}$ 型交叉项即 Coriolis 耦合，使振动角动量 $l$ 与整体转动发生相干混合，在简并弯曲振动中产生 $l$ 型倍频分裂。其三，$\hat{\pi}_{\alpha}\hat{\pi}_{\beta}$ 项给出振动的非谐与角动量耦合修正。零级近似下能级即 $E(v,J)=E_{\mathrm{vib}}(v)+B_vJ(J+1)$，而

$$B_v=\left\langle\chi_v\left|\frac{\hbar^{2}}{2\mu R^{2}}\right|\chi_v\right\rangle
=B_e-\alpha_e\left(v+\tfrac12\right)+\cdots$$

把振动平均转动常数与振转耦合常数 $\alpha_e$ 联系起来。
$D_v=4B_v^{3}/\omega_v^{2}$ 的推导可作微扰方法的范本：把 $\hat{H}_{\mathrm{rot}}=B(\hat{R})\hat{J}^{2}$ 写成 $B_v\hat{J}^{2}+\bigl[B(\hat{R})-B_v\bigr]\hat{J}^{2}$，后者对角部分仅重整化 $B_v$，非对角部分以振动矩阵元 $\langle v\pm1|B(\hat{R})-B_v|v\rangle\approx\mp\alpha_e\sqrt{(v+1)/2}$ 与转动能差 $\approx\pm\omega_v$ 代入二级微扰，逐项求和后恰整理为 $-4B_v^{3}\omega_v^{-2}\,J^{2}(J+1)^{2}$。同一逻辑给出多原子离心畸变常数与 Coriolis 型分裂：线性分子简并弯曲态的 $l$ 型倍频源于振动角动量项与 $K$ 结构的联合作用，量级为 $q\sim B^{2}/\omega$。Watson 展开的有效性判据是 $D_v J^{2}(J+1)^{2}\ll B_v$ 与 $|\zeta|J\ll1$；超转子与 floppy 体系必须放弃微扰、返回精确动能算符做全变分处理。

**物理含义。** Watson 哈密顿量是高分辨转动光谱拟合的标准模型（分子常数 $A,B,C,D_J,D_{JK},\ldots,\zeta$ 的物理载体）。离心畸变度量化学键的转动软化：转得越快，键被离心力拉伸，$B$ 下降，能级相对刚转子逐级下压；Coriolis 耦合则是转动能级内振动角动量再分配的通道，是红外–微波双共振与 $l$ 型倍频光谱的核心机制。对范德华络合物与其它 floppy 体系，微扰展开失效，须回到精确动能算符做变分处理——这正是 DVR 与密耦方法的价值所在。

**数值陷阱。** 其一，度量约定：FGH/DVR 求解器返回的本征矢系数在格点 Kronecker 度量下归一化，而 quadrature 型积分函数（如 `calc_rotational_constants_bv`）期望物理波函数 $\psi=z/\sqrt{\Delta x}$；直接以原始系数调用会得到缩小 $\Delta x$ 倍的 $B_v$（下述实验实测 $0.6726$ 对 $60.18\ \mathrm{cm^{-1}}$）。其二，Watson 伪势 $-\hbar^{2}\sum_\alpha\mu_{\alpha\alpha}/8$ 虽小，对轻氢化物可达 $0.1\ \mathrm{cm^{-1}}$ 量级，高分辨拟合中不可忽略。其三，A reduction 在近球形陀螺极限退化，应改用 S reduction；ζ 常数的符号约定在 Wilson–Decius–Cross 与 Watson 两套文献中相反。其四，微扰离心畸变常数在低频模体系（范德华模、弯曲模）失效，误差可达一个量级。

**GeneralModule 实现映射。** 源码 `src/mod_rovibrational.f90`：`build_rovibrational_hamiltonian`（零级转振能级组装）、`calc_rotational_constants_bv`（振动平均转动常数 $\langle\chi_v|\hbar^{2}/(2\mu R^{2})|\chi_v\rangle$）、`calc_vibrational_dipole_matrix` 与 `calc_franck_condon_factors`；`src/mod_dvr_grid.f90`：`dvr_sinc_init` 与 `fgh_solve_bound_states`（提供 $\chi_v$）。文献依据：E. B. Wilson, Jr., J. C. Decius 与 P. C. Cross, *Molecular Vibrations*, McGraw-Hill, New York (1955)；J. K. G. Watson, Mol. Phys. 15, 479 (1968)；Bunker 与 Jensen (1998)。GitHub 直链：[src/mod_rovibrational.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rovibrational.f90)。二维码：`qr/src__mod_rovibrational.f90.png`。

![转振模块二维码](qr/src__mod_rovibrational.f90.png)

**小型数值实验。** 对 Morse 型 $\mathrm{H}_2$（$D_e=0.1744\ E_h$，$r_e=1.401\ a_0$，$\beta=1.044\ a_0^{-1}$，$\mu=0.5\ \mathrm{amu}$），在 $[0.4,6.0]\ a_0$ 网格（$n=500$）上以 Sinc-DVR FGH 求得振动本征态：基本带间隔 $E_1-E_0=4219.92\ \mathrm{cm^{-1}}$；盒内束缚能级 15 条（半经典估计 17 条，最高两条的外转折点越出盒界，演示盒子误差）。以 $\psi=z/\sqrt{\Delta x}$ 正确归一化后调用 `calc_rotational_constants_bv` 得 $B_0=60.18$、$B_1=57.71$、$B_2=55.12$、$B_3=52.39$、$B_4=49.52\ \mathrm{cm^{-1}}$，呈单调下降，拟合得 $\alpha_e=B_0-B_1=2.46\ \mathrm{cm^{-1}}$（实验 $\mathrm{H}_2$ 为 $3.06\ \mathrm{cm^{-1}}$，差异源于 Morse 参数与真实势能面的偏离）；离心畸变估计 $D_J\approx4B_0^{3}/(E_1-E_0)^{2}=0.04895\ \mathrm{cm^{-1}}$，与实验值 $0.047\ \mathrm{cm^{-1}}$ 偏差约 $4\%$，验证了 Wilson–Decius–Cross 微扰公式的定量可靠性。若跳过 $\sqrt{\Delta x}$ 归一化直接调用，$B_0$ 输出为 $0.6726\ \mathrm{cm^{-1}}$，恰为正确值乘以 $\Delta x=0.011196$，即 1.1 节所述度量陷阱的实例。

**练习。** (1) 对 Morse 振子推导 $D_v=4B_v^{3}/\omega_v^{2}$ 并指出二级微扰的适用条件。(2) 证明 Eckart 条件使动能中 $\omega\cdot\sum_i m_i\mathbf{r}_{i,e}\times\mathbf{a}_i$ 型线性 Coriolis 项消失。(3) 以本库 $B_v$ 序列拟合 $B_v=B_e-\alpha_e(v+\frac12)+\gamma_e(v+\frac12)^{2}$，报告 $B_e,\alpha_e,\gamma_e$ 并讨论三参数拟合的残差结构。

---

## 第四章 单位制与本征求解的数值分析

理论正确而数值失真，是计算物理中最隐蔽的一类失败。本章讨论三层数值基础：单位制与量纲一致性（4.1 节）、实对称本征问题的 Householder 三对角化与 QL 隐式位移迭代（4.2、4.3 节）、以及贯穿一切算法的条件数、误差传播与收敛阶理论（4.4 节）。GeneralModule 的设计原则——零外部库依赖、统一 `real(dp)` 强类型、纯函数契约——正是在这三层基础上确立的。
与前三章不同，本章的主题不依赖具体物理体系，而是普适的数值工程准则；但其全部结论仍以物理算例呈现：单位制的自洽以组合常数交叉验证，本征求解器以解析谱验收，误差理论以随机矩阵与差分实验定标。这种以已知答案校验未知计算的方法论，与 1.3 节的谱定理验收、2.2 节的正交性校验一脉相承，构成本库持续集成测试体系的理论基础。

### 4.1 单位制、原子单位与量纲一致性

**定义。** 计算量子力学并行使用国际单位制（SI）与原子单位制（a.u.）。原子单位定义为 $\hbar=e=m_e=4\pi\varepsilon_0=1$，由此导出能量单位 Hartree

$$E_h=\frac{m_e e^{4}}{(4\pi\varepsilon_0)^{2}\hbar^{2}}=4.3597447\times10^{-18}\ \mathrm{J}
=27.211386\ \mathrm{eV}=219474.63\ \mathrm{cm^{-1}},$$

长度单位 Bohr 半径 $a_0=0.529177\ \text{\AA}$，时间单位 $\hbar/E_h=2.418884\times10^{-17}\ \mathrm{s}$。质量以电子质量计，$1\ \mathrm{amu}=1822.888486\ m_e$（库内 `AMU2AU`）。

**公式与推导。** 原子单位下电子的薛定谔方程化为无量纲形式 $[-\frac12\nabla^{2}-\frac{Z}{r}]\psi=E\psi$，数值量级集中于 $\mathcal{O}(1)$，是浮点运算最友好的区间。单位换算的内部自洽性可由组合常数交叉验证：以库内常数 $A_{\mathrm{eV}}=$ `AU2EV` 与 $A_{\mathrm{cm}}=$ `AU2CM` 相除，得

$$\frac{A_{\mathrm{eV}}}{A_{\mathrm{cm}}}=1.239842\times10^{-4}\ \mathrm{eV\cdot cm},$$

这正是组合 $hc=1.239841984\times10^{-4}\ \mathrm{eV\cdot cm}$（CODATA）的数值，两个独立存储的换算因子在八位有效数字内自洽。同理，光强换算因子 `AU2W_CM2`= $3.5094452\times10^{16}\ \mathrm{W/cm^{2}}$ 由 $I=\frac12\varepsilon_0 c E_0^{2}$ 以 $E_0=$ `AU2VM` $=5.14221\times10^{11}\ \mathrm{V/m}$ 代入即得，磁场换算 `AU2TESLA`= $2.3505176\times10^{5}\ \mathrm{T}$ 由玻尔磁子对应关系导出。温度以能量计的约定：库内 `K2AU` 把以开尔文计的温度值转换为 $k_B T$ 的 Hartree 数（$1\ \mathrm{K}\leftrightarrow3.1578\times10^{-5}\ k_B\cdot\mathrm{K}/E_h$ 量级），故 Boltzmann 因子直接写为 $\exp(-E/(T\cdot\mathrm{K2AU}))$。
量纲分析还能给出超越直接计算的解析标度。以超冷散射为例，范德华势 $-C_6/r^{6}$ 下唯一可构造的长度量纲组合是 $(2\mu C_6/\hbar^{2})^{1/4}$，Gribakin–Flambaum 的半经典分析进一步给出平均散射长度

$$\bar{a}=\frac{2\pi}{\Gamma(1/4)^{2}}\left(\frac{2\mu C_6}{\hbar^{2}}\right)^{1/4}\approx0.4779888\left(\frac{2\mu C_6}{\hbar^{2}}\right)^{1/4},$$

其中数值常数仅由无量纲的 WKB 作用量相位决定。这类结果的启示是：凡在单位换算后出现反常的大数或小数，首先应怀疑某个隐藏的特征尺度未被约化。本库各模块以 a.u. 为统一内部单位，正是为了让此类标度检验可以直接进行。

**物理含义。** 单位制选择即物理问题的尺度选择：原子与分子过程以 Hartree、Bohr、飞秒为自然尺度；光谱学惯用波数 $\mathrm{cm^{-1}}$；动力学速率用 $\mathrm{cm^{3}\,mol^{-1}\,s^{-1}}$ 或 $\mathrm{cm^{3}\,s^{-1}}$。本库统一以 a.u. 为内部计算单位，以字符串参数 `to_au`/`from_au` 在边界处换算，使核心算法与单位策略解耦。
单位一致性的工程价值在于其可检验性：任何以不同推导路径得到的同一物理量必须在机器精度内重合，例如由能级差换算的波数与直接以波数输入的跃迁频率、由速率常数积分还原的截面与直接计算的截面。本库的若干模块内置了此类双向校验（如 Breit–Rabi 解析式与数值对角化的自动比对），把单位问题从程序员纪律提升为可自动执行的断言。建议读者在自建计算流程中遵循同一原则：每个涉及换算的接口都应配一条往返恒等式测试。

**数值陷阱。** 其一，单位混用是静默错误的头号来源：eV 与 Hartree 相差 27.2 倍、$\mathrm{cm^{-1}}$ 与 Hartree 相差 $2.2\times10^{5}$ 倍，程序照常运行但物理结论全非；对策是在每个模块入口强制换算并保留单位注释。其二，字符串换算函数的默认分支：`to_au` 对未识别的单位串不报错而原值返回，若拼写笔误（如 `cm-1` 写成 `cm^-1` 以外的变体）将被静默吞掉；调用后宜以量级断言抽查。其三，温度换算的双重含义：$T\cdot$`K2AU` 是 $k_BT$ 而非温度本身的能量数值，混用于 Maxwell–Boltzmann 权重时差一个 $k_B$。其四，约化质量须用核质量而非原子量整数，$\mu_{\mathrm{H_2}}=0.5\times1.0078=0.5039\ \mathrm{amu}$ 与 $0.5\ \mathrm{amu}$ 在高分辨比较时差 $0.8\%$。

**GeneralModule 实现映射。** 源码 `src/mod_constants.f90`：CODATA 基础常数（`C_LIGHT`、`HBAR`、`M_E`、`CHARGE_E`、`EPS0`、`KB`、`AMU2AU`）、数学常数（`PI`、`SQRTPI`、`EYE`）、约二十组双向换算因子（`AU2EV`/`EV2AU`、`AU2CM`/`CM2AU`、`AU2FS`/`FS2AU`、`AU2DEBYE`、`AU2VM`、`AU2W_CM2`、`AU2TESLA` 等）与纯函数 `to_au`、`from_au`。GitHub 直链：[src/mod_constants.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90)。二维码：`qr/src__mod_constants.f90.png`。

![常数模块二维码](qr/src__mod_constants.f90.png)

**小型数值实验。** 两项自洽性检验：(a) 计算 `AU2EV/AU2CM`，得 $1.239842\times10^{-4}\ \mathrm{eV\cdot cm}$，与 CODATA 的 $hc$ 在第八位有效数字一致；(b) 以 `from_au` 将 $4351.6\ \mathrm{cm^{-1}}$ 换算为 a.u. 再以 `to_au` 换回，得 $4351.6000000000004$，往返误差 $4\times10^{-13}\ \mathrm{cm^{-1}}$，纯属双精度舍入。此类零成本检验建议纳入持续集成。

**练习。** (1) 由 $\hbar,e,m_e,\varepsilon_0$ 的量纲推导 $E_h$、$a_0$ 与 a.u. 磁场单位 $2.35\times10^{5}\ \mathrm{T}$ 的表达式。(2) 证明 `AU2W_CM2` 可由 `AU2VM`、`EPS0`、`C_LIGHT` 完整重构。(3) 以本库常数计算 $1\ \mathrm{K}$ 对应的波数 $\mathrm{(k_B/hc)}=0.695\ \mathrm{cm^{-1}}$ 并与 CODATA 比对。

### 4.2 实对称本征问题与 Householder 三对角化

**定义。** 实对称本征问题 $A z=\lambda z$（$A=A^{T}$）的现代标准解法是两阶段方法：第一阶段以正交相似变换把 $A$ 化为三对角形式 $T=Q^{T}AQ$；第二阶段以带隐式位移的 QL（或 QR）迭代求 $T$ 的谱。Householder 反射定义为

$$P=I-\frac{2\,vv^{T}}{v^{T}v},\qquad P=P^{T}=P^{-1},$$

即对称、正交、对合的变换矩阵。

**公式与推导。** 正交性可直接验证：$P^{T}P=\left(I-\frac{2vv^{T}}{v^{T}v}\right)^{2}=I-\frac{4vv^{T}}{v^{T}v}+\frac{4v(v^{T}v)v^{T}}{(v^{T}v)^{2}}=I$。给定列向量 $u$（欲消元的目标），取 $v=u+\mathrm{sign}(u_1)\|u\|\,\mathbf{e}_1$，则

$$Pu=-\mathrm{sign}(u_1)\|u\|\,\mathbf{e}_1+0\cdot\mathbf{e}_2+\cdots+0\cdot\mathbf{e}_n,$$

即一次反射即可把一列的 $n-1$ 个非零元素清零。对 $A$ 逐列（本库 `tred2` 自第 $n$ 列向第 2 列推进）施加嵌入式的反射 $P_i=\mathrm{diag}(I_{n-i},\,P')$（子块反射不动已完成的部分），$n-2$ 步后得三对角矩阵 $T$；反射的符号选择 $g=-\mathrm{sign}(f)\sqrt{h}$（$f$ 为对角元、$h$ 为模方）避免了 $v$ 构造中的相消。累积全部反射得正交变换矩阵 $Q=P_1P_2\cdots P_{n-2}$，使 $Z=QZ_T$ 直接作用于三对角阶段的特征矢量。运算量：三对角化 $\frac{4}{3}n^{3}$ flops（含矢量累积 $\mathcal{O}(n^{3})$ 附加），数值稳定性由每一步的正交性保证，无增长因子问题。
与备选方案的比较有助于理解该选择的必然性。经典 Gram–Schmidt 正交化在近线性相关情形迅速丧失正交性，改进格式仅部分缓解；Jacobi 旋转法直接对全矩阵对角化，误差性质优美（渐近二次收敛、正交性接近机器精度），但运算量常数约十倍于两阶段方法，仅在结构特殊情形有竞争力。Householder 路线的另一工程优势是内存访问模式规整：三对角化阶段以列为主的对称秩二更新具有良好的缓存局部性，与本库 Fortran 列主序的数组布局天然匹配。对只需能级不需波函数的场景（如大规模扫描中的中间步骤），跳过反射累积可再节省约一半运算量。

**物理含义。** DVR、有限基组、密耦通道等一切厄米哈密顿量离散化后的束缚态问题都归结于此：本征值即能级，本征矢即波函数系数。两阶段结构还允许按需取舍——只要能级时只做三对角化加无矢量 QL（$\mathcal{O}(n^{2})$ 加速一个量级），需要波函数时再回补累积变换。

**数值陷阱。** 其一，输入必须严格对称：组装循环应显式镜像，或在调用前以 $\max_{ij}|A_{ij}-A_{ji}|$ 断言。其二，本库 `tred2` 以绝对阈值 $10^{-35}$ 判断列向量是否为零向量（scale 判据），对元素量级悬殊（例如混合 a.u. 与 $\mathrm{cm^{-1}}$ 单位后）的矩阵，阈值语义随范数缩放失效——换算单位必须在组装之前完成。其三，正交性随 $n$ 缓慢退化：$Z^{T}Z$ 对角偏差约 $n\varepsilon$ 量级，$n=10^{4}$ 时达 $10^{-12}$，对超高维问题应改用分块、随机化或迭代精化方案。

**GeneralModule 实现映射。** 源码 `src/mod_linear_algebra.f90`：`tred2`（EISPACK TRED2 的现代化移植，含反射累积）、`pythag`（防溢出的 $\sqrt{a^{2}+b^{2}}$ 安全计算）、`diag_symmetric_matrix`（统一入口）。算法源流：B. T. Smith 等, *Matrix Eigensystem Routines — EISPACK Guide*, 2nd ed., Springer, Berlin (1976)；G. H. Golub 与 C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, Baltimore (2013)。GitHub 直链：[src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90)。二维码：`qr/src__mod_linear_algebra.f90.png`。

**小型数值实验。** 见 1.3 节实验：$n=300$ 随机实对称矩阵，最大特征残余 $4.4\times10^{-14}$、正交性偏差 $2.2\times10^{-14}$，均符合 $n\varepsilon\|A\|$ 的舍入标度，证实 `tred2` 累积变换的正交性达到后向稳定标准。

**练习。** (1) 对 $3\times3$ 实对称矩阵手工执行一次 Householder 反射并验证三对角化第一步。(2) 证明嵌入式反射 $P_i$ 保持已完成列的零结构。(3) 统计 `tred2` 的乘加次数，验证 $\frac{4}{3}n^{3}$ 的量级。

### 4.3 QL 隐式位移迭代、Wilkinson 位移与收敛理论

**定义。** 对三对角矩阵 $T$，QL 迭代取正交阵 $Q_k$ 使 $T_{k+1}=Q_k^{T}T_kQ_k$ 且 $Q_k^{T}$ 的首列为 $e_1$ 方向的平面旋转乘积（QL 与 QR 分别自矩阵的下端与上端消元，二者等价地收敛于同一谱）。隐式位移技巧把位移 $\mu$ 融入旋转序列而不显式形成 $T-\mu I$；Wilkinson 位移取尾端 $2\times2$ 子矩阵靠近 $d_n$ 的那个本征值：

$$\mu=\mathrm{eig}\begin{pmatrix} d_{n-1} & e_{n-1} \\ e_{n-1} & d_n \end{pmatrix}_{\text{closer to } d_n}.$$

**公式与推导。** 对称三对角矩阵的 QR/QL 迭代保持对称性与三对角性，且次对角元以 Wilkkinson 位移达到三阶收敛：每次迭代后 $|e_{n-1}|=\mathcal{O}(|e_{n-1}|^{3}/\mathrm{gap})$。当 $|e_m|\le\varepsilon\cdot\mathrm{tst1}$（$\mathrm{tst1}$ 为全局范数标尺）时该次对角元置零，矩阵分裂为两个独立三对角块（deflation），各块独立迭代。全过程仅含平面旋转，数值上正交性精确保持，本征值误差不超过 $\mathcal{O}(\varepsilon\|A\|)$ 的后向稳定界。
Wilkinson 位移的选择有其确定论依据：位移取尾端 $2\times2$ 子阵靠近 $d_n$ 的本征值时，可证明对称情形下 $|e_{n-1}|$ 一次迭代后收缩至其三次幂量级，且位移序列单调锁定于某条本征值附近，不存在停滞。与之并行的另一条路线是 Sturm 序列二分法：利用三对角矩阵顺序主子式的符号变化计数区间内本征值个数，能以 $\mathcal{O}(n)$ 代价二分定位任意目标谱段，特别适合只需要部分能级的大规模问题；求出能级后再以逆迭代补充本征矢。本库选择 QL 全谱方案，取其实现紧凑、谱与波函数一体输出，且对分子哈密顿的典型规模（$n\le10^{4}$）性能完全足够。总代价：仅本征值 $\mathcal{O}(n^{2})$，含本征矢 $\mathcal{O}(n^{3})$。本库 `tql2` 按 EISPACK TQL2 实现，逐条本征值最多迭代 60 次后置非零错误码 `ierr`，最终以插入排序输出升序谱。

**物理含义。** 能级按升序输出直接对应光谱学的基态到激发态排列；简并能级由迭代自动给出正交本征矢；谱的排序稳定性使相邻能级间隔（转动常数、非谐常数）的差分运算可靠。三阶收敛意味着每条能级一般仅需二至三次迭代，$n=10^{4}$ 的全谱对角化在现代单机上为分钟量级。

**数值陷阱。** 其一，错误码检查：`stat`/`ierr` 非零表示某条本征值 60 次迭代未收敛，静默忽略将把未收敛值当作能级使用；本库统一在返回参数中传出，调用方必须显式判别。其二，位移计算中的除法 $p=(d_{l+1}-d_l)/(2e_l)$ 在 $e_l\to0$ 时上溢，须配合 `pythag` 的防溢出缩放。其三，deflation 判据 $\varepsilon\cdot\mathrm{tst1}$ 以全局范数为标尺，对谱跨度悬殊（如同时含深束缚与高激发能级）的哈密顿量，小幅本征值的相对精度可能受限，必要时应分块处理。其四，末端插入排序为 $\mathcal{O}(n^{2})$，超大规模下可替换为归并排序，但须保持本征矢与本征值的同步交换。

**GeneralModule 实现映射。** 源码 `src/mod_linear_algebra.f90`：`tql2`（QL 隐式位移、Wilkinson 位移、deflation、升序排序），供 `diag_symmetric_matrix` 调用；测试基线见 `tests/test_dvr_grid.f90` 与 `tests/test_propagators.f90`。算法源流：J. H. Wilkinson, *The Algebraic Eigenvalue Problem*, Clarendon Press, Oxford (1965)；Golub 与 Van Loan (2013)。GitHub 直链：[src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90)。二维码：`qr/src__mod_linear_algebra.f90.png`。

**小型数值实验。** 一维 Laplacian 验收：取 $h=0.01$、$n=200$ 的三对角矩阵 $A=\mathrm{tridiag}(-1,2,-1)/h^{2}$，其精确谱为 $\lambda_k=2\bigl(1-\cos(k\pi/201)\bigr)/h^{2}$。调用 `diag_symmetric_matrix` 得 `stat=0`，全谱最大偏差 $2.5\times10^{-11}$，即 $6\times10^{-16}\|A\|$（$\|A\|=4/h^{2}=4\times10^{4}$），一次到位地演示了后向稳定界 $\mathcal{O}(\varepsilon\|A\|)$ 的实际含义：绝对误差随矩阵范数放大，相对精度始终为机器水平。

**练习。** (1) 证明一次显式 QL 步保持三对角性与对称性（描画 bulge 追赶过程）。(2) 按 Wilkinson 的分析说明为何对称三对角 QL 迭代在 Wilkinson 位移下三阶收敛。(3) 修改实验为 $n=1000$，验证最大偏差仍满足 $c\,\varepsilon\|A\|$ 并估计常数 $c$。

### 4.4 条件数、误差传播与收敛阶

**定义。** 函数求值的（相对）条件数定义为 $\kappa(x)=\left|\dfrac{x f'(x)}{f(x)}\right|$，度量输入扰动到输出的放大倍数。矩阵问题的条件数 $\kappa(A)=\|A\|\,\|A^{-1}\|$；对称本征值问题的绝对条件数则由 Weyl 定理给出最优上界：对扰动 $A\to A+\delta A$，

$$\max_i|\delta\lambda_i|\le\|\delta A\|_{2},$$

即对称谱的绝对条件数为 1，这是谱方法稳定性的理论根基。本征矢的条件数依赖能隙（Davis–Kahan $\sin\theta$ 定理）：$\sin\theta\le\|\delta A\|_{2}/\mathrm{gap}$。算法的收敛阶 $p$ 由 $\lim\|e_{k+1}\|/\|e_k\|^{p}=C>0$ 定义：Newton 法 $p=2$，Wilkinson 位移的三对角 QL 迭代 $p=3$，$m$ 阶数值求积与差分格式的截断误差 $\mathcal{O}(h^{m})$。

**公式与推导。** 误差传播的一阶展开：$\delta f=\sum_i\frac{\partial f}{\partial x_i}\delta x_i$，相对形式为 $\frac{\delta f}{f}=\sum_i\kappa_i\frac{\delta x_i}{x_i}$，多个 $\kappa_i\gg1$ 的环节串联时误差指数式放大。数值微分的经典分析给出截断与舍入的权衡：中心差分

$$f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{f'''(x)}{6}h^{2}+\mathcal{O}(h^{4})$$

的总误差约为 $\frac{|f'''|}{6}h^{2}+\frac{\varepsilon|f|}{h}$，极小化得最优步长

$$h^{*}=\left(\frac{3\varepsilon|f|}{|f'''|}\right)^{1/3},\qquad
\delta f_{\min}\sim\varepsilon^{2/3},$$

对双精度 $\varepsilon=2.2\times10^{-16}$ 有 $\varepsilon^{2/3}\approx3.7\times10^{-11}$。Richardson 外推可把相邻两个步长的结果按误差主项系数组合，将收敛阶提高一倍（$2\to4$）。

**物理含义。** 条件数告诉研究者哪些物理量本质上可精确预测：对称哈密顿的非简并能级（绝对条件数 1，稳健）；而近简并子空间的混合角、共振宽度、谱的高阶导数（条件数大，脆弱）则须以简并微扰论或联合拟合处理。收敛阶则是网格与步长选择的定量依据：从两次计算的误差比可实测收敛阶 $p=\ln(e(h)/e(h/10))/\ln 10$，据此外推至 $h\to0$ 的极限值。
补充两个定量参照。其一，病态问题的典型标本是 Hilbert 矩阵 $H_{ij}=(i+j-1)^{-1}$：其条件数随维数指数增长，$n=12$ 时已达 $10^{16}$ 量级，双精度求逆完全失真。凡涉及病态最小二乘的计算（如由谱线位置拟合分子常数），应改用奇异值分解或正则化方案，并以残差 $\|A\,x-b\|$ 监控。其二，后向稳定与向前精度是两个概念：`tql2` 保证计算谱是某个邻近矩阵 $A+\delta A$ 的精确谱（后向稳定），但单条能级的向前误差仍受条件数控制；对称问题因绝对条件数为 1 而两者相当，非对称广义问题（复标度哈密顿、含重叠矩阵的广义本征问题）则可能出现后向稳定而向前失真的局面，判别手段是 4.3 节的残差复核。

**数值陷阱。** 其一，病态矩阵求逆：Gauss–Jordan 全主元（本库 `inv_real_matrix`）只在 $\kappa(A)\ll1/\varepsilon$ 时可靠，$\kappa$ 接近 $10^{16}$ 时返回的逆完全失真，且 `stat` 未必报警；对策是以 $\|A A^{-1}-I\|$ 残差复核。其二，差分步长不能想当然取小：$h=10^{-8}$ 处中心差分误差已开始反转上升，最优点在 $10^{-5}\sim10^{-6}$ 量级（见实验）。其三，测得的收敛阶若持续低于设计阶（如二阶格式测得 $p\approx1.6$），优先怀疑边界处理或奇点污染，而非盲目加密网格。其四，简并或近简并谱的本征矢比较须用子空间投影（1.3 节），任何以单本征矢为对象的判据在能隙趋于零时失效。

**GeneralModule 实现映射。** 源码 `src/mod_linear_algebra.f90`：`inv_real_matrix`、`inv_complex_matrix`（全主元 Gauss–Jordan，奇异检测 `stat=-1`）、`pythag`（溢出安全）；`src/mod_interpolation.f90` 提供解析一阶/二阶导数以替代数值差分。理论依据：Golub 与 Van Loan (2013)，第 7 章；Wilkinson (1965)。GitHub 直链：[src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90)。二维码：`qr/src__mod_linear_algebra.f90.png`。

**小型数值实验。** 两项定标实验。(a) Weyl 界实测：$n=200$ 随机对称矩阵 $A$ 加扰动 $E$（元素均匀随机、幅度 $0.05$），$\|E\|_{F}=5.79$，对角化 $A$ 与 $A+E$ 比对谱得 $\max_i|\delta\lambda_i|=0.122$，比值 $0.021\le1$，与 Weyl 界 $\max_i|\delta\lambda_i|\le\|E\|_2\le\|E\|_F$ 相容（实际远小于 Frobenius 界，因 $\|E\|_2\approx0.12$ 体现随机矩阵的谱集中）。(b) 收敛阶与最优步长：对 $f(x)=\cos x$ 在 $x=1$ 处的中心差分，在截断区实测收敛阶 $p=2.0000$，误差在 $h^{*}=1.0\times10^{-5}$ 处取得极小 $1.1\times10^{-11}$，与理论 $h^{*}=(3\varepsilon|f|/|f'''|)^{1/3}\approx9\times10^{-6}$ 及 $\varepsilon^{2/3}\approx3.7\times10^{-11}$ 两者吻合；继续缩小 $h$ 至 $10^{-8}$ 以下误差反转上升，演示了舍入主导区。

**练习。** (1) 证明 Weyl 定理在 $2\times2$ 情形，并用 Rayleigh 商极小极大原理推广。(2) 推导中心差分的总误差函数并求最优步长。(3) 以 Richardson 外推由 $h$ 与 $h/2$ 两点外推 $f'(1)$，与精确值 $-\sin1$ 比较并实测四阶收敛。

---

## 参考文献

1. D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific, Singapore (1988). DOI: 10.1142/0270.
2. J. J. Sakurai and J. Napolitano, *Modern Quantum Mechanics*, 3rd ed., Cambridge University Press, Cambridge (2020).
3. C. Cohen-Tannoudji, B. Diu, and F. Laloë, *Quantum Mechanics*, Vols. I–II, Wiley, New York (1977).
4. M. Reed and B. Simon, *Methods of Modern Mathematical Physics I: Functional Analysis*, Academic Press, New York (1972).
5. P. R. Bunker and P. Jensen, *Molecular Symmetry and Spectroscopy*, 2nd ed., NRC Research Press, Ottawa (1998).
6. E. B. Wilson, Jr., J. C. Decius, and P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra*, McGraw-Hill, New York (1955).
7. J. K. G. Watson, "The vibration–rotation Hamiltonian of polyatomic molecules", Mol. Phys. 15, 479 (1968).
8. C. Eckart, "Some studies concerning rotating axes and polyatomic molecules", Phys. Rev. 47, 552 (1935). DOI: 10.1103/PhysRev.47.552.
9. G. Herzberg, *Molecular Spectra and Molecular Structure I: Spectra of Diatomic Molecules*, 2nd ed., Van Nostrand, New York (1950).
10. D. T. Colbert and W. H. Miller, J. Chem. Phys. 96, 1982 (1992). DOI: 10.1063/1.462125.
11. C. C. Marston and G. G. Balint-Kurti, J. Chem. Phys. 91, 3571 (1989). DOI: 10.1063/1.456888.
12. M. D. Feit, J. A. Fleck, Jr., and A. Steiger, J. Comput. Phys. 47, 412 (1982). DOI: 10.1016/0021-9991(82)90091-2.
13. R. Kosloff, J. Phys. Chem. 92, 2087 (1988). DOI: 10.1021/j100319a003.
14. J. C. Light and T. Carrington, "Discrete-variable representations and their utilization", Adv. Chem. Phys. 114, 263 (2000).
15. B. R. Johnson, J. Chem. Phys. 73, 5051 (1980). DOI: 10.1063/1.440058.
16. R. T. Pack and G. A. Parker, J. Chem. Phys. 87, 3888 (1987). DOI: 10.1063/1.452944.
17. D. G. Truhlar and C. J. Horowitz, J. Chem. Phys. 68, 2466 (1978). DOI: 10.1063/1.436019.
18. S. Sato, J. Chem. Phys. 23, 592 (1955). DOI: 10.1063/1.1742050.
19. M. V. Berry, Proc. R. Soc. Lond. A 392, 45 (1984). DOI: 10.1098/rspa.1984.0023.
20. G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, Baltimore (2013).
21. J. H. Wilkinson, *The Algebraic Eigenvalue Problem*, Clarendon Press, Oxford (1965).
22. B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema, and C. B. Moler, *Matrix Eigensystem Routines — EISPACK Guide*, 2nd ed., Lecture Notes in Computer Science 6, Springer, Berlin (1976).

---

*本章节所有源码引用基于仓库 main 分支当前版本；数值实验环境为 GNU Fortran 11.4（`-O2 -std=f2008`），双精度 `real64`。项目文件在本章节撰写过程中未作任何修改。*


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

以下二维码均指向 GeneralModule 仓库中的对应文件。二维码图像保存在 `qr/`，文件名规则是把源码路径中的 `/` 替换为 `__`。

| 内容 | 路径 | 二维码 |
|---|---|---|
| 仓库主页 | <https://github.com/l1Ha/QuantumGeneralModule> | `qr/repo_root.png` |
| README | [README.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/README.md) | `qr/README.md.png` |
| 配置指南 | [CONFIG_GUIDE.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/CONFIG_GUIDE.md) | `qr/CONFIG_GUIDE.md.png` |
| 文献映射 | [LITERATURE.md](https://github.com/l1Ha/QuantumGeneralModule/blob/main/LITERATURE.md) | `qr/LITERATURE.md.png` |
| 常数与单位 | [src/mod_constants.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90) | `qr/src__mod_constants.f90.png` |
| 特殊函数 | [src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90) | `qr/src__mod_special_functions.f90.png` |
| 线性代数 | [src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90) | `qr/src__mod_linear_algebra.f90.png` |
| DVR | [src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90) | `qr/src__mod_dvr_grid.f90.png` |
| 波包传播 | [src/mod_wavepacket_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_wavepacket_propagator.f90) | `qr/src__mod_wavepacket_propagator.f90.png` |
| 吸收边界 | [src/mod_absorbing_boundary.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_absorbing_boundary.f90) | `qr/src__mod_absorbing_boundary.f90.png` |
| Chebyshev | [src/mod_chebyshev_propagator.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_chebyshev_propagator.f90) | `qr/src__mod_chebyshev_propagator.f90.png` |
| Lindblad | [src/mod_open_quantum.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_open_quantum.f90) | `qr/src__mod_open_quantum.f90.png` |
| Krotov | [src/mod_optimal_control.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_optimal_control.f90) | `qr/src__mod_optimal_control.f90.png` |
| 定态散射 | [src/mod_ti_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_ti_scattering.f90) | `qr/src__mod_ti_scattering.f90.png` |
| 含时散射 | [src/mod_td_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_td_scattering.f90) | `qr/src__mod_td_scattering.f90.png` |
| 外场散射 | [src/mod_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_field_scattering.f90) | `qr/src__mod_field_scattering.f90.png` |
| 转振光谱 | [src/mod_rovibrational.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rovibrational.f90) | `qr/src__mod_rovibrational.f90.png` |
| 光碎片 | [src/mod_photofragment_flux.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_photofragment_flux.f90) | `qr/src__mod_photofragment_flux.f90.png` |
| FSSH | [src/mod_surface_hopping_fssh.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_hopping_fssh.f90) | `qr/src__mod_surface_hopping_fssh.f90.png` |
| Jacobi 坐标 | [src/mod_triatomic_geometry.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_triatomic_geometry.f90) | `qr/src__mod_triatomic_geometry.f90.png` |
| 超球反应 | [src/mod_hyperspherical_reactive.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_hyperspherical_reactive.f90) | `qr/src__mod_hyperspherical_reactive.f90.png` |
| 三体复合 | [src/mod_three_body_recombination.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_three_body_recombination.f90) | `qr/src__mod_three_body_recombination.f90.png` |
| CIR | [src/mod_confined_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_confined_scattering.f90) | `qr/src__mod_confined_scattering.f90.png` |
| Fano/CCR | [src/mod_autoionization_fano.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_autoionization_fano.f90) | `qr/src__mod_autoionization_fano.f90.png` |
| 交叉场 | [src/mod_crossed_field_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_crossed_field_scattering.f90) | `qr/src__mod_crossed_field_scattering.f90.png` |
| 旋量 BEC | [src/mod_spinor_bec.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_spinor_bec.f90) | `qr/src__mod_spinor_bec.f90.png` |
| 表面散射 | [src/mod_surface_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_scattering.f90) | `qr/src__mod_surface_scattering.f90.png` |
| Eley–Rideal | [src/mod_surface_reaction_er.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_reaction_er.f90) | `qr/src__mod_surface_reaction_er.f90.png` |
| 电子摩擦 | [src/mod_surface_electronic_friction.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_surface_electronic_friction.f90) | `qr/src__mod_surface_electronic_friction.f90.png` |
| 相对论原子 | [src/mod_relativistic_atomic.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_relativistic_atomic.f90) | `qr/src__mod_relativistic_atomic.f90.png` |
| RIXS | [src/mod_resonant_xray_scattering.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_resonant_xray_scattering.f90) | `qr/src__mod_resonant_xray_scattering.f90.png` |
| Penning | [src/mod_penning_associative_ionization.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_penning_associative_ionization.f90) | `qr/src__mod_penning_associative_ionization.f90.png` |
| 强场 NSDI | [src/mod_strong_field_nsdi.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_strong_field_nsdi.f90) | `qr/src__mod_strong_field_nsdi.f90.png` |
| 里德堡阻塞 | [src/mod_rydberg_blockade.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_rydberg_blockade.f90) | `qr/src__mod_rydberg_blockade.f90.png` |
| FGH 示例 | [examples/ex01_fgh_diatomic_bound_states.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex01_fgh_diatomic_bound_states.f90) | `qr/examples__ex01_fgh_diatomic_bound_states.f90.png` |
| Split 示例 | [examples/ex03_split_operator_1d.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex03_split_operator_1d.f90) | `qr/examples__ex03_split_operator_1d.f90.png` |
| TI/TD 示例 | [examples/ex07_scattering_wavefunctions_ti_td.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex07_scattering_wavefunctions_ti_td.f90) | `qr/examples__ex07_scattering_wavefunctions_ti_td.f90.png` |
| Feshbach 示例 | [examples/ex08_ultracold_feshbach_segmented.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex08_ultracold_feshbach_segmented.f90) | `qr/examples__ex08_ultracold_feshbach_segmented.f90.png` |
| FSSH 示例 | [examples/ex30_tully_surface_hopping.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/examples/ex30_tully_surface_hopping.f90) | `qr/examples__ex30_tully_surface_hopping.f90.png` |
| Python 可视化 | [python/pygenmod/visualizer.py](https://github.com/l1Ha/QuantumGeneralModule/blob/main/python/pygenmod/visualizer.py) | `qr/python__pygenmod__visualizer.py.png` |

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
# 第十二部分　基础代码能力与缺口

## 12.1　能力分层

GeneralModule 的基础代码能力可以分成四层：

| 层级 | 内容 | 当前状态 |
|---|---|---|
| L0 | 类型、常数、单位转换、基础 I/O | 已具备 |
| L1 | 特殊函数、稠密本征、FFT、插值 | 已具备 |
| L2 | DVR、波包传播、散射、开放系统 | 已具备 |
| L3 | 稀疏方法、并行化、高维张量、自动微分 | 部分或未具备 |
| L4 | 通用测试框架、基准套件、持续性能监控 | 部分具备 |

这一章只讨论 L0/L1 和少量 L2 的基础能力。目的不是罗列 API，而是说明哪些能力已经稳定、哪些还不完整、在科研工作流中应该怎样使用。

## 12.2　数值类型与常数模块

`mod_constants` 提供统一的双精度类型别名和物理常数：

```fortran
use, intrinsic :: iso_fortran_env, only: dp => real64, int32, int64
use mod_constants, only: PI, EYE, HBAR, AU2EV, EV2AU
```

核心约定是：

- 所有实数使用 `real(dp)`；
- 所有复数使用 `complex(dp)`；
- 避免混用 `real(8)`、`double precision` 和编译器私有种类；
- 单位转换集中在 `to_au` 和 `from_au`。

### 已有能力

- 物理常数：$\hbar$、$c$、$m_e$、$e$、$\varepsilon_0$、$k_B$；
- 能量单位：Hartree、eV、$\mathrm{cm^{-1}}$、K、J；
- 长度单位：Bohr、Å、nm、m；
- 时间单位：原子单位、fs、ps、s；
- 电场与磁场：V/m、MV/cm、Tesla、Gauss。

### 已知限制

- 单位接口使用字符串分派，未知单位会静默返回原值；
- 没有编译期单位类型；
- 没有向量/张量级别的单位检查；
- 没有自动收集当前使用单位的元数据。

### 建议的改进

引入一个显式状态码版本：

```fortran
pure subroutine to_au_checked(value, unit_name, value_au, stat)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: unit_name
    real(dp), intent(out) :: value_au
    integer, intent(out) :: stat
end subroutine to_au_checked
```

并增加单元测试：

$$
\left|
\operatorname{from\_au}(\operatorname{to\_au}(x,u),u)-x
\right|
\le
\epsilon |x|.
$$

### 对应源码

- [src/mod_constants.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_constants.f90)
- 二维码：

![](qr/src__mod_constants.f90.png){width=2.0cm}

## 12.3　特殊函数与角动量代数

`mod_special_functions` 提供 Legendre 多项式、缔合 Legendre 函数、Wigner $3j$、Clebsch–Gordan、$6j$、$9j$ 和转动偶极矩阵元。

### 已有能力

- Legendre 递推；
- Wigner $3j$ 和 Clebsch–Gordan 系数；
- 半整数角动量；
- $6j$ 和 $9j$ 符号；
- 基本选择定则。

### 已知限制

- 没有统一的缓存表；
- 没有批量接口一次性返回整张耦合系数矩阵；
- 没有显式的对称性检验模块；
- 大角动量时的递推稳定性需要更多参数化测试。

### 建议的改进

增加批量接口：

```fortran
pure subroutine cg_table(j1, j2, jmin, jmax, coeffs, stat)
    integer, intent(in) :: j1, j2, jmin, jmax
    real(dp), intent(out) :: coeffs(:, :, :)
    integer, intent(out) :: stat
end subroutine cg_table
```

并增加对称性测试：

$$
C^{JM}_{j_1m_1j_2m_2}
=
(-1)^{j_1-j_2+M}
C^{J,-M}_{j_1,-m_1j_2,-m_2}.
$$

### 对应源码

- [src/mod_special_functions.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_special_functions.f90)
- 二维码：

![](qr/src__mod_special_functions.f90.png){width=2.0cm}

## 12.4　稠密线性代数

`mod_linear_algebra` 提供四类基础能力：

| 能力 | 接口 | 说明 |
|---|---|---|
| 实对称本征求解 | `diag_symmetric_matrix` | Householder 三对角化 + QL 迭代 |
| 实矩阵求逆 | `inv_real_matrix` | Gauss–Jordan 求逆 |
| 复矩阵求逆 | `inv_complex_matrix` | Gauss–Jordan 求逆 |
| FFT | `fft_1d`, `fft_2d` | Cooley–Tukey + 慢 DFT 回退 |

### 已有能力

- 实对称矩阵本征值和本征向量；
- 实矩阵和复矩阵求逆；
- 一维和二维 FFT；
- 非 2 的幂长度的慢 DFT 回退。

### 已知限制

- 没有稀疏本征求解器；
- 没有复 Hermitian 本征求解器；
- 没有 QR、SVD、Cholesky、最小二乘；
- 没有块算法和多线程；
- FFT 只支持特定长度，缺少通用混合基数优化；
- 没有统一的矩阵条件数估计接口。

### 建议的改进

优先级最高的三个扩展是：

1. `diag_hermitian_matrix`；
2. `eig_sparse_symmetric` 或 Lanczos 接口；
3. `svd_real_matrix` 和 `svd_complex_matrix`。

例如复 Hermitian 本征求解可以先转成实对称块矩阵：

$$
H=A+iB,\qquad H=H^{\dagger},
$$

则

$$
A=A^T,\qquad B=-B^T.
$$

可构造

$$
\mathcal H
=
\begin{pmatrix}
A & -B\\
B & A
\end{pmatrix},
$$

其本征值与 $H$ 的本征值一一对应。这样可以在不引入 LAPACK 的情况下复用现有实对称本征求解器。

### 对应源码

- [src/mod_linear_algebra.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_linear_algebra.f90)
- 二维码：

![](qr/src__mod_linear_algebra.f90.png){width=2.0cm}

## 12.5　FFT 与谱计算

FFT 是 Split-Operator 和 FGH 的基础。当前 `fft_1d` 已支持常见 2 的幂长度，并对非 2 的幂使用慢 DFT。

### 已有能力

- 一维 FFT；
- 二维 FFT；
- 逆变换；
- 非 2 的幂回退。

### 已知限制

- 没有原地/非原地统一接口；
- 没有实数专用 FFT；
- 没有批处理接口；
- 没有多线程；
- 没有显式的归一化约定说明；
- 缺少与 NumPy FFT 的大规模对比测试。

### 建议的测试

对随机复数场验证 Parseval 定理：

$$
\sum_n |f_n|^2
=
\frac1N
\sum_k |\tilde f_k|^2.
$$

并验证循环卷积：

$$
(f*g)_n
=
\sum_m f_m g_{n-m}
=
\mathcal F^{-1}
\left[
\tilde f_k\tilde g_k
\right]_n.
$$

## 12.6　DVR 与网格能力

`mod_dvr_grid` 提供一维 Sinc-DVR、Legendre-DVR 和 FGH。

### 已有能力

- `dvr_sinc_init`；
- `dvr_legendre_init`；
- `fgh_solve_bound_states`；
- `dvr_expectation_value`；
- `dvr_matrix_element`。

### 已知限制

- 只有一维 Sinc-DVR 和 Legendre-DVR；
- 没有通用正交多项式工厂；
- 没有多维张量积网格容器；
- 没有非均匀网格；
- 没有自适应网格；
- 没有 DVR 基组与有限基组之间的显式变换接口。

### 建议的改进

引入统一网格接口：

```fortran
type, abstract :: grid_t
    integer :: n_points
contains
    procedure(grid_weights_if), deferred :: weights
    procedure(grid_kinetic_if), deferred :: kinetic
end type grid_t
```

然后派生：

- `sinc_grid_t`
- `legendre_grid_t`
- `hermite_grid_t`
- `laguerre_grid_t`
- `fourier_grid_t`

这样可以在不重写物理模块的情况下扩展更多 DVR。

### 对应源码

- [src/mod_dvr_grid.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_dvr_grid.f90)
- 二维码：

![](qr/src__mod_dvr_grid.f90.png){width=2.0cm}

## 12.7　插值与势能面预处理

`mod_interpolation` 提供三次样条、导数计算和势能面插值。

### 已有能力

- 自然边界三次样条；
- 固定一阶导数边界；
- 一阶导数和二阶导数；
- 势能面插值到网格；
- 长程 $C_6/r^6$ 外推。

### 已知限制

- 只支持一维；
- 没有二维或三维张量样条；
- 没有单调性约束；
- 没有自动外推不确定性估计；
- 没有对剧烈势能变化的局部加密。

### 建议的改进

增加二维双三次样条接口：

```fortran
type :: spline_2d_t
    integer :: nx, ny
    real(dp), allocatable :: x(:), y(:)
    real(dp), allocatable :: coeffs(:, :, :, :)
end type spline_2d_t
```

并增加误差诊断：

$$
\epsilon_{\mathrm{interp}}
=
\max_i
\left|
V_{\mathrm{interp}}(x_i)-V_{\mathrm{ref}}(x_i)
\right|.
$$

### 对应源码

- [src/mod_interpolation.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_interpolation.f90)
- 二维码：

![](qr/src__mod_interpolation.f90.png){width=2.0cm}

## 12.8　数据 I/O

`mod_io_utils` 提供文本表和矩阵输出。

### 已有能力

- `save_data_table_1d`；
- `save_data_table_2d`；
- `save_matrix_dat`；
- 进度条和横幅。

### 已知限制

- 只有文本格式；
- 没有二进制格式；
- 没有 HDF5/NetCDF；
- 没有元数据标准；
- 没有自动单位记录；
- 没有数据版本号；
- 没有异常安全的文件打开/关闭封装。

### 建议的改进

定义统一实验记录结构：

```fortran
type :: experiment_record_t
    character(len=:), allocatable :: title
    character(len=:), allocatable :: commit
    character(len=:), allocatable :: compiler
    character(len=:), allocatable :: units
    real(dp), allocatable :: parameters(:)
end type experiment_record_t
```

并输出 JSON 或 TOML 头，便于 Python 自动读取。

### 对应源码

- [src/mod_io_utils.f90](https://github.com/l1Ha/QuantumGeneralModule/blob/main/src/mod_io_utils.f90)
- 二维码：

![](qr/src__mod_io_utils.f90.png){width=2.0cm}

## 12.9　测试与构建能力

项目当前有自实现测试和三套构建入口。

### 已有能力

- 41 个 Fortran 测试程序；
- 36 个物理示例；
- `Makefile`、`fpm.toml`、`CMakeLists.txt`；
- GitHub Actions CI；
- Python 交叉验证。

### 已知限制

- 没有统一断言库；
- 没有测试覆盖率统计；
- 没有基准测试；
- 没有数值回归阈值管理；
- CI 缺少多维编译器矩阵；
- 缺少自动生成 API 文档。

### 建议的断言接口

```fortran
subroutine expect_close(actual, expected, tol, name, n_pass, n_total)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: name
    integer, intent(inout) :: n_pass, n_total

    n_total = n_total + 1
    if (abs(actual - expected) <= tol) then
        n_pass = n_pass + 1
    else
        print '(A,A,ES12.4,A,ES12.4)', 'FAIL: ', name, actual, expected
    end if
end subroutine expect_close
```

### 对应源码

- [tests/run_all_tests.sh](https://github.com/l1Ha/QuantumGeneralModule/blob/main/tests/run_all_tests.sh)
- 二维码：

![](qr/tests__run_all_tests.sh.png){width=2.0cm}

- [fpm.toml](https://github.com/l1Ha/QuantumGeneralModule/blob/main/fpm.toml)

- 二维码：

![](qr/fpm.toml.png){width=2.0cm}

## 12.10　基础能力缺口汇总

| 领域 | 当前状态 | 优先级 | 建议目标 |
|---|---|---|---|
| 复 Hermitian 本征 | 未具备 | 高 | 支持复哈密顿量 |
| 稀疏本征 | 未具备 | 高 | Lanczos/Arnoldi |
| SVD/QR/Cholesky | 未具备 | 中 | 支持拟合和正则化 |
| 多维插值 | 未具备 | 高 | 二维/三维势能面 |
| 自适应网格 | 未具备 | 中 | 局部加密 |
| 并行化 | 基本未具备 | 高 | OpenMP/MPI |
| 单位类型 | 弱 | 中 | 显式状态码 |
| 数据格式 | 文本为主 | 高 | JSON/TOML/HDF5 |
| 通用测试库 | 自实现 | 中 | 统一断言库 |
| 基准测试 | 缺失 | 中 | 性能回归 |
| 自动微分 | 未具备 | 低 | 支持灵敏度和控制 |

## 12.11　最小扩展路线

1. 增加 `diag_hermitian_matrix`；
2. 增加稀疏对称本征求解；
3. 增加二维样条；
4. 增加 JSON/TOML 实验记录；
5. 增加统一断言库；
6. 增加 OpenMP 版本的 FFT 和矩阵乘法；
7. 增加基准脚本；
8. 为每个基础模块增加收敛性测试。

## 12.12　参考文献

1. G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, Baltimore, 2013.
2. Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, Philadelphia, 2003.
3. C. R. Harris et al., “Array programming with NumPy”, *Nature* **585**, 357 (2020). DOI: 10.1038/s41586-020-2649-2.
4. GeneralModule 源码与文档：<https://github.com/l1Ha/QuantumGeneralModule>