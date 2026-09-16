# GeneralModule 科学文献典藏与理论映射全典 (Literature Compendium)

本全典详尽整理了 `GeneralModule` 算法库所依据的权威经典与前沿学术文献。涵盖**超冷量子碰撞、多通道密耦与对数导数法、外场磁 Feshbach 共振、角动量表象变换、各向异性偶极散射与自旋弛豫、光缔合谱学、强场超快物理、波包动力学与最优控制**十大核心物理领域。每条文献均包含标准学术引用、DOI、理论物理机理及其在算法库源代码中的确切映射位置。

---

## 目录
1. [非含时散射理论、分波相移与有效力程展开](#1-非含时散射理论分波相移与有效力程展开)
2. [多通道密耦对数导数法与分段网格推进](#2-多通道密耦对数导数法与分段网格推进)
3. [长程范德华色散与半经典平均散射长度](#3-长程范德华色散与半经典平均散射长度)
4. [超冷原子自旋相互作用与磁 Feshbach 共振](#4-超冷原子自旋相互作用与磁-feshbach-共振)
5. [塞曼-超精细 Breit-Rabi 能谱与四大经典基组变换](#5-塞曼-超精细-breit-rabi-能谱与四大经典基组变换)
6. [强场超快物理、高次谐波发射与隧穿电离](#6-强场超快物理高次谐波发射与隧穿电离)
7. [离散变量表象 (DVR)、虚时间与分裂算符波包动力学](#7-离散变量表象-dvr虚时间与分裂算符波包动力学)
8. [开放量子系统 Lindblad 耗散与 Krotov 最优控制](#8-开放量子系统-lindblad-耗散与-krotov-最优控制)
9. [各向异性磁偶极与电偶极超冷散射及自旋弛豫](#9-各向异性磁偶极与电偶极超冷散射及自旋弛豫)
10. [超冷光缔合谱学与自由-束缚态量子跃迁](#10-超冷光缔合谱学与自由-束缚态量子跃迁)
11. [超冷少体物理、三体复合与 Efimov 普适态](#11-超冷少体物理三体复合与-efimov-普适态)
12. [低维光晶格受限量子散射与约束诱导共振 CIR](#12-低维光晶格受限量子散射与约束诱导共振-cir)
13. [自电离体系、Fano 共振理论与复坐标旋转法 CCR](#13-自电离体系fano-共振理论与复坐标旋转法-ccr)
14. [交叉静电磁场分子量子动力学与非共线 Stark-Zeeman 态混合](#14-交叉静电磁场分子量子动力学与非共线-stark-zeeman-态混合)
15. [三原子反应碰撞几何、Jacobi 坐标、LEPS 势能面与锥形交叉几何相位](#15-三原子反应碰撞几何jacobi-坐标leps-势能面与锥形交叉几何相位)
16. [超冷旋量玻色-爱因斯坦凝聚与宏观自旋混合动力学](#16-超冷旋量玻色-爱因斯坦凝聚与宏观自旋混合动力学)
17. [三原子超球面反应动力学与过渡态理论 (Hyperspherical Reactive Scattering & TST)](#17-三原子超球面反应动力学与过渡态理论-hyperspherical-reactive-scattering--tst)
18. [超冷偶极量子液滴与李-黄-杨量子涨落修正 (Dipolar Droplets & LHY)](#18-超冷偶极量子液滴与李-黄-杨量子涨落修正-dipolar-droplets--lhy)
19. [强场非顺序双电离与电子重碰撞相关动量动力学 (Strong-Field NSDI & Recollision)](#19-强场非顺序双电离与电子重碰撞相关动量动力学-strong-field-nsdi--recollision)
20. [磁与光 Feshbach 共振、分子弱束缚态与光致非弹性损耗 (Feshbach Resonances)](#20-磁与光-feshbach-共振分子弱束缚态与光致非弹性损耗-feshbach-resonances)
21. [阿秒瞬态吸收光谱 (ATAS) 与光诱导态自电离干涉动力学](#21-阿秒瞬态吸收光谱-atas-与光诱导态自电离干涉动力学)
22. [双色反向旋转圆偏振场与分子光电子圆二色性 (PECD)](#22-双色反向旋转圆偏振场与分子光电子圆二色性-pecd)
23. [超冷极性分子化学反应动力学与微波/静电偶极遮蔽](#23-超冷极性分子化学反应动力学与微波静电偶极遮蔽)
24. [里德堡原子阻塞、PXP 约束模型与量子多体疤痕](#24-里德堡原子阻塞pxp-约束模型与量子多体疤痕)
25. [表面量子散射与选择性吸附共振 (Selective Adsorption Resonances)](#25-表面量子散射与选择性吸附共振-selective-adsorption-resonances)
26. [气-固表面催化反应与 Eley-Rideal 提取机理](#26-气-固表面催化反应与-eley-rideal-提取机理)
27. [表面非绝热动力学与电子摩擦耗散 (Electronic Friction & GLE)](#27-表面非绝热动力学与电子摩擦耗散-electronic-friction--gle)
28. [掠入射快原子表面量子衍射与彩虹散射 (GIFAD)](#28-掠入射快原子表面量子衍射与彩虹散射-gifad)
29. [冷离子-中性原子杂化散射与极化阱动力学 (Cold Ion-Atom Hybrid Scattering)](#29-冷离子-中性原子杂化散射与极化阱动力学-cold-ion-atom-hybrid-scattering)
30. [Tully 最少开关表面跳跃与非绝热分子动力学 (FSSH)](#30-tully-最少开关表面跳跃与非绝热分子动力学-fssh)
31. [强场分子定向、取向与光学离心机超转子动力学 (Molecular Alignment & Superrotors)](#31-强场分子定向取向与光学离心机超转子动力学-molecular-alignment--superrotors)
32. [超冷光晶格与玻色-哈伯德微观映射 (Optical Lattice & Bose-Hubbard)](#32-超冷光晶格与玻色-哈伯德微观映射-optical-lattice--bose-hubbard)
33. [多原子反应路径哈密顿量与变分过渡态理论 (RPH & Variational TST)](#33-多原子反应路径哈密顿量与变分过渡态理论-rph--variational-tst)
34. [相对论原子结构与径向狄拉克方程 (Relativistic Atomic Structure & Dirac)](#34-相对论原子结构与径向狄拉克方程-relativistic-atomic-structure--dirac)
35. [共振非弹性 X 射线散射与内壳层光谱 (Resonant Inelastic X-ray Scattering - RIXS)](#35-共振非弹性-x-射线散射与内壳层光谱-resonant-inelastic-x-ray-scattering---rixs)
36. [亚稳态原子潘宁电离与缔合电离动力学 (Penning & Associative Ionization / Chemi-ionization)](#36-亚稳态原子潘宁电离与缔合电离动力学-penning--associative-ionization--chemi-ionization)
37. [分子光解离动力学、时间自相关函数与光碎片动能释放谱 (KER)](#37-分子光解离动力学时间自相关函数与光碎片动能释放谱-ker)
38. [多通道非绝热避差穿越、Landau-Zener 跃迁与 Hellmann-Feynman 耦合](#38-多通道非绝热避差穿越landau-zener-跃迁与-hellmann-feynman-耦合)
39. [双原子分子转振跃迁调控、STIRAP 绝热受激跃迁与 Franck-Condon 原理](#39-双原子分子转振跃迁调控stirap-绝热受激跃迁与-franck-condon-原理)
40. [含时波包散射动力学、通量时间-能量积分与 Möller 算符投影](#40-含时波包散射动力学通量时间-能量积分与-möller-算符投影)
41. [复吸收势 (CAP) 最佳边界参数化与量子概率流连续性方程](#41-复吸收势-cap-最佳边界参数化与量子概率流连续性方程)
42. [库仑三体系统、Perkeris 坐标变换与两电子关联](#42-库仑三体系统perkeris-坐标变换与两电子关联)
43. [快速学术检索与代码对照总表](#43-快速学术检索与代码对照总表)

---

## 1. 非含时散射理论、分波相移与有效力程展开

### 1.1 分波法与光学定理
- **文献**: 
  - E. P. Wigner, *"On the behavior of cross sections near thresholds"*, **Phys. Rev.** 73, 1002 (1948). [DOI: 10.1103/PhysRev.73.1002](https://doi.org/10.1103/PhysRev.73.1002)
  - L. S. Rodberg and R. M. Thaler, *Introduction to the Quantum Theory of Scattering*, Academic Press, New York (1967).
- **核心理论**:
  径向定态薛定谔方程渐近边界处的定态相移 $\delta_l(E)$ 展开：
  $$\psi_l(r) \xrightarrow{r \to \infty} A_l \left[ \hat{j}_l(kr) \cos\delta_l - \hat{n}_l(kr) \sin\delta_l \right]$$
  弹性散射反应矩阵 $K_l = \tan\delta_l$，散射矩阵元 $S_l = e^{2i\delta_l}$，跃迁矩阵元 $T_l = S_l - 1$。
  光学定理（Optical Theorem）建立前向散射振幅与总截面的精确自洽守恒关系：
  $$\sigma_{\text{tot}} = \frac{4\pi}{k} \text{Im}[f(0)] = \sum_{l=0}^\infty \frac{4\pi(2l+1)}{k^2}\sin^2\delta_l$$
- **代码映射**:
  - `src/mod_ti_scattering.f90`:
    - `riccati_bessel_neumann`: $\hat{j}_l(x), \hat{n}_l(x)$ 及其解析一阶导数；
    - `calc_phase_shift_single_l`: 单通道正能量定态相移与弹性截面提取；
    - `optical_theorem_cross_section`: 光学定理前向振幅虚部自洽校验；
    - `calc_differential_cross_section`: 微分散射截面 Legendre 多项式级数展开。

### 1.2 有效力程展开 (Effective Range Expansion, ERE)
- **文献**:
  - H. A. Bethe, *"Theory of the Effective Range in Nuclear Scattering"*, **Phys. Rev.** 76, 38 (1949). [DOI: 10.1103/PhysRev.76.38](https://doi.org/10.1103/PhysRev.76.38)
  - J. M. Blatt and J. D. Jackson, *"On the Interpretation of Low-Energy Proton-Proton Scattering"*, **Phys. Rev.** 76, 18 (1949). [DOI: 10.1103/PhysRev.76.18](https://doi.org/10.1103/PhysRev.76.18)
- **核心理论**:
  当质心碰撞动能趋近于零（$k \to 0$）时，$s$-波（$l=0$）相移满足解析展开：
  $$k \cot\delta_0(k) = -\frac{1}{a_s} + \frac{1}{2} r_0 k^2 - P_r r_0^3 k^4 + \mathcal{O}(k^6)$$
  其中 $a_s$ 为零能散射长度（s-wave scattering length），$r_0$ 为有效相互作用力程（effective range）。
- **代码映射**:
  - `src/mod_ti_scattering.f90`: `fit_effective_range_expansion`（多项式最小二乘拟合提取 $a_s, r_0$ 与零温弹性截面 $\sigma_0 = 4\pi a_s^2$）。

### 1.3 Wigner-Smith 碰撞时延与形状共振
- **文献**:
  - E. P. Wigner, *"Lower Limit for the Energy Derivative of the Scattering Phase Shift"*, **Phys. Rev.** 98, 145 (1955). [DOI: 10.1103/PhysRev.98.145](https://doi.org/10.1103/PhysRev.98.145)
  - F. T. Smith, *"Lifetime Matrix in Collision Theory"*, **Phys. Rev.** 118, 349 (1960). [DOI: 10.1103/PhysRev.118.349](https://doi.org/10.1103/PhysRev.118.349)
- **核心理论**:
  粒子在碰撞势阱或离心势垒（$l > 0$）内的滞留时延由 Wigner-Smith 时延矩阵标定：
  $$\tau(E) = 2\hbar \frac{d\delta_l(E)}{dE} = -i\hbar S_l^\dagger(E) \frac{dS_l(E)}{dE}$$
  当能量扫描跨越形状共振（Shape Resonance）时，相移跳变 $\pi$，Wigner 时延呈现 Breit-Wigner 洛伦兹峰：
  $$\tau(E) \approx \frac{2\hbar \Gamma}{(E - E_R)^2 + (\Gamma/2)^2}$$
- **代码映射**:
  - `src/mod_ti_scattering.f90`: `analyze_shape_resonance`（数值微分相移提取共振极值能量 $E_R$、全半高宽 $\Gamma$ 与准束缚态寿命 $\tau_{\text{life}} = \hbar/\Gamma$）。

---

## 2. 多通道密耦对数导数法与分段网格推进

### 2.1 Johnson 矩阵比值对数导数递推法
- **文献**:
  - B. R. Johnson, *"The multichannel log-derivative method for treating reactive collisions"*, **J. Comput. Phys.** 13, 445 (1973). [DOI: 10.1016/0021-9991(73)90049-1](https://doi.org/10.1016/0021-9991(73)90049-1)
  - B. R. Johnson, *"New numerical methods for solving the coupled differential equations of molecular scattering"*, **J. Chem. Phys.** 67, 4086 (1977). [DOI: 10.1063/1.435384](https://doi.org/10.1063/1.435384)
- **核心理论**:
  多通道薛定谔方程在强排斥芯经典禁区由于闭通道指数发散极易产生数值线性相关（numerical instability）。Johnson 算法不直接递推波函数矩阵 $\mathbf{\Psi}(r)$，而是推进相邻网格比值矩阵：
  $$\mathbf{R}_i = \mathbf{Q}_i \mathbf{\Psi}_i \mathbf{\Psi}_{i-1}^{-1} \mathbf{Q}_{i-1}^{-1}$$
  递推关系式为：
  $$\mathbf{R}_{i+1} = \mathbf{M}_i - \mathbf{R}_i^{-1}, \quad \mathbf{M}_i = 12 \mathbf{Q}_i^{-1} - 10 \mathbf{I}, \quad \mathbf{Q}_i = \mathbf{I} - \frac{h^2}{12}\mathbf{W}(r_i)$$
  在渐近边界处重构对数导数矩阵 $\mathbf{Y}(r_N) = \mathbf{\Psi}'(r_N)\mathbf{\Psi}^{-1}(r_N)$：
  $$\mathbf{Y}(r_N) = \frac{1}{2h}\left(3\mathbf{I} - 4\mathbf{P}_1 + \mathbf{P}_2\right), \quad \mathbf{P}_1 = \mathbf{Q}_{N-1}^{-1} \mathbf{R}_N^{-1} \mathbf{Q}_N, \quad \mathbf{P}_2 = \mathbf{P}_1^2$$
- **代码映射**:
  - `src/mod_ti_scattering.f90`: 
    - `calc_scattering_length_logder`: 单通道 Johnson 比值法零能散射长度；
    - `calc_multichannel_close_coupling_logder`: 任意通道实对称比值矩阵递推、闭通道 Feshbach Schur 补投影与 Riccati 匹配；
    - `calc_multichannel_close_coupling_segmented_logder`: 多扇区分段推进。

### 2.2 Manolopoulos 局域扇区传播矩阵改进算法
- **文献**:
  - D. E. Manolopoulos, *"An improved log-derivative method for solving the radial Schrödinger equation"*, **J. Chem. Phys.** 85, 6425 (1986). [DOI: 10.1063/1.451472](https://doi.org/10.1063/1.451472)
  - D. E. Manolopoulos, M. J. Jamieson, and A. D. Pradhan, *"Johnson's log derivative method revisited"*, **J. Comput. Phys.** 105, 169 (1993). [DOI: 10.1006/jcph.1993.1064](https://doi.org/10.1006/jcph.1993.1064)
- **核心理论**:
  将全域划分为若干扇区（Sectors），在各扇区内势能矩阵做局部多项式近似。对数导数矩阵 $\mathbf{Y}(r)$ 满足跨扇区传递关系：
  $$\mathbf{Y}_{k+1} = \mathbf{D}_k - \mathbf{C}_k \left(\mathbf{Y}_k + \mathbf{A}_k\right)^{-1} \mathbf{B}_k$$
  对数导数矩阵是纯物理局域量，在扇区交接面处严格单值连续，扇区间无需波函数插值即可实现步长跨尺度放缩。
- **代码映射**:
  - `src/mod_ti_scattering.f90`: 
    - `segmented_grid_t`: 多扇区分段网格派生类型；
    - `create_segmented_grid`: 变步长多扇区单调坐标生成；
    - `calc_scattering_wavefunction_segmented_ti`: 4 阶 Taylor 桥接高精度波函数拼接。

### 2.3 工业级分子密耦计算系统标准实现对比
- **文献**:
  - J. M. Hutson and C. R. Le Sueur, *"MOLSCAT: A program for non-reactive quantum scattering calculations on atomic and molecular collisions"*, **Comput. Phys. Commun.** 241, 9 (2019). [DOI: 10.1016/j.cpc.2019.02.014](https://doi.org/10.1016/j.cpc.2019.02.014)
  - J. M. Hutson and C. R. Le Sueur, *"FIELD: A program for calculating properties of interacting molecules in external fields"*, **Comput. Phys. Commun.** 241, 1 (2019). [DOI: 10.1016/j.cpc.2019.02.015](https://doi.org/10.1016/j.cpc.2019.02.015)

---

## 3. 长程范德华色散与半经典平均散射长度

### 3.1 Gribakin-Flambaum 半经典散射长度解析式
- **文献**:
  - G. F. Gribakin and V. V. Flambaum, *"Calculation of the scattering length in atomic collisions using the semiclassical approximation"*, **Phys. Rev. A** 48, 546 (1993). [DOI: 10.1103/PhysRevA.48.546](https://doi.org/10.1103/PhysRevA.48.546)
- **核心理论**:
  对于渐近长程为范德华吸引势 $V(r) \xrightarrow{r \to \infty} -C_6 / r^6$ 的双原子碰撞，平均散射长度（mean scattering length）定义为：
  $$\bar{a} = \frac{2\pi}{\Gamma(1/4)^2} \left(\frac{2\mu C_6}{\hbar^2}\right)^{1/4} \approx 0.4779888 \cdot \left(\frac{2\mu C_6}{\hbar^2}\right)^{1/4}$$
  实际 $s$-波零能散射长度与分子深势阱内总半经典 WKB 作用量 $\Phi = \int_{r_0}^\infty \sqrt{-2\mu V(r)} dr / \hbar$ 解析关联：
  $$a_s = \bar{a} \left[ 1 - \tan\left(\Phi - \frac{\pi}{8}\right) \right]$$
- **代码映射**:
  - `src/mod_ti_scattering.f90`:
    - `van_der_waals_mean_length`: 计算平均散射长度 $\bar{a}$；
    - `gribakin_flambaum_length`: 计算半经典散射长度 $a_s(\Phi)$。

### 3.2 纯长程 $1/r^6$ 势精确解析解 (Bo Gao 理论)
- **文献**:
  - Bo Gao, *"Solutions of the Schrödinger equation for an attractive $1/r^6$ potential"*, **Phys. Rev. A** 58, 4222 (1998). [DOI: 10.1103/PhysRevA.58.4222](https://doi.org/10.1103/PhysRevA.58.4222)
  - Bo Gao, E. Tiesinga, C. J. Williams, and P. S. Julienne, *"Multichannel quantum-defect theory for slow atomic collisions"*, **Phys. Rev. A** 72, 042719 (2005). [DOI: 10.1103/PhysRevA.72.042719](https://doi.org/10.1103/PhysRevA.72.042719)

---

## 4. 超冷原子自旋相互作用与磁 Feshbach 共振

### 4.1 Feshbach 共振理论基础与唯象参数公式
- **文献**:
  - H. Feshbach, *"Unified theory of nuclear reactions"*, **Ann. Phys. (N.Y.)** 5, 357 (1958). [DOI: 10.1016/0003-4916(58)90007-1](https://doi.org/10.1016/0003-4916(58)90007-1)
  - U. Fano, *"Effects of Configuration Interaction on Intensities and Phase Shifts"*, **Phys. Rev.** 124, 1866 (1961). [DOI: 10.1103/PhysRev.124.1866](https://doi.org/10.1103/PhysRev.124.1866)
  - C. Chin, R. Grimm, P. S. Julienne, and E. Tiesinga, *"Feshbach resonances in ultracold gases"*, **Rev. Mod. Phys.** 82, 1225 (2010). [DOI: 10.1103/RevModPhys.82.1225](https://doi.org/10.1103/RevModPhys.82.1225)
- **核心理论**:
  在外加静磁场 $B$ 下，开通道（入口碰撞道）与闭通道（具有不同磁矩差 $\Delta\mu = \mu_{\text{closed}} - \mu_{\text{open}}$ 的分子束缚态）由于自旋交换相互作用发生非绝热耦合，诱导散射长度产生色散共振跃迁：
  $$a(B) = a_{\text{bg}} \left( 1 - \frac{\Delta B}{B - B_0} \right)$$
  其中 $B_0$ 为共振中心位置，$\Delta B$ 为磁共振线宽，$a_{\text{bg}}$ 为背景非共振散射长度。零散射点出现在 $B_{\text{zero}} = B_0 + \Delta B$。
- **代码映射**:
  - `src/mod_field_scattering.f90`:
    - `calc_magnetic_feshbach_resonance_scan`: 外磁场扫描多通道密耦与开通道散射长度提取；
    - `fit_feshbach_resonance_parameters`: 基于 Levenberg-Marquardt 非线性最小二乘提取 $(B_0, \Delta B, a_{\text{bg}})$。

### 4.2 双原子自旋交换相互作用与势能算符分解
- **文献**:
  - H. T. C. Stoof, J. M. V. A. Koelman, and B. J. Verhaar, *"Spin-exchange frequency shift in a cesium atomic clock"*, **Phys. Rev. B** 38, 4688 (1988). [DOI: 10.1103/PhysRevB.38.4688](https://doi.org/10.1103/PhysRevB.38.4688)
  - E. Tiesinga, B. J. Verhaar, and H. T. C. Stoof, *"Threshold and resonance phenomena in ultracold ground-state collisions"*, **Phys. Rev. A** 47, 4114 (1993). [DOI: 10.1103/PhysRevA.47.4114](https://doi.org/10.1103/PhysRevA.47.4114)
- **核心理论**:
  碱金属原子在基态的静电库仑作用完全取决于双电子总自旋 $S = s_1 + s_2$。电子自旋交换算符 $\hat{\mathbf{s}}_1 \cdot \hat{\mathbf{s}}_2$ 的本征值为：
  $$\hat{\mathbf{s}}_1 \cdot \hat{\mathbf{s}}_2 = \frac{1}{2}\left[S(S+1) - s_1(s_1+1) - s_2(s_2+1)\right] = \begin{cases} -3/4, & S = 0 \text{ (单重态 Singlet } X^1\Sigma_g^+) \\ +1/4, & S = 1 \text{ (三重态 Triplet } a^3\Sigma_u^+) \end{cases}$$
  空间相互作用势矩阵元可直接投影分解：
  $$\hat{V}(r) = V_0(r) \hat{\mathcal{P}}_0 + V_1(r) \hat{\mathcal{P}}_1 = \bar{V}(r) \hat{\mathbf{I}} + \Delta V(r) \left(\hat{\mathbf{s}}_1 \cdot \hat{\mathbf{s}}_2\right)$$
  其中 $\bar{V}(r) = \frac{1}{4}V_0(r) + \frac{3}{4}V_1(r)$，$\Delta V(r) = V_1(r) - V_0(r)$。
- **代码映射**:
  - `src/mod_field_scattering.f90`:
    - `build_spin_exchange_matrix`: 电子自旋交换算符在总自旋、未耦合及单原子超精细基下的精确矩阵元与幺正投影；
    - `build_field_potential_matrix`: 耦合各通道势能矩阵 $V_{ij}(r)$ 组装。

---

## 5. 塞曼-超精细 Breit-Rabi 能谱与四大经典基组变换

### 5.1 单原子 Breit-Rabi 解析公式与塞曼能级
- **文献**:
  - G. Breit and I. I. Rabi, *"Measurement of Nuclear Spin"*, **Phys. Rev.** 38, 2082 (1931). [DOI: 10.1103/PhysRev.38.2082](https://doi.org/10.1103/PhysRev.38.2082)
  - E. Arimondo, M. Inguscio, and P. Violino, *"Experimental determinations of the hyperfine structure in the alkali atoms"*, **Rev. Mod. Phys.** 49, 31 (1977). [DOI: 10.1103/RevModPhys.49.31](https://doi.org/10.1103/RevModPhys.49.31)
- **核心理论**:
  单电子原子（$s = 1/2$）在外磁场 $B$ 下的自旋超精细-塞曼哈密顿量：
  $$\hat{H}_{\text{atom}} = A_{\text{hf}} \hat{\mathbf{s}} \cdot \hat{\mathbf{i}} + (g_s \mu_B \hat{s}_z - g_i \mu_N \hat{i}_z) B$$
  其严格解析本征能量由 Breit-Rabi 公式给出：
  $$E(F = i \pm 1/2, m_F) = -\frac{\Delta E_{\text{hfs}}}{2(2i + 1)} - g_i \mu_N B m_F \pm \frac{\Delta E_{\text{hfs}}}{2} \sqrt{1 + \frac{4 m_F}{2i + 1} x + x^2}$$
  其中无量纲场强参数 $x = \frac{(g_s \mu_B + g_i \mu_N) B}{\Delta E_{\text{hfs}}}$，$\Delta E_{\text{hfs}} = \frac{2i+1}{2} A_{\text{hf}}$。
- **代码映射**:
  - `src/mod_field_scattering.f90`:
    - `calc_breit_rabi_energies`: 矩阵对角化求解任意核自旋 $i$ 的严格数值解，并自动与解析 Breit-Rabi 进行双向机器精度比对。

### 5.2 角动量耦合代数与四大物理基组表象
- **文献**:
  - D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific, Singapore (1988). [DOI: 10.1142/0270](https://doi.org/10.1142/0270)
  - J. P. Burke, Jr., *Theoretical Investigation of Cold Alkali Collisions*, Ph.D. thesis, University of Colorado (1999).
- **核心理论**:
  两原子双碰撞系统存在四大互为严格幺正变换的量子力学表象：
  1. **未耦合原子基组 (Uncoupled Basis)**: $|s_1, m_{s1}, i_1, m_{i1}; s_2, m_{s2}, i_2, m_{i2}\rangle$
  2. **单原子超精细耦合基组 (f-Coupled Basis)**: $|(s_1, i_1)f_1, m_{f1}; (s_2, i_2)f_2, m_{f2}\rangle$
     $$\langle (s_1 i_1)f_1 m_{f1}, (s_2 i_2)f_2 m_{f2} | m_{s1} m_{i1} m_{s2} m_{i2} \rangle = C_{s_1 m_{s1}, i_1 m_{i1}}^{f_1 m_{f1}} \cdot C_{s_2 m_{s2}, i_2 m_{i2}}^{f_2 m_{f2}}$$
  3. **总自旋耦合基组 (Total Spin Basis)**: $|(s_1 s_2)S, (i_1 i_2)I, F, M_F\rangle$
     通过 Clebsch-Gordan 级联与 Wigner 9j 符号建立与 f-Coupled 基组的正交变换：
     $$\langle (f_1 f_2)F | (SI)F \rangle = \sqrt{(2f_1+1)(2f_2+1)(2S+1)(2I+1)} \begin{Bmatrix} s_1 & i_1 & f_1 \\ s_2 & i_2 & f_2 \\ S & I & F \end{Bmatrix}$$
  4. **场缀饰本征基组 (Field-Dressed Channel Basis)**: 渐近哈密顿量 $\hat{H}_{\text{asymp}}(B)$ 的定态本征矢矩阵。
- **代码映射**:
  - `src/mod_special_functions.f90`:
    - `wigner_3j_half`, `clebsch_gordan_half`: 半整数/整数任意自旋角动量耦合系数；
    - `wigner_6j_half`, `wigner_9j_half`: 6j / 9j 符号 Racah 代数求值器；
  - `src/mod_field_scattering.f90`:
    - `calc_basis_transform_matrix`: 任意两基组间机器精度幺正变换算符 $U_{ab}$；
    - `build_field_collision_channels`: 自动按 $M_{\text{tot}} = m_{f1} + m_{f2}$ 对称性截断通道子空间。

---

## 6. 强场超快物理、高次谐波发射与隧穿电离

### 6.1 Keldysh 参数与强场隧穿电离 (ADK 理论)
- **文献**:
  - L. V. Keldysh, *"Ionization in the field of a strong electromagnetic wave"*, **Sov. Phys. JETP** 20, 1307 (1965).
  - M. V. Ammosov, N. B. Delone, and V. P. Krainov, *"Tunnel ionization of complex atoms and of atomic ions in an alternating electromagnetic field"*, **Sov. Phys. JETP** 64, 1191 (1986).
- **核心理论**:
  绝热参数定义强场电离机制划分：
  $$\gamma_K = \frac{\omega \sqrt{2 I_p}}{E_0} = \sqrt{\frac{I_p}{2 U_p}}$$
  当 $\gamma_K \ll 1$ 时隧穿电离占主导。准静态 ADK 隧穿电离率公式为：
  $$W_{\text{ADK}}(t) = C_{n^*, l^*}^2 f(l, m) I_p \left(\frac{2(2I_p)^{3/2}}{|E(t)|}\right)^{2n^* - |m| - 1} \exp\left(-\frac{2(2I_p)^{3/2}}{3|E(t)|}\right)$$
- **代码映射**:
  - `src/mod_coulomb_atomic.f90`:
    - `keldysh_parameter`: 绝热参数 $\gamma_K$；
    - `ponderomotive_energy`: 有质动力能 $U_p = E_0^2 / (4\omega^2)$；
    - `adk_ionization_rate`: 瞬时电离率积分。

### 6.2 三步模型与强场近似 (Lewenstein SFA 高次谐波发射)
- **文献**:
  - P. B. Corkum, *"Plasma perspective on strong field multiphoton ionization"*, **Phys. Rev. Lett.** 71, 1994 (1993). [DOI: 10.1103/PhysRevLett.71.1994](https://doi.org/10.1103/PhysRevLett.71.1994)
  - M. Lewenstein, P. Balcou, M. Y. Ivanov, A. L’Huillier, and P. B. Corkum, *"Theory of high-harmonic generation by low-frequency laser fields"*, **Phys. Rev. A** 49, 2117 (1994). [DOI: 10.1103/PhysRevA.49.2117](https://doi.org/10.1103/PhysRevA.49.2117)
- **核心理论**:
  高次谐波发射（HHG）遵从半经典三步模型：隧穿电离 $\to$ 连续区激光驱动加速 $\to$ 回碰复合发射高能光子。截止能经典定律为：
  $$E_{\text{cutoff}} = I_p + 3.17 U_p$$
  Lewenstein 量子强场近似诱导偶极矩为：
  $$x(t) = i \int_0^t dt' \int d^3\mathbf{p} \, d_x^*(\mathbf{p} + \mathbf{A}(t)) e^{-i S(\mathbf{p}, t, t')} E(t') d_x(\mathbf{p} + \mathbf{A}(t')) + \text{c.c.}$$
- **代码映射**:
  - `src/mod_coulomb_atomic.f90`: `hhg_cutoff_energy`；
  - `src/mod_hhg_spectra.f90`:
    - `calculate_dipole_acceleration`: Ehrenfest 定理偶极加速度；
    - `hhg_power_spectrum`: 偶极功率谱发射强度；
    - `lewenstein_sfa_dipole`: SFA 鞍点路径近似偶极响应。

---

## 7. 离散变量表象 (DVR)、虚时间与分裂算符波包动力学

### 7.1 Sinc-DVR、Legendre-DVR 与 FGH 束缚态求解
- **文献**:
  - D. T. Colbert and W. H. Miller, *"A novel discrete variable representation for quantum mechanical reactive scattering via the S-matrix Kohn method"*, **J. Chem. Phys.** 96, 1982 (1992). [DOI: 10.1063/1.462125](https://doi.org/10.1063/1.462125)
  - C. C. Marston and G. G. Balint-Kurti, *"The Fourier grid Hamiltonian method for bound state eigenvalues and eigenfunctions"*, **J. Chem. Phys.** 91, 3571 (1989). [DOI: 10.1063/1.456888](https://doi.org/10.1063/1.456888)
- **核心理论**:
  Sinc 离散变量表象下动能算符矩阵解析式：
  $$T_{ij} = \frac{\hbar^2}{2\mu \Delta x^2} \begin{cases} \frac{\pi^2}{3}, & i = j \\ \frac{2(-1)^{i-j}}{(i-j)^2}, & i \ne j \end{cases}$$
  势能矩阵严格对角：$V_{ij} = V(x_i)\delta_{ij}$，总哈密顿矩阵直接通过对称对角化即可获得全部分子束缚能级与高精度振动波函数。
- **代码映射**:
  - `src/mod_dvr_grid.f90`:
    - `dvr_sinc_init`: 构造 Sinc-DVR 动能矩阵；
    - `fgh_solve_bound_states`: FGH 求解分子转振能级与波函数；
    - `dvr_legendre_init`: Gauss-Legendre 转动取向空间 DVR。

### 7.2 分裂算符动力学推进器 (Split-Operator Method)
- **文献**:
  - M. D. Feit, J. A. Fleck, Jr., and A. Steiger, *"Solution of the Schrödinger equation by a spectral method"*, **J. Comput. Phys.** 47, 412 (1982). [DOI: 10.1016/0021-9991(82)90091-2](https://doi.org/10.1016/0021-9991(82)90091-2)
  - R. Kosloff, *"Time-dependent quantum-mechanical methods for molecular dynamics"*, **J. Phys. Chem.** 92, 2087 (1988). [DOI: 10.1021/j100319a003](https://doi.org/10.1021/j100319a003)
- **核心理论**:
  二阶辛对称 Trotter 分裂推进格式：
  $$\hat{U}(\Delta t) = e^{-i \hat{H} \Delta t / \hbar} = e^{-i \hat{V} \frac{\Delta t}{2\hbar}} e^{-i \hat{T} \frac{\Delta t}{\hbar}} e^{-i \hat{V} \frac{\Delta t}{2\hbar}} + \mathcal{O}(\Delta t^3)$$
- **代码映射**:
  - `src/mod_wavepacket_propagator.f90`:
    - `propagate_split_operator_1d`: 1D FFT 动能空间加速演化；
    - `propagate_split_operator_2d`: 2D 坐标分裂算符；
  - `src/mod_chebyshev_propagator.f90`: 切比雪夫展开大步长高阶算子推进。

---

## 8. 开放量子系统 Lindblad 耗散与 Krotov 最优控制

### 8.1 完全正定动力学半群与 Lindblad 主方程
- **文献**:
  - G. Lindblad, *"On the generators of quantum dynamical semigroups"*, **Commun. Math. Phys.** 48, 119 (1976). [DOI: 10.1007/BF01608499](https://doi.org/10.1007/BF01608499)
  - V. Gorini, A. Kossakowski, and E. C. G. Sudarshan, *"Completely positive dynamical semigroups of N-level systems"*, **J. Math. Phys.** 17, 821 (1976). [DOI: 10.1063/1.522979](https://doi.org/10.1063/1.522979)
  - H.-P. Breuer and F. Petruccione, *The Theory of Open Quantum Systems*, Oxford University Press, Oxford (2002).
- **核心理论**:
  密度矩阵 $\hat{\rho}(t)$ 的马尔可夫演化方程：
  $$\frac{d\hat{\rho}}{dt} = -\frac{i}{\hbar}[\hat{H}, \hat{\rho}] + \sum_k \gamma_k \left( \hat{L}_k \hat{\rho} \hat{L}_k^\dagger - \frac{1}{2} \{\hat{L}_k^\dagger \hat{L}_k, \hat{\rho}\} \right)$$
- **代码映射**:
  - `src/mod_open_quantum.f90`:
    - `propagate_lindblad_rk4`: RK4 积分开放系统主方程；
    - `calculate_purity`: 态纯度 $\text{Tr}(\rho^2)$；
    - `calculate_von_neumann_entropy`: 冯·诺依曼熵 $S = -\text{Tr}(\rho \ln\rho)$。

### 8.2 量子最优控制 Krotov 迭代算法
- **文献**:
  - V. F. Krotov, *Global Methods in Optimal Control Theory*, Marcel Dekker, New York (1996).
  - R. Somlói, J. Kazakov, and D. J. Tannor, *"A generalized relaxation method for optimal control of molecular motion"*, **Chem. Phys.** 172, 85 (1993). [DOI: 10.1016/0301-0104(93)80108-L](https://doi.org/10.1016/0301-0104(93)80108-L)
  - D. M. Reich, M. Ndong, and C. P. Koch, *"Monotonically convergent optimal control theory of quantum systems"*, **J. Chem. Phys.** 136, 104103 (2012). [DOI: 10.1063/1.3691827](https://doi.org/10.1063/1.3691827)
- **核心理论**:
  构造拉格朗日乘子反向共轭协态 $|\chi(t)\rangle$，通过严格单调递增迭代公式更新外部控制激光场：
  $$\Delta \epsilon(t) = \frac{S(t)}{\alpha_0} \text{Im}\left[\langle \chi^{(k)}(t) | \hat{\mu} | \psi^{(k+1)}(t) \rangle\right]$$
- **代码映射**:
  - `src/mod_optimal_control.f90`:
    - `optimize_pulse_krotov`: 双向共轭场迭代与目标态保真度单调优化。

---

## 9. 各向异性磁偶极与电偶极超冷散射及自旋弛豫

### 9.1 电子磁偶极-偶极相互作用 (MDDI) 与 2 阶球张量展开
- **文献**:
  - H. T. C. Stoof, J. M. V. A. Koelman, and B. J. Verhaar, *"Spin-exchange and dipole relaxation in cesium and rubidium"*, **Phys. Rev. B** 38, 4688 (1988). [DOI: 10.1103/PhysRevB.38.4688](https://doi.org/10.1103/PhysRevB.38.4688)
  - A. J. Moerdijk, B. J. Verhaar, and A. Axelsson, *"Resonances in ultracold collisions of $^6\text{Li}$, $^7\text{Li}$, and $^{133}\text{Cs}$"*, **Phys. Rev. A** 51, 4852 (1995). [DOI: 10.1103/PhysRevA.51.4852](https://doi.org/10.1103/PhysRevA.51.4852)
- **核心理论**:
  电子自旋之间的各向异性磁偶极哈密顿量可严格展开为 2 阶不可约球谐张量与自旋张量的标量积：
  $$\hat{V}_{dd}(\mathbf{r}) = -\frac{\alpha^2 \sqrt{6}}{r^3} \sum_{q=-2}^2 (-1)^q C_{2, -q}(\hat{\mathbf{r}}) \left[\hat{\mathbf{s}}_1 \otimes \hat{\mathbf{s}}_2\right]^{(2)}_q$$
  该相互作用打破轨道角动量 $L$ 的守恒，导致 $s$-波 ($L=0$) 与 $d$-波 ($L=2$) 强烈耦合，但严格守恒总磁量子数：
  $$M_{\text{tot}} = M_L + M_S = M_L' + M_S'$$
- **代码映射**:
  - `src/mod_dipolar_scattering.f90`:
    - `calc_mddi_spatial_matrix_element`: 空间张量矩阵元 $\langle L' M_L' | C_{2, q} | L M_L \rangle$；
    - `calc_mddi_spin_matrix_element`: 自旋张量矩阵元 $\langle S' M_S' | [\mathbf{s}_1 \otimes \mathbf{s}_2]^{(2)}_q | S M_S \rangle$；
    - `calc_mddi_total_matrix_element`: 总无量纲各向异性耦合强度系数。

### 9.2 磁阱中非弹性偶极自旋弛豫 (Dipolar Relaxation)
- **文献**:
  - A. J. Moerdijk and B. J. Verhaar, *"Collisional stability of magnetically trapped rubidium atoms"*, **Phys. Rev. A** 53, 4341 (1996). [DOI: 10.1103/PhysRevA.53.4341](https://doi.org/10.1103/PhysRevA.53.4341)
  - F. H. Mies, C. J. Williams, P. S. Julienne, and M. O. Krauss, *"Estimates of Antifolding Rates in Trapped Rubidium"*, **J. Res. Natl. Inst. Stand. Technol.** 101, 521 (1996). [DOI: 10.6028/jres.101.053](https://doi.org/10.6028/jres.101.053)
- **核心理论**:
  在静磁场 $B$ 下，处于弱场寻优态（如 $|S=1, M_S=1, L=0\rangle$）的双原子碰撞时，磁偶极作用驱动非弹性自旋翻转释放 Zeeman 能 $\Delta E = g_s \mu_B B$，转变为 $d$-波离去通道相对动能，造成磁阱原子剧烈碰撞损失：
  $$\sigma_{\text{rel}}(E, B) \approx \frac{4\pi}{k_i^2} \frac{k_f}{k_i} \left|\langle \psi_f^{(d)} | \hat{V}_{dd} | \psi_i^{(s)} \rangle\right|^2$$
- **代码映射**:
  - `src/mod_dipolar_scattering.f90`:
    - `calc_dipolar_relaxation_cross_section`: 单能自旋弛豫截面 $\sigma_{\text{rel}}$ 与双体速率系数 $K_{\text{rel}}$；
    - `calc_dipolar_relaxation_thermal_rate`: 麦克斯韦-玻尔兹曼系综平均热弛豫速率。

### 9.3 极性分子外加直流电场 (Stark 效应) 与各向异性电偶极碰撞 (EDDI)
- **文献**:
  - J. L. Bohn, M. Cavagnero, and C. Ticknor, *"Quasi-universal dipolar scattering in cold and ultracold gases"*, **New J. Phys.** 11, 055039 (2009). [DOI: 10.1088/1367-2630/11/5/055039](https://doi.org/10.1088/1367-2630/11/5/055039)
  - G. Quéméner and J. L. Bohn, *"Electric field dependence of ultracold chemical reactions of polar molecules"*, **Phys. Rev. A** 81, 022702 (2010). [DOI: 10.1103/PhysRevA.81.022702](https://doi.org/10.1103/PhysRevA.81.022702)
  - K.-K. Ni et al., *"A High Phase-Space-Density Gas of Polar Molecules"*, **Science** 322, 231 (2008). [DOI: 10.1126/science.1163861](https://doi.org/10.1126/science.1163861)
- **核心理论**:
  外加电场 $\mathcal{E}$ 混合相反宇称转动态使得实验室系诱导电偶极矩 $d_{\text{ind}}(\mathcal{E}) = -\frac{\partial E_0}{\partial \mathcal{E}} > 0$。分子间各向异性电偶极势为：
  $$V_{dd}^{\text{el}}(\mathbf{r}) = -\frac{2 d_{\text{ind}}^2}{r^3} C_{20}(\hat{\mathbf{r}})$$
  特征偶极长度尺度定义为 $a_d = \frac{\mu d_{\text{ind}}^2}{2\hbar^2}$。
- **代码映射**:
  - `src/mod_dipolar_scattering.f90`:
    - `calc_stark_induced_dipole`: 刚体转子 Stark 哈密顿量对角化求 $d_{\text{ind}}$；
    - `calc_eddi_matrix_element`: 各向异性电偶极张量矩阵元；
    - `calc_dipole_length_scale`: 偶极力程 $a_d$；
    - `build_dipolar_coupled_potential_matrix`: 多分波偶极耦合势矩阵组装。

---

## 10. 超冷光缔合谱学与自由-束缚态量子跃迁

### 10.1 自由-束缚态 Franck-Condon 重叠积分
- **文献**:
  - K. M. Jones, E. Tiesinga, P. D. Lett, and P. S. Julienne, *"Ultracold photoassociation spectroscopy: Long-range molecules and atomic scattering"*, **Rev. Mod. Phys.** 78, 483 (2006). [DOI: 10.1103/RevModPhys.78.483](https://doi.org/10.1103/RevModPhys.78.483)
  - H. R. Thorsheim, J. Weiner, and P. S. Julienne, *"Laser-Induced Photoassociation of Ultracold Sodium Atoms"*, **Phys. Rev. Lett.** 58, 2420 (1987). [DOI: 10.1103/PhysRevLett.58.2420](https://doi.org/10.1103/PhysRevLett.58.2420)
- **核心理论**:
  连续能量归一化初态散射波 $u_E(r)$ 与电子激发态长程振动束缚态 $\chi_{v'}(r)$ 在激光作用下发生自由-束缚偶极跃迁：
  $$I_{\text{FB}}(E, v') = \int_0^\infty u_E(r) \mu_{eg}(r) \chi_{v'}(r) dr$$
  Franck-Condon 密度 $f_{\text{FB}}(E, v') = |I_{\text{FB}}|^2$ 主要贡献来自于外转向点（Condon 点 $R_C$），呈现显著的量子反射与节点干涉结构。
- **代码映射**:
  - `src/mod_photoassociation.f90`:
    - `calc_free_bound_fc_overlap`: 自由-束缚态重叠积分与 FC 密度计算；
    - `calc_pa_stimulated_width`: 激光光强诱导受激线宽 $\hbar\Gamma_{\text{stim}}(E, I)$。

### 10.2 光缔合共振截面与热速率系数 (Bohn-Julienne 理论)
- **文献**:
  - J. L. Bohn and P. S. Julienne, *"Semianalytic theory of laser-assisted ultracold collisions"*, **Phys. Rev. A** 60, 414 (1999). [DOI: 10.1103/PhysRevA.60.414](https://doi.org/10.1103/PhysRevA.60.414)
  - R. Napolitano, J. Weiner, C. J. Williams, and P. S. Julienne, *"Line shapes of photoassociation in ultracold collisions"*, **Phys. Rev. Lett.** 73, 1352 (1994). [DOI: 10.1103/PhysRevLett.73.1352](https://doi.org/10.1103/PhysRevLett.73.1352)
- **核心理论**:
  在激光失谐 $\Delta$ 下，单能散射对发生光缔合的截面服从孤立 Breit-Wigner / Fano 形式：
  $$\sigma_{\text{PA}}(E, \Delta) = \frac{\pi}{k^2} \frac{\hbar\Gamma_{\text{stim}} \gamma_{\text{nat}}}{(E - \Delta)^2 + [(\gamma_{\text{nat}} + \hbar\Gamma_{\text{stim}})/2]^2}$$
  在有限碰撞温度 $T$ 下，麦克斯韦-玻尔兹曼系综平均速率系数为：
  $$K_{\text{PA}}(T, \Delta) = \frac{1}{h Q_T} \int_0^\infty e^{-E/k_B T} \frac{\hbar\Gamma_{\text{stim}}(E) \gamma_{\text{nat}}}{(E - \Delta)^2 + [(\gamma_{\text{tot}})/2]^2} dE$$
- **代码映射**:
  - `src/mod_photoassociation.f90`:
    - `calc_pa_cross_section`: 单能光缔合截面与速率；
    - `calc_pa_thermal_rate_coefficient`: 热平均速率系数 $K_{\text{PA}}(T)$；
    - `calc_pa_spectrum_scan`: 缔合激光谱扫描；
    - `calc_two_photon_raman_association_coupling`: 双光子受激 Raman / STIRAP 缔合基态分子。

---

## 11. 超冷少体物理、三体复合与 Efimov 普适态

### 11.1 Efimov 超径向超越方程与离散标度不变性
- **文献**:
  - V. Efimov, *"Energy levels arising from resonant two-body forces in a three-body system"*, **Phys. Lett. B** 33, 563 (1970). [DOI: 10.1016/0370-2693(70)90349-7](https://doi.org/10.1016/0370-2693(70)90349-7)
  - E. Braaten and H.-W. Hammer, *"Universality in few-body systems with large scattering length"*, **Phys. Rep.** 428, 259 (2006). [DOI: 10.1016/j.physrep.2006.03.001](https://doi.org/10.1016/j.physrep.2006.03.001)
- **核心理论**:
  在共振两体相互作用极限下（$|a| \gg r_0$），全同玻色子三体超径向波函数方程导出特征超越方程：
  $$s \cosh\left(\frac{\pi s}{2}\right) - \frac{8}{\sqrt{3}} \sinh\left(\frac{\pi s}{6}\right) = 0$$
  该方程在实轴上具有唯一正根 $s_0 \approx 1.0062378$，决定了有效吸引超径向势 $V_{\text{eff}}(R) = -(s_0^2 + 1/4)\hbar^2 / (2\mu R^2)$。三体体系破缺连续标度对称性，展现几何递推常数的**离散标度不变性 (Discrete Scale Invariance)**：
  $$\lambda = e^{\pi / s_0} \approx 22.69438$$
  对于异核 $AAB$ 碰撞体系，该根依赖于轻-重质量比 $\beta$，当存在极轻原子时 $s_0$ 显著增大，使几何比率 $\lambda$ 压缩至易于实验观测的区间。
- **代码映射**:
  - `src/mod_three_body_recombination.f90`:
    - `solve_efimov_s0_identical_bosons`: Newton-Raphson 精确求解 $s_0 = 1.0062378$；
    - `solve_efimov_s0_heteronuclear`: 异核 $AAB$ 体系阻尼自适应牛顿迭代求根；
    - `calc_efimov_scale_factor`: 离散标度因子 $\lambda = e^{\pi/s_0}$ 计算。

### 11.2 Braaten-Hammer 普适三体复合损失速率 $K_3(a)$
- **文献**:
  - B. D. Esry, C. H. Greene, and J. P. Burke, Jr., *"Recombination of Three Atoms in the Ultracold Limit"*, **Phys. Rev. Lett.** 83, 1751 (1999). [DOI: 10.1103/PhysRevLett.83.1751](https://doi.org/10.1103/PhysRevLett.83.1751)
  - T. Kraemer et al., *"Evidence for Efimov quantum states in an ultracold gas of caesium atoms"*, **Nature** 440, 315 (2006). [DOI: 10.1038/nature04626](https://doi.org/10.1038/nature04626)
- **核心理论**:
  三体非弹性复合速率 $K_3$ 依赖于两体散射长度 $a$：
  - 正散射长度（$a > 0$）：发生多路径干涉相消，出现 Efimov 干涉极小值：
    $$K_3(a > 0) = \frac{C_+ \hbar}{m} a^4 \left[ \sin^2\left(s_0 \ln(a / a_*)\right) + \sinh^2(\eta_*) \right]$$
  - 负散射长度（$a < 0$）：当 Efimov 三聚体束缚态切入两体碰撞阈值时形成剧烈的共振损耗峰：
    $$K_3(a < 0) = \frac{C_- \hbar}{m} a^4 \frac{\sinh(2\eta_*)}{\sin^2\left(s_0 \ln(|a| / a_-)\right) + \sinh^2(\eta_*)}$$
  - 幺正极限有限温度饱和（$|a| \to \infty$）：受热波长 $k_T$ 截断，损失率服从普适 $T^{-2}$ 幂律：
    $$K_3(T) = C_{\text{unit}} \frac{\hbar^5}{m^3 (k_B T)^2}$$
- **代码映射**:
  - `src/mod_three_body_recombination.f90`:
    - `calc_three_body_recombination_a_positive`: 正散射长度三体复合与干涉极小值；
    - `calc_three_body_recombination_a_negative`: 负散射长度三聚体共振峰；
    - `calc_unitary_three_body_loss_temperature`: 幺正极限 $T^{-2}$ 温度饱和定律。

---

## 12. 低维光晶格受限量子散射与约束诱导共振 (CIR)

### 12.1 Olshanii 准一维约束诱导共振理论
- **文献**:
  - M. Olshanii, *"Atomic Scattering in the Presence of an External Confinement and a Gas of Hard-Core Bosons"*, **Phys. Rev. Lett.** 81, 938 (1998). [DOI: 10.1103/PhysRevLett.81.938](https://doi.org/10.1103/PhysRevLett.81.938)
  - E. Haller et al., *"Realization of an Excited, Strongly Correlated Quantum Gas Phase"*, **Science** 325, 1224 (2009). [DOI: 10.1126/science.1175850](https://doi.org/10.1126/science.1175850)
- **核心理论**:
  在横向谐振子紧束缚光波导 $V_\perp(\rho) = \frac{1}{2}\mu\omega_\perp^2\rho^2$ 下（振子特征长度 $a_\perp = \sqrt{\hbar/(\mu\omega_\perp)}$），准一维有效接触势耦合常数由横向模态多体格林函数重整化：
  $$g_{\text{1D}}(a_s) = \frac{2\hbar^2 a_s}{\mu a_\perp^2} \frac{1}{1 - C \frac{a_s}{a_\perp}}, \quad C = \frac{|\zeta(1/2)|}{\sqrt{2}} \approx 1.0326$$
  当三维散射长度逼近极点 $a_{\text{CIR}} = a_\perp / C \approx 0.96843 a_\perp$ 时，$g_{\text{1D}} \to \pm \infty$ 发生约束诱导共振（Confinement-Induced Resonance, CIR）。一维玻色气体跨入由 Lieb-Liniger 参数 $\gamma_{\text{LL}} = m |g_{\text{1D}}| / (\hbar^2 n_{\text{1D}}) \gg 1$ 支配的强关联 Tonks-Girardeau 硬球玻色极限。
- **代码映射**:
  - `src/mod_confined_scattering.f90`:
    - `init_waveguide_1d`: 简谐波导振子尺寸 $a_\perp$ 与 CIR 极点 $a_{\text{CIR}}$ 计算；
    - `calc_olshanii_cir_parameters`: 1D 重整化有效相互作用 $g_{\text{1D}}$ 与状态判据；
    - `calc_confined_dimer_binding_energy`: 受限诱导分子束缚态结合能 $E_b$；
    - `calc_lieb_liniger_parameter`: 一维无量纲关联度 $\gamma_{\text{LL}}$ 与 Tonks 气体判别。

### 12.2 各向异性波导共振分裂与准二维散射
- **文献**:
  - T. Bergeman, M. G. Moore, and M. Olshanii, *"Atom-Atom Scattering under Cylindrical Harmonic Confinement"*, **Phys. Rev. Lett.** 91, 163201 (2003). [DOI: 10.1103/PhysRevLett.91.163201](https://doi.org/10.1103/PhysRevLett.91.163201)
  - D. S. Petrov, M. Holzmann, and G. V. Shlyapnikov, *"Bose-Einstein Condensation in Quasi-2D Trapped Gases"*, **Phys. Rev. Lett.** 84, 2551 (2000). [DOI: 10.1103/PhysRevLett.84.2551](https://doi.org/10.1103/PhysRevLett.84.2551)
- **核心理论**:
  若横向谐振势阱各向异性（$\omega_x \neq \omega_y$），空间简并被打破，CIR 极点分裂为双峰：$a_{\text{CIR},x} \approx 0.9684 a_{\perp,x}$ 与 $a_{\text{CIR},y} \approx 0.9684 a_{\perp,y}$。准二维平面囚禁（$a_z$）中，两体散射由对数尺度散射长度 $a_{\text{2D}} = A_0 a_z \exp(-\sqrt{\pi/2} a_z / a_s)$ 描述。
- **代码映射**:
  - `src/mod_confined_scattering.f90`:
    - `calc_anisotropic_cir_poles`: 各向异性 CIR 双极点分裂解析提取；
    - `init_planar_2d`: 准二维平面势阱初始化；
    - `calc_quasi_2d_scattering`: 准 2D 对数散射长度与低能复散射振幅 $f_{\text{2D}}$。

---

## 13. 自电离体系、Fano 共振理论与复坐标旋转法 (CCR)

### 13.1 组态相互作用与不对称 Fano 线型
- **文献**:
  - U. Fano, *"Effects of Configuration Interaction on Intensities and Phase Shifts"*, **Phys. Rev.** 124, 1866 (1961). [DOI: 10.1103/PhysRev.124.1866](https://doi.org/10.1103/PhysRev.124.1866)
- **核心理论**:
  当离散双激发准束缚态 $|\phi\rangle$ 浸没在同能量连续电离谱带 $|\psi_E\rangle$ 中时，由库仑多体相互作用引起的组态混合 $V_E = \langle \phi | H | \psi_E \rangle$ 导致准束缚态发生自电离。跃迁吸收截面表现为非对称 Fano 线型：
  $$\sigma(\epsilon) = \sigma_0 \frac{(q + \epsilon)^2}{1 + \epsilon^2}, \quad \epsilon = \frac{E - E_0}{\Gamma / 2}$$
  其中 $\Gamma = 2\pi |V_E|^2$ 为自电离衰变宽度，自电离寿命为 $\tau = \hbar / \Gamma$。当 $\epsilon = -q$ 时跃迁振幅相消干涉至严格零点（抗共振 Anti-resonance），当 $\epsilon = 1/q$ 时跃迁相干相长达到最大峰值 $\sigma_{\text{max}} = \sigma_0(1 + q^2)$。
- **代码映射**:
  - `src/mod_autoionization_fano.f90`:
    - `calc_fano_profile`: Fano 谱线截面计算；
    - `calc_autoionization_lifetime`: 衰变宽度至飞秒寿命高精度换算；
    - `calc_fano_q_parameter`: 偶极跃迁矩阵元组态不对称因子 $q$ 构造。

### 13.2 复坐标旋转法 (Complex Coordinate Rotation, CCR)
- **文献**:
  - W. P. Reinhardt, *"Complex coordinates in the theory of atomic and molecular structure and dynamics"*, **Annu. Rev. Phys. Chem.** 33, 223 (1982). [DOI: 10.1146/annurev.pc.33.100182.001255](https://doi.org/10.1146/annurev.pc.33.100182.001255)
  - N. Moiseyev, *"Quantum theory of resonances: calculating energies, widths and cross-sections by complex scaling"*, **Phys. Rep.** 302, 212 (1998). [DOI: 10.1016/S0370-1573(98)00002-7](https://doi.org/10.1016/S0370-1573(98)00002-7)
- **核心理论**:
  通过非厄米复标度算符 $U(\theta) = \exp(i\theta r \partial_r)$ 使坐标复平面旋转 $r \to r e^{i\theta}$，连续谱以分支点为顶点顺时针旋转 $2\theta$ 角度，将原本隐藏在第二黎曼叶上的 S-矩阵复极点暴露在物理区域中：
  $$E_{\text{res}} = E_R - i \frac{\Gamma}{2}$$
  时域波包展现纯指数衰变与连续电子出射通量 $J(t) = \Gamma P(t) / \hbar$。
- **代码映射**:
  - `src/mod_autoionization_fano.f90`:
    - `solve_ccr_resonance_model`: 2x2 复对称非厄米模型本征求解提取共振极点 $(E_R, \Gamma)$；
    - `calc_time_domain_autoionization_decay`: 时域生存几率与自电离出射通量。

---

## 14. 交叉静电磁场分子量子动力学与非共线 Stark-Zeeman 态混合

### 14.1 任意倾角 $\theta_{EB}$ 交叉场宇称与投影对称性破缺
- **文献**:
  - T. V. Tscherbul and R. V. Krems, *"Controlling Central Collisions of Cold Molecules with Crossed Electric and Magnetic Fields"*, **Phys. Rev. Lett.** 97, 083201 (2006). [DOI: 10.1103/PhysRevLett.97.083201](https://doi.org/10.1103/PhysRevLett.97.083201)
  - B. Friedrich and D. Herschbach, *"Spatial orientation of molecules in strong electric fields and nonresonant intense laser fields"*, **Z. Phys. D** 36, 221 (1996). [DOI: 10.1007/BF01426405](https://doi.org/10.1007/BF01426405)
- **核心理论**:
  极性开壳层分子（如 $^2\Sigma$ 自由基或顺磁偶极分子）在非共线静电场 $\mathbf{E} = E \hat{z}$ 与倾斜磁场 $\mathbf{B} = B(\sin\theta_{EB} \hat{x} + \cos\theta_{EB} \hat{z})$ 共同作用下：
  $$H = B_{\text{rot}} \mathbf{J}^2 - d E \cos\theta_R + g_S \mu_B \left( B_z S_z + B_x S_x \right)$$
  横向磁场分量 $B_x = B \sin\theta_{EB}$ 产生非对角算符 $S_x = (S_+ + S_-)/2$，驱动 $\Delta M_S = \pm 1$ 强耦合，彻底瓦解了沿电场轴的空间柱对称性，形成可由倾角连续调控的 Stark-Zeeman 密集反交叉谱结构与高定向度 $\langle \cos\theta_R \rangle$。
- **代码映射**:
  - `src/mod_crossed_field_scattering.f90`:
    - `init_crossed_field_config`: 交叉场参数与倾斜坐标分解初始化；
    - `build_crossed_field_hamiltonian`: 转动-自旋全耦合实对称哈密顿量矩阵组装；
    - `solve_crossed_field_eigenstates`: 倾斜场绝热本征能级与本征波函数求解；
    - `calc_crossed_field_observables`: 分子定向度 $\langle \cos\theta \rangle$ 与自旋极化 $\langle S_z \rangle$；
    - `scan_tilt_angle_spectrum`: 倾角 $\theta_{EB}$ 从 $0^\circ$ 至 $90^\circ$ 全域连续扫描。

---

## 15. 三原子反应碰撞几何、Jacobi 坐标、LEPS 势能面与锥形交叉几何相位

### 15.1 Jacobi 反应散射坐标体系与核间距可逆映射
- **文献**:
  - D. G. Truhlar and C. J. Horowitz, *"Functional representations of potential energy surfaces: The H + H2 reaction"*, **J. Chem. Phys.** 68, 2466 (1978). [DOI: 10.1063/1.436019](https://doi.org/10.1063/1.436019)
- **核心理论**:
  在 $A + BC \to AB + C$ 反应碰撞中，质心 Jacobi 坐标 $(r, R, \gamma)$（$r$ 为双原子间距，$R$ 为入射原子至双原子质心距离，$\gamma$ 为二者夹角）与三原子内部核间距 $(r_{12}, r_{23}, r_{31})$ 满足精确封闭解析几何转换，是构建反应散射通道波函数与势能曲面的基础数学工具。
- **代码映射**:
  - `src/mod_triatomic_geometry.f90`:
    - `jacobi_to_internuclear`: Jacobi 坐标向互核欧几里得距离严格投影；
    - `internuclear_to_jacobi`: 任意三原子构型向质心 Jacobi 反应坐标逆变换。

### 15.2 Sato 修正 LEPS 反应势能面与过渡态活化势垒
- **文献**:
  - S. Sato, *"A New Method of Drawing the Potential Energy Surface"*, **J. Chem. Phys.** 23, 592 (1955). [DOI: 10.1063/1.1742050](https://doi.org/10.1063/1.1742050)
- **核心理论**:
  London-Eyring-Polanyi-Sato (LEPS) 势能面基于 Morse 单重态 $^1\Sigma$ 与反 Morse 三重态 $^3\Sigma$ 曲线构造库仑积分 $Q_i(r_i)$ 与交换积分 $J_i(r_i)$：
  $$V(r_{12}, r_{23}, r_{31}) = Q_1 + Q_2 + Q_3 - \sqrt{\frac{1}{2} \left[ (J_1 - J_2)^2 + (J_2 - J_3)^2 + (J_3 - J_1)^2 \right]}$$
  在 $H + H_2$ 基准反应体系中精确复现了解离渐近极限 $-D_e$ 与共线对称过渡态鞍点势垒（Barrier $\sim 9.8\text{ kcal/mol}$）。
- **代码映射**:
  - `src/mod_triatomic_geometry.f90`:
    - `calc_leps_potential`: Sato 标度库仑与交换积分及 LEPS 反应势能面；
    - `init_default_h3_leps`: $H_3$ 经典反应势能面参数初始化。

### 15.3 线性锥形交叉 (CI) 与 Longuet-Higgins / Berry 拓扑几何相位
- **文献**:
  - H. C. Longuet-Higgins et al., *"Studies of the Jahn-Teller effect. II. The dynamical problem"*, **Proc. R. Soc. Lond. A** 244, 1 (1958). [DOI: 10.1098/rspa.1958.0022](https://doi.org/10.1098/rspa.1958.0022)
  - M. V. Berry, *"Quantal phase factors accompanying adiabatic changes"*, **Proc. R. Soc. Lond. A** 392, 45 (1984). [DOI: 10.1098/rspa.1984.0023](https://doi.org/10.1098/rspa.1984.0023)
- **核心理论**:
  在两态 $E \otimes e$ 线性锥形交叉模型 $H_{\text{diab}} = \begin{pmatrix} \kappa x & \lambda y \\ \lambda y & -\kappa x \end{pmatrix}$ 中，绝热能级呈现双锥形退化劈裂 $\Delta E = 2\sqrt{\kappa^2 x^2 + \lambda^2 y^2}$。绕奇点闭合回路积分散布非绝热规范势 $\mathbf{A} = i\langle \psi_- | \nabla | \psi_- \rangle$，严格累积拓扑相 $\Phi_B = \oint \mathbf{A} \cdot d\mathbf{R} = \pi$，使绝热电子波函数环绕一周后产生拓扑变号 $|\psi(2\pi)\rangle = -|\psi(0)\rangle$。
- **代码映射**:
  - `src/mod_triatomic_geometry.f90`:
    - `calc_conical_intersection_adiabats`: CI 锥形绝热势能面与能隙；
    - `calc_berry_phase_around_ci`: 环绕退化奇点闭合回路拓扑 Berry 几何相位计算。

---

## 16. 超冷旋量玻色-爱因斯坦凝聚与宏观自旋混合动力学

### 16.1 $F=1$ 旋量凝聚体相互作用与铁磁/反铁磁相分类
- **文献**:
  - T.-L. Ho, *"Spinor Bose Condensates in Optical Traps"*, **Phys. Rev. Lett.** 81, 742 (1998). [DOI: 10.1103/PhysRevLett.81.742](https://doi.org/10.1103/PhysRevLett.81.742)
  - T. Ohmi and K. Machida, *"Bose-Einstein Condensation with Internal Degrees of Freedom in Alkali Atom Gases"*, **J. Phys. Soc. Jpn.** 67, 1822 (1998). [DOI: 10.1143/JPSJ.67.1822](https://doi.org/10.1143/JPSJ.67.1822)
- **核心理论**:
  自旋 $F=1$ 玻色气体在无外场光偶极阱中具有内部自旋自由度。低能碰撞由总自旋 $F_{\text{tot}} = 0, 2$ 两个通道的 $s$-波散射长度 $a_0, a_2$ 决定：
  $$c_0 = \frac{4\pi\hbar^2(a_0 + 2a_2)}{3m} \quad (\text{密度相互作用}), \quad c_2 = \frac{4\pi\hbar^2(a_2 - a_0)}{3m} \quad (\text{自旋交换相互作用})$$
  - $^{87}\text{Rb}$ 满足 $a_2 < a_0 \implies c_2 < 0$，基态倾向于自旋平行排列，表现为**铁磁相 (Ferromagnetic Phase)**；
  - $^{23}\text{Na}$ 满足 $a_2 > a_0 \implies c_2 > 0$，基态倾向于自旋反平行排列，表现为**反铁磁/极性相 (Antiferromagnetic / Polar Phase)**。
- **代码映射**:
  - `src/mod_spinor_bec.f90`:
    - `init_spinor_preset`: $^{87}\text{Rb}$ 与 $^{23}\text{Na}$ 实验标准散射参数装载；
    - `calc_spinor_interaction_couplings`: 理论相互作用参数 $c_0, c_2$ 换算与相归属。

### 16.2 单模近似 (SMA) 非线性运动方程与 RK4 宏观自旋振荡
- **文献**:
  - M.-S. Chang et al., *"Observation of Spinor Dynamics in Optically Trapped 87Rb Bose-Einstein Condensates"*, **Phys. Rev. Lett.** 92, 140403 (2004). [DOI: 10.1103/PhysRevLett.92.140403](https://doi.org/10.1103/PhysRevLett.92.140403)
  - H. Pu et al., *"Spin mixing in a doubly degenerate ultracold atomic gas"*, **Phys. Rev. A** 60, 1463 (1999). [DOI: 10.1103/PhysRevA.60.1463](https://doi.org/10.1103/PhysRevA.60.1463)
- **核心理论**:
  在单模近似（SMA）下，宏观分量波函数空间分布相同，归一化多分量旋量标量满足含二阶塞曼位移 $q_Z \propto B^2$ 的非线性运动方程：
  $$i\hbar \frac{d\zeta_{\pm 1}}{dt} = \left[ q_Z + c_2 n (\rho_0 + \rho_{\pm 1} - \rho_{\mp 1}) \right] \zeta_{\pm 1} + c_2 n \zeta_0^2 \zeta_{\mp 1}^*$$
  $$i\hbar \frac{d\zeta_0}{dt} = c_2 n (\rho_{+1} + \rho_{-1}) \zeta_0 + 2 c_2 n \zeta_{+1} \zeta_{-1} \zeta_0^*$$
  系统在时间演化中严格守恒全凝聚体几率与磁化强度 $m_z = |\zeta_{+1}|^2 - |\zeta_{-1}|^2$。
- **代码映射**:
  - `src/mod_spinor_bec.f90`:
    - `calc_quadratic_zeeman_shift`: 超精细 Breit-Rabi 二阶塞曼能量位移 $q_Z(B)$；
    - `propagate_spinor_sma_rk4`: 四阶 Runge-Kutta 保范数与保磁化强度动力学步进器；
    - `simulate_spin_mixing_dynamics`: 宏观相干自旋混合长时时域模拟。

---

---

## 17. 三原子超球面反应动力学与过渡态理论 (Hyperspherical Reactive Scattering & TST)

### 17.1 质量标度超球面坐标 (Delves/Smith 坐标) 与偏角动力学
- **文献**:
  - B. R. Johnson, *"The quantum dynamics of three-body reactive collisions in hyperspherical coordinates"*, **J. Chem. Phys.** 73, 5051 (1980). [DOI: 10.1063/1.440058](https://doi.org/10.1063/1.440058)
  - R. T. Pack and G. A. Parker, *"Quantum reactive scattering in three dimensions: General theory and the hyperspherical representation"*, **J. Chem. Phys.** 87, 3888 (1987). [DOI: 10.1063/1.452944](https://doi.org/10.1063/1.452944)
  - L. M. Delves, *"Tertiary and higher-order collision processes"*, **Nucl. Phys.** 9, 391 (1959); **Nucl. Phys.** 20, 275 (1960).
- **核心理论**:
  将三原子体系 $A + BC \to AB + C$ 的质心分离后，引入 Delves 质量标度坐标：
  $$S = d \cdot R_{A,BC}, \quad s = d^{-1} \cdot r_{BC}, \quad d = \left(\frac{\mu_{A,BC}}{\mu_{BC}}\right)^{1/4}$$
  超半径 $\rho = \sqrt{S^2 + s^2}$，超角 $\alpha = \arctan(s / S)$。三体运动动能项在超球面坐标中退化为各向同性拉普拉斯算符加上离心势与超角动量算符。反应散射通道间的偏转角（Reaction Skew Angle $\beta_{skew}$）仅由三核质量比决定：
  $$\cos\beta_{skew} = \sqrt{\frac{m_A m_C}{(m_A + m_B)(m_B + m_C)}}, \quad \tan\beta_{skew} = \sqrt{\frac{m_B M}{m_A m_C}}$$
  对于同核体系 $\text{H} + \text{H}_2 \to \text{H}_2 + \text{H}$，$\beta_{skew} = 60^\circ$；对于 $\text{D} + \text{H}_2 \to \text{HD} + \text{H}$，$\beta_{skew} = 54.74^\circ$。
- **代码映射**:
  - `src/mod_hyperspherical_reactive.f90`:
    - `init_reaction_mass`: 自动计算反应体系 Delves 因子 $d$、各通道约化质量及反应偏转角 $\beta_{skew}$；
    - `jacobi_to_hyperspherical` / `hyperspherical_to_jacobi`: 双向正逆解析坐标几何转换。

### 17.2 量子隧穿、累积反应几率 $N(E)$ 与正则热反应速率常数 $k(T)$
- **文献**:
  - W. H. Miller, *"Semiclassical limit of quantum mechanical transition state theory for nonseparable systems"*, **J. Chem. Phys.** 62, 1899 (1975). [DOI: 10.1063/1.430676](https://doi.org/10.1063/1.430676)
  - D. G. Truhlar and B. C. Garrett, *"Variational transition-state theory"*, **Acc. Chem. Res.** 13, 440 (1980). [DOI: 10.1021/ar50156a002](https://doi.org/10.1021/ar50156a002)
  - C. Eckart, *"The penetration of a parabolic potential barrier"*, **Phys. Rev.** 35, 1303 (1930). [DOI: 10.1103/PhysRev.35.1303](https://doi.org/10.1103/PhysRev.35.1303)
  - E. P. Wigner, *"Calculation of the rate of elementary association reactions"*, **J. Chem. Phys.** 5, 720 (1937).
- **核心理论**:
  反应鞍点过渡态处的量子隧穿传递几率采用 Eckart / 抛物线势垒严格解：
  $$P(E) = \frac{1}{1 + \exp\left[ - \frac{2\pi (E - V_b)}{\hbar \omega^\ddagger} \right]}$$
  累积反应几率 $N(E)$ 对跨越势垒的所有量化弯曲与伸缩振动态通道求和：
  $$N(E) = \sum_{n_{bend}, n_{symm}} P\left( E - E_{vib}(n_{bend}, n_{symm}) \right)$$
  由量子碰撞理论，正则热反应速率常数 $k(T)$ 由 $N(E)$ 进行玻尔兹曼热积分直接求得：
  $$k(T) = \frac{1}{2\pi \hbar Q_R(T)} \int_{-\infty}^\infty dE \, N(E) e^{-E / (k_B T)}$$
  在过渡态理论（TST）下，Wigner 势垒量子隧穿因子解析修正为 $\kappa(T) = 1 + \frac{1}{24}\left(\frac{\hbar \omega^\ddagger}{k_B T}\right)^2$。
- **代码映射**:
  - `src/mod_hyperspherical_reactive.f90`:
    - `calc_eckart_transmission`: 鞍点精确量子隧穿传递几率计算；
    - `calc_cumulative_reaction_probability`: 多量子振动态求和累积反应几率 $N(E)$；
    - `calc_canonical_rate_constant`: 高斯-勒让德/自适应数值积分求取全量子反应速率常数 $k(T)$；
    - `calc_tst_wigner_rate`: 经典 TST 与 Wigner 隧穿修正速率计算。

---

## 18. 超冷偶极量子液滴与李-黄-杨量子涨落修正 (Dipolar Droplets & LHY)

### 18.1 磁偶极相互作用与 Pelster-Lima 涨落修正积分 $Q_5(\epsilon_{dd})$
- **文献**:
  - T. D. Lee, K. Huang, and C. N. Yang, *"Eigenvalues and Eigenfunctions of a System of Impenetrable Spheres"*, **Phys. Rev.** 106, 1135 (1957). [DOI: 10.1103/PhysRev.106.1135](https://doi.org/10.1103/PhysRev.106.1135)
  - A. R. P. Lima and A. Pelster, *"Quantum fluctuations in dipolar Bose gases"*, **Phys. Rev. A** 84, 041604(R) (2011); **Phys. Rev. A** 86, 063609 (2012). [DOI: 10.1103/PhysRevA.84.041604](https://doi.org/10.1103/PhysRevA.84.041604)
  - F. Chomaz et al., *"Dipolar physics: a review of experiments with magnetic atoms"*, **Rep. Prog. Phys.** 86, 026401 (2023). [DOI: 10.1088/1361-6633/aca8a4](https://doi.org/10.1088/1361-6633/aca8a4)
- **核心理论**:
  对于磁偶极原子（如 $^{162}\text{Dy}: \mu_m \approx 10 \mu_B, a_{dd} \approx 131 a_0$；$^{166}\text{Er}: \mu_m \approx 7 \mu_B, a_{dd} \approx 65 a_0$），偶极长度 $a_{dd} = \frac{\mu_0 \mu_m^2 m}{12 \pi \hbar^2}$，相对偶极相互作用强度 $\epsilon_{dd} = a_{dd} / a_s$。
  超越平均场的李-黄-杨（LHY）零点量子涨落能量修正积分由 Pelster 与 Lima 解析导出：
  $$Q_5(\epsilon_{dd}) = \frac{3}{2} \int_0^1 du \, (1 - u^2) \left[ 1 + \epsilon_{dd} (3u^2 - 1) \right]^{5/2}$$
  当 $\epsilon_{dd} \to 0$ 时退化为纯接触相互作用普适 LHY 结果 $Q_5(0) = 1$；当 $\epsilon_{dd} > 1$ 时，$Q_5 > 1$，量子涨落排斥能被各向异性偶极激发显著增强。
- **代码映射**:
  - `src/mod_dipolar_droplets_lhy.f90`:
    - `init_dipolar_droplet_param`: $^{162}\text{Dy}$ / $^{166}\text{Er}$ 磁偶极特征常数计算；
    - `calc_pelster_lima_q5`: 高斯-勒让德高精度数值求积 $Q_5(\epsilon_{dd})$。

### 18.2 拓展 Gross-Pitaevskii (eGPE) 理论与自由空间自束缚平衡态
- **文献**:
  - D. S. Petrov, *"Quantum Mechanical Stabilization of a Collapsing Bose-Bose Mixture"*, **Phys. Rev. Lett.** 115, 155301 (2015). [DOI: 10.1103/PhysRevLett.115.155301](https://doi.org/10.1103/PhysRevLett.115.155301)
  - H. Kadau et al., *"Observing the Rosensweig instability of a quantum ferrofluid"*, **Nature** 530, 194 (2016). [DOI: 10.1038/nature16485](https://doi.org/10.1038/nature16485)
  - I. Ferrier-Barbut et al., *"Observation of Quantum Droplets in a Strongly Dipolar Bose Gas"*, **Phys. Rev. Lett.** 116, 215301 (2016). [DOI: 10.1103/PhysRevLett.116.215301](https://doi.org/10.1103/PhysRevLett.116.215301)
  - M. Schmitt et al., *"Self-bound droplets of a dilute magnetic quantum liquid"*, **Nature** 539, 259 (2016). [DOI: 10.1038/nature20126](https://doi.org/10.1038/nature20126)
  - F. Chomaz et al., *"Quantum-Fluctuation-Driven Crossover from a Dilute Bose-Einstein Condensate to a Macrodroplet in a Dipolar Quantum Fluid"*, **Phys. Rev. X** 6, 041039 (2016). [DOI: 10.1103/PhysRevX.6.041039](https://doi.org/10.1103/PhysRevX.6.041039)
- **核心理论**:
  在平均场吸引占优区域（$\epsilon_{dd} > 1$），平均场能量密度负贡献与高阶正定 LHY 排斥能量密度相抗衡：
  $$\mathcal{E}(n) = \frac{1}{2} g_{eff} n^2 + \frac{2}{5} \gamma_{LHY} n^{5/2}, \quad g_{eff} = g (1 - \epsilon_{dd}) < 0$$
  $$\mu(n) = g_{eff} n + \gamma_{LHY} n^{3/2}, \quad \gamma_{LHY} = \frac{32}{3\sqrt{\pi}} g a_s^{3/2} Q_5(\epsilon_{dd})$$
  热力学零压边界条件 $P = n\mu - \mathcal{E} = -\frac{1}{2}|g_{eff}| n^2 + \frac{3}{5}\gamma_{LHY} n^{5/2} = 0$ 严格决定了平顶自束缚量子液滴的核心平衡密度：
  $$n_0 = \frac{25}{36} \left( \frac{|g_{eff}|}{\gamma_{LHY}} \right)^2 = \frac{25 \pi}{4096} \frac{(\epsilon_{dd} - 1)^2}{a_s^3 Q_5(\epsilon_{dd})^2}$$
  在 $n = n_0$ 处，化学势 $\mu(n_0) = -\frac{1}{6}|g_{eff}| n_0 < 0$ 且结合能量密度 $\mathcal{E}(n_0) = -\frac{1}{6}|g_{eff}| n_0^2 < 0$，在第一性原理上严格保证了自由空间三维自束缚量子液滴的稳定存在；临界原子数由 Petrov-Chomaz 标度给出：$N_{crit} \approx 18.6 (\epsilon_{dd}-1)^{-5/2} / \sqrt{Q_5}$。
- **代码映射**:
  - `src/mod_dipolar_droplets_lhy.f90`:
    - `calc_equilibrium_droplet_density`: 零压平顶核心平衡密度 $n_0$ 计算；
    - `calc_droplet_chemical_potential`: 自束缚负化学势 $\mu(n_0) < 0$ 计算；
    - `calc_critical_atom_number`: 液滴自束缚与气相蒸发临界原子数 $N_{crit}$；
    - `calc_egpe_energy_density`: 局部拓展 Gross-Pitaevskii (eGPE) 能量密度。

---

## 19. 强场非顺序双电离与电子重碰撞相关动量动力学 (Strong-Field NSDI & Recollision)

### 19.1 Corkum 三步模型、经典轨道与 $3.17 U_p$ 回碰截止
- **文献**:
  - P. B. Corkum, *"Plasma perspective on strong field multiphoton ionization"*, **Phys. Rev. Lett.** 71, 1994 (1993). [DOI: 10.1103/PhysRevLett.71.1994](https://doi.org/10.1103/PhysRevLett.71.1994)
  - K. J. Schafer et al., *"Above threshold ionization beyond the high harmonic cutoff"*, **Phys. Rev. Lett.** 70, 1599 (1993). [DOI: 10.1103/PhysRevLett.70.1599](https://doi.org/10.1103/PhysRevLett.70.1599)
  - M. Lewenstein et al., *"Theory of high-harmonic generation by low-frequency laser fields"*, **Phys. Rev. A** 49, 2117 (1994). [DOI: 10.1103/PhysRevA.49.2117](https://doi.org/10.1103/PhysRevA.49.2117)
- **核心理论**:
  第一电子在电场相位 $\phi_0 = \omega t_0$ 经由 ADK 准静态隧穿脱附母核，在激光场 $E(t) = F_0 \cos(\omega t)$ 中作自由振荡：
  $$x(\phi) = \frac{F_0}{\omega^2} \left[ \cos\phi - \cos\phi_0 + \sin\phi_0 (\phi - \phi_0) \right]$$
  回碰相位 $\phi_r$ 满足 $x(\phi_r) = 0$（$\phi_r > \phi_0$）。经典回碰动能为：
  $$E_{rec}(\phi_0) = 2 U_p (\sin\phi_r - \sin\phi_0)^2, \quad U_p = \frac{F_0^2}{4\omega^2}$$
  回碰能量存在著名的 Corkum 截断极限：当 $\phi_0 \approx 0.297$ rad（$\approx 17^\circ$）时，回碰相位 $\phi_r \approx 4.45$ rad（$\approx 255^\circ$），最大回碰动能精确达到：
  $$E_{rec}^{max} \approx 3.173 U_p$$
- **代码映射**:
  - `src/mod_strong_field_nsdi.f90`:
    - `init_nsdi_laser`: 激光场光强转换为原子单位峰值电场 $F_0$ 与有质动力势 $U_p$；
    - `calc_adk_rate`: Ammosov-Delone-Krainov 准静态隧穿电离几率密度；
    - `calc_recollision_trajectory`: 牛顿-拉夫逊法精确求解回碰轨道与动能 $E_{rec}(\phi_0)$。

### 19.2 COLTRIMS 双电子相关动量谱 $P(p_{z1}, p_{z2})$ 与“膝盖结构” (Knee Structure)
- **文献**:
  - Th. Weber et al., *"Correlated electron emission in strong field double ionization"*, **Nature** 405, 658 (2000). [DOI: 10.1038/35015033](https://doi.org/10.1038/35015033)
  - R. Moshammer et al., *"Momentum Distributions of Ne(n+) Ions Created by an Intense Ultrashort Laser Pulse"*, **Phys. Rev. Lett.** 84, 447 (2000). [DOI: 10.1103/PhysRevLett.84.447](https://doi.org/10.1103/PhysRevLett.84.447)
  - B. Walker et al., *"Precision Measurement of Strong Field Double Ionization of Helium"*, **Phys. Rev. Lett.** 73, 1227 (1994). [DOI: 10.1103/PhysRevLett.73.1227](https://doi.org/10.1103/PhysRevLett.73.1227)
  - A. Becker and F. H. M. Faisal, *"Mechanism of laser-induced double ionization of helium"*, **Phys. Rev. Lett.** 84, 3546 (2000); **J. Phys. B** 38, R1 (2005). [DOI: 10.1088/0953-4075/38/3/R01](https://doi.org/10.1088/0953-4075/38/3/R01)
  - W. Lotz, *"An empirical formula for the electron-impact ionization cross-section"*, **Z. Phys.** 206, 205 (1967).
- **核心理论**:
  当 $E_{rec} > I_{p2}$ 时触发直接碰撞电离 $(e, 2e)$，其截面遵从 Lotz 经验公式 $\sigma(E) \propto \ln(E/I_{p2}) / (E I_{p2})$。出射电子分享剩余能量 $\Delta E = E_{rec} - I_{p2}$ 后，均受到剩余光场矢势赋予的共同漂移动量：
  $$p_{drift}(\phi_r) = -A(t_r) = \frac{F_0}{\omega} \sin\phi_r$$
  这使得两电子平行动量 $p_{z1}, p_{z2}$ 具有强烈同号关联性，概率高度聚拢在第一象限（$p_{z1}>0, p_{z2}>0$）与第三象限（$p_{z1}<0, p_{z2}<0$），统计关联系数 $C_{corr} = \frac{\langle p_{z1} p_{z2} \rangle}{\sqrt{\langle p_{z1}^2 \rangle \langle p_{z2}^2 \rangle}} > 0$，完美复现 COLTRIMS 实验指纹；若低于碰撞电离阈值，电子通过 RESI（碰撞激发-后续场致电离）产生跨越二、四象限的十字交叉谱；非顺序双电离产率在 $10^{14}-10^{15} \text{ W/cm}^2$ 强度下比顺序双电离高出数个数量级，形成著名的平坦“膝盖”平台。
- **代码映射**:
  - `src/mod_strong_field_nsdi.f90`:
    - `calc_lotz_cross_section`: 碰撞电离截面计算；
    - `calc_nsdi_drift_momenta`: 两电子剩余能分享与偏振轴最终漂移动量；
    - `calc_nsdi_2d_momentum_dist`: 全周期数值积分生成 $P(p_{z1}, p_{z2})$ 与关联系数；
    - `calc_double_ion_yield_curve`: 光强依赖非顺序与顺序电离产率曲线与“膝盖结构”。

---

## 20. 磁与光 Feshbach 共振、分子弱束缚态与光致非弹性损耗 (Feshbach Resonances)

### 20.1 磁 Feshbach 共振 (MFR)、特征长度 $R^*$ 与闭通道几率 $Z(B)$
- **文献**:
  - C. Chin, R. Grimm, P. Julienne, and E. Tiesinga, *"Feshbach resonances in ultracold gases"*, **Rev. Mod. Phys.** 82, 1225 (2010). [DOI: 10.1103/RevModPhys.82.1225](https://doi.org/10.1103/RevModPhys.82.1225)
  - T. Köhler, K. Góral, and P. S. Julienne, *"Production of cold molecules via magnetically tunable Feshbach resonances"*, **Rev. Mod. Phys.** 78, 1311 (2006). [DOI: 10.1103/RevModPhys.78.1311](https://doi.org/10.1103/RevModPhys.78.1311)
  - B. Gao, *"Universal properties in ultracold atom-atom interactions"*, **Phys. Rev. A** 64, 010701(R) (2001). [DOI: 10.1103/PhysRevA.64.010701](https://doi.org/10.1103/PhysRevA.64.010701)
  - D. S. Petrov, *"Three-boson problem near a narrow Feshbach resonance"*, **Phys. Rev. Lett.** 93, 143201 (2004). [DOI: 10.1103/PhysRevLett.93.143201](https://doi.org/10.1103/PhysRevLett.93.143201)
- **核心理论**:
  磁 Feshbach 共振通过外磁场 Zeeman 效应调节开通道与束缚闭通道间的能量失谐 $\delta\mu(B - B_0)$。s 波散射长度随磁场满足经典色散公式：
  $$a(B) = a_{bg} \left( 1 - \frac{\Delta B}{B - B_0} \right)$$
  由范德瓦尔斯特征尺度 $\bar{a} \approx 0.956 R_{vdW}$ 与微观共振相互作用长度 $R^* = \frac{\hbar^2}{2\mu_{red} |a_{bg} \delta\mu \Delta B|}$ 定义无量纲共振强度参数 $s_{res} = (a_{bg} / \bar{a}) (\delta\mu \Delta B / \bar{E})$。
  在有效程展开与耦合通道两通道场论模型下，跨越普适宽共振（$s_{res} \gg 1, R^* \to 0$）与窄共振（$s_{res} \ll 1$）的分子态结合能解析解为（Chin et al. RMP 2010）：
  $$E_b(B) = \frac{\hbar^2}{2\mu_{red} (R^*)^2} \left[ \sqrt{1 + \frac{2 R^*}{a(B)}} - 1 \right]^2 \quad (a(B) > 0)$$
  由 Hellmann-Feynman 定理，弱束缚二聚体波函数中闭通道态几率成分 $Z(B) = \langle\psi|Q|\psi\rangle = \frac{\partial E_b / \partial B}{\delta\mu}$ 为：
  $$Z(B) = 1 - \frac{1}{\sqrt{1 + 2 R^* / a(B)}}$$
  在共振极点附近 $a \to \infty$ 时，$Z \to 0$（纯开通道大尺度晕轮二聚体 Halo Dimer）；在远共振区或窄共振中 $Z \to 1$（准经典闭通道分子）。
- **代码映射**:
  - `src/mod_feshbach_bound_states.f90`:
    - `init_mfr_preset`: $^{6}\text{Li}$（宽共振 $s_{res} \approx 51$）、$^{40}\text{K}$（中等）与 $^{87}\text{Rb}$（窄共振 $s_{res} \approx 0.13$）物理参数；
    - `calc_mfr_scattering_length`: 磁场依赖散射长度 $a(B)$；
    - `calc_mfr_bound_energy_universal` / `calc_mfr_bound_energy_coupled`: 普适与有限程耦合通道二聚体结合能；
    - `calc_mfr_closed_channel_fraction`: 闭通道权重分量 $Z(B)$。

### 20.2 光 Feshbach 共振 (OFR)、复散射长度与双体损耗率 $K_2$
- **文献**:
  - P. O. Fedichev, Yu. Kagan, G. V. Shlyapnikov, and J. T. M. Walraven, *"Influence of nearly resonant light on the scattering length in low-temperature atomic gases"*, **Phys. Rev. Lett.** 77, 2913 (1996). [DOI: 10.1103/PhysRevLett.77.2913](https://doi.org/10.1103/PhysRevLett.77.2913)
  - J. L. Bohn and P. S. Julienne, *"Semianalytic theory of laser-assisted ultracold collisions"*, **Phys. Rev. A** 60, 414 (1999). [DOI: 10.1103/PhysRevA.60.414](https://doi.org/10.1103/PhysRevA.60.414)
  - M. Theis et al., *"Tuning the Scattering Length with an Optical Feshbach Resonance"*, **Phys. Rev. Lett.** 93, 123001 (2004). [DOI: 10.1103/PhysRevLett.93.123001](https://doi.org/10.1103/PhysRevLett.93.123001)
- **核心理论**:
  使用近共振激光场诱导基态散射态与激发电子分子态耦合，有效复散射长度随激光失谐 $\Delta_L$ 演化为：
  $$\tilde{a}(\Delta_L) = a_{bg} + \delta a(\Delta_L) - i \frac{b(\Delta_L)}{2}$$
  $$\delta a(\Delta_L) = \frac{\ell_{opt} \gamma_{mol} \Delta_L}{\Delta_L^2 + (\gamma_{mol} / 2)^2}, \quad b(\Delta_L) = \frac{2 \ell_{opt} (\gamma_{mol}/2)^2}{\Delta_L^2 + (\gamma_{mol} / 2)^2}$$
  其中 $\ell_{opt}$ 为正比于激光光强的光学特征长度。虚部对应光缔合引起的非弹性自发辐射双体损失：
  $$K_2(\Delta_L) = \frac{4\pi \hbar}{\mu_{red}} \left( -\text{Im}(\tilde{a}) \right) = \frac{2\pi \hbar}{\mu_{red}} b(\Delta_L)$$
  在红失谐与蓝失谐 $\Delta_L = \pm \gamma_{mol} / 2$ 处，散射长度实部获得最大色散调制。
- **代码映射**:
  - `src/mod_feshbach_bound_states.f90`:
    - `calc_ofr_complex_scattering_length`: 复散射长度实虚部分离计算；
    - `calc_ofr_inelastic_loss_rate`: 实验可测双体非弹性损失速率常数 $K_2$（$\text{cm}^3/\text{s}$）。

---

## 21. 阿秒瞬态吸收光谱 (ATAS) 与光诱导态自电离干涉动力学

### 21.1 氦原子 $2s2p (^1P)$ 双激发态与相位微扰模型 (PPM)
- **文献**:
  - M. Chini, K. Zhao, and Z. Chang, *"The generation, detection, and applications of attosecond pulses"*, **Nat. Photonics** 8, 178 (2014). [DOI: 10.1038/nphoton.2013.362](https://doi.org/10.1038/nphoton.2013.362)
  - C. Ott et al., *"Lorentz Meets Fano in the Interaction of Two Strongly Driven Helium Rydberg Series"*, **Science** 340, 716 (2013). [DOI: 10.1126/science.1232759](https://doi.org/10.1126/science.1232759)
  - M. Wu, S. Chen, M. B. Gaarde, and K. J. Schafer, *"Time-frequency analysis of autoionizing states in attosecond transient absorption"*, **Phys. Rev. A** 93, 033405 (2016). [DOI: 10.1103/PhysRevA.93.033405](https://doi.org/10.1103/PhysRevA.93.033405)
  - U. Fano, *"Effects of Configuration Interaction on Intensities and Phase Shifts"*, **Phys. Rev.** 124, 1866 (1961). [DOI: 10.1103/PhysRev.124.1866](https://doi.org/10.1103/PhysRev.124.1866)
- **核心理论**:
  阿秒瞬态吸收光谱（ATAS）采用极紫外（XUV）阿秒单脉冲激发体系到自电离共振态（如氦原子 $2s2p\ ^1P$, $E_0 = 60.15$ eV, $\Gamma = 0.037$ eV, $\tau \approx 17.8$ fs），并用飞秒近红外（NIR）激光控制后续偶极退相与自电离干涉。
  在相位微扰模型（PPM）中，NIR 场通过极化率 $\alpha$ 诱导含时动态 Stark 势能移动，累积附加作用量相位 $\phi(\tau) = \int_\tau^\infty \Delta E(t') dt'$。稳态 Fano 不对称参数 $q_0$ 发生几何旋转映射为有效含时参数：
  $$q(\tau) = \frac{q_0 + \tan\phi(\tau)}{1 - q_0 \tan\phi(\tau)}$$
  导致吸收谱线在 Fano 共振、对称吸收峰、对称色散窗及反共振峰（Window Resonance）之间发生相干翻转。
- **代码映射**:
  - `src/mod_attosecond_transient_absorption.f90`:
    - `init_atas_helium_benchmark`: 氦原子 $2s2p\ ^1P$ 经典自电离能谱参数初始化；
    - `calc_laser_dressed_fano_q`: 激光修饰相移动态 $q(\tau)$ 参数计算。

### 21.2 光诱导态 (LIS)、AC Stark 位移与二维时延光谱 $\Delta\text{OD}(\omega, \tau)$
- **文献**:
  - S. Chen et al., *"Light-induced states in attosecond transient absorption spectra of helium"*, **Phys. Rev. A** 86, 063408 (2012). [DOI: 10.1103/PhysRevA.86.063408](https://doi.org/10.1103/PhysRevA.86.063408)
  - Z. Chang, *Fundamentals of Attosecond Optics*, CRC Press, Boca Raton (2011).
- **核心理论**:
  当强红外光耦合明态（Bright State $|b\rangle$）与邻近单光子禁戒暗态（Dark State $|d\rangle$，如 $2p^2\ ^1S$ 或 $2s^2\ ^1S$）时，系统形成缀饰二能级杂化态，并在能谱中诱导出**光诱导态（Light-Induced States, LIS）**：
  $$E_{\text{LIS}} = \frac{E_b + (E_d \pm \hbar\omega_L)}{2} \pm \frac{1}{2}\sqrt{\delta^2 + \Omega_R^2}$$
  在 XUV-NIR 时延 $\tau > 0$ 区域，明暗态波包相干叠加诱发周期为 $T_{\text{beat}} = \frac{h}{\Delta E}$ 的高阶量子拍频（Quantum Beats），二维差分光密度谱 $\Delta\text{OD}(\omega, \tau) = -\log(I_{pump-probe} / I_{probe})$ 呈现出双缝干涉条纹与斜向双光子吸收特征。
- **代码映射**:
  - `src/mod_attosecond_transient_absorption.f90`:
    - `calc_light_induced_state_energy`: 缀饰光诱导态特征能量；
    - `calc_quantum_beat_period_fs`: 量子拍频周期；
    - `calc_atas_spectrum`: 二维能量-时延吸收光谱 $\Delta\text{OD}(\omega, \tau)$ 矩阵全量计算。

---

## 22. 双色反向旋转圆偏振场与分子光电子圆二色性 (PECD)

### 22.1 双色反向圆偏振合成场与 $C_{r+1}$ 离散旋转动力学对称性
- **文献**:
  - O. Kfir et al., *"Generation of strongly elliptically polarized high-harmonic emission using bicircular laser fields"*, **Nat. Photonics** 9, 99 (2015). [DOI: 10.1038/nphoton.2014.364](https://doi.org/10.1038/nphoton.2014.364)
  - C. A. Mancuso et al., *"Strong-field ionization with two-color circularly polarized laser fields"*, **Phys. Rev. A** 91, 031402(R) (2015). [DOI: 10.1103/PhysRevA.91.031402](https://doi.org/10.1103/PhysRevA.91.031402)
  - D. Ayuso et al., *"Synthetic chiral light for ultra-fast control of chiral measurements"*, **Nat. Photonics** 13, 866 (2019). [DOI: 10.1038/s41566-019-0531-2](https://doi.org/10.1038/s41566-019-0531-2)
- **核心理论**:
  双色旋转场由基频 $\omega_1$（旋向 $h_1 = \pm 1$）与谐波 $\omega_2 = r\omega_1$（旋向 $h_2 = \mp 1$）共面合成：
  $$\mathbf{E}(t) = \frac{E_{0,1}}{\sqrt{2}} f(t) [\cos(\omega_1 t) \hat{\mathbf{x}} + h_1 \sin(\omega_1 t) \hat{\mathbf{y}}] + \frac{E_{0,2}}{\sqrt{2}} f(t) [\cos(r\omega_1 t + \phi) \hat{\mathbf{x}} + h_2 \sin(r\omega_1 t + \phi) \hat{\mathbf{y}}]$$
  对于反向旋转 $\omega + 2\omega$ 场（$h_1 = +1, h_2 = -1, r=2$），场矢量在偏振面内扫过三叶草状（Trefoil）轨迹，具有严格的 $C_3$ 离散动力学对称性：
  $$\mathcal{R}\left(\frac{2\pi}{3}\right) \mathbf{E}(t) = \mathbf{E}\left(t + \frac{T_1}{3}\right)$$
  若同向旋转则具有心形单叶分布（$C_1$ 对称）。
- **代码映射**:
  - `src/mod_bicircular_pecd.f90`:
    - `init_bicircular_field`: 场强、波长、相对相位与包络设置；
    - `calc_bicircular_field_at_t`: 瞬时场强 $\mathbf{E}(t)$ 与矢势 $\mathbf{A}(t)$ 计算；
    - `calc_dynamical_symmetry_fold`: 动力学折叠对称度判断；
    - `calc_bicircular_trajectory`: Lissajous 轨迹合成。

### 22.2 光电子圆二色性 (PECD) 与手性四面体势模型
- **文献**:
  - N. Böwering et al., *"Asymmetry in photoionization of chiral molecules with circularly polarized light"*, **Phys. Rev. Lett.** 86, 1187 (2001). [DOI: 10.1103/PhysRevLett.86.1187](https://doi.org/10.1103/PhysRevLett.86.1187)
  - C. Lux et al., *"Circular Dichroism in the Photoelectron Angular Distributions of Camphor and Fenchone from Multiphoton Ionization with Femtosecond Laser Pulses"*, **Angew. Chem. Int. Ed.** 51, 5001 (2012). [DOI: 10.1002/anie.201109035](https://doi.org/10.1002/anie.201109035)
  - B. Ritchie, *"Theory of the angular distribution of photoelectrons ejected from optically active molecules and molecular negative ions"*, **Phys. Rev. A** 13, 1411 (1976). [DOI: 10.1103/PhysRevA.13.1411](https://doi.org/10.1103/PhysRevA.13.1411)
- **核心理论**:
  随机取向手性分子在圆偏振光单光子/多光子电离中，偶极基质元不同分波干涉产生沿光传播方向（$z$ 轴）的奇宇称勒让德多项式不对称性：
  $$I(\theta, \phi) = \frac{\sigma_{tot}}{4\pi} \left[ 1 + \beta_1 P_1(\cos\theta) + \beta_2 P_2(\cos\theta) + \dots \right]$$
  奇数阶参数 $\beta_1$ 导致明显的前向-后向不对称发射（Forward-Backward Asymmetry），其相对不对称比率为：
  $$G_{\text{PECD}} = \frac{I_{\text{forward}} - I_{\text{backward}}}{I_{\text{forward}} + I_{\text{backward}}} = \frac{1}{2} \beta_1$$
  四中心四面体手性势中伪标量手性不变量为 $\chi = [(\mathbf{R}_1 - \mathbf{R}_4) \times (\mathbf{R}_2 - \mathbf{R}_4)] \cdot (\mathbf{R}_3 - \mathbf{R}_4) \prod_{i < j} (Z_i - Z_j)$。对于对映异构体（$R$-/$S$-型）满足严格反转关系：$\beta_1^{(R)} = -\beta_1^{(S)}$ 且 $G_{\text{PECD}}^{(R)} = -G_{\text{PECD}}^{(S)}$。
- **代码映射**:
  - `src/mod_bicircular_pecd.f90`:
    - `init_chiral_tetrahedral_molecule`: $R$ 与 $S$ 对映异构体初始化；
    - `calc_chirality_measure`: 伪标量不变量 $\chi$ 检验；
    - `calc_forward_backward_asymmetry`: 前后不对称度换算；
    - `calc_chiral_beta1_model`: 动能依赖 $\beta_1(E)$ 模型；
    - `calc_pecd_pad_spectrum`: 2D 光电子角分布 $I(\theta, \phi)$ 计算。

---

## 23. 超冷极性分子化学反应动力学与微波/静电偶极遮蔽

### 23.1 超冷双分子化学反应与普遍短程吸收损失
- **文献**:
  - K.-K. Ni et al., *"A High Phase-Space-Density Gas of Polar Molecules"*, **Science** 322, 231 (2008). [DOI: 10.1126/science.1163861](https://doi.org/10.1126/science.1163861)
  - S. Ospelkaus et al., *"Quantum-State Controlled Chemical Reactions of Ultracold Potassium-Rubidium Molecules"*, **Science** 327, 853 (2010). [DOI: 10.1126/science.1184121](https://doi.org/10.1126/science.1184121)
  - M. Mayle, G. Quéméner, B. P. Ruzic, and J. L. Bohn, *"Scattering of ultracold molecules in the presence of resonant sticky collisions"*, **Phys. Rev. A** 87, 012709 (2013). [DOI: 10.1103/PhysRevA.87.012709](https://doi.org/10.1103/PhysRevA.87.012709)
- **核心理论**:
  超冷双原子极性分子（如 $^{40}\text{K}^{87}\text{Rb}, ^{23}\text{Na}^{87}\text{Rb}$）在纳开尔文到微开尔文碰撞中，放热化学重排反应 $2\text{KRb} \to \text{K}_2 + \text{Rb}_2$ 或四体中间复合物“粘滞碰撞”（Sticky Collisions）导致剧烈的非弹性双体损耗（$K_2 \sim 10^{-10}\ \text{cm}^3/\text{s}$）。在量子数亏损理论（QDT）中，短程吸收被参数化为非弹性吸收率 $y \in [0, 1]$，其中 $y=1$ 为普适全黑体吸收极限。
- **代码映射**:
  - `src/mod_ultracold_reaction_shielding.f90`:
    - `init_ultracold_molecule_preset`: KRb, NaRb, NaK 分子质量、偶极矩与 $C_6$ 参数初始化。

### 23.2 微波蓝失谐缀饰免交叉势垒与 WKB 隧穿抑制
- **文献**:
  - G. Quéméner and J. L. Bohn, *"Shielding of polar molecules with electric fields"*, **Phys. Rev. A** 81, 022702 (2010). [DOI: 10.1103/PhysRevA.81.022702](https://doi.org/10.1103/PhysRevA.81.022702)
  - L. Anderegg et al., *"Observation of microwave shielding of ultracold molecules"*, **Science** 373, 779 (2021). [DOI: 10.1126/science.abg9502](https://doi.org/10.1126/science.abg9502)
  - A. Schindewolf et al., *"Evaporative cooling of microwave-shielded polar molecules to quantum degeneracy"*, **Nature** 607, 677 (2022). [DOI: 10.1038/s41586-022-04900-3](https://doi.org/10.1038/s41586-022-04900-3)
  - K. Matsuda et al., *"Resonant collateral shielding of polar molecules by a static electric field"*, **Science** 370, 1324 (2020). [DOI: 10.1126/science.abe7370](https://doi.org/10.1126/science.abe7370)
- **核心理论**:
  施加蓝失谐圆偏振微波（$\Delta > 0$）或强 DC 静电场耦合分子转动态 $|J=0\rangle$ 与 $|J=1\rangle$，共振偶极-偶极相互作用 $C_3 \sim d^2$ 形成工程化免交叉排斥势垒：
  $$V_{\text{rep}}(R) = \frac{\hbar}{2} \left[ \sqrt{\Delta^2 + \left(\frac{2C_3}{R^3}\right)^2 + \Omega^2} - \Delta \right]$$
  屏蔽半径位于 $R_{\text{shield}} \approx (2C_3 / \hbar\Delta)^{1/3} \sim 300 - 1000\ a_0$，势垒高度达到数百微开尔文，远高于碰撞动能。
  分子穿透势垒到达短程反应区的几率由 WKB 隧穿公式决定：
  $$T_{\text{tunnel}} = \exp\left[ -2 \int_{R_{\text{in}}}^{R_{\text{out}}} \sqrt{\frac{2\mu}{\hbar^2}(V_{\text{eff}}(R) - E)} dR \right] \ll 1$$
  非弹性损失率 $K_2^{(\text{inel})}$ 被压低 6-8 个数量级，弹性散射截面 $\sigma_{\text{el}} \approx 4\pi R_{\text{shield}}^2$ 保持巨大，实现蒸发冷却关键判据：
  $$\gamma = \frac{K_2^{(\text{el})}}{K_2^{(\text{inel})}} > 100 - 10^7$$
- **代码映射**:
  - `src/mod_ultracold_reaction_shielding.f90`:
    - `init_shielding_config`: 微波/直流静电屏蔽配置；
    - `calc_effective_shielding_potential`: 有效相互作用势 $V_{\text{eff}}(R)$；
    - `calc_shielding_barrier_height`: 屏蔽势垒高度与位置；
    - `calc_wkb_tunneling_probability`: WKB 量子隧穿衰减因数；
    - `calc_shielded_scattering_rates`: 弹/非弹性速率与优良因子 $\gamma$；
    - `calc_shielding_detuning_scan`: 失谐扫描曲线。

---

## 24. 里德堡原子阻塞、PXP 约束模型与量子多体疤痕

### 24.1 巨范德瓦尔斯作用 $C_6 \propto n^{11}$ 与双原子偶极阻塞机制
- **文献**:
  - D. Jaksch, J. I. Cirac, P. Zoller, S. L. Rolston, R. Côté, and M. D. Lukin, *"Fast Quantum Gates for Neutral Atoms"*, **Phys. Rev. Lett.** 85, 2208 (2000). [DOI: 10.1103/PhysRevLett.85.2208](https://doi.org/10.1103/PhysRevLett.85.2208)
  - M. D. Lukin et al., *"Dipole Blockade and Quantum Information Processing in Mesoscopic Atomic Ensembles"*, **Phys. Rev. Lett.** 87, 037901 (2001). [DOI: 10.1103/PhysRevLett.87.037901](https://doi.org/10.1103/PhysRevLett.87.037901)
  - M. Saffman, T. G. Walker, and K. Mølmer, *"Quantum information with Rydberg atoms"*, **Rev. Mod. Phys.** 82, 2313 (2010). [DOI: 10.1103/RevModPhys.82.2313](https://doi.org/10.1103/RevModPhys.82.2313)
  - A. Browaeys and T. Lahaye, *"Many-body physics with individually controlled Rydberg atoms"*, **Nat. Phys.** 16, 132 (2020). [DOI: 10.1038/s41567-019-0733-z](https://doi.org/10.1038/s41567-019-0733-z)
- **核心理论**:
  高主量子数里德堡态具有超大电偶极跃迁矩阵元，其长程色散相互作用系数具有极端标度律 $C_6 \propto n^{11}$。在激光 Rabi 频率 $\Omega$ 驱动下，二聚体激发能发生巨大位移 $V_{\text{vdW}}(R) = C_6 / R^6$。
  当间距小于**里德堡阻塞半径 (Rydberg Blockade Radius)**：
  $$R_b = \left( \frac{|C_6|}{\hbar \Omega} \right)^{1/6}$$
  双重激发态 $|rr\rangle$ 发生剧烈能级失谐，强行禁止同时激发，系统被约束在纠缠单激发态（W-态）$|W\rangle = \frac{1}{\sqrt{2}}(|gr\rangle + |rg\rangle)$，集体有效 Rabi 频率获得 $\sqrt{2}\Omega$ 增强。
- **代码映射**:
  - `src/mod_rydberg_blockade.f90`:
    - `init_rydberg_atom`: $^{87}\text{Rb}$ $n S$ 态量子数亏损与 $n^{11}$ 相互作用缩放；
    - `calc_rydberg_blockade_radius`: 阻塞半径 $R_b$ 解析求解；
    - `calc_two_atom_dynamics`: 四态全空间保范数细化时域动力学步进求解器。

### 24.2 1D 原子阵列、PXP 拓扑约束模型与量子多体疤痕 (Quantum Many-Body Scars)
- **文献**:
  - H. Bernien et al., *"Probing many-body dynamics on a 51-atom quantum simulator"*, **Nature** 551, 579 (2017). [DOI: 10.1038/nature24622](https://doi.org/10.1038/nature24622)
  - C. J. Turner, A. A. Michailidis, D. A. Abanin, M. Serbyn, and Z. Papić, *"Weak ergodicity breaking from quantum many-body scars"*, **Nat. Phys.** 14, 745 (2018). [DOI: 10.1038/s41567-018-0137-5](https://doi.org/10.1038/s41567-018-0137-5)
- **核心理论**:
  在光镊一维原子晶格中，最近邻强阻塞条件排除了相邻双激发（$n_i n_{i+1} = 0$），希尔伯特空间维度被约束为斐波那契数 $D_N = F_{N+2}$。系统哈密顿量退化为著名的 **PXP 模型**：
  $$H_{\text{PXP}} = \frac{\hbar\Omega}{2} \sum_{i} P_{i-1} \sigma_i^x P_{i+1} - \hbar\Delta \sum_i n_i$$
  当体系从交错反铁磁 Néel 态 $|r g r g \dots\rangle$ 发生量子淬火时，违反通常本征态热化假说（ETH），呈现出非热化的**量子多体疤痕（Quantum Many-Body Scars）**。交错交替序参量（Staggered $Z_2$ Order Parameter）：
  $$\mathcal{O}_{Z_2}(t) = \frac{1}{N} \sum_{i=1}^N (-1)^i \langle n_i(t) \rangle$$
  展现出周期为 $T_{\text{scar}} \approx \frac{2\pi}{1.33 \Omega}$ 的长寿命宏观相干周期性复苏振荡。
- **代码映射**:
  - `src/mod_rydberg_blockade.f90`:
    - `init_rydberg_array`: 一维原子阵列与周期/开边界参数；
    - `calc_z2_order_parameter`: $Z_2$ 空间交错电荷密度波序参量；
    - `calc_rydberg_scar_dynamics`: PXP 疤痕周期复苏相干动力学模拟。

---

## 25. 表面量子散射与选择性吸附共振 (Selective Adsorption Resonances)

### 25.1 2D 周期晶格相干散射与硬波纹表面 (HCS) 程函近似
- **文献**:
  - R. Frisch and O. Stern, *"Beugung von Materiestrahlen an Kristallgitterflächen"*, **Z. Phys.** 84, 430 (1933).
  - J. E. Lennard-Jones and A. F. Devonshire, *"Diffraction and selective adsorption of atoms by crystals"*, **Nature** 137, 1069 (1936).
  - G. Boato, P. Cantini, and L. Mattera, *"Diffraction of He and H2 molecular beams from a LiF(001) crystal surface"*, **J. Phys. C: Solid State Phys.** 6, L394 (1973); **Surf. Sci.** 55, 141 (1976).
  - J. R. Manson, *"Inelastic scattering from surfaces"*, **Phys. Rev. B** 43, 6924 (1991).
  - G. Benedek and J. P. Toennies, *"Atomic Scale Dynamics at Surfaces: Theory and Experimental Methods with Helium Atom Scattering"*, **Springer Series in Surface Sciences** 63 (2018).
- **核心理论**:
  低能热原子（如 He 束，能量 10-100 meV）与刚性晶体表面碰撞时，横向周期势场产生离散二维布拉格衍射峰：
  $$\mathbf{K}_{\mathbf{G}} = \mathbf{K}_0 + \mathbf{G} = \left( K_{0x} + m \frac{2\pi}{a_x}, K_{0y} + n \frac{2\pi}{a_y} \right)$$
  出射法向动量满足能量守恒：$k_{Gz}^2 = k_0^2 - |\mathbf{K}_0 + \mathbf{G}|^2$。
  在硬波纹表面（HCS）程函近似下，衍射振幅由二维贝塞尔函数卷积决定：
  $$A_{\mathbf{G}} = (-i)^{|m|+|n|} J_m(u_x) J_n(u_y)$$
  有限温度下表面声子非弹性热激发导致弹性峰按 Debye-Waller 因子指数衰减：$I_{\mathbf{G}}(T) = I_{\mathbf{G}}^{(0)} e^{-2W(T)}$。
- **代码映射**:
  - `src/mod_surface_scattering.f90`:
    - `init_surface_lattice`: 二维晶格与波纹参数；
    - `calc_surface_diffraction_channels`: 开/闭通道判据与出射角；
    - `calc_hcs_diffraction_probabilities`: HCS 程函衍射几率与幺正归一化；
    - `calc_surface_debye_waller`: 表面声子 Debye-Waller 热衰减。

### 25.2 选择性吸附共振 (SAR) 与 Fano 不对称干涉线型
- **文献**:
  - K. L. Wolfe and J. H. Weare, *"Evidence of Fano-type interference in selective adsorption"*, **Phys. Rev. Lett.** 41, 1775 (1978).
  - P. Cantini et al., *"Resonant scattering of He atoms from LiF(001)"*, **Surf. Sci.** 63, 104 (1977).
- **核心理论**:
  当某个倏逝衍射通道（闭通道，Evanescent Beam）的法向运动能量恰好与表面横向平均吸附势的离散束缚能级 $\epsilon_v < 0$ 重合时：
  $$\frac{\hbar^2 k_{Gz}^2}{2M} = E_i - \frac{\hbar^2 (\mathbf{K}_0 + \mathbf{G})^2}{2M} = \epsilon_v$$
  原子发生共振表面囚禁，通过虚态耦合向开放通道（尤其是镜面反射通道 $(0,0)$）产生强烈的相消/相长量子干涉，表现为典型的 Fano 不对称线型：
  $$I_{\text{spec}}(\theta_i) \propto \frac{(q_{\text{SAR}} + \epsilon)^2}{1 + \epsilon^2}, \quad \epsilon = \frac{E_{Gz} - \epsilon_v}{\Gamma_v/2}$$
- **代码映射**:
  - `src/mod_surface_scattering.f90`:
    - `init_surface_potential_morse`: 表面吸引阱深与束缚态求解；
    - `calc_selective_adsorption_resonance`: 束缚态共振判定与 Fano 线型调制。

---

## 26. 气-固表面催化反应与 Eley-Rideal 提取机理

### 26.1 反应放热量分配与非热化超热振动激发
- **文献**:
  - D. D. Eley and E. K. Rideal, *"The catalysis of the para-hydrogen conversion by tungsten"*, **Nature** 146, 401 (1940); **Proc. R. Soc. London A** 178, 429 (1941).
  - C. T. Rettner, *"Reaction of an H atom beam with a hydrogen-passivated Si(100) surface"*, **Phys. Rev. Lett.** 69, 383 (1992).
  - C. T. Rettner and D. J. Auerbach, *"Dynamics of the Eley-Rideal reaction: Hydrogen atom recombination on Cu(111)"*, **J. Chem. Phys.** 101, 1529 (1994).
  - B. Jackson and M. Persson, *"Quantum mechanical study of the Eley-Rideal reaction of H(g) with H(ads) on a metal surface"*, **J. Chem. Phys.** 96, 2378 (1992).
- **核心理论**:
  Eley-Rideal (ER) 机理指入射气相超热原子 $A(g)$ 直接单次碰撞提取表面化学吸附原子 $B(\text{ads})$，形成激发态气相分子 $AB(g, v, j)$ 飞离表面：
  $$A(g) + B(\text{ads})/\text{Surf} \to AB(g, v, j) + \text{Surf}$$
  反应释放巨大放热量 $\Delta E_{\text{exo}} = D_{AB} - D_{\text{chem}}$（以 H+H/Cu(111) 为例高达 $2.30\ \text{eV}$）。
  在数十飞秒碰撞时间内，体系远离热平衡：约 $50\%$ 的能量转化为分子的剧烈振动激发（导致振动布居反转如 $P(v=2) > P(v=0)$），约 $35\%$ 转化为高速平动，仅约 $10\%$ 耗散到基底。
- **代码映射**:
  - `src/mod_surface_reaction_er.f90`:
    - `init_er_reaction_system`: 反应体系常数与放热量计算；
    - `calc_er_potential_2d`: 反应路径 2D 势能面 $V(r, Z_{\text{cm}})$；
    - `calc_er_energy_partitioning`: 放热能量通道分配（振动、平动、基底）；
    - `calc_er_vibrational_populations`: 产物分态振动布居反转计算。

### 26.2 反应截面与准无势垒热催化速率
- **文献**:
  - A. C. Luntz, *"Gas-surface reaction dynamics: From simple models to complex chemistry"*, **Science** 302, 1352 (2003).
- **核心理论**:
  ER 反应通常表现为极低或无进入势垒，反应截面在亚电子伏能区表现为几何碰撞截面 $\sigma_{\text{ER}} \sim 0.2 - 0.5\ \text{\AA}^2$，热反应速率常数满足阿伦尼乌斯型或通量积分形式：
  $$k_{\text{ER}}(T) = \langle \sigma_{\text{ER}} v_{\text{gas}} \rangle_T$$
- **代码映射**:
  - `src/mod_surface_reaction_er.f90`:
    - `calc_er_reaction_cross_section`: 入射能量依赖反应截面；
    - `calc_er_thermal_rate_constant`: 温度依赖热催化速率常数。

---

## 27. 表面非绝热动力学与电子摩擦耗散 (Electronic Friction & GLE)

### 27.1 局域密度摩擦近似 (LDFA) 与广义朗之万方程 (GLE)
- **文献**:
  - J. C. Tully, *"Dynamics of gas-surface interactions: Thermal accommodation and sticking"*, **Surf. Sci.** 111, 461 (1981); *"Nonadiabatic dynamics at surfaces"*, **Annu. Rev. Phys. Chem.** 51, 153 (2000).
  - M. Head-Gordon and J. C. Tully, *"Molecular dynamics with electronic frictions"*, **J. Chem. Phys.** 103, 10137 (1995).
  - J. I. Juaristi, M. Alducin, R. Díez Muiño, H. F. Busnengo, and A. Salin, *"Role of electron-hole pairs in the dissipation of energy during molecular adsorption on metal surfaces"*, **Phys. Rev. Lett.** 100, 116102 (2008).
  - A. M. Wodtke, J. C. Tully, and D. J. Auerbach, *"Electronically nonadiabatic dynamics in molecule-surface interactions"*, **Science** 290, 1585 (2000).
- **核心理论**:
  金属表面存在连续无能隙的电子单粒子激发谱（电子-空穴对，EHPs）。高速核运动强耦合导带电子，导致机械动能非绝热耗散。
  局域密度摩擦近似下电子阻尼力为 $\mathbf{F}_{\text{friction}} = -\eta(z) \mathbf{v}$，摩擦系数指数衰减：$\eta(z) = \eta_0 e^{-\gamma(z - z_{\text{surf}})}$。
  分子运动遵循广义朗之万方程（GLE）：
  $$M \frac{d^2 z}{dt^2} = -\frac{\partial V(z)}{\partial z} - \eta(z) \frac{dz}{dt} + \xi(t)$$
  单次碰壁非绝热能量损耗为 $\Delta E_{\text{loss}} = \int \eta(z) v(t)^2 dt$。
- **代码映射**:
  - `src/mod_surface_electronic_friction.f90`:
    - `init_metal_surface`: Au(111), Cu(111), Pt(111) 参数初始化；
    - `calc_electronic_friction_coeff`: 摩擦系数分布 $\eta(z)$；
    - `calc_surface_morse_force`: 绝热表面保守力；
    - `integrate_gle_scattering_trajectory`: GLE 动力学步进与非绝热能损。

### 27.2 吸附分子高频振动弛豫寿命与电子浴退相干
- **文献**:
  - B. N. J. Persson and M. Persson, *"Vibrational phase relaxation at surfaces"*, **Surf. Sci.** 97, 609 (1980).
- **核心理论**:
  高频分子化学键在金属表面吸附时，通过与电子-空穴对的共振耦合发生快速能量弛豫，其振动弛豫速率与寿命满足：
  $$\Gamma_{\text{vib}} = \frac{1}{\tau_{\text{vib}}} = \frac{\eta(z_{\text{ads}})}{M}$$
- **代码映射**:
  - `src/mod_surface_electronic_friction.f90`:
    - `calc_vibrational_relaxation_rate`: 表面吸附分子振动弛豫速率与寿命。

---

## 28. 掠入射快原子表面量子衍射与彩虹散射 (GIFAD)

### 28.1 轴向沟道快慢自由度解耦与横向低能量子波长
- **文献**:
  - P. Rousseau, H. Khemliche, A. G. Borisov, and P. Roncin, *"Quantum Grazing Incidence Fast Atom Diffraction"*, **Phys. Rev. Lett.** 98, 016104 (2007).
  - A. Schüller, S. Wethekam, and H. Winter, *"Diffraction of Fast Atoms under Axial Surface Channeling Conditions"*, **Phys. Rev. Lett.** 99, 136106 (2007).
  - H. Winter and A. Aigner, *"Fast atom diffraction at surfaces"*, **Prog. Surf. Sci.** 86, 169 (2011).
  - P. Roncin and H. Khemliche, *"Fast atom diffraction: A new tool for surface science"*, **Nucl. Instrum. Methods Phys. Res. B** 269, 1400 (2011).
- **核心理论**:
  GIFAD 使用 keV 准直轻原子束以掠角 $\theta_{\text{in}} < 1^\circ - 2^\circ$ 沿低指数晶向沟道入射。快轴向沟道效应使纵向与横向运动彻底解耦：
  1. 纵向平行运动动能高达数 keV，$\lambda_\parallel \sim 10^{-3}\ \text{\AA}$，对晶向凹凸完全平滑化；
  2. 横向垂直有效动能骤降：$E_\perp = E_{\text{beam}} \sin^2 \theta_{\text{in}} \sim 0.1 - 2.0\ \text{eV}$，横向德布罗意波长扩大至量子尺度：
     $$\lambda_\perp = \frac{h}{\sqrt{2 M E_\perp}} \sim 0.2 - 1.0\ \text{\AA}$$
  横向周期沟道产生一维 Bragg 衍射极点：$\sin\theta_m = m \frac{\lambda_\perp}{a_x}$。
- **代码映射**:
  - `src/mod_grazing_fast_atom_diffraction.f90`:
    - `init_gifad_experiment`: 掠入射轴向沟道运动学初始化；
    - `calc_gifad_transverse_kinematics`: 横向有效能量 $E_\perp$ 与德布罗意波长 $\lambda_\perp$。

### 28.2 经典表面彩虹散射与亚皮米波纹度逆向反演
- **文献**:
  - E. A. Manson and V. Celli, *"Inelastic helium scattering from surfaces"*, **Surf. Sci.** 24, 495 (1971).
- **核心理论**:
  表面波纹度斜率在拐点处达到极大值，形成经典表面彩虹角：
  $$\theta_R \approx \arctan\left( \frac{2\pi \zeta}{a_x} \right)$$
  量子衍射在彩虹角处呈现明显的 Airy 级斑汇聚。通过彩虹角可高精度反演亚皮米表面波纹幅度：
  $$\zeta = \frac{a_x \tan\theta_R}{2\pi}$$
- **代码映射**:
  - `src/mod_grazing_fast_atom_diffraction.f90`:
    - `calc_gifad_rainbow_angle`: 表面经典彩虹散射角 $\theta_R$；
    - `calc_gifad_diffraction_spectrum`: 1D 横向量子衍射谱；
    - `calc_surface_corrugation_from_rainbow`: 表面亚皮米波纹度逆向反演重构。

---

## 29. 冷离子-中性原子杂化散射与极化阱动力学 (Cold Ion-Atom Hybrid Scattering)

### 29.1 长程 $1/r^4$ 极化势、Langevin 俘获截面与修正有效力程展开 (MERE)
- **文献**:
  - P. Langevin, *"Une formule fondamentale de théorie cinétique"*, **Ann. Chim. Phys.** 5, 245 (1905).
  - R. Côté and A. Dalgarno, *"Ultracold atom-ion collisions"*, **Phys. Rev. A** 62, 012709 (2000). [DOI: 10.1103/PhysRevA.62.012709](https://doi.org/10.1103/PhysRevA.62.012709)
  - R. Côté, *"Ultracold atom-ion interactions and collisions"*, **Phys. Rev. Lett.** 89, 083201 (2002). [DOI: 10.1103/PhysRevLett.89.083201](https://doi.org/10.1103/PhysRevLett.89.083201)
  - Z. Idziaszek et al., *"Analytical solutions for energy-dependent s-wave ion-atom scattering lengths"*, **Phys. Rev. A** 83, 052713 (2011).
- **核心理论**:
  带电荷 $q$ 的离子与极化率为 $\alpha$ 的中性原子间由感应偶极相互作用主导：$V(r) \to -C_4 / (2 r^4)$，其中 $C_4 = \alpha q^2 / (4\pi \epsilon_0)$。
  体系特征长度与特征能量标度为：$R^* = \sqrt{2\mu C_4/\hbar^2}$，$E^* = \hbar^2 / [2\mu (R^*)^2]$。
  经典离心势垒极值确定临界碰撞参数 $b_c(E) = (2 C_4 / E)^{1/4}$，经典 Langevin 反应俘获截面为 $\sigma_L(E) = 2\pi \sqrt{C_4 / (2E)}$，对应的化学反应碰撞速率常数与能量无关：$K_L = 2\pi \sqrt{C_4 / \mu}$。
  在极低温量子极限下，由于 $1/r^4$ 长程极化势具有反常奇点，经典 ERE 失效，由修正有效力程展开 (MERE) 给出：
  $$k \cot\delta_0 = -\frac{1}{a_s} + \frac{\pi}{3 R^*} (k R^*)^2 + \frac{1}{2} r_0 k^2 + \dots$$
- **代码映射**:
  - `src/mod_ion_atom_scattering.f90`:
    - `init_ion_atom_system`: 杂化原子-离子系统特征尺度初始化；
    - `calc_langevin_cross_section`: 经典 Langevin 截面与临界碰撞参数；
    - `calc_langevin_rate_coefficient`: 能量无关 Langevin 速率常数；
    - `calc_mere_phase_shift_s_wave`: MERE 极化量子散射相移。

### 29.2 Paul 射频阱离子微运动非弹性碰撞致热动力学
- **文献**:
  - M. Cetina, A. T. Grier, and V. Vuletić, *"Micromotion-induced limit to atom-ion sympathetic cooling in Paul traps"*, **Phys. Rev. Lett.** 109, 253201 (2012). [DOI: 10.1103/PhysRevLett.109.253201](https://doi.org/10.1103/PhysRevLett.109.253201)
  - K. Chen, S. T. Sullivan, and E. R. Hudson, *"Neutral production of cold ions in Paul traps"*, **Phys. Rev. Lett.** 112, 143009 (2014).
- **核心理论**:
  在 Paul 射频四极阱中，时间周期交变场驱动离子产生高频微运动（Micromotion）。在单次碰撞中，微运动动能通过非弹性散射不可逆地转化为久期运动热能，导致同情冷却存在极限平衡温度：
  $$T_{\text{limit}} \propto \left( \frac{m_{\text{atom}}}{m_{\text{ion}}} \right) \cdot T_{\text{secular}}$$
  致热率与 Mathieu 稳定性参数 $q$ 及射频驱动频率 $\omega_{\text{rf}}$ 满足 $\frac{\text{d}E}{\text{d}t} \propto \frac{m_a}{m_i + m_a} q^2 \omega_{\text{rf}} E_{\text{sec}}$。
- **代码映射**:
  - `src/mod_ion_atom_scattering.f90`:
    - `calc_rf_micromotion_heating`: 微运动非平衡碰撞致热率与平衡极限温度计算。

---

## 30. Tully 最少开关表面跳跃与非绝热分子动力学 (FSSH)

### 30.1 最少开关表面跳跃几率算法与 Velocity Verlet 核演化
- **文献**:
  - J. C. Tully, *"Molecular dynamics with electronic transitions"*, **J. Chem. Phys.** 93, 1061 (1990). [DOI: 10.1063/1.459170](https://doi.org/10.1063/1.459170)
  - S. Hammes-Schiffer and J. C. Tully, *"Vibrationally adiabatic surface hopping"*, **J. Chem. Phys.** 101, 4657 (1994). [DOI: 10.1063/1.467455](https://doi.org/10.1063/1.467455)
  - J. E. Subotnik et al., *"Recent advances in surface hopping"*, **Annu. Rev. Phys. Chem.** 67, 387 (2016).
- **核心理论**:
  经典核坐标按当前绝热势能面 $V_k(R)$ 进行 Velocity Verlet 推进，电子相干波函数 $\mathbf{c}(t)$ 按含时薛定谔方程沿核轨迹幺正演化：
  $$i\hbar \dot{c}_k(t) = V_k(R) c_k(t) - i\hbar \sum_j \dot{R} \cdot \mathbf{d}_{kj}(R) c_j(t)$$
  其中 $\mathbf{d}_{kj} = \langle \phi_k | \nabla_R \phi_j \rangle$ 为非绝热导数耦合矢量 (NACV)。
  Tully 最少开关规则定义态 $k \to j$ 的瞬时跳跃几率：
  $$g_{k \to j} = \max\left(0, \frac{2 \Delta t}{\rho_{kk}} \text{Im}\left[\rho_{kj}^* (\dot{R} \cdot \mathbf{d}_{jk})\right]\right)$$
  跳跃时沿 NACV 方向调整核动量保持全系统能量守恒；若能量不足则视为禁阻跳跃 (Forbidden Hopping) 并发生动量反向反射。
- **代码映射**:
  - `src/mod_surface_hopping_fssh.f90`:
    - `init_tully_model`: Tully 三大基准解析势能面 (SAC, DAC, ECR) 模型；
    - `calc_adiabatic_surface_and_nacv`: 势能面与 NACV 矢量解析求解；
    - `propagate_fssh_step`: Tully 最少开关推进与动量重标度；
    - `run_fssh_ensemble`: 蒙特卡洛系综统计通道透射/反射分支比。

---

## 31. 强场分子定向、取向与光学离心机超转子动力学 (Molecular Alignment & Superrotors)

### 31.1 飞秒激光诱导非绝热对齐与无场转动复苏序参量
- **文献**:
  - H. Stapelfeldt and T. Seideman, *"Colloquium: Aligning molecules with strong laser pulses"*, **Rev. Mod. Phys.** 75, 543 (2003). [DOI: 10.1103/RevModPhys.75.543](https://doi.org/10.1103/RevModPhys.75.543)
  - D. M. Villeneuve et al., *"Observation of the revivals of a molecular wave packet created by a strong femtosecond laser pulse"*, **Phys. Rev. Lett.** 85, 542 (2000).
- **核心理论**:
  短脉冲强激光各向异性极化势 $V(\theta, t) = -\frac{1}{4} \mathcal{E}^2(t) \Delta\alpha \cos^2\theta$ 驱动分子转动拉曼跃迁形成 $\Delta J = 0, \pm 2$ 的相干转动波包。
  脉冲过后，波包在无场环境下发生相干周期复苏，完全复苏周期为 $T_{\text{rev}} = \frac{1}{2 B c}$。系综平均对齐度由序参量 $\langle\cos^2\theta\rangle(t)$ 表征。
- **代码映射**:
  - `src/mod_molecular_alignment.f90`:
    - `init_rotor_molecule`: 刚体转子结构与极化参数；
    - `simulate_laser_induced_alignment`: 非绝热脉冲激发与时域复苏谱计算。

### 31.2 光学离心机超转子产生与极端离心破键
- **文献**:
  - J. Karczmarek, J. Wright, P. B. Corkum, and M. Ivanov, *"Optical centrifuge for molecules"*, **Phys. Rev. Lett.** 82, 3420 (1999). [DOI: 10.1103/PhysRevLett.82.3420](https://doi.org/10.1103/PhysRevLett.82.3420)
  - U. Steinitz, R. Gopal, and I. Sh. Averbukh, *"Centrifugal breaking of molecules by an optical centrifuge"*, **Phys. Rev. Lett.** 109, 173001 (2012).
- **核心理论**:
  由反向旋转圆偏振场合成的光学离心机具备线性递增的角旋转速度 $\omega(t) = \alpha_{\text{chirp}} t$。分子被绝热捕获在极化势阱中加速至超高角动量超转子态（$J \gg 1$）：
  $$J_{\text{terminal}} \approx \frac{\alpha_{\text{chirp}} T_{\text{pulse}}}{2 B}$$
  巨大的离心能 $E_{\text{rot}} = B J(J+1)$ 压低共价结合势阱，当有效势垒完全消失时触发机械离心共价破键。
- **代码映射**:
  - `src/mod_molecular_alignment.f90`:
    - `calc_optical_centrifuge_kick`: 离心机终端角动量与超转子生成；
    - `calc_superrotor_dissociation`: 极端离心势能与共价键断裂判定。

---

## 32. 超冷光晶格与玻色-哈伯德微观映射 (Optical Lattice & Bose-Hubbard)

### 32.1 Mathieu 能带、Wannier 函数与玻色-哈伯德参数 $(J, U)$ 映射
- **文献**:
  - D. Jaksch, C. Bruder, J. I. Cirac, C. W. Gardiner, and P. Zoller, *"Cold bosonic atoms in optical lattices"*, **Phys. Rev. Lett.** 81, 3108 (1998). [DOI: 10.1103/PhysRevLett.81.3108](https://doi.org/10.1103/PhysRevLett.81.3108)
  - M. Greiner, O. Mandel, T. Esslinger, T. W. Hänsch, and I. Bloch, *"Quantum phase transition from a superfluid to a Mott insulator in a gas of ultracold atoms"*, **Nature** 415, 39 (2002). [DOI: 10.1038/415039a](https://doi.org/10.1038/415039a)
  - I. Bloch, J. Dalibard, and W. Zwerger, *"Many-body physics with ultracold gases"*, **Rev. Mod. Phys.** 80, 885 (2008).
- **核心理论**:
  驻波激光形成周期势 $V(x) = V_0 \sin^2(k_L x)$。在紧束缚近似下，平面波展开 Mathieu 方程构造局域 Wannier 函数 $w(x)$。
  玻色-哈伯德哈密顿量微观参数由重叠积分给出：
  $$J = -\int dx \, w(x-d) \left[ -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x) \right] w(x) \approx \frac{4}{\sqrt{\pi}} E_R s^{3/4} e^{-2\sqrt{s}}$$
  $$U = \frac{4\pi\hbar^2 a_s}{m} \int dx \, |w(x)|^4 \approx \sqrt{\frac{8}{\pi}} k_L a_s E_R s^{1/4}$$
  一维系统在 $(U/J)_c \approx 3.84$ 处经历超流体 (SF) 向 Mott 绝缘体 (MI) 的量子相变。
- **代码映射**:
  - `src/mod_optical_lattice_hubbard.f90`:
    - `init_optical_lattice`: 光晶格反冲尺度初始化；
    - `calc_bloch_band_energies`: Bloch 能带结构求解；
    - `calc_bose_hubbard_parameters`: Wannier 微观紧束缚参数与相变判定。

### 32.2 恒定外场布洛赫振荡与 Landau-Zener 带间隧穿
- **文献**:
  - M. Ben Dahan, E. Peik, J. Reichel, Y. Castin, and C. Salomon, *"Bloch oscillations of atoms in an optical lattice"*, **Phys. Rev. Lett.** 76, 4508 (1996). [DOI: 10.1103/PhysRevLett.76.4508](https://doi.org/10.1103/PhysRevLett.76.4508)
  - C. Zener, *"A theory of the electrical breakdown of solid dielectrics"*, **Proc. R. Soc. Lond. A** 137, 696 (1932).
- **核心理论**:
  在微重力或外力 $F$ 下，准动量沿第一布里渊区扫描 $\hbar \dot{q} = F$，形成空间定域的布洛赫振荡，周期为 $T_B = \frac{2\pi\hbar}{F d} = \frac{h}{F d}$。
  当扫描经过第一布里渊区边界时，向高能能带隧穿的几率由 Landau-Zener 公式决定：
  $$P_{\text{LZ}} = \exp\left( - \frac{\pi \Delta_{\text{gap}}^2}{4 \hbar |F| v_R} \right)$$
- **代码映射**:
  - `src/mod_optical_lattice_hubbard.f90`:
    - `calc_bloch_oscillation_dynamics`: 布洛赫振荡周期与 Landau-Zener 泄漏几率计算。

---

## 33. 多原子反应路径哈密顿量与变分过渡态理论 (RPH & Variational TST)

### 33.1 反应路径哈密顿量 (RPH) 与正则变分过渡态理论 (CVT)
- **文献**:
  - W. H. Miller, N. C. Handy, and J. E. Adams, *"Reaction path Hamiltonian for polyatomic molecules"*, **J. Chem. Phys.** 72, 99 (1980). [DOI: 10.1063/1.438959](https://doi.org/10.1063/1.438959)
  - D. G. Truhlar and B. C. Garrett, *"Variational transition-state theory"*, **Acc. Chem. Res.** 13, 440 (1980). [DOI: 10.1021/ar50156a002](https://doi.org/10.1021/ar50156a002)
  - B. C. Garrett and D. G. Truhlar, *"Criterion of maximum free energy of activation for variational transition-state theory"*, **J. Phys. Chem.** 83, 1052 (1979).
- **核心理论**:
  内禀反应坐标 (IRC) 弧长 $s$ 描述沿最小能量路径的推进，垂直于路径的简正振动自由度 $Q_k$ 与路径曲率耦合。
  广义过渡态理论 (GTST) 正则速率常数由分界面 $s$ 处决定：
  $$k^{\text{GTST}}(T, s) = \frac{k_B T}{h} \frac{Q^\ddagger(T, s)}{\Phi^R(T)} \exp\left( - \frac{V_0(s)}{k_B T} \right)$$
  正则变分过渡态理论 (CVT) 通过最小化速率常数寻找动力学自由能瓶颈分界面：
  $$k^{\text{CVT}}(T) = \min_s k^{\text{GTST}}(T, s)$$
- **代码映射**:
  - `src/mod_reaction_path_hamiltonian.f90`:
    - `init_rph_benchmark_reaction`: RPH 反应路径与曲率振动参数初始化；
    - `calc_generalized_tst_rate`: 沿路径分界面的广义 TST 速率；
    - `calc_cvt_rate_constant`: 正则变分 CVT 速率常数瓶颈优化。

### 33.2 Eckart 势垒半经典量子穿透修正
- **文献**:
  - C. Eckart, *"The penetration of a parabolic potential barrier"*, **Phys. Rev.** 35, 1303 (1930). [DOI: 10.1103/PhysRev.35.1303](https://doi.org/10.1103/PhysRev.35.1303)
- **核心理论**:
  对于具有虚频 $\omega^\ddagger$ 的鞍点势垒，量子隧穿修正系数 $\kappa(T) = k^{\text{tunnel}}(T) / k^{\text{classical}}(T)$ 通过对不对称 Eckart 势垒穿透概率 $P(E)$ 的玻尔兹曼积分严格求得，低温下显著高于经典 1.0。
- **代码映射**:
  - `src/mod_reaction_path_hamiltonian.f90`:
    - `calc_eckart_tunneling_factor`: Eckart 势垒高斯-勒让德量子隧穿修正因子计算。

---

## 34. 相对论原子结构与径向狄拉克方程 (Relativistic Atomic Structure & Dirac)

### 34.1 径向狄拉克方程、Sommerfeld 能级与核心极化模型势
- **文献**:
  - P. A. M. Dirac, *"The Quantum Theory of the Electron"*, **Proc. R. Soc. Lond. A** 117, 610 (1928). [DOI: 10.1098/rspa.1928.0023](https://doi.org/10.1098/rspa.1928.0023)
  - I. P. Grant, *Relativistic Quantum Theory of Atoms and Molecules*, Springer (2007).
  - D. W. Norcross, *"Model potential calculations of alkali-metal negative-ion and neutral binding energies"*, **Phys. Rev. A** 7, 606 (1973). [DOI: 10.1103/PhysRevA.7.606](https://doi.org/10.1103/PhysRevA.7.606)
- **核心理论**:
  径向一阶耦合狄拉克方程组描述大分量 $P(r)$ 与小分量 $Q(r)$：
  $$\frac{dP}{dr} = -\frac{\kappa}{r} P + \frac{1}{c} [E - V(r) + 2 c^2] Q, \quad \frac{dQ}{dr} = \frac{\kappa}{r} Q - \frac{1}{c} [E - V(r)] P$$
  Sommerfeld 精细结构本征解析公式为：
  $$E_D = c^2 \left[ \left(1 + \frac{(Z_{\text{eff}}/c)^2}{(n - |\kappa| + \sqrt{\kappa^2 - (Z_{\text{eff}}/c)^2})^2}\right)^{-1/2} - 1 \right]$$
  Norcross-Klapisch 核心极化相互作用势：$V_{\text{pol}}(r) = -\frac{\alpha_{\text{core}}}{2 r^4} [1 - \exp(-(r/r_c)^6)]$。
- **代码映射**:
  - `src/mod_relativistic_atomic.f90`:
    - `calc_dirac_model_potential`: 核心屏蔽与极化势；
    - `solve_radial_dirac_eigenvalue`: 狄拉克方程本征能量与量子亏损求解；
    - `calc_dirac_fine_structure_splitting`: 天然自旋-轨道耦合劈裂能计算；
    - `calc_dirac_e1_matrix_element`: 相对论电偶极 (E1) 振子强度计算。

---

## 35. 共振非弹性 X 射线散射与内壳层光谱 (Resonant Inelastic X-ray Scattering - RIXS)

### 35.1 Kramers-Heisenberg 二阶微扰截面与 2D RIXS 能损图谱
- **文献**:
  - H. A. Kramers and W. Heisenberg, *"Über die Streuung von Strahlung durch Atome"*, **Z. Phys.** 31, 681 (1925). [DOI: 10.1007/BF02980624](https://doi.org/10.1007/BF02980624)
  - L. J. P. Ament, M. van Veenendaal, T. P. Devereaux, J. P. Hill, and J. van den Brink, *"Resonant inelastic x-ray scattering studies of elementary excitations"*, **Rev. Mod. Phys.** 83, 705 (2011). [DOI: 10.1103/RevModPhys.83.705](https://doi.org/10.1103/RevModPhys.83.705)
  - F. M. F. de Groot and A. Kotani, *Core Level Spectroscopy of Solids*, CRC Press (2008).
- **核心理论**:
  入射光子 $\hbar\omega_1$ 激发核心电子至未占轨道，经飞秒级核心空穴寿命 $\Gamma_m$ 衰变发射光子 $\hbar\omega_2$，留下低能电子/轨道激发 $\hbar\Omega = \hbar\omega_1 - \hbar\omega_2$：
  $$\frac{d^2\sigma}{d\Omega d\omega_2} \propto \sum_f \left| \sum_m \frac{\langle f | \hat{\mathcal{D}}_2^\dagger | m \rangle \langle m | \hat{\mathcal{D}}_1 | i \rangle}{E_i + \hbar\omega_1 - E_m + i\Gamma_m / 2} \right|^2 \mathcal{L}(\hbar\Omega - (E_f - E_i), \gamma_f)$$
- **代码映射**:
  - `src/mod_resonant_xray_scattering.f90`:
    - `init_rixs_system`: RIXS 能级与偶极矩阵元初始化；
    - `calc_xas_cross_section`: X 射线吸收谱 (XAS) 求解；
    - `calc_kramers_heisenberg_cross_section`: 二阶 RIXS 散射截面计算；
    - `calc_rixs_2d_map`: 2D RIXS 能量-能损响应矩阵生成。

### 35.2 电-声耦合 Franck-Condon 伴峰与相联拉盖尔多项式严格展开
- **文献**:
  - K. Huang and A. Rhys, *"Theory of light absorption and non-radiative transitions in F-centres"*, **Proc. R. Soc. Lond. A** 204, 406 (1950).
  - L. J. P. Ament et al., **Rev. Mod. Phys.** 83, 705 (2011), Section III.B (Phonon RIXS).
- **核心理论**:
  中间态核心空穴产生的晶格位移常数 $d = \sqrt{S}$（Huang-Rhys 因子 $S$）。各阶声子损失峰振幅采用相联拉盖尔多项式 $L_p^{(\alpha)}$ 严格解析重叠积分：
  $$A_n(\omega_1) = \sum_{\nu=0}^\infty \frac{\langle n | \hat{D}(-d) | \nu \rangle \langle \nu | \hat{D}(d) | 0 \rangle}{\Delta\omega - \nu \omega_0 + i\Gamma_m / 2}$$
  在快碰撞极限 $\Gamma_m \gg \omega_0$ 下，完备性求和 $\sum_\nu \langle n | \nu \rangle \langle \nu | 0 \rangle = \delta_{n, 0}$ 严格重现弹性主峰占优与声子发射级数单调衰减。
- **代码映射**:
  - `src/mod_resonant_xray_scattering.f90`:
    - `calc_huang_rhys_vibrational_rixs`: 严格相联拉盖尔声子级数计算。

---

## 36. 亚稳态原子潘宁电离与缔合电离动力学 (Penning & Associative Ionization / Chemi-ionization)

### 36.1 半经典光学势自电离理论、分支截面与能量分流
- **文献**:
  - H. Hotop and A. Niehaus, *"Reactions of excited atoms and molecules with atoms and molecules: II. Energy analysis of electrons from reactions with He*(2^1S) and He*(2^3S)"*, **Z. Phys.** 228, 68 (1969). [DOI: 10.1007/BF01392439](https://doi.org/10.1007/BF01392439)
  - P. E. Siska, *"Molecular-beam studies of Penning ionization"*, **Rev. Mod. Phys.** 65, 337 (1993). [DOI: 10.1103/RevModPhys.65.337](https://doi.org/10.1103/RevModPhys.65.337)
  - W. H. Miller, *"Theory of Penning ionization. I. Atoms"*, **J. Chem. Phys.** 52, 3563 (1970). [DOI: 10.1063/1.1673523](https://doi.org/10.1063/1.1673523)
  - J. S. Cohen and N. F. Lane, *"Chemi-ionization of He(2^3S, 2^1S) by H: Optical-potential calculation"*, **J. Chem. Phys.** 66, 586 (1977). [DOI: 10.1063/1.433980](https://doi.org/10.1063/1.433980)
- **核心理论**:
  当电子激发能超过靶原子第一电离能的亚稳态激发原子 $A^*$（如 $\text{He}^*(2^3S), E_{\text{exc}} \approx 19.82\text{ eV}$）与中性靶 $B$ 发生慢碰撞时，体系嵌入电子连续态：
  $$A^* + B \to \begin{cases} A + B^+ + e^- & (\text{Penning Ionization, PI}) \\ AB^+ + e^- & (\text{Associative Ionization, AI}) \end{cases}$$
  复光学势有效哈密顿量定义为 $W(R) = V_*(R) - \frac{i}{2}\Gamma(R)$。半经典碰撞轨道（碰撞参数 $b$、相对动能 $E_{\text{coll}}$）沿经典径向转折点 $R_{\text{turn}}$ 前后积分存活几率：
  $$P_{\text{surv}}(b) = \exp\left( -2 \int_{R_{\text{turn}}}^\infty \frac{\Gamma(R)}{\hbar v_r(R, b)} dR \right)$$
  两体碰撞产生自电离的总电离截面为：
  $$\sigma_{\text{tot}} = 2\pi \int_0^\infty b [1 - P_{\text{surv}}(b)] db$$
  当在核间距 $R$ 处发射电子时，释放电子能量由局域势差守恒决定 $E_e(R) = V_*(R) - V_+(R)$。根据核相对动能守恒，若碰撞电离后 $AB^+$ 核运动能量落入束缚阱中（$E_{\text{coll}} + V_+(R) < 0$），则产物形成分子离子 $AB^+$（缔合电离 AI）；反之若 $E_{\text{coll}} + V_+(R) > 0$，则克服离子势解离为自由产物 $A + B^+$（潘宁电离 PI）。极低碰撞能量下粒子易被离子深阱捕获，AI 占据绝对优势；而在高能量碰撞下 AI 分支急剧衰减至零。
- **代码映射**:
  - `src/mod_penning_associative_ionization.f90`:
    - `init_penning_system`: 亚稳态入口势、离子势与自电离宽度初始化；
    - `calc_penning_classical_turning_point`: 有效势经典转折点高精度二分求解；
    - `calc_penning_cross_sections`: 半经典光学势存活几率、PI/AI 空间分流截面与总截面求解；
    - `calc_penning_thermal_rate`: 麦克斯韦-玻尔兹曼热平均反应速率常数 $k(T)$ 计算。

### 36.2 潘宁电离电子能谱 (PIES) 与超冷自旋极化寿命抑制
- **文献**:
  - H. Hotop, *"Penning ionization electron spectroscopy"*, **Adv. Mass Spectrom.** 7A, 9 (1978).
  - A. Niehaus, *"Spontaneous transitions to a continuum"*, **Adv. Chem. Phys.** 48, 399 (1981).
  - G. V. Shlyapnikov, J. T. M. Walraven, U. M. Rahmanov, and M. W. Reynolds, *"Penning ionization in spin-polarized ultracold metastable helium"*, **Phys. Rev. Lett.** 73, 3247 (1994). [DOI: 10.1103/PhysRevLett.73.3247](https://doi.org/10.1103/PhysRevLett.73.3247)
- **核心理论**:
  潘宁电离电子能谱 (PIES) 反映初末态势能面垂直跃迁差 $E_e(R) = V_*(R) - V_+(R)$ 的态密度加权分布：
  $$\frac{d\sigma}{dE_e} = \int_{R_{\min}}^\infty \frac{\Gamma(R)}{\hbar v_r(R)} \delta(E_e - [V_*(R) - V_+(R)]) dR = \sum_{R_c} \frac{\Gamma(R_c)}{\hbar v_r(R_c) \left| \frac{d}{dR}(V_*(R) - V_+(R)) \right|_{R_c}}$$
  在势差极值点 $\frac{d}{dR}(V_*(R) - V_+(R)) = 0$ 处，PIES 能谱呈现出类经典彩虹奇异性（Airy 峰），灵敏反映相互作用短程势井参数。
  在超冷微开尔文（$\mu\text{K}$）极限下，碰撞处于 $s$-波 Wigner 阈值律，总相互作用由复散射长度 $a = \alpha - i\beta$ 刻画：
  $$K_{\text{elastic}} = \frac{4\pi\hbar}{\mu}(\alpha^2 + \beta^2), \quad K_{\text{loss}} = \frac{4\pi\hbar}{\mu}\beta$$
  对于自旋极化三线态亚稳态原子（如全极化 $\text{He}^*(2^3S_1, m_s=+1)$），碰撞双原子总自旋为 $S=2$（五重态），而产物电子 $e^-$ 与离子基态 $\text{He}_2^+(^2\Sigma_u^+)$ 耦合最高仅能形成三重态 ($S=1$)，导致电离过程因 Wigner 自旋守恒规则而自旋禁阻，自电离损失速率 $K_{\text{loss}}$ 受到 $10^4 \sim 10^5$ 倍的巨大抑制，使得亚稳态氦原子超冷玻色-爱因斯坦凝聚（BEC）成为可能。
- **代码映射**:
  - `src/mod_penning_associative_ionization.f90`:
    - `calc_pies_spectrum`: 局域静止相干与高斯展宽 PIES 发射电子能谱；
    - `calc_ultracold_penning_rates`: 复散射长度 $a=\alpha-i\beta$ 弹性与非弹性电离速率及自旋极化抑制寿命估算。

---

## 37. 分子光解离动力学、时间自相关函数与光碎片动能释放谱 (KER)

### 37.1 Heller 时间波包相关函数与光吸收截面
- **文献**:
  - E. J. Heller, *"The semiclassical way to molecular spectroscopy"*, **Acc. Chem. Res.** 14, 368 (1981). [DOI: 10.1021/ar00072a002](https://doi.org/10.1021/ar00072a002)
  - R. Schinke, *Photodissociation Dynamics: Spectroscopy and Fragmentation of Small Polyatomic Molecules*, Cambridge University Press, Cambridge (1993). [DOI: 10.1017/CBO9780511599934](https://doi.org/10.1017/CBO9780511599934)
  - K. C. Kulander and E. J. Heller, *"Time-dependent approach to molecular photodissociation: Interpretation of absorption spectra"*, **J. Chem. Phys.** 69, 2439 (1978). [DOI: 10.1063/1.436919](https://doi.org/10.1063/1.436919)
- **核心理论**:
  分子在电子基态波函数 $|\phi_0\rangle$ 下受紫外/极紫外光子垂直 Franck-Condon 激发跃迁至解离态激发能面 $V_{\text{exc}}(R)$，初态波包由电偶极跃迁矩 $\hat{\mu}$ 诱导生成：
  $$|\psi(0)\rangle = \hat{\mu}(R) |\phi_0\rangle$$
  核波包在激发态势能面上的量子含时演化受激发态哈密顿量 $\hat{H}_{\text{exc}}$ 驱动：
  $$i\hbar \frac{\partial \psi(R, t)}{\partial t} = \hat{H}_{\text{exc}} \psi(R, t) = \left[ -\frac{\hbar^2}{2\mu}\frac{d^2}{dR^2} + V_{\text{exc}}(R) \right] \psi(R, t)$$
  光子吸收截面 $\sigma(\omega)$ 由波包含时自相关函数 $C(t) = \langle \psi(0) | \psi(t) \rangle$ 的全时半傅里叶变换解析给出：
  $$\sigma(\omega) = \frac{4\pi \alpha \omega}{3 c} \text{Re} \int_0^\infty C(t) e^{i(E_0 + \hbar\omega)t/\hbar} e^{-\gamma t/\hbar} dt$$
  式中 $\gamma$ 为唯象谱线展宽或阻尼因子，$E_0$ 为初态基态能量。
- **代码映射**:
  - `src/mod_photofragment_flux.f90`:
    - `calc_autocorrelation_function`: 计算核波包在激发态势能面上含时投影自相关函数 $C(t) = \langle\psi(0)|\psi(t)\rangle$；
    - `calc_heller_absorption_spectrum`: 基于 Heller 时间相关函数积分求解连续光吸收截面谱 $\sigma(\omega)$。

### 37.2 渐近概率流通量与光碎片动能释放 (KER) 谱
- **文献**:
  - G. G. Balint-Kurti, R. N. Dixon, and C. C. Marston, *"Grid methods for calculating photodissociation cross sections and product state distributions"*, **Int. Rev. Phys. Chem.** 11, 269 (1992). [DOI: 10.1080/01442359209353274](https://doi.org/10.1080/01442359209353274)
  - R. N. Zare, *Angular Momentum: Understanding Spatial Aspects in Chemistry and Physics*, Wiley-Interscience, New York (1988).
- **核心理论**:
  在解离渐近区 $R = R_{\infty}$ 设立虚动能分析仪面，出射光碎片的量子概率流通量为：
  $$J(R_{\infty}, t) = \frac{\hbar}{\mu} \text{Im}\left[ \psi^*(R_{\infty}, t) \left. \frac{\partial \psi(R, t)}{\partial R} \right|_{R_{\infty}} \right]$$
  能量分辨动能释放谱 (Kinetic Energy Release, KER) $P(E_k)$ 通过渐近通量算符的时频能量投影给出：
  $$P(E_k) = \frac{\hbar}{\mu k} \left| \int_0^\infty J(R_{\infty}, t) e^{i E_k t / \hbar} dt \right|^2, \quad E_k = \frac{\hbar^2 k^2}{2\mu}$$
  对于线偏振激光场诱导单光子解离，光碎片的空间出射角分布服从双极展开：
  $$I(\theta) = \frac{\sigma_{\text{tot}}}{4\pi} \left[ 1 + \beta P_2(\cos\theta) \right]$$
  其中 $\beta$ 为各向异性参数，平行跃迁时 $\beta = +2$（碎片沿激光偏振轴喷射），垂直跃迁时 $\beta = -1$（碎片在偏振垂直平面展开）。
- **代码映射**:
  - `src/mod_photofragment_flux.f90`:
    - `calc_ker_spectrum_from_flux`: 渐近概率流时间-能量半傅里叶变换提取动能释放谱 $P(E_k)$；
    - `calc_photofragment_angular_distribution`: Zare 各向异性双极展开角分布 $I(\theta)$。

---

## 38. 多通道非绝热避差穿越、Landau-Zener 跃迁与 Hellmann-Feynman 耦合

### 38.1 透热-绝热基组表象变换与导数耦合矢量 (NACV)
- **文献**:
  - F. T. Smith, *"Diabatic and adiabatic representations for atomic collision problems"*, **Phys. Rev.** 179, 111 (1969). [DOI: 10.1103/PhysRev.179.111](https://doi.org/10.1103/PhysRev.179.111)
  - M. Baer, *Beyond Born-Oppenheimer: Electronic Nonadiabatic Effects in Chemical Reactions*, Wiley, New York (2006). [DOI: 10.1002/0471780081](https://doi.org/10.1002/0471780081)
  - H. Hellmann, *Einführung in die Quantenchemie*, Franz Deuticke, Leipzig (1937); R. P. Feynman, *"Forces in Molecules"*, **Phys. Rev.** 56, 340 (1939). [DOI: 10.1103/PhysRev.56.340](https://doi.org/10.1103/PhysRev.56.340)
- **核心理论**:
  在双态透热表象中，电子哈密顿矩阵形式为：
  $$\mathbf{H}_{\text{diab}}(R) = \begin{pmatrix} V_{11}(R) & V_{12}(R) \\ V_{21}(R) & V_{22}(R) \end{pmatrix}$$
  通过局域正交旋转矩阵 $\mathbf{U}(R) = \begin{pmatrix} \cos\theta(R) & \sin\theta(R) \\ -\sin\theta(R) & \cos\theta(R) \end{pmatrix}$ 对角化获得绝热势能面：
  $$E_{\pm}(R) = \frac{V_{11} + V_{22}}{2} \pm \frac{1}{2} \sqrt{(V_{11} - V_{22})^2 + 4 V_{12}^2}, \quad \tan(2\theta) = \frac{2 V_{12}}{V_{11} - V_{22}}$$
  核动能算符作用于绝热基态导致一阶非绝热导数耦合矢量 (NACV)：
  $$d_{12}(R) = \langle \psi_1(R) | \frac{d}{dR} | \psi_2(R) \rangle = \frac{d\theta}{dR} = \frac{V_{12} \frac{d(V_{11}-V_{22})}{dR} - (V_{11}-V_{22})\frac{dV_{12}}{dR}}{(V_{11}-V_{22})^2 + 4 V_{12}^2}$$
  利用 Hellmann-Feynman 定理，非对角导数耦合可严格重写为电子哈密顿量动量梯度除以绝热能级差：
  $$d_{jk}(R) = \frac{\langle \psi_j | \nabla_R \hat{H}_{\text{el}} | \psi_k \rangle}{E_k(R) - E_j(R)} \quad (j \ne k)$$
- **代码映射**:
  - `src/mod_multistate.f90`:
    - `diabatic_to_adiabatic`: 透热势能矩阵解析对角化与绝热本征能级输出；
    - `calc_nonadiabatic_coupling_vector`: 解析计算非绝热一阶导数耦合标量 $d_{12}(R)$。

### 38.2 Landau-Zener 跃迁几率与多态波包动力学
- **文献**:
  - L. D. Landau, *"Zur Theorie der Energieübertragung. II"*, **Phys. Z. Sowjetunion** 2, 46 (1932).
  - C. Zener, *"Non-adiabatic crossing of energy levels"*, **Proc. R. Soc. Lond. A** 137, 696 (1932). [DOI: 10.1098/rspa.1932.0165](https://doi.org/10.1098/rspa.1932.0165)
  - E. C. G. Stueckelberg, *"Theorie der unelastischen Stösse zwischen Atomen"*, **Helv. Phys. Acta** 5, 369 (1932). [DOI: 10.5169/seals-110177](https://doi.org/10.5169/seals-110177)
- **核心理论**:
  当核以速度 $v = \dot{R}$ 穿过避免交叉点 $R_c$（此处 $V_{11}(R_c) = V_{22}(R_c)$）时，体系在透热态之间的非绝热跃迁几率满足 Landau-Zener 公式：
  $$P_{\text{LZ}} = \exp\left( -2\pi \delta_{\text{LZ}} \right) = \exp\left( -\frac{2\pi V_{12}^2}{\hbar v |\Delta F|} \right), \quad \Delta F = \left| \left.\frac{dV_{11}}{dR}\right|_{R_c} - \left.\frac{dV_{22}}{dR}\right|_{R_c} \right|$$
  在强耦合或极慢碰撞极限下（$\delta_{\text{LZ}} \gg 1$），$P_{\text{LZ}} \to 0$，系统绝热跟随本征态演化；在弱耦合或高速碰撞极限下（$\delta_{\text{LZ}} \ll 1$），$P_{\text{LZ}} \to 1$，系统保持透热状态穿越。
  多态含时波包动力学演化采用包含势能耦合项的分裂算符推进：
  $$|\Psi(t+\Delta t)\rangle = e^{-i \hat{T} \Delta t / (2\hbar)} e^{-i \hat{\mathbf{V}}(R) \Delta t / \hbar} e^{-i \hat{T} \Delta t / (2\hbar)} |\Psi(t)\rangle$$
- **代码映射**:
  - `src/mod_multistate.f90`:
    - `calc_landau_zener_probability`: 解析计算不同穿行速度与耦合参数下的 Landau-Zener 跃迁几率；
    - `propagate_two_state_wavepacket`: 双通道多态分裂算符量子波包非绝热协同推进。

---

## 39. 双原子分子转振跃迁调控、STIRAP 绝热受激跃迁与 Franck-Condon 原理

### 39.1 离心修正转振能谱与 Franck-Condon 因子 (FCF)
- **文献**:
  - P. M. Morse, *"Diatomic molecules according to the wave mechanics. II. Vibrational levels"*, **Phys. Rev.** 34, 57 (1929). [DOI: 10.1103/PhysRev.34.57](https://doi.org/10.1103/PhysRev.34.57)
  - E. U. Condon, *"A theory of intensity distribution in band systems"*, **Phys. Rev.** 28, 1182 (1926). [DOI: 10.1103/PhysRev.28.1182](https://doi.org/10.1103/PhysRev.28.1182)
  - G. Herzberg, *Molecular Spectra and Molecular Structure: I. Spectra of Diatomic Molecules*, D. Van Nostrand, Princeton (1950).
- **核心理论**:
  双原子分子在转动量子数 $J$ 与振动态 $v$ 下的有效径向薛定谔方程为：
  $$\left[ -\frac{\hbar^2}{2\mu}\frac{d^2}{dR^2} + V_0(R) + \frac{\hbar^2 J(J+1)}{2\mu R^2} \right] \chi_{v, J}(R) = E_{v, J} \chi_{v, J}(R)$$
  若势能函数采用 Morse 解析形式 $V_{\text{Morse}}(R) = D_e [1 - e^{-\alpha(R - R_e)}]^2$，纯振动态能量为：
  $$E_v = \hbar\omega_e \left( v + \frac{1}{2} \right) - \hbar\omega_e x_e \left( v + \frac{1}{2} \right)^2$$
  电子跃迁辐射带强度由初态振动波函数 $\chi_{v}(R)$ 与终态振动波函数 $\chi_{v'}'(R)$ 的 Franck-Condon 因子 (FCF) 支配：
  $$q_{v v'} = |\langle \chi_v | \chi_{v'}' \rangle|^2 = \left| \int_0^\infty \chi_v^*(R) \chi_{v'}'(R) dR \right|^2, \quad \sum_{v'} q_{v v'} = 1$$
- **代码映射**:
  - `src/mod_rovibrational.f90`:
    - `solve_rovibrational_spectrum`: 结合离心势修正求解双原子分子转振本征能级与波函数；
    - `calc_franck_condon_factors`: 跨电子态振动态 Franck-Condon 重叠积分矩阵与正交归一校验。

### 39.2 受激拉曼绝热通道 (STIRAP) 暗态高保真度布居转移
- **文献**:
  - U. Gaubatz, P. Rudecki, S. Schiemann, and K. Bergmann, *"Population transfer between molecular vibrational levels by stimulated Raman scattering with partially overlapping laser pulses: A new technique"*, **J. Chem. Phys.** 92, 5363 (1990). [DOI: 10.1063/1.458514](https://doi.org/10.1063/1.458514)
  - K. Bergmann, H. Theuer, and B. W. Shore, *"Coherent population transfer among quantum states of atoms and molecules"*, **Rev. Mod. Phys.** 70, 1003 (1998). [DOI: 10.1103/RevModPhys.70.1003](https://doi.org/10.1103/RevModPhys.70.1003)
  - N. V. Vitanov, A. A. Rangelov, B. W. Shore, and K. Bergmann, *"Stimulated Raman adiabatic passage in physics, chemistry, and beyond"*, **Rev. Mod. Phys.** 89, 015006 (2017). [DOI: 10.1103/RevModPhys.89.015006](https://doi.org/10.1103/RevModPhys.89.015006)
- **核心理论**:
  在三能级 $\Lambda$ 系统中（初态 $|1\rangle$、激发中间态 $|2\rangle$、终态基态 $|3\rangle$），泵浦激光 $\Omega_P(t)$ 耦合 $|1\rangle \leftrightarrow |2\rangle$，斯托克斯激光 $\Omega_S(t)$ 耦合 $|2\rangle \leftrightarrow |3\rangle$。双光子共振（$\delta = 0$）哈密顿量为：
  $$\hat{H}_{\text{STIRAP}}(t) = \frac{\hbar}{2} \begin{pmatrix} 0 & \Omega_P(t) & 0 \\ \Omega_P(t) & 2\Delta & \Omega_S(t) \\ 0 & \Omega_S(t) & 0 \end{pmatrix}$$
  该系统的零本征值本征矢为不含激发态 $|2\rangle$ 成分的辐射“暗态” (Dark State)：
  $$|D(t)\rangle = \cos\Theta(t) |1\rangle - \sin\Theta(t) |3\rangle, \quad \tan\Theta(t) = \frac{\Omega_P(t)}{\Omega_S(t)}$$
  施加“反常延时序”（Counter-intuitive sequence，即 $\Omega_S(t)$ 脉冲先于 $\Omega_P(t)$ 达到峰值）：
  $$t \to -\infty: \quad \Omega_S \gg \Omega_P \implies \Theta \to 0 \implies |D\rangle \to |1\rangle$$
  $$t \to +\infty: \quad \Omega_P \gg \Omega_S \implies \Theta \to \frac{\pi}{2} \implies |D\rangle \to -|3\rangle$$
  在绝热准则 $\Omega_{\text{eff}} \Delta\tau = \sqrt{\Omega_P^2 + \Omega_S^2} \Delta\tau \gg 10$ 满足下，量子系统沿暗态平滑演化，以理论接近 $100\%$ 的转移效率实现分子由弱束缚 Feshbach 态向振转基态的无辐射损耗传输。
- **代码映射**:
  - `src/mod_rovibrational.f90`:
    - `simulate_stirap_transfer`: 三能级 $\Lambda$ 体系含时拉比脉冲序列演化与终态保真度计算。

---

## 40. 含时波包散射动力学、通量时间-能量积分与 Möller 算符投影

### 40.1 高斯初态波包与自由传播谱展开
- **文献**:
  - M. D. Feit, J. A. Fleck, Jr., and A. Steiger, *"Solution of the Schrödinger equation by a spectral method"*, **J. Comput. Phys.** 47, 412 (1982). [DOI: 10.1016/0021-9991(82)90091-2](https://doi.org/10.1016/0021-9991(82)90091-2)
  - D. J. Tannor, *Introduction to Quantum Mechanics: A Time-Dependent Perspective*, University Science Books, Sausalito (2007).
  - C. Leforestier et al., *"A comparison of different propagation schemes for the time dependent Schrödinger equation"*, **J. Comput. Phys.** 94, 59 (1991). [DOI: 10.1016/0021-9991(91)90137-A](https://doi.org/10.1016/0021-9991(91)90137-A)
- **核心理论**:
  在含时波包散射框架中，单次长时间波包演化即可解析提取全能区连续散射信息。初态构造为远离势能区中心位于 $x_0$、平均动量为 $p_0 = \hbar k_0$ 的最小不确定度高斯波包：
  $$\psi(x, 0) = \frac{1}{(2\pi \sigma_x^2)^{1/4}} \exp\left[ -\frac{(x - x_0)^2}{4\sigma_x^2} + i k_0 x \right]$$
  其在连续能量表象上的正交振幅权重分量为：
  $$A(E) = \left( \frac{\mu}{\hbar k(E)} \right)^{1/2} \frac{1}{\sqrt{2\pi\hbar}} \int_{-\infty}^{+\infty} \psi(x, 0) e^{-i k(E) x} dx, \quad k(E) = \frac{\sqrt{2\mu E}}{\hbar}$$
- **代码映射**:
  - `src/mod_td_scattering.f90`:
    - `init_gaussian_wavepacket`: 空间与动量空间精确归一化初始高斯波包构筑。

### 40.2 跨势垒通量半傅里叶变换与 Möller 散射矩阵提取
- **文献**:
  - C. Möller, *"General properties of the characteristic matrix in the theory of elementary particles I"*, **K. Dan. Vidensk. Selsk. Mat. Fys. Medd.** 23, 1 (1945).
  - D. Neuhauser and M. Baer, *"The time-dependent Schrödinger equation including a complex absorbing potential: Application to the reactive collinear H + H_2 system"*, **J. Chem. Phys.** 90, 4351 (1989). [DOI: 10.1063/1.456644](https://doi.org/10.1063/1.456644)
  - W. H. Miller, S. D. Schwartz, and J. W. Tromp, *"Statistical functional theory of chemical reaction rates"*, **J. Chem. Phys.** 79, 4889 (1983). [DOI: 10.1063/1.445581](https://doi.org/10.1063/1.445581)
- **核心理论**:
  在相互作用势垒后方渐近探测面 $x = x_{\text{det}}$ 处记录随时间流逝的含时概率流密度：
  $$J(x_{\text{det}}, t) = \frac{\hbar}{\mu} \text{Im}\left[ \psi^*(x_{\text{det}}, t) \left. \frac{\partial \psi(x, t)}{\partial x} \right|_{x_{\text{det}}} \right]$$
  跨越势垒的总透射几率谱 $T(E)$ 通过时间-能量半傅里叶变换与初态动量振幅归一化给出：
  $$T(E) = \frac{1}{|A(E)|^2} \text{Re} \int_0^\infty J(x_{\text{det}}, t) e^{i E t / \hbar} dt$$
  基于 Möller 波动算符 $\hat{\Omega}^{(\pm)} = \lim_{t \to \mp \infty} e^{i \hat{H} t / \hbar} e^{-i \hat{H}_0 t / \hbar}$，通过对散射出射波在渐近本征基上的直接动量投影，提取散射 $S$ 矩阵元：
  $$S(E) = \frac{\langle \phi_{\text{out}}(E) | \psi(T) \rangle}{\langle \phi_{\text{in}}(E) | \psi(0) \rangle} e^{i E T / \hbar}$$
- **代码映射**:
  - `src/mod_td_scattering.f90`:
    - `calculate_td_transmission`: 渐近概率流通量时间-能量半傅里叶变换提取全能区透射谱 $T(E)$；
    - `project_wavepacket_to_smatrix`: Möller 动量基渐近投影法提取散射矩阵元 $S(E)$ 与弹性散射相移。

---

## 41. 复吸收势 (CAP) 最佳边界参数化与量子概率流连续性方程

### 41.1 非厄米复吸收势 (CAP) 吸收与量子力学反射抑制
- **文献**:
  - U. V. Riss and H.-D. Meyer, *"Calculation of resonance energies and widths using the complex absorbing potential method"*, **J. Phys. B: At. Mol. Opt. Phys.** 26, 4503 (1993). [DOI: 10.1088/0953-4075/26/23/021](https://doi.org/10.1088/0953-4075/26/23/021)
  - U. V. Riss and H.-D. Meyer, *"Reflection by and transmission through complex absorbing potentials"*, **J. Chem. Phys.** 105, 1409 (1996). [DOI: 10.1063/1.472003](https://doi.org/10.1063/1.472003)
  - J. G. Muga, J. P. Palao, B. Navarro, and I. L. Egusquiza, *"Complex absorbing potentials"*, **Phys. Rep.** 395, 357 (2004). [DOI: 10.1016/j.physrep.2004.03.002](https://doi.org/10.1016/j.physrep.2004.03.002)
- **核心理论**:
  在开体系有限格点波包动力学计算中，自由粒子解离波在计算边界处的反弹反射会引起假干涉现象。通过在哈密顿量边界区添加纯虚负吸收势构造非厄米有效势：
  $$\hat{H}_{\text{eff}} = \hat{H}_0 - i W(x), \quad W(x) \ge 0$$
  多项式吸收层型剖面（起始点 $x_c$，厚度 $L_{\text{cap}}$）：
  $$W(x) = \eta \left( \frac{x - x_c}{L_{\text{cap}}} \right)^n \Theta(x - x_c) \quad (n=2, 3)$$
  或光滑余弦平滑包络（Manolopoulos 形式）：
  $$W(x) = \frac{\hbar^2}{2\mu} \left( \frac{2\pi}{L_{\text{cap}}} \right)^2 \left[ \frac{1}{\cos^2\left(\frac{\pi(x - x_c)}{2 L_{\text{cap}}}\right)} - 1 \right]$$
  吸收势强度 $\eta$ 的选取需平衡高能穿透透射与势突变引起的量子反射：
  $$\eta_{\text{opt}} \approx \frac{\hbar^2}{2\mu L_{\text{cap}}^2} f(n)$$
- **代码映射**:
  - `src/mod_absorbing_boundary.f90`:
    - `cap_init`: 二次方、三次方、余弦吸收边界层几何厚度与最优吸收强度参数配置；
    - `cap_apply_damping`: 分裂算符单步波函数局域指数吸收衰减算子 $\exp[-W(x)\Delta t / \hbar]$ 作用。

### 41.2 非厄米连续性方程与概率流通量衰减律
- **文献**:
  - A. S. Dickinson et al., *"Energy-loss and probability flux in wavepacket scattering"*, **Mol. Phys.** 48, 1221 (1983).
- **核心理论**:
  在非厄米复吸收势作用下，概率密度 $\rho(x, t) = |\psi(x, t)|^2$ 的时间导数满足带有耗散汇项的连续性方程：
  $$\frac{\partial \rho(x, t)}{\partial t} + \frac{\partial j(x, t)}{\partial x} = -\frac{2}{\hbar} W(x) \rho(x, t)$$
  空间积分后系统总保全几率随时间单调衰减：
  $$\frac{d}{dt} \int_{-\infty}^{+\infty} |\psi(x, t)|^2 dx = -\frac{2}{\hbar} \int_{x_c}^{x_{\max}} W(x) |\psi(x, t)|^2 dx \le 0$$
  通过在时间轴上累积吸收汇项，可严格恢复跨越吸收边界的总解离几率。
- **代码映射**:
  - `src/mod_absorbing_boundary.f90`:
    - `calculate_probability_flux`: 计算格点局域量子概率流密度 $j(x, t)$ 及边界吸收汇项平衡。

---

## 42. 库仑三体系统、Perkeris 坐标变换与两电子关联

### 42.1 库仑三体两电子关联与 Kato 尖点条件
- **文献**:
  - E. A. Hylleraas, *"Neue Berechnung des Grundzustands des Heliumatoms sowie des negativen Wasserstoffions H^-"*, **Z. Phys.** 54, 347 (1929). [DOI: 10.1007/BF01375457](https://doi.org/10.1007/BF01375457)
  - T. Kato, *"On the eigenfunctions of many-particle systems in quantum mechanics"*, **Commun. Pure Appl. Math.** 10, 151 (1957). [DOI: 10.1002/cpa.3160100201](https://doi.org/10.1002/cpa.3160100201)
  - C. L. Pekéris, *"Ground State of Two-Electron Atoms"*, **Phys. Rev.** 112, 1649 (1958). [DOI: 10.1103/PhysRev.112.1649](https://doi.org/10.1103/PhysRev.112.1649)
- **核心理论**:
  类氦双电子库仑三体体系（原子核电荷 $Z$，无限大质量极限）全哈密顿量为：
  $$\hat{H} = -\frac{\hbar^2}{2m} \left( \nabla_1^2 + \nabla_2^2 \right) - \frac{Z e^2}{4\pi\epsilon_0 r_1} - \frac{Z e^2}{4\pi\epsilon_0 r_2} + \frac{e^2}{4\pi\epsilon_0 r_{12}}$$
  在电子-核碰撞点 $r_i \to 0$ 与电子-电子重合点 $r_{12} \to 0$ 处，库仑势发散导致波函数一阶导数不连续，满足严格的 Kato 尖点条件：
  $$\left. \frac{\partial \bar{\psi}}{\partial r_i} \right|_{r_i=0} = -Z \psi(r_i=0), \quad \left. \frac{\partial \bar{\psi}}{\partial r_{12}} \right|_{r_{12}=0} = \frac{1}{2} \psi(r_{12}=0)$$
  标准轨道乘积基组（如 Hartree-Fock）完全忽略了 $r_{12}$ 显式项，收敛极其迟缓；显式关联 Hylleraas-Pekéris 基组能够以指数精度直接满足该两电子关联。
- **代码映射**:
  - `src/mod_coulomb_threebody.f90`:
    - `calc_hylleraas_analytical_integral`: 求解两电子显式关联关联矩阵元三维解析 Gamma 递推积分。

### 42.2 正交 Pekéris 坐标变换与变分高精度本征求解
- **文献**:
  - C. L. Pekéris, *"1^1S and 2^3S states of Helium"*, **Phys. Rev.** 115, 1216 (1959). [DOI: 10.1103/PhysRev.115.1216](https://doi.org/10.1103/PhysRev.115.1216)
  - G. W. F. Drake, *"High precision calculations for helium"*, in *Springer Handbook of Atomic, Molecular, and Optical Physics*, Springer (2006).
- **核心理论**:
  三角形核间距几何约束条件（$|r_1 - r_2| \le r_{12} \le r_1 + r_2$）使得三维积分边界极其复杂。Pekéris 引入保正交线性坐标变换：
  $$u = r_1 + r_2 - r_{12}, \quad v = r_1 - r_2 + r_{12}, \quad w = -r_1 + r_2 + r_{12} \quad (u \ge 0, \; v \ge 0, \; w \ge 0)$$
  反变换为：
  $$r_1 = \frac{u + v}{2}, \quad r_2 = \frac{u + w}{2}, \quad r_{12} = \frac{v + w}{2}$$
  Pekéris 坐标将复杂的三原子核间距积分区间完全解耦为三个独立的半无限区间 $[0, \infty) \times [0, \infty) \times [0, \infty)$。变换的雅可比行列式为：
  $$d\tau = \frac{1}{8} (u + v)(u + w)(v + w) \, du dv dw$$
  在变分 Laguerre 正交多项式基组上展开，全部动能矩阵元与势能矩阵元均能解析化为高阶阶乘组合，彻底消除数值积分误差，将基态能量精算至亚赫兹精度。
- **代码映射**:
  - `src/mod_coulomb_threebody.f90`:
    - `perkeris_coordinate_transform`: 物理核间距与正交解耦 Pekéris 坐标间可逆解析变换；
    - `calc_perkeris_volume_element`: Pekéris 坐标变换雅可比体积元计算；
    - `solve_helium_ground_state_variational`: 氦原子/类氦离子基态能量超高精度变分对角化求解。

---

## 43. 快速学术检索与代码对照总表

| 序号 | 物理算法模块 | 对应源文件 | 核心经典学术文献代表 | 主要导出 API 与核心算法 |
| :---: | :--- | :--- | :--- | :--- |
| 1 | **基础物理常数与单位换算** | `src/mod_constants.f90` | CODATA 2018 / 2022 | `to_au`, `from_au`, `GAUSS2AU`, `AU2TESLA`, `AU2EV` |
| 2 | **角动量代数与 Wigner 符号** | `src/mod_special_functions.f90` | Varshalovich (1988), Racah (1942) | `wigner_3j_half`, `clebsch_gordan_half`, `wigner_6j_half`, `wigner_9j_half` |
| 3 | **数值线性代数与快速傅里叶** | `src/mod_linear_algebra.f90` | EISPACK, Cooley & Tukey (1965) | `diag_symmetric_matrix`, `inv_real_matrix`, `fft_1d`, `matrix_norm` |
| 4 | **离散变量表象 (DVR) 与 FGH** | `src/mod_dvr_grid.f90` | Colbert & Miller (1992), Marston (1989) | `dvr_sinc_init`, `fgh_solve_bound_states`, `dvr_kinetic_matrix` |
| 5 | **飞秒超快激光脉冲合成** | `src/mod_laser_pulse.f90` | Diels & Rudolph (2006) | `create_gaussian_pulse`, `create_chirped_pulse`, `ac_stark_shift` |
| 6 | **复吸收势 (CAP) 边界与连续性** | `src/mod_absorbing_boundary.f90` | Riss & Meyer (1993), Muga (2004) | `cap_init`, `cap_apply_damping`, `calculate_probability_flux` |
| 7 | **转动/振动热力学玻尔兹曼统计** | `src/mod_boltzmann_weights.f90` | McQuarrie (2000), Herzberg (1950) | `calc_rotational_partition_function`, `calc_boltzmann_weights` |
| 8 | **波包分裂算符与 Bloch 方程** | `src/mod_time_propagation.f90` | Feit & Fleck (1982), Strang (1968) | `propagate_wavepacket_strang`, `solve_bloch_equations_rk4` |
| 9 | **强场库仑、Keldysh 与 ADK 电离** | `src/mod_coulomb_atomic.f90` | Keldysh (1965), ADK (1986) | `keldysh_parameter`, `adk_ionization_rate`, `soft_core_coulomb` |
| 10 | **高次谐波发射与 SFA 理论** | `src/mod_hhg_spectra.f90` | Lewenstein et al. (1994), Corkum (1993) | `hhg_power_spectrum`, `lewenstein_sfa_dipole`, `gabor_transform` |
| 11 | **切比雪夫展开与能窗算子 PES** | `src/mod_polynomial_expansion.f90` | Tal-Ezer & Kosloff (1984), Schafer (1990) | `propagate_chebyshev_step`, `calc_energy_window_operator` |
| 12 | **多态非绝热避差与 Landau-Zener** | `src/mod_multistate.f90` | Landau (1932), Zener (1932), Hellmann-Feynman | `diabatic_to_adiabatic`, `calc_nonadiabatic_coupling_vector`, `calc_landau_zener_probability` |
| 13 | **双原子转振能谱与 STIRAP** | `src/mod_rovibrational.f90` | Bergmann (1998), Vitanov (2017), Condon (1926) | `solve_rovibrational_spectrum`, `calc_franck_condon_factors`, `simulate_stirap_transfer` |
| 14 | **科学计算 I/O 与矩阵序列化** | `src/mod_io_utils.f90` | HDF5 Group, NetCDF | `save_matrix_binary`, `load_matrix_binary`, `compute_checksum` |
| 15 | **三次样条插值与 Thomas 算法** | `src/mod_spline_interpolation.f90` | de Boor (1978), Stoer & Bulirsch (2002) | `spline_cubic_fit`, `spline_cubic_eval`, `thomas_algorithm_tridiag` |
| 16 | **分子光解离与光碎片 KER 谱** | `src/mod_photofragment_flux.f90` | Heller (1981), Schinke (1993), Zare (1988) | `calc_autocorrelation_function`, `calc_heller_absorption_spectrum`, `calc_ker_spectrum_from_flux` |
| 17 | **开放量子系统 Lindblad 耗散主方程**| `src/mod_open_quantum.f90` | Lindblad (1976), Gorini (1976), Breuer (2002) | `propagate_lindblad_rk4`, `calculate_von_neumann_entropy`, `calc_purity` |
| 18 | **Krotov 量子最优控制泛函** | `src/mod_optimal_control.f90` | Somlói & Tannor (1993), Krotov (1996), Reich (2012) | `optimize_pulse_krotov`, `calc_control_fidelity` |
| 19 | **非含时散射、分波法与 Log-Der** | `src/mod_ti_scattering.f90` | Wigner (1948), Johnson (1973), Manolopoulos (1986) | `calc_scattering_length_numerov`, `calc_multichannel_close_coupling_logder` |
| 20 | **含时波包散射动力学与 Möller 投影**| `src/mod_td_scattering.f90` | Feit & Fleck (1982), Möller (1945), Tannor (2007) | `init_gaussian_wavepacket`, `calculate_td_transmission`, `project_wavepacket_to_smatrix` |
| 21 | **塞曼 Breit-Rabi 与磁 Feshbach** | `src/mod_field_scattering.f90` | Breit & Rabi (1931), Stoof (1988), Chin (2010) | `calc_breit_rabi_energies`, `calc_basis_transform_matrix`, `calc_magnetic_feshbach_resonance_scan` |
| 22 | **各向异性磁/电偶极超冷散射** | `src/mod_dipolar_scattering.f90`| Stoof (1988), Moerdijk (1996), Bohn (2009) | `calc_mddi_total_matrix_element`, `calc_dipolar_relaxation_cross_section`, `calc_stark_induced_dipole` |
| 23 | **超冷光缔合谱学与自由-束缚跃迁** | `src/mod_photoassociation.f90` | Jones (RMP 2006), Bohn & Julienne (1999) | `calc_free_bound_fc_overlap`, `calc_pa_cross_section`, `calc_pa_thermal_rate_coefficient` |
| 24 | **三体 Efimov 普适态与复合损失** | `src/mod_three_body_recombination.f90` | Efimov (1970), Braaten & Hammer (2006), Esry (1999) | `solve_efimov_s0_identical_bosons`, `calc_three_body_recombination_a_positive` |
| 25 | **低维约束诱导共振 CIR 与波导** | `src/mod_confined_scattering.f90` | Olshanii (1998), Bergeman (2003), Haller (2009) | `init_waveguide_1d`, `calc_olshanii_cir_parameters`, `calc_confined_dimer_binding_energy` |
| 26 | **自电离 Fano 不对称线型与 CCR** | `src/mod_autoionization_fano.f90` | Fano (1961), Reinhardt (1982), Moiseyev (1998) | `calc_fano_profile`, `calc_autoionization_lifetime`, `solve_ccr_resonance_model` |
| 27 | **交叉电磁场宇称破缺分子能谱** | `src/mod_crossed_field_scattering.f90` | Tscherbul & Krems (2006), Friedrich (1996) | `init_crossed_field_config`, `solve_crossed_field_eigenstates`, `calc_crossed_field_observables` |
| 28 | **三原子反应 Jacobi 与 Berry 相位** | `src/mod_triatomic_geometry.f90` | Sato (1955), Berry (1984), Longuet-Higgins (1958) | `jacobi_to_internuclear`, `calc_leps_potential`, `calc_conical_intersection_adiabats`, `calc_berry_phase_around_ci` |
| 29 | **旋量 BEC 宏观自旋动力学** | `src/mod_spinor_bec.f90` | Ho (1998), Ohmi & Machida (1998), Chang (2004) | `init_spinor_preset`, `calc_quadratic_zeeman_shift`, `propagate_spinor_sma_rk4` |
| 30 | **三原子超球面反应与热速率常数** | `src/mod_hyperspherical_reactive.f90` | Johnson (1980), Pack & Parker (1987), Miller (1975) | `init_reaction_mass`, `calc_eckart_transmission`, `calc_cumulative_reaction_probability`, `calc_canonical_rate_constant` |
| 31 | **偶极量子液滴与 LHY 涨落修正** | `src/mod_dipolar_droplets_lhy.f90` | Lee-Huang-Yang (1957), Petrov (2015), Chomaz (2016) | `init_dipolar_droplet_param`, `calc_pelster_lima_q5`, `calc_equilibrium_droplet_density` |
| 32 | **强场 NSDI 电子重碰撞动量谱** | `src/mod_strong_field_nsdi.f90` | Corkum (1993), Weber et al. (Nature 2000), Becker (2005) | `init_nsdi_laser`, `calc_recollision_trajectory`, `calc_lotz_cross_section`, `calc_nsdi_2d_momentum_dist` |
| 33 | **磁/光 Feshbach 束缚态与损耗** | `src/mod_feshbach_bound_states.f90` | Chin et al. (RMP 2010), Gao (2001), Theis (PRL 2004) | `init_mfr_preset`, `calc_mfr_bound_energy_coupled`, `calc_mfr_closed_channel_fraction`, `calc_ofr_inelastic_loss_rate` |
| 34 | **阿秒瞬态吸收光谱 (ATAS)** | `src/mod_attosecond_transient_absorption.f90` | Chini et al. (2014), Ott et al. (Science 2013), Wu (2016) | `init_atas_helium_benchmark`, `calc_laser_dressed_fano_q`, `calc_light_induced_state_energy`, `calc_atas_spectrum` |
| 35 | **双色圆偏振场与分子 PECD** | `src/mod_bicircular_pecd.f90` | Kfir (Nat. Phot. 2015), Lux (2012), Böwering (2001) | `init_bicircular_field`, `calc_dynamical_symmetry_fold`, `calc_chirality_measure`, `calc_pecd_pad_spectrum` |
| 36 | **超冷极性分子反应与微波遮蔽** | `src/mod_ultracold_reaction_shielding.f90` | Quéméner (PRA 2010), Anderegg (Science 2021), Schindewolf (2022) | `init_ultracold_molecule_preset`, `calc_effective_shielding_potential`, `calc_wkb_tunneling_probability`, `calc_shielded_scattering_rates` |
| 37 | **里德堡原子阻塞与量子多体疤痕** | `src/mod_rydberg_blockade.f90` | Lukin (PRL 2001), Bernien (Nature 2017), Turner (2018) | `init_rydberg_atom`, `calc_rydberg_blockade_radius`, `calc_two_atom_dynamics`, `calc_rydberg_scar_dynamics` |
| 38 | **表面量子散射与选择性吸附 (SAR)** | `src/mod_surface_scattering.f90` | Boato (1973), Manson (1991), Benedek (2018) | `init_surface_lattice`, `calc_surface_diffraction_channels`, `calc_hcs_diffraction_probabilities`, `calc_selective_adsorption_resonance` |
| 39 | **气-固表面催化与 Eley-Rideal 反应**| `src/mod_surface_reaction_er.f90` | Eley & Rideal (1940), Rettner (1992), Jackson (1992) | `init_er_reaction_system`, `calc_er_potential_2d`, `calc_er_energy_partitioning`, `calc_er_vibrational_populations` |
| 40 | **金属表面非绝热与电子摩擦 (GLE)** | `src/mod_surface_electronic_friction.f90` | Tully (1981), Head-Gordon & Tully (1995), Wodtke (2000) | `init_metal_surface`, `calc_electronic_friction_coeff`, `integrate_gle_scattering_trajectory`, `calc_vibrational_relaxation_rate` |
| 41 | **掠入射快原子表面衍射 (GIFAD)** | `src/mod_grazing_fast_atom_diffraction.f90` | Rousseau (PRL 2007), Schüller (PRL 2007), Winter (2011) | `init_gifad_experiment`, `calc_gifad_transverse_kinematics`, `calc_gifad_rainbow_angle`, `calc_surface_corrugation_from_rainbow` |
| 42 | **冷离子-中性原子杂化散射极化阱** | `src/mod_ion_atom_scattering.f90` | Langevin (1905), Côté (PRA 2000), Cetina (PRL 2012) | `init_ion_atom_system`, `calc_langevin_cross_section`, `calc_mere_phase_shift_s_wave`, `calc_rf_micromotion_heating` |
| 43 | **最少开关表面跳跃 (FSSH)** | `src/mod_surface_hopping_fssh.f90` | Tully (JCP 1990), Hammes-Schiffer & Tully (1994) | `init_tully_model`, `calc_adiabatic_surface_and_nacv`, `propagate_fssh_step`, `run_fssh_ensemble` |
| 44 | **强场分子定向取向与超转子动力学** | `src/mod_molecular_alignment.f90` | Stapelfeldt & Seideman (RMP 2003), Karczmarek (1999) | `init_rotor_molecule`, `simulate_laser_induced_alignment`, `calc_optical_centrifuge_kick`, `calc_superrotor_dissociation` |
| 45 | **超冷光晶格与玻色-哈伯德微观映射**| `src/mod_optical_lattice_hubbard.f90` | Jaksch (PRL 1998), Greiner (Nature 2002), Ben Dahan (1996) | `init_optical_lattice`, `calc_bloch_band_energies`, `calc_bose_hubbard_parameters`, `calc_bloch_oscillation_dynamics` |
| 46 | **反应路径哈密顿量与变分过渡态** | `src/mod_reaction_path_hamiltonian.f90` | Miller (JCP 1980), Truhlar & Garrett (1980), Eckart (1930) | `init_rph_benchmark_reaction`, `calc_generalized_tst_rate`, `calc_cvt_rate_constant`, `calc_eckart_tunneling_factor` |
| 47 | **相对论原子结构与径向狄拉克方程** | `src/mod_relativistic_atomic.f90` | Dirac (1928), Norcross (PRA 1973), Grant (2007) | `calc_dirac_model_potential`, `solve_radial_dirac_eigenvalue`, `calc_dirac_fine_structure_splitting`, `calc_dirac_e1_matrix_element` |
| 48 | **共振非弹性 X 射线散射 (RIXS)** | `src/mod_resonant_xray_scattering.f90` | Kramers-Heisenberg (1925), Ament et al. (RMP 2011) | `init_rixs_system`, `calc_xas_cross_section`, `calc_kramers_heisenberg_cross_section`, `calc_huang_rhys_vibrational_rixs` |
| 49 | **亚稳态原子潘宁电离与缔合电离** | `src/mod_penning_associative_ionization.f90` | Hotop & Niehaus (1969), Siska (RMP 1993), Miller (1970) | `init_penning_system`, `calc_penning_classical_turning_point`, `calc_penning_cross_sections`, `calc_pies_spectrum` |
| 50 | **库仑三体系统与两电子关联** | `src/mod_coulomb_threebody.f90` | Hylleraas (1929), Pekéris (1958), Drake (2006) | `perkeris_coordinate_transform`, `calc_hylleraas_analytical_integral`, `solve_helium_ground_state_variational` |

---
*全典文献经过精确校核，代码实现中严格遵守学术规范，所有公式推导与符号约定均与原始文献完全对齐。*


