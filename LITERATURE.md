# GeneralModule 科学文献典藏与理论映射全典 (Literature Compendium)

本全典详尽整理了 `GeneralModule` 算法库所依据的权威经典与前沿学术文献。涵盖**超冷量子碰撞、多通道密耦与对数导数法、外场磁 Feshbach 共振、角动量表象变换、强场超快物理、波包动力学与最优控制**八大核心物理领域。每条文献均包含标准学术引用、DOI、理论物理机理及其在算法库源代码中的确切映射位置。

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
9. [快速学术检索与代码对照总表](#9-快速学术检索与代码对照总表)

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

## 9. 快速学术检索与代码对照总表

| 物理模块 | 对应源文件 | 核心经典文献代表 | 主要导出 API 与算法 |
| :--- | :--- | :--- | :--- |
| **基础常数与单位** | `mod_constants.f90` | CODATA 2018 / 2022 | `to_au`, `from_au` |
| **角动量代数** | `mod_special_functions.f90` | Varshalovich (1988) | `wigner_3j_half`, `clebsch_gordan_half`, `wigner_9j_half` |
| **数值线性代数** | `mod_linear_algebra.f90` | EISPACK / Cooley-Tukey | `diag_symmetric_matrix`, `inv_real_matrix`, `fft_1d` |
| **DVR 与格点谱方法**| `mod_dvr_grid.f90` | Colbert & Miller (1992) | `dvr_sinc_init`, `fgh_solve_bound_states` |
| **激光脉冲合成** | `mod_laser_pulse.f90` | Diels & Rudolph (2006) | `create_gaussian_pulse`, `create_chirped_pulse` |
| **复吸收势边界** | `mod_absorbing_boundary.f90` | Riss & Meyer (1993) | `cap_init`, `calculate_probability_flux` |
| **强场超快物理** | `mod_coulomb_atomic.f90` | Keldysh (1965), Corkum (1993) | `keldysh_parameter`, `adk_ionization_rate` |
| **高次谐波发射** | `mod_hhg_spectra.f90` | Lewenstein et al. (1994) | `hhg_power_spectrum`, `lewenstein_sfa_dipole` |
| **分子转振控制** | `mod_rovibrational.f90` | Bergmann, Theuer & Shore (1998)| `calc_franck_condon_factors`, `solve_rovibrational_spectrum` |
| **非含时单/多通道散射**| `mod_ti_scattering.f90` | Wigner (1948), Johnson (1973), Manolopoulos (1986) | `calc_scattering_length_numerov`, `calc_multichannel_close_coupling_logder`, `create_segmented_grid` |
| **含时波包散射** | `mod_td_scattering.f90` | Feit & Fleck (1982), Möller (1945) | `calculate_td_transmission`, `project_wavepacket_to_smatrix` |
| **外场多基组超冷散射**| `mod_field_scattering.f90` | Breit & Rabi (1931), Stoof (1988), Chin et al. (2010) | `calc_breit_rabi_energies`, `calc_basis_transform_matrix`, `calc_magnetic_feshbach_resonance_scan` |
| **开放量子系统** | `mod_open_quantum.f90` | Lindblad (1976), Gorini (1976) | `propagate_lindblad_rk4`, `calculate_von_neumann_entropy` |
| **量子最优控制** | `mod_optimal_control.f90` | Somlói & Tannor (1993), Krotov (1996) | `optimize_pulse_krotov` |

---
*全典文献经过精确校核，代码实现中严格遵守学术规范，所有公式推导与符号约定均与原始文献完全对齐。*
