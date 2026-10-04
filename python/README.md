# pygenmod

`pygenmod` 是 [GeneralModule](../README.md) 的 Python 伴侣包：把 Fortran 库中常用的
前处理计算（原子参数、脉冲、DVR 初值、Feshbach 拟合）与出版级可视化封装成
NumPy 友好的接口，便于在 Jupyter、脚本或 CI 流水线中直接调用。

## 环境要求

- Python ≥ 3.8
- NumPy ≥ 1.20
- SciPy ≥ 1.5（QCT 模块的 EBK 根查找与作用量积分）
- Matplotlib ≥ 3.3

在无显示环境（CI 容器、SSH 会话、批处理节点）下，包在导入时会自动把
Matplotlib 切到 `Agg` 后端，不需要手动设置 `MPLBACKEND`。

## 安装

作为源码树的一部分直接使用：

```bash
cd python
python3 test_pygenmod.py      # 运行 12 组单元测试
```

或以标准 Python 包安装：

```bash
pip install ./python
```

## 模块一览

| 模块 | 内容 |
| --- | --- |
| `constants.py` | 原子单位制常量与 `to_au` / `from_au` 单位换算 |
| `pulse.py` | 激光脉冲包络 `pulse_envelope`、电场 `pulse_electric_field` 与 Stark 位移 |
| `dvr.py` | Sinc-DVR 网格构造与 FGH 束缚态求解 |
| `coulomb.py` | 原子势、软核库仑、Keldysh 参数、HHG 截止能与 ADK 电离率 |
| `hhg.py` | 偶极加速度与高次谐波功率谱 |
| `multistate.py` | Landau-Zener 跃迁概率与布居演化 |
| `rovibrational.py` | 转振态索引映射、Franck-Condon 因子、转动常数与跃迁偶极矩阵 |
| `scattering.py` | Numerov 散射长度、多通道紧密耦合、散射波函数与微分截面 |
| `field_scattering.py` | Breit-Rabi 能级、场致碰撞通道、自旋交换矩阵与 Feshbach 拟合 |
| `qct.py` | QCT 反应散射：LEPS 势能面、EBK 作用量初条件、Velocity-Verlet、不透明度函数、截面与热速率（教科书第 17 章） |
| `sop_hamiltonian.py` | 和积（SOP）哈密顿量的矩阵自由作用与 POTFIT 分解（第 18.1 节） |
| `tensor_train.py` | 张量列车：TT-SVD、舍入、内积、TT 算符应用（第 18.3 节） |
| `mctdh_core.py` | MCTDH 核心：A 系数与单粒子函数（SPF）联合传播（第 18.2 节） |
| `sparse_grid.py` | Smolyak 稀疏网格与嵌套 Clenshaw–Curtis 规则（第 18.4 节） |

可运行示例：`python3 examples/demo_qct_highdim.py` 在数秒内跑通第 17—18 章
的主流程（QCT 系综、SOP 作用、POTFIT、TT 舍入、Smolyak 积分、MCTDH 传播）。
| `visualizer.py` | 发表级绘图样式与波函数、脉冲、取向对齐动力学出图 |

## 快速示例

```python
import numpy as np
from pygenmod import calc_scattering_length_numerov, plot_scattering_length_wavefunction

r = np.linspace(0.01, 12.0, 600)
v_pot = np.where(r <= 2.0, -1.0, 0.0)      # 单位 Hartree 的吸引方势阱
a_s, u_wf = calc_scattering_length_numerov(r, v_pot, mass=1.0)

plot_scattering_length_wavefunction(r, v_pot, u_wf, a_s, save_path="scattering_length.png")
```

各 API 的完整参数说明与物理背景见仓库根目录的
[README](../README.md#-python-辅助分析套件-pygenmod) 与
[CONFIG_GUIDE.md](../CONFIG_GUIDE.md)。

## 测试

```bash
python3 test_pygenmod.py
```

该脚本使用标准库 `unittest`，不依赖 pytest；全部 28 项测试覆盖常数换算、
脉冲、DVR、HHG、多通道散射、冷原子场致散射、绘图冒烟测试，以及第 17—18 章
参考实现的验证层次：QCT 守恒律与统计一致性、SOP 作用对照稠密 Kronecker 和、
TT 舍入往返、Smolyak 两种构造逐点一致、MCTDH 全空间精确性/规范与守恒/SPF 收敛。
`docs/audit_textbook.py` 会对同一套不变量做独立复核。

## 性能与服务器多核利用

- **QCT 向量化与多进程扩展**：全部轨迹合并为 `(N, 3, 3)` 数组推进，
  Velocity-Verlet 每步仅做一次势能梯度求值（力复用），300 条轨迹单核约 0.26 s。
  在多核服务器上，可通过设置 `QCTConfig(..., n_workers=-1)` 自动打满全部 CPU 核心，
  或显式指定 `n_workers=16/32`，各个工作进程使用独立的种子流并行演化。
- **EBK 轨道表缓存**：初始条件的径向轨道按 `(v, j)` 缓存，整个系综只积分一次一维轨道；
- **MCTDH 矩阵共享**：SPF 基矩阵 eta 在 A 方程与平均场之间共享，每步只构造一次；
- **稀疏网格缓存**：Clenshaw–Curtis 规则带 LRU 缓存。
- 更大规模的 QCT/MCTDH 计算建议直接在服务器或批处理集群上执行，
  可先用小系综确认物理收敛再将轨迹数放大至 $10^4\sim 10^6$。运行 `python3 examples/demo_qct_highdim.py`
  可直接观察单核 vs 多核服务器加速效果。
