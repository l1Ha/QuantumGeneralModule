# pygenmod

`pygenmod` 是 [GeneralModule](../README.md) 的 Python 伴侣包：把 Fortran 库中常用的
前处理计算（原子参数、脉冲、DVR 初值、Feshbach 拟合）与出版级可视化封装成
NumPy 友好的接口，便于在 Jupyter、脚本或 CI 流水线中直接调用。

## 环境要求

- Python ≥ 3.8
- NumPy ≥ 1.20
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

该脚本使用标准库 `unittest`，不依赖 pytest；全部 12 组测试覆盖常数换算、
脉冲、DVR、HHG、多通道散射、冷原子场致散射与绘图冒烟测试。
