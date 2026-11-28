# GPU 实测结果（NVIDIA A100-SXM4-40GB）

日期：2026-11-28 · 执行节点：192.168.6.24（192 核 Xeon，755 GB 内存）

| 项目 | 版本/参数 |
|---|---|
| GPU | NVIDIA A100-SXM4-40GB（sm_80/cc80） |
| 驱动 | 610.57.04 |
| GPU 编译器 | NVHPC 25.5 `nvfortran`（`-O2 -acc=gpu -gpu=cc80`），用户级安装 `~/nvhpc/inst129` |
| CPU 编译器 | gfortran 14.3.1（`-O2 -fopenmp`，192 线程） |

## 基准数字（`benchmarks/bench_gpu.f90`，5 次平均）

| 核 | CPU | GPU (A100) | 加速比 |
|---|---|---|---|
| `sop_apply_3d`（48³ 网格，M=8 项，双精度） | 749.8 ms | **8.01 ms** | **≈ 93.6×** |
| `fft_1d`（2²² 点） | 1154.6 ms | 1121.9 ms（CPU 路径，未卸载） | — |

## 结论与说明

1. **`sop_apply_3d` 是 GPU 高效目标**：三维张量与 1D 矩阵的收缩为带宽受限的稠密运算，
   A100 上 93.6× 加速。这对 MCTDH/SOP 类多维量子动力学的"作用算符"路径意义直接：
   时间演化步的核心开销可整体卸载。
2. **`fft_1d` 蝶形循环不做逐 m 卸载**（首次尝试在 2²² 规模下因每步隐式整组拷贝挂死，
   指令已移除并留注释）。FFT 的 GPU 路径应经由 `fft_2d` 的**列级并行**
   （每列独立 1D 变换，天然可并行），留作后续项。
3. **OpenACC 深拷贝陷阱（已修复并固化）**：派生类型组件数组（`sop%h1(m)%mat`）不会被
   隐式数据子句拷贝到设备，直接引用触发 `CUDA_ERROR_ILLEGAL_ADDRESS`。
   修复：每项先拷入连续局部数组（`hm1/hm2/hm3`），并以
   `!$acc data copyin(in_vec, hm1, hm2, hm3) copy(tmp1, tmp2, tmp3)` 显式包裹。
   gfortran 默认构建完全不受影响（指令即注释），全量 604 断言回归绿。

## 验收与构建

- CPU 基准与 GPU 基准为同一源码两编译器构建；`make GPU=1`（gfortran 回退路径）与
  `nvfortran -acc=gpu` 在 CPU 回退语义下行为一致。
- OpenMP 与 OpenACC 指令互斥启用（双指令互为注释），两种构建均通过全量单元测试。
