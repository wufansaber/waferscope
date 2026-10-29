# S4 独立运行与光栅验证记录

日期：2026-10-22。

## 结论

**S4 主线已经开始实际运行，不再只是选型。零级界面和一维裸光栅通过检查，三维颗粒超胞尚未收敛。**

本轮不将试验后端接入生产 GUI，不输出经过精度认证的三维颗粒相机图或 SNR。grcwa 保留为辅助对照，FEM 保留为独立离散方法；没有删除或替换旧模型。

## 运行环境

- S4 官方源码固定提交：`5b89848cfd15719ccaa04c853fe883ecf8d6ac15`。这是固定历史提交，不是宣称存在的官方发布标签。
- Lua 5.2.4 前端，项目私有 w64devkit GCC 16.2.0 工具链。
- 使用 S4 官方内置 RNP 特征值求解器，BLAS 使用现有私有 OpenBLAS DLL；不使用本机测试失败的旧 LAPACK 绑定路径。
- 程序：`build/s4_runtime/bin/S4.exe`。
- 构建脚本：`tools/build_s4_runtime.py`；输入校验和、编译命令和二进制校验和：`build/s4_runtime/build_manifest.json`。
- 源码编译兼容修复：`s4_prototype/patches/pre_api_gcc_compat.patch`，将厚度指针与零比较改为检查指针及其数值。不修改 RCWA 方程。
- 核心库和 Lua 前端按版本及特征值后端放在独立对象目录，避免旧 Makefile 漏依赖导致跨版本链接。
- 没有安装系统 WSL、改全局 PATH，也没有修改 NIST2/FEM 的 Python 依赖。

首次下载文件、解压工具链后，可重新执行：

```powershell
& build/rcwa_env/Scripts/python.exe -m tools.build_s4_runtime
```

此命令不是联网安装器；需要已有的三个校验和匹配下载文件和私有 OpenBLAS。构建脚本中的 `INPUTS` 是实际检查规则。

源文件来自 [S4 固定提交](https://github.com/victorliu/S4/tree/5b89848cfd15719ccaa04c853fe883ecf8d6ac15)、[Lua 5.2.4](https://www.lua.org/ftp/lua-5.2.4.tar.gz) 和 [w64devkit v2.10.1](https://github.com/skeeto/w64devkit/releases/tag/v2.10.1)。`build/` 不属于本轮 Git 源码交付；换机需要另行复制私有运行库或按校验规则重建。

## 验证顺序与结果

### 1. 零级与平面界面

覆盖真空、无损玻璃、吸收硅，s/p 偏振，0°/45°/75° 入射，以及一维旋转坐标与二维全局坐标。

检查复反射系数、功率、总场近场相位。无损情况核对 R+T=1；吸收硅额外核对有限厚度内的解析衰减，不能把底部测量面的功率直接当成界面透射率。

增加零对比度颗粒不改变背景、弱对比度信号连续趋零的检查。零对比度图形直接按均匀材料表示，避免为没有物理界面的材料生成退化图形特征值问题；球的分层仍保留于该退化测试中。

S4 检查 41 项通过，连同既有 grcwa 回归共 69 项通过。原始测试记录为 `output/s4_validation_20261022/tests.xml`。

### 2. 与 NIST 同截断的裸光栅

输入：266 nm，75°，p 偏振；光栅周期 500 nm、线宽 250 nm、深度 100 nm。使用 `n=5.57, k=0.387` 的旧硅输入。这个数值仅用于匹配，不是 266 nm 真实硅光学常数的标定。

S4 原生一维默认展开采用 731 模式，即 -365 到 +365；NIST RCW 同为 ±365 阶。

| 量 | 结果 |
|---|---:|
| NIST 总反射率 | 0.27780809342204854 |
| S4 总反射率 | 0.2778080934480334 |
| 总反射率绝对差 | 2.60e-11 |
| 总反射率相对差 | 9.35e-11 |
| 功率加权复电场相对差 | 5.29e-11 |
| 四个传播级最大功率相对差 | 2.54e-10 |

逐级比较包括弱衍射级 m=-2，而不仅是总功率。没有幅度拟合、相位拟合或强制左右对称。

**这证明了该固定输入与相同截断下两个实现的一致性，不等于无限阶理论误差为 1e-10，也不认证有颗粒、有限光栅、真实材料或成像。**

### 3. 三维超胞试算

同样入射与材料；颗粒直径 100 nm，颗粒 `n=1.59, k=0`，离晶圆间隙 10 nm，球体 8 个圆形截面层。超胞 0.5×1.0 μm，即一个横向周期内重复一个球，背景含两个 500 nm 周期。它是周期颗粒阵列，不是孤立球。

sample/reference 必须同输入、同基组、同因子化、同程序，才允许相减。

- 默认展开实际 121→319 模式时，颗粒差分复场的功率加权变化为 **17.37%**，未达到 2% 目标。
- 默认三维背景反射率：121 模式 0.24117476；319 模式 0.25525305；641 模式 0.23548830。加阶未单调趋近 NIST 的 0.27780809，不能宣称收敛。
- 当前旧版内置 normal、vector、Jones 基底也没有使 319 模式背景通过验收；保留为诊断选项，不将名字等同于精度保证。
- 已计算 NA=0.6、倾斜 +45°/-45°/0° 的三路离散级次收集量，保存在报告中，均标为未收敛试算。
- 分别保存 `|E_sample-E_reference|²` 对应的差分场功率与真实强度差。二者不同，后者包含相干干涉项，允许为负。

尚未完成：三维背景的因子化/各向异性基组收敛、球分层收敛、超胞间距收敛、平面周期球对照、同几何 FEM/grcwa 幅相交叉验证及定量成像。有限光栅/AFMM 也没有因此得到认证。

## 复现命令

在项目根目录运行。输出目录必须为空或尚不存在，避免覆盖原始结果。

```powershell
# 同 NIST ±365 阶裸光栅
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_new/bare_731 --one-dimensional --modes 731 --timeout 180

# 同三维基组的 sample/reference 对
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_new/reference_321 --modes 321 --timeout 240
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_new/particle_321 --modes 321 --particle --timeout 240

# 重建本轮报告与对比图
& build/rcwa_env/Scripts/python.exe -m tools.s4_validation_report --plot

# 回归检查
& build/rcwa_env/Scripts/python.exe -m pytest tests/test_s4_prototype.py tests/test_rcwa_prototype.py tests/test_rcwa_grating_fourier.py -q -p no:cacheprovider
```

默认监控进程树占用不超过 8 GiB，保留 3 GiB 系统余量，并限制时间。这是采样式监控，不是操作系统硬内存配额，也不是误差控制。731 模式一维基准本机约 58.5 秒，进程树峰值统计约 0.99 GiB；319 模式颗粒试算约 38 秒、0.58 GiB。不是当前内存不够导致未通过。

## 原始数据

本轮目录：`output/s4_validation_20261022/`。

每个成功算例含 `inputs.json`、`input.lua`、`stdout.txt`、`stderr.txt`、`modes.npz`、`report.json`；外侧有 `.resources.json`。失败尝试也保留输入与日志，但不生成成功结果。

```python
import numpy as np

with np.load('output/s4_validation_20261022/particle_default_321/modes.npz',
             allow_pickle=False) as data:
    orders = data['orders']
    electric = data['electric_reflected_z0']  # complex128, shape (N, 3)
    efficiency = data['reflectance_by_order']
    valid = data['propagating']
```

关键字段：

| 字段 | 含义 |
|---|---|
| `orders` | 全局超胞整数倒格矢；一维结果转换为对应全局级次 |
| `electric_reflected_z0` | 反射 E/E0 复矢量，参考面为晶圆 ridge 顶部 z=0 |
| `directions` | 全局 +Z 向上、+KY 为入射传播投影；非传播级为 NaN |
| `propagating` | 传播级掩码 |
| `reflectance_by_order` | 各级功率/该超胞入射功率，无量纲，不是 W 或连续 DSC |
| `raw_orders` | S4 原生整数级次，未做全局映射 |
| `raw_forward_hxy`, `raw_backward_hxy` | S4 顶层原生模态振幅，在均匀空气中按 Hx/Hy 解析 |
| `raw_flux_by_order` | S4 原生每级正/反向实部及虚部通量 |
| `raw_top_flux`, `raw_bottom_flux` | S4 原生总通量，未做功率归一化 |
| `near_electric`（请求时） | 总近场 E/E0，不是单独的反射场 |
| `near_magnetic_impedance_scaled`（请求时） | 全局 Z0H/E0，总近场 |

近场可通过 Python 的 `solve(..., near_points=[(x,y,z), ...])` 传入任意点，长度单位 μm；不是本轮颗粒算例默认计算的完整近场网格。原始 S4 坐标的近场复数也单独保留。

总报告 `summary.json` 包含 NIST 参考校验和、逐级幅相误差、三路量、截断变化和资源记录；`VALIDATION.md` 是其可读摘要，`nist_s4_validation.png` 是对比图。

## 构建失败记录与后续顺序

本机最初 master 提交构建未通过真空透射检查，保存在 `build/s4_runtime/master_7fd00a2_failed_smoke/`；随后历史提交的外置 LAPACK 路径在图案化特征值求解中访问异常。未把它们视为有效求解器，也不由此断言官方所有平台都存在同样问题。当前采用实际通过前述检查的历史提交与内置 RNP 路径。

下一轮先解决    二维超胞裸光栅的材料因子化与方向性截断，要求每个传播级及复场都达标；再扫描球分层与间距，并做独立交叉验证。只有这些误差可控后才进入 sample/reference 成像与相机 SNR。不把继续加内存作为唯一方案。

官方接口依据：[S4 Lua API](https://web.stanford.edu/group/fan/S4/lua_api.html)、[S4 官方源码](https://github.com/victorliu/S4)。本文件中的数值来自保存的本机算例，不是文档给出的理论保证。
