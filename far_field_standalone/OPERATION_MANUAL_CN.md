# FarField Simulation 中文操作手册

适用版本：`2026.10.10-roughness`，Windows x64 分文件发行版。更新日期：2026-10-10。

本手册按当前界面的英文按钮和字段名称说明操作。SCATMECH、局部 FEM、整光斑 SPM 与混合模型的结果含义不同，请先选对计算入口。本版 FEM/混合功能仍属实验模型，不是完整晶圆的解析解。

## 目录

- [1. 安装与启动](#1-安装与启动)
- [2. 选择计算入口](#2-选择计算入口)
- [3. 坐标和单位](#3-坐标和单位)
- [4. 第一次运行：10 nm 颗粒](#4-第一次运行10-nm-颗粒)
- [5. 照明设置](#5-照明设置)
- [6. 粗糙面设置](#6-粗糙面设置)
- [7. 三路孔径与成像](#7-三路孔径与成像)
- [8. 整光斑 haze 和混合计算](#8-整光斑-haze-和混合计算)
- [9. 多表面平均与 NIST 对照](#9-多表面平均与-nist-对照)
- [10. 收敛与模型检查](#10-收敛与模型检查)
- [11. SCATMECH 页操作](#11-scatmech-页操作)
- [12. 结果保存和重新成像](#12-结果保存和重新成像)
- [13. 常见问题](#13-常见问题)
- [14. 正式结果检查清单](#14-正式结果检查清单)

## 1. 安装与启动

1. 从 [GitHub Release](https://github.com/wufansaber/aptk/releases/tag/v2026.10.10-roughness) 下载 `FarField_Simulation-2026.10.10-roughness-win64.zip`。
2. 完整解压到本地可访问的目录，不要直接在 ZIP 内运行。
3. 双击 `FarField_Simulation.exe`。

主要文件结构如下，实际目录还包含许可、示例和验证摘要：

```text
FarField_Simulation/
  FarField_Simulation.exe
  _internal/
  fem/
    FarField_FEM.exe
    _internal/
  docs/
  verification/
  README_CN.md
  build_manifest.json
```

不要只移动一个 EXE，也不要删除任意 `_internal`。发行版附带界面和 FEM 两套运行时，无需另装 Anaconda、Python、SCATMECH 或 NGSolve。

本手册作为现有 Release 的独立 MD 附件补充上传；此前下载的 ZIP 不会自动更新。后续打包会自动将其收入 `docs/far_field_standalone`。

可用 PowerShell 核对 ZIP 的 SHA-256，与 Release 中同名 `.sha256` 文件比较：

```powershell
Get-FileHash .\FarField_Simulation-2026.10.10-roughness-win64.zip -Algorithm SHA256
```

### 运行环境

界面建议使用至少 1280×800 的显示区域。左侧参数可滚动，也可拖动中间分隔条调整宽度；结果页右侧箭头可滚动标签，方框图标用于适配图片。

FEM 页默认 `Memory limit = 10 GiB`、`Time limit = 30 min`、2 个 FEM 线程。内存限制是作业预算，不代表计算机总内存；保护程序还默认保留约 3 GiB 系统可用内存，且是周期采样保护，不是操作系统硬上限。

复杂混合/收敛任务建议 32 GB RAM、足够磁盘空间，并关闭占用内存较大的应用。不能仅因计算慢就无限提高线程或放开内存限制；时间和内存需求随网格、阶数、表面范围及样本数变化。

运行中可点击 `Cancel` 请求停止，等待程序确认子进程退出。触发取消、超时或内存保护后，不要将目录内已生成的中间图片视为成功结果。

### 从源码启动（开发环境）

在项目根目录，可通过当前电脑的 NIST2 解释器直接启动界面，不依赖 PowerShell 的 conda 激活脚本：

```powershell
& 'D:\anaconda\envs\NIST2\python.exe' .\far_field_standalone\main.py
```

源码界面依赖 NIST2 中的 PyQt5/pySCATMECH，FEM 求解器使用另一个包含 NGSolve 的隔离环境，当前路径为 `build/ngsolve_env/Scripts/python.exe`。迁移电脑时需要重新配置这些环境和求解器路径；发行版 EXE 不需要这些开发环境。

## 2. 选择计算入口

顶部有 `SCATMECH` 和 `Full-wave FEM (experimental)` 两页。FEM 页的不同按钮不是同一模型的显示选项。

| 需求 | 入口 | 计算范围与注意事项 |
| --- | --- | --- |
| 球形颗粒远场、镀膜、多层粗糙度、粒径扫描 | SCATMECH 页 | NIST 模型；不是近场 FEM |
| 颗粒与局部粗糙斑块的近场、远场、像面 | FEM 页 `Calculate` | 局部 FEM 域；斑块外仍是平整晶圆 |
| 平滑晶圆上的颗粒与局部纯 haze 独立对照 | `Compare particle / haze` | 先加载完整局部结果，必要时补算平滑晶圆颗粒源 |
| 多张局部随机粗糙面的 FEM 强度平均 | `Average haze` | 无颗粒；不是整光斑全局 FEM |
| 整个有限光斑覆盖粗糙面的纯 haze | `Whole-beam haze (SPM)` | 相干照明、一阶表面微扰，无颗粒 |
| 颗粒与整光斑粗糙面共同计算 | `Particle + whole haze (hybrid)` | 整面 SPM + 共享局部 FEM − 相同局部 SPM |
| 使用保存的混合复场重新成像 | `Re-image hybrid` | 不改变原始表面、颗粒或照明形状 |

如果目标是“颗粒 + 整个光斑下的粗糙背景”，应选择 hybrid，而不是认为 `Calculate` 已把局部粗糙斑块扩展到了整个照明区域。

## 3. 坐标和单位

- 晶圆名义表面为 XY 平面，+Z 指向空气，入射角相对表面法线定义。
- 光源位于 −KY 侧，入射传播方向在表面上的投影为 +KY。入射方向可写成 `(0, sinθ, −cosθ)`。
- 三路物镜轴在 XZ 平面倾斜；正倾角朝 +KX，负倾角朝 −KX，0° 沿 +Z。
- 远场二维全局坐标为 `kx/k`、`ky/k`，圆盘半径对应 `sin(散射极角)`，不是晶圆位置或像面位置。
- 同一各向同性、平滑或统计对称模型下，±X 两路可能对称；单张随机粗糙面不必逐像素对称。

| 量 | 单位/含义 |
| --- | --- |
| 波长、颗粒直径、RMS、像差 | nm |
| 光斑、局部粗糙面宽度、相关长度、像面坐标 | μm；界面写作 `um` |
| 物镜焦距 | mm |
| 激光功率 | W |
| 远场 DSC | μm²/sr；半球或孔径积分后为 μm² |
| 像面照度 | W/μm² |
| 收集或像面总功率 | W，图中可能换算为 pW/nW |
| 近场图 | 相对参考场的 `|E/E_inc|²`，不是绝对功率图 |

`n` 为折射率实部，`k` 为消光系数，复折射率是 `n+ik`；这里的 `k` 不是介电常数的虚部。更改波长不会自动取得该波长的材料色散，必须填写对应的 n/k。

## 4. 第一次运行：10 nm 颗粒

### 4.1 先验证纯颗粒流程

1. 打开 FEM 页，选择 `Parameter set = 10 nm PSL preview`，点击 `Apply preset`。
2. 确认 `Particle diameter = 10 nm`。预设的部分关键条件见下表。
3. 首次检查可将 `Roughness RMS = 0 nm`，先不启用 `Camera response`，以排除粗糙背景和探测器噪声。
4. 保持程序识别的 `fem/FarField_FEM.exe` 求解器路径，点击 `Calculate`。
5. 等待完成提示，查看 `Near field`、`Far field`、`Three image planes` 和 `Result parameters`。
6. 需要局部粗糙背景时，再把 RMS 改为目标值并重新 `Calculate`；需要整光斑背景时按第 8 节运行 hybrid。

| 参数 | 当前 10 nm 预设/默认条件 |
| --- | --- |
| Wavelength | 355 nm |
| Particle n / k | 1.59 / 0 |
| Wafer n / k | 5.57 / 0.387 |
| Incidence / normal | 75°，p 偏振 |
| Sphere / wafer gap | 10 nm |
| Illumination / Finite beam definition | `gaussian` / `wafer_gaussian` |
| Laser power / wafer spot | 1 W / 2×2 μm 目标强度直径 |
| Rough surface model / RMS / PSD length | Gaussian / 2 nm / 0.05 μm；纯颗粒检查手动改 RMS=0 |
| Rough surface nodes / seed | 9 点/轴 / 7 |
| Channel 1/2/3 NA | 0.5 / 0.5 / 0.5 |
| Channel 1/2/3 tilt | +45° / −45° / 0° |
| Magnification / objective focal length | 20 / 20 mm |
| FEM profile / order / direct solver | preview / 1 / umfpack |
| Global / local mesh scale | 1.50 / 1.00 |
| Particle mesh / source quadrature | 0.00150 μm / 7 |
| Pupil samples / beam quadrature | 128 / 31 |

这是低内存预览，不是已对任意参数证明收敛的高精度配置。手动修改后预设显示 `Custom`；恢复预设会覆盖相关输入，请先记录自定义参数。

### 4.2 改为 266 nm、单路 NA=0.96

先恢复并记录预设，再设 `Wavelength = 266 nm`，填写该波长的颗粒和晶圆 n/k。将不用的两路 NA 设为 0，将保留通道 NA 设为 0.96、倾角设为所需方向，然后重新计算。

高 NA 改变收集范围和像面采样要求，应检查光瞳采样、像场大小和镜面方向遮挡。程序的散射图没有自动加入真实暗场挡光器，也不能用固定 355 nm 材料预设直接声称得到真实 266 nm 实验结果。

## 5. 照明设置

`Illumination = gaussian` 是兼容命名，当前也用于矩形/实测有限包络；实际形状由 `Finite beam definition` 决定。`plane` 是平面波，与有限 2×2 μm 光束不是同一条件。

| Finite beam definition | 长宽定义 | 操作要点 |
| --- | --- | --- |
| `wafer_gaussian` | 晶圆面目标强度 1/e² 直径；length 沿 Y、width 沿 X | 需要晶圆上的 2×2 μm 光斑时选此项 |
| `beam_cross_section` | 垂直光轴的高斯截面 1/e² 直径 | 晶圆面几何 Y 长度约为 length/cosθ，不是仍保持 2 μm |
| `wafer_rectangle` | 晶圆面的完全相干均匀矩形目标包络 | 传播谱截断后不再是严格硬边矩形 |
| `measured_wafer` | 晶圆面实测强度和可选包络相位 | 用文件选择按钮选 NPZ，检查缺失相位策略 |

四种有限模式均设 `Beam waist = 0`，不要同时填写旧的非零 waist。激光总功率独立设置，默认 1 W。

有限单频照明使用相干角谱。掠入射时目标谱的一部分可能位于传播圆盘外，实际生成光斑会改变；查看 `Illumination` 中的实际分布、功率和覆盖报告。排除的标量谱权重不是激光功率损失率，传播场仍按设定总功率归一化。

### 实测光斑文件

NPZ 必须包含以下数组，示意代码只说明格式，变量应由真实测量提供：

```python
import numpy as np

np.savez("beam.npz",
         x_um=x, y_um=y,
         intensity=I,
         phase_envelope_rad=phase)
```

- `x_um`、`y_um` 是均匀递增的像素中心坐标，单位 μm，每轴 4～256 个点。
- `intensity` 形状为 `[Y, X]`，必须有限、非负且不全为零；允许相对强度，功率由 Laser power 指定。
- `phase_envelope_rad` 可省略，单位 rad，表示不含入射 +KY 载波的包络相位。
- 缺少相位时，`flat` 明确假设平坦包络相位；`reject` 拒绝运行。强度不能唯一确定相位。
- 测量范围外按零处理，不自动补出未知光斑尾部；文件改变后必须重新求解。

矩形/实测输入通常需要较高 `Beam quadrature order`，界面切换时至少设为 61。局部场积分检查不通过时逐级增加到 81、101、151，并检查报告；局部阶数通过不保证整光斑积分收敛。

## 6. 粗糙面设置

| 字段 | 含义及注意事项 |
| --- | --- |
| Rough surface model | 新计算优先用 Gaussian 或 ABC PSD；legacy 是历史高度生成器 |
| Roughness RMS | 平稳表面高度 RMS，单位 nm；局部边缘渐隐后实际 RMS 可不同 |
| PSD length B / L | Gaussian 为相关长度 L，满足 `C(r)=σ² exp[-(r/L)²]`；ABC 的 B 不是同一参数 |
| ABC exponent C | 高频衰减；必须大于 2。C≤4 时无限带宽 RMS 斜率不有限，必须关注物理高频截止 |
| Rough patch width | 局部 FEM 粗糙斑块边长，不等于光斑长宽或整面规划宽度 |
| Rough surface nodes / axis | 局部表面节点数；Gaussian 至少需要每个 L 两个网格间隔，再做加密检查 |
| Roughness seed | 固定种子便于复现；修改种子得到另一张表面，但须保持其他条件相同 |
| Edge taper | 物理边缘渐隐宽度，单位 μm；需覆盖至少一个表面网格间隔 |

RMS=0 可用于平滑晶圆颗粒检查。二维 haze 图含散斑、偏斜或局部热点，不代表网格一定不对称；统计 PSD 描述随机表面的集合，不要求每个 realization 与 NIST 平均图完全相同。

### 实测粗糙输入

在 `Rough surface model` 选择 `measured_height` 或 `measured_psd`，然后选择对应测量文件。RMS 由文件读取；仪器传递函数和高频补偿不会自动完成。

| 输入 | 必需字段或格式 |
| --- | --- |
| 高度 NPZ | 二维 `height_nm[Y,X]`，正数 `dx_um`、`dy_um` |
| 高度 CSV | 无标题二维 nm 高度矩阵；同名 JSON 提供 `dx_um`、`dy_um` |
| 二维 PSD NPZ | `fx_per_um`、`fy_per_um`、`psd_um4[Y,X]`；双边、含零频的均匀 DFT 网格，满足反演对称 |
| 径向 PSD CSV | 两列 `f_per_um,S2_um4`；频率递增、非负；不要混用一维 PSD |

PSD 频率单位是 cycles/μm，不是 rad/μm；二维 PSD 单位是 μm⁴。高度可选 `none/mean/plane` 去趋势；缺失点可选拒绝或受限的 nearest 填充。

当前 `Whole-beam haze (SPM)` 和 hybrid 支持 Gaussian/ABC PSD，不支持把实测高度或 PSD 直接作为整面输入。实测高度是一张确定性表面，不能当作多张独立样本重复计数；实测 PSD 可以用于局部 `Average haze`。

## 7. 三路孔径与成像

每一路独立设置 `Channel N NA` 和 `Channel N tilt`。NA=0 关闭该路；关闭通道后，若仍模拟相机，暗电流/读出噪声不一定为零。倾角不是 NA，也不是扫描像面的位移。

默认 NA=0.5，对空气物镜对应约 30° 光锥半角。三路默认倾角 +45°/−45°/0°，分别传播到独立像面，不相干合束。

| 参数 | 操作建议 |
| --- | --- |
| Magnification / Objective focal length | 默认 20×、20 mm，对应理想管镜焦距 400 mm |
| Pupil samples | 每轴光瞳采样；逐步增加检查三路功率和图像，不是 FEM 网格参数 |
| Image padding | 改善像面离散采样；不能替代增加 Pupil samples 来扩大物理计算像场 |
| Image defocus、Astigmatism、Coma、Spherical | 像面离焦/像差；可对保存源重新成像 |
| Power transmission | 成像功率透过率，不是额外激光功率 |
| Image archive | compact 保存重成像所需源及必要结果；full 额外保留像面复 E/H，文件更大 |
| Camera response | 启用电子数及噪声读出；需配置像素、ROI、曝光、量子效率和噪声 |

三路总功率应读取报告，不要仅依据亮度判断。FEM 三路总场像面共用绝对色标；某些分解/对照列使用不同的明确标注色标。SCATMECH 三路成像也应按各自色标和数字比较，不能把颜色相同解释为功率相同。

`Coherent decomposition` 的颗粒变化场是同一粗糙面有/无颗粒的复场差，包含耦合变化和干涉项；它不是平滑晶圆上的独立颗粒实验。`Compare particle / haze` 才是“平滑晶圆颗粒”与“保存粗糙面无颗粒”的独立对照；其功率不能相加来重建颗粒位于粗糙面上的总结果。

相机的 `Ideal SNR` 是指定噪声假设下的指标，不代表实际设备检测率。额外均匀相机背景不等同于物理 haze。注意满阱饱和、ROI 越界和实际波长处的量子效率。

## 8. 整光斑 haze 和混合计算

### 8.1 纯 haze：Whole-beam haze (SPM)

1. 设置有限照明定义、功率、波长、入射角、偏振和晶圆 n/k，保持旧 waist=0。
2. 选择 Gaussian 或 ABC PSD，设置非零 RMS、相关长度/指数、种子。
3. 设置三路 NA/倾角、倍率和光瞳采样。
4. `Surfaces` 可先设 1 作流程检查；统计分析再选固定的较大样本数。
5. 保持 `Whole-beam power coverage = 0.99`、`Whole-beam PSD variance fraction = 0.90` 作为起始规划条件，检查 Edge taper。
6. 点击 `Whole-beam haze (SPM)`。完成后查看整面高度/照明、`Mean far field`、`Mean image planes` 和 `Model / coverage / power`。

颗粒尺寸不参与此纯 haze 模式。一阶近似要求 RMS/波长足够小；程序限制 RMS/波长<0.05 不等于满足全部小斜率及数值收敛条件。

光斑 2×2 μm 不代表粗糙面规划也应是 2×2 μm。软件根据实际传播光束、照明尾部和边缘渐隐选择更大的计算窗口，读取报告中的实际宽度，不要把过去的 7.5 μm 当作固定规则。

### 8.2 颗粒 + 整面：Particle + whole haze (hybrid)

1. 在上述照明和 PSD 条件上，填写颗粒直径、n/k、间隙及局部 FEM 网格。
2. 设置 `Hybrid mother range factor`，默认 8，可选 2/4/8；先确认内存预算。
3. 点击 `Particle + whole haze (hybrid)`，等待共享母表面、有/无颗粒局部 FEM、整面辐射和成像全部完成。
4. 检查 `Illumination`、局部 XY/XZ/YZ、远场、三路像面、相干分解以及保存参数。

混合在复振幅层面组合并减去重复局部粗糙散射。扩大母表面因子会改变新的随机实现，需要新求解局部 FEM，不能把旧 FEM 源直接拼接到新母表面。

hybrid 近场仍仅覆盖局部 FEM 域，不是整片粗糙晶圆的近场。外部粗糙面与颗粒之间的全局高阶往返耦合尚未求解；不要把这一路径称为整个光斑的全波 FEM。

## 9. 多表面平均与 NIST 对照

### 9.1 局部 FEM 平均：Average haze

1. 选择 Gaussian/ABC 或实测 PSD，RMS>0，设局部宽度、节点数和种子。
2. 若需同时显示 NIST 平面波 PSD 对照，设 `Illumination = plane`、waist=0。
3. 设置 `Surfaces`、`Edge taper` 和 `Compare area`。0.2×0.2 μm 斑块的名义面积是 0.04 μm²；比较面积不是有限激光光斑面积，也不是自动复原渐隐损失的系数。
4. 点击 `Average haze`，完成后查看 `Haze ensemble`、统计指标和有限斑块参考。

有限照明也可用于局部 FEM 平均，但该模式不会自动生成平面波 NIST 对照。不同表面平均的是强度/功率，不是复振幅。增加样本只能降低统计波动，不能消除小斑块、渐隐、网格、PSD 截断或模型误差。

### 9.2 整面 SPM 扩展到 128 张

1. 先运行或用 `Load result` 加载完整的 whole-beam haze 结果。
2. 打开 `Model checks`，选择 `Extend whole haze to 128 surfaces`。
3. 软件锁定已保存的物理输入和原规划，复用原样本并追加种子，写入新结果目录；左侧临时编辑值不会替换旧源。
4. 查看逐指标 95% 区间、同时区间和目标是否满足。已有 128 张时无需重复扩展。

128 不是所有输入都能达到 5% 精度的保证。发布版验证算例的三路同时区间半宽约为 3.24%、3.80%、3.18%，仅适用于那组 haze-only SPM 参数；不包括表面范围、FEM 或高阶耦合误差。正式统计应预先确定样本数，不因反复观察区间而提前停止。

### 9.3 SCATMECH 页的输入对齐

`Match displayed FEM haze result` 读取已完成且当前显示的局部 FEM 保存结果，不读取可能被修改的表单。用于平面波、非零 Gaussian/ABC 粗糙度、均匀晶圆条件；需要先移除 SCATMECH 镀膜。

点击后仍须在 SCATMECH 页重新运行远场。此操作不会把有限斑块 realization、有限相干光斑或 hybrid 转换成相同的无限统计 PSD 模型；不应要求两张图逐像素一致。

## 10. 收敛与模型检查

### Model checks

先加载完整 hybrid 源，才可启用相应菜单；仅有 PNG/报告的发布示例不能做需要复场的检查。

| 菜单项 | 含义 |
| --- | --- |
| Whole-surface range | 同母表面扩大外部积分范围，检查复场、功率和积分误差 |
| Same-surface grid / taper | 同一频带表面改变空间采样和渐隐；不自动增加 PSD 高频模式 |
| Low-height absolute power / energy balance | 缩小粗糙高度，比较 FEM/SPM 绝对功率及控制体内总场通量/吸收 |
| Low-height FEM / SPM (P1) | 一阶局部 FEM 的低高度极限检查 |
| Thin-limit FEM (P2, sparse) | 独立 P2 检查，可能需要更多内存 |

这些操作另存结果，不覆盖生产源。范围审计通过也不会自动将大窗口替换到原像面。

`Convergence` 用于局部 haze 的多档参数研究。对话框可填写逗号分隔的比较档位；空行表示跳过该维度。设置每工况样本数、数值容差及统计区间目标后运行，查看配对误差、功率、区间和资源状态。

**99% 入射功率覆盖率 ≠ 散射结果误差<1%；低求解残差 ≠ 物理模型正确；增加 FFT padding ≠ 表面/FEM 网格加密。** 正式结论应同时考虑表面范围与频带、FEM 网格/阶数、PML、角度积分、光瞳和统计误差。

## 11. SCATMECH 页操作

### 远场与三路物镜

1. 设置 Wafer material、晶圆 n/k；有膜层时设置 coating stack。表格第一行靠近基底，最后一行靠近空气。
2. 设置 Particle material、Sphere diameter、颗粒 n/k、激光波长、入射角和偏振。
3. 要纯颗粒时关闭 `Include roughness background`；要纯 PSD haze 时开启它并勾选 `Haze only (no particle)`。
4. 配置三路独立 NA/倾角，点击 `Run far-field simulation`。
5. 点击 `Run ideal objective example` 查看物镜光锥和出瞳；颗粒或照明参数改变后须先重新运行远场。
6. 平滑晶圆颗粒成像使用 `Calculate three image planes`；当前这一路径须关闭粗糙度，不能用它替代 FEM/hybrid haze 成像。

`Volume` 是光锥在单位上半球中覆盖的几何体积占比，`Energy` 是按实际 DSC 加权的能量占比。三路合并占比按方向掩膜并集积分，重叠只计一次；独立探测器各自的功率不能当作这个去重并集占比。

SCATMECH 颗粒+PSD 背景按非相干强度合成，不包含同一 realization 的颗粒/粗糙面相干耦合。SCATMECH 可计算膜层；当前 FEM/hybrid 的晶圆仍是均匀半空间，两页不能因输入看起来相近就直接比较多层膜结果。

### 粒径序列

设置 `Start diameter`、`End diameter`、`Diameter step`，点击 `Run particle-size sweep`。例如 30～200 nm、步长 10 nm，可查看三路积分 DSC 和能量占比曲线；NA=0 的通道保持关闭。

粒径扫描使用序列中的直径，不使用单个 Sphere diameter。它是 SCATMECH 扫描，不是批量 FEM/hybrid 全流程扫描。启用粗糙背景时，小颗粒曲线可能出现背景平台，应同时检查纯颗粒结果。

## 12. 结果保存和重新成像

发行版 FEM 结果自动保存在：

```text
%LOCALAPPDATA%\FarField_Simulation\results\时间戳目录
```

可把该路径粘贴到资源管理器地址栏，或按界面日志中的输出路径定位。源码运行通常保存在项目的 `build/fullwave_ui`。若设置 `FARFIELD_RESULTS_DIR`，以该环境变量指定的路径为准。

| 文件/目录 | 用途 |
| --- | --- |
| request.json、resources.json | 实际启动命令、作业内存/时间及完成状态 |
| scattering/report.json、paired_fields.npz | 局部保存源的参数和复场；NPZ 位于 scattering 目录 |
| imaging_report.json、image_fields.npz | 局部成像参数、各路总功率及像面数组 |
| near_field_xz/xy/yz.png | 局部近场；XY 位于颗粒中心高度，XZ 固定 Y=0，YZ 固定 X=0 |
| image_channel_1/2/3.png | 三路独立像面 |
| hybrid_source_report.json、hybrid_fields.npz、hybrid_report.json | 混合源、来源签名和混合成像报告 |
| whole_haze_report.json、whole_haze_fields.npz | 整面 SPM 统计、实际规划、功率和数组 |
| haze_ensemble_report.json | 局部 FEM 平均及参考比较 |

文件随计算模式不同，不会每个目录都包含上表全部文件。近场范围是局部求解范围；二维 PNG 没有完整三维复相位，不能代替散射源。

用 `Load result` 选择结果目录。普通完整局部结果可用 `Image current result` 或 `Re-image source`；混合结果使用 `Re-image hybrid`。whole-beam SPM 参数修改后通过原计算入口另算，不应选普通局部重成像按钮。

| 修改内容 | 后续操作 |
| --- | --- |
| 颗粒尺寸/n/k/间隙、晶圆材料、波长、入射角/偏振 | 重新运行相应物理求解 |
| 粗糙 PSD/RMS/相关长度/种子/测量文件、母域/局部几何 | 重新运行相应物理求解 |
| 光斑定义、长宽、实测强度/相位 | 重新求解；旧复场不是新光斑下的散射源 |
| NA、倾角、倍率、离焦、像差、透过率、相机设置 | 对支持的保存源重新成像 |
| 仅激光功率，保持照明形状及物理条件不变 | 对支持的保存源重新成像按功率缩放；不是重算材料非线性 |
| 仅改变色标或视图 | 不改变原始 DSC 和积分功率 |

搬运可重算的结果时复制整个目录及其实际依赖源，尤其不要删除局部 paired_fields、混合复场或测量输入。部分报告引用原绝对路径，移机后需要保留对应依赖或重新求解；只有图像的示例可用于查看，但不能补算缺失近场。

SCATMECH 页用 `Export figures and data` 导出 PNG、NPZ 和参数；粒径扫描还导出 CSV。每个工况使用独立目录，避免覆盖同名文件。

## 13. 常见问题

| 现象 | 检查与处理 |
| --- | --- |
| 改尺寸后结果不变 | 是否只查看旧图/重成像旧源？重新物理求解，检查 Result parameters 中已保存的直径 |
| 点 Calculate 后没结果 | 查看底部日志及 resources.json；确认求解器路径、输入检查、内存和超时，未完成任务不当作有效结果 |
| Gaussian footprint/legacy waist 报错 | 有限光斑定义下 Beam waist 必须为 0，使用长宽配置光斑 |
| 高阶光束积分未通过 | 逐级提高 Beam quadrature order；检查实际 Illumination，不要忽略角谱或覆盖积分警告 |
| Gaussian 粗糙面无法运行 | 检查节点间距与 L、局部宽度是否适配 FEM 域、求解器及预算；增加节点可能显著增加内存 |
| FEM 和 NIST 图不同 | 区分单表面/集合均值、有限/无限斑块、有限光束/平面波、实际渐隐面积和绝对色标；先统一保存的输入 |
| 中路暗或像环形 | 检查向量偏振、物镜方向、色标及光瞳/FEM 收敛；不是所有颗粒像都必须为实心圆，矩形缺块仍需检查采样 |
| 小功率时图很暗 | 比较色标和总功率而非颜色；功率按线性模型缩放，不应为显示而改变物理归一化 |
| Model checks 灰色 | 先加载对应完整 hybrid 或 whole-beam 源；不同源只启用支持的动作 |
| ROI 越界 | 减小 ROI/像素尺寸或增大 Pupil samples；padding 不能替代扩大物理像场 |
| 触发资源保护/取消 | 保存日志，降低预览负载或增加合理预算，重新运行；不要把残留图片当作完成结果 |
| DLL 缺失 | 重新完整解压，检查主程序和 fem 两套 _internal，确认文件未被隔离 |

## 14. 正式结果检查清单

1. 记录软件版本、计算入口和完整参数，保存原始测量文件与标定单位。
2. 确认波长处的材料 n/k、光斑定义、入射方向、功率、三路 NA/倾角和单位。
3. 先做纯颗粒/纯 haze 检查，区分独立实验、相干变化场和总场。
4. 确认资源报告完成，核对 warnings、生产窗口、表面频带及局部近场范围。
5. 独立检查 FEM/PML、表面范围/网格/渐隐、角度积分和光瞳采样，不能只看残差。
6. 统计分析预先确定样本数，报告逐指标及同时区间；不把统计区间当成全部模型误差。
7. 从保存报告读取绝对功率；DSC 不能直接乘激光总功率来得到 W，应使用对应参考辐照度和程序归一化。
8. 保存可复用复场及依赖源。高阶全局耦合、多层 FEM、镜面泄漏和未知测量相位等未建模因素须注明。

更多说明：[模型实现与本次校验](../fullwave_prototype/ROUGHNESS_CAMPAIGN_CN.md)、[混合模型](../fullwave_prototype/HYBRID_MODEL_CN.md)、[整光斑 haze](../fullwave_prototype/WHOLE_BEAM_HAZE_CN.md)、[原有详细指南](USER_GUIDE_CN.md)。仓库内这些链接相对有效；单独下载本手册时，可在 GitHub 同目录查看对应文件。
