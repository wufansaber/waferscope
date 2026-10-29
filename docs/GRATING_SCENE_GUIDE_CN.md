# S4 周期光栅场景操作与模型说明

本模块将同一光栅的无颗粒参考场和含颗粒样品场，连接到三平面近场、离散远场、三路独立收光、理想相干成像及模拟面阵相机。它使用真正的 S4 原生程序，不以 grcwa 输出替代 S4。

**这是周期超胞模型，不是已验证的孤立缺陷或有限光栅模型。Images 是理想晶圆共轭 Jones 中继；Camera DN 是该模型经过像元积分和噪声后的模拟灰度，不是实机测量。** 原有 FEM 和 SCATMECH 页面保留。

本页电磁计算只包含确定性的光栅和球，不自动加入统计粗糙度/haze 或多层膜。相机计算可单独加入传感器噪声，但不会把均匀电子背景伪装成晶圆 haze。

## 运行与操作

在项目根目录启动现有程序，选择新增的 `S4 Grating (periodic)` 页签：

```powershell
& D:/anaconda/envs/NIST2/python.exe far_field_standalone/main.py
```

Qt 界面运行在 NIST2；原生求解由独立的 `build/rcwa_env/Scripts/python.exe` 调用 `build/s4_runtime/bin_stripes/S4.exe`。两个路径可在页面上指定，不需要向 NIST2 安装另一套原生求解库。

1. 设置波长、入射角、光栅周期/线宽/深度、超胞尺寸及材料 n、k。光栅沿 X 不变，沿 Y 周期；光源来自 -KY，传播投影沿 +KY。
2. 设置颗粒直径和底部间隙。间隙默认 **0 nm**，球底接触脊顶。球体采用垂直分层近似，切片数仍需收敛；默认材料数值不是经 266 nm 标定的 PSL/晶圆数据。
3. 各路分别设置 NA、倾角；NA=0 关闭该路。倾角绕 Y 轴变化，正/负值分别朝向正/负 KX。
4. `Calculate` 重新求解参考和样品。`Preview` 是快速预览，不是已通过精度验收的配置。`Higher cutoff` 增加阶数和切片，也不自动代表收敛。
5. 在结果中查看 Geometry、Near XY/XZ/YZ、Far field 2D/3D、Images 和 Collected powers。较小窗口可滚动图像；图中数据不按峰值归一化，`Log color scale` 仅改变显示色标。
6. 只改 NA、倾角、倍率、分析器、离焦、透过率或功率归一化时，可用 `Update optics`。改变几何、材料、偏振或近场观察面必须重新求解；程序会拒绝复用不匹配的场。
7. `Load` 加载已有完整结果；`Export raw bundle` 导出完整 ZIP，包含原始 S4 输入、输出和复场，不只是截图。`Cancel` 终止本次求解拥有的工作进程和原生子进程，保留之前加载的结果。
8. 侧栏下部设置面阵相机参数，点击工具栏的显示器图标，从当前**已保存**光学结果计算相机，不重新求解 S4。Camera DN 显示三路 sample/reference/差值，Camera metrics 显示原始 ADC 帧计算的 ROI SNR。计算相机前，表单中的光学修改必须先用 Update optics 应用。
9. Diagnostics 列出 R、T、1-R-T、横向场、X 镜像复场误差和两侧功率不对称度。1-R-T 是通量账本，不是独立体积分得到的吸收验证；镜像诊断不强制修改场。

2026-10-28 分文件版在页面顶部提供三组 **sphere + grating** 示例。`Apply parameters` 仅填入参数，
`Open example` 加载完整的原生近/远/像面或相机记录，无需等待重算。EXE 自带 `fem` 工作进程与 `s4` 原生程序，
不依赖开发电脑上的 Anaconda；请保留整个文件夹。混合场景的独立实现对照见随包
`verification/grating_hybrid_validation/HYBRID_VALIDATION.md`。当前示例仍标记为实验性、未完全收敛。

加载结果后，工具栏的详细检查图标可运行当前结果的 X/Y 截断、切片和联合加密对照。结果列在 `Validation` 页，保留原场景，不会把局部检查通过改称为孤立极限通过。此入口针对已保存结果，不采用可能已改动的表单物理参数；时间上限按各对照算例计算。

示例输入为 [接触球预览配置](../examples/grating_contact_preview.json)。命令行等价调用：

```powershell
& build/rcwa_env/Scripts/python.exe -m grating_scene.run --config examples/grating_contact_preview.json --output output/my_grating_scene --timeout 1800
```

输出目录必须为空。计算中断的目录没有完整 `scene_report.json`，不能作为成功结果加载。资源保护保留 3 GiB 系统可用内存，进程树统计默认上限 8 GiB；这是采样监控，不是操作系统硬限制。

`Sample native near planes` 默认开启。只检查远场/成像收敛时可以关闭，报告将明确标记未采样近场，三个近场页签也会禁用；不会拿旧近场充当新结果。

## 功率与照明定义

当前电磁激励是完全相干的 **平面波**。参考面积上的功率输入只用来确定幅度，不会产生有限矩形或高斯光斑，也没有改变原先暂缓的高斯光束建模。

设参考晶圆面面积为 `Aref=Lref*Wref`，其上的入射功率为 `Pref`，入射角为 theta，超胞面积为 `Acell=Lx*Ly`：

```text
I_beam_normal = Pref / (Aref * cos(theta))             [W/um^2]
P_incident_cell = I_beam_normal * Acell * cos(theta)   [W]
                = Pref * Acell / Aref
P_collected = P_incident_cell * sum_selected(R_g)
```

默认 1 W、2x2 um 是上述参考晶圆面归一化，不是仿真了边缘有限的 2x2 um 光斑。S4 单位平面波的斜入射功率带有 cos(theta) 因子，见 [S4 单位与功率约定](https://web.stanford.edu/group/fan/S4/units.html)。

三路收集的是各自光锥内的传播衍射级，每路可以有重叠。联合收光用级次集合的并集，只计一次；不能把三路的和当作去重后的总功率。

## 近场和远场

XY 平面高度由 `XY observation z` 指定；XZ 在 y=0，YZ 在 x=0。坐标全局 z 向上，z=0 为未刻蚀脊顶。近场来自 S4 的原生 `GetFields`，包含入射和散射的总 E、H；参见 [S4 原生接口](https://web.stanford.edu/group/fan/S4/lua_api.html)。

近场图显示 `|E/E0|^2`。差场图是 `|(Es-Er)/E0|^2`，不是局部电磁能量密度，也不是近场能流之差；介质中的能量密度还涉及材料和色散。磁场保存为阻抗缩放的复幅。

周期远场只有离散衍射级，图中的点不会被插值成连续 DSC 热图。全局坐标为 kx/k、ky/k、kz/k，颜色为每级反射效率。3D 图中的光锥表示收集方向，并非人为归一化后的能量体积。

## 相干成像定义

像面采用理想晶圆共轭、平移不变的 Jones 中继。倾角决定收集光锥和偏振输运，坐标按晶圆面统一倍率 M 标定；这不等于未作几何校正的斜物镜相机平面，没有自动加入 Scheimpflug 校正、实际镜组或像差。

在被收集的传播级中，将偏振无损输运到对应通道的横向基底，得到二分量场 J_g。对每级应用：

```text
C_g = sqrt(T * n_g,z) / M * Analyzer(J_g) * exp(i*defocus_phase_g)
Eimage(x,y) = sum_g C_g * exp(-2*pi*i*(gx*x/(M*Lx)+gy*y/(M*Ly)))
Iimage = I_beam_normal * sum_components(|Eimage|^2)
```

保留衍射级之间的原始复相位；只省略对强度无影响的共同 Bloch 载波和共同相位。以 NumPy FFT 实现 Fourier 求和，并用独立直接求和测试核对。相干成像必须先叠加复幅再取强度，相关基础见 [TU Delft 相干成像频率响应](https://qiweb.tudelft.nl/aoi/coherentimaging/coherentimaging/)。这里的具体 Jones 输运和功率权重属于本项目明确选择的理想模型，不由该文献直接认证为实际斜物镜。

每个放大超胞内按 Parseval 检查像面总功率与传输后的模态功率一致。倍率改变面积和强度密度，不凭空改变总功率。这个检验只验证所定义中继的代数和归一化，不验证真实镜组。

必须区别：

```text
Delta I = I_sample - I_reference
I_delta_field = |E_sample - E_reference|^2 * I_beam_normal
Delta I = I_delta_field + interference
```

`Images` 第三列为真实强度相减，可为负。中路接收光栅背景衍射时，干涉项可能很显著；不能用差场强度替换。相机可计算条件性灰度/SNR，但光学模型没有因此获得定量认证。

若边路只收集一个 kx 级次族，像面可能沿 X 呈条带。这是周期模型和离散收光的结果，不是孤立球必然成条纹，也不能靠插值或按通道增亮改成点像。需要增大超胞或验证有限边界模型。

## 面阵相机和灰度 SNR

相机不从截图或经过显示拉伸的像面取样，而是对上节的复级次计算强度 Fourier 自相关。对方形有效像元积分时，每个强度频率乘以相应的 sinc 因子。有效边长为 `pixel_pitch*sqrt(fill_factor)`；二维面积积分与独立 Gauss 积分、一个完整放大超胞的绝对功率均有测试。这避免粗糙图像网格造成的相机功率误差。

像面是周期函数。相机 ROI 超过 `M*Lx`、`M*Ly` 会接收到**模型中真实重复的颗粒和光栅**，不是增加了孤立缺陷视野。Camera metrics 的 Periodic copies X/Y 明确列出覆盖的周期数，不能把多胞 SNR 当作单颗粒 SNR。

像元功率按曝光和光子能量转为光子数，经 QE、光子/暗电流 Poisson 噪声、读出噪声、固定 PRNU/DSNU/DCNU、满阱和线性 ADC，输出 uint16 容器内的原始 DN。复用现有相机核心，定义参考 [EMVA 1288 Release 4.0 Linear](https://www.emva.org/wp-content/uploads/EMVA1288Linear_4.0Release.pdf)；这是模型计算，不是 EMVA 合规测量或认证。没有 TDI 累加。

参考帧是同一平面波下的**无颗粒光栅**，不是 haze。每组独立采集一个 sample 和 N 个 reference，计算原始 DN 的 ROI 和：

```text
D_j = sum_ROI(sample_DN_j - mean_N(reference_DN_j))
SNR_raw = abs(mean_j(D_j)) / std_j(D_j, ddof=1)
```

`Raw SNR pairs=0` 不采集；否则支持 100~20000 组，受像元总数保护。存在满阱、ADC 触轨或不可分辨噪声时明确标记无效，不输出可信 SNR。这里是同一静态模型条件下的时间统计，不包含表面变化、运动、配准误差、真实 haze 或测量标定误差。零 NA 仍可有相机暗噪声，但没有光学信号。

Gemi 预设沿用已有手册输入：266 nm QE 0.639、读出噪声 12.1 e-、满阱 16000 e-；5 um 像元、10 bit。增益 16 e-/DN、黑电平 64 DN、暗电流 1000 e-/s 是现有示例假设，不是手册保证。其他波长应输入该波长的实测 QE。功率过高导致饱和时，应明确设置透过率或曝光，程序不会偷偷归一化或自动调暗。

`Camera DN` 每幅图标注原始 DN 的最小值、最大值和平均值。三个通道的 Sample/Reference 使用共同色标；`Observed DN (shared)` 使用这些图的共同最小/最大值，便于观察低亮度分布，`Full ADC range` 使用 ADC 全范围。两者只改变显示窗口，仍以真实 DN 标注，不做通道峰值归一化；差值使用共同的正负 DN 色标。移动鼠标或点击像素可在图上方读取该像素的坐标和原始 DN；Reference 平均和差值保留小数。

`Camera brightness` 表格给出逐幅图的最小值、最大值、均值、ROI DN 总和、处于 ADC 上限的像素数和 Sample 实际满阱像素数。平均 Reference 不包含逐帧的满阱掩膜，因此显示 N/A，而不是推算其满阱率。10 bit 上限是 1023 DN，与 uint16 存储容器无关。ADC MAX / FULL WELL 红字提示真实饱和；全白且所有像素为 1023 时，显示拉伸不能恢复丢失的亮度或形状，必须降低明确设置的曝光或透过率后重新模拟相机。旧相机目录也可直接显示这些统计；新报告另存 `brightness` 摘要，不改变原始 NPZ/TIFF。

相机结果保存到独立目录：`camera_data.npz`、`camera_settings.json`、`camera_report.json`、`ch1/sample_raw.tiff`、`ch1/reference_first_raw.tiff` 等，并完整复制 `optical_scene/` 原始光学场。Load 支持相机目录；Export raw bundle 导出当前相机及其原始光学数据。TIFF 是单次原始曝光，不是显示图；参考平均可能是小数，保存在 NPZ。

NPZ 中的 `ch1_digital_dn`、`ch1_reference_mean_dn`、`ch1_difference_dn` 为显示帧；`ch1_gray_roi_series_dn` 列为 sample ROI 和、平均 reference ROI 和及两者之差。采集量不超过现有原始栈上限时，还保存 `ch1_gray_sample_raw_dn` 和 `ch1_gray_reference_raw_dn`；超过时仍保存每组 ROI DN 序列和首组原始图。CH2/CH3 同理。

示例相机设置见 [Gemi 面阵近似](../examples/grating_camera_gemi_area.json)，对应光学中继设置见 [40 倍、透过率 1e-6 的明确测试条件](../examples/grating_contact_camera_optics.json)。这只是可复现的测试条件，不是推荐实机曝光：

```powershell
& build/rcwa_env/Scripts/python.exe -m grating_scene.camera --source output/my_grating_scene --output output/my_grating_camera --settings examples/grating_camera_gemi_area.json --timeout 600
```

用 `grating_scene.camera.load_camera(path)` 可同时验证原始文件哈希并读取报告与数组。相机报告的 `optical_model_validated=false` 是有意保留的边界。

## 原始数据读取

使用 `numpy.load(path, allow_pickle=False)` 读取所有 NPZ。复数以 NumPy 原生 complex128 保存，无需手工拼接实虚部。

| 文件 | 内容 |
|---|---|
| `config.json` | 本次完整输入及物理单位 |
| `reference/input.lua`, `sample/input.lua` | 实际执行的 S4 输入 |
| `reference/stdout.txt`, `sample/stdout.txt` | 原生带标签输出，包括复模态和近场 |
| `reference/modes.npz`, `sample/modes.npz` | 原始复反射级次、原生 E/H、功率、坐标 |
| `near_fields.npz` | XY/XZ/YZ 的坐标和 reference/sample/delta 复 E/H |
| `far_field.npz` | 级次、方向、传播标记、三路掩膜、并集、复反射场及效率分解 |
| `images.npz` | 像面坐标、三路 Jones 场、强度、真实差值及干涉项 |
| `scene_report.json` | 模型范围、检查状态、逐路功率和原始数据 SHA256 清单 |

`near_fields.npz` 使用 `xy_points_um`、`xy_reference_electric` 等键，数组为 `(N,N,3)`；XZ/YZ 同理。近场坐标单位 um，E 为相对单位入射 E0 的复幅。

`far_field.npz` 的 `orders` 为 `(modes,2)`；`directions` 为 `(modes,3)`。**先用 `propagating` 筛选**，倏逝级的方向为 NaN。`efficiency_delta` 是差场效率，`efficiency_difference` 是真实功率差，二者不同。

`images.npz` 使用 `ch1_electric_sample`、`ch1_sample_w_um2`、`ch1_reference_w_um2`、`ch1_difference_w_um2`、`ch1_delta_w_um2`、`ch1_interference_w_um2`；CH2/CH3 同理。Jones 场为 `(Ny,Nx,2)`，强度为 `(Ny,Nx)`，单位 W/um^2。像面 x/y 坐标已包含倍率。

```python
import numpy as np
from grating_scene.pipeline import load_scene

config, report, arrays = load_scene('output/my_grating_scene')  # 同时验证清单哈希
sample = arrays['images']['ch3_sample_w_um2']
reference = arrays['images']['ch3_reference_w_um2']
difference = arrays['images']['ch3_difference_w_um2']
assert np.allclose(difference, sample-reference)
```

## 验证范围

软件检查包括：原生 S4 幅相/功率、零对比度、FFT 与直接 Fourier 求和、Parseval、倍率、分析器、相干交叉项、解析像元积分与独立积分、原始 DN 的 SNR 重算、原始文件完整性、Qt 异步调用与取消。它们不能替代空间截断、球切片、超胞间距、有限边界、材料标定和真实相机验证。

可运行固定超胞原生收敛对照：

```powershell
& build/rcwa_env/Scripts/python.exe -m tools.grating_scene_validate --baseline output/my_grating_scene --output output/my_grating_controls --timeout 1800
```

该工具分别增加 X 截断、Y 截断、球切片，再检查联合加密，逐项保存原始场和资源记录；默认省去这些对照中重复的近场提取，`--with-near-fields` 可一并计算。被内存预检拦下的算例明确记为 blocked，而不是视作通过。

### 内存预检的可选标定

默认采用通用保守预估。开发者可用 `tools.grating_memory_calibration` 对同一个 S4 二进制已完成任务的实测进程树峰值进行核对，再给 `grating_scene.run` 指定 `--memory-profile`。标定必须匹配二进制哈希、保留至少 15% 的已测余量，并限制矩阵尺寸和层数的外推；不是把内存守护阈值调低。

UI 的 Memory profile (optional) 可填入该 JSON，应用于原生计算和截断/切片对照。留空仍使用通用预估。镜像一致性检查采用 2% 的项目门槛；未通过时局部步也不标记为通过，不对原始复场做对称化。

无论是否使用标定预估，3 GiB 系统余量、进程树上限和超时都继续生效。迁移环境或更换 S4 二进制后应重新标定，不能把旧预估当成硬内存保证。

新增空气区近远场检查使用所有反射级（包含倏逝级）和入射波重建 E/H，与高于整个结构的原生近场逐点比较。它检查相位、坐标和两类导出的衔接；即使误差接近舍入精度，也不认证有限基组或孤立缺陷。
