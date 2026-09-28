# Far-field 远场仿真程序用户使用指南

## 全波 FEM 实验页

运行本目录 `main.py` 后，顶部新增 `Full-wave FEM (experimental)`，
可查看近场、远场、三路独立像面及相干干涉分解。原有功能仍在 `SCATMECH` 页。
该页使用单独的 NGSolve 环境，默认三阶 `balanced` 与 128 点光瞳，
可切换二阶 `preview`；带取消、内存/时间保护与结果载入。
新增 `PML order (0=same)`、`PML full-order fraction` 和 `Local mesh scale`。
默认 PML order=0，保持原来的同阶吸收层；局部尺寸倍率默认 1。
选择 PML order=2 时，full-order fraction=0.5 让吸收层内半部保持物理区阶数，
只降低外半部单元的阶数；局部倍率 0.9 将颗粒及粗糙区目标单元尺寸减少 10%。
这是可选的内存/精度折中，不是默认高精度方案。整个 PML 均降为二阶的试验
在独立 Mie 校验中产生约 2.88% 消光误差，不推荐用于高精度结果。
保留内半部后该基准误差约 0.77%，仍大于原全三阶的约 0.18%。
详见 [低内存三阶加密报告](../fullwave_prototype/LOW_MEMORY_REFINEMENT_CN.md)。
详细参数、已完成基准与限制见
[后续验证与实验界面说明](../fullwave_prototype/FOLLOWUP_VALIDATION_CN.md)。
最新公式与数值结果见 [理论模型](../fullwave_prototype/THEORETICAL_MODEL_CN.md)、
[理论交叉校验结果](../fullwave_prototype/THEORY_VALIDATION_CN.md)。
2026-09-28 分文件发行版已包含 FEM 界面和独立求解器，无需安装 Anaconda。
完整解压后启动 `FarField_Simulation.exe`，保留同级 `fem` 和 `_internal` 目录。
发行版 FEM 结果默认位于 `%LOCALAPPDATA%\FarField_Simulation\results`。

## 三路独立 20× 成像

在左侧 `Independent image planes` 中设置成像参数，点击
`Calculate three image planes`，在 `Three image planes` 页查看三路独立像面。
默认物镜焦距为 20 mm、倍率为 20×，管镜焦距自动计算为 400 mm。
三路沿用各自的 NA 和倾角；NA=0 时该路关闭。

`Incident irradiance` 是垂直于入射光束截面的照度，单位 mW/μm²。
`Image defocus` 是探测器相对于理想像面的轴向离焦，单位 μm。
`Power transmission` 是功率透过率；`Pupil analyzer` 作用于各路准直光瞳的局部偏振。
以上成像设置共同应用于三路，但三路复场和图像分别计算，不进行相干合束。

图像显示绝对照度 W/μm²，各路使用独立色标，颜色采用幂律拉伸以显示弱信号，
截线为线性强度。总功率单位为 W，不做峰值归一化。
`Pupil samples / axis` 控制光瞳采样，`Multipole order` 控制球谐展开阶数。
应逐步增加两者检查收敛；图中很小的 `Flux error` 仅表示离散传播的功率守恒，
不代表物理模型已经具有同等精度。

该模式目前计算平滑多层晶圆上的单颗粒散射像，采用理想无像差物镜和管镜。
不包含完整近场、粗糙度散斑和镜面泄漏，因此运行时需关闭粗糙度。
当收集光锥包含镜面反射方向时会提示需要暗场挡光，但不会自动模拟实际挡光器。
颗粒和膜层材料仍使用界面中的 n、k；更改波长不会自动更新材料色散。

导出包括三路图像 PNG、保存复 Jones 光瞳及像面 E/H 模态的 NPZ 和完整参数 JSON。
详细公式、场单位和模型边界见 [IMAGING_MODEL.md](IMAGING_MODEL.md)。

### 三路探测器

勾选 `Detectors (three independent ROIs)` 后，点击同一个成像计算按钮，
在 `Detectors` 页查看三路的预期颗粒电子数和带噪声读出图。相机参数共同应用于三路，
但每路按自己的像场独立积分，噪声也独立生成，不相干叠加三路光场。

- `Pixel pitch`：像素边长；`ROI columns/rows`：以像面原点为中心的矩形像素区域。
- `Exposure`：曝光时间；`Quantum efficiency`：量子效率，需填写激光波长处的实际值。
- `Read noise`：读出噪声 RMS，单位电子；`Dark current / pixel`：单像素暗电流。
- `Background / pixel`：额外均匀背景电子率，不代表已计算了晶圆粗糙度或镜面泄漏。
- `Full well / pixel`：满阱电子数；`Noise seed`：固定后可复现同一参数下的随机读出。

程序先按像素实际面积积分照度，再根据光子能量和量子效率换算电子数；
对信号、暗电流和背景生成泊松噪声，按满阱截断后再加入高斯读出噪声。
读出可能出现负值或略高于满阱值，这是读出噪声的结果，不是负光子数。
默认相机参数仅为示例，不对应具体型号。

图中 `ROI` 表示该矩形像素区域接收到的功率占本路计算像面总功率的比例，
不是相对于半球的收集率。`Ideal SNR` 是背景均值已知时、饱和前的理论信噪比。
出现 `saturated` 提示时应缩短曝光或降低照度，不能直接把此 SNR 当作实际读出性能。
NA=0 时颗粒信号为零，但开启的探测器仍有暗电流和读出噪声。

ROI 超出计算像场会报错，请增加 `Pupil samples / axis` 或减小 ROI/像素尺寸。
导出新增探测器图、各像素光子数/电子数、随机读出、饱和掩码和参数。

## 1. 软件介绍

Far-field Simulation 是一个面向晶圆表面颗粒散射分析的独立桌面程序。
程序使用 pySCATMECH 建立球形颗粒位于晶圆表面的散射模型，计算上半空间的
远场强度和偏振分布，并分析三路理想物镜能够收集到的远场能量。

程序主要提供以下功能：

- 设置晶圆和球形颗粒的复折射率、厚度及颗粒尺寸。
- 设置从晶圆到表面的多层镀膜结构。
- 使用 SCATMECH ABC PSD 模型加入表面或多层界面粗糙度。
- 使用矩形照明区域的长和宽确定粗糙度背景的有效面积。
- 设置激光波长、入射角和偏振态。
- 显示远场强度、偏振方向及入射几何关系。
- 分别设置三路物镜的 NA 和倾斜角度。
- 显示三路物镜的三维收集光锥和二维出瞳能量分布。
- 分别计算每一路及三路去重后的几何体积占比和能量收集占比。
- 批量扫描一系列颗粒尺寸，并绘制三路积分收集强度和能量占比曲线。
- 导出图像、数值数据和仿真参数。

## 2. 启动程序

打开发布目录：

```text
FarField_Simulation\
|-- FarField_Simulation.exe
`-- _internal\
```

双击 `FarField_Simulation.exe` 启动程序。

注意：本程序采用目录式打包。`FarField_Simulation.exe` 必须和 `_internal`
文件夹保持在同一目录中。移动或发送程序时，请复制整个
`FarField_Simulation` 文件夹，不要只复制 EXE。

## 3. 界面组成

窗口左侧是参数和操作区域，右侧是绘图区域。绘图区域包含三个标签页：

- `3-D far-field / pupil`：三维远场或三维物镜收集光锥。
- `2-D far-field / pupil`：二维远场或三路物镜出瞳能量分布。
- `Particle-size sweep`：颗粒尺寸扫描的三路强度与能量占比曲线。

左侧底部显示仿真进度和运行状态。

## 4. 晶圆参数

### Wafer material

选择晶圆材料。程序提供 Silicon、Fused silica、Sapphire、Silicon carbide、
Gallium arsenide 和 Custom 等预设值。

材料预设值用于快速填写折射率。选择材料后，仍可手动修改 `n` 和 `k`。

### Refractive index n

晶圆复折射率的实部 `n`。

### Extinction coefficient k

晶圆复折射率的消光系数 `k`。无吸收材料可设置为 0。

### Thickness

晶圆厚度，单位为 `um`，可设置范围为 0.001 至 5000 um。

对于光学厚度很大的吸收晶圆，程序会自动使用半无限基底模型；对于满足条件的
薄层，会使用有限厚度相干膜层模型。实际采用的方式会显示在运行状态中。

### Wafer coating stack

膜层表格用于设置晶圆表面的多层镀膜。表格顺序固定为：

```text
第 1 行：最靠近晶圆基底
最后一行：最靠近空气的顶层
```

每一层可以设置材料名称、复折射率 `n + ik` 和厚度，厚度单位为 `nm`。
`Add`、`Remove`、`Up`、`Down` 分别用于增加、删除和调整膜层顺序。

所有膜层均通过 SCATMECH `FilmStack/Stack_StackModel` 传入
`Bobbert_Vlieger_BRDF_Model`，不是使用程序自定义的薄膜近似。

### Surface and interface roughness

选中 `Include roughness background` 后，程序会把颗粒散射与粗糙度散射合并。
界面模式包括：

- `Top exposed interface`：只计算最外侧暴露表面的粗糙度。
- `All interfaces correlated`：所有膜层界面具有相同且完全相关的粗糙度。
- `All interfaces uncorrelated`：所有膜层界面的粗糙度互不相关。

程序使用 SCATMECH Rayleigh-Rice 光滑表面微扰模型。粗糙度采用二维 ABC PSD：

```text
S2(f) = A / (1 + B^2 f^2)^(C/2)
A = RMS^2 B^2 (C - 2) / (2 pi)
```

其中：

- `RMS height`：积分均方根高度，单位为 `nm`。
- `ABC length B`：PSD 转折长度参数，单位为 `um`。
- `ABC exponent C`：高频衰减指数，必须大于 2。

无膜层时使用 `Microroughness_BRDF_Model`；有膜层时根据界面模式使用
`Roughness_Stack_BRDF_Model`、`Correlated_Roughness_Stack_BRDF_Model`
或 `Uncorrelated_Roughness_Stack_BRDF_Model`。

### Rectangular illumination footprint

`Footprint length` 和 `Footprint width` 是晶圆表面上的矩形照明区域尺寸，
单位为 `um`。程序自动计算：

```text
Effective area = Footprint length x Footprint width
```

当前组合模型假设该矩形照明区域内有一个指定 PSL 颗粒，因此 Bobbert-Vlieger
颗粒面密度设置为 `1 / Effective area`。SCATMECH
`Two_Source_BRDF_Model` 以 `correlation=0` 将颗粒与粗糙度 Mueller BRDF
按非相干强度相加，再按同一面积和入射/散射余弦因子换算回有效 DSC。

矩形长或宽增大时，颗粒本身的单颗粒 DSC 不变，但光斑覆盖的粗糙表面积增加，
因此粗糙度背景会相应增加。关闭粗糙度后，长和宽不参与计算，结果保持原来的
单颗粒 Bobbert-Vlieger DSC。

## 5. 颗粒参数

### Particle material

选择颗粒材料。程序提供 PSL、Silica、Silicon、Alumina、Aluminum、Gold 和
Custom 等预设值。

### Refractive index n / Extinction coefficient k

颗粒复折射率的实部和消光系数。材料切换后可继续手动修改。

### Sphere diameter

球形颗粒直径，单位为 `nm`，可设置范围为 1 至 10000 nm。

修改颗粒尺寸后，需要重新执行 `Run far-field simulation`。只有新的远场仿真
完成后，物镜收集计算才会使用更新后的颗粒远场结果。

## 6. 激光与入射参数

### Wavelength

真空波长，单位为 `nm`，可设置范围为 100 至 2000 nm。

### Incidence angle

入射光与晶圆表面法线之间的夹角，单位为度，可设置范围为 0 至 89.9 度。

### In-plane projection

软件采用以下坐标约定：晶圆表面为 `KX-KY` 平面，表面法线为 `+KZ`，
入射面为 `KY-KZ` 平面，光源位于负 `KY` 方向：

- 光源位于 `-KY` 一侧。
- 斜入射传播矢量为 `(0, +sin(theta_i), -cos(theta_i))`。
- 入射光在晶圆平面上的传播投影指向 `+KY`。
- 法向入射时面内投影为零，但仍沿用同一坐标基准。

pySCATMECH 使用相对于入射面的散射方位角，程序把其中的 `phi=0` 映射到
全局 `+KY`。计算结果的位置坐标和偏振基矢均绕 `+KZ` 旋转 90 度，物镜采样时
执行对应的逆变换，因此远场图、偏振方向和收集能量使用同一全局坐标。

### Polarization

- `P polarized`：P 偏振。
- `S polarized`：S 偏振。
- `Unpolarized`：非偏振光。

## 7. 数值与显示参数

### Angular step

远场角度采样步长，可设置范围为 2 至 15 度。

- 较小步长：角度分辨率更高，计算时间更长。
- 较大步长：计算速度更快，适合快速预览。

计算量大致随采样步长减小而快速增加。建议先使用 8 至 12 度检查参数，再使用
2 至 4 度生成最终结果。

### Display range

图像显示动态范围，可设置为 10 至 80 dB。该参数只改变颜色显示范围，不改变
原始远场能量和物镜收集占比。

## 8. 运行远场仿真

1. 设置晶圆、颗粒、激光和数值参数。
2. 单击 `Run far-field simulation`。
3. 等待进度条达到 100%。
4. 在右侧两个标签页查看三维和二维结果。

仿真在后台线程运行，计算过程中不要关闭程序。完成后，状态栏会显示 DSC
数值范围和晶圆模型说明。

### 三维远场图

三维表面表示上半空间方向半球：

- 坐标为归一化方向分量 `kx/k`、`ky/k` 和 `kz/k`。
- 表面颜色表示相对 DSC，单位为 dB。
- 短线表示局部线偏振方向，长度与线偏振度相关。
- 图中同时显示晶圆平面、表面法线和激光入射方向。
- 光源标在 `-KY` 一侧，橙色面内箭头指向 `+KY`。

### 二维远场图

二维图是方向半球在 `kx/k`、`ky/k` 平面上的投影：

- 圆周对应接近 90 度的掠射方向。
- 中心对应晶圆表面法线方向。
- 颜色表示用于观察分布形状的归一化散射强度。
- 短线表示投影后的偏振方向。
- 图左侧箭头表示入射传播方向的面内投影 `+KY`。

二维远场颜色用于显示分布，不应直接当作绝对功率读数。物镜能量占比使用未被
二维显示归一化修改的原始远场强度计算。

## 9. 三路物镜收集设置

物镜被视为焦距为 `f`、无像差的理想单透镜。每一路均可单独设置 NA 和倾角。

### Focal length

三路物镜共用焦距，单位为 `mm`，默认值为 20 mm。

### Channel 1/2/3 NA

每一路的数值孔径可分别设置为 0 至 0.99。

- 光锥半角为 `asin(NA)`。
- 理想出瞳半径为 `f x NA`。
- `NA=0` 表示关闭该路。关闭通道的能量为 0，不显示光锥，也不参与合并计算。

### Channel 1/2/3 tilt

每一路物镜光轴相对晶圆表面法线的倾角，可设置范围为 -89.9 至 89.9 度。

- `0 deg`：物镜光轴沿晶圆表面法线。
- 正角度：光轴朝 `+kx` 方向倾斜。
- 负角度：光轴朝 `-kx` 方向倾斜。

默认三路角度为 `+45 deg`、`-45 deg` 和 `0 deg`。

## 10. 运行物镜收集分析

推荐操作顺序：

1. 先设置颗粒等参数并完成一次 `Run far-field simulation`。
2. 设置三路物镜的 NA、焦距和倾角。
3. 单击 `Run ideal objective example`。
4. 查看三维光锥、二维出瞳、几何体积占比和能量占比。

如果尚未运行远场仿真，物镜分析会使用内置的 synthetic dark-field 示例分布，
状态栏会明确显示这一点。如果已经完成远场仿真，则使用最近一次 SCATMECH
结果。修改颗粒尺寸、波长、材料或入射条件后，应先重新运行远场仿真，再运行
物镜收集分析。

### 三维物镜收集图

- 半球颜色显示远场能量分布。
- CH1、CH2、CH3 光锥表示三路物镜实际覆盖的方向范围。
- 高亮区域表示至少被一路物镜收集的远场区域。
- 标题分别显示三路去重后的合并几何体积占比和能量占比。

### 二维三路出瞳图

- 三个子图分别对应 CH1、CH2 和 CH3。
- 坐标单位为 mm，显示物镜出瞳上的能量分布。
- 三个通道使用同一个固定 dB 色标，不对每一路单独归一化，因此可以比较强弱。
- 标题分别显示该路 NA、倾角、几何体积占比、能量占比和平均 DoLP。
- 白色短线表示出瞳上的线偏振方向。
- `NA=0` 的通道显示为 `Disabled`。

## 11. 几何体积与能量占比

程序将 `Volume` 和 `Energy` 作为两套独立指标计算和显示。

### 几何体积占比 Volume

几何体积计算完全不使用远场强度，只由每路物镜的 NA、倾角以及光锥之间的
重叠关系决定。程序在单位半球上计算光锥覆盖的立体角。由于径向长度相同，
单位半球中的光锥体积比等于立体角占比：

```text
几何体积占比 = 光锥覆盖区域内的积分(sin(theta))
               / 整个上半球的积分(sin(theta))
```

因此，改变颗粒尺寸、材料、激光波长、偏振态或远场强弱不会改变 `Volume`。
改变 NA 或倾角会改变该指标。

### 能量占比 Energy

能量占比根据当前远场能量分布积分计算：

对上半空间方向，程序按以下形式计算：

```text
能量占比 = 光锥覆盖区域内的积分(I(theta, phi) sin(theta))
           / 整个上半球的积分(I(theta, phi) sin(theta))
```

其中 `I(theta, phi)` 是远场强度，`sin(theta)` 是立体角积分权重。因此：

- 相同大小的光锥位于不同方向时，几何体积可以相同，但能量占比可能不同。
- NA 增大时覆盖范围增大，但新增区域的实际贡献取决于该处远场强度。
- 两套指标的三路合并结果都使用三路光锥掩膜的并集。
- 几何体积和能量计算中的多路重叠区域都只计算一次，不会重复累加。

每一路标题显示该路单独的 `Volume` 和 `Energy`；图像总标题显示三路去重后的
`union volume` 和 `union energy`。

## 12. 颗粒尺寸扫描

`Particle-size sweep` 参数组用于连续计算一系列球形颗粒尺寸：

- `Start diameter`：扫描起始直径，单位为 nm。
- `End diameter`：扫描终止直径，单位为 nm。
- `Diameter step`：相邻计算点的直径间隔，单位为 nm。

扫描至少需要两个点，最多允许 200 个点。如果终止直径不能被步长精确到达，
程序计算不超过终止值的最后一个点，并在完成状态中显示实际范围。

### 运行步骤

1. 设置晶圆、颗粒材料、波长、入射角、偏振和角度步长。
2. 设置三路物镜的 NA 和倾角。`NA=0` 的通道在扫描中保持关闭。
3. 设置起始直径、终止直径和步长。
4. 单击 `Run particle-size sweep`。
5. 等待进度达到 100%，程序会自动切换到 `Particle-size sweep` 标签页。

扫描不使用单个 `Sphere diameter` 的值，而是为扫描序列中的每个尺寸重新建立
SCATMECH 模型。每个点都使用当前界面的其他参数，不需要先单独运行远场仿真。

如果启用了粗糙度，扫描曲线表示“每个粒径的颗粒散射 + 同一矩形照明区域的
粗糙度背景”。粗糙度背景不会随颗粒直径变化，因此小颗粒区域可能出现由表面
背景决定的平台。关闭粗糙度可得到纯单颗粒尺寸响应。

### 结果曲线

上图显示三路物镜分别收集的积分 DSC：

```text
Collected DSC = 光锥内积分(DSC(theta, phi) sin(theta))
```

单位为 `um^2`。纵轴使用对数刻度，适合比较跨越多个数量级的颗粒散射强度。

下图显示每一路积分 DSC 相对于整个上半球积分 DSC 的能量占比，单位为 `%`。
三条曲线分别对应 CH1、CH2 和 CH3，图例同时显示该路 NA 和倾角。

扫描点数和较小的 `Angular step` 都会增加计算时间。建议先用较大角度步长进行
少量点预览，再用目标角度步长运行最终扫描。

## 13. 导出结果

完成远场仿真或物镜收集分析后，单击 `Export figures and data`，选择输出目录。

程序会生成：

| 文件 | 内容 |
| --- | --- |
| `far_field_3d_intensity_with_polarization.png` | 当前三维结果图 |
| `far_field_2d_intensity_with_polarization.png` | 当前二维结果图 |
| `far_field_data.npz` | NumPy 压缩数组数据 |
| `simulation_parameters.json` | 本次仿真参数和入射坐标方向约定 |

`NPZ` 文件可使用 NumPy 读取：

```python
import numpy as np

data = np.load("far_field_data.npz")
print(data.files)
intensity = data["intensity"]
```

建议每次导出到独立目录，以免覆盖之前同名的结果文件。

粒径扫描结果会单独导出以下文件：

| 文件 | 内容 |
| --- | --- |
| `particle_size_sweep.png` | 三路扫描曲线图 |
| `particle_size_sweep.csv` | 每个粒径点的三路积分 DSC 和能量占比 |
| `particle_size_sweep_data.npz` | NumPy 数组数据 |
| `particle_size_sweep_parameters.json` | 仿真、物镜、扫描参数和入射坐标方向约定 |

## 14. 常见问题

### 单击物镜分析后，颗粒尺寸变化没有反映出来

修改颗粒尺寸后，先单击 `Run far-field simulation` 生成新的远场结果，再单击
`Run ideal objective example`。物镜分析使用最近一次成功完成的远场结果。

### 某一路不需要使用

将该路 NA 设置为 0。不要只依靠倾角把光锥移出显示范围。

### 计算速度较慢

增大 `Angular step`。快速检查可使用 10 至 15 度，最终计算再减小到 2 至 4 度。
粒径扫描还可通过增大 `Diameter step` 或缩小扫描范围来减少计算点数。

### EXE 无法启动或提示缺少 DLL

确认 EXE 旁边存在完整的 `_internal` 文件夹。不要从发布目录中单独取出 EXE。
若文件由压缩包传输，请先完整解压后再运行。

### 远场仿真按钮不可用

这通常表示当前运行环境缺少 pySCATMECH 或 SCATPY。使用发布目录中的完整打包
版本；从源码运行时，应使用项目配置的 Anaconda `NIST2` 环境。

### 显示范围改变后，能量占比是否会改变

不会。`Display range` 仅改变图像色标范围，积分使用原始远场强度。

### 颗粒参数改变后，几何体积为何没有变化

这是正常结果。几何体积只由 NA、倾角和光锥重叠决定；颗粒参数只会影响远场
能量分布和 `Energy` 指标。

### 粒径扫描曲线中某一路为零

检查该路 NA 是否设置为 0。`NA=0` 表示关闭该通道，积分 DSC 和能量占比均为 0。

## 15. 模型适用范围

- 当前原生模型针对位于平面晶圆表面的球形颗粒。
- 颗粒、膜系、粗糙度及其 Mueller BRDF 合成均调用 SCATMECH 模型。
- 多层镀膜按相干 `FilmStack` 计算，材料预设仅是参考值，实际计算使用表格中的
  `n`、`k` 和厚度。
- 粗糙度采用 Rayleigh-Rice 光滑表面一阶微扰模型。RMS 较大或表面斜率很大时，
  应检查该模型的适用性；程序会在光滑表面参数明显偏大时给出运行状态警告。
- 粗糙度与颗粒默认按非相干强度相加，不包含二者之间的相位相关或多重散射。
- 矩形照明区域是晶圆表面 footprint，不是入射光束垂直截面。
- 物镜模型为无像差理想单透镜，不包含实际镜头像差、遮挡、透过率和镀膜损耗。
- 几何体积占比使用单位上半球内光锥体积作为基准，与绝对焦距和远场强度无关。
- 能量占比以仿真得到的上半空间远场为分母，不包含下半空间传播能量。
- 材料预设值是接近 355 nm 的参考值。波长变化较大时，应填写对应波长下的实际
  复折射率。
- 最终工程设计应结合真实物镜透过率、偏振响应、机械遮挡和实验标定结果。

SCATMECH 官方模型说明：

- [Bobbert_Vlieger_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)
- [StackModel](https://pages.nist.gov/SCATMECH/docs/stackmodel.htm)
- [Microroughness_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/microroughbrdf.htm)
- [Correlated_Roughness_Stack_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/correlated_stack.htm)
- [Two_Source_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/two_source.htm)
