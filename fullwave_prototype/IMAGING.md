# 近场、远场、三路独立像面

最新公式与检查：[理论模型](THEORETICAL_MODEL_CN.md)、[理论交叉校验](THEORY_VALIDATION_CN.md)。
源码 UI 的 FEM 页现在默认使用 Gaussian PSD 粗糙面、9x9 高度节点、
一阶 `preview`、128 点光瞳和 12 GiB 记账内存限制；
命令行省略 `--accuracy` 时仍为二阶 `preview`，以保留旧命令的含义。

重要：后续[物理校验](PHYSICAL_AUDIT_CN.md)发现原二阶配置的能量闭合与 PML 敏感性
尚未达到 1% 筛查门槛。出瞳到像面的离散功率一致性不能代替物理精度验证。

新增入口 `fullwave_prototype.imaging_run` 完成一条可运行的实验成像链路：

```text
颗粒 + 粗糙表面 Maxwell 全波求解
  -> 复数近场 E/H
  -> 半空间体源辐射积分，得到复数远场
  -> 三路各自的 NA / 倾角出瞳
  -> 理想无像差物镜 + 管镜矢量传播
  -> 三路独立像面上的法向能流密度
```

不是将远场强度直接做 FFT，也没有假设所有角度具有相同相位。
出瞳直接由保存的体极化源计算，不从 5 度远场显示图插值。
复用的矢量传播算法已移到根目录 `optical_propagation.py`；原 SCATMECH
成像模块仍通过原接口调用同一算法，不要求新环境安装 PyQt 或 SCATMECH。

## 一次运行看到三个位置

在仓库根目录运行：

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.imaging_run --output build/fullwave_imaging/default
```

生成文件：

- `near_far_image.png`：近场 YZ 截面、上半球远场与收集边界、三路像面。
- `image_decomposition.png`：每路有颗粒、无颗粒背景、场变化和干涉项。
- `image_fields.npz`：每路复数出瞳、复数像面 E/H、能流密度及坐标。
- `imaging_report.json`：成像参数、收集截面、像面功率、功率一致性和周期边缘指标。
- `scattering/`：本次全波近远场结果与可重复使用的体源积分数据。

默认示例保持上一阶段的 355 nm、100 nm 颗粒、2 nm RMS 粗糙面、
球底高于平均面 10 nm、75 度 p 偏振照明。默认三路 NA 均为 0.5，
倾角分别 +45、-45、0 度，绕全局 Y 轴倾斜，三路独立计算、不相干涉。
物镜焦距 20 mm，倍率 20，管镜焦距 400 mm；默认无偏振片、透过率 1。
默认入射总功率 1 W，晶圆表面有效照明长宽分别为 2 um、2 um。
这些默认值同时适用于 FEM 界面与不指定照明参数的命令行运行。
界面显示的光束法向辐照度按 `I = P / (L × W × cos θ_i)` 换算；
默认 75 度入射时为约 0.965926 W/um^2。这里的 2 × 2 um 是晶圆表面
**有效功率归一化面积**，当前 FEM 默认仍是局部均匀平面波，未把它建成
真实的方形截断光斑，也未将有限粗糙面扩展到整个 2 × 2 um。
远场 DSC 与功率无关；图中半球功率和像面功率/能量按该辐照度缩放。

## 不重新解全波，调整成像参数

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.imaging_run --source build/fullwave_imaging/default/scattering --output build/fullwave_imaging/na04 --na 0.4 0.4 0.4 --tilts 45 -45 0
```

`--na` 和 `--tilts` 各给出三个数；NA=0 关闭对应通道。
其他选项为 `--magnification`、`--objective-focal-mm`、`--defocus-um`、
`--illumination-power-w`、`--illumination-length-um`、`--illumination-width-um`、
`--irradiance-w-um2`、`--transmission`、`--analyzer none/x/y/45/135`、
`--pupil-samples 32/64/96/128` 和 `--padding 2/4/8`。
功率和两个光斑尺寸必须一起指定。旧流程仍可只传 `--irradiance-w-um2`
直接指定光束法向辐照度；若同时传入两种口径，数值必须一致。
偏振片方向是各自物镜局部出瞳的 X/Y，不是所有通道统一的全局 X/Y。

`--source` 模式只改变成像参数，使用源报告中的材料、波长、形状与全波结果。
首次求解可用 `--order`、`--mesh-scale`、`--rms-nm`、`--gap-nm`。
FEM 页的 `Rough surface model` 可切换旧的逐样本归一化表面与
`gaussian` 统计相关表面；后者还可设置相关长度、粗糙斑块宽度和高度节点数。
选择 `gaussian` 时界面自动改为一阶 FEM 预览以控制内存；高阶计算需手动选择并检查资源。
只有重新运行 `Calculate` 才会改变散射源；`Re-image source` 不会重新生成粗糙面。
其他散射参数可调用 `rough_run.run(Parameters(...), Roughness(...), path)` 设置，
随后将这个目录传给 `--source`。旧版结果没有保存体源，会明确报错要求重新求解。

## 三个位置的量不相同

1. 近场：显示总电场的 |E/E_inc|^2，包含入射、反射和透射背景。
2. 远场：颜色表示 DSC，单位 um^2/sr；半球投影半径是 sin(theta)，
   不是距离。图上的三个边界表示三个物镜的收集范围。
   星形标记指向样品看向入射光源的方向：来自 -KY，默认 75 度入射；
   该标记不是被画入上半球散射场的入射光强度。
3. 像面：显示沿该像面法线的能流密度，单位 W/um^2；使用复数 E/H
   计算 Re(E_x H_y* - E_y H_x*)，不简单将三个电场分量的模平方相加。

远场和像面仅包含相对于平面参考的散射光，不包含无限平面波镜面光束。
近场与像面的原始复电场不是同一种归一化：近场为 E/E_inc，成像内部为
通量归一化模式，不能把输出数组直接当作 SI V/m。

总览的三路像面共用绝对色标，使用平方根颜色刻度，不按各路峰值归一化。
分解图每一列在三路间共用色标，但背景和干涉列使用各自标出的绝对范围。
像面图仅放大显示中心区域；标题中的功率是完整 FFT 平面的积分，不是显示窗口积分。

## 光学假设

采用理想无像差、满足正弦条件的空气物镜，矢量旋转保留偏振，
通过 1/sqrt(cos(theta_local)) 将固体角能量转换为出瞳面积采样权重。
物像方向余弦满足 k_image,xy = -k_object,xy / M，对应倒像。
倍率固定时同时按比例改变物镜和管镜焦距，不改变这套理想模型的像面分布；
物理出瞳坐标则随焦距变化。

每路的物方参考焦平面经过 (0,0,0)，并垂直于该路光轴。
因此两个斜视通道和正视通道不是在观察同一个方向的二维切片。
颗粒高度和粗糙面形貌的相位已经包括在体源辐射中，不能再次人为加一遍离焦。
`--defocus-um` 是相对于理想共轭像面的额外像方轴向位移。

三路是独立观察结果。NA 锥有重叠时不能将三路功率相加当作实际总收集率；
本模块不对重叠角度分光，也不包含多物镜机械遮挡或分光器损耗。
`admits_specular` 会提示收集锥是否包含镜面方向，但模型没有自动补入镜面背景。

## 相干背景分解与验证

分别传播有颗粒和无颗粒的复数场；场变化 dE/dH 包含颗粒对粗糙面响应的改变。
像面保持如下分解：

```text
S_with = S_without + S_change + S_interference
S_interference = Re(dE x conj(H_without) + E_without x conj(dH))_z
```

相减两张强度图得到的是 S_with-S_without，而不是 S_change。
测试覆盖自由空间偶极子的解析收集截面、位移倒像、倍率缩放、离焦、
偏振片功率分解、NA=0 关闭通道及像面相干分解。
出瞳功率与像面积分的一致性是离散传播自洽检查，不是 FEM 物理精度保证。

已完成部分 FEM 边界、近场、能量平衡及独立成像基准，结果不代表全部工况已经收敛。
FFT 图像具有周期性；提高出瞳采样数扩大无混叠视场，增加 padding 仅使像面采样更密，
不能替代出瞳采样收敛。硬边 NA 有衍射长尾，有限显示区域不会包含全部能量。
可选有限角谱光束、Zernike 像差、相机像素和噪声已经接入；默认关闭以保留原来的理想模型。
FEM 仍不包含多层膜，也没有真实物镜或相机的实验标定。
这是完整的理想光学计算链路，不是已经标定验证的真实仪器数字模型。
源码 UI 和 2026-09-28 发布包已接入全波实验页；新增功能是否包含在 EXE 中以对应发布记录为准。

## 数值配置更新

新计算保留默认 `--accuracy preview` 二阶设置。显式 `--accuracy refined` 使用
三阶有限元、静态凝聚、原线性径向 PML、七阶体源积分，需要明显更多内存。
保存的体源仍使用其原始
计算参数，不会被自动升级。高阶配置不代表全部工况具有 1% 精度保证，
实验依据与限制见 [SOLVER_DECISION_CN.md](SOLVER_DECISION_CN.md)。

新增 `--accuracy balanced` 保留三阶与静态凝聚，PML 强度设为 1、PML 目标网格比例为
1.25，并使用 sparsecholesky。它已完成本轮三阶粗糙面配对成像，仍需按目标工况做收敛检查。
`--pml-strength`、`--pml-mesh-scale` 和 `--integration-bonus` 可显式覆盖数值控制；
在 `--source` 模式下这些选项会被拒绝，以免旧源被误标为新数值配置。

## 独立积分与线程控制

新增 `--threads` 和源码 UI 的 `FEM threads`，显式设置 NGSolve 工作线程。
资源保护保持 BLAS/OpenMP 单线程，内存与超时限制仍然生效。
复用 `--source` 时无需进行 FEM 求解，因此拒绝附加新的 `--threads` 参数。
同一 FEM 出瞳的直接矢量积分、独立 BV 极坐标成像及本机线程实测见
[IMAGING_THREAD_AUDIT_CN.md](IMAGING_THREAD_AUDIT_CN.md)。

## 新增光束、像差、相机与缓存

详见 [MODEL_UPGRADE_CN.md](MODEL_UPGRADE_CN.md)。CLI 默认输出 compact 存档，省略像面
E/H 四个大数组，但保留复出瞳、全部强度和相干分解。需要像面复场时选择 `--archive-mode full`。
API `calculate_images` 默认仍保留复场，兼容原有分析脚本。

`Re-image source` 选择含 `report.json` 与 `paired_fields.npz` 的 scattering 目录，
或选择其上一级完整结果目录。此操作只应用成像和相机参数，物理参数继续使用原始散射源。
改变颗粒、材料、粗糙度、照明方向或光束形状请使用 `Calculate`，不能用重成像代替重新求解。
仅改变总功率或有效照明面积可用 `Re-image source` 重新缩放功率和像面结果；
它不会改变归一化近场或远场 DSC，也不会生成实际的有限光斑场。
界面自动使用按源内容校验的出瞳缓存，CLI 通过 `--pupil-cache DIRECTORY` 启用。

## 三个近场截面

新计算直接从同一次 FEM 解采样三个 101×101 截面：XZ（Y=0）、YZ（X=0），
以及 XY（Z=颗粒中心高度，即粒径/2+间隙）。XY 并非晶圆表面投影；截面位置在图标题中标注。
横纵坐标保持等比例；三个截面共用 `|E/E_inc|²` 色标，不按各自峰值归一化。
平方根仅用于色彩映射。这里是包含照明背景的总电场强度指标，并非近场能流密度。
青色实线标出球截面和真实三角粗糙面的交线，竖直截面中的虚线为 Z=0 参考面。

界面 `Near field` 内可选 `Overview`、`XZ`、`XY`、`YZ`。
输出总览 `near_field.png` 和三张独立图 `near_field_xz.png`、`near_field_xy.png`、
`near_field_yz.png`。`scattering/paired_fields.npz` 中保存每个截面的复电场、复磁场
（有颗粒、无颗粒和两者差场）、坐标与位置：`near_<xz|xy|yz>_e_<with|without|change>`、
`near_<plane>_hz0_<state>`、`near_<plane>_u_um/v_um/position_um`。
磁场单位为 `H*Z0/E_inc`，数组维度为 `[纵坐标, 横坐标, XYZ分量]`。
旧的 YZ 键仍然保留，报告中的 `near_field_planes` 记录各平面信息。

XZ/YZ 沿用原来的显示范围；XY 方形窗口必要时缩小，保证采样点位于物理域内而非 PML。
增加截面不增加 FEM 求解次数，仅增加采样和导出工作。
旧结果只含 YZ，仍然可以读取或重成像，但无法据此生成缺失的 XZ/XY；请重新计算以获得三个截面。

## 独立颗粒与粗糙面对照

在 FEM 页完成计算或载入含可复用 `scattering/paired_fields.npz` 的完整结果后，
可分别设置 Ch1/Ch2/Ch3 的 NA 和倾角，再点击 `Compare particle / haze`。
任一路 NA=0 表示关闭；例如仅开 Ch2，或开 Ch1+Ch3，三路 NA 都可独立取值。
载入结果时先恢复原结果的成像设置，之后可直接修改任意一路。
它按载入结果保存的物理参数再求解一次平整晶圆上的颗粒，
并从原结果提取同一粗糙面、无颗粒的 haze 场；不是把“有颗粒减无颗粒”的差场误称为孤立颗粒。
平整面源保存在结果目录的 `particle_only_scattering/`，后续只更改 NA、倾角等成像参数时
直接复用两份 FEM 散射源，重新计算出瞳和像面，不再重新求解 FEM。

完成后 `Particle / haze` 页提供 `Far field`、`Far log`、`Image plane`、
`Image log` 和 `Totals`。颗粒与 haze 图采用各自带单位的数值色标，使弱颗粒信号仍可看清；
同一分量的三路像面共用色标。不同分量不可只凭颜色比较强弱，需看色条数值和总功率。
对数页改变颜色映射，不改动原始 DSC、辐照度或积分值。
远场总量为上半球 DSC 积分（um²）及按入射辐照度换算的功率；像面给出三路各自
完整像面的积分功率（pW）和按相机曝光时间换算的能量（fJ）。
三路物镜可能收集相同方向，不能把三路功率简单相加作为唯一总收集功率。

这两个输出属于**独立实验条件**：颗粒位于平整晶圆上、粗糙晶圆没有颗粒。
它们的功率或强度不能简单相加还原“颗粒在粗糙晶圆上”的相干结果。
原有 `Coherent decomposition` 页继续显示实际粗糙面上颗粒引起的差场与干涉。
UI 只在比较成功后启用新页；旧结果可正常打开。命令行入口为
`python -m fullwave_prototype.component_compare --result RESULT_DIRECTORY`，
通过界面启动的计算受内存和时间保护；直接运行上述命令不使用界面的保护设置。
EXE 需要重新打包才能包含新增代码。
