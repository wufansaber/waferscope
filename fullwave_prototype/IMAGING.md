# 近场、远场、三路独立像面

最新公式与检查：[理论模型](THEORETICAL_MODEL_CN.md)、[理论交叉校验](THEORY_VALIDATION_CN.md)。
源码 UI 的 FEM 页现在默认使用三阶 `balanced`、128 点光瞳和 12 GiB 记账内存限制；
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
入射辐照度为 1e-6 W/um^2，定义在垂直于入射传播方向的面积上。

## 不重新解全波，调整成像参数

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.imaging_run --source build/fullwave_imaging/default/scattering --output build/fullwave_imaging/na04 --na 0.4 0.4 0.4 --tilts 45 -45 0
```

`--na` 和 `--tilts` 各给出三个数；NA=0 关闭对应通道。
其他选项为 `--magnification`、`--objective-focal-mm`、`--defocus-um`、
`--irradiance-w-um2`、`--transmission`、`--analyzer none/x/y/45/135`、
`--pupil-samples 32/64/96/128` 和 `--padding 2/4/8`。
偏振片方向是各自物镜局部出瞳的 X/Y，不是所有通道统一的全局 X/Y。

`--source` 模式只改变成像参数，使用源报告中的材料、波长、形状与全波结果。
首次求解可用 `--order`、`--mesh-scale`、`--rms-nm`、`--gap-nm`。
其他散射参数可调用 `rough_run.run(Parameters(...), Roughness(...), path)` 设置，
随后将这个目录传给 `--source`。旧版结果没有保存体源，会明确报错要求重新求解。

## 三个位置的量不相同

1. 近场：显示总电场的 |E/E_inc|^2，包含入射、反射和透射背景。
2. 远场：颜色表示 DSC，单位 um^2/sr；半球投影半径是 sin(theta)，
   不是距离。图上的三个边界表示三个物镜的收集范围。
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
目前不包含镜头像差、多层膜、有限照明光斑、相机像素和噪声。
这是完整的理想光学计算链路，不是已经标定验证的真实仪器数字模型。
源码 UI 已接入全波实验页；现有 EXE 尚未包含这一更新。

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
