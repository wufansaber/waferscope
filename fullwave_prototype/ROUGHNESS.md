# 粗糙表面与颗粒的相干全波原型

整光斑独立 haze 已新增 [Whole-beam haze (SPM)](WHOLE_BEAM_HAZE_CN.md)。
它自动按照明功率覆盖与 PSD 带宽选择大面积随机表面，输出相干远场、
三路独立像面及统计区间，不把局部 FEM 面积放大。SPM 是一阶近似，
原局部 FEM 与颗粒耦合求解的范围、物理含义保持不变。

有限完全相干 Gaussian 光斑已接入实际背景场、配对求解、随机表面平均和
收敛计划；详情见 [有限相干照明](COHERENT_ILLUMINATION_CN.md)。该模式不按
有效面积外推局部 haze，不与无限平面波 NIST PSD 结果直接比较。界面可切回
`plane` 使用原有 NIST 对照流程。粗糙斑块之外仍是平整晶圆。

这是独立的实验 FEM 模块，已接入程序的 FEM 页；没有修改 NIST2，也没有重新打包 EXE。
它是有限元数值解，不是任意粗糙几何的闭式解析解。

## 本次实现

- 用真实三角形高度面构造局部凸起和凹陷，区域外仍为平面半无限基底。
- 颗粒与表面在同一个 Maxwell 方程中求解，包含相互作用，不是两个强度图相加。
- 同一网格分别计算有颗粒、无颗粒，位置、入射相位和粗糙面实现保持相同。
- 输出复数近场 E、H，半球远场复振幅及散射截面，并分离背景、场变化和干涉。
- 颗粒完全位于空气中。本版本要求球底高于整块表面的最高点；尚不处理粗糙面贴合、穿透或嵌入。

## 运行

在仓库根目录执行：

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.rough_run --output build/fem_roughness/default --angular-step 5
```

参数入口在 `rough_run.py`。可设置 `--rms-nm`、`--width-um`、`--seed`、
`--grid-points`、`--smoothing-cells`、`--model gaussian`、
`--correlation-length-um`、`--gap-nm`、`--polarization`、
`--order`、`--mesh-scale`、`--pml-thickness`、`--angular-step`。
其他光学参数可以通过 Python 调用 `run(Parameters(...), Roughness(...), output)` 设置。

默认验证示例是 355 nm、100 nm 球、n=1.59/k=0、基底 n=5.57/k=0.387、
75 度 p 偏振，入射从负 KY 侧朝正 KY/负 Z 传播。
粗糙区域为 200 x 200 nm，RMS=2 nm，5 x 5 高度节点，随机种子 7。
球底位于 z=0 参考平面上方 10 nm，不是距局部表面恒定 10 nm；
旧模式的样本平均面正好在 z=0，高斯模式的单次样本平均面则可能偏移。
这些粗糙度和间隙是本模块的测试条件，**不是正式 UI 的默认粗糙度/接触条件**。
材料常数是参数，不是可直接用于任意波长的色散数据库。

`--rms-nm 0` 退回平面网格，无颗粒散射为零。
`--mesh-scale 0.8` 细化 FEM 网格但保持同一高度面。
改变 `--grid-points` 会改变表面实现，不能当作同一个表面的单纯网格收敛检查。

## 粗糙面的定义

默认 `legacy` 模式保持旧结果含义：随机高度按网格单元平滑、边缘渐隐，
并按照分片线性三角形的投影面积将**每块样本**均值设为零、RMS 设为指定值。
这会抑制有限样本的低空间频率，不能当作有指定 PSD 的晶圆统计模型。

`gaussian` 模式以 RMS `sigma` 和横向相关长度 `L` 定义未截断高度场：

```text
C(r) = sigma^2 exp[-(r/L)^2]
S(fx,fy) = pi sigma^2 L^2 exp[-(pi L)^2(fx^2+fy^2)]
```

程序按该协方差生成相关高斯高度，然后在有限斑块边缘渐隐以接回平面晶圆。
输入 RMS 是**渐隐前平稳随机表面**的统计 RMS，不强制每块样本均值为零，
因此单个 FEM 斑块的平均高度、RMS 和低频起伏会自然变化，实际 RMS 通常低于输入值。
边缘渐隐改变了
最终有限斑块的 PSD；它不是无限平稳表面的精确 PSD 实现。要得到可置信的
统计 haze，需要多随机种子平均并做斑块大小、表面采样和 FEM 收敛检查。

`abc` 模式增加了与 NIST `ABC_PSD_Function` 一致的高频幂律尾部：

```text
S_ABC(f) = A / [1 + (B |f|)^2]^(C/2)
A = sigma^2 B^2 (C-2) / (2 pi),  C > 2
```

这里 `f` 的单位是 1/um，`A` 的单位是 um4，`B` 是 PSD 转折长度，
**不是** Gaussian 协方差的 `L`。此归一化满足二维积分
`integral S_ABC(f) d^2f = sigma^2`。高度节点由该 PSD 的二维逆 Fourier
变换协方差生成；令 `a=(C-2)/2`、`z=2 pi r/B`，相关函数是
`2^(1-a) z^a K_a(z)/Gamma(a)`（在 `r=0` 处取 1）。`C=3` 时它精确简化为
`exp(-2 pi r/B)`。每张有限表面仍用相同的边缘渐隐和三角网格表示。

ABC 的高频尾部对采样更敏感。报告里的 `nyquist_variance_fraction` 是
无限 ABC PSD 落在高度网格 Nyquist 内切圆内的方差比例下界；低于 0.9 会警告，
但它**不代表**剩余高频被正确忽略，离散节点会产生混叠。`C<=4` 时无限
ABC PSD 的 RMS 斜率不收敛，有限 FEM 网格提供了人为高频截止，因此斜率和
散射结果必须做网格收敛检查。对 266 nm、2 nm RMS，`sigma/lambda=0.0075`，
可作为一阶光滑面参考的探索条件；当 `sigma/lambda>=0.05` 时程序给出警告，
不是把这个近似界线当成严格判据。

一个可重现的 ABC 探索例子为 266 nm、75 度 p 入射、2 nm RMS、
`B=0.2 um`、`C=4.5`、0.2 x 0.2 um 斑块、9 x 9 高度节点、每边
0.025 um 渐隐、16 张表面。半球 DSC 为 FEM `1.845e-5`、有限斑块解析
期望 `1.824e-5`、NIST 无限 ABC `1.945e-5` um2；FEM/有限参考
`1.012`，FEM/NIST `0.949`，FEM/NIST 角分布余弦相似度 `0.987`。
Nyquist 内切圆覆盖 97.1% 的无限 PSD 方差，平均实际斑块 RMS 仍只有
1.55 nm；这**不是**网格、PML、斑块尺寸和样本数收敛证明。

ABC PSD 的形式和参数定义见
[NIST SCATMECH PSD_Function](https://pages.nist.gov/SCATMECH/docs/psd.htm)；
一阶矢量微粗糙度模型见
[Microroughness_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/microroughbrdf.htm)。

### Haze-only 随机表面集合

FEM 页的 `Average haze` 只计算无颗粒的 Gaussian 或 ABC 粗糙面，默认对 16 张独立表面
分别求解，再对远场 DSC **强度**求平均；不会平均复振幅。`Surfaces`、
`Edge taper`（每边渐隐宽度，um）和 `Compare area`（FEM/NIST 共同参考面积，um2）
可单独设置。NIST 曲线由相同 RMS、PSD 长度、ABC 指数（适用时）、入射条件
及参考面积的对应 `Gaussian_PSD_Function` 或 `ABC_PSD_Function` 生成。
程序还计算实际 FEM 三角网格和边缘渐隐对应的
**有限斑块解析期望 PSD**。输出的 FEM 均值、有限斑块参考、NIST PSD 图采用相同坐标和
绝对 DSC 色标，并附有逐角标准误差图、半球 DSC、总量比和角分布相似度。
原始每张表面的 DSC 仍单独保存。
纯 haze 路径不要求粗糙峰值低于虚拟颗粒底部；有颗粒求解仍执行该间隙检查。

设未渐隐的节点高度向量为 `h`，其协方差 `C_ij` 来自所选 PSD 的
二维逆 Fourier 变换（Gaussian 为
`C_ij = sigma^2 exp[-(|xi-xj|/L)^2]`）；`w_i` 为节点渐隐权重，
`b_i(f)` 是三角形分片线性节点基函数的二维 Fourier 积分。程序直接计算

```text
S_finite(f) = E[|sum_i w_i h_i b_i(f)|^2] / A_eff
            = sum_ij w_i w_j C_ij b_i(f) conj(b_j(f)) / A_eff
A_eff = integral (sum_i w_i phi_i(x,y))^2 dx dy
```

其中 `A_eff` 是渐隐窗平方的积分。有限斑块参考把此 PSD 乘以 NIST
`Unit_PSD_Function` 的一阶矢量散射响应；因此它仍是 Rayleigh-Rice **一阶**
参考，不是第二次 FEM 求解。它显式计入形貌节点间距、三角形线性插值、有限宽度
和边缘渐隐，但不包含 FEM 网格/PML 误差或高阶粗糙面散射。平均后的 FEM DSC
同样先按 `A_eff` 换算到 `Compare area`；这不是调整高度或归一化颜色。

此前 Gaussian 例子采用 266 nm、75 度 p 入射、2 nm RMS、50 nm 相关长度、0.2 x 0.2 um
斑块的 16 张示例中，半球 DSC 分别为 FEM `2.006e-5`、有限斑块解析期望
`1.907e-5`、无限 NIST PSD `2.343e-5` um2。FEM/有限参考为 `1.052`，
有限参考/NIST 为 `0.814`，FEM/NIST 为 `0.856`；FEM/NIST 角分布相似度
为 `0.988`。这说明该例的总量差异主要来自**有限形貌表示**，但尚
**不构成严格一致或物理收敛的证明**。该斑块仅为 4 个相关
长度，渐隐后的平均实际 RMS 为 1.57 nm，网格、PML、斑块宽度和样本数
都需要继续收敛检验。仅 1-2 张表面的图像可能明显偏斜，不宜用来判断
统计模型是否吻合。

例如使用 200 nm 宽斑块、50 nm 相关长度、每轴 9 个高度节点：

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.rough_run --model gaussian --correlation-length-um 0.05 --grid-points 9 --output build/fem_roughness/gaussian
```

GUI 中选择 `Rough surface model = gaussian` 会把默认高度节点数提高到 9，
并明确把 FEM 档位设为 `preview`、阶次设为 `1`，以适配当前电脑的内存。
这只是可运行的预览计算，不代表一阶网格已达到物理收敛；用户可手动提高阶次并重新校验资源与结果。
至少要求每个相关长度有两个网格间隔；更多节点会显著提高 FEM 耗时和内存。
`report.json` 的 `surface_statistics` 记录目标模型、采样间距和该次样本的实际均值/RMS。
`metrics.haze_angular_diagnostics` 另记录纯粗糙面散射峰值与中心光线反射方向
`(theta=入射角, phi=90 度)` 的偏离和强度比；它不把镜面反射加入 haze。
高斯 PSD 的定义与 NIST SCATMECH 的
[Gaussian_PSD_Function](https://pages.nist.gov/SCATMECH/docs/psd.htm) 一致，
但这里是 FEM 有限斑块的形貌生成器，不是直接调用 SCATMECH 的 Rayleigh-Rice BRDF。
选择 `abc` 同样采用可调整的高度节点数；界面中的长度框此时表示 NIST 的 `B`，
另有指数 `C`。当前 PC 上高 `grid_points` 与高 FEM 阶次可能很耗内存，
先看 Nyquist 方差比例与有限斑块参考，再逐步加密网格。

`smoothing_cells` 是生成器的平滑尺度，不是经过标定的物理相关长度。
旧模式是有限区域的一次实现，不是无限平稳粗糙面，也不匹配现有 SCATMECH
ABC PSD 的统计模型。只有 RMS 相同不足以比较两种粗糙面的散射结果。
粗糙区域宽度不是有效照明宽度，当前激励仍为无限平面波。

## 方程与远场

时间约定 exp(-i omega t)，参考介质为 z>=0 的空气、z<0 的平面基底。
相对参考介质的介电常数差为：

| 区域 | delta_eps |
| --- | --- |
| 空气中的球 | eps_particle - 1；无颗粒计算时为 0 |
| z>0 的基底凸起 | eps_substrate - 1 |
| z<0 的空气凹陷 | 1 - eps_substrate |

```text
curl curl E_s - k0^2 eps_actual E_s = k0^2 delta_eps E_background
E_total = E_background + E_s
A_alpha(k) = k0^2/(4*pi) integral delta_eps E_total(r) dot E_recip,alpha(r) dV
DSC(k) = |A_s|^2 + |A_p|^2
```

体积分使用曲面四面体的原生积分点及 Jacobian，不用规则空间采样去近似纳米薄层。
互易投影不对电场取共轭。空气中的互易波是直达波和 Fresnel 反射波；
基底中的互易波是从空气入射的透射波，包含吸收衰减。
这允许凹陷内的体源辐射到上半球，不能将它当成空气中的偶极子。

具体地，对观测方向 k=(kx,ky,w)，q=sqrt(eps_substrate-kx^2-ky^2)，
取被动介质分支；s=(-ky,kx,0)/rho，p+=k cross s。
基底内的互易波相位为 exp[-i*k0*(kx*x+ky*y+q*z)]，偏振系数为：

```text
t_s = (1+r_s) s
t_p,xy = (1-r_p) p+,xy
t_p,z = -2 (kx,ky) dot p+,xy / (eps_substrate*w+q)
```

无折射率反差时直接退回自由空间表达式，避免掠射极限 0/0。
已验证界面处切向 E 连续、法向 eps*E 连续以及自由空间平移相位。
这些单元测试不是完整的粗糙面物理准确性验证。

体极化形式可参考 [Martin & Piller (1998)](https://www.epfl.ch/labs/nam/wp-content/uploads/2023/05/031.pdf)。
分层背景近远场转换的互易方法可参考
[Yang, Hugonin & Lalanne (2016)](https://arxiv.org/abs/1510.06344)。
本模块使用半空间体源重叠积分，不是复现该论文的全部算法或导模功能。

## 成对结果的含义

设 A1 是有颗粒振幅，A0 是相同表面无颗粒振幅，dA=A1-A0：

```text
|A1|^2 = |A0|^2 + |dA|^2 + 2 Re[dA dot conj(A0)]
```

- `total`：有颗粒时的散射截面。
- `background`：无颗粒时的粗糙面散射截面。
- `change_field`：两次复电场之差的平方，包含颗粒对表面响应的改变。
- `interference`：带正负号的干涉项。
- `intensity_difference`：半球内 total-background，不等于 change_field。

这里的远场 total 仍是**相对于平面参考背景的散射场**，不包含无限平面波的
入射束和镜面反射束。近场 total 则包含参考入射/反射/透射电场。
结果没有逐图最大值归一化；近场除以入射场幅值，色标使用平方根刻度增强可读性。
DSC 单位 um^2/sr，半球积分单位 um^2，后者不是能量百分比。
给定垂直于传播方向的入射辐照度 I 后，积分截面乘 I 才是对应散射功率。

## 输出和边界

- `rough_comparison.png`：高度面、近场、远场切线与积分分解。
- `paired_fields.npz`：复振幅、复 E/H、高度三角形和四项 DSC。
- `report.json`：完整参数、网格、残差、体源积分检查和半球积分。

结果存放于 Git 忽略的 build 目录，源码与说明不会被忽略。
加密网格的直接矩阵求解可能占用数 GB 到十余 GB 内存，DOF 上限不是内存保证。
当前不包含薄膜、多层基底、大面积统计平均或有限光斑。
近场界面值、尖角峰值、PML/域大小、总能量平衡和独立粗糙面基准仍需验证。
三路独立理想物镜与 20x 成像链路现已作为实验扩展接入，见 [IMAGING.md](IMAGING.md)。
该扩展不代表上述物理收敛检查已经完成。

## 自动收敛扫描

FEM 页新增 `Convergence` 按钮，读取当前物理参数、三路 NA/倾角、
边缘渐隐宽度和比较面积。弹窗中每行填写按比较顺序排列的数值，
逗号分隔；清空某行可省略该扫描。建议首先使用 4 至 8 张表面，
确认资源和数值变化后再增加样本数。GUI 的时间限制是整个任务的限制。

扫描项包括 FEM 网格比例、PML 厚度与强度、物理域半径、高度节点数、
斑块宽度、渐隐宽度、远场角采样、体源积分阶数及 `fem_order` 多项式阶次。
FEM 网格、PML、域半径、阶次和积分阶数的比较使用相同种子与完全相同的高度面，
报告对每个样本能量差的成对 95% Student t 区间。斑块宽度扫描保持
物理渐隐宽度和不增大的高度网格间距，因此节点数可能增加；
它和渐隐扫描报告均值、误差条和敏感度，不采用相同表面的成对判据。
Gaussian/ABC 高度网格扫描现在从所有等级节点的联合协方差生成一张共享
母表面：公共节点完全一致，新节点保持目标协方差，不是对粗网格简单插值。
报告的 `paired_same_realization` 区别于 `paired_same_surface`；三角形表示仍会变化，
但可以对同一母表面的离散误差计算成对区间。实测 PSD 高度扫描暂不采用这个判据。
联合节点最多 2048 个；细网格建议使用嵌套等级，例如 9、17、33。

最后两个扫描等级用于检查：

```text
(abs(mean paired energy difference) + paired 95% half-width)
    / abs(last-level mean energy) <= requested tolerance
```

数值检查通过仅表示指定等级、指定输入和指定孔径范围内满足容差。
随机表面至少需要三个样本才允许该检查通过；单个实测高度面属于确定性
数值比较，不产生随机表面误差条。当前不自动证明多参数联合作用或能量守恒。
角采样扫描复用保存的体源，不重复 FEM 求解，但阶跃孔径掩膜的积分仍需
通过角采样加密确认。三路独立积分，NA=0 的一路为零，不求三路并集。

每个样本独立受内存和时间保护，串行计算。缓存以物理/数值参数、种子、
实测文件 SHA256 和求解模块源码指纹识别；相同计划可续算完成的样本。
结果保存为 `haze_convergence_report.json`、`haze_convergence.csv` 和
`haze_convergence.png`，包含半球/三路积分、95% 区间、自由度、耗时和峰值内存。
异常时保留已完成样本及失败报告。某等级超出资源限制或网格生成失败时，
整组等级排除，不能仅平均成功的表面造成选择偏差；继续其余扫描并标注
`complete_with_limits`。它代表扫描过程结束，不代表所有等级成功或物理收敛。
可在 UI 的 `Load result` 中加载扫描目录。在本机的一阶验证中，
`sparsecholesky` 比 UMFPACK 显著节省内存；需要时在 FEM 输入的
Direct solver 中选择它，再运行扫描。

命令行使用 JSON 计划：

```powershell
& 'build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.haze_convergence --input 'build\my_plan.json' --output 'build\my_study' --max-memory-gib 10 --sample-timeout-seconds 180
```

计划包含 `parameters`、`roughness`、`samples`、`numerical_apertures`、
`tilt_degrees`、`comparison_area_um2`、`tolerance`、`target_relative_ci`
及 `scans`；可从 UI 产生的 `convergence_plan.json` 修改。

`Average haze` 现在同时保存 `haze_sample_convergence.png`，在 2、4、8、16...
个样本前缀上报告半球与三路能量的点态 95% Student t 区间。
这是有限样本近似，依赖独立实现，不是三路联合置信保证。
一个样本不产生虚构标准误差，也不能声明统计收敛。

## 实验高度图与 PSD

FEM 的 `Rough surface model` 新增 `measured_height` 和 `measured_psd`。
选择文件后，RMS 来自测量数据；原有 RMS/B/C 不对实测振幅重新缩放。
高度图的行方向为全局 Y，列方向为全局 X，入射传播的面内方向为 +KY。

| 输入 | 文件内容 | 单位 |
| --- | --- | --- |
| 高度 NPZ | `height_nm` 二维矩阵、`dx_um`、`dy_um` 标量 | nm、um |
| 高度 CSV | 无表头的二维高度矩阵；同名 JSON 必须有 `dx_um`、`dy_um` | nm、um |
| 二维 PSD NPZ | `fx_per_um`、`fy_per_um` 一维均匀、含零的 fftshift DFT 频率轴，`psd_um4` 矩阵按 Y/X 排列 | cycles/um、um4 |
| 径向 PSD CSV | 两列 `f_per_um,S2_um4`，频率非负且严格递增；表头用 `#` 注释 | cycles/um、um4 |

高度处理可选 `none`、`mean`、`plane`。缺失点默认拒绝，显式选择
`nearest` 才允许填补，且缺失比例不能超过 5%；正负无穷仍拒绝。
高度图可以是长方形；FEM 从中心截取指定宽度的正方形，双线性插值到
指定高度节点数，再施加明确的边缘渐隐，不能超出测量范围。
导入后的几何仍受当前每轴 4..41 个节点的限制。
高度图和二维 PSD 输入均支持每轴最多 1024 个测量点；FEM 的节点数限制
是另一项设置。复杂表面的 OCC/Netgen 布尔操作与建网格可能失败，
失败时保持原始测量不变并检查任务日志，不把失败的表面混入成功样本平均。

PSD 必须非负并满足实高度的 `S(-fx,-fy)=S(fx,fy)`。
积分采用 DFT 频率单元求和，得到已测带宽内的方差；未测频段设为零，
不外推，也不补齐到任意指定 RMS。径向 CSV 会离散成各向同性二维频谱；
二维 NPZ 保留各向异性。PSD 合成采用高斯随机高度假设，不能由 PSD
推断高度概率分布或恢复原始实测表面。

同名 JSON 可记录仪器、校准、测量带宽和传递函数信息；导入报告保存
完整元数据、像素尺寸、去趋势系数、缺失点处理、频率范围和 SHA256。
当前仪器传递函数只记录，不执行未知响应的自动反卷积；应导入经校准的 PSD。

实测高度直接进入近场/远场/成像流程，也可用 `samples=1` 执行收敛扫描。
它不能通过换种子变成多张独立测量。实测 PSD 可以参与随机表面平均，
并使用 NIST Unit PSD 的电磁响应乘上该频谱，以及相同三角形节点协方差
的有限区域参考。对于单张实测高度，收敛报告保留实际积分截面，
不按假定无限平稳面面积放大。

输入审计可单独运行：

```powershell
& 'D:\anaconda\envs\NIST2\python.exe' -m fullwave_prototype.experimental_input --input 'D:\measurements\wafer.npz' --kind height --detrend plane --output 'build\wafer_input_audit'
```

输出处理后的高度、二维 PSD、输入图和 JSON 审计。
高度到 PSD 使用无额外窗的双边 periodogram，验证
`sum(PSD)*dfx*dfy = mean((height_nm/1000)^2)`，报告 Parseval 误差。
有限测量窗会造成频谱泄漏；该 PSD 不是无限表面的精确统计真值。
测试脚本 `tools/prepare_haze_validation.py` 产生的输入是明确标注的合成数据，
用于导入流程验证，不是实测晶圆证据。

## 几何与孔径积分的数值改进

粗糙面 CAD 构造与布尔运算现在在 nm 尺度进行，之后将线性网格和 CAD
几何恢复到 um，再构造高阶曲面映射。原始测量高度、PSD、光学波长和
Maxwell 求解单位不变。不要只缩放网格顶点而保留旧 CAD 几何，否则曲面映射
会不一致。实现已独立到 `rough_geometry.py`。
当前验证修复了 0.25 um ABC 斑块的边界重叠及一个含颗粒的规则周期高度样例。
这不保证任意复杂表面都能成功；资源保护和失败样本整组排除仍然保留。

几何体积审计需要先禁用 PML；PML 的人工拉伸不是实际物理体积。
可独立验证保存的参数：

```powershell
& 'build\ngsolve_env\Scripts\python.exe' -m tools.rough_mesh_audit --input 'build\sample\input.json' --output 'build\mesh_audit'
```

三路统计能量现在使用圆锥坐标下的 Gauss 立体角积分，越过半球地平线的
部分准确截除。地图数据使用周期双线性 DSC 插值，阶数为 20；仍需检查
角采样收敛，插值不产生原地图缺少的物理细节。NA=0 的一路严格为零。
不再用离散地图点的硬孔径掩膜做这些统计积分，避免圆锥边界台阶误差。
半球积分仍按原有球坐标采样做数值积分。

另有 `tools.aperture_convergence_audit` 直接对保存的复 FEM 体源在孔径积分
节点重新计算辐射场，用作独立角积分参考，不需要重做 Maxwell 求解。
它不是独立电磁算法，因此不能证明 FEM 场本身已准确。
收敛 JSON 还提供 `sample_budget_estimate`：基于先导样本波动估计达到指定
区间半宽需要多少独立表面。它不是收敛保证，也不是反复查看普通 t 区间
即可停止的规则；当前单次任务仍限制最多 128 张。

## 二阶局部与边界收敛

收敛计划支持三个额外的成对扫描轴：

- `rough_mesh_scale`：只缩放粗糙凸起/凹坑的 FEM 网格尺寸，不改变高度节点或 PSD。
- `pml_mesh_scale`：只缩放吸收层的网格尺寸，不改变物理区域的尺寸设置。
- `pml_order`：检查降低吸收层阶次是否引入偏差，取 1..物理区域阶次的整数。

吸收层降阶只是资源优化，不能默认认定没有数值误差。厚度、强度、
计算域半径和 FEM 阶次仍可分别扫描。局部网格尺寸会通过网格生成器
的过渡影响邻近单元；它不是固定其他单元后的严格嵌套网格。

使用 `--sample-cache` 可以把验证报告放在新目录，同时复用旧目录内的
样本。复用要求输入和求解源码指纹完全匹配，且资源记录为成功完成、
退出码为零、无存活的任务子进程；未完成的样本不会被接受。
报告的 `sample_folders` 保留每张表面的实际数据位置。

```powershell
& 'build\ngsolve_env\Scripts\python.exe' -m tools.prepare_second_order_study --input 'build\model_optimization_20261008\optimization_plan.json' --output 'build\model_optimization_20261009\local_plan.json' --stage local
& 'build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.haze_convergence --input 'build\model_optimization_20261009\local_plan.json' --output 'build\model_optimization_20261009\local' --sample-cache 'build\model_optimization_20261008\study' --max-memory-gib 10 --sample-timeout-seconds 480
```

`--stage boundary` 使用完整二阶 PML，生成强度、厚度和域半径计划。上述路径是本次
先导验证的示例；工具需要已保存的输入计划，不会自动选取某张实测晶圆。

新增 `haze_paired_errors.png`，黑点表示均值变化，条形表示成对差值的
点态 95% 误差上界，虚线表示设定容差。绿色只表示最后两个指定层级
满足当前判据，不表示全部离散参数或实际晶圆已经收敛。缺失或资源
超限的层级明确标记，不会画成零误差；旧结果无此图时 UI 隐藏该页。

旧的单表面 `haze_crosscheck` 也已改用同一圆锥积分；其半球积分和
NIST 电磁模型没有更换。旧报告中的网格掩膜数值不会自动改写，需重新
生成报告。这样可避免新旧对比工具的孔径边界误差被误认为模型差异。
