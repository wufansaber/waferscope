# 粗糙表面与颗粒的相干全波原型

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

可选 `gaussian` 模式以 RMS `sigma` 和横向相关长度 `L` 定义未截断高度场：

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
