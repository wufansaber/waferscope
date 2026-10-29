# PSL 球形标准颗粒：验证范围、证据与参考文献

更新日期：2026-10-27。

## 1. 可以声明什么

本阶段将正确性验证限定为：**单个、均匀、各向同性、无颗粒包覆层的理想球形 PSL 颗粒，在空气中接触平整硅半空间，受单色平面波照明时的电磁散射。** 首先建立 100 nm 级基准；10～100 nm 系列仍区分“理想材料参数的数值验证”和“真实标准颗粒的实验验证”。

几何接触定义为 `gap=delta=0`、球心高度 `z_c=d/2`。改变直径时保持球底接触。这里不声称真实聚合物与晶圆之间没有形变、残留物或界面薄层。

主目标波长为 266 nm；442/532/633 nm 只作为已发表 PSL-on-wafer 实验的复现条件，不能将可见光实验直接当作 266 nm 的精度证明。

本阶段不外推到随机粗糙面/haze、图形晶圆、光栅、凹坑、十字/线条、非球形、团聚体、多颗粒耦合或未经测量的镀膜。已知氧化层/包覆层是后续另列参数的扩展，不混入无膜基准。

首先验收到绝对复远场、角分布和三路收集功率。近场需要独立检查局部场与界面条件；成像、实际相机灰度和 SNR 另立验证链，不能因散射通过而自动认为它们也通过。

**当前不能声明：本程序已在 266 nm、整个 10～100 nm PSL 标准颗粒范围内获得实验认证。**

## 2. 不依赖人员判断的证据链

不是用“SCATMECH 来自 NIST”“FEM 是成熟方法”或者“两个结果接近”代替证明，而是预先固定输入与判据，由外部解析解、物理定律、独立算法和可追溯样品共同约束。

| 层级 | 要回答的问题 | 证据 | 不足以证明的内容 |
|---|---|---|---|
| 理论 | 是否求解同一个物理边值问题？ | Maxwell 方程、界面/辐射条件、BV 推导 | 不保证某次程序调用正确 |
| 实现 | 方程和归一化是否实现正确？ | 独立 Mie/Bessel、偶极子、Fresnel 极限及完整复场对照 | 不保证具体大网格已收敛 |
| 数值 | 当前结果离极限解多远？ | 多极阶数、网格阶数/密度、计算域、PML 和积分收敛 | 小矩阵残差不能替代这些检查 |
| 实验 | 是否代表实际 PSL 标准样品？ | 有效批次证书、波长相关 n/k、晶圆状态、绝对散射测量和不确定度 | 粒径证书不认证散射、成像或 SNR |

严格数学误差界与工程数值验证也不同。当前尚无覆盖本项目接触几何、材料和 PML 的完整严格误差界，不能将收敛步长变化直接称作已证明的全局误差上界。

### 2.1 SCATMECH 为什么有理论依据

BV 模型用球的多极响应以及平面基底的反射/平移耦合描述球与基底的相互作用。[1] 其结构可示意为：

\[
\mathbf a=\mathbf T(\mathbf b+\mathbf R\mathbf a),\qquad
(\mathbf I-\mathbf T\mathbf R)\mathbf a=\mathbf T\mathbf b.
\]

其中 T 表示球的散射响应，R 包含基底反射与相关展开变换。它不是仅把自由空间强度乘上一个 Fresnel 系数。

NIST 文档明确要求检查收敛；`delta=0` 表示接触，完整相互作用使用 `order=-1`，准确计算应关闭 `Norm_Inc_Approx`，并检查 `lmax` 与数值稳定性。[2] 这说明理论适用性，不是对所有输入及本软件封装的自动认证。

### 2.2 FEM 为什么可以独立对照

在非磁性、线性材料及 `exp(-i omega t)` 约定下，背景场 Eb 已满足平面晶圆问题时，散射场满足：

\[
\nabla\times\nabla\times\mathbf E_s-k_0^2\varepsilon_r\mathbf E_s
=k_0^2(\varepsilon_r-\varepsilon_{r,b})\mathbf E_b.
\]

同时必须满足材料界面条件和出射辐射条件。FEM 使用 H(curl) 空间离散这个边值问题，PML 则是有限计算域对开放空间的数值处理。相关数学框架见 Monk。[5]

**采用这个弱式，不等于当前网格已正确解析该解。** 代数残差小只表示离散方程解得准确；还必须排除几何、网格、积分及 PML 误差。

## 3. PSL 范围内的验证顺序

| 编号 | 基准 | 比较量与限制 | 当前证据状态 |
|---|---|---|---|
| A | 同一颗粒去掉基底折射率差，退化为均匀介质球 | 独立 Mie 系列的散射/消光截面、角分布、近场与相位；不仅看 d^6 | 已有直接 Bessel 与独立 miepython 组件测试；不等于全部 PSL 尺寸已逐项验收 |
| B | 无颗粒平整晶圆 | 复 Fresnel 反射/透射及背景场；无颗粒散射增量应为零 | 无颗粒零散射路径已检查；当前材料色散仍待追溯 |
| C | PSL 接触球/平整晶圆 | BV 与独立 FEM 的绝对 DSC、复远场、三路积分功率；只消除无意义的整体常数相位 | 已有 266 nm 同参数结果；FEM 计算域/PML 收敛仍未通过 |
| D | 物理收支与边界 | 检查物理控制面能流、吸收和界面条件；数值 PML 吸收不算样品吸收 | 不以 FFT/Parseval 或矩阵残差代替整个系统的能量闭合 |
| E | 已发表 PSL-on-wafer 实验 | 先复现文献的粒径、波长、偏振、角度、材料和晶圆状态 | 已核对全文及输入，已复算五组文献条件曲线；尚无逐点实验误差 |
| F | 实际批次 266 nm PSL 样品 | 实际有效证书、n/k、晶圆/氧化层状态和独立标定的绝对测量 | 尚未完成，不能以可见光论文或软件一致性替代 |

去除基底折射率差才是均匀介质 Mie 极限；不能只把球移远就忽略仍存在的反射照明和反射出射波。

对于无吸收背景中的孤立球，Mie 理论给出：

\[
\sigma_{sca}=\frac{2\pi}{k^2}\sum_{\ell=1}^{\infty}(2\ell+1)
(|a_\ell|^2+|b_\ell|^2),\qquad
\sigma_{ext}=\frac{2\pi}{k^2}\sum_{\ell=1}^{\infty}(2\ell+1)
\operatorname{Re}(a_\ell+b_\ell).
\]

这类解析量可用于排查绝对归一化；相关计算方法见 Wiscombe。[4] 无损孤立球应满足消光等于散射；**吸收硅晶圆体系不能只积分上半球就宣称能量闭合**，还要处理基底中的通量、吸收和背景干涉。

近场比较选择有明确位置、离界面距离和分辨率的点或积分量，不把接触点附近未经收敛的最大像素值当作认证指标。

## 4. 真实 PSL 标准颗粒的输入要求

### 粒径证书不等于光学响应证书

NIST SRM 1963a 是一个可追溯粒径标准的历史例子：本次获取的证书给出气溶胶颗粒的众数直径 101.8 nm、扩展不确定度 ±1.1 nm，而不是“每颗球都精确 100 nm”。[6]

该证书版本标明有效期截至 2022-03-31。这里只用它说明被测量和溯源方式，不把它作为 2026 年当前可用批次的有效证明，也不据此判断供应状态。实际验证必须使用手中批次的有效证书。

应记录证书的被测量是众数、数均值、迁移直径还是其他定义；不能用名义值替代。标准颗粒分布、单颗粒差异和沉积后的状态要分开处理。也不能把证书中的粒径分布宽度当作众数认证不确定度。

### 266 nm 材料光学常数必须独立追溯

本轮旧输入为 PSL `n=1.59,k=0`、Si `n=5.57,k=0.387`。它们目前只是固定的对照输入；两个求解器使用同样数值而一致，不证明这些数值适用于真实样品的 266 nm 响应。

French 等提供聚苯乙烯 UV/VUV 光学性质的材料文献入口；Aspnes 与 Studna 提供包含 266 nm 能量区间的硅光学性质研究。[7,8] 但材料形态、样品制备与批次不同，不能未经核对就当作特定 PSL 乳胶颗粒和实际晶圆的精确输入。

特别是不能仅从可见光折射率外推 UV 的 n，也不能未经证据将 k 固定为 0。本次不据摘要或图片估读强行填入新的 n/k，更不通过调整 n/k 来让 FEM 拟合 SCATMECH。

还需记录实际晶圆氧化层、颗粒沉积残留、球形度、团聚、是否发生光损伤，以及功率/光斑/偏振的测量。若它们与无膜理想球假设不符，应报告模型失配，不能悄悄改变模型后仍称为同一基准。

## 5. 本项目已有证据和明确缺口

2026-10-27 新增[外部参考基准实际过程与结果](../output/psl_reference_20261027/VERIFICATION_CN.md)：在 100 nm、266 nm 的**均匀空气**问题上，独立 Mie 与当前 BV 的退化极限一致；本次新算的细网格 FEM 总散射误差 0.927%，能量不闭合仍为 1.324%。这不替代下述接触硅晶圆验证。

- 2026-10-27 在 NGSolve 环境复跑 `test_theory_validation.py` 与 `test_physical_audit.py`：11 项通过，其中包含不依赖 SCATMECH 的 Bessel/Mie 交叉核对。它们主要是理论组件测试，**没有因此完成一次新的 FEM 网格收敛扫描**，也不是全在 266 nm/真实 PSL 参数下完成。
- NIST2 环境同组测试为 9 项通过、2 项因缺少 miepython 跳过；两项独立 Mie 测试在上述 NGSolve 环境实际执行并通过，不把跳过说成通过。
- 接触几何、球心随粒径变化、相位原点已独立检查；只有整体常数相位可以消除，不拟合强度比例、位置或角度相关相位。
- SCATMECH 多极阶数、球冠积分和光瞳采样已有部分收敛对照，覆盖范围见原报告；不能据端点测试声称整个参数空间都已证明。
- 某个 100 nm、75° P 的细化 FEM 算例与同输入 BV 的三路功率差约 1.00%～1.17%，**这不是已认证误差上界**。范围与 PML 的变化仍影响结果，进一步 PML 细化触发内存保护。
- 材料、批次和 266 nm 实验闭环仍缺失。单次接近 NIST、测试项数多、曲线单调或矩阵残差很小，都不能填补这些缺口。

当前预设的“连续两级功率变化小于 1%、复场差异小于 2%”只是项目工程目标，不是文献赋予的通用标准，也不是严格误差界。正式验收应根据应用允许误差和粒径、材料、测量不确定度预先冻结预算，不能看到结果后再放宽。

原始数据和数值限制：[尺寸/亮度报告](../output/particle_brightness_contact_266_20261026/BRIGHTNESS_VALIDATION_CN.md)、[数值审查](../output/particle_brightness_contact_266_20261026/NUMERICAL_REVIEW_CN.md)。本次理论测试记录：[NGSolve 环境](../build/contact_phase_release_20261026/psl_native_theory_tests_20261027.xml)、[NIST2 环境](../build/contact_phase_release_20261026/psl_theory_tests_20261027.xml)。

## 6. 参考文献及证明边界

以下已核对作者/出版信息与可访问的原始摘要、官方文档或证书；**不声称已逐图复现每篇论文**。

1. P. A. Bobbert, J. Vlieger (1986). *Light scattering by a sphere on a substrate*. Physica A **137**, 209–242. [原论文记录](https://www.sciencedirect.com/science/article/pii/0378437186900725)，DOI: `10.1016/0378-4371(86)90072-5`。支持球/平面基底模型的理论来源，不认证本软件的具体输入和数值精度。

2. NIST. [SCATMECH: Bobbert_Vlieger_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)。支持接触参数、完整相互作用/近似开关和收敛检查要求；不是所有场景的无条件金标准。

3. J. H. Kim, S. H. Ehrman, G. W. Mulholland, T. A. Germer (2002). *Polarized light scattering by dielectric and metallic spheres on silicon wafers*. Applied Optics **41**(25), 5405-5412. [作者机构 NIST 记录](https://www.nist.gov/publications/polarized-light-scattering-dielectric-and-metallic-spheres-silicon-wafers-0)，[DOI](https://doi.org/10.1364/AO.41.005405)。PSL 粒径为 101/155/218 nm；92 nm 属于铜颗粒。照明为 442/532/633 nm，实验比较含 1.5 nm 氧化层。它是相关的实验依据，但不能直接证明 266 nm 或 10 nm PSL 的精度；文中关于金属颗粒的具体结论也不能原封不动搬到 PSL。

4. W. J. Wiscombe (1980). *Improved Mie scattering algorithms*. Applied Optics **19**(9), 1505–1509. [原论文重印本](https://indico.ictp.it/event/a02060/contribution/2/material/0/0.pdf)，DOI: `10.1364/AO.19.001505`。支持孤立球的独立系列计算与数值参考，不直接解决接触晶圆或成像问题。

5. P. Monk (2003). *Finite Element Methods for Maxwell's Equations*. Oxford University Press. [出版社：散射理论章节](https://academic.oup.com/book/41656/chapter-abstract/353725057)，DOI: `10.1093/acprof:oso/9780198508885.001.0001`。支持 Maxwell/FEM 的数学与散射框架；书中的收敛理论有相应假设，不是本项目所有几何/PML 配置的自动保证。

6. NIST. [SRM 1963a 官方证书](https://tsapps.nist.gov/srmext/certificates/1963a.pdf)，本次获取版本签发于 2019-05-01、有效期至 2022-03-31。支持历史粒径认证数值与溯源定义，不认证 266 nm n/k、沉积后形状或散射强度。

7. R. H. French et al. (2007). *Optical Properties and van der Waals–London Dispersion Interactions of Polystyrene Determined by Vacuum Ultraviolet Spectroscopy and Spectroscopic Ellipsometry*. Australian Journal of Chemistry **60**(4), 251–263. [出版社原文入口](https://www.publish.csiro.au/CH/CH06222)，DOI: `10.1071/CH06222`。提供聚苯乙烯光学性质的材料参考；不是特定 PSL 标准颗粒批次的 n/k 证书。

8. D. E. Aspnes, A. A. Studna (1983). *Dielectric functions and optical parameters of Si, Ge, GaP, GaAs, GaSb, InP, InAs, and InSb from 1.5 to 6.0 eV*. Physical Review B **27**, 985. [出版社原文入口](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.27.985)，DOI: `10.1103/PhysRevB.27.985`。提供硅光学性质的参考，仍须匹配实际晶圆状态；不证明旧程序常数适用于 266 nm。
