# 整光斑粗糙面与相干 haze

## 入口与范围

源码 FEM 页新增 `Whole-beam haze (SPM)`。它覆盖主要照明区域，不再将
0.2 um 粗糙斑块以外的大部分受照晶圆当成平面，也不把局部 FEM 的
功率乘面积比例来冒充整光斑结果。原 `Calculate` 和 `Average haze`
仍是局部几何的全波 FEM；新入口是**一阶矢量表面微扰 SPM**，两者分开标注。

输入沿用同一组波长、基底 n/k、入射角、偏振、Gaussian 光斑长宽、
总功率、Gaussian/ABC PSD、RMS、相关长度 B/L、ABC 指数、随机种子，
以及三路独立 NA/倾角、20 倍成像参数。颗粒参数不参与这个 haze-only 模型。
默认总功率 1 W、光斑目标直径 2×2 um、覆盖目标 0.99、
Nyquist 圆内 PSD 方差目标 0.90；默认平均 16 张独立表面。

需要同时加入颗粒时，使用新的 `Particle + whole haze (hybrid)`，
见 [相干混合模型](HYBRID_MODEL_CN.md)。它使用一张共享母表面与局部 FEM 修正，
不是本入口的 16 张 haze-only 平均，也不是整片晶圆的严格全波解。

界面左侧可设置 `Whole-beam power coverage` 与
`Whole-beam PSD variance fraction`。0.99 表示 99%，不是 0.99%。
Edge taper 控制边缘渐隐的目标宽度，实际值至少跨两格，并写入结果。
软件自动选择粗糙面范围和 FFT 网格，不使用局部 FEM 的斑块宽度/域半径
作为整光斑尺寸。当前每轴最多 2048 点；超过上限时明确拒绝计算，
不会暗中降低 PSD 分辨率或按样本 RMS 补回未采样频段。

目前支持 Gaussian 和 ABC 统计表面及 Gaussian 有限照明。实测高度/PSD、
平面波无穷照明以及相机噪声不在新入口支持范围内，错误会明确报告；
这些参数在原 FEM/SCATMECH 流程中的已有功能保持不变。

## 覆盖不是几何面积比

Gaussian 光斑没有有限的严格边界，有限范围不能宣称覆盖 100%。先在
规划网格上扩展表面，再在最终细网格上检查：

```text
eta = integral [-Re(Einc cross conj(Hinc*Z0))_z] w(x,y)^2 dxdy / F_beam
```

`w` 为粗糙高度的边缘渐隐函数，内部等于 1，外边缘连续降到 0。
采用 `w²` 检查，防止只统计几何白框却忽略渐隐区对高度方差的削弱。
程序要求 eta 达到用户设定目标；结果同时保存未乘渐隐权重的几何覆盖。
它们都按真实光场能量分布积分，不是窗口面积比，也不是散射收集率。
窗口外剩余照明尾部仍未建模，不外推、不重归一成 100%。

照明通过固定相位的复数角谱相干叠加，规划时在整个表面检查角谱
阶次变化。必要时使用 41、61、81、101、151 阶；达不到 0.1% 相邻
复场变化目标时停止。全部复场仍以整束 1 W 的 Poynting/Parseval
积分归一，而不是按扩大后的表面面积重新定义入射辐照度。

## 随机面与电磁核

二维双边 PSD 采用 cycles/um、um4。白高斯噪声在频域乘
`sqrt(S)/(sqrt(dx*dy))` 后逆变换，产生统计平稳随机高度。
在双倍宽度的父域上合成并截取中央区域，避开周期拼接边界。
高度的 DC 与每张表面的统计 RMS 自然波动，不逐张强制零均值/固定 RMS。
有限父域、有限高频截止和最终渐隐仍需要数值收敛检查。

表面位移的一阶辐射公式为：

```text
T(k_parallel) = integral h(x,y) w(x,y) Ebg(x,y,0+) exp[-i k0 k_parallel·(x,y)] dxdy
A_alpha(k) = k0² (eps_s-1)/(4 pi) b_alpha(k) · T(k_parallel)
b_s = (1+r_s) s
b_p,parallel = (p_up + r_p p_down)_parallel
b_p,z = (p_up + r_p p_down)_z / eps_s
DSC = |A_s|² + |A_p|²
```

`Ebg` 是有限照明的入射加 Fresnel 反射总背景。
法向项的 `/eps_s` 来自法向 D 连续，不能把全部 E 都视为界面连续。
先对整个粗糙面求复振幅，再取平方；不同位置不是独立强度相加。
不同随机实现只在求完各自的散射/像面强度后进行统计平均。

同类理论来源是 NIST 的
[Microroughness_BRDF_Model](https://pages.nist.gov/SCATMECH/docs/microroughbrdf.htm)，
即单界面光滑极限的一阶矢量 Rayleigh-Rice 微扰。
新核在平面波、相同 n/k 和单位 PSD 下，与 NIST 绝对 DSC 的
s/p 偏振交叉校验最大相对误差约 1.6×10^-14。
这验证了单位、偏振与核，不代表任意粗糙面的一阶近似都成立，
也不意味着有限 Gaussian 光束等同于无限平面波的 NIST 统计图。

当前入口要求 RMS/波长 < 0.05；该门槛只是拒绝明显不适合的输入，
不是理论误差保证。ABC C<=4 的无限带宽斜率方差不收敛，会明确警告；
格点高频截止仍会影响结果。模型无多层膜、颗粒耦合、粗糙面高阶
多次散射或镜面光束，不能标为整片晶圆的严格 Maxwell FEM 解。

## 远场、成像与误差

总功率采用直接表面谱的 Gauss 立体角积分，半球和三路圆锥分别计算，
不是对 2.5° 显示图的像素或插值求和。积分检查 48/64/96/128 阶，
相邻变化目标为 1%；到上限仍未达到时保留明确警告。
每个随机样本还在选定方向比较 FFT 插值与独立空间 Fourier 求和，
复振幅相对 L2 超过 1% 时拒绝输出成功结果。选定点检查不是全角度误差保证。

三路复数出瞳直接来自整面辐射谱，使用现有无像差矢量传播到独立像面。
像面显示绝对 W/um²，三路共用色标，不逐路归一；平均三路图保持各自
功率与随机实现的非相干平均。NA=0 关闭该路。
远场收集功率与带透过率/偏振片的像面功率含义不同，不能混用。
三路重叠方向也不能简单相加成为唯一收集功率。
有限 FFT 像面是周期图；视场不够时会警告，增加 padding 不能代替出瞳采样。

输出 `whole_haze_report.json`、`whole_haze_fields.npz`、
`whole_surface.png`、`whole_haze_far.png`、`whole_haze_images.png`
及 `whole_haze_statistics.png`，可用 Load result 再次打开。
统计页提供点态 95% Student t 区间；达到相邻数值积分目标不等于
达到随机表面统计目标，更不等于通过物理/实验精度认证。

## 本机运行证据

2026-10-09，266 nm、75°、p 偏振、1 W、2×2 um 目标光斑，
Gaussian RMS 2 nm、L=0.05 um，三路 NA=0.5、倾角 +45/-45/0°，
128 点出瞳、20 倍成像，平均 16 张随机表面：

| 检查 | 结果 |
| --- | --- |
| 自动粗糙面 | 7.5×7.5 um、256×256 点 |
| 网格间隔 | 0.029296875 um |
| 实际渐隐宽度 | 0.05859375 um |
| 渐隐权重后的照明覆盖 | 99.2938% |
| Nyquist 圆内目标 PSD 方差 | 99.9700% |
| 整面复照明角谱变化 | 0.00924% |
| 平均半球散射功率 | 2.13285 mW |
| Ch1/Ch2/Ch3 远场收集功率 | 0.38940 / 0.36717 / 0.47641 mW |
| Ch1/Ch2/Ch3 像面平均功率 | 0.39069 / 0.36770 / 0.47712 mW |
| 单样本直接角积分检查 | 16/16 达到相邻 1% 目标 |
| 记账内存峰值 / 总耗时 | 约 0.25 GiB / 12.75 s |
| 统计 95% 区间相对半宽 | 半球 3.56%，三路 8.03% / 5.90% / 6.59% |

结果目录：`build/whole_haze_20261009/direct/result`。
半球已达到当前 5% 统计半宽目标，三路尚未全部达到，软件不会宣称
整体统计收敛。扩大样本数仍需检查表面范围、PSD 带宽和成像采样。
材料 n/k 在该示例中沿用输入参数，不是 266 nm 已标定色散。

运行命令：

```powershell
& 'build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.guarded_run --report build/whole_haze_20261009/direct/resources.json --max-memory-gib 10 --reserve-gib 3 --timeout-seconds 600 whole_haze --input build/whole_haze_20261009/direct/whole_haze_input.json --output build/whole_haze_20261009/direct/result
```

上述 input JSON 由 `tools.whole_haze_smoke` 从实际 UI 控件生成。
NIST2 完整回归 195 项通过、16 项跳过；九项新物理/数值测试包含
NIST 核、PSD 绝对振幅、相干网格、范围覆盖及独立 Fourier 积分。
独立求解环境中八项新测试通过，NIST 核测试在该环境跳过，已由 NIST2 验证。
ABC B=0.2 um、C=4.5 的两样本整光斑流程也在本机通过，输出位于
`build/whole_haze_20261009/abc/result`，两样本不构成统计收敛证据。
源码已实现，既有 EXE 不会自动包含新增功能，需另行打包。
