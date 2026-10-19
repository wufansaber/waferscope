# 远场到像面模型正确性校验

校验日期：2026-10-12。生产代码基线：`7def4aa`。

## 结论与边界

**在理想无像差、满足正弦条件、物像两侧为空气的模型假设下，现有复远场到像面的矢量传播实现通过本次校验。** 未发现偏振旋转、放大倍率、像面反转、离焦符号或绝对功率归一化的实现错误。本次增加校验脚本和回归测试，没有修改生产传播算法、UI 或已有 EXE。

不能把这个结论扩展为“任何参数下均已达到 1% 精度”或“已经与真实设备完全一致”。主要限制如下：

1. 32 点光瞳不适合高精度定量。NA=0.96 的本次基准像面剖面误差达到约 3%～4%。
2. 96/128 点的本次基准误差明显减小，但硬圆形孔径的离散边界导致误差非单调，128 不保证每一项都优于 96。
3. SCATMECH 成像入口使用原生振幅焦点参考；FEM/混合模型以晶圆坐标原点为焦点参考。不能直接比较两个入口的零离焦图而忽略坐标转换。
4. 平整晶圆上的单球有独立 BV 参考；粗糙面与有限光斑案例本次验证的是给定复远场的传播正确性，不是重新证明 FEM/SPM 散射源和高阶颗粒-粗糙面耦合全部正确。
5. 没有验证真实镜头像差、装调误差、遮挡/分光、杂散光或实测相机。三路仍是三个独立理想接收器。

## 物理公式与代码

生产内核为 [`optical_propagation.py`](../optical_propagation.py)，SCATMECH 和 FEM 两个成像入口共用它。

设物镜局部方向为 `(u,v,w_o)`，其中 `w_o=sqrt(1-u²-v²)`；空气中 `u²+v²<NA²`。倾斜物镜先旋转到全局坐标，再剔除全局 `k_z<=0` 的方向。积分测度为：

```text
dOmega = du dv / w_o
P_collected = I_ref * integral ||F_far(s,p)||^2 dOmega
```

这里 `F_far` 是复振幅，不能用只有强度的远场图反推相干像面。s/p 分量先还原为笛卡尔矢量，再旋转到物镜坐标并平行输运到准直光瞳。光瞳振幅除以 `sqrt(w_o)`，保留上式功率。

倍率 `M` 将光线映射为像方方向：

```text
q = (-u/M, -v/M, w_i),  w_i = sqrt(1-(u²+v²)/M²)
E(x,y,z) = 1/(lambda*M) integral
    sqrt(T) R_image R_object F_far / sqrt(w_o*w_i)
    * exp{ i*(2*pi/lambda)*[-(u*x+v*y)/M + z*(w_i-1)] } du dv
H*Z0: 对每个角谱方向以 q cross E 构造，再独立传播
I_image = I_ref * Re[E_x*(H*Z0)_y^* - E_y*(H*Z0)_x^*]
```

分析器作用在准直光瞳局部 x/y 轴上；像差是光瞳相位。省略的公共轴向载波不影响强度或同参考下的干涉。20 mm 物镜、20 倍对应 400 mm 管镜；20 倍时 NA=0.96 的像方 NA 为 0.048，而不是 0.96。倍率加大时像宽增大、强度密度下降，总功率不应变化。

以上采用矢量角谱/Richards–Wolf 思路；方向余弦积分测度和矢量高 NA 要求可参见开源 [PSF Generator 理论文档](https://psf-generator.readthedocs.io/src/theory/basic.html)。本报告中的归一化与误差数字是对本项目代码的独立推导和实际复算，不是该文档对本项目的背书。

## 独立传播与颗粒参考

使用三个不同层次，避免只重复生产 FFT：

- 独立 Rodrigues 旋转与直接平面波积分，不调用生产光瞳几何、透镜旋转或 FFT。
- 均匀 x 偏振圆光瞳，用解析方位角积分得到 Bessel/Hankel 参考；径向 Gauss 积分提高阶数检查参考自身收敛。
- 平整晶圆单球，用 NIST Bobbert–Vlieger Jones 振幅加连续极坐标积分，与已有 FEM/FFT 像面比较；另外提高 BV 多极阶数。

| 项目 | 本次结果 | 验证范围 |
| --- | --- | --- |
| 独立直接积分对 Bessel/Hankel 复 E/H | 最大相对 L2 `9.52e-15` | NA=0.8，M=1/20，含 M=1 离焦 0.4 μm |
| 已有局部粗糙面三路像面，独立复 E/H 和强度复算 | 最大相对 L2 `4.91e-16` | 每路中心 ±12 μm，729 个像点；同一离散复光瞳 |
| 整光斑混合模型三路强度复算 | 最大相对 L2 `4.29e-16` | 每路中心 ±12 μm，81 个像点；精简存档不含像面复 E/H |
| BV 连续积分加密 | 最大强度相对 L2 `3.04e-15` | 径向/方位 24×96 对 48×192 |
| BV 多极阶数提高 | 最大强度相对 L2 `7.72e-11` | lmax=14 对 18 |
| 离散光瞳到完整 FFT 像面的功率闭合 | 最大相对误差 `6.20e-14` | 本次偶极基准矩阵 |

局部粗糙面源为 `build/near_pml_order/pipeline_outer`。整光斑源为 `build/model_campaign_20261010/large_hybrid`：266 nm、100 nm 颗粒、1 W、晶圆面 2×2 μm 高斯照明、NA=0.5、三路倾斜 +45/-45/0 度。上述源来自已有存档，本次没有重新求解 FEM 网格；复远场的直接积分与 BV 参考计算均重新运行。

BV/FEM 比较条件：355 nm、100 nm 球、颗粒 n=1.59/k=0、基底 n=5.57/k=0.387、间隙 10 nm、入射 75 度、NA=0.5、20 倍、128 点光瞳，ROI 为像面中心 ±12 μm，不做峰值归一化、位置配准或拟合相位。

| 入射偏振 | +45° 路强度相对 L2 | -45° 路强度相对 L2 | 0° 路强度相对 L2 |
| --- | --- | --- | --- |
| p | 0.1623% | 0.1691% | 0.4916% |
| s | 0.5259% | 0.3307% | 0.4238% |

这些数包含已有 FEM 散射误差和 FFT 光瞳离散误差。不能将其当成所有颗粒尺寸、材料、粗糙面或波长的统一误差上界。

## 光瞳采样与高 NA

新增 [`imaging_convergence_audit.py`](imaging_convergence_audit.py) 用自由空间复偶极辐射和均匀圆光瞳检查连续积分误差。偶极矩 `(0.3,0.2i,1)` 是数学基准，不是本次实测 PSL；266 nm、20 倍，报告绝对幅度下的两条中心线剖面误差，不冒充整个二维像面的误差上界。

高 NA 倾斜孔径的参考积分解析处理半球裁切边界；边界处径向替换平滑根号变化，避免通过不连续网格掩膜定义“真值”。两套 64×128、96×192 积分的本次最大剖面差异 `1.71e-14`。

| NA / 倾斜 | 光瞳点数 | 收集功率误差 | 中心剖面强度相对 L2 |
| --- | --- | --- | --- |
| 0.5 / 0° | 32 | 1.497% | 1.931% |
| 0.5 / 0° | 64 | 0.537% | 0.679% |
| 0.5 / 0° | 96 | 0.136% | 0.195% |
| 0.5 / 0° | 128 | 0.293% | 0.373% |
| 0.96 / 0° | 32 | 3.075% | 3.271% |
| 0.96 / 0° | 64 | 1.198% | 1.162% |
| 0.96 / 0° | 96 | 0.324% | 0.462% |
| 0.96 / 0° | 128 | 0.649% | 0.660% |
| 0.96 / 45° | 32 | 3.495% | 4.000% |
| 0.96 / 45° | 64 | 0.968% | 0.773% |
| 0.96 / 45° | 96 | 0.036% | 0.264% |
| 0.96 / 45° | 128 | 0.402% | 0.305% |

另用 z 向偶极的严格收集公式验证绝对测度：`C=2*pi*(1-c)^2*(c+2)/3`，`c=sqrt(1-NA²)`。这不是“按孔径几何面积算能量”，而是积分偶极的 `sin²(theta)` 辐射分布。

额外将倍率降到 M=1，验证真正高像方 NA=0.96 的矢量传播；128 点均匀光瞳在焦点的强度剖面误差为 0.393%，离焦 0.4 μm 时为 0.714%。这些测试超过默认 20 倍的像方角度，但仍受光瞳离散限制。

**使用建议：定量结果至少比较 64、96、128 点的总功率、中心强度和所需 ROI 能量；不要只看 Parseval。FFT padding 增大只细化像面采样，不提高散射角/光瞳积分精度。** 极低 NA、大离焦、近地平线窄峰、很大颗粒仍需另做参数相关收敛。

## 焦点参考与粗糙面含义

SCATMECH 成像入口直接使用原生 JonesDSC，保留其振幅参考；FEM 的焦点平面经过晶圆原点。对于本次无镀膜基准，高度 `h=半径+间隙`，将球心原点的出射振幅换到晶圆原点时需乘 `exp(-i*k*h*k_z_global)`。这保持每个方向的 DSC 不变，却改变像面相位、位置和离焦。本次 BV 独立参考已做该转换，没有对 UI 结果拟合。NIST 的 [BV 定义](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm) 和 [原生源码](https://pages.nist.gov/SCATMECH/code/bobvlieg1.cpp) 可用于检查其原始场约定。

因此“远场强度一样，像面不一样”不必然意味着传播错误。必须先统一复场相位原点、焦点、偏振、孔径和透过率。

只有 PSD/平均 Mueller BRDF 不能唯一确定相干像面。SCATMECH 平均粗糙背景的直接相干成像被现有入口拒绝；FEM/整光斑 SPM 是先对具体随机表面算复场，再传播，统计时对强度求平均。这个传播流程可校验，但本次没有重新验证所有随机表面的物理统计与高阶耦合误差。

三路各自成像；重叠角域在真实仪器中不能不经分光就被多个接收器同时完整接收。三路功率之和不是去重后的联合收集率。接纳镜面方向的提示也不代表已模拟真实挡光器或镜面泄漏。

## 可复算命令

在仓库根目录、NIST2 Python 下执行：

```powershell
& D:/anaconda/envs/NIST2/python.exe -m fullwave_prototype.imaging_convergence_audit `
  --output build/imaging_validation_20261012/convergence `
  --saved build/near_pml_order/pipeline_outer

& D:/anaconda/envs/NIST2/python.exe -m fullwave_prototype.imaging_audit `
  --source build/model_campaign_20261010/large_hybrid `
  --flat-source build/imaging_independent/flat128 `
  --reference-python D:/anaconda/envs/NIST2/python.exe `
  --output build/imaging_validation_20261012/independent

& D:/anaconda/envs/NIST2/python.exe -m pytest `
  tests/test_imaging.py tests/test_fullwave_imaging.py tests/test_direct_imaging.py `
  tests/test_theory_validation.py tests/test_imaging_entry.py tests/test_aperture_integrals.py `
  tests/test_imaging_convergence_audit.py tests/test_detector.py `
  -q --disable-warnings -p no:cacheprovider `
  --junitxml=build/imaging_validation_20261012/focused_tests.xml
```

`--saved` 可以省略，此时纯解析/连续积分基准不依赖历史 FEM 文件。第二条命令和存档复算需要本机现有结果文件，克隆代码并不会自动带上被忽略的 `build` 目录。

本次成像、探测器与物理参考相关测试：**59 passed、2 skipped**；跳过的是 NIST2 中未安装 `miepython` 的两个独立 Mie 校验，未将其计为通过。40 个警告来自现有 Qt/SIP 废弃 API。数值结果为 `build/imaging_validation_20261012/{independent,convergence}/report.json`，收敛图为 `convergence/convergence.png`，测试记录为 `focused_tests.xml`。
