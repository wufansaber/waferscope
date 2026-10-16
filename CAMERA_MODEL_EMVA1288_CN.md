# 相机模型与 EMVA 1288 校验

更新：2026-10-12。SCATMECH 和 FEM 共用同一套相机模型；三路分别成像、分别生成时间噪声和固定像素图样。新增的是可配置单色 CMOS 线性读出模型，**不是某个相机型号的实测标定，也不是 EMVA 认证**。已有 EXE 未重新打包。

## 噪声定义

采用 [EMVA 1288 Release 4.0 Linear 官方标准](https://www.emva.org/wp-content/uploads/EMVA1288Linear_4.0Release.pdf) 的时间噪声、空间不均匀和增益约定；对应核心公式为第 2.4、2.5、8.3、8.5 节。标准的系统增益 K 为 DN/e⁻，程序输入 `gain_e_dn` 为 e⁻/DN，二者互为倒数。

| 参数 | 模型含义 | 单位/填写方式 |
| --- | --- | --- |
| quantum_efficiency | 当前波长有效 QE；需由相机数据或标定给出 | 0～1；示例为 0.6，不代表 266 nm 实测值 |
| read_noise_e | 最低曝光的电子域时间暗噪声，不含 ADC 量化 | e⁻ RMS；曝光相关暗电流散粒噪声另外加入 |
| dark_current_e_s | 平均暗电流，温度和读出模式在一次计算中固定 | e⁻/pixel/s |
| dsnu_e | 固定暗偏置的像素间 RMS；同一 ROI/传感器多帧不重新生成 | e⁻；不是 e⁻/s |
| prnu_fraction | 固定光响应相对 RMS | `0.01` = 1% |
| dark_current_std_e_s | 额外暗电流率不均匀 DCNU | e⁻/s；属于扩展，不混入 DSNU 定义 |
| gain_e_dn | 电子数到灰度码的换算 | e⁻/DN；默认 0.5，即 K=2 DN/e⁻ |
| offset_dn / bit_depth | 黑电平 / ADC 位深 | DN / 8～16 bit |
| full_well_e | 光电荷积累阶段的满阱限制 | e⁻；可能先满阱，也可能先 ADC 饱和 |

PRNU 使用均值为 1 的固定对数正态响应图，再将有效 QE 限制在 0～1；DSNU 使用固定零均值高斯偏置图。分布是合成建模假设，不是标准强制指定的分布。QE 接近 1、PRNU 较大时裁切会改变实现后的响应 RMS，应检查导出的 `qe_map`。

## 从光场到原始图像

```text
pixel_power_W = 对像面 W/um^2 积分像素敏感面积
mean_photons = pixel_power_W * exposure_s / (h*c/lambda)
S = mean_photons * effective_QE
B = exposure_s * (dark_current + extra_background)
charge ~ Poisson(S+B)
electrons = min(charge, full_well) + Normal(0, read_noise_e) + fixed_DSNU_offset
raw_DN = clip(round(electrons/gain_e_dn + offset_dn), 0, 2^bit_depth-1)
```

像素不是像面中心点采样，而是面积积分。可选 `fill_factor` 用居中的正方形敏感区描述，默认 1；若有效 QE 已包含填充因子，保持 1，避免重复扣光。图像积分采用离散像场的分片常数近似，还需要光瞳/像面采样收敛。

量化由实际 round/clip 实现，不额外加一次随机量化噪声。理论时间方差用电子域 `S+B+sigma_read²+gain²/12` 预测；空间加时间预测再加 `DSNU²+(PRNU*S)²`，开启 DCNU 时加其曝光相关项。量化不足或饱和区域用 `digital_snr_valid` 标识无效。低于噪声分辨条件时，1/12 DN² 近似不应被视为精确的暗图方差。

`digital_ideal_snr` 是给定像素映射、已知背景条件下的时间 SNR；`emva_total_snr` 使用平均 QE，预测未校正空间不均匀时的总 SNR。ROI 预测假设像素空间差异彼此独立，不涵盖相关条纹或背景估计误差。输出有 `noise_definition` 和标准链接说明适用范围。

光学 haze 应先进入像面光场，再经过同一相机链路。`background_e_s` 只是额外均匀电子率，不能替代 PSD/FEM haze。当前纯整光斑 SPM 平均像面入口仍限制为光学结果，未改变其相机开关限制。

## 验证结果

本次验证为合成平场与解析积分；未复现标准的实测照明均匀性校正、高通滤波和行/列分量测量流程，因此参数恢复结果不是认证结果。

50 档平场，每档两帧，通过差分去除固定图样；另对暗/亮场各平均 64 帧并扣除残留时间方差。反推结果：

| 项目 | 输入 | 恢复值 | 相对误差 |
| --- | --- | --- | --- |
| K | 2 DN/e⁻ | 1.999949 DN/e⁻ | 0.0026% |
| QE | 0.6 | 0.600032 | 0.0053% |
| 时间暗噪声 | 2 e⁻ RMS | 1.991579 e⁻ | 0.4211% |
| DSNU | 3 e⁻ RMS | 3.005438 e⁻ | 0.1813% |
| PRNU | 2% | 1.994731% | 0.2634% |

另用偏心高斯像面的解析 erf 像素积分检查“W → 光子 → 电子”转换，相对 L2 误差为 `2.659e-6`（0.000266%）。测试还覆盖固定图样跨帧稳定、不同通道独立、零信号仍有黑电平、ADC 先于满阱饱和、满阱后加读出噪声和 16 位 TIFF 无损往返。全量测试 **275 passed、16 skipped**；NIST2 缺少可选原生求解/参考依赖的测试保持跳过，并非全部物理案例通过。记录为 `build/camera_validation_20261012/all_tests.xml`。

## 三路颗粒图像示例

生成 `build/camera_validation_20261012/particle_10nm/`：10 nm PSL、266 nm、p 偏振、75° 入射、NA=0.5、倾斜 +45/−45/0°、20 倍、128 点光瞳。曝光 0.01 ms，6.5 μm 像素，16×16 ROI，16 bit，0.5 e⁻/DN，黑电平 100 DN，PRNU=1%，DSNU=0.5 e⁻。

光学照明按局部均匀平面波处理，参考强度来自 1 W、晶圆面 2×2 μm 矩形面积及入射投影；没有在这个 SCATMECH 示例里求解有限 2×2 μm 光束。输入材料 n/k 和 QE 都是示例值，不是该波长的标定数据。该图不含 haze，不证明真实系统能检测 10 nm 颗粒。

两侧光学 ROI 功率均约 `8.553e-10 W`，中路约 `2.475e-10 W`。固定 PRNU 图不同，使相机电子数有小幅差异，不代表入射方向不对称。示例均未满阱或 ADC 饱和。

输出 `detector_pixels.png` 为可读的对比图；`camera_ch1.tiff`～`camera_ch3.tiff` 为 16 位容器内的原始 ADC 灰度码，没有每路归一化。`coherent_imaging.npz` 保存像素光功率、电子数、传感器图样、灰度码和有效性掩膜。10 nm 颗粒的相机图像受衍射和像素积分决定，不是直径 10 nm 的几何圆。

FEM 额外验证复用了已有 355 nm、100 nm 颗粒及局部粗糙面存档，接入相机并输出灰度 TIFF；1 ms 曝光下两侧各 4 个像素满阱，中路未满阱。未重新求解 FEM，也未把饱和后的图像解释为线性 SNR。

## 操作与复算

SCATMECH 页启用 `Detectors`；FEM 页启用 `Camera response`。调整增益、ADC、PRNU/DSNU 和曝光后重新成像。参数随 JSON 保存；旧结果缺少的新字段按默认值加载。源码运行入口不变，2026.10.16 分文件发布版已接入相机与灰度 SNR；更早的 EXE 不包含这些更新。

```powershell
& D:/anaconda/envs/NIST2/python.exe -m camera_validation `
  --output build/camera_validation_20261012
```

本机有 FEM 像面存档时，可以额外加 `--saved-image build/camera_validation_20261012/fem_rough`，仅复算相机电子学并输出 `fem_camera/`，不再算散射。

不包括滚动快门运动、温度依赖标定、EMCCD 倍增、像素串扰、blooming、坏点、相关行列噪声、颜色滤片、非线性响应或实际平场/暗场文件导入。使用具体硬件时应填写其当前波长、温度、曝光/读出模式下的标定参数，不能把示例参数称为实测相机。
