# 相机端缺陷、haze 与噪声检测 SNR

更新：2026-10-14。适用于 FEM 和相干混合模型的配对像面；本次没有重新打包 EXE。

## UI 操作

1. 在 FEM 页设置颗粒、照明、晶圆、粗糙面和三路孔径。非零粗糙 RMS 才产生物理 haze。
2. 勾选 `Camera response`，配置像素、ROI、曝光、QE 和 EMVA 1288 相机参数。
3. `Camera reference frames` 是同一无缺陷粗糙面的独立噪声参考帧数，默认 1，可平均 1～128 帧。
4. `Camera haze mismatch fraction` 是扣除参考后的额外 haze RMS/均值，默认 0 表示同一静态表面；`0.1` 表示 10%。这是 SNR 预算的用户假设，不是从一张 FEM 表面推断出的 PSD 统计精度，也不在原始相机帧上人为叠加新的光场。
5. 计算完整流程，或对保存源重新成像。结果的 `Camera` 页新增四个子页：

| 子页 | 内容 |
| --- | --- |
| Raw camera | 含缺陷的原始 DN、无缺陷 haze 参考平均 DN、两者差值 |
| Defect / haze / noise | 缺陷电子数对比、haze 电子数、检测噪声 RMS；三路共用对应绝对色标 |
| Pixel detection SNR | 带符号像素 SNR；饱和、ADC 截断或量化不足的像素用灰色屏蔽 |
| Detection metrics | 三路缺陷/haze/噪声预算、ROI 求和 SNR、已知模板 SNR及有效状态 |

启用新配对相机结果时，Camera 子页默认显示指标表。旧结果没有配对信息时仅显示旧相机图，不使用残留的上次 SNR。

## 信号定义

使用同一粗糙面、同一照明和焦点参考的两套复光场：有颗粒的 `with` 与无颗粒的 `without`。先传播，再按相机敏感面积积分并转换为电子数。

```text
mu_with[j] = 含颗粒及粗糙背景的平均光电子数
mu_haze[j] = 同一粗糙面无颗粒的平均光电子数
delta[j] = mu_with[j] - mu_haze[j]
```

`delta` 可以为负，保留颗粒和 haze 的相干干涉、耦合影响。它既不是总亮度，也不是简单取差分复场的自能量，更不能把“独立平滑晶圆颗粒功率”和“独立 haze 功率”相加代替总场。

## 噪声与 SNR

两套原始相机读出使用同一传感器固定 QE/DSNU/DCNU 图，时间散粒/读出噪声独立；参考帧逐帧完成满阱及 ADC 操作后再平均，不先平均光子就绕过饱和。

令 N 为参考帧数，b 是平均暗电流及额外均匀电子背景，sigma_r 是电子域时间读出 RMS，g 是 e⁻/DN，f 是附加 haze 不匹配比例。未饱和、量化充分、像素独立时：

```text
V[j] = mu_with[j] + mu_haze[j]/N
       + (b_with[j] + sigma_r^2 + g^2/12)
       + (b_haze[j] + sigma_r^2 + g^2/12)/N
       + (f * mu_haze[j])^2

signed_pixel_SNR[j] = delta[j] / sqrt(V[j])
ROI_sum_SNR = abs(sum(delta)) / sqrt(sum(V))
known_template_SNR = sqrt(sum(delta[j]^2 / V[j]))
```

含颗粒帧的光子方差 `mu_with` 已含 haze 光子，不另加一次完整 haze 方差。参考 haze 光子噪声是额外的 `mu_haze/N`。原始噪声按 [EMVA 1288 相机模型](CAMERA_MODEL_EMVA1288_CN.md) 生成。

同一相机同曝光下，固定 DSNU/DCNU 偏置在参考扣除中抵消；不能把它们重新生成成两套独立的时间噪声。PRNU 保留在像素实际 QE 和缺陷响应中。这里的检测 SNR 条件于固定像素图和已知缺陷模板，不等于旧的单幅未校正空间不均匀总 SNR。真实暗场/平场校正残差、相关行列噪声仍需要实测输入或进一步模型。

ROI 求和会使正负干涉对比相互抵消；已知模板 SNR 使用图样信息，在独立高斯噪声近似下更能描述一个已知模板的可辨识度。它不是未知缺陷的实测检出率；低光子计数下也不能直接把该数当作精确的假警概率。

全 ROI 有任意无效像素时，主 ROI/模板 SNR 显示 `Invalid`，JSON 保存 `null`，不报出虚假的完整 ROI 检测能力。内部保存的 `valid_pixel_template_snr` 仅描述剩余有效像素，不作为主指标。

**一个随机表面只能给出该表面的 haze 图与条件 SNR，不能给出独立表面 speckle 的统计误报率。** 默认 f=0 依赖“同一静态表面参考”假设；参考来自不同位置/不同表面时需评估残差或做多表面统计，不能把 f=0 当成实际 wafer 扫描背景已被完全消除。

## 本次验证

全量测试：**293 passed、16 skipped**，记录 `build/detection_validation_20261014/tests.xml`。跳过依赖缺失的原生求解/独立参考案例，不表示所有物理参数下都已证明正确。

新增校验包含解析电子噪声预算、纯 haze 时缺陷 SNR 必须为零、破坏性干涉保留负对比、固定缺陷下 haze 增大使 SNR 下降、参考平均噪声按 1/N 减小、DSNU 固定偏置抵消、Monte Carlo 差分帧方差、关闭孔径、饱和屏蔽、FEM 端到端接线及 UI 旧结果清理/缺失文件原子校验。

实际存档复算：355 nm、100 nm 颗粒、局部粗糙面 RMS=2 nm、NA=0.5、三路 +45/-45/0 度、20 倍、16×16 ROI、6.5 μm 像素、输入参考照度 1e-6 W/μm²。未重新求解 FEM；像面沿用已验证存档。

| 通道 | 曝光 0.01 ms：缺陷对比 e⁻ | haze e⁻ | ROI 求和 SNR | 已知模板 SNR |
| --- | --- | --- | --- | --- |
| +45° | 1818.71 | 3.318 | 32.615 | 40.432 |
| -45° | 2179.02 | 27.639 | 36.831 | 44.193 |
| 0° | 576.82 | 14.373 | 13.298 | 21.099 |

参考平均 4 帧，附加 haze 不匹配为 0。曝光改为 1 ms 后，两侧各 4 个满阱像素；两侧的主 ROI/模板 SNR 正确变为无效。中路在这套输入下仍有效。这些数不是当前默认 1 W 高斯照明或 10 nm 颗粒的检测结果。

复算相机端而不重新求解散射：

```powershell
& D:/anaconda/envs/NIST2/python.exe -m fullwave_prototype.detection_replay `
  --source build/camera_validation_20261012/fem_rough `
  --output build/detection_validation_20261014/unsaturated `
  --exposure-ms 0.01 --reference-frames 4
```

输出的 `imaging_report.json` 和 `image_fields.npz` 可直接由 UI 加载；NPZ 增加相机缺陷对比、haze、方差预算、参考平均 DN、差值 DN、像素 SNR 和有效性掩膜。`detector_components.png`、`detector_snr.png` 为新增诊断图；原始有缺陷帧继续导出三路灰度 TIFF。该脚本输出目录必须与原始存档目录分离，不修改原始光学结果。

SCATMECH 纯颗粒相干成像页仍不将平均 Mueller haze BRDF 当成确定复场；完整颗粒+haze 相机检测功能在 FEM/混合模型页使用。2026.10.16 分文件发布版已接入该功能及原始灰度帧 SNR，操作见随包手册第 15 节；更早的 EXE 不包含这些更新。
