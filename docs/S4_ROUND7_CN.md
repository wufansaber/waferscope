# 第七轮：固定超胞的 Fourier 截断检查

日期：2026-10-27。承接 [第六轮](S4_ROUND6_CN.md)。本轮是 grcwa 诊断计算，S4 主线、GUI 和 FEM 后端不变。

## 1. 问题和固定计划

第六轮 X 胞长从 2/3 增至 5/6 um 后，连续边路/中路收光变化仅 -0.5163%/-1.6405%，但复场仍变化 3.7410%/4.5653%。当时固定的是物理最大 Fourier 频率，改变胞长仍会改变频率采样密度，因此不能将差异直接认定为周期颗粒之间的物理耦合。

本轮固定 Lx=2/3 um、Ly=1 um，把 X 截断从 Mx=4 提至 5，保持 My=40。分别为 729 和 891 模态，X 最大空间频率从 6 提至 7.5 cycles/um；不扩大几何域，不改变照明、材料或颗粒。

固定输入仍是 266 nm、75 度 P 平面波、光源来自 -KY；周期/线宽/深度=500/250/100 nm；100 nm 球、球底与脊顶间隙 10 nm；受控旧材料值 Si=5.57+0.387i、球=1.59+0i。沿用 8 球切片、每片 8 点 z 积分，三路独立 NA=0.6、倾角 +45/-45/0 度。

间隙和材料沿用仅为控制变量，不代表经 266 nm 标定的接触 PSL 颗粒。未经收敛的周期电流结果不用于定量成像或相机 SNR。

## 2. 本轮状态

固定计划见 [plan.json](../output/s4_round7_20261027/plan.json)。新的 891 模态颗粒解、裸参考及两级连续孔径积分均已正常完成。

颗粒解耗时 3023.33 秒，峰值进程树统计内存 2.66785 GiB，监测最低系统可用内存 3.44619 GiB；保护阈值仍为进程树统计 8 GiB、系统余量 3 GiB。正常退出，无遗留子进程。保护是采样监控，不是操作系统硬内存上限。

30 个共同传播级的同基组互易复场重建相对差 `3.10667e-12`，差场功率相对差 `7.06324e-13`。这是系数、参考面与互易表达式的一致性，不是基组误差界。

固定胞长后，内部复 E 的 L2 变化 **2.09509%**，相干对比极化矩变化 **0.360391%**。径向/方位角 20x80、32x128 两级实空间积分相对解析圆盘结果的 L2 差异不超过 `6.04e-13`，仅核对存储截断场的积分方法。记录见 [独立源场积分核对](../output/s4_round7_20261027/source_overlap_audit.json)。

比较必须使用相同连续角度、固角权重和反向背景截断；保留绝对复幅，不拟合整体相位。先通过各自的同基组复场重建，再比较内部源场 L2、相干对比极化矩和三路连续差场截面。

一次单轴加密即使小于局部门槛，也不能认证 Y 截断、球切片、z 积分、超胞孤立性或有限光栅端点。

沿用 8 层中点球切片的体积比解析球体积高 **0.78125%**；横向圆盘 Fourier 系数解析不等于连续球体的垂直几何也无误差。本轮不改变切片，避免将该误差混入 X 截断对照。

## 3. 比较检查和原始数据约定

新增 [固定超胞比较器](../rcwa_prototype/cutoff_check.py) 与 [报告工具](../tools/rcwa_cutoff_report.py)：

- 仅允许一个 Fourier 半截断增加，拒绝胞长、材料、间隙、分层或另一轴同时变化。
- 验证颗粒/裸参考 NPZ 哈希及各自输入，同基组复辐射重建必须先通过。
- 核对求解器版本、无损传播分支、全局坐标及相位参考；不混用旧分支或重新拟合相位。
- 读取原始复 F、方向和固角权重重新计算每路截面，拒绝角度、权重或保存积分不一致。
- 内部场使用物理坐标中的圆盘解析重叠；不同基组不能按数组索引直接相减。

源场重建约定为 `E(r,z_j)=sum_G E_G(z_j)*exp(i*kxy_G dot r)`；`source_kxy_per_um` 的单位为 rad/um，坐标以脊顶中心的全局 z=0 为参考，时间约定 `exp(-i*omega*t)`。z 采样与权重见 `source_z_um`、`source_weight_um`；内部积分半径按 `source_radius_um`，并非用球外整个超胞体积代替。

连续辐射保存 `far_electric_um`、`directions`、`solid_angle_weights` 及各点 DSC；收光截面为 `sum(w*|F|^2)`，单位 um^2。它不是实际相机灰度或 `|E_sample|^2-|E_reference|^2`；三路分开保存，不按峰值归一化或强制左右对称。

比较器测试还覆盖纯相位变化、原始参考文件完整性、错误版本/参考面拒绝，以及报告不覆盖旧结果。基准输入、旧源场与旧辐射分别位于第五轮，裸参考位于第四轮；复现时应一并保留这些依赖目录。

相关回归 **180 项通过**，见 [测试报告](../output/s4_round7_20261027/tests.xml)。不是整个应用全量回归或所有物理场景的验证。

## 4. 固定超胞的实际结果

两组均使用相同的 24x72 连续孔径积分、反向光栅半截断 40。单位 um^2；三路由原始复 F 分别积分，没有峰值归一化、整体相位拟合或左右取均值。

| 通道 | Mx=4, My=40 截面 | Mx=5, My=40 截面 | 截面变化 | 加权复场变化 |
|---|---:|---:|---:|---:|
| CH1 | 4.71746109e-4 | 4.74636576e-4 | +0.612717% | 0.558377% |
| CH2 | 4.71746109e-4 | 4.74636576e-4 | +0.612717% | 0.558377% |
| CH3 | 4.36232546e-4 | 4.37015664e-4 | +0.179519% | 0.343011% |

复场变化定义仍是 `sqrt(sum(w*|F_refined-F_baseline|^2)/sum(w*|F_baseline|^2))`，不去掉整体相位。裸光栅参考的受激子块控制符合 X 不变性：反射率均为 `0.2760783016950492`，对应加权复反射场变化为本次计算中的 `0.0`；这来自同一个精确背景子块，不是两套独立物理方法达到零误差。

新源电流的角度积分从 16x48 加密到 24x72 后，边路截面变化 -0.215711%、中路 +0.151234%；总积分耗时 997.125 秒，峰值进程树统计内存 0.085743 GiB，正常退出、无遗留子进程。相邻加密变化不作为绝对误差界；固定同一角度网格比较两个源也不能消除所有积分误差。

**判断：这一次固定超胞 X 向加密通过功率 <5%、复场 <2% 的局部门槛，但不是全模型收敛。** 在已测的 2/3 um 超胞中，本次 X 截断增量引起的远场变化比第六轮相邻胞长变化小；不能由一次加密推断所有后续截断误差均小，也不能把两种百分比相减推算周期复制耦合。

内部场 L2 仍变化约 2.10%，比集成极化矩和孔径辐射变化大；远场积分趋稳不等于近场每一点也达到同样精度。

![固定超胞的 X 截断对照](../output/s4_round7_20261027/x_cutoff_comparison/cutoff_sensitivity.png)

下一步优先固定同一超胞及 Mx，单独增加 My（例如 40 -> 48）检查 Y 截断。之后再安排 Y 向颗粒间距、球切片/z 积分以及更高联合截断，均需先通过实际内存预检。保留 S4 匹配复核和 FEM 独立方法验证；本轮不替换 GUI，不发布定量孤立颗粒成像或相机 SNR 结论。

## 5. 文件与复现

- [机器可读比较和数据哈希](../output/s4_round7_20261027/x_cutoff_comparison/comparison.json) 与 [完整表格](../output/s4_round7_20261027/x_cutoff_comparison/COMPARISON.md)
- [新增原始模态与内部复场](../output/s4_round7_20261027/particle_lx0667_mx5_current/modes.npz)
- [新增连续方向、权重及原始复 F](../output/s4_round7_20261027/radiation_lx0667_mx5/cone_24_72_ref40.npz)
- [互易重建检查](../output/s4_round7_20261027/radiation_lx0667_mx5/gate.json) 与 [独立实空间积分核对](../output/s4_round7_20261027/source_overlap_audit.json)
- [颗粒资源记录](../output/s4_round7_20261027/particle_lx0667_mx5_current.resources.json)、[裸参考资源记录](../output/s4_round7_20261027/reference_lx0667_mx5.resources.json)、[连续积分资源记录](../output/s4_round7_20261027/radiation_lx0667_mx5.resources.json)

旧基准辐射没有内嵌 NPZ 哈希，报告明确记为 `legacy_no_embedded_radiation_hash`，同时记录本轮读取时的实际哈希；不将它补称为当时已有完整的哈希链。新源场、参考和辐射文件的记录哈希均已验证。

实际调用参数如下；输出目录必须为空，重跑应更换目录并相应调整引用，不能覆盖本轮记录。

```powershell
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_trial --output output/s4_round7_20261027/particle_lx0667_mx5_current --harmonics 5 40 --lx 0.6666666666666666 --particle --particle-method analytic --grating-method analytic_li --current-z-order 8 --timeout 3600
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_trial --output output/s4_round7_20261027/reference_lx0667_mx5 --harmonics 5 40 --lx 0.6666666666666666 --grating-method analytic_li --timeout 300
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_current_radiation --source output/s4_round7_20261027/particle_lx0667_mx5_current --reference output/s4_round7_20261027/reference_lx0667_mx5 --output output/s4_round7_20261027/radiation_lx0667_mx5 --level 16 48 40 --level 24 72 40 --timeout 1800
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_cutoff_report --baseline output/s4_round5_20261023/radiation_lx0667 --refined output/s4_round7_20261027/radiation_lx0667_mx5 --output output/s4_round7_20261027/x_cutoff_comparison
```
