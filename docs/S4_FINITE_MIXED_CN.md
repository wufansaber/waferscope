# S4 有限光栅与单颗粒：保留的程序和计算方案

日期：2026-10-29。状态：**程序已保存，本组原生 S4 仿真尚未运行完成，没有新的连续远场或像面结果。用户选择稍后再算，本轮不继续启动求解。**

换机安装、原生 S4 运行包恢复和启动步骤见 [项目迁移手册](PROJECT_MIGRATION_CN.md)。

## 选择依据与已有记录

继续采用 S4 主线，承接 [求解器路线](GRATING_SOLVER_ROADMAP_CN.md)、[周期混合场景](GRATING_SCENE_GUIDE_CN.md) 和 [混合场景对照摘要](migration_20261029/PERIODIC_VALIDATION_20261028.md)。FEM 继续作为独立方法对照，不在此次选择中替换 S4。轻量迁移版不附历史原始算例。

必须区分三个版本：

1. 原有 S4 UI：无限长直线光栅加周期复制颗粒，直接输出离散衍射级与周期像面。
2. 2026-10-28 连续角度诊断：使用 grcwa 内部颗粒电流，不是新的 S4 连续全混合场结果。
3. 此次保留的新程序：真正调用原生 S4，构造有实际端部的有限光栅沟槽、中心颗粒及四周平坦缓冲区，提取沟槽和颗粒的完整对比电流，然后计算连续角分布与相干像面。

**新程序仍是 S4 周期超胞近似。** 一个超胞中只放一个有限片和一个颗粒，但整个超胞仍会被周期复制。有限沟槽的端部和侧壁确实进入几何与源积分；是否已逼近一个孤立有限片，必须另做超胞尺寸收敛。它不是已实现 AFMM/PML 的非周期解。

## 保留的物理与光学设置

- 266 nm，75° P 偏振入射，从 -KY 方向入射。
- 两条有限沟槽，周期 500 nm、沟槽宽度 250 nm、深度 100 nm、长度 0.8 um；外侧恢复为平坦晶圆。
- 中央脊顶为 100 nm 球，默认间隙 0 nm，采用 8 层中点圆盘表示。
- 初始超胞 1.5 × 2 um，Mx=6、My=8，共 221 个 Fourier 模态。这只是准备的预览设置，未经精度验收。
- Si=5.57+0.387i、颗粒=1.59：沿用受控旧常数，不称为 266 nm 标定材料。
- 平面波照明，光强归一化为 1e-6 W/um²。没有把参考面积当作有限光斑。
- 三路独立空气物镜，NA=0.6，绕 Y 倾角 +45°、-45°、0°，20 倍，理想无像差相干成像。

连续远场是相对于平坦半空间的散射场，包含有限沟槽和颗粒的共同辐射。平坦晶圆的无限平面波镜面反射不加入连续 DSC；上述三路为暗场收集。既不将颗粒差场当作混合总场，也不将离散衍射级插值成连续分布。

## 计算链

```text
S4：相同有限片，分别求无颗粒参考和含颗粒样品
  -> 在有限沟槽体积与分层球体积内采样原生复 E
  -> 对比极化源 (epsilon - epsilon_flat_reference) * E * dV
  -> 平坦半空间格林函数重建每个观察方向的复散射振幅
  -> 保存 reference、sample、delta 及强度差/干涉项
  -> 在每路物镜孔径上重新计算复振幅
  -> 偏振输运、孔径衍射与像面复场叠加
  -> 保存有颗粒、无颗粒、强度差、差场强度及干涉项
```

沟槽采用硅基底中的有限 Air 矩形；有限端部破坏沿槽不变性，因此新入口禁止使用旧的一维解析 stripe shortcut。原有周期场景的默认参数和调用方式不变。

像面保留 `E_sample = E_reference + deltaE` 的完整干涉，不能将 `|deltaE|²` 代替 `I_sample - I_reference`。本次只做光学像面，不加入相机 DN、读出噪声或 UI 接口。

## 稍后运行

在仓库根目录执行；输出目录必须为空：

```powershell
& build/rcwa_env/Scripts/python.exe -m tools.s4_finite_mixed --output output/s4_finite_mixed_20261029/preview_221 --cell-x 1.5 --cell-periods-y 4 --mx 6 --my 8 --slices 8 --length 0.8 --trenches 2 --gap-nm 0 --source-xy 10 --source-z 8 --polar-samples 24 --azimuthal-samples 72 --na 0.6 --magnification 20 --irradiance 1e-6 --pupil-samples 64 --extraction-tolerance 0.02 --timeout 1800
```

资源保护仍保留 3 GiB 系统余量、默认 8 GiB 进程树上限；先做矩阵内存规划。首次尝试在启动守护检查时因可用内存不足而退出，未启动原生 S4。未降低保护阈值，未关闭用户应用。

## 预期输出与完成判定

- `reference/`、`sample/`：原生 S4 Lua 输入、标准输出、近场 E/H、离散模态和构建/哈希信息。
- `source_quadrature.npz`、`*_sources.npz`：实际源坐标、体积权重与复对比极化源。
- `modal_extraction_diagnostic.npz`、`extraction_diagnostic.json`：用连续源重建相同传播级的颗粒差场，与原生 S4 差场比较，不拟合幅相。默认加权复场相对差门槛 2%；不通过时停止像面发布，留下诊断。
- `continuous_far_field.npz/csv/png`：上半球逐方向复场及 reference/sample/delta/强度差/干涉项，DSC 单位 um²/sr。
- `image_fields.npz`、`mixed_images.png`、`imaging_report.json`：三路复像面场、W/um² 强度和收光功率。
- `geometry.png`：有限沟槽和中心颗粒示意。
- `finite_mixed_report.json`：最后写入的完成标记，含完整设置、误差诊断和文件哈希。没有这个文件就不能声称该组结果已完成。
- 输出目录旁 `.admission.json` 与 `.resources.json`：内存预检及实际运行资源记录。

已有相关测试 57 项通过，覆盖有限源体积/矩、缓冲区和 stripe shortcut 拒绝，以及原有 S4 接口回归。这不是新物理模型的全波验收，也不是新的像面实算结果。添加了运行时门槛和预检后，依用户“稍后再算”的选择，没有继续发起求解。

待实际运行后再安排：源体积分辨率、完整二维 Fourier 截断、球切片、超胞复制耦合、连续角积分与孔径采样的收敛。只有局部一致性门槛通过，不能将结果升级为孤立有限光栅的定量认证。
