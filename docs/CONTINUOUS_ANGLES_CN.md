# 连续角分布：计算、查看与原始数据

2026-10-28。本功能已加入源码界面的 **S4 Grating (periodic)** 页，工具栏按钮为 **Continuous angles**，结果页签同名。现有已打包的 EXE 尚未重建，运行源码才能使用新增界面。

## 物理含义

采用已有 `contrast_current.project_ray` 的互易辐射表达式，从单个超胞中颗粒内部的复电场，在裸周期光栅/晶圆背景上逐方向计算连续复振幅 F。没有把离散衍射级插值成连续强度。DSC 为 `sum(abs(F)**2)`，单位 um²/sr；F 的单位为 um。

结果是单胞颗粒对比电流的连续辐射。原始电流仍受周期复制颗粒影响，尚不能认证为孤立颗粒的定量解。它不是整个样品的总远场，也不是样品与裸晶圆强度之差。有限光斑、有限光栅端点和相机灰度没有在本功能中建模。

每次计算先验证原始数据哈希、颗粒/裸参考的几何与求解器元数据，并通过同基组复场重建门槛：加权复场相对差不大于 1e-7，功率相对差绝对值不大于 2e-7。这个门槛核对表达式、系数和参考面的一致性，不替代 Fourier 截断、分层或孤立性收敛。

## 源码界面操作

1. 用项目原有启动方式运行源码界面，进入 **S4 Grating (periodic)**。
2. 查看已生成的数据时，点击 **Load**，选择包含 `continuous_report.json` 的文件夹，例如 `output/continuous_angles_20261028/hemisphere_32x96`。
3. **Continuous angles** 页签显示 theta/phi 二维图和四条邻近采样方位角截线。勾选 **3D direction surface** 切换三维方向表面，颜色表示绝对 DSC，半径表示方向而非强度。工具栏 **Log color scale** 控制色标和截线纵轴。
4. **Export raw bundle** 导出完整 ZIP，含原始颗粒和裸参考模式、报告、复场、CSV 和图片；导出结果可移动后重新加载。结果目录必须为空，计算不会覆盖已有结果。
5. 重新计算时点击 **Continuous angles** 按钮，选择匹配的颗粒电流与裸参考目录，设置 polar/azimuthal 采样数和反向背景截断。运行时使用左侧 **Worker Python / EXE**，开发目录默认 `build/rcwa_env/Scripts/python.exe`。该环境必须具有项目使用的 grcwa；原有 FEM 环境没有 grcwa 时不能计算此功能。支持内存、时间保护和 Cancel。

S4 普通场景的 `modes.npz` 未导出颗粒内部电流，不能直接充当本功能的源。已有可复用的配对输入为：

- `output/s4_round7_20261027/particle_lx0667_mx5_current`
- `output/s4_round7_20261027/reference_lx0667_mx5`

这组受控基准采用 266 nm、75° P 入射、100 nm 球、10 nm 间隙、500/250/100 nm 光栅，源基组 Mx=5/My=40。材料沿用旧值 Si=5.57+0.387i、球=1.59，不能解释为标定后的接触 PSL 算例。加载连续结果时界面显示其自身的波长、球径和间隙；左侧 S4 场景参数不会隐式替换已存的电流。

## 命令行复现

在仓库根目录执行，输出目录须为空：

```powershell
& build/rcwa_env/Scripts/python.exe -m rcwa_prototype.angular_distribution --source output/s4_round7_20261027/particle_lx0667_mx5_current --reference output/s4_round7_20261027/reference_lx0667_mx5 --output output/continuous_angles_new --polar-samples 32 --azimuthal-samples 96 --reference-halfwidth 40 --timeout 3600
```

新电流可由 `tools.rcwa_trial --particle --particle-method analytic --current-z-order 8` 导出，裸参考应使用同一几何、基组和求解器设置。详见 [第七轮记录](S4_ROUND7_CN.md)。本功能不自动发起新的高成本颗粒全波求解。

## 坐标、积分与文件

theta 从 +Z 法线计，0° 到 90°；phi 从 +X 向 +Y 计，0° 到 360°。在 `mu=cos(theta)` 上采用 Gauss-Legendre 求积，方位角使用均匀中点采样，严格避开精确掠射方向。角度网格的固角权重和为 2π，积分为 `sum(weight*DSC)`，单位 um²。网格是连续函数的数值采样，不含“每个像素一个衍射级”的含义；没有峰值归一化、左右取均值或全局相位拟合。

- `angular_distribution.npz`：theta、phi、方向、固角权重、原始三分量复 F、DSC。
- `angular_distribution.csv`：每行一个方向，含坐标、权重、复场实虚部与 DSC，表头带单位。
- `continuous_angles.png` / `continuous_3d.png`：二维角度图、截线、三维方向图。
- `modal_gate.npz`：同基组重建与直接差场。
- `source/`、`reference/`：源与裸参考的原始 `modes.npz` 和 `report.json`。
- `continuous_report.json`：完成标记、积分、源记录、采样设置、适用范围及所有文件哈希，最后写入。
- 输出目录旁的 `.resources.json`：守护进程的时间和内存记录。

没有完成标记的取消/失败结果不会被界面作为有效结果加载。角度加密对照只检查当前源与反向背景截断下的角积分敏感性，不认证整个物理模型。

## 本轮实算与验证

原始结果保存在 `output/continuous_angles_20261028/`。同一源与反向背景半截断 40 下，上半球截面从 24×72 网格的 0.00208818071937 um² 变为 32×96 网格的 0.00209509369861 um²，相对变化 +0.3310527%。这不是点态误差界，也不是所有后续加密误差的上限。同基组复场重建相对差为 3.10667e-12，功率相对差为 7.06324e-13。

两次守护计算均正常完成，无遗留子进程；耗时分别约 197 秒和 352 秒，进程树峰值统计内存分别约 0.121 和 0.124 GiB。原始资源记录、比较 JSON、数值测试 XML、界面测试 XML 与实际结果加载截图均在上述目录。

数值相关 57 项测试通过，界面/打包路径/亮度显示相关 21 项测试通过。覆盖固角和为 2π、半球方向一阶矩、解析偶极子 4π/3 积分、复相位保留、失败门槛阻止输出、损坏数据拒绝、二维/三维切换和 ZIP 移动后的完整加载。相关回归共 78 项，不代表整个应用的全量物理验证。
