# FarField Simulation Windows 分文件版

## 启动

完整解压 ZIP 后，双击 `FarField_Simulation.exe`。不要从压缩包内部运行，
也不要单独移动 EXE；请保留 `_internal` 和 `fem` 文件夹。
本版附带两个独立运行时，不需要安装 Anaconda、Python 或 SCATMECH/NGSolve。

- `SCATMECH`：颗粒远场、三路收集、独立像面和探测器。
- `Full-wave FEM (experimental)`：近场、远场、三路独立像面和相干分解。
- `fem/FarField_FEM.exe`：主界面自动调用的求解器，通常不用手动启动。

FEM 结果默认保存在 `%LOCALAPPDATA%\FarField_Simulation\results`。
默认使用三阶 balanced、128 点光瞳；粗糙面算例内存需求较大，本轮验收使用 32 GB RAM，
并保留足够空闲内存。界面保留内存、时间限制和取消按钮；触发保护不是有效计算结果。

## 示例与说明

在 FEM 页选择 `Load result`，载入 `examples/validated_355nm_100nm`，可以直接查看
已计算的近场、远场和三路像面。示例是 355 nm、100 nm 球、三路 NA=0.5、
倾角 +45/-45/0 度、20 倍；它不是 266 nm 参数，也不是全三阶默认配置的结果。
完整计算参数显示在 `Result parameters` 页。示例仅包含图像和报告，不包含大体积复场缓存。

- [用户指南](docs/far_field_standalone/USER_GUIDE_CN.md)
- [完整操作手册 PDF](docs/FarField_Simulation_中文操作手册.pdf)
- [软件功能与案例 PDF](docs/FarField_Simulation_功能与案例手册.pdf)
- [理论模型](docs/fullwave_prototype/THEORETICAL_MODEL_CN.md)
- [验证与限制](docs/fullwave_prototype/LOW_MEMORY_REFINEMENT_CN.md)

FEM 仍是实验模型，不是完备解析解：有限粗糙区域、几何尖角、材料色散、网格与吸收边界
仍存在误差限制。文档中指向源码 `build` 的历史原始数据链接只在开发工作区可用。
未经代码签名的程序可能受到 Windows SmartScreen 提示；请核对来源与发布的 SHA-256。

第三方依赖的许可和版本见 `THIRD_PARTY_NOTICES.md` 与 `licenses` 文件夹。

## 本次更新

2026-10-19 重新打包，加入横线、竖线与十字有限缺陷、计算域自动扩展，以及
适配缺陷尺寸的近场与像面显示范围。图案仍需做网格和尖角收敛验证，不能将预览图视为准确形貌预测。
附带 21 页中文操作手册和 10 页功能与案例手册。案例混合远场已使用同一保存源，
按 401×401 方向网格重算，消除旧 5° 采样的扇区色块，保留单张粗糙面的真实相干散斑。
该加密是案例图更新，不代表主程序所有远场默认采样均已改为 401×401，也不是全模型收敛证明。

2026-10-16 增加 Gemi 面阵近似预设、原始灰度帧采样和独立灰度 ROI SNR 指标、
相机单独重算、低起伏 266 nm 预设，以及保存的 NIST/对称性报告查看。
相机不考虑 TDI 累积；转换增益暂估，不能将模拟 SNR 当作实拍标定。
透过率支持 0.0001，不会再因小数精度被界面舍入为零。
操作步骤见 [操作手册第15节](docs/far_field_standalone/OPERATION_MANUAL_CN.md)。
顶部文档图标可选择 `verification/camera_20261016` 查看保存的同参数 NIST 对照。
该低 RMS 例子的统计 haze 功率差约 0.2%，不代表全模型误差已经收敛到 1%。

2026-10-10 的历史版本增加整光斑粗糙面 SPM、颗粒/局部 FEM 的共享母表面和复场混合、
固定 128 张随机表面的统计扩展，以及表面范围、网格、渐隐、低 RMS 和功率校验。
支持晶圆面高斯、光束截面高斯、均匀矩形和实测 NPZ 四种照明定义。
实测光斑缺少相位时必须明确选择平相位假设或拒绝；示范 NPZ 不代表真实测量。

2026-10-10 校验的三路同时置信区间半宽为 3.24%、3.80%、3.18%（仅 haze SPM）。
同一小粗糙面的低 RMS FEM/SPM 半球功率差为 0.419%，能量闭合误差约 0.088%。
这些数值只适用于已记录的校验参数，不是任意模型误差保证。
90.7 μm 大范围仅用于收敛审计，未自动替换生产成像窗口；整光斑高阶耦合仍待验证。
详见 [本次粗糙面校验](docs/fullwave_prototype/ROUGHNESS_CAMPAIGN_CN.md)。
随包 `verification` 保存本次摘要，`historical_20260928` 是旧版验证，不能视为本次 EXE 验收。

保留独立成像校验、FEM 线程设置、有限角谱光束、Zernike 像差与相机响应。
`Re-image source` 复用已有散射源，只应用成像和相机设置；修改颗粒或材料需重新 `Calculate`。
默认 compact 存档并自动缓存复出瞳；需要像面复 E/H 时选择 full。
本机验收使用 32 GB RAM，默认保留 2 个 FEM 线程；提高线程数会增加内存需求。
模型说明和未完成研究项见 [本次模型扩展](docs/fullwave_prototype/MODEL_UPGRADE_CN.md)。
独立校验和线程实测见 [验收报告](docs/fullwave_prototype/IMAGING_THREAD_AUDIT_CN.md)。
