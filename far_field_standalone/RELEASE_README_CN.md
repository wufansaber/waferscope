# FarField Simulation Windows 分文件版

## 启动

完整解压 ZIP 后，双击 `FarField_Simulation.exe`。不要从压缩包内部运行，
也不要单独移动 EXE；请保留 `_internal` 和 `fem` 文件夹。
本版附带两个独立运行时，不需要安装 Anaconda、Python 或 SCATMECH/NGSolve。

- `SCATMECH`：颗粒远场、三路收集、独立像面和探测器。
- `Full-wave FEM (experimental)`：近场、远场、三路独立像面和相干分解。
- `fem/FarField_FEM.exe`：主界面自动调用的求解器，通常不用手动启动。

FEM 结果默认保存在 `%LOCALAPPDATA%\FarField_Simulation\results`。
默认使用三阶 balanced、128 点光瞳；粗糙面算例内存需求较大，建议至少 16 GB RAM，
并保留足够空闲内存。界面保留内存、时间限制和取消按钮；触发保护不是有效计算结果。

## 示例与说明

在 FEM 页选择 `Load result`，载入 `examples/validated_355nm_100nm`，可以直接查看
已计算的近场、远场和三路像面。示例是 355 nm、100 nm 球、三路 NA=0.5、
倾角 +45/-45/0 度、20 倍；它不是 266 nm 参数，也不是全三阶默认配置的结果。
完整计算参数显示在 `Result parameters` 页。示例仅包含图像和报告，不包含大体积复场缓存。

- [用户指南](docs/far_field_standalone/USER_GUIDE_CN.md)
- [理论模型](docs/fullwave_prototype/THEORETICAL_MODEL_CN.md)
- [验证与限制](docs/fullwave_prototype/LOW_MEMORY_REFINEMENT_CN.md)

FEM 仍是实验模型，不是完备解析解：有限粗糙区域、几何尖角、材料色散、网格与吸收边界
仍存在误差限制。文档中指向源码 `build` 的历史原始数据链接只在开发工作区可用。
未经代码签名的程序可能受到 Windows SmartScreen 提示；请核对来源与发布的 SHA-256。

第三方依赖的许可和版本见 `THIRD_PARTY_NOTICES.md` 与 `licenses` 文件夹。
