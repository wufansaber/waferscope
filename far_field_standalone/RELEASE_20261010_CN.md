# FarField Simulation 2026-10-10 粗糙面与相干照明更新

Windows x64 分文件版。完整解压 ZIP 后运行 `FarField_Simulation.exe`，保留
`_internal` 和 `fem` 文件夹；无需另装 Python、Anaconda 或 NGSolve。

## 更新

- 整光斑 PSD 粗糙面 SPM 与颗粒局部 FEM 使用共享母表面、照明和复场相位，减去重复的局部粗糙散射。
- 增加同母表面的范围、采样、边缘渐隐校验，低 RMS FEM/SPM 绝对功率比较与能量闭合检查。
- 支持固定 128 张 haze 表面统计扩展及三路同时置信区间。
- 支持晶圆面高斯、光束截面高斯、均匀矩形及实测 NPZ 光斑；明确缺失相位的假设。
- 输出近场 XY/XZ/YZ、远场及独立三路像面，更新模型和校验文档。

## 验证与限制

已记录算例：低 RMS 半球 FEM/SPM 差 0.419%，能量闭合误差约 0.088%；
固定 128 张 haze SPM 表面的三路同时置信区间半宽为 3.24%、3.80%、3.18%。
上述结果不是所有输入参数下的误差保证。

FEM/混合模型仍为实验功能：90.7 μm 范围仅用于收敛审计，未自动替换生产成像窗口；
整光斑高阶颗粒/粗糙面耦合仍需独立验证。请阅读随包模型文档和 `verification` 摘要。

ZIP 与 SHA-256 校验文件在本 Release 提供。程序未作代码签名。

## 操作手册

[新版中文操作手册（MD）](https://github.com/wufansaber/aptk/releases/download/v2026.10.10-roughness/OPERATION_MANUAL_CN.md) 包含参数设置、10 nm 示例、照明定义、粗糙面统计、三路成像、结果保存和常见问题。
手册是独立补充附件，原 ZIP 和 EXE 未修改。
