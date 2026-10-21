# WaferScope

晶圆缺陷暗场光学仿真：从散射近场、远场到三路独立成像和相机灰度信噪比。

项目原名为 DFI Simulation，更早使用 APTK / FarField Simulation 名称。现有计算接口、结果目录和旧版本文件保持兼容；更名不改变物理模型或数值精度。

## 运行

在包含 PyQt5 和 pySCATMECH 的 NIST2 环境中，从项目根目录运行：

```powershell
conda run -n NIST2 python far_field_standalone\main.py
```

`main.py` 是保留的综合平台入口，两个入口均显示 WaferScope。FEM 使用独立的 `build/ngsolve_env` 环境；RCWA 试验使用 `build/rcwa_env`，尚未接入 UI。

## 模型与边界

- SCATMECH：球形颗粒、膜层、PSD 粗糙散射、孔径收集、成像与相机。
- FEM：实验性 Maxwell 近场、远场及成像；有限光栅仍需网格、PML 和提取方法收敛。
- RCWA：实验性 grcwa 三维超胞，用于独立对照；周期颗粒不等于孤立缺陷，当前结果尚未收敛。

## 文档

- [操作手册](far_field_standalone/OPERATION_MANUAL_CN.md)
- [迁移与继续开发](docs/PROJECT_MIGRATION_CN.md)
- [独立界面与分文件打包](far_field_standalone/README.md)
- [FEM 模型与验证](fullwave_prototype/README.md)

仓库：<https://github.com/wufansaber/waferscope>。当前源码的构建目标为 `WaferScope.exe` 和 `fem/WaferScope_FEM.exe`；历史 Releases 中的 EXE 不会因仓库更名自动更新。
