# WaferScope

晶圆缺陷暗场光学仿真：从散射近场、远场到三路独立成像和相机灰度信噪比。

**换机入口（2026-10-29）：** [迁移手册](docs/PROJECT_MIGRATION_CN.md)；
[轻量源码与修改版 S4 下载](https://github.com/wufansaber/waferscope/releases/tag/waferscope-migration-lite-2026-10-29)。
有限光栅＋单颗粒的连续远场/像面程序已保存，首次原生求解因内存预检未启动，仍待新电脑运行。

项目原名为 DFI Simulation，更早使用 APTK / FarField Simulation 名称。现有计算接口、结果目录和旧版本文件保持兼容；更名不改变物理模型或数值精度。

## 运行

在包含 PyQt5 和 pySCATMECH 的 NIST2 环境中，从项目根目录运行：

```powershell
conda run -n NIST2 python far_field_standalone\main.py
```

`main.py` 是保留的综合平台入口，两个入口均显示 WaferScope。FEM 使用独立的 `build/ngsolve_env` 环境。新增 `S4 Grating (periodic)` 页面使用 `build/rcwa_env` 调用私有 S4 stripes 原生程序，输出三平面近场、离散远场、理想周期相干像面，以及像元积分后的模拟面阵相机原始 DN 和时间统计 SNR。

## 模型与边界

- SCATMECH：球形颗粒、膜层、PSD 粗糙散射、孔径收集、成像与相机。
- FEM：实验性 Maxwell 近场、远场及成像；有限光栅仍需网格、PML 和提取方法收敛。
- RCWA：S4 周期场景为主，grcwa 用于辅助核对；保留原始复场和收敛检查。周期颗粒不等于孤立缺陷，理想中继像面不等于实际相机数据，尚未完成孤立极限认证。

## 文档

- [操作手册](far_field_standalone/OPERATION_MANUAL_CN.md)
- [迁移与继续开发](docs/PROJECT_MIGRATION_CN.md)
- [独立界面与分文件打包](far_field_standalone/README.md)
- [FEM 模型与验证](fullwave_prototype/README.md)
- [S4 周期光栅页面与原始数据](docs/GRATING_SCENE_GUIDE_CN.md)
- [光栅场景本轮开发与实算验收](docs/GRATING_SCENE_DEVELOPMENT_CN.md)
- [颗粒电流的连续角分布与原始复场导出](docs/CONTINUOUS_ANGLES_CN.md)
- [S4 有限光栅＋颗粒的连续远场与像面方案](docs/S4_FINITE_MIXED_CN.md)

仓库：<https://github.com/wufansaber/waferscope>。当前源码的构建目标为 `WaferScope.exe` 和 `fem/WaferScope_FEM.exe`；历史 Releases 中的 EXE 不会因仓库更名自动更新。
