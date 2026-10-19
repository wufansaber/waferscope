# Windows 分文件打包验证

日期：2026-09-28。平台：Windows x64。

## 发布结构

界面使用 NIST2 / Python 3.8 / PyInstaller 6.9.0，求解器使用独立 Python 3.11 /
PyInstaller 6.16.0。两套 Python DLL、NumPy 和依赖分别存放，不混用环境。
NGSolve 6.2.2607、Netgen、OCCT、OpenBLAS 及运行库随 FEM EXE 附带。
打包时补齐了冻结环境所需的 libffi 和 setuptools 的 backports 模块。

主程序自动选取 `fem/FarField_FEM.exe`；源码模式仍支持原独立 Python。
结果写入用户目录，启动目录不必是源码目录。数据不会写入 `_internal`。
FEM 保留原有时间/内存保护、取消和子进程清理机制。所有物理默认值未因打包而改变。

## 历史版本 fullwave-preview 测试

1. NIST2 测试集：90 项，79 通过、11 项因原生依赖条件跳过。
2. 发布版 GUI 完整计算：清除 PATH 中的 Anaconda，清空 PYTHONPATH/PYTHONHOME，
   从临时目录启动 EXE，经 QProcess 启动独立 FEM EXE，再启动受保护工作进程。
3. 用 preview、RMS=0、32 点光瞳进行配对计算，成功输出近场、远场和三路像面，
   GUI 自动加载结果并截图。耗时 58.4 秒，计量峰值 3.352 GiB，无残留子进程。
   三路功率约 159.3 / 159.0 / 41.76 pW；这是启动回归算例，不是高精度粗糙面验证值。
4. 发布版 SCATMECH `--headless-test` 成功计算并导出图像。

原始本机记录：`build/package_checks/gui/package_smoke.json`、
`build/package_checks/gui/results/*/resources.json`、
`build/package_checks/gui/packaged_gui.png`、`build/package_checks/scatmech`。

测试是在当前 Windows 机器上隔离路径运行，不等同于在全新 Windows 虚拟机上安装测试。
程序未经商业代码签名，不能保证所有杀毒软件和所有 Windows 版本的兼容性。
安装包中的验证示例采用此前完整粗糙面计算，不以低精度启动测试代替物理校验。

## 发布边界

源码和构建说明提交 Git；分文件 ZIP 与 SHA-256 作为 GitHub Release 附件上传。
`build`、`dist`、本地 `.codex` 配置、Python 环境和完整仿真缓存不提交。
FEM 仍标记 experimental，发布标记为预发布，不暗示所有模型误差均已解决。

## 当前版本 2026.09.28-imaging-upgrade

- NIST2：109 项测试，98 通过、11 项因原生依赖条件跳过；NGSolve 环境补跑的 49 项全部通过。
- 新增像差与独立直接积分对照、有限束界面连续性、缓存失效、精简存档和相机回归检查。
- 冻结 GUI 在不含 Anaconda 的 PATH、清空 PYTHONPATH/PYTHONHOME、临时工作目录下通过完整测试。
- GUI 完整计算使用 preview、RMS=0、有限束腰 2 um、3 阶角谱、32 点光瞳、4 个 FEM 线程、
  15 nm 彗差和相机。耗时 33.984 秒、计量峰值 3.815 GiB；这是启动回归，不是高精度基准。
- 同一 GUI 紧接着选择已保存的源、移除彗差重新成像，三路全部命中缓存，耗时 3.922 秒、
  计量峰值 0.176 GiB。相机图、compact 存档和结果自动加载均通过，无遗留受管子进程。
- 独立冻结求解器重成像也通过；与未共享传播函数的直接积分比较，三路强度 L2 差异最高 3.38e-16。
- SCATMECH 冻结版计算、二维/三维图导出通过。已人工查看 FEM 界面、相机图及 SCATMECH 图像。

原始记录：`build/release_upgrade/packaged/package_smoke.json`、其 results 下的 resources.json、
`build/release_upgrade/native_tests.json`、`build/release_upgrade/frozen_reimage_resources.json`。
发布包 verification 文件夹附有关键 JSON 汇总，docs 包含本次理论、数值边界和使用说明。
build_manifest 记录源提交及工作区状态；ZIP 的 SHA-256 随包提供。仍非全新虚拟机或实验仪器标定测试。
