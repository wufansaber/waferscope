# WaferScope 轻量换机迁移手册

更新：2026-10-29。平台：Windows x64。按用户选择，**本次不发布历史仿真原始数据，通用依赖使用官方下载链接**。本轮仅迁移代码和方案，没有启动新的电磁求解。

## 1. 获取代码

[轻量迁移发布](https://github.com/wufansaber/waferscope/releases/tag/waferscope-migration-lite-2026-10-29)

- 源码标签：`waferscope-migration-lite-2026-10-29`
- 分支：`codex/migration-s4-lite-20261029`
- `WaferScope_Source_20261029.zip`：源码、测试、手册、参数示例、环境版本快照；无历史近场/远场/像面原始数据。
- `WaferScope_S4_Modified_20261029.zip`：项目修改版 S4 原生程序、对应源码、许可证和补丁；无 OpenBLAS DLL、编译器和 Python 环境。
- 两包附 `.sha256` 和逐文件 manifest。文档中的 PDF 手册是文档，不是历史原始场数据。

推荐克隆固定标签并限制历史下载：

```powershell
git clone --depth 1 --branch waferscope-migration-lite-2026-10-29 https://github.com/wufansaber/waferscope.git D:\code\waferscope
Set-Location D:\code\waferscope
git switch -c codex/continue-finite-s4
```

也可解压源码 ZIP，进入包含 `main.py`、`tools`、`s4_prototype` 的目录。将修改版 S4 包内部根文件夹的**内容合并到项目根目录**，得到：

```text
build/s4_runtime/bin_stripes/S4.exe
build/s4_runtime/bin_stripes/build_manifest.json
build/s4_runtime/licenses/
build/s4_runtime/downloads/                 # 仅 S4/Lua 原始源码归档
build/s4_runtime/victorliu-S4-5b89848/       # 修改后的对应源码
```

解压前可用 `Get-FileHash 下载路径 -Algorithm SHA256`，与同名 `.sha256` 比较。此包不是新版 WaferScope GUI EXE。

## 2. 依赖下载地址

| 依赖 | 官方地址 / 安装入口 | 用途 |
| --- | --- | --- |
| Git for Windows | [官方下载](https://git-scm.com/install/windows) | 克隆项目及安装固定 Git 依赖 |
| Python 3.11 x64 | [Windows 版本列表](https://www.python.org/downloads/windows/) | S4 调度、grcwa 辅助核对和 FEM 独立环境；选择 3.11 系列 |
| NumPy、SciPy、Matplotlib、psutil、pytest、autograd | [PyPI](https://pypi.org/)；按 `rcwa_prototype/requirements.txt` 安装固定版本 | 通用 Python 计算依赖 |
| grcwa | [固定源码提交](https://github.com/weiliangjinca/grcwa/tree/c52f10a138f4b2d5c897b7670967813c36af8d9b) | 辅助核对及源诊断，实际 S4 求解仍调用原生程序 |
| OpenBLAS 0.3.33 Windows x64 | [ngsolve-openblas 固定发行文件](https://pypi.org/project/ngsolve-openblas/0.3.33/#files) | S4 所需 `openblas.dll`，下节脚本仅下载并取出 DLL |
| Conda（可选） | [Miniconda 官方说明](https://www.anaconda.com/docs/getting-started/miniconda/main) | 恢复 GUI/NIST 环境 |
| pySCATMECH 0.1.9（可选） | [PyPI](https://pypi.org/project/pySCATMECH/0.1.9/) | 原有 NIST/SCATMECH 页面 |
| NGSolve 6.2.2607（可选） | [PyPI](https://pypi.org/project/ngsolve/6.2.2607/) | FEM 独立求解，不是此次 S4 的必装依赖 |
| w64devkit 2.10.1（仅重编 S4） | [固定 x64 编译器下载](https://github.com/skeeto/w64devkit/releases/download/v2.10.1/w64devkit-x64-2.10.1.7z.exe) | 默认直接运行现有 S4 时不需要 |
| S4 原始源码（仅重编） | [固定提交归档](https://api.github.com/repos/victorliu/S4/tarball/5b89848cfd15719ccaa04c853fe883ecf8d6ac15) | 对应源码已随修改版 S4 包提供 |
| Lua 5.2.4（仅重编） | [官方源码](https://www.lua.org/ftp/lua-5.2.4.tar.gz) | 对应归档已随修改版 S4 包提供 |

**项目修改版 S4 无法用同名 pip 包或未经修改的上游 S4 替代。** 这里增加了显式矩形 Fourier 基和 stripe 接口；有限沟槽使用完整二维几何，不启用一维 stripe shortcut。

## 3. 恢复 S4 计算环境

安装 Git、64 位 Python 3.11 后，在项目根目录执行：

```powershell
py -3.11 -m venv build\rcwa_env
& .\build\rcwa_env\Scripts\python.exe -m pip install -r rcwa_prototype\requirements.txt
& .\build\rcwa_env\Scripts\python.exe -m tools.restore_s4_dependencies --download
& .\build\rcwa_env\Scripts\python.exe -m tools.s4_finite_mixed --help
```

先解压修改版 S4 包，再执行 `restore_s4_dependencies`。该脚本从官方 PyPI 下载 `ngsolve-openblas==0.3.33` 的 Windows x64 wheel，**不安装 NGSolve**，只将其中的 `openblas.dll` 恢复到 S4 目录，并检查它是否与原生构建清单的 SHA-256 一致。

也可手动从上表下载 Windows x64 wheel，再运行：

```powershell
& .\build\rcwa_env\Scripts\python.exe -m tools.restore_s4_dependencies --wheel C:\Downloads\ngsolve_openblas-0.3.33-py3-none-win_amd64.whl
```

文件名以实际下载为准；不使用 ARM/Linux wheel。`S4.exe` 与 `openblas.dll` 必须同目录。

```text
S4.exe SHA256       1da0f38040f8f8b168269280fd7c013c536a84a252bcb7f2fb08b898f2271741
openblas.dll SHA256 ea5706aaa7e88ea76f350faf4181c5e3393817a1ce1cd800501bf3ee21f205e8
```

不复制原机 venv/Anaconda 目录。原机三套环境版本保存在 `docs/migration_20261029/ENVIRONMENT_SNAPSHOTS.json`，新机解释器路径可以不同。若没有 Python Launcher，用新机的 Python 3.11 解释器替换 `py -3.11`。

## 4. 待继续的计算

目前保存了三条不同范围的程序：

1. `grating_scene/`、`grating_ui/`：S4 无限周期光栅＋周期球，原生离散级次、周期像面与相机模拟。
2. `rcwa_prototype/angular_distribution.py`：grcwa 的颗粒电流连续角诊断，不能标作 S4 全混合场。
3. **`tools/s4_finite_mixed.py`**：此次目标——有限沟槽片＋单颗粒、边缘散射、连续远场及三路相干像面。**第一次原生求解因旧机内存预检未通过而未启动，没有该组完成结果。**

新入口仍采用周期超胞内的有限片近似；四周平坦缓冲区和实际沟槽端部进入模型，但胞间耦合与孤立极限尚未验证。程序由完整沟槽和颗粒复对比源重建连续辐射，不对离散级次强度做插值。FEM 保留独立对照；MATLAB 独立程序在 `matlab/independent_farfield/abgle.m`，用户自行验证，不与 UI 连接。

内存充足后，从空输出目录开始：

```powershell
& .\build\rcwa_env\Scripts\python.exe -c "import psutil; print('available GiB=',psutil.virtual_memory().available/2**30)"
& .\build\rcwa_env\Scripts\python.exe -m tools.s4_finite_mixed --output output/s4_finite_mixed_20261029/preview_221 --cell-x 1.5 --cell-periods-y 4 --mx 6 --my 8 --slices 8 --length 0.8 --trenches 2 --gap-nm 0 --source-xy 10 --source-z 8 --polar-samples 24 --azimuthal-samples 72 --na 0.6 --magnification 20 --irradiance 1e-6 --pupil-samples 64 --extraction-tolerance 0.02 --timeout 1800
```

参数：266 nm、75° P、500 nm 周期、250 nm 槽宽、100 nm 槽深、两条长 0.8 um 的有限槽、100 nm 球接触中央脊顶、1.5 × 2 um 超胞、8 球切片、221 模态。三路 NA=0.6，倾角 +45/-45/0 度，20 倍理想相干成像。这是预览设置，未收敛。Si=5.57+0.387i、球=1.59 是沿用对照常数，未完成实验材料标定。

默认保留 3 GiB 系统内存余量、进程树预算 8 GiB、1800 秒超时；规划或监控失败就退出，不通过降低保护强行求解。监控不是操作系统硬内存上限。

完成条件：`.resources.json` 为 `completed`，`extraction_diagnostic.json` 中 `passed=true`，且最后写入 `finite_mixed_report.json`。原生差场/体源重建相对复场门槛为 2%，失败先检查源积分和 Fourier 表示；局部门槛通过不等于孤立物理模型验收。保留原始复场、参考/样品、差场强度、真实强度差与干涉项。

远场是相对平坦晶圆的散射场，不加入无限平面波镜面反射。详细设置、输出和后续收敛项目见 [S4_FINITE_MIXED_CN.md](S4_FINITE_MIXED_CN.md)。

## 5. 按需恢复 GUI / FEM

仅运行独立 S4 程序无需 GUI 或 NGSolve。需要原界面时：

```powershell
conda env create -n NIST2 -f environment.yml
conda run -n NIST2 python far_field_standalone\main.py
```

GUI 采用独立 Python 3.8 / NumPy 1.24 / Matplotlib 3.7 / PyQt5 / pySCATMECH 环境。若 pySCATMECH 没有匹配 Windows wheel，安装可能需要 Microsoft C++ 工具。不要把 NumPy 2.x 或 FEM 原生库混入旧 NIST2。重新设置 UI 求解器路径、IDE 解释器与 Qt 首选项。因为本次不附历史数据，旧结果图库需自行复制原结果或重新计算；JSON 参数示例仍在 `examples/`。

FEM 独立环境：

```powershell
py -3.11 -m venv build\ngsolve_env
& .\build\ngsolve_env\Scripts\python.exe -m pip install --only-binary=:all: -r fullwave_prototype\requirements-validation.txt
& .\build\ngsolve_env\Scripts\python.exe -m pip install pytest
```

如需重编 S4：从上表下载编译器并改名 `w64devkit-api.7z.exe`；其 SHA-256 必须为 `add1edf984ce4d473cfeb7c52e0ab38db73a855bba6c7f564d5e16648eb2896a`。在单独的 `build/s4_rebuild` 准备下载归档和解压后的 `w64devkit/bin`，避免覆盖已保存运行库：

```powershell
& .\build\rcwa_env\Scripts\python.exe -m tools.build_s4_runtime --root build/s4_rebuild --variant stripes --eigensolver builtin --blas build/s4_runtime/bin_stripes/openblas.dll
```

运行新编程序时显式传 `--executable build/s4_rebuild/bin_stripes/S4.exe`；构建脚本会核对所有输入归档。

## 6. 交接和旧机保留

本次轻量提交不新增历史原始数据；源码树只保留文档 PDF。压缩包由已提交 Git 树生成，含版本快照和文件校验。`MIGRATION_MANIFEST.json` 的 `omitted_local_outputs` 仅列出旧机未随包提供文件的路径、大小及 SHA-256，不包含这些原始场文件。仓库更早的历史记录不属于本次新上传内容，使用上面的 `--depth 1` 获取本版即可。

如仍需旧结果，清理旧机前自行备份完整 `output/` 和需要的 `build/` 到移动硬盘；新机验收后再处理旧副本。FEM 检查点需同时保留 `checkpoint.json` 和 `field.ngs`，并使用匹配 NGSolve。GitHub 登录、IDE 私人设置及软件安装需在新机恢复。

给新电脑的后续任务：

> 继续 WaferScope S4 主线。先读 docs/PROJECT_MIGRATION_CN.md、docs/S4_FINITE_MIXED_CN.md。恢复修改版 S4 和 build/rcwa_env，从官方 PyPI 下载匹配 OpenBLAS，再检查内存，运行 tools.s4_finite_mixed 的 preview_221。用户选择有限尺寸光栅＋单颗粒，需要边缘散射、连续远场与三路相干像面。第一次原生求解尚未启动，没有该组完成结果；保持 3 GiB 系统余量，不执行 MATLAB 验证。保留完整复场，继续检查 Fourier、球切片、源积分、角度、孔径与超胞收敛。此次不带历史原始数据。