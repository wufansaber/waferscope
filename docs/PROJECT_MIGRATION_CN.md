# DFI Simulation 项目迁移与继续开发

更新日期：2026-10-21。目标平台：Windows x64。

本说明用于把当前暗场光学仿真项目迁移到新电脑并继续开发，不是软件精度认证。迁移包包含源码、操作文档、环境版本快照和精选验证数据；不包含可直接搬走使用的 Python 环境，也不是新打包的 EXE。

## 1. 当前状态

| 模块 | 状态与入口 |
| --- | --- |
| SCATMECH 远场、收集、成像和相机 | 原有主程序路径保留，入口 `far_field_standalone/main.py` |
| FEM 近场、远场、成像 | `fullwave_prototype/` 与 `fullwave_ui/`，实验模型 |
| 有限三维光栅 | `tools/finite_grating_run.py`，有限长度沟槽；当前此入口不含颗粒 |
| 周期光栅独立验证 | `tools/grating_reference.py`、`tools/periodic_grating_audit.py`、`tools/vector_grating_audit.py` |
| RCWA 三维颗粒超胞 | `rcwa_prototype/`、`tools/rcwa_trial.py`，独立 grcwa 后端，尚未接入 UI |

S4 和 AFMM 尚未完成实算。当前机器未检测到 S4 编译工具链或可用 WSL，因此本轮先使用 grcwa；不要把 grcwa 结果描述为 S4 或 AFMM 结果。

### 必须保留的物理约定

- 激光来自负 KY 一侧，入射传播方向在晶圆面上的投影沿正 KY；入射角相对表面法线。
- 三路为独立光路，常用倾角为 +45、-45、0 度；不要用总和替代各路独立成像。
- 新光栅验证常用 266 nm、75 度、p 偏振、500 nm 周期、250 nm 线宽、100 nm 深度。参数以每个 `report.json` 为准，不与旧 355 nm 案例混用。
- RCWA 颗粒试算为 100 nm 球、10 nm 间隙、0.5 × 1 μm 超胞。球会周期重复，不是孤立颗粒。它也不是有限两沟槽 FEM 的同一几何。
- Si 的 `n=5.57, k=0.387` 和颗粒 `n=1.59, k=0` 是沿用的数值对照输入，不能视为经过 266 nm 标定的材料常数。
- `|E_sample-E_reference|²` 与 `|E_sample|²-|E_reference|²` 不相同；后者包含相干干涉。
- 高斯光束物理建模方式的修改仍按用户要求暂缓，迁移时不要顺带修改。

## 2. 硬件与内存

原电脑物理内存约 31.50 GiB，属于标称 32 GB 配置。以下为同一 RCWA 三维颗粒场景的实测值，总模式数不是单轴傅里叶阶数。

| 总模式数 | 峰值进程计入内存 | 完整试算耗时 |
| --- | --- | --- |
| 117 | 0.106 GiB | 2.36 秒 |
| 319 | 0.371 GiB | 33.61 秒 |
| 633 | 1.373 GiB | 375.97 秒 |

最大算例运行时，系统至少还有 12.70 GiB 可用。当前 32 GB 机器足够继续小规模验证。

同层数、同后端下按密集矩阵近似平方增长外推：1200、2400、4000 总模式分别约需 4.93、19.74、54.82 GiB 进程内存。这些不是实测值，必须另留系统和工作区余量。保持相同空间频率分辨率时，超胞长宽均翻倍，模式数约增至四倍，密集矩阵内存约增至十六倍。

若迁移主要为扩大三维超胞，64 GB 可作为中等规模开发起点；128 GB 给更大规模验证留出余量，两者均不是最终收敛保证。优先保留足够 SSD 空间给 Python 环境、网格和场缓存。当前后端没有自动使用 GPU，购买显卡不能直接解决本轮内存或收敛问题。

有限光栅 FEM 三阶算例曾使用约 11.33 GiB 计入内存，不能用 RCWA 的 1.37 GiB 代替 FEM 内存估算。先保持 FEM 线程为 2，逐次测量，不要直接把线程数设成全部逻辑核心。

## 3. 获取项目

GitHub 仓库：<https://github.com/wufansaber/aptk>，私有仓库需使用有权限的账号。

迁移发布标签：`migration-2026-10-21`。优先下载该发布中的 `DFI_Simulation_Migration_20261021.zip`；它比 GitHub 自动生成的源码 ZIP 多包含环境快照及逐文件 SHA-256 清单。

也可以安装 Git 后克隆标签：

```powershell
git clone --branch migration-2026-10-21 https://github.com/wufansaber/aptk.git D:\code\aptk
Set-Location D:\code\aptk
git switch -c codex/migration-next
```

从 ZIP 获取时，完整解压到可写目录，在包含 `main.py`、`fullwave_prototype` 的目录运行以下命令。源码 ZIP 不包含 `.git` 历史；继续开发并同步时优先使用 Git 克隆。

```powershell
Get-FileHash .\DFI_Simulation_Migration_20261021.zip -Algorithm SHA256
```

将结果与同名 `.sha256` 文件比较。压缩包内 `MIGRATION_MANIFEST.json` 记录源码提交、文件大小、逐文件摘要和主动排除项；`ENVIRONMENT_SNAPSHOTS.json` 是原机环境版本记录，不包含凭据。

## 4. 重建三个独立环境

不要复制原电脑 `build/ngsolve_env`、`build/rcwa_env` 或 Anaconda 环境目录。venv 含有原机基础 Python 路径，复制后不能保证运行。不要把 NumPy 2.x 原生 FEM 依赖安装进旧 NIST2 环境。

### 4.1 原有 GUI 与 NIST 基准

仓库的 `environment.yml` 定义旧版 Python 3.8 / SCATMECH 环境。先按现有方案恢复，兼容性升级另行验证：

```powershell
conda env create -n NIST2 -f environment.yml
conda run -n NIST2 python -c "import sys, numpy, matplotlib, PyQt5, pySCATMECH; print(sys.version); print(numpy.__version__)"
conda run -n NIST2 python far_field_standalone\main.py
```

如 pySCATMECH 没有匹配的 Windows 二进制发行包，安装可能需要 C++ 编译工具；以安装器输出为准，不要混装不同 Python/NumPy ABI 的 DLL。仅运行已经发布的 EXE 时，应下载旧版完整分文件程序；该需求与源码开发环境不同。

### 4.2 FEM 原生环境

安装 64 位 Python 3.11 后，从项目根目录执行：

```powershell
py -3.11 -m venv build\ngsolve_env
& .\build\ngsolve_env\Scripts\python.exe -m pip install --only-binary=:all: -r fullwave_prototype\requirements-validation.txt
& .\build\ngsolve_env\Scripts\python.exe -m pip install pytest
& .\build\ngsolve_env\Scripts\python.exe -c "import ngsolve, netgen, numpy; print(ngsolve.__version__); print(numpy.__version__)"
```

当前测试使用 NGSolve 6.2.2607。GUI 默认会从项目内 `build/ngsolve_env/Scripts/python.exe` 查找求解器，入口规则见 `fullwave_ui/runtime_paths.py`。保持该相对目录即可，不需要复用旧机 `D:\anaconda`。

### 4.3 RCWA 独立环境

```powershell
py -3.11 -m venv build\rcwa_env
& .\build\rcwa_env\Scripts\python.exe -m pip install -r rcwa_prototype\requirements.txt
& .\build\rcwa_env\Scripts\python.exe -m pytest tests\test_rcwa_prototype.py -q -p no:cacheprovider
```

本轮上述测试为 17 项通过。`grcwa` 必须使用 requirements 中的官方 Git 提交 `c52f10a138f4b2d5c897b7670967813c36af8d9b`；不能仅执行 `pip install grcwa==0.1.3`，当前测试安装时 PyPI 没有该版本。安装 Git 引用需要 Git 和可用网络，不能用同名 `S4` 包随意替换后端。

本轮以 Laurent 因子化作为试算默认值；`--factorization pol` 保留作对照，不代表它在该高折射率算例中更准确。

### 4.4 PowerShell 和机器路径

若原电脑的 conda-hook.ps1 被执行策略阻止，以上 `conda run` 和直接调用解释器方式不需要激活脚本。可以在 Anaconda Prompt 执行 conda 命令；不要为了迁移全局禁用 PowerShell 安全策略。

`--reference-python`、旧 `.vscode` 解释器路径、Qt 保存的求解器路径、结果中的历史绝对路径均可能指向旧机。迁移包不携带个人 `.vscode` 配置。VS Code 中重新选择新机解释器；不要将报告里的原机路径当作必须存在的目录。

## 5. 按顺序验收

### 5.1 先检查测试和界面

```powershell
& .\build\ngsolve_env\Scripts\python.exe -m unittest discover -s tests -p test_fullwave.py -v
& .\build\ngsolve_env\Scripts\python.exe -m pytest tests\test_finite_grating.py tests\test_patterned_wafer.py -q -p no:cacheprovider
conda run -n NIST2 python far_field_standalone\main.py
```

确认界面可启动、FEM 求解器解释器存在，再做大算例。当前试验 RCWA 未接入 UI，不应在界面寻找一个已完成的 S4/RCWA 成像按钮。

### 5.2 运行小规模 RCWA 对照

```powershell
& .\build\rcwa_env\Scripts\python.exe -m tools.rcwa_trial --output build\migration_check\reference_121 --modes 121 --near
& .\build\rcwa_env\Scripts\python.exe -m tools.rcwa_trial --output build\migration_check\particle_121 --modes 121 --particle --near
& .\build\rcwa_env\Scripts\python.exe -m tools.rcwa_trial_report --root build\migration_check
```

原始参数快照及离散衍射级数据在 `output/rcwa_trial_20261021/`。如需重建 NIST 对照：

```powershell
conda run -n NIST2 python -m tools.grating_reference --output build\migration_check\nist245.npz --order 245
& .\build\rcwa_env\Scripts\python.exe -m tools.rcwa_trial_report --root build\migration_check
```

默认每任务预算 8 GiB，保留系统可用内存 3 GiB。启动前有规划检查，运行中约每 0.25 秒采样内存；它不是操作系统硬分配上限。状态必须为 `completed`，但完成本身不等于物理收敛。

### 5.3 再恢复有限光栅 FEM

先做几何预览，不立即启动高内存三阶求解：

```powershell
& .\build\ngsolve_env\Scripts\python.exe -m tools.finite_grating_run --output build\migration_check\finite_preview --preview-only
& .\build\ngsolve_env\Scripts\python.exe -m tools.finite_grating_run --output build\migration_check\finite_probe --order 2 --factor 1.8 --probe
```

确认自由度及空闲内存后，按保存的 `.resources.json` 中 `command` 重建参数，替换旧解释器和输出路径。不要覆盖历史结果。三阶原算例的预算为 16 GiB；32 GB 机器运行前应关闭自己的高内存任务，不要在仍有大量占用时直接放宽保护。

## 6. 验证结果如何解读

### RCWA 仍未收敛

依据 `output/rcwa_trial_20261021/summary.json`：

- 一维无颗粒 Laurent 算例的总反射率相对 NIST 差约 0.693%，但弱衍射级 -2 仍差约 7.63%。总反射率接近不能代替逐级验证。
- 三维模式网格中的裸光栅在 633 总模式下，总反射率相对 NIST 差约 -19.66%，尚不能视为合格周期背景。
- 319 到 633 总模式时，颗粒差分复场的功率加权相对变化约 15.71%。这是相邻截断变化，不是已知真实解的误差界限。
- 尚未证明超胞扩大后的孤立颗粒极限，离散级次效率也不是连续 DSC；不能直接送入原来的连续远场相机积分冒充相同数据。

### 有限光栅 FEM 仍未认证

依据 `output/finite_grating_20261020/comparison/comparison.json`，同几何和网格控制下，二阶与三阶选定 CH3 收集功率变化约 18.38%，复场变化约 29.82%。阶数变化也会改变默认 PML 阶数，需拆开检查。

三阶体源与闭合曲面提取的未加权复场差约 14.98% 至 15.67%；两包围面的复场差约 2.27%。小代数残差、较小积分功率差或两面接近均不能单独认证近场和远场。

周期无限光栅、有限沟槽、平整晶圆颗粒的验证必须分开。NIST 平整晶圆颗粒功率基准不直接适用于有图形的散射场，二维无限周期基准也不能直接认证有限三维端点加颗粒。

## 7. 新电脑上继续工作的优先顺序

1. 重建三套环境，通过第 5 节小算例；与保存的参数和功率比较，不做幅相拟合或强制对称。
2. 固定裸周期光栅几何，检查逐衍射级复场和功率、傅里叶截断、介电采样网格及因子化；先解决三维周期背景不一致。
3. 在合格背景上加入颗粒，分别做模式数、球体分层、横向采样收敛；再扩大超胞，监控周期副本作用。
4. 保留有限 FEM 的独立网格、物理域、PML、边角和远场提取验证。需要时评估 S4 或真正三维 AFMM，但先提供可重复的小基准。
5. 统一已验证散射复场的坐标、功率、相位和光瞳传播，再比较 sample/reference/真实强度差相机图及灰度 SNR。
6. 实验对比前补充 266 nm 材料光学常数、真实光栅高度/截面、照明与相机标定。SEM 灰度不能直接解释为三维高度。

## 8. 随包数据和另行备份

精选数据保留在源码原有 `output/` 相对路径，包含 RCWA 模式 NPZ、报告、资源记录，以及有限/周期光栅验证摘要和必要图像。包内清单以 `MIGRATION_MANIFEST.json` 为准。

原始 `volume_sources.npz`、`field.ngs` 不进入轻量 Git 迁移提交；它们留在旧电脑，清单另行记录大小及 SHA-256。需要直接复用旧 FEM 检查点时，再单独复制同目录的 `field.ngs` 和 `checkpoint.json`，并保持相同的 NGSolve 版本。当前 `field.ngs` 同时保存网格与场，只复制其中一个检查点文件不能安全恢复。只加载自己信任的检查点，不从未知存档反序列化网格。

`build/` 的历史粗糙面、相机校验及已打包 EXE 不随源码迁移；部分历史文档链接因此仅在旧工作区可用。旧版完整 EXE 可从历史 Releases 另行获取。删除旧电脑项目前，应先验收新机，并额外备份需要保留的 `build/`、用户结果目录和大型 FEM 检查点。

迁移包不会自动恢复 Qt 用户首选项、GitHub 凭据、浏览器登录、IDE 配置，也不会改变新电脑的系统策略。需要重新生成分文件 EXE 时，参照 `far_field_standalone/build_release.py` 和 `PACKAGING_VALIDATION_CN.md`，先完成新机环境与数值验收；本迁移包不声称当前试验 RCWA 已随 EXE 发布。
