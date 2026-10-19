# 开源全波求解路线与本轮改进

2026-09-27：此前未完成的三阶平面晶圆成像基准已完成；粗糙面完整链路及资源限制见
[FOLLOWUP_VALIDATION_CN.md](FOLLOWUP_VALIDATION_CN.md)。本文保留前一轮实验记录。

## 决定

保留 NGSolve 的三维频域矢量 Maxwell 有限元路线，继续完善，不把现有原型替换为 FDTD。
这是针对当前 Windows、本地开源、单波长、颗粒加有损晶圆及有限粗糙面的选择，
并不是说 FDTD 不好，或有限元天然更准确。

FDTD 在时间域推进 Maxwell 方程，频域有限元直接求解同一方程的稳态复振幅。
对本项目的单色激光，两者都能描述多次散射、干涉和近场；数值精度必须通过
网格、边界、守恒和独立基准判断，不能由方法名称判断。任意粗糙面与颗粒的
一般组合不能保证存在闭式解析解；当前路线是数值全波解加解析背景/远场传播。

## 调研依据

| 方案 | 已核实的信息 | 本项目判断 |
|---|---|---|
| Meep | 官方不支持原生 Windows；Windows 安装建议 WSL | 不符合此前本机原生环境的约束 |
| openEMS | 开源 EC-FDTD，有原生 Windows 包、色散材料和复数远场输出 | 可作为未来独立交叉验证器，但不是更好的直接替换 |
| NGSolve | 已在本机运行曲面矢量有限元、吸收材料、PML 和复场 | 保留，优先降低已量化的离散/边界误差 |

openEMS 的标准 NF2FF 要求包围辐射结构，且包围面不能穿过有损材料。
本项目的有损半无限晶圆穿过有限计算区域，因此不能直接照搬其标准闭合盒
NF2FF 得到所需上半球场，需要另做分层介质/半空间远场变换和激励处理。
这并不代表 openEMS 无法模拟晶圆，而是替换求解器仍然不能免除这些物理实现与验证。
其官方 Windows 页面还提示已知性能问题，不能仅凭有安装包判定更适合。

本轮没有安装 openEMS，也没有做它与 FEM 的同参数速度或精度实测，
所以不声称已通过横向跑分证明 FEM 胜过所有 FDTD 软件。

一手来源（本轮核查）：

- [Meep 安装说明](https://meep.readthedocs.io/en/latest/Installation/)
- [openEMS Windows 包及限制](https://docs.openems.de/en/latest/install/package.html)
- [openEMS 近远场变换的适用条件](https://docs.openems.de/en/latest/concepts/nf2ff.html)
- [openEMS 色散材料](https://docs.openems.de/en/latest/concepts/dispersive.html)
- [NGSolve Maxwell 示例](https://ngsolve.org/ngsolve/docs/i-tutorials/wta/maxwell.html)
- [NGSolve 静态凝聚](https://docu.ngsolve.org/latest/i-tutorials/unit-1.4-staticcond/staticcond.html)

## 已实现的改进

1. 新增吸收层独立网格比例 `pml_mesh_scale`，不必同时细化颗粒网格。
   几何网格器可能将细化传播到邻接物理区，不能理解为物理区网格绝对不变。
2. 新增可复现的径向渐变吸收 `pml_power=1/2/3`，共同使用同一外边界复位移：

   `r_tilde = r + i * strength * T * ((r-R)/T)^power`

   `power>=2` 在吸收层入口的变换雅可比连续为单位阵。但这不是精度保证，
   本轮二次渐变实验未胜出，故默认仍为线性径向 PML。
3. 新增高阶静态凝聚：消去单元内部未知量，再完整恢复内部电场。
   不降阶、不裁剪复场，近场、磁场和体源仍使用恢复后的完整解。
   已用同一网格比较未凝聚/凝聚的所有场系数及 E/H。
   凝聚计算的 `algebraic_residual` 明确标注为 Schur 补系统残差，
   不冒充完整未凝聚矩阵残差，更不冒充物理误差。
4. 新增 `preview` 和 `refined` 数值配置。精度模式 `refined` 使用
   三阶场、静态凝聚、原线性 PML、七阶体源积分。默认仍为旧二阶 `preview`，
   因为本轮高阶大模型复核曾明显拖慢本机，不宜强制提高默认资源需求。
   底层 `Parameters()` 和旧报告的默认值不变，旧结果不会被自动“升级”。
5. 已保存体源的重新成像禁止同时传入新的数值配置或吸收层参数，
   避免用户以为改了求解设置、实际仍复用旧光场。

## 自由球实测筛选

无基底自由球，355 nm、直径 100 nm、折射率实部 1.59。所有误差均相对原始
绝对值计算，没有通过峰值归一化或拟合振幅减小差异。

| 配置 | 散射截面相对 Mie 误差，k=0 | 能量闭合偏差，k=0 |
|---|---:|---:|
| 原二阶、线性 PML | 0.5073% | 4.7797% |
| 二阶、二次渐变 PML | 1.2486% | 4.5338% |
| 二阶、仅 PML 网格尺度乘 0.65 | 0.2464% | 1.9067% |
| 三阶、线性 PML、静态凝聚 | 0.0490% | 0.8087% |

三阶有吸收测试球 `k=0.1` 的散射误差为 0.00871%，消光误差为 0.3693%，
能量闭合偏差为 0.3595%。`k=0.1` 只是守恒测试，不是声称 PSL 的实际吸收系数。
三阶两种球的近场积分散射通量与远场积分相差约 0.81% 至 0.83%。
这些限定算例通过列出的 1% 筛查项，但不构成所有工况的误差界。

凝聚前后三阶无吸收球的散射结果在浮点精度内一致；本轮凝聚计算耗时约 79 秒，
此前未凝聚同阶约 102 秒（非隔离环境跑分，仅供参考）。实际参与分解的未知量
为 185580，完整场自由度仍为 247080。

原始实验：`build/boundary_upgrade/graded2/`、`refined_pml/`、`order3_condensed/`；
旧基线：`build/physical_audit/`。

## 粗糙晶圆边界检查

固定同一粗糙面、同一张三阶曲面网格，分别设置 PML 强度为 1、2、3。
仅计算有颗粒场，以强度 2 为参照，不拟合振幅或相位。

| 指标 | 原二阶最大变化 | 三阶凝聚最大变化 |
|---|---:|---:|
| 上半球散射积分 | 1.6077% | 0.5340% |
| 四个近场探针的复电场 L2 | 1.0612% | 0.3239% |

三阶强度 2 的半球散射截面为 `0.001393650984 um^2`，使用 Gauss 球面求积，
不是几何体积占比。三次求解均为 298136 个完整场自由度，实际分解 228360 个；
每次 Maxwell 求解约 178 至 185 秒，后处理耗时另计。
近场探针不在接触奇点；该结果不是整个近场区域的逐点误差界。
这里只改变吸收强度，不能替代吸收层厚度、计算域大小和粗糙几何的系统收敛分析。

## 资源限制与未完成项

高阶自由球与粗糙面 PML 扫描已经完成。随后三阶平面晶圆成像检查的
UMFPACK 计算造成明显资源压力，主动终止；单线程 BLAS 配合
`sparsecholesky` 的大网格重试也未完成，已停止。没有把中断任务记为通过，
也没有用旧成像图冒充本轮更新后的图像。完整三阶“近场-远场-三路像面”
输出本轮未重新生成，原有二阶输出仍可用。

新增可选的对称直接求解器 `direct_solver='sparsecholesky'`，已通过有损晶圆
小网格的完整场系数对照，但尚未通过本轮大网格物理基准，因此不是默认求解器。
任一直解器出现非有限场或代数相对残差大于 `1e-8` 都会拒绝导出。
该数值残差门槛不是物理误差门槛。

## 复现

在项目根目录运行，解释器仍为独立环境，不改动 NIST2：

```powershell
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.physical_audit --case mie --order 3 --static-condensation --output build/boundary_upgrade/order3_condensed
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.physical_audit --case pml --order 3 --static-condensation --output build/boundary_upgrade/order3_condensed
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.physical_audit --case flat --order 3 --static-condensation --output build/boundary_upgrade/order3_condensed
& '.\build\ngsolve_env\Scripts\python.exe' -m fullwave_prototype.imaging_run --accuracy refined --output build/fullwave_imaging/refined
```

快速预览添加 `--accuracy preview`。`--order` 是显式覆盖，实际设置写入输出报告。
高阶粗糙面验证中曾观测到约 18.5 GB 的进程工作集，这不是内存峰值上限；
内存不足时使用预览配置。自由度限制只能挡住规模过大的网格，不能保证内存足够。
三路仍相互独立，默认倾角 +45/-45/0 度、NA 0.5、物镜焦距 20 mm、倍率 20x、
tube 焦距 400 mm。默认验证光学参数为 355 nm、100 nm 球、10 nm 间隙、
有限粗糙斑 RMS 2 nm；不是 266 nm 结果。折射率是测试输入，不是材料测量数据库。

## 尚不能声称的能力

不能称为“已全局最优”或“全部工况误差小于 1%”。尚需其他波长/粒径/材料、
接触点局部近场、粗糙统计系综以及独立粗糙面求解器的基准验证。
当前全波原型仍没有多层膜、任意颗粒形状、有限光斑激励和相机噪声的完整链路。
它仍独立于旧 UI，本轮没有把原型接入 UI 或重新打包 EXE。
后续加入膜层必须同时升级背景场和半空间辐射算子，不能只改变网格材料。
