# 低内存三阶局部加密与验证

日期：2026-09-27。接续 [近场峰值诊断](NEAR_CONVERGENCE_CN.md)。

## 1. 本轮实际改进

新增分区阶数的 H(curl) 空间。颗粒、粗糙面、空气和晶圆物理区保持三阶，只将
选定的 PML 单元及其非保护边/面降为二阶。物理区与 PML 共享的边/面仍保持三阶，
保留切向场连续性。几何曲边阶数不降低，Maxwell 弱形式、材料、照明和实际粗糙面都不改动。

这样降低的是外部吸收层的离散成本，不是删掉物理区的高阶场分量。
`field_space.py` 将空间构建独立出来，求解仍使用已经验证的直接求解器和静态凝聚。
单元测试逐个核对物理区单元自由度数及其边/面阶数，与原三阶空间一致。

此前超出内存的三阶局部尺寸倍率 0.9 算例，现在已经成功；倍率 0.85 也已成功。
但 PML 降阶引入额外误差，所以**它是可选低内存模式，不替换原默认配置**。

## 2. 界面和命令行

FEM 页新增三个设置：

- `PML order (0=same)`：0 保留原设置；2 表示 PML 内部使用二阶。
- `PML full-order fraction`：界面默认 0.5。仅在 PML 降阶时生效；按单元顶点
  平均位置判断其半径，在 `R + fraction*T` 内的单元保持原阶数。它是按单元
  选择的保护区，不是在该半径处重新划分几何。受保护单元共享的边/面也保持原阶数。
- `Local mesh scale`：默认 1。0.9 将颗粒与粗糙区目标单元尺寸减少 10%，
  0.85 减少 15%；不会改变空气、晶圆、PML 的目标尺寸或粗糙面的随机几何。

下方早期试验使用 `balanced`、PML order=2、full-order fraction=0、Local mesh scale=0.9。
这可运行，但独立 Mie 检查发现其精度代价较大，不能推荐为高精度方案。
0.85 单次近场求解接近内存上限，不建议直接据此认定完整配对流程也有足够余量。
不同尺寸、波长、材料和电脑需重新检查资源与收敛。

阶数和保护比例写入结果参数，局部倍率体现在颗粒及粗糙区网格尺寸中。
`--source` 复用旧源时拒绝同时修改这些设置，防止把旧结果冒充新网格结果。
读取历史报告时，缺失的 `pml_order` 和 `pml_inner_fraction` 均按 0 处理。
命令行保护比例默认 0，复算新折中方案时必须显式传入 `--pml-inner-fraction .5`。
`--near-diagnostics` 在完整联算的同一个有颗粒有限元解上导出密集近场，不能用于旧源重成像。

## 3. 内存与计算时间

均为同一粗糙面、355 nm、100 nm 球，参考间隙 10 nm、p 偏振 75 度。
除最后一行明确标注的失败配对试验外，本表是有颗粒的单次近场求解加密集采样，
时间含导出，不是完整配对成像耗时。

| 物理/PML 阶数 | 局部尺寸倍率 | 总自由度 | 计量内存峰值 GiB | 耗时秒 | 状态 |
|---|---:|---:|---:|---:|---|
| 3/3 | 1.00 | 261432 | 10.171 | 136.1 | 上轮成功基线 |
| 3/2 | 1.00 | 195365 | 8.074 | 106.1 | 成功 |
| 3/3 | 0.90 | 287172 | 11.818，停止时 | 36.9，停止时 | 上轮触发系统内存保留保护 |
| 3/2 | 0.90 | 220616 | 9.477 | 136.0 | 成功 |
| 3/2 | 0.85 | 236069 | 11.823 | 171.0 | 成功 |
| 3/内3外2 | 0.90 | 248183 | 12.296，停止时 | 69.8，停止时 | 本轮配对流程第一解触发内存保护 |

均保留 12 GiB 抽样上限和 3 GiB 系统可用内存保留，没有为了运行而放宽保护。
超限停止任务不是完成的物理解，不能拿它导出光场。抽样保护可能在两次采样间越过
阈值，因此报告中的峰值可以略大于配置上限；它不是操作系统级硬限额。
直接分解的内存还取决于网格连接和矩阵填充，不能用自由度数简单线性外推。

## 4. 三阶局部网格收敛

对比 0.90 与 0.85，均使用物理区三阶、PML 二阶，完全相同的粗糙面和采样坐标。
下表是相对 0.85 解的离散差异，**不是相对未知真解的误差上限**。

| 区域 | 复 E 相对 L2 | 强度相对 L2 | 采样峰值差异 |
|---|---:|---:|---:|
| 全 YZ 切片，401x401 | 0.45325% | 0.53128% | 0.07934% |
| 距粗糙面 >10 nm 且距球面 >5 nm | 0.22528% | 0.16602% | 0.16875% |
| 局部间隙，界面距离 >=1 nm | 0.35742% | 0.44522% | 0.02066% |
| 局部间隙，界面距离 >=2 nm | 0.34672% | 0.50312% | 0.07135% |

![三阶局部网格对比](../build/near_pml_order/h_convergence/near_convergence.png)

这次比较不再混用二阶细网格与三阶粗网格。另一方面，0.90 到 0.85 只是一次局部加密，
物理区其他网格未加密；不能仅凭一对网格接近就宣称严格渐近收敛。

## 5. 降低 PML 阶数的代价

固定局部倍率为 1，只把 PML 从三阶改为二阶：

| 区域 | 复 E 相对变化 | 采样峰值变化 |
|---|---:|---:|
| 全 YZ 切片 | 0.49807% | 0.91501% |
| 局部间隙 >=1 nm | 0.76034% | 0.91271% |
| 局部间隙 >=2 nm | 0.75977% | 0.91591% |

这说明低内存并非没有代价；上节小于 1% 的网格差异不能抹掉本节的边界敏感性。
两种变化同时引入后，原三阶/PML三阶/局部倍率1 与新三阶/PML二阶/倍率0.85 对比，
>=2 nm 区域复场相差 0.93848%，但峰值相差 1.69673%。不能把这个指标隐藏起来。

同一有限元解的取样加密，仍会让未限定界面距离的峰值向尖角靠近；本轮未声称解决连续
三维场最大值的奇异性问题，也没有人为平滑、改圆角或拟合幅度。

## 6. 独立 Mie 校验：保留内层 PML 精度

另用无衬底、无吸收的 100 nm 球，355 nm、n=1.59，与独立球谐 Mie 解比较。
不拟合幅度、相位或位置。物理区、PML 厚度和几何网格相同，只改变 PML 阶数。

| PML 配置 | 散射截面误差 | 消光截面误差 | 能量闭合差异 | 近场 E 最大相对 L2 | 近场 H 最大相对 L2 |
|---|---:|---:|---:|---:|---:|
| 全三阶，原配置 | 0.04097% | 0.17811% | 0.21869% | 0.09842% | 0.10352% |
| 全二阶，初始试验 | 0.21148% | 2.88423% | 2.59783% | 0.43705% | 0.26955% |
| 内半部三阶、外半部二阶 | 0.12570% | 0.76706% | 0.88596% | 0.16695% | 0.12852% |

E/H 指四个固定半径采样球壳中最大的相对 L2，而非连续空间最大点误差。
外半部降阶方案在最外采样球壳上，弱散射 E 的相对误差仍为 1.10656%，
不能用包含入射场的总 E 误差较小来掩盖这个指标。

结论：保留内层三阶能明显减轻全部降阶的边界误差，但原全三阶精度仍最好。
因此默认不变，新方案仅供内存受限时做局部加密和敏感性分析。
该自由空间基准不能替代粗糙衬底的独立解，也不能证明任意波长/颗粒均小于 1%。

- [全部 PML 降阶的 Mie 原始报告](../build/near_pml_order/mie/report.json)
- [只降外半部的 Mie 原始报告](../build/near_pml_order/mie_outer/report.json)
- [原全三阶的 Mie 原始报告](../build/theory_audit/o3_compact/report.json)

## 7. 未采用的方法

试验过两类迭代求解思路，但没有进入交付代码：

- BDDC 预条件 GMRES 在粗网格有损晶圆例子上收敛缓慢，未通过所要求的残差检查。
- 普通稀疏 ILU 预条件在该例子中出现数值溢出；虽然迭代器给出误导性的零预条件残差，
  独立原始方程残差仍为 1，因此被拒绝，没有导出为有效场。

没有保留无法通过验证的求解器选项。残差校验仍基于原始离散方程和有限场值，
不以预条件残差或迭代器退出状态代替它。
BDDC 的一般结构参考 [NGSolve 官方文档](https://docu.ngsolve.org/v6.2.2203/i-tutorials/unit-2.1.3-bddc/bddc.html)，
本例试验失败不表示该方法在所有电磁问题上不适用。

## 8. 复算与原始数据

```powershell
# 单次近场的两档三阶网格。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_pml_order/resources_h90.json --max-memory-gib 12 --timeout-seconds 600 near_audit --scale .9 --pml-order 2 --output build/near_pml_order/h90
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_pml_order/resources_h85.json --max-memory-gib 12 --timeout-seconds 600 near_audit --scale .85 --pml-order 2 --output build/near_pml_order/h85
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.near_comparison --base build/near_pml_order/h90 --refined build/near_pml_order/h85 --output build/near_pml_order/h_convergence

# 初始全部 PML 降阶试验，不作为高精度推荐配置。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_pml_order/resources_pipeline90.json --max-memory-gib 12 --timeout-seconds 900 imaging_run --accuracy balanced --pml-order 2 --local-mesh-scale .9 --pupil-samples 128 --output build/near_pml_order/pipeline90

# 保留内半部三阶，原局部网格，并导出同一个解的密集近场诊断。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_pml_order/resources_pipeline_outer.json --max-memory-gib 12 --timeout-seconds 900 imaging_run --accuracy balanced --pml-order 2 --pml-inner-fraction .5 --local-mesh-scale 1 --near-diagnostics --pupil-samples 128 --output build/near_pml_order/pipeline_outer
```

- [网格变化详细指标](../build/near_pml_order/h_convergence/comparison.json)
- [PML 阶数敏感性](../build/near_pml_order/pml_sensitivity/comparison.json)
- [与原配置同时改变后的差异](../build/near_pml_order/combined_change/comparison.json)

所有输出保存在本机 `build` 中，默认不会随 Git 上传。本轮只更新源码，不重新打包 EXE。

## 9. 最终完成的近场、远场与三路成像

成功配置：物理区三阶、内半部 PML 三阶、外半部二阶、局部尺寸倍率 1。
同一几何网格分别计算有颗粒/无颗粒复场，再用保存的体积源计算远场和三路独立像面。
总自由度 222265，凝聚后 181887；完整流程 313.7 秒，计量峰值 9.181 GiB，
系统可用内存最低 6.167 GiB，任务结束没有残留子进程。

物理参数：355 nm、球径 100 nm、颗粒 n=1.59/k=0、晶圆 n=5.57/k=0.387，
p 偏振、相对法线 75 度入射、间隙参考高度 10 nm；粗糙区域宽 0.2 um、RMS 2 nm、
5x5 高度节点、随机种子 7。照度为 1e-6 W/um^2，三路 NA=0.5，
倍率 20，物镜焦距 20 mm，管镜焦距 400 mm，光瞳采样 128x128。

| 独立通道 | 全三阶基线总像面功率 pW | 只降低外半部 PML 后 pW |
|---|---:|---:|
| +45 度 | 173.76203 | 173.24015 |
| -45 度 | 212.58837 | 211.98986 |
| 0 度 | 57.45877 | 57.40772 |

这里是含颗粒和粗糙背景的散射总光场功率，不含镜面照明光，不等同于背景扣除后的
颗粒信号。三路互相独立，不合束、不扣除跨通道重叠。离散传播功率误差约 2.54e-14，
仅说明传播积分一致，不能据此宣称电磁解达到该精度。

与相同网格的全三阶近场相比：

| 区域 | 复 E 相对 L2 | 强度相对 L2 | 采样峰值差异 |
|---|---:|---:|---:|
| 全 YZ 切片 | 0.13316% | 0.23769% | 0.36947% |
| 远离界面的区域 | 0.12450% | 0.13853% | 0.31253% |
| 间隙，界面距离 >=1 nm | 0.19154% | 0.36701% | 0.36933% |
| 间隙，界面距离 >=2 nm | 0.19145% | 0.36539% | 0.36994% |

近场比较以外半部降阶解为分母，与前文工具约定一致；不是严格真解误差。
全切片 101/201/401 点采样峰值仍随采样加密变化。连续三维最大值、粗糙面尖角、
更细几何、PML 截断、材料色散、真实照明和像差仍需分别研究，不能认定全部问题已解决。

![完成的近场、远场及三路像面](../build/near_pml_order/pipeline_outer/near_far_image.png)

- [完整像面参数与功率](../build/near_pml_order/pipeline_outer/imaging_report.json)
- [同一有限元解的密集近场](../build/near_pml_order/pipeline_outer/scattering/near_audit/report.json)
- [与全三阶近场的比较](../build/near_pml_order/outer_sensitivity/comparison.json)
- [完整流程资源记录](../build/near_pml_order/resources_pipeline_outer.json)
- [加密组合的内存保护记录](../build/near_pml_order/resources_pipeline_outer90.json)

## 10. 回归检查

- NIST2 环境：87 项测试，76 通过、11 因缺少原生依赖跳过。
- NGSolve 原生环境：全波、阶数、吸收边界、近场诊断、理论与资源保护共 56 项，
  52 通过、4 个 Qt 界面测试跳过。物理区及受保护 PML 单元自由度和共享边/面阶数
  均与原三阶空间逐一核对。
- 曾在 NGSolve 环境直接发现全部测试：4 个 SCATMECH/界面模块因该环境没有 PyQt5
  导入失败。没有掩盖该失败，也未向独立求解环境安装界面依赖；这些模块在 NIST2 全套中通过。
- 完成 1200x800、1680x980 窗口截图检查，新参数文字可见，三路像面正常显示。
- 原全三阶默认求解、旧报告默认值、旧源复用限制和用户已有改动均保留。
- NIST2 源码编译检查及 Git 已跟踪改动的空白检查通过。
