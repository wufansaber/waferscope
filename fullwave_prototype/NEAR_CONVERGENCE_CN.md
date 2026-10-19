# 近场峰值问题的专项诊断

日期：2026-09-27。对象是上一份报告中的“粗糙面采样峰值变化约 2.45%”。

后续进展：[低内存三阶局部加密](LOW_MEMORY_REFINEMENT_CN.md) 已完成此前内存不足的
三阶加密算例；本页保留当时的诊断和失败记录。

## 1. 结论和纠正

**可以修复统计口径并进一步检查数值收敛，但目前不能声称原来的全图峰值已经完全解决。**
2.45% 是两张 101x101 二维 YZ 切片的采样最大值之差，不是真实三维最大场的误差上限。
旧程序的 L2 指标排除了界面带区，而全图峰值没有排除；旧说明没有充分区分这两种指标。
现在明确分开输出全图采样峰值、排除带区后的峰值和固定界面距离的间隙指标。

新的距离计算使用实际粗糙三角面和球面，不把参考平面 z=0 当作真实粗糙面。
所有计算保持物理粗糙面、材料、入射和颗粒不变，没有平滑场值或拟合归一化。
保留原 balanced 默认设置，没有把未收敛的新配置直接替换进主程序。

## 2. 三种误差必须分开

- 有限元离散：修改场的阶数或单元尺寸，重新求解 Maxwell 方程。
- 取样离散：同一个已求出的场，用不同密度的显示点取值，不改变有限元解。
- 物理几何：三角面顶点、实际圆角、最小间隙改变时，解的问题本身也改变。

在理想边角处，Maxwell 界面问题可能存在奇异场，不能预先保证未限定距离的点峰值有
稳定的有限极限。依据：[Costabel、Dauge、Nicaise，Singularities of Maxwell interface problems，1999](https://numdam.org/item/M2AN_1999__33_3_627_0/)。
这项一般理论不是对本例特定有损介质角点已经发散的证明；本轮实际证明的是取样敏感性。

因此，不通过随意添加圆角改变几何来“制造收敛”。真实圆角或测量平均尺度需要有物理依据。

## 3. 同一三阶解，只加密取样

仍为 355 nm、100 nm 球、2 nm RMS 粗糙面、10 nm 参考间隙、p 偏振 75 度。
balanced 三阶重新求解后，原 101 点结果复现至数值舍入精度。

| YZ 切片取样 | 峰值 `|E/E_inc|^2` | 峰值 z 坐标 nm |
|---|---:|---:|
| 101x101 | 13.931184 | 5.76180 |
| 201x201 | 15.236358 | 4.60790 |
| 401x401 | 16.147641 | 4.03095 |

三者峰值都在 x=y=0；该位置的粗糙面顶点 z 约为 3.79974 nm。
加密取样时采样点越来越靠近这个顶点。不能把上表的增长解释为激光能量变大，
也不能把 401 点的值自动当作真实最大值。

## 4. 固定距离的局部检查

另取 x=0、y 在 [-15,15] nm 的局部间隙切片。按离实际粗糙面和球面的最短距离同时
排除 1 nm 或 2 nm 带区。对这个**明确限定的区域**，比较 61、121、241 点嵌套取样：

| 最小界面距离 | 三阶 61 点峰值 | 121 点峰值 | 241 点峰值 |
|---|---:|---:|---:|
| 1 nm | 14.883423 | 14.945881 | 14.946433 |
| 2 nm | 13.867471 | 13.867471 | 13.887834 |

121 到 241 点的峰值变化分别约 0.0037%、0.1466%。这是取样敏感性，不是有限元误差。
排除区域不属于“已经解决的全图峰值”，程序继续保留未排除的值，避免隐藏尖角问题。

## 5. 不同离散方式的交叉检查

内存可完成的一组为：二阶，球和粗糙区目标单元尺寸分别从 15/25 nm 减到 10.5/17.5 nm。
空气、晶圆、PML 尺寸及强度保持不变，与原 balanced 三阶在完全相同的采样点上比较。
阶数变化也改变球面曲边近似阶数，因此这是 hp 交叉检查，不是单变量的三阶 h 收敛。
两种结果都不是解析真值，以下数字是相对三阶结果的差异。

| 区域 | 复 E 相对 L2 | 强度相对 L2 | 区域采样峰值差异 |
|---|---:|---:|---:|
| 全 YZ 切片，401x401 | 1.0385% | 1.3507% | 3.2010% |
| 距真实粗糙面 >10 nm 且距球面 >5 nm | 0.7276% | 0.5204% | 0.3128% |
| 局部间隙，界面距离 >=1 nm | 1.3235% | 1.3335% | 0.9231% |
| 局部间隙，界面距离 >=2 nm | 1.3065% | 1.4422% | 0.8405% |

峰值指标小于 1% 并不意味着整片复场小于 1%；机器可读报告的 `gap_under_one_percent`
同时检查复场和峰值，所以该项为 **false**。不选择性地只展示较小指标。

![近场取样与网格交叉检查](../build/near_convergence/comparison/near_convergence.png)

### 再加密一级后的结果

二阶局部倍率进一步从 0.7 减到 0.6，颗粒/粗糙区的目标尺寸变为 9/15 nm。
相对原三阶 balanced，得到：

| 区域 | 复 E 相对 L2 | 强度相对 L2 | 区域采样峰值差异 |
|---|---:|---:|---:|
| 全 YZ 切片，401x401 | 1.0676% | 1.2546% | 2.7399% |
| 距真实粗糙面 >10 nm 且距球面 >5 nm | 0.6512% | 0.3936% | 0.2836% |
| 局部间隙，界面距离 >=1 nm | 1.3784% | 1.1931% | 0.7978% |
| 局部间隙，界面距离 >=2 nm | 1.1174% | 1.2300% | 0.1073% |

更细的二阶解不保证每项指标都单调接近三阶解；三阶解也并非精确参照。
其中 >=1 nm 区域的复场差异略增，不能隐藏这个变化。
两张二阶网格 0.7 与 0.6 直接比较时，>=2 nm 区域复场差异为 0.9253%，峰值差异为
0.7325%；>=1 nm 区域复场仍相差 1.3714%。因此仍没有通过整体 1% 的检查。

最新图和数据：[二阶 0.6 与三阶对比图](../build/near_convergence/comparison60/near_convergence.png)、
[完整指标](../build/near_convergence/comparison60/comparison.json)、
[两张二阶网格的收敛检查](../build/near_convergence/h_convergence/comparison.json)。

## 6. 资源限制和失败记录

本轮继续使用 12 GiB 抽样内存阈值和 3 GiB 系统可用内存保留，不解除保护。

- 原三阶求解成功，约 136 秒，计量峰值 10.17 GiB。
- 三阶局部网格倍率 0.9、0.95 均触发系统内存保留保护。
- 三阶 0.9 加粗外部 PML 网格的试算触发计量上限；这组没有结果，不能据此声称 PML 改动无影响。
- 二阶局部倍率 0.5 仍触发内存保留保护。
- 二阶局部倍率 0.7 成功，约 60 秒，计量峰值 5.64 GiB。
- 二阶局部倍率 0.6 成功，约 112 秒，计量峰值 9.72 GiB。

所有这些任务均已清理自有子进程。抽样保护允许短时超调，不是操作系统硬配额。
直接分解的内存受网格连接关系影响，不能只凭自由度数或目标单元尺寸推断峰值内存。
失败的资源记录保留在 `build/near_convergence/resources_*.json`，没有伪造对应场结果。

## 7. 如何复算

在项目根目录使用独立 FEM 环境，不需要改动 NIST2 的依赖：

```powershell
# balanced 三阶，附加密集取样诊断，只算有颗粒近场。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_convergence/resources_base.json --max-memory-gib 12 --timeout-seconds 600 near_audit --scale 1 --output build/near_convergence/base

# 不同离散方式的独立数值交叉检查。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_convergence/resources_h70p2.json --max-memory-gib 12 --timeout-seconds 600 near_audit --scale .7 --order 2 --output build/near_convergence/h70p2

# 比较固定坐标、材料和三角面；输出原始指标及图像。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.near_comparison

# 再加密一级，并保留原三阶作为交叉比较，而非精确参照。
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run --report build/near_convergence/resources_h60p2.json --max-memory-gib 12 --timeout-seconds 600 near_audit --scale .6 --order 2 --output build/near_convergence/h60p2
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.near_comparison --base build/near_convergence/h60p2 --refined build/near_convergence/base --output build/near_convergence/comparison60
.\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.near_comparison --base build/near_convergence/h70p2 --refined build/near_convergence/h60p2 --output build/near_convergence/h_convergence
```

完整的数组、峰值坐标、不同区域的样本数量见
[比较 JSON](../build/near_convergence/comparison/comparison.json)。
新增距离与掩膜回归测试 4 项通过；NIST2 全部测试 84 项，74 通过、10 个依赖相关项跳过。

## 8. 尚未解决的部分

没有证明任意靠近三角尖角的连续三维峰值已经收敛，也没有把有限元分布与本例的未知解析解比较。
要继续降低这部分数值不确定性，需要更省内存的局部高阶/迭代求解策略，或更多可用内存完成
三阶网格收敛序列；要预测真实表面最大增强，还需真实几何的圆角和空间分辨率。
目前不能把全图峰值标成“误差小于 1%”，也不能宣称已经通过本轮修复将它变成精确量。
