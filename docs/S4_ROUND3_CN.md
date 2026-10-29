# S4 第三轮：分轴截断与球体分层

日期：2026-10-23。承接 [第二轮](S4_ROUND2_CN.md)。

## 本轮结论

**三个单变量加密步骤均通过局部容差检查，但尚不能宣布完整模型收敛。**

| 单独变化的量 | 其余设置 | 模式数 | 功率加权复差场变化 | CH1 / CH2 / CH3 差场功率变化 |
|---|---|---:|---:|---|
| Mx：3→4 | My=40，8 层 | 567→729 | 0.43183% | -0.32072% / -0.32070% / +0.08289% |
| My：40→48 | Mx=3，8 层 | 567→679 | 0.27596% | -0.06242% / -0.06163% / +0.10617% |
| 球切片：8→12 | Mx=3，My=40 | 567→567 | 0.99487% | +1.38413% / +1.38413% / +1.05864% |

本轮门槛为复差场变化小于 2%、每路差场功率变化小于 5%。功率变化按 `coarse/fine - 1` 计算；负号表示细化后功率增加。**这是相邻计算的变化，不是真实解误差上界，也不是能量守恒的独立验证。**

![收敛检查](../output/s4_round3_20261023/convergence.png)

第二轮在更低截断下的 X/Y 变化仍为 2.6574%/3.9592%，保留原记录，不用本轮数值覆盖。本轮三项从共同基点分别出发，尚未计算 `Mx=4,My=48,12 层` 的联合加密解。

## 输入与求解器

- 真空波长 266 nm，入射角 75°，p 偏振，光源来自 -KY，传播投影沿 +KY；全局 +Z 朝晶圆外。
- 光栅沿 X 不变、沿 Y 周期变化：周期 500 nm，线宽 250 nm，深度 100 nm。
- 超胞为 0.5×1.0 μm，每胞一颗直径 100 nm 的球，球底离脊顶 10 nm。
- 晶圆 `n=5.57,k=0.387`，球 `n=1.59,k=0`，均为沿用的受控对照输入，**未经 266 nm 材料标定**。
- 仍使用固定提交 `5b89848cfd15719ccaa04c853fe883ecf8d6ac15` 的 S4、项目条纹层 Li 适配器和矩形基组；球形切片保留原 S4 材料路径。
- 二进制为 `build/s4_runtime/bin_stripes/S4.exe`，SHA256 为 `1da0f38040f8f8b168269280fd7c013c536a84a252bcb7f2fb08b898f2271741`。本轮未修改原生数值核心或重建二进制。

所有算例仍是完全耦合的二维横向模态展开，不把有球结构缩减为一维。它表示周期球阵列，不是有限光栅上的孤立颗粒。高斯照明修改仍暂缓。

## 比较方式

先以同一设置计算有球场和无球参考场，在相同全局脊顶参考面上构造复矢量差：

```text
delta_E_g = E_sample,g - E_reference,g
w_g = (kz_g/k0) / cos(theta_incident)
relative_field_change = sqrt(
    sum_g w_g * |delta_E_coarse,g - delta_E_fine,g|^2
    / sum_g w_g * |delta_E_fine,g|^2
)
```

这里空气侧 `k0=2*pi/lambda`，`kz_g/k0` 为出射方向的 Z 分量。求和覆盖全部共同传播级，包含三分量复场。本轮各计算均保留同样的 20 个传播级。不做振幅拟合、全局相位对齐或左右平均；相同功率而不同相位也会被识别为不一致。

新增 `tools/s4_convergence_report.py` 使用显式比较计划。每一步只允许一个声明的数值参数增加；材料、入射、几何、其他截断参数、材料规则或求解器二进制不同会拒绝比较。输入 Lua、原始 NPZ 和 NIST 基准均检查 SHA256。

12 层球算例复用相同截断下的裸背景原始文件，**没有把旧参考改名为新计算**。无球输入中的球顶至脊顶间只有一层均匀空气，厚度由球直径和间隙决定，和球切片数量无关。工具仅允许这种 `slices` 差异；参考面、材料、超胞或截断不同不允许复用。来源和复用标志写入报告。

## 三路收光

三路独立，均为 NA=0.6，光轴倾角 +45°、-45°、0°；方向分别为 `(sin(t),0,cos(t))`。下面给出 **12 层、Mx=3、My=40** 的结果，不将其称为联合最高精度解。

| 通道 | 有球样品效率 | 无球参考效率 | 复差场效率 | 实际样品减参考效率 |
|---|---:|---:|---:|---:|
| CH1 | 0.00261871505 | 约 0 | 0.00261871505 | 0.00261871505 |
| CH2 | 0.00261871505 | 约 0 | 0.00261871505 | 0.00261871505 |
| CH3 | 0.10663045922 | 0.10356025658 | 0.00371106254 | 0.00307020264 |

效率是该超胞收集功率除以该超胞平面波入射功率，不是 W、几何体积占比或孤立颗粒 DSC。CH1/CH2 各收集 4 个传播级，CH3 收集 8 个。

CH3 的裸光栅衍射背景不为零，因此

```text
|E_sample|^2 - |E_reference|^2
  = |delta_E|^2 + 2*Re(E_reference* conjugate(delta_E))
```

其中矢量项取点积。CH3 本例的相干交叉项效率约为 `-0.00064085990`，不能把差场效率 `0.00371106` 当作真实强度差 `0.00307020`。此差异不是 haze 或相机噪声，本轮没有加入这两项。

12 层结果的 CH1/CH2 相对不对称度约 `1.92e-9`，分母为两路均值；本轮所有算例的这一数值小于 `7.86e-6`。没有强制左右相等。

## 背景与球几何检查

| 裸背景设置 | 对 NIST 总反射率相对差 | 功率加权复场差 |
|---|---:|---:|
| My=40，Mx=3 或 4 | -0.62266% | 0.68704% |
| My=48，Mx=3 | -0.37778% | 0.43296% |

NIST 基准仍为同几何裸直线光栅的 ±365 阶结果。这里是不同截断的敏感性检查，不能宣称它认证了三维颗粒。第二轮 `Mx=1,My=120` 的高阶背景结果仍保留，其误差更小，不和本轮较低 My 的结果混用。

球体中点切片的体积相对偏差由 8 层的 +0.78125% 降为 12 层的 +0.34722%。没有修改球折射率、半径或输出振幅来补偿体积。体积接近也不等于场已收敛。

## 性能与未采用的优化

| 颗粒算例 | 监测总时间 | 进程树峰值统计 |
|---|---:|---:|
| Mx=4，My=40，8 层 | 766 s | 1.811 GiB |
| Mx=3，My=48，8 层 | 609 s | 1.609 GiB |
| Mx=3，My=40，12 层 | 418 s | 1.459 GiB |

本轮逐个运行主要算例，限制 8 GiB、系统保留 3 GiB；统计是采样监测，不是操作系统硬内存上限。最小可用物理内存仍约 8.72 GiB，当前主要限制是密集特征值求解时间。

另外试验了 `AddLayerCopy` 复用球体上下相同切片。远场很接近，但 5 层回归的近场 E/H 范数分别出现最高约 5.44%/2.90% 差异。原因尚未分离，不能断言所有 S4 版本有同样问题。**这项提速没有通过回归，已从活动代码撤回，本轮正式算例未使用。** 原始成对输入、复场及说明保留在 `output/s4_round3_20261023/rejected_layer_copy_audit/`，不混入已接受结果。

## 原始数据与复现

主目录：`output/s4_round3_20261023/`。

- `plan.json`：四组样品/参考的真实来源和三个单变量步骤。
- `convergence.json`、`CONVERGENCE.md`、`convergence.png`：完整量化结果及图表。
- 各算例下的 `input.lua`、`inputs.json`、`stdout.txt`、`stderr.txt`、`modes.npz`、`report.json`：原始输入、模态复场、逐级功率和构建清单。
- 算例同级的 `*.resources.json`：时间、峰值内存、退出码及进程清理状态。
- `tests.xml`：最终 **88 项测试通过**。旧的 `layer_reuse_tests.xml` 是被拒绝的提速试验记录，不是最终通过清单。

本报告还引用 `output/s4_round2_20261022/` 的基点和参考，以及 `output/rcwa_grating_20261022/nist365.npz`；迁移时需要一起保留，不能只复制第三轮文件夹。NPZ 的 S4 字段规则见 [第一轮记录](S4_VALIDATION_CN.md)，不要直接套用 SCATMECH 导出字段名。

```powershell
# 使用新路径；求解器拒绝覆盖非空的算例文件夹
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce/reference_rect_4_40 --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 4 40 --timeout 500
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce/particle_rect_4_40 --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 4 40 --particle --timeout 1500
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce/reference_rect_3_48 --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 3 48 --timeout 500
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce/particle_rect_3_48 --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 3 48 --particle --timeout 1500
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce/particle_rect_3_40_s12 --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 3 40 --slices 12 --particle --timeout 1500

# 从已保存原始结果重建本轮报告，不重新求解
& build/rcwa_env/Scripts/python.exe -m tools.s4_convergence_report --config output/s4_round3_20261023/plan.json --output output/s4_round3_20261023 --plot

& build/rcwa_env/Scripts/python.exe -m pytest tests/test_s4_prototype.py tests/test_s4_stripes.py tests/test_s4_convergence_report.py tests/test_rcwa_prototype.py tests/test_rcwa_grating_fourier.py -q -p no:cacheprovider --basetemp=build/pytest_s4_round3_rerun
```

## 下一步边界

1. 做联合高截断/细分层计算，再增加一档，排查不同数值误差的交叉或抵消。
2. 扩大超胞时按物理空间频率保持截断密度，比较每颗粒的量和收集积分；不同超胞传播级不同，不能直接逐级相减或比较入射功率分母变化后的效率。
3. 平整晶圆颗粒与适用的 Bobbert-Vlieger 基准比较，再做同几何 S4/grcwa 复场交叉检查及 FEM 独立方法对照。
4. 本轮未验证局部近场的截断收敛，也未建立有限光栅端点或孤立颗粒结论。以上完成后再接定量成像及相机 SNR。

现有 GUI/FEM 没有切换到试验后端；本轮未打包、提交或上传。
