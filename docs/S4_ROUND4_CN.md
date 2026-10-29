# S4 第四轮：联合加密、实现交叉核对和超胞先导

日期：2026-10-23。承接 [第三轮](S4_ROUND3_CN.md)。

## 最重要的结论

**截断与分层的联合局部加密已通过，但超胞尺寸仍明显影响收光，不能把当前周期颗粒阵列结果当作已收敛的孤立缺陷结果。**

| 检查 | 结果 | 含义 |
|---|---|---|
| S4 联合加密到 873 模、10 层 | 复差场变化 0.66399%，每路差场功率变化不超过 0.85% | 通过本轮局部门槛，不是误差上界 |
| 同截断 S4/grcwa 颗粒复差场 | 相差 0.32640% | 实现交叉核对，不是 FEM 独立方法认证 |
| grcwa 超胞 X 尺寸 0.5→0.6667 μm | 每颗粒周期差场截面 CH1/CH2 增加约 60%，CH3 减少约 4% | 孤立极限和孔径角采样未收敛 |
| 最终相关回归 | 111 项通过 | 包括原有回归、解析几何、比较约束和原始数据防覆盖 |

本轮保留 S4 主线、grcwa 辅助、FEM 独立对照的安排。现有 GUI 没有切换后端，没有输出已认证的相机图像或 SNR。

## 固定物理输入

266 nm、75°、p 偏振，光源来自 -KY，传播投影沿 +KY，全局 +Z 向上。光栅沿 X 不变、沿 Y 周期变化；周期 500 nm、线宽 250 nm、深度 100 nm。球直径 100 nm，球底距脊顶 10 nm。

晶圆 `n=5.57,k=0.387`，球 `n=1.59,k=0`，均为受控比较输入，**不是经过 266 nm 标定的真实材料数据**。采用单位振幅平面波，不把它直接换称为有限光斑的 1 W 功率；高斯照明修改仍暂缓。

三路独立，NA=0.6，倾角 +45°/-45°/0°。收光按实际模态功率求和，不按光锥几何体积计算。

## 联合加密

超胞仍为 0.5×1 μm；S4 继续使用固定版本和第二轮条纹层适配器，没有重建或修改原生求解核心。

| 起点 Mx/My/球切片数 | 联合终点 | 复差场变化 | CH1 / CH2 / CH3 差场功率变化 |
|---|---|---:|---|
| 3/40/8 | 4/48/10 | 0.66399% | +0.54966% / +0.54940% / +0.84930% |
| 4/40/8 | 4/48/10 | 0.71380% | +0.87318% / +0.87291% / +0.76578% |
| 3/48/8 | 4/48/10 | 0.59702% | +0.61246% / +0.61141% / +0.74235% |

功率变化仍按 `coarse/fine - 1`；正号表示联合加密后该量略小。复场使用同一全局 z=0 参考面、全部 20 个传播级和功率权重，不拟合振幅或相位。门槛为复差场变化小于 2%、每路差场功率变化小于 5%。

![联合加密](../output/s4_round4_20261023/convergence.png)

新增报告逻辑显式区分单变量和联合步骤。联合步骤只允许增加声明的数值参数，不能降低某一参数或同时修改材料、颗粒、照明、超胞等物理输入；原有单变量严格检查没有放宽。

这些结果只说明本次联合步没有出现超门槛漂移，不证明数值误差完全没有抵消。**终点是 10 层，不是第三轮单独检查过的 12 层，也不是 4/48/12 的计算结果。**

### 三路联合结果

| 通道 | 样品效率 | 无球参考效率 | 复差场效率 | 实际样品减参考效率 |
|---|---:|---:|---:|---:|
| CH1 | 0.00264044801 | 约 0 | 0.00264044801 | 0.00264044801 |
| CH2 | 0.00264045477 | 约 0 | 0.00264045477 | 0.00264045477 |
| CH3 | 0.10673839748 | 0.10361688285 | 0.00371876568 | 0.00312151463 |

效率的分母为该超胞平面波入射功率。`|E_sample-E_reference|²` 不等于 `|E_sample|²-|E_reference|²`；CH3 的交叉项效率约 `-0.00059725105`，没有丢弃。

背景相对 NIST 裸光栅 ±365 阶结果：总反射率差 -0.37778%，功率加权复场差 0.43296%，最差传播级功率相对差约 0.94713%。该比较仍仅认证裸周期背景的对应量，不认证球或有限光栅。球体 10 层中点切片的体积偏差为 +0.5%，未以振幅补偿。

## S4 与 grcwa 的匹配交叉核对

本轮使用实际安装的 **grcwa 0.1.3**，不是另一个 Python `rcwa` 库的 1.0.48 版本。新增功能在项目模块 `rcwa_prototype/sphere_fourier.py` 中，没有修改安装包。

以前 grcwa 的圆形颗粒边界由像素网格表示，S4 用解析圆形系数。现在可以显式选择 `--particle-method analytic`，使两者采用相同的球体中点切片和圆形边界，避免把像素化误差误认为求解器差异。默认 `fft` 路径仍保留。

对半径 r、面积 A=Lx*Ly 的圆盘，介电卷积矩阵为

```text
q = 2*pi*(Delta_m/Lx, Delta_n/Ly)
epsilon_GG' = delta_GG' + (epsilon_particle-1) * pi*r^2/A
              * 2*J1(|q|*r)/(|q|*r)
              * exp[-i*pi*(Delta_m+Delta_n)]
```

零频处 `2*J1(u)/u` 取极限 1。相位项将 grcwa 的胞角坐标中心转换为全局球心 `(0,0)`。圆盘层采用 Laurent 材料规则，直线光栅层采用既有 Li 规则；特征模和层间散射仍由 grcwa 求解，没有把有颗粒的整个结构缩减成一维。

测试包括独立极坐标数值积分、无损/有损颗粒、中心平移相位、与像素网格无关、与 S4 切片体积相同，以及独立体积公式 `V_slice/V_sphere=1+1/(2*Nslice²)`。**圆盘系数解析，不代表整个三维问题是无误差的解析解。**

实际交叉核对使用 `Mx=3,My=40,8 层`，两套实现均保留 567 模式：

| 比较量 | 相对复场差，分母为 S4 范数 |
|---|---:|
| 总样品场 | 0.0843306% |
| 无球参考场 | 约 9.88e-12% |
| 颗粒复差场 | 0.326404% |

三路差场功率的 `grcwa/S4-1` 分别为 -0.14830%、-0.14787%、-0.02927%。CH3 的实际功率差相差约 -0.75774%，再次说明不能用差场功率替代实际强度差。

![实现交叉核对](../output/s4_round4_20261023/crosscheck.png)

裸背景一致不能掩盖颗粒复差场仍有的差别。剩余偏差没有拟合掉，也没有被证明仅来自舍入误差。两者同属 Fourier 模态法，且共享本项目对材料规则的约定；必须继续保留不同离散方法的独立验证。

## 超胞先导揭示的限制

**本节后端是 grcwa，不是 S4。** 选择同一输入的 8 层球，保持实际频率截断 `Mx/Lx=6 cycles/μm`、`My/Ly=40 cycles/μm`：

- 小胞：Lx=0.5 μm，Ly=1 μm，Mx=3，My=40。
- 扩大 X：Lx=2/3 μm，Ly=1 μm，Mx=4，My=40。

X 沿直线光栅不变方向，因此这里不要求 Lx 是光栅 Y 周期的整数倍；Y 仍为两个原始周期。

不同超胞的效率分母随入射面积变化，不能直接比较原始效率。这里换算每胞一颗粒的周期差场截面：

```text
sigma_delta = Lx*Ly*cos(theta_i)*eta_delta
```

它仍是周期阵列的量，不提前称为孤立颗粒 DSC。

| 量 | 小胞截面 μm² | 扩大 X 后截面 μm² | 大胞/小胞-1 | 孔径内传播级数 |
|---|---:|---:|---:|---|
| CH1 差场 | 0.00034306777 | 0.00054840320 | +59.8527% | 4→8 |
| CH2 差场 | 0.00034306924 | 0.00054819826 | +59.7923% | 4→8 |
| CH3 差场 | 0.00048518884 | 0.00046583185 | -3.98958% | 8→12 |
| 半球差场 | 0.00229609869 | 0.00214698194 | -6.49435% | 全部传播级 |

![超胞敏感性](../output/s4_round4_20261023/supercell.png)

CH3 实际强度差对应的截面从 `0.00038951606` 变到 `0.00027682514 μm²`，其变化也不同于差场截面。

这个差异已经去除了简单的面积分母变化，但**仍同时包含周期复制相互作用与离散出射角/孔径采样变化**。两张超胞不足以把两者分离，不能把约 60% 全部解释成近邻耦合。改变胞长后，同一整数级次对应的实际波矢也不同，因此工具没有把它们当成相同方向直接相减复场。

当前结果不能用于认证孤立缺陷的收光率、定量成像或 SNR。仅提高同一小胞的截断，无法消除这个限制。

## 资源与计算边界

S4 的 873 模、10 层颗粒算例约 1528 s，进程树峰值统计 2.768 GiB；匹配背景约 181 s、1.289 GiB。grcwa 的 567 模/8 层颗粒算例约 127 s、1.044 GiB；扩大超胞后的 729 模/8 层约 389 s、1.681 GiB。不同规模和负载的时间不直接当作性能倍数。

最初请求的 `4/48/12` 在当时两次预检中被拦下，**没有启动，也没有这组原始计算数据**。其保守规划量约 4.7155 GiB，加系统保留 3 GiB，需要启动时可用内存约 7.72 GiB。`4/48/10` 是明确记录的中间联合算例，不是假称完成 12 层。

这些限制指当时可用内存，不等于整机总内存不足。没有降低 3 GiB 系统保留量、关闭用户其他应用或启用上一轮未通过的层复制优化。监测属于采样保护，不是操作系统硬上限。所有正式计算均正常完成，资源记录中的残留所属进程列表为空。

## 文件与复现

`output/s4_round4_20261023/` 保留：

- `plan.json`、`convergence.json`、`CONVERGENCE.md`、`convergence.png`：S4 联合检查。
- `crosscheck_plan.json`、`crosscheck.json`、`CROSSCHECK.md`、`crosscheck.png`：同截断实现核对。
- `supercell_plan.json`、`supercell.json`、`supercell.png`：grcwa 超胞先导。
- `admission_notes.json`：未启动请求及内存保护说明；最终状态以各实际算例资源记录为准。
- 各算例的 `modes.npz`、`report.json` 和同级 `*.resources.json`；S4 还保留 `input.lua`、`inputs.json`、标准输出与错误输出。
- `tests.xml`：最终 111 项测试通过；早期 `circle_tests.xml` 为解析圆形子集检查。

S4 与 grcwa 的 NPZ 都保留原始模态复电场，不只有绘图结果。公共字段为 `orders`、`electric_reflected_z0`、`directions`、`propagating`、`reflectance_by_order`；使用 `numpy.load(...,allow_pickle=False)` 读取。S4 的额外原生模态字段见 [第一轮读取说明](S4_VALIDATION_CN.md)。不能将这些字段名称直接冒充 SCATMECH 导出格式。

本报告还引用第二、三轮原始算例及 NIST 基准目录，迁移时需要一起保留。新的 grcwa 命令拒绝覆盖非空算例目录。

```powershell
# 求解示例须使用新的空路径
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce4/particle --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 4 48 --slices 10 --particle --timeout 3000
& build/rcwa_env/Scripts/python.exe -m tools.s4_trial --output output/s4_reproduce4/reference --executable build/s4_runtime/bin_stripes/S4.exe --stripe-rule --harmonics 4 48 --slices 10 --timeout 700
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_trial --output output/grcwa_reproduce4/particle --harmonics 3 40 --particle --particle-method analytic --grating-method analytic_li --timeout 900
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_trial --output output/grcwa_reproduce4/large_cell --harmonics 4 40 --lx 0.6666666666666666 --particle --particle-method analytic --grating-method analytic_li --timeout 1200

# 仅从已有数据重建报告
& build/rcwa_env/Scripts/python.exe -m tools.s4_convergence_report --config output/s4_round4_20261023/plan.json --output output/s4_round4_20261023 --plot
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_crosscheck_report --config output/s4_round4_20261023/crosscheck_plan.json --output output/s4_round4_20261023
& build/rcwa_env/Scripts/python.exe -m tools.rcwa_supercell_report --config output/s4_round4_20261023/supercell_plan.json --output output/s4_round4_20261023
```

## 后续优先级

1. 优先解决超胞与孔径角采样敏感性：扩大 X/Y，保持物理频率截断，再交叉增加截断密度，检查每颗粒量与收集积分。不要用插值平滑或强制对称替代缺失的物理计算。
2. 条件允许时完成 S4 `4/48/12`，再做更高联合步；将同几何 grcwa 交叉核对扩展到联合终点及扩大后的超胞。
3. 在匹配边界的条件下做平整晶圆颗粒基准及 FEM 独立离散方法对照。Bobbert-Vlieger 不能直接认证有光栅的颗粒结构。
4. 继续区分周期阵列、孤立缺陷和有限光栅端点。接入定量成像之前，仍需验证局部近场收敛、功率及真实 sample/reference 相干差；本轮没有跳过这些项目。

本轮未打包、提交或上传，也未改变现有 UI 的默认物理模型。
