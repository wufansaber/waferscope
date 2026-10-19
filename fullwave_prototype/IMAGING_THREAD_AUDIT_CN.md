# 独立成像校验与本机线程测试

## 校验范围

本轮针对近场到远场计算后的理想三路成像，以及现有电脑的 FEM 线程配置。
不修改材料常数、颗粒形状、照明方向、物镜 NA 或真实仪器的假设。
三路保持独立，不将三个像面相干叠加。

采用三层检查，避免把同一套传播程序的两次调用当作独立验证：

1. 独立直接积分与圆孔 Bessel/Hankel 解比较，含高像方 NA、离焦和 20 倍成像。
2. 读取已有粗糙面计算的复数远场，用独立坐标旋转和逐光线积分复算三路像面。
3. 对平面晶圆，使用 SCATMECH 的 Bobbert–Vlieger 复振幅和独立极坐标积分，
   对照已有 FEM 像面；分别测试 p、s 入射。检查积分阶数和多极截断收敛。

第二层保留同一 FEM 出瞳，因此验证传播实现，不独立验证 FEM 源。
第三层同时更换散射求解方法及传播积分方法，但只适用于平面晶圆上的球形颗粒。
不能将平面参考的误差当作随机粗糙表面、镀膜或任意颗粒的误差保证。

## 2026-09-28 实测结论

工况为 355 nm、直径 100 nm 球、颗粒 n=1.59/k=0、晶圆 n=5.57/k=0.387、
球底间隙 10 nm、入射角 75 度；三路 NA=0.5、倾角 +45/-45/0 度、20 倍成像。
这些是固定的模型参数，不是本轮测量得到的材料色散。

### 成像独立校验

- 直接矢量积分对 Bessel/Hankel 解：复电场/磁场相对 L2 误差最高 7.08e-15。
- 已有粗糙面三路复场对独立直接积分：最高 4.45e-16；无颗粒背景也逐路检查。
- 独立 BV 积分的径向/方位采样从 24x96 加密至 48x192：强度相对变化最高 2.62e-15。
- BV 多极截断从 14 增至 18：像面强度相对变化最高 7.72e-11。

FEM 的 128 点瞳面像面与独立 BV 极坐标积分的强度相对 L2 误差：

| 入射偏振 | +45 度 | -45 度 | 0 度 |
|---|---:|---:|---:|
| p | 0.1623% | 0.1691% | 0.4916% |
| s | 0.5259% | 0.3307% | 0.4238% |

误差在像面 x/y 均约 ±12 um 的共同采样点计算，没有归一化或平移配准。
这同时包含 FEM、瞳面采样及传播误差，不是所有像素的最大相对误差。
独立积分和 FEM/FFT 不共用传播函数。
此前保存的 32 点瞳面结果最高误差为 **2.1076%**，超过本轮预先设定的 2% 门槛；
保留失败记录在 `build/imaging_independent/current`，重新计算 128 点结果后通过，
而非放宽验收标准。生产 UI 原本就默认 128 点，本轮升级的是验证基准的采样。
仅将 BV 用 FFT128 传播对比连续积分，残余瞳面采样误差约 0.29%–0.39%。

粗糙面中路在像面 (0, 4.4375) um 的直接积分强度为 5.67575e-15 W/um²，
是该区域峰值的 0.8845%。平面晶圆的独立 BV 模型也在 p 偏振下出现相似暗区，
细网格上暗区约位于 (0, 4.5) um；s 偏振形态不同。
这排除了本工况中 FFT 漏点/显示缺块的解释，支持偏振与相干干涉的模型解释，
**不等于粗糙面绝对强度已经通过实验验证**，也不应人为填成完整圆斑。

通过的报告和图：

- `build/imaging_independent/refined/report.json`
- `build/imaging_independent/refined/rough_independent.png`
- `build/imaging_independent/refined/flat_bv_independent.png`
- 高采样 FEM/BV 原始基准：`build/imaging_independent/flat128/flat_p`、`flat_s`

### 线程基准

本机 Core Ultra 5 225H / 32 GB；NGSolve 6.2.2607。
使用 222265 自由度、静态凝聚后 181887 自由度的粗糙面算例；
物理区域三阶，PML 内半段保持三阶，外半段二阶，sparsecholesky 求解器。
每组重复两次，第二轮反序执行。

| NGSolve 线程 | 求解时间中位数 | 两次范围 | 峰值记账内存 |
|---:|---:|---:|---:|
| 1 | 337.16 s | 335.74–338.59 s | 8.15 GiB |
| 2 | 118.33 s | 95.44–141.21 s | 9.04 GiB |
| 4 | 67.42 s | 57.18–77.65 s | 11.00 GiB |
| 8 | 50.39 s | 45.88–54.90 s | 11.44 GiB |

在已测试的线程数中，8 线程最快；相对 2 线程中位数为 **2.35 倍求解速度**，
时间减少约 **57.4%**。这是单次有颗粒 FEM 求解时间，不是整个成像流水线的倍数。
复电场、磁场和远场振幅探针最大相对变化 2.97e-14，所有代数残差通过。
原始数据在 `build/thread_audit/rough/summary.json`。
未测试所有 1–14 线程组合，不宣称全局最佳；电源模式、后台程序和温度会影响时间。

8 线程已接近本算例 12 GiB 保护上限，建议在有足够空闲内存且使用同样数值配置时使用；
需要更多余量可选 4 或 2。不要据此直接增加网格规模或关闭内存保护。

### 完整流程与回归

另以 8 线程完整重算有/无颗粒、近场诊断、远场、128 点出瞳和三路成像：

- 总耗时 173.91 s，峰值记账内存 11.54 GiB，最低系统可用内存 4.20 GiB。
- 无颗粒/有颗粒 FEM 求解分别为 43.11 / 45.80 s；后续远场与成像并非同样比例提速。
- 与原 2 线程结果相比，除线程数外物理、网格和成像设置相同。
- 全部比较的像面复场/强度/干涉项最高相对差异 2.15e-14；近场/远场及场变化最高 5.92e-14。
- 三路像面积分功率为 173.2401 / 211.9899 / 57.4077 pW，保持独立计算。

完整结果目录：`build/thread_audit/pipeline8`；逐项对比：`thread_comparison.json`。
资源保护记录：`build/thread_audit/pipeline8_resources.json`。
可在界面的 `Load result` 中打开该目录查看三位置结果；加载结果不会自动替换当前输入设置。

回归测试：NIST2 运行 100 项，89 项通过、11 项因环境缺少对应可选后端跳过；
原生 NGSolve 环境的相关 41 项全部通过（两组存在交集，不应直接相加）。
新增测试还覆盖错误复场/非有限值拒绝、无效采样提前拒绝、CLI/UI 线程参数传递、
资源限制及不推荐错误/失败的快速试算。Qt 离屏截图已检查线程行的显示。

## 独立传播公式

令物镜局部方向为 $(u,v,w_o)$，其中 $w_o=\sqrt{1-u^2-v^2}$；
倍率 $M$ 下像方光线方向为

$$\mathbf q=(-u/M,-v/M,w_i),\qquad w_i=\sqrt{1-(u^2+v^2)/M^2}.$$

将远场 s/p 复振幅转为全局三维矢量，转入物镜坐标后，用 Rodrigues 公式实现
收集光线到光轴、以及光轴到像方光线的两次最小旋转。
偏振片作用于准直区的横向场；功率透射率为 $T$，场乘以 $\sqrt T$。
记旋转及偏振片之后的场为 $\mathbf F(u,v)$，则本程序使用的通量标度为

$$\mathbf E(x,y,z)=\frac{\sqrt T}{\lambda M}
\iint_{\mathrm{pupil}}\frac{\mathbf F(u,v)}{\sqrt{w_o w_i}}
\exp\!\left[\frac{2\pi i}{\lambda}
\left(-\frac{ux+vy}{M}+z(w_i-1)\right)\right]\,du\,dv.$$

这里 $\mathbf F$ 不包含透射率；计算磁场时将被积场替换为
$\mathbf q\times\mathbf F$。电磁场采用通量标度，不是 SI 电场/磁场单位。
像面强度为

$$I=I_{\mathrm{inc}}\operatorname{Re}(E_xH_y^*-E_yH_x^*).$$

注意积分变量是投影瞳面 $du\,dv$，不是直接以其代替立体角。
$d\Omega=du\,dv/w_o$，上述权重包含消球差物镜的振幅变换。
极坐标检查使用 $du\,dv=r\,dr\,d\phi$。
不做峰值归一化，不拟合平移或强度系数，不用 FFT 生成独立参考图。
BV 的球心原点通过已知高度相位移到晶圆表面原点，而非拟合相位曲面。

理论依据：[Richards–Wolf 矢量衍射与 Cartesian / spherical 参数化](https://psf-generator.readthedocs.io/src/theory/basic.html)。
独立实现为 `direct_imaging.py`；生产传播仍为 `optical_propagation.py`。

## 线程测试方法

固定同一组物理参数、粗糙面种子和单线程网格生成过程，分别改变 NGSolve 工作线程。
各任务依次运行，不并行启动多个高内存 FEM 任务。
记录装配、分解、总求解时间、进程树峰值内存，以及复数近场/远场探针差异。
不同线程的电场、磁场、散射振幅相对差异必须小于 $10^{-8}$。
至少重复两次才给出自动建议；速度差距在 5% 内时优先少线程。

内存保护维持 12 GiB 上限、3 GiB 系统可用内存预留，每 0.25 秒采样。
这是采样保护，不是操作系统硬内存配额，瞬时尖峰可能高于采样记录。
BLAS/OpenMP 始终限制为 1；调整的是 NGSolve 的 `SetNumThreads`，避免嵌套并行。
依据：[NGSolve 共享内存并行文档](https://docu.ngsolve.org/latest/how_to/howto_parallel.html)。

## 复现入口

在项目根目录用现有隔离环境运行：

```powershell
& .\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.thread_benchmark `
  --source build/near_pml_order/pipeline_outer/scattering/report.json `
  --output build/thread_audit/new_run --threads 1 2 4 8 --repeats 2
```

线程测试的输出目录必须为空，防止混合不同批次。
成像检查入口为 `python -m fullwave_prototype.imaging_audit --help`。
需要已有 FEM 结果及可导入 pySCATMECH 的 NIST2 解释器。

成像校验命令如下；不要与计时测试同时运行：

```powershell
& .\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.imaging_audit `
  --source build/near_pml_order/pipeline_outer `
  --flat-source build/imaging_independent/flat128 `
  --reference-python D:/anaconda/envs/NIST2/python.exe `
  --output build/imaging_independent/new_run
```

缺少高采样平面参考时，先执行：

```powershell
& .\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run `
  --report build/imaging_independent/flat128_resources.json --max-memory-gib 12 `
  physical_audit --case flat --order 3 --threads 8 --static-condensation `
  --direct-solver sparsecholesky --pml-strength 1 --pml-mesh-scale 1.25 `
  --pupil-samples 128 --output build/imaging_independent/flat128
```

本轮通过完整验证的加速配置：

```powershell
& .\build\ngsolve_env\Scripts\python.exe -m fullwave_prototype.guarded_run `
  --report build/thread_audit/new_pipeline_resources.json --max-memory-gib 12 `
  imaging_run --output build/thread_audit/new_pipeline --accuracy balanced `
  --pml-order 2 --pml-inner-fraction .5 --threads 8 --pupil-samples 128 --near-diagnostics --archive-mode full
```

源码界面中的 `FEM threads` 可选择线程数，求解报告会记录实际设置。
原有通用默认 2 线程暂不自动调高：不同 PML 阶数、网格、颗粒和粗糙面
会改变矩阵分解的内存需求，不能把单一基准的最快设置强制用于所有工况。
已有旧 EXE 不会因修改源码自动更新；线程控制已经纳入后续的 2026.09.28-imaging-upgrade 分文件包。

## 限制

这是理想光学模型和选定工况的数值验证，不是真实镜头、相机或实验样品的标定。
有限显示区域不包含全部衍射尾部；全像面功率守恒测试不能替代逐点场检查。
网格、计算域、PML、粗糙度统计、有限光斑及像差等问题仍需单独收敛或实验验证。
线程提速也不改变这些物理误差。
