# SCATMECH 原始数据读取规范

适用版本：WaferScope `2026.10.21-raw-data`；内嵌归档格式 `waferscope-raw-result`，版本 `1`。本说明仅涉及 SCATMECH 页面，不涉及 FEM。规则根据当前导出代码及打包版实际输出核对。

## 1. 文件与读取入口

| 计算类型 | 数值文件 | 配套说明 |
| --- | --- | --- |
| 原生远场、三路物镜收集 | `far_field_data.npz` | `simulation_parameters.json` |
| 相干成像及可选相机 | `coherent_imaging.npz` | `imaging_parameters.json`、启用时的 `camera_ch1.tiff` 等 |
| 粒径扫描 | `particle_size_sweep_data.npz` | `particle_size_sweep_parameters.json`、CSV |

新版各 NPZ 均有同名 `.manifest.json`，例如 `far_field_data.manifest.json`。PNG 只是显示结果，不适合作为原始数值读取入口。

同一个 `far_field_data.npz` 文件名可能表示原生远场，也可能表示物镜收集；不能仅凭文件名判断。重建后 `result['plot_kind'] == 'objective_collection'` 表示物镜，`coherent_imaging` 表示成像，`particle_size_sweep` 表示扫描。原生远场当前没有 `plot_kind`，可按是否存在原生 `mueller_raw`、`theta` 等字段识别。

## 2. 容器与跨语言数据类型

NPZ 是 ZIP 容器，每个数组是一个 `字段名.npy`。读取库须支持 NumPy NPY，而不能把解压后的文件当作纯二进制浮点数。

- NPY 头部给出 `descr`（dtype）、`shape` 和 `fortran_order`，必须读取这些信息，不硬编码维度、字节序或内存顺序。
- 常见 dtype：`<f8` 为小端 float64；`<c16` 为小端 complex128；`<u2` 为小端 uint16；`|b1` 为 bool。以实际头部为准。
- complex128 每个元素包含 float64 实部和虚部。保留两者，不要只读取实部、模值或强度。
- Python 数组通常采用 C 顺序，最后一个下标变化最快。MATLAB 等列优先语言必须正确处理维度顺序；不要为配合显示而直接转置物理数组。
- NPY Unicode 标量不是 UTF-8 字节数组；按 NPY 的 Unicode dtype 解码后，才能解析里面的 JSON。单独的 `.json` 文件则是 UTF-8。
- 不使用 pickle，不需要反序列化 Python 类。建议 Python 始终使用 `allow_pickle=False`。
- 数值无损压缩，保留复数、bool、NaN。孔径外 NaN 代表无效区域，不代表零信号；统计时同时使用掩膜。

数组尺寸随角分辨率、瞳面采样、FFT padding、相机 ROI 和粒径数量变化。本文的 `Nt`、`Np`、`N`、`L`、`Ny`、`Nx` 都是变量。

## 3. 内嵌元数据：必须优先采用的映射规则

NPZ 中保留键 `_archive_metadata_json`，其值是零维 Unicode 数组；读取其字符串并解析 JSON。其结构如下：

```json
{
  "format": "waferscope-raw-result",
  "version": 1,
  "result": {"kind": "dict", "items": {}},
  "source": null,
  "arrays": {"intensity": {"shape": [9, 31], "dtype": "<f8"}},
  "conventions": {}
}
```

上例只演示结构，不是某个完整结果。`result` 是计算结果的树，`source` 是可选的原生远场源树，`arrays` 是实际数值数组清单。参数、通道积分值等普通标量通常存储在树中，不一定有同名 NPY。

各节点的解码规则：

| `kind` | 解码方式 |
| --- | --- |
| `array` | 读取 `key` 指向的 NPZ 数组 |
| `dict` | 递归解码 `items` 的各个键值 |
| `list`、`tuple` | 按 `items` 原顺序递归解码；列表下标从 0 开始 |
| `scalar` | 直接使用 `value`，包括 JSON null、布尔值、数字和字符串 |
| `complex` | 由 `real`、`imag` 构造复数标量 |
| `nonfinite` | 将 `value` 的 `nan`、`inf`、`-inf` 还原为浮点特殊值 |

**不要依靠拆分 `__` 来重建结构。** 导出为了兼容旧字段和减少重复数组，会让不同逻辑字段引用同一个实际 NPZ 键。例如 CH1 数组可能仍存为 `intensity`，后两路存为 `channel_pupils__1__...`；源数组也可能引用 `far_x`，而不另存 `source__x`。唯一可靠映射是节点中的 `key`。

无元数据的历史文件，只能读取现有的平面数组和外部 JSON。新版读取函数不支持直接重建旧文件的嵌套结构，旧版未保存的数据也不能从图片补回。

## 4. 坐标、角度和单位

单位以后缀和字段定义为准，不能仅使用归档里概括性的 angles 提示：

- 网格 `theta`、`phi`、`global_phi`、物镜 `aolp` 使用弧度。
- 输入参数 `incidence_deg`、`tilt_deg`、`tilt_degrees` 使用度。`tilt_degrees` 虽不以 `_deg` 结尾，也仍然是度。
- 波长、粒径参数通常以 `_nm` 结尾；像面长度为 um；物镜出瞳长度为 mm。
- 原生远场 `x/y/z` 是方向余弦 `kx/k, ky/k, kz/k`，不是空间位置；物镜的 `x/y` 则是出瞳坐标 mm，不能混用。

晶圆在 XY 平面，+Z 为晶圆上方。照明来自 -Y 侧，在 YZ 入射面中沿 `(0, sin(i), -cos(i))` 传播；`i` 为相对法线的入射角。光源方向是传播方向的反向。

原生远场中 `phi == scatmech_phi`，是 SCATMECH 局部方位角。全局坐标相对它绕 +Z 旋转 90 度：

```text
global_phi = (phi + pi/2) mod (2*pi)
x = -sin(theta)*sin(phi)
y =  sin(theta)*cos(phi)
z =  cos(theta)
```

物镜通道中的 `phi = atan2(ky, kx)` 已是全局方位角，可能为负，不能再加 90 度。三路倾斜绕 +Y 轴，物镜轴为 `(sin(t), 0, cos(t))`。

二维图建议用真实坐标绘制，并设置 `origin='lower'`；数组第 0 轴是行/纵向，第 1 轴是列/横向，禁止从 PNG 的上下方向推断数据顺序。

## 5. 原生远场字段

令 `S = (Nt, Np)`。原生数组第 0 轴为 theta，第 1 轴为 SCATMECH phi。

| 字段 | 形状 | 含义/单位 |
| --- | --- | --- |
| `theta`, `phi`, `scatmech_phi`, `global_phi` | S | 角度网格，rad |
| `x`, `y`, `z` | S | 全局方向余弦，无量纲 |
| `intensity` | S | 绝对 DSC 或由 BRDF 转换的等效 DSC，um^2/sr |
| `intensity_db` | S | 相对当前峰值的 dB 显示值，非绝对强度 |
| `q`, `u`, `v` | S | 已乘转换系数的 Stokes Q/U/V，um^2/sr |
| `mueller_raw` | S + (4,4) | SCATMECH 原始 Mueller 矩阵，单位见 `mueller_raw_units` |
| `incident_stokes` | (4,) | 入射 Stokes，顺序 I/Q/U/V |
| `stokes_raw` | S + (4,) | 原始 Mueller × 入射 Stokes，尚未乘 DSC 转换系数 |
| `mueller_to_dsc_scale` | S | 每个角点的 DSC 转换系数 |
| `jones_dsc_um` | S + (2,2) | 存在时为复数 JonesDSC，振幅单位按 DSC 平方根约定记为 um |
| `dolp`, `docp` | S | 线偏振度、带符号圆偏振度，无量纲 |
| `e_x`, `e_y`, `e_z` | S | 从 Stokes 导出的线偏振方向，不是复电场 |

Jones 的末两轴依次为输出 s/p 分量、输入 s/p 分量；使用 SCATMECH `coords='psps'` 约定。Q/U/V 的正负含义沿用这个基底，不应另套相反的圆偏振符号约定。

若 `scatmech_output_kind == 'dsc'`，Mueller 原始单位是 um^2/sr，转换系数为 1。若为 BRDF 路径，原始单位是 1/sr，程序使用：

```text
scale = effective_area_um2 * abs(cos(incidence)) * abs(cos(theta))
stokes_raw = mueller_raw @ incident_stokes
intensity = max(stokes_raw[..., 0] * scale, 0)
```

Q/U/V 也乘同一 scale，但不作非负截断。直接采用已保存的 `mueller_to_dsc_scale`，不要再重复乘面积或余弦。颗粒加粗糙度的统计模型不保证存在可代表全部背景的相干 Jones 场；`jones_dsc_um` 缺失时不应自行从 Mueller 构造唯一相位。

`theta[:,0]` 和 `phi[0,:]` 可取为两个轴。方位角当前包含 0 和 2*pi 两个相同方向的端点；theta 通常止于 89 度，而非严格的 90 度，最后一段可能不是常步长。读取实际轴，不要用固定步长乘数组总和。

## 6. 三路物镜收集字段

重建 `result['channel_pupils']`，下标 0/1/2 对应 CH1/CH2/CH3。每一路的二维数组形状为 `(N,N)`，参数保存在 `result['parameters']`，原始远场及其材料参数在返回的 `source` 树中。

| 每路字段 | 含义/单位 |
| --- | --- |
| `x`, `y` | 物镜局部出瞳坐标，mm |
| `kx`, `ky`, `kz` | 全局方向余弦 |
| `local_kx`, `local_ky` | 沿物镜局部坐标的方向余弦 |
| `theta`, `phi` | 全局散射极角、方位角，rad |
| `sampled_dsc_um2_sr` | 从原始远场插值得到的 DSC，尚未按孔径裁切 |
| `aperture_mask` | bool，有效收光区域；NA=0 时全 False |
| `pupil_energy_density` | DSC / local_kz；孔径外为 NaN，按局部方向余弦平面积分可得截面 |
| `intensity`, `intensity_db` | 归一化/裁切的显示值，不用于绝对功率计算 |
| `q`, `u`, `dolp`, `aolp`, `pol_x`, `pol_y` | 采样 Stokes 与投影偏振显示；`aolp` 为 rad |
| `numerical_aperture`, `tilt_deg`, `pupil_radius_mm` | 各路 NA、倾角（度）和出瞳半径 |

顶层 `intensity`、`pupil_energy_density`、`x/y` 等只代表 CH1，顶层 `z` 甚至是显示强度别名而不是方向余弦。顶层 `q/u/v/e_x/e_y/e_z` 等兼容占位值不能当作物镜实际偏振。`far_intensity` 同样是归一化显示值；绝对原生远场应读取 `source['intensity']`。

`source_kind` 必须为 `scatmech` 才表示真实 SCATMECH 源；`synthetic` 是示意光场，不具有标定后的绝对 DSC 意义。

能量比例使用 `energy_channel_fractions`、`energy_collection_fraction`；几何体积比例使用 `geometric_channel_fractions`、`geometric_collection_fraction`，均是 0~1 而非百分数。三路总值来自孔径并集，重叠只计一次，不可把三路比例直接相加。

**重要实现约定：** 当前物镜 `collected_energy`、`hemisphere_energy`、`channel_collected_energies` 是 `sum(DSC*sin(theta))` 的离散加权和，未乘真实角度步长，不是 J、W，也不能直接当作 um^2。保存的能量比例是该离散求和的比值；重新按真实轴积分可能与它略有差异。`relative_collected_power` 是显示参考比，也不是绝对功率。

## 7. 相干成像与相机字段

`result['channels'][0..2]` 为三路独立成像；`parameters` 为材料/照明输入，`objective_parameters` 为物镜输入，`imaging_parameters` 为成像及相机输入。当前这一路径是颗粒相干成像，要求关闭粗糙度，不提供近场，也没有 FEM 那种 sample/reference 差分场。

| 每路字段 | 形状 | 含义/单位 |
| --- | --- | --- |
| `pupil_mask` | (N,N) | 有效孔径 |
| `pupil_axis_mm` | (N,) | 出瞳局部轴，mm |
| `jones_dsc_um` | (N,N,2,2) | 输出 s/p、输入 s/p 的复数 JonesDSC |
| `pupil_modes` | (N,N,2,2) | 出瞳横向 X/Y 分量、输入 s/p 模式 |
| `image_axis_um` | (L,) | 像面 X/Y 公用物理坐标轴，um |
| `image_e_modes`, `image_hz0_modes` | (L,L,3,2) | XYZ 分量、输入 s/p 模式的复场，通量归一化，非 SI V/m 或 A/m |
| `irradiance_w_um2` | (L,L) | 绝对像面强度，W/um^2 |
| `collected_dsc_um2`, `transmitted_dsc_um2`, `image_integral_dsc_um2` | 标量 | 收集、透过、像面积分等效截面，um^2 |
| `power_w` | 标量 | 整个计算像面窗口的功率，W |

`L = pupil_samples * padding`。像面坐标已经包含倍率，不要再次乘倍率；除以倍率得到的是各路局部物方平面的坐标，不自动等于晶圆 XY 坐标。两个输入模式的强度权重：s 光 `(1,0)`、p 光 `(0,1)`、非偏振光 `(0.5,0.5)`；非偏振输入不能先把两个复场相加再求模平方。像面强度按 E 与 H 的轴向能流计算，不建议简单用 `abs(E)**2` 替代。

启用相机时，每路增加 `detector`：

| 相机字段 | 形状/单位 | 说明 |
| --- | --- | --- |
| `pixel_x_um`, `pixel_y_um` | (Nx,)、(Ny,)，um | 像面像素中心 |
| `pixel_power_w` | (Ny,Nx)，W | 每像素有效感光面积内积分功率，不是强度密度 |
| `mean_photons` | (Ny,Nx)，光子数 | 曝光内均值 |
| `mean_signal_e`, `mean_background_e`, `mean_total_e` | (Ny,Nx)，电子数 | 信号、暗/背景、总电子数均值 |
| `measured_e` | (Ny,Nx)，电子数 | 一次噪声采样，含满阱截断及读出/固定偏移 |
| `digital_dn` | (Ny,Nx)，uint16 | 实际 ADC 灰度，0~2^bit_depth-1，无显示拉伸 |
| `nominal_dn` | (Ny,Nx)，float64 | 模型标称灰度，不是一次实际采样 |
| `digital_snr_valid`, `saturated`, `adc_saturated` 等 | (Ny,Nx)，bool | SNR 有效性和不同截断标记，不互相替代 |
| `ideal_snr`, `digital_ideal_snr`, `emva_total_snr` | (Ny,Nx)，无量纲 | 模型 SNR，不是从单张图实测得到的时域 SNR |

`camera_ch1.tiff` 等与 `digital_dn` 对应，用无损 16 位容器存储，即使有效位深是 8/10/12 位。读 TIFF 时不要自动缩放到 0~255。必须使用实际保存的 gain（e/DN）、offset（DN）、QE、曝光和截断标志，不能把灰度当作 W。

## 8. 粒径扫描字段

令 D 为粒径数量，`diameters_nm` 给出实际扫描点，不要只根据起止值推测点数。

| 字段 | 形状/单位 |
| --- | --- |
| `diameters_nm` | (D,)，nm |
| `channel_collected_dsc_um2` | (D,3)，um^2；列依次 CH1/CH2/CH3 |
| `hemisphere_dsc_um2` | (D,)，um^2；实际计算范围按保存的 theta 轴 |
| `channel_energy_fractions` | (D,3)，无量纲；有重叠时不能直接求和当作总比例 |
| `far_field_samples` | D 个原生远场结果，保留各粒径材料参数及 Mueller 等数组 |
| `sampled_dsc_um2_sr` | (D,Nt2,Np2)，收集积分使用的插值 DSC |
| `collection_sampling` | 共同的 theta/phi 轴、方向网格、sin_theta 和三路 bool 掩膜 |

扫描的积分采用保存轴上的梯形积分；它与第 6 节物镜显示结果的历史离散加权和不是同一种数值实现，不应混用两个“energy”字段。

## 9. 独立 Python 读取示例

只需 NumPy 和随说明附带的 `result_archive.py`，不需要 Qt、pySCATMECH 或原软件安装。使用普通 Python 的文件路径语法，路径按实际位置替换。

也附带 `read_scatmech_data.py`，可先核对文件哈希并导出字段清单：

```powershell
python read_scatmech_data.py "D:\results\far_field_data.npz" --inventory "D:\results\fields.json"
```

脚本与 `result_archive.py` 放在同一目录即可。无配套 manifest 时返回 `manifest_verified: null`，不是校验通过。清单只含实际数组的名字、形状和 dtype；完整层级及标量仍以 NPZ 内嵌元数据为准。

```python
from pathlib import Path
import numpy as np
from result_archive import load_result_archive

path = Path(r'D:\results\far_field_data.npz')
result, source = load_result_archive(path)

if result.get('plot_kind') == 'objective_collection':
    for index, channel in enumerate(result['channel_pupils'], 1):
        mask = channel['aperture_mask']
        dsc = np.where(mask, channel['sampled_dsc_um2_sr'], np.nan)
        print(index, dsc.shape, channel['numerical_aperture'])
    native = source
else:
    native = result

if native is not None and 'mueller_raw' in native:
    print(native['parameters'])
    theta_axis = native['theta'][:, 0]
    phi_axis = native['phi'][0, :]
    dsc = native['intensity']
    integrate = getattr(np, 'trapezoid', None) or np.trapz
    cross_section_um2 = integrate(
        integrate(dsc * np.sin(theta_axis[:, None]), x=phi_axis, axis=1),
        x=theta_axis, axis=0)
    print('integrated DSC (um^2):', cross_section_um2)
```

这里的积分仅覆盖已采样的 theta 范围，并非补全了未采样的最后 1 度。DSC 是截面，不是功率；只有乘以定义一致的入射辐照度 W/um^2，才得到 W。入射到晶圆平面的平均辐照度与光束法向辐照度在斜入射时不同，不要额外乱乘或漏乘 cos(i)。

只需要平面数组或读取旧版时：

```python
import numpy as np
with np.load(r'D:\results\far_field_data.npz', allow_pickle=False) as data:
    for name in data.files:
        array = data[name]
        print(name, array.shape, array.dtype)
```

## 10. 完整性与第三方程序验收

1. 对整个 NPZ 文件计算 SHA256，与同名 `.manifest.json` 的 `sha256` 比较。哈希是对文件字节计算，不是对解压后的数组求和。
2. 验证内嵌 `format` 和 `version`；不认识的版本应拒绝或提示，不要猜字段。
3. 数组清单应等于 `metadata['arrays']` 加元数据键；检查每个实际 shape/dtype。数值 NaN 本身不代表文件损坏。
4. 按树节点的 `key` 解析别名，三路顺序与原始参数一起保存。不要从标题“左/右”推导通道顺序。
5. 检查 `stokes_raw == mueller_raw @ incident_stokes`，以及第 5 节 DSC 转换关系。允许正常浮点舍入误差。
6. 成像模式可检查 `sum(irradiance_w_um2)*dx*dy` 与 `power_w` 一致；相机 ROI 功率可以小于整幅像面功率，不应强制相等。
7. 第三方绘图保持实际坐标、共同色标和原始绝对值；若另作归一化，单独输出新字段，绝不覆盖原数组。

`load_result_archive` 会验证归档内数组 shape/dtype 并还原结构，但当前不会自动核对外部 manifest 的 SHA256；需要独立核对，示例脚本已实现该步骤。归档没有逐数组的单位字典，单位应按本文及配套参数解释，不能仅凭 dtype 或名字猜测。
