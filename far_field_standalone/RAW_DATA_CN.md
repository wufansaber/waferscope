# 原始数据保存与迁移

本版 UI 默认完整保存数值结果，不只保存图片或归一化显示数据。

## SCATMECH 导出

- `far_field_data.npz`：远场网格、原始 Mueller/Stokes、DSC、偏振、输入参数；物镜计算保留三路瞳面、孔径掩膜及原始远场源。
- `coherent_imaging.npz`：各通道 Jones、瞳面、像面复场、强度及启用相机时的数值结果。
- `particle_size_sweep_data.npz`：积分曲线及每个粒径的原始远场、采样 DSC 与收集掩膜。
- 同名 `.manifest.json`：NPZ 的 SHA256、数组尺寸、数据类型及约定。原有 JSON、CSV、TIFF 和 PNG 仍保留。

读取所有版本的普通数值数组：

```python
import numpy as np
with np.load(r'D:\results\far_field_data.npz', allow_pickle=False) as data:
    print(data.files)
    for name in data.files:
        print(name, data[name].shape, data[name].dtype)
```

新版嵌套结构可使用随包 `docs/result_archive.py` 中的读取函数；在该模块可导入的 Python 环境运行：

```python
from result_archive import load_result_archive
result, source = load_result_archive(r'D:\results\far_field_data.npz')
for channel in result.get('channel_pupils', []):
    print(channel['sampled_dsc_um2_sr'], channel['aperture_mask'])
if source is not None:
    print(source['intensity'])
```

原生远场 `intensity` 单位为 DSC（um^2/sr）；物镜结果的顶层 `intensity` 是 CH1 的显示别名，不能作为原始绝对强度。物镜原始值使用各路 `sampled_dsc_um2_sr` 和 `pupil_energy_density`，后者包含出瞳投影 Jacobian；原始远场从 `source` 读取。角度默认弧度，明确以 `_deg` 命名的量除外。复数相位和孔径外 NaN 保留。`e_x/e_y/e_z` 是偏振方向，不是复电场。

## FEM 完整结果

UI 使用 `full` 存档。`image_fields.npz` 保存像面复场、强度及启用时的相机数值；`scattering/paired_fields.npz` 保存近场、远场和可复用体源，`report.json` 保存散射输入。重新成像会将原始散射源复制到新结果目录。混合模型同时保留共享表面、局部 FEM 和外部 SPM 所需的源文件。

迁移时复制整个结果目录，不要只复制 PNG 或某一个 NPZ。压缩为无损压缩，不压缩有效数字；文件较大是正常现象。旧版曾省略的原始量不能从图片恢复，需要重新计算。CLI 显式 `--archive-mode compact` 仍为旧脚本保留，但不等价于完整存档。此处原始数据是已输出的物理场与输入源，不包括求解器 LU 等内部工作缓存。
