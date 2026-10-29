"""Write the quantitative audit note from actual saved results, without fitted scaling."""
from pathlib import Path


def phase_section(rows):
    if not rows:
        return []
    lines = ['## 相位原点独立对照', '',
        '以下复光瞳相对 L2 差异使用同一 FEM 场。只消除一个不影响强度的整体常数相位，不能拟合幅值、位置或角度相关相位。', '',
        '| 直径 nm | 通道 | 转换前复场差异 % | 晶圆面相位原点下差异 % | 转换引起的光瞳功率变化 % |',
        '|---:|---:|---:|---:|---:|']
    for row in rows:
        lines.append('| {} | CH{} | {:.6g} | {:.6g} | {:.6g} |'.format(
            row['diameter_nm'], row['channel'], row['old_relative_l2']*100,
            row['wafer_reference_relative_l2']*100, row['pupil_power_relative_change']*100))
    return lines+['', '相位转换不改变收集功率。100 nm 下剩余复场差异不能靠相位转换消除，仍须检查 FEM 的离散与边界收敛。', '']


def mesh_section(rows):
    if not rows:
        return []
    lines = ['## 同规则网格细化', '',
        '固定 100 nm 接触球、0° P、半径 0.22 μm、二阶物理/PML、PML 幂次 1/厚度 0.16 μm。每个细化序列内的 PML 强度不变。',
        '空气/基底/PML 按 mesh_scale 同时细化，颗粒网格固定 0.01 μm。表中是相对前一级 FEM 的变化，不是相对 NIST 的误差。', '',
        '| 新一级目录 | 基底网格 μm | PML 强度 | 通道 | 功率变化 % | 复光瞳 L2 差异 % | 单步满足功率 1% / 场 2% 目标 |',
        '|---|---:|---:|---:|---:|---:|---|']
    for row in rows:
        lines.append('| {} | {:g} | {:g} | CH{} | {:+.5g} | {:.5g} | {} |'.format(
            row['current'], row['substrate_mesh_um'], row['pml_strength'], row['channel'],
            row['relative_power_change']*100, row['pupil_relative_l2']*100,
            '是' if row['within_step_targets'] else '否'))
    return lines+['', '单步通过仍不等于全部收敛，尚需连续两级及独立计算域/PML 检查；复场只消除整体常数相位。', '']


def symmetry_section(rows):
    if not rows:
        return []
    lines = ['## 边路对称性', '',
        '使用同一绝对量的相对差 `2|CH1-CH2|/(CH1+CH2)`，不分别归一化图像。下表为 10～100 nm 中的最大值。', '',
        '| 入射角 / 偏振 | 球冠积分功率差异 % | 离散光瞳功率差异 % | 相机期望净峰值差异 % |',
        '|---|---:|---:|---:|']
    for row in rows:
        lines.append('| {:g}° / {} | {:.6g} | {:.6g} | {:.6g} |'.format(
            row['incidence_deg'], row['polarization'].upper(),
            row['max_relative_cap_dsc_um2_asymmetry']*100,
            row['max_relative_pupil_dsc_um2_asymmetry']*100,
            row['max_relative_expected_peak_signal_dn_asymmetry']*100))
    return lines+['', '总功率与单像元峰值不能混用。峰值还受有限像面采样、像素积分和亚像素位置影响；随机 raw DN 帧会另外包含独立噪声。这是平整晶圆结果，不证明每张随机粗糙面的两路也相同。', '']


def write_note(base, nist, verification, nist_name):
    rows = nist['rows']
    c = nist['controls']
    p, s, detector = c['parameters'], c['settings'], c['settings']['detector']
    lines = ['# 266 nm 接触球粒径与亮度验证', '', '数据采集：2026-10-26；验证报告更新：2026-10-27。', '',
        '## 结论边界', '',
        '球底始终接触平整晶圆：gap=delta=0，球心高度随尺寸变化为 d/2。未使用固定球心或悬空间隙。',
        'SCATMECH 是理想球/平面基底的独立散射理论参考，不是任意缺陷、粗糙面或实际相机的无条件金标准。',
        '验证声明限定为理想球形 PSL/平整硅基底；材料未标定、无有效批次及实验对照时，不能称为真实 PSL 标准颗粒的已认证结果。详见 [PSL 验证范围与参考文献](../../docs/PSL_SPHERE_VALIDATION_CN.md)。',
        '复场已经从球心相位原点转换到晶圆顶面；不拟合或归一化强度。原生 JonesDSC 与施加的相位分别保存。',
        '无 haze、无镜面漏光、无 TDI；正常入射的中路会在真实系统收进镜面反射，本报告未模拟它，不能将该组 SNR 当作真实暗场探测性能。', '',
        '追溯用历史记录为 line_x、200×50×50 nm、0°入射、gap=10 nm，并不是一个 100 nm 接触球；本轮单独重建球体，未将该旧结果改名当作球形数据。', '',
        '## 固定输入', '',
        '- λ=266 nm；球直径 10～100 nm，步长 10 nm；接触平整晶圆。',
        '- 颗粒 n/k={:g}/{:g}；晶圆 n/k={:g}/{:g}。这些是继承的受控输入，不代表已标定的 266 nm PSL/Si 色散。'.format(
            p['particle_n'], p['particle_k'], p['substrate_n'], p['substrate_k']),
        '- 三路 NA={}，倾角={}°；{} 倍，物镜焦距 {} mm。'.format(
            s['numerical_apertures'], s['tilt_degrees'], s['magnification'], s['focal_length_mm']),
        '- 平面波入射辐照度保持 {:.6g} W/μm²；透过率 {:.6g}。两入射角控制组固定辐照度，不声称两者的晶圆面 2×2 μm 光斑均为 1 W。'.format(
            s['irradiance_w_um2'], s['transmission']),
        '- 相机 {}×{} 像元，像元 {} μm，曝光 {} ms，QE={}，增益 {} e-/DN，黑电平 {} DN，读噪声 {} e-。'.format(
            detector['pixels_x'], detector['pixels_y'], detector['pixel_pitch_um'], detector['exposure_ms'],
            detector['quantum_efficiency'], detector['gain_e_dn'], detector['offset_dn'], detector['read_noise_e']),
        '- 每种输入保存 100 对模拟 sample/reference 原始灰度帧，三路固定像元响应与噪声种子均记录。', '',
        '## 相机净信号峰值', '',
        '以下是无噪声期望信号的像素积分峰值 DN，不含黑电平；不是在噪声帧中寻找的最大灰度。所有尺寸使用同一功率、曝光和增益，没有逐图归一化。', '']
    for incidence in (0., 75.):
        lines += ['### P 偏振，入射角 {:g}°'.format(incidence), '',
                  '| 直径 nm | CH1 净峰值 DN | CH2 净峰值 DN | CH3 净峰值 DN |', '|---:|---:|---:|---:|']
        for diameter in range(10, 101, 10):
            selected = [next(r for r in rows if r['diameter_nm'] == diameter and r['incidence_deg'] == incidence
                             and r['polarization'] == 'p' and r['channel'] == ch) for ch in (1, 2, 3)]
            lines.append('| {} | {} |'.format(diameter, ' | '.join('{:.6g}'.format(r['expected_peak_signal_dn']) for r in selected)))
        lines.append('')
    monotonic = all(r['monotonic_increasing'] for r in verification['size_statistics'])
    lines += ['## 粒径规律与收敛', '',
        '- 0°/75°、P/S、三路收集功率均随尺寸增大而单调上升：{}。'.format('是' if monotonic else '否，见 verification.json'),
        '- d^6 是自由空间 Rayleigh 渐近规律，不能强制应用到整个 10～100 nm 球/晶圆体系。100 nm/266 nm 的尺寸参数约为 1.18；球心高度、局域反射场和多极响应均随粒径变化。',
        '- 中路并非普遍最亮；应按入射角、偏振及收集功率/像素峰值分别判断。',
        '- 计算相机期望峰值与实际随机帧 SNR 分开。微弱尺寸的 raw DN 主要是黑电平和噪声，不能据一个亮像素判断检测到颗粒。']
    conv = verification['convergence']
    for label, chosen, key in (
            ('多极阶数 12→16，积分 DSC 最大变化', [r for r in conv if r['lmax'] == 16], 'cap_change'),
            ('球冠积分阶数 24→32 最大变化', [r for r in conv if r['cap_order'] == 32], 'cap_change'),
            ('光瞳 64→96，相机净峰值最大变化', [r for r in conv if r['pupil_samples'] == 96], 'camera_peak_change')):
        if chosen:
            lines.append('- {}：{:.6g}%。仅覆盖已测试的 10/100 nm 端点。'.format(label, max(abs(r[key]) for r in chosen)*100))
    lines += ['- 光瞳积分与独立球冠积分最大差异 {:.5g}%；像面/光瞳 Parseval 通量残差 {:.5g}。后者不是整个晶圆/颗粒系统的入射-散射-吸收能量闭合证明。'.format(
        verification['max_pupil_vs_cap_relative_change']*100, verification['max_flux_relative_error']), '']
    lines += phase_section(verification.get('phase_origin_comparisons', []))
    lines += symmetry_section(verification.get('side_symmetry', []))
    lines += [
        '## FEM 对照', '', '| 算例 | 粒径 nm / 入射角° | 三路功率相对差异 % | 三路相机净峰值相对差异 % |', '|---|---:|---|---|']
    groups = {}
    for row in verification['fem_comparison']:
        groups.setdefault(row['source'], []).append(row)
    for source, group in groups.items():
        name = '/'.join(Path(source).parts[-3:-1])
        lines.append('| {} | {:g} / {:g} | {} | {} |'.format(name, group[0]['diameter_nm'], group[0]['incidence_deg'],
            ', '.join('{:+.4f}'.format(r['relative_power_change']*100) for r in group),
            ', '.join('{:+.4f}'.format(r['relative_peak_change']*100) for r in group)))
    lines += ['', '只列出完成并保存成像报告的 FEM；触发内存保护的尝试不填 0、不当作有效结果。',
        '最初约 22% 的偏差出现在这轮内存受限的一阶 PML 算例，不是所有 UI 默认配置的统一误差，也不能将偏差全部归因于 PML 阶数。生产界面的 PML order=0 表示跟随物理区域阶数。',
        '不同阶数、范围、网格、PML 改动必须结合看，不能从某一组最接近 NIST 的数字宣布收敛。绝不使用乘法系数强行对齐。', '',
    ]
    lines += mesh_section(verification.get('mesh_convergence', []))
    lines += [
        '## 结果文件', '',
        '- [亮度曲线](brightness_curves.png)：0°及75°入射的绝对收集功率、相机净峰值和局部幂指数。',
        '- [中路/边路比值](channel_ratios.png)：P/S 及两入射角分开。',
        '- [模拟原始灰度对比](camera_size_montage.png)：0° P、六个代表尺寸、共同显示窗口。',
        '- [全部尺寸数值]('+nist_name+'/brightness.csv)：含功率、期望 DN、原始灰度 ROI SNR；ROI 是完整相机窗口。',
        '- [数值校验](verification.json)。',
        '- [数值收敛与求解器审查](NUMERICAL_REVIEW_CN.md)：列明未通过的验证和下一阶段验收条件。',
        '- 原始复 Jones、光瞳、像面 E/HZ0、角积分采样、相机期望电子数及 100 对 raw DN 在 `'+nist_name+'/cases/*/raw_fields.npz`。E/HZ0 是通量归一化场，不直接等于 SI V/m。',
        '- 历史输入与校验摘要保留 SHA-256；早期 gap=10 nm 和未转换相位的扫描未覆盖，用于追溯，不用于当前结论。', '',
        '参考：[NIST 模型定义](https://pages.nist.gov/SCATMECH/docs/bobvlieg.htm)。', '']
    (base/'BRIGHTNESS_VALIDATION_CN.md').write_text('\n'.join(lines), encoding='utf-8')
