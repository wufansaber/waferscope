"""Summarize original mixed-scene checks without promoting failed gates."""
import argparse
import json
from pathlib import Path

import numpy as np

from grating_scene.collection import collect
from grating_scene.config import SceneConfig
from grating_scene.pipeline import digest, write_json
from rcwa_prototype.comparison import field_change
from tools.rcwa_crosscheck_report import checked_load


def summarize(root):
    root = Path(root)
    rows, cases = [], {}
    for modes in (441, 585, 729):
        folder = root/('crosscheck_'+str(modes))
        report = json.loads((folder/'crosscheck.json').read_text(encoding='utf-8'))
        config = SceneConfig.from_dict(json.loads((folder/'config.json').read_text(encoding='utf-8')))
        reference_info, reference = checked_load(folder/'grcwa_reference')
        sample_info, sample = checked_load(folder/'grcwa_sample')
        if not sample_info['particle'] or reference_info['particle'] or config.scene.gap_um != 0:
            raise ValueError('Expected a genuine contact-particle plus grating comparison')
        far, _ = collect(reference, sample, config)
        delta = dict(sample, electric_reflected_z0=far['reflected_electric_delta'])
        cases[modes] = dict(config=config, delta=delta, masks=far['channel_masks'], report=report)
        power = []
        for a, b in zip(report['channels']['s4']['channels'], report['channels']['grcwa']['channels']):
            power.append(b['difference_field_efficiency']/a['difference_field_efficiency']-1)
        mirror = {backend: report[backend+'_diagnostics']['fields']['delta']['mirror_relative_l2']
                  for backend in ('s4', 'grcwa')}
        field = report['fields']['delta']['power_weighted_relative_complex_field_difference']
        rows.append(dict(modes=modes, harmonics=[config.scene.harmonics_x, config.scene.harmonics_y],
            slices=config.scene.slices, sample_complex_field_relative_error=
                report['fields']['sample']['power_weighted_relative_complex_field_difference'],
            delta_complex_field_relative_error=field, channel_delta_power_relative_errors=power,
            sample_near_xy=report.get('sample_near_xy'), delta_mirror_relative_errors=mirror,
            consistency_gates_passed=bool(field < .02 and max(abs(p) for p in power) < .05
                and all(value < .02 for value in mirror.values())),
            report_sha256=digest(folder/'crosscheck.json'), quantitative_hybrid_certified=False))
    convergence = []
    for old, new in ((441, 585), (585, 729)):
        a, b = cases[old], cases[new]
        channels = []
        for i in range(3):
            left, right = dict(a['delta'], propagating=a['masks'][i]), dict(b['delta'], propagating=b['masks'][i])
            metric = field_change(left, right, a['config'].scene.incidence_deg)
            pa, pb = [case['report']['channels']['grcwa']['channels'][i]['difference_field_efficiency']
                      for case in (a, b)]
            channels.append(dict(channel=i+1, complex_field_change=metric, power_change=pb/pa-1))
        convergence.append(dict(old_modes=old, new_modes=new, channels=channels,
            scope='Fixed cell and slices; My changes only. Field norm denominator is the refined grcwa field.'))
    blocked = root/'crosscheck_715.admission.json'
    return dict(status='mixed_scene_audited_not_quantitatively_certified', cases=rows,
        grcwa_y_convergence=convergence,
        blocked_715=json.loads(blocked.read_text(encoding='utf-8')) if blocked.is_file() else None,
        quantitative_hybrid_certified=False,
        conclusion='441-mode implementation and coincident near-field check agrees well. Higher-cutoff S4 '
            'mixed fields retain a mirror defect. Do not replace mixed validation with a bare-grating benchmark.',
        limitations=['S4/grcwa are two Fourier-modal implementations, not two independent physical methods.',
            'Fixed-cell agreement cannot certify isolated particles, finite grating ends, or real material inputs.',
            'No phase/amplitude fitting or post-hoc symmetry enforcement was used.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.root)
    write_json(args.root/'hybrid_summary.json', result)
    text = ['# 接触颗粒加光栅：单独验证结果', '', '日期：2026-10-28。此处验证的不是裸光栅。', '',
        '固定输入：266 nm、75 度 p 偏振、来自 -KY；500/250/100 nm 光栅；100 nm 接触球，gap=0；',
        '0.5x1 um 周期超胞，三路 NA=0.6、+45/-45/0 度。n/k 沿用受控比较输入，不是材料标定。', '',
        '| 模态 / 切片 | 样品复场差 | 颗粒差分复场差 | CH1 / CH2 / CH3 差场功率差 | S4 镜像复场误差 | 综合一致性门槛 |',
        '|---|---:|---:|---|---:|---|']
    for row in result['cases']:
        power = ' / '.join('{:+.6g}%'.format(100*v) for v in row['channel_delta_power_relative_errors'])
        text.append('| {} / {} | {:.6g}% | {:.6g}% | {} | {:.6g}% | {} |'.format(row['modes'], row['slices'],
            100*row['sample_complex_field_relative_error'], 100*row['delta_complex_field_relative_error'],
            power, 100*row['delta_mirror_relative_errors']['s4'], '通过这一实现核对' if row['consistency_gates_passed'] else '未通过'))
    near = result['cases'][0]['sample_near_xy']
    text += ['', '441 模态下，在 {} 个完全重合 XY 点上比较样品总 E，复场相对 L2 差为 {:.6g}；没有插值或幅相拟合。'.format(
        near['points'], near['relative_complex_e_l2']), '', '## grcwa 自身 Y 截断检查', '',
        '| 模态步 | CH1 / CH2 / CH3 复场变化 | CH1 / CH2 / CH3 差场功率变化 |', '|---|---|---|']
    for row in result['grcwa_y_convergence']:
        field = ' / '.join('{:.6g}%'.format(100*c['complex_field_change']['relative_power_weighted_complex_field_change'])
                           for c in row['channels'])
        power = ' / '.join('{:+.6g}%'.format(100*c['power_change']) for c in row['channels'])
        text.append('| {} -> {} | {} | {} |'.format(row['old_modes'], row['new_modes'], field, power))
    text += ['', '## 结论与限制', '',
        '- 裸光栅的 NIST RCW 对齐不能证明颗粒加光栅。此次额外计算的是完整周期混合样品及其无颗粒参考。',
        '- 一组良好实现一致性不等于数值收敛；较高基组的 S4 镜像复场仍有问题，未强行平均左右通道。',
        '- 715 模态、16 切片的独立对照被保留 3 GiB 系统余量的内存预检阻止，未运行、未记为通过。',
        '- 两种求解器都属于 Fourier modal/RCWA；仍缺少同问题有限边界 FEM 等独立方法认证。',
        '- 周期复制、球切片、联合截断、有限照明、实机镜组、材料标定仍限制定量结论。',
        '- 因此本版示例标记为实验性周期混合场景，不是已认证的孤立颗粒/有限光栅/实机 SNR 模型。', '',
        '同目录 raw/ 中保存各组原生输入、模态 NPZ、复场对比、近场及资源记录；源码工作区对应 crosscheck_* 目录。', '']
    (args.root/'HYBRID_VALIDATION.md').write_text('\n'.join(text), encoding='utf-8')


if __name__ == '__main__':
    main()
