"""Collect all benchmark runs without cherry-picking the best agreement."""
import argparse
import json
from pathlib import Path


def summarize(directory):
    root = Path(directory)
    records = []
    lines = ['# Native FEM benchmark runs', '',
             'Experimental results, not a certification of near-field accuracy.', '',
             '| Run | Gap nm | Pol | Order | DOFs | FEM upper DSC um^2 | Reference um^2 | Difference |',
             '| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |']
    for path in sorted(root.glob('*/report.json')):
        record = json.loads(path.read_text(encoding='utf-8'))
        record['run'] = path.parent.name
        records.append(record)
        p, m = record['parameters'], record['metrics']
        reference = m.get('reference_hemisphere_dsc_um2')
        difference = m.get('hemisphere_relative_error')
        lines.append('| {} | {} | {} | {} | {} | {:.8g} | {} | {} |'.format(
            record['run'], p['gap_nm'], p['polarization'], p['order'], m['dofs'],
            m['hemisphere_dsc_um2'], '{:.8g}'.format(reference) if reference is not None else '-',
            '{:.3%}'.format(difference) if difference is not None else '-'))
    lines.extend(['', 'Mesh, PML and order can all affect these results. Read each report.json',
                  'for the complete configuration; rows with different gaps are different physics.',
                  'SCATMECH phase alignment, energy balance and contact-field convergence remain open.'])
    (root/'summary.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    (root/'SUMMARY.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', nargs='?', default='build/fem_validation')
    summarize(parser.parse_args().directory)
