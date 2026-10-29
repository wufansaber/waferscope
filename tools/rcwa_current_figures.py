"""Plot preserved frozen-current data; no smoothing, rescaling, or fresh solves."""
import argparse
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def figures(base, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    refined = base/'radiation_lx0500_refined'
    report = read(refined/'report.json')
    lattice = read(base/'frozen_lattice'/'report.json')
    output.mkdir(parents=True, exist_ok=True)
    final = report['levels'][-1]
    with np.load(refined/final['raw'], allow_pickle=False) as raw:
        data = {key: raw[key] for key in raw.files}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), layout='constrained')
    x = np.arange(3)
    values = [row['cross_section_um2'] for row in lattice['lattices']]+[final['cross_section_um2']]
    labels = ['Lattice Lx=0.500 um', 'Lattice Lx=0.667 um', 'Continuous cone integral']
    for offset, row, label, color in zip((-.25, 0., .25), values, labels, ('#4878A8', '#CF8050', '#3A9275')):
        axes[0].bar(x+offset, row, .23, label=label, color=color)
    axes[0].set(xticks=x, xticklabels=['CH1 (+45 deg)', 'CH2 (-45 deg)', 'CH3 (0 deg)'],
                ylabel='Collected cross section (um^2)', title='Identical frozen current, different angular sampling')
    axes[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    axes[0].legend(frameon=False, fontsize=8)
    levels = read(base/'radiation_lx0500'/'report.json')['levels']+report['levels']
    target = np.asarray(final['cross_section_um2'])
    for index, (label, color) in enumerate(zip(('CH1', 'CH2', 'CH3'), ('#4878A8', '#CF8050', '#3A9275'))):
        axes[1].plot(np.arange(len(levels)), [100*(row['cross_section_um2'][index]/target[index]-1)
                                            for row in levels], 'o-', label=label, color=color)
    axes[1].set(xticks=np.arange(len(levels)),
                xticklabels=['{}x{}\nref {}'.format(row['radial'], row['azimuthal'], row['reference_halfwidth'])
                             for row in levels], ylabel='Change relative to finest diagnostic (%)',
                title='Angular and reference-kernel refinement')
    axes[1].axhline(0., color='0.5', lw=.8)
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.savefig(output/'sampling_comparison.png', dpi=180)
    plt.close(fig)
    radial, azimuthal = final['radial'], final['azimuthal']
    mu, _ = np.polynomial.legendre.leggauss(radial)
    lower = np.sqrt(1-report['na']**2)
    mu = lower+(mu+1)*(1-lower)/2
    # Edges only frame the unsmoothed Gauss-node samples; integrals use saved weights.
    edges = np.r_[lower, .5*(mu[1:]+mu[:-1]), 1.]
    radii = np.sqrt(np.maximum(0., 1-edges**2))[::-1]
    phi = np.linspace(0, 2*np.pi, azimuthal+1)
    xx, yy = np.outer(radii, np.cos(phi)), np.outer(radii, np.sin(phi))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), layout='constrained')
    norm = Normalize(0., float(np.max(data['dsc_um2_per_sr'])))
    for index, ax in enumerate(axes):
        values = data['dsc_um2_per_sr'][index].reshape(radial, azimuthal)[::-1]
        mesh = ax.pcolormesh(xx, yy, values, shading='flat', norm=norm, cmap='viridis', rasterized=True)
        ax.set(aspect='equal', xlim=(-.62, .62), ylim=(-.62, .62),
               xlabel='Local transverse direction u', ylabel='Local transverse direction v',
               title='CH{} | tilt {:+g} deg\nIntegral = {:.6g} um^2'.format(
                   index+1, report['tilts_deg'][index], final['cross_section_um2'][index]))
    fig.colorbar(mesh, ax=axes, label='Frozen-current DSC (um^2 / sr)', shrink=.8)
    fig.suptitle('NA=0.6 | one shared absolute color scale | not an isolated-defect solution')
    fig.savefig(output/'aperture_radiation.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    figures(args.base, args.output)


if __name__ == '__main__':
    main()
