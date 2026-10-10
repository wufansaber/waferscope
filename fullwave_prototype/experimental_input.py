"""Prepare and audit experimental height maps or two-sided PSDs."""
import argparse
import json
from pathlib import Path

import numpy as np

from .experimental_surface import height_periodogram, load_height, load_psd


def prepare(source, output, kind='height', detrend='plane', missing='reject'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    height = None
    if kind == 'height':
        height, report = load_height(source, detrend, missing)
        fx, fy, density = height_periodogram(height, report['dx_um'], report['dy_um'])
        np.savez_compressed(output/'processed_height.npz', height_nm=height,
                            dx_um=report['dx_um'], dy_um=report['dy_um'])
    elif kind == 'psd':
        fx, fy, density, report = load_psd(source)
    else:
        raise ValueError('Measurement kind must be height or psd')
    variance = float(density.sum()*(fx[1]-fx[0])*(fy[1]-fy[0]))
    report['psd_band_rms_nm'] = float(np.sqrt(variance)*1000)
    if height is not None:
        report['parseval_relative_error'] = float(abs(variance-
            np.mean((height/1000)**2))/max(variance, 1e-300))
    np.savez_compressed(output/'measured_psd.npz', fx_per_um=fx, fy_per_um=fy, psd_um4=density)
    for name in ('measured_psd.json', 'processed_height.json'):
        if name.startswith('processed') and height is None:
            continue
        (output/name).write_text(json.dumps(dict(provenance=report), indent=2), encoding='utf-8')
    (output/'measurement_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    fig, axes = plt.subplots(1, 2 if height is not None else 1, figsize=(10 if height is not None else 6, 4.5),
                             layout='constrained', squeeze=False)
    ax = axes[0, -1]
    maximum = max(float(density.max()), 1e-30)
    image = ax.pcolormesh(fx, fy, np.maximum(density, maximum*1e-8), shading='auto',
        norm=LogNorm(vmin=maximum*1e-8, vmax=maximum), cmap='viridis')
    ax.set(xlabel='fx (cycles/um)', ylabel='fy (cycles/um)', title='Measured-band 2D PSD')
    fig.colorbar(image, ax=ax, label='PSD (um4)')
    if height is not None:
        image = axes[0, 0].imshow(height, origin='lower', cmap='coolwarm',
            extent=(0, report['footprint_xy_um'][0], 0, report['footprint_xy_um'][1]))
        axes[0, 0].set(xlabel='X (um)', ylabel='Y (um)', title='Processed measured height')
        fig.colorbar(image, ax=axes[0, 0], label='Height (nm)')
    fig.savefig(output/'measurement.png', dpi=150)
    plt.close(fig)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--kind', choices=('height', 'psd'), default='height')
    parser.add_argument('--detrend', choices=('none', 'mean', 'plane'), default='plane')
    parser.add_argument('--missing-policy', choices=('reject', 'nearest'), default='reject')
    args = parser.parse_args()
    print(json.dumps(prepare(args.input, args.output, args.kind, args.detrend, args.missing_policy), indent=2))


if __name__ == '__main__':
    main()
