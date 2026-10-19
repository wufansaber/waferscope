"""Generate reproducible brochure figures from real solvers and saved fields."""
import csv
import hashlib
import json
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/showcase_20261017'
ASSETS = OUT / 'assets'
LOW = ROOT / 'build/channel_symmetry_20261016/low_haze_02nm'
SMALL = ROOT / 'build/fullwave_ui/20261017_164144_979822'
LARGE = ROOT / 'build/cross_2um_20261017/background_refined'
COLORS = ['#007f86', '#bb3971', '#c58c16']


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(fig, name):
    fig.savefig(ASSETS / (name+'.png'), dpi=190, facecolor='white')
    plt.close(fig)


def disk(ax, theta, phi, values, title, vmax=None):
    # The saved FEM azimuth is already global; do not rotate it a second time.
    def edges(axis):
        axis = np.asarray(axis)
        return np.r_[axis[0]-(axis[1]-axis[0])/2, (axis[:-1]+axis[1:])/2,
                     axis[-1]+(axis[-1]-axis[-2])/2]
    ph, th = np.meshgrid(edges(np.unwrap(phi)), np.clip(edges(theta), 0, np.pi/2))
    art = ax.pcolormesh(np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), values,
        shading='flat', cmap='viridis', vmin=0, vmax=vmax, rasterized=True)
    ax.add_patch(plt.Circle((0, 0), 1, fill=False, color='#333333', lw=.6))
    ax.set(xlim=(-1.03, 1.03), ylim=(-1.03, 1.03), aspect='equal', title=title,
           xlabel='kx/k', ylabel='ky/k')
    return art


def load_dense_far(asset_dir=ASSETS, source=LOW):
    report = read(asset_dir/'case01_dense_far_report.json')
    signature = hashlib.sha256((source/'hybrid_source_report.json').read_bytes()).hexdigest()
    if report['source_metadata_sha256'] != signature:
        raise ValueError('Dense far field belongs to a different physical source; re-evaluate it')
    with np.load(asset_dir/'case01_dense_far.npz', allow_pickle=False) as stored:
        axis, values, mask = (stored[k] for k in ('axis', 'dsc_total', 'mask'))
    if (len(axis) < 401 or values.shape != (len(axis), len(axis)) or mask.shape != values.shape
            or not np.all(np.diff(axis)>0) or not np.isfinite(values[mask]).all()
            or np.any(values[mask]<0) or np.isfinite(values[~mask]).any()):
        raise ValueError('Invalid dense physical far-field array')
    return axis, values, report


def cartesian_disk(ax, axis, values, title, vmax=None):
    art = ax.pcolormesh(axis, axis, values, shading='auto', cmap='viridis', vmin=0, vmax=vmax,
                       rasterized=True)
    ax.add_patch(plt.Circle((0,0), 1, fill=False, color='#333333', lw=.7))
    ax.set(xlim=(-1.03,1.03), ylim=(-1.03,1.03), aspect='equal', title=title,
           xlabel='kx/k', ylabel='ky/k')
    return art


def workflow_case():
    dense_axis, dense_values, dense_report = load_dense_far()
    with np.load(LOW/'local_scattering/paired_fields.npz') as near, np.load(LOW/'hybrid_far_fields.npz') as far, np.load(LOW/'image_fields.npz') as images:
        fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.1), constrained_layout=True)
        for ax, plane, labels in zip(axes[0], ('xy', 'xz', 'yz'), (('x', 'y'), ('x', 'z'), ('y', 'z'))):
            e = near['near_'+plane+'_e_with']
            intensity = np.sum(abs(e)**2, axis=-1)
            art = ax.pcolormesh(near['near_'+plane+'_u_um'], near['near_'+plane+'_v_um'], intensity,
                                shading='auto', cmap='magma')
            ax.set(title='Local FEM | '+plane.upper(), xlabel=labels[0]+' (um)', ylabel=labels[1]+' (um)', aspect='equal')
            fig.colorbar(art, ax=ax, label='|E/E_inc|^2', shrink=.75)
        art = cartesian_disk(axes[1, 0], dense_axis, dense_values, 'Hybrid far field')
        fig.colorbar(art, ax=axes[1, 0], label='DSC (um^2/sr)', shrink=.75)
        vmax = max(np.max(images['ch%d_image_with_w_um2'%i]) for i in (1, 2))
        for ax, i in zip(axes[1, 1:], (1, 2)):
            axis = images['ch%d_image_axis_um'%i]
            art = ax.pcolormesh(axis, axis, images['ch%d_image_with_w_um2'%i], shading='auto', cmap='magma', vmin=0, vmax=vmax)
            ax.set(title='Image CH%d | 20x'%i, xlim=(-30, 30), ylim=(-30, 30), aspect='equal', xlabel='Image x (um)', ylabel='Image y (um)')
            fig.colorbar(art, ax=ax, label='W/um^2', shrink=.75)
        save(fig, 'case01_pipeline')
        fig, ax = plt.subplots(figsize=(5.6,5.0), constrained_layout=True)
        art = cartesian_disk(ax, dense_axis, dense_values, 'Hybrid far field | dense angular evaluation')
        ax.plot(0,-np.sin(np.deg2rad(75)), marker='*', color='white', markeredgecolor='#333333', markersize=9)
        ax.annotate('Laser from -KY\n75 deg incidence', xy=(0,-.9659), xytext=(-.87,-.66),
                    fontsize=9, color='#202020', arrowprops={'arrowstyle':'->','color':'#333333'},
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':3})
        fig.colorbar(art, ax=ax, label='DSC (um^2/sr)', shrink=.85)
        save(fig, 'case01_far_field_dense')
        fig, axes = plt.subplots(1,2,figsize=(10,4.6),constrained_layout=True)
        vmax = max(float(np.nanmax(dense_values)),float(np.max(far['dsc_total'])))
        disk(axes[0],far['theta'],far['phi'],far['dsc_total'],'Before | 5 deg angular grid',vmax)
        art = cartesian_disk(axes[1],dense_axis,dense_values,'After | 401 x 401 direction disk',vmax)
        fig.colorbar(art,ax=list(axes),label='DSC (um^2/sr) | shared absolute scale',shrink=.8)
        save(fig,'case01_sampling_comparison')
    return dict(source=str(LOW.relative_to(ROOT)), physical=read(LOW/'imaging_report.json')['scattering_parameters'],
                far_field_angular_evaluation=dense_report,
                scope='Local FEM near fields; coherent hybrid far/image fields, not a global full-wave wafer solve.')


def sweep_and_apertures():
    from scatmech_ui.parameters import SimulationParameters, ObjectiveCollectionParameters
    from scatmech_ui.native_far_field import calculate_far_field
    from scatmech_ui.size_sweep import _collection_sampling_grid, _integrate_dsc
    from scatmech_ui.objective_math import sample_far_field, objective_far_field_view
    p = SimulationParameters(wavelength_nm=266, angular_step_deg=4, illumination_length_um=2, illumination_width_um=2)
    objective = ObjectiveCollectionParameters(channel_numerical_apertures=(.4, .4, .4))
    sampling = _collection_sampling_grid(objective)
    diameters = np.arange(30, 201, 10)
    rows, reference = [], None
    for d in diameters:
        result = calculate_far_field(replace(p, particle_diameter_nm=float(d)))
        intensity = sample_far_field(result, sampling['kx'], sampling['ky'], sampling['kz'])['intensity']
        powers = [_integrate_dsc(intensity, mask, sampling) for mask in sampling['channel_masks']]
        rows.append([int(d)]+powers)
        if d == 100:
            reference = result
        print('SCATMECH diameter %d nm completed'%d, flush=True)
    with (OUT/'particle_size_sweep.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['diameter_nm', 'ch1_collected_dsc_um2', 'ch2_collected_dsc_um2', 'ch3_collected_dsc_um2'])
        writer.writerows(rows)
    a = np.array(rows)
    fig, ax = plt.subplots(figsize=(10.6, 4.6), constrained_layout=True)
    for i in range(3):
        ax.semilogy(a[:, 0], a[:, i+1], color=COLORS[i], marker=('o', 's', '^')[i],
                    fillstyle='none' if i == 1 else 'full', label='CH%d | NA 0.4 | %+d deg'%(i+1, (45, -45, 0)[i]))
    ax.set(xlabel='Particle diameter (nm)', ylabel='Collected DSC (um^2)', title='Sphere on wafer | 266 nm | 75 deg | p polarization')
    ax.grid(alpha=.18); ax.legend()
    save(fig, 'case02_size_sweep')
    configs = [(.2, .2, .2), (.4, .4, .4), (.5, .3, 0.)]
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.2), constrained_layout=True)
    aperture_rows = []
    for ax, nas in zip(axes, configs):
        o = ObjectiveCollectionParameters(channel_numerical_apertures=nas)
        view = objective_far_field_view(o, reference)
        art = disk(ax, reference['theta'][:, 0], reference['global_phi'][0], reference['intensity'], 'NA '+ '/'.join('%g'%n for n in nas))
        for i, boundary in enumerate(view['cone_boundaries']):
            if nas[i] > 0:
                b = np.asarray(boundary)
                ax.plot(b[0], b[1], color=COLORS[i], lw=1.7)
        ax.annotate('Source -KY', xy=(0, -.97), xytext=(-.9, -.73), fontsize=7,
                    arrowprops={'arrowstyle': '->', 'color': '#333333'})
        energy, geometry = view['energy_collection_fraction'], view['geometric_collection_fraction']
        ax.text(.5, -.30, 'Energy union: %.2f%%\nGeometry union: %.2f%%'%(100*energy, 100*geometry),
                transform=ax.transAxes, ha='center', fontsize=10)
        aperture_rows.append(dict(na=nas, energy_union_fraction=energy, geometric_union_fraction=geometry,
                                  energy_channel_fractions=view['energy_channel_fractions']))
    fig.colorbar(art, ax=list(axes), label='DSC (um^2/sr)', shrink=.65)
    save(fig, 'case03_apertures')
    return dict(parameters=asdict(p), objective=asdict(objective), rows=rows, aperture_cases=aperture_rows,
                scope='New SCATMECH Bobbert-Vlieger calculations. 4-degree source grid; interpolated collection quadrature. No gap parameter in this SCATMECH setup. Illustrative n/k held fixed, not calibrated material dispersion.')


def haze_case():
    base = ROOT/'build/channel_symmetry_20261016'
    r = read(base/'nist_haze_report.json')
    with np.load(base/'nist_haze_maps.npz') as fields:
        fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.0), constrained_layout=True)
        peak = max(fields['nist_plane_w_sr'].max(), fields['finite_w_sr'].max())*1e9
        for ax, key, title in zip(axes[:2], ('nist_plane_w_sr', 'finite_w_sr'), ('NIST plane-wave PSD', 'Finite coherent PSD reference')):
            art = disk(ax, fields['theta'], fields['phi'], fields[key]*1e9, title, peak)
        fig.colorbar(art, ax=list(axes[:2]), label='nW/sr', shrink=.7)
        ch = np.arange(3); width = .34
        axes[2].bar(ch-width/2, [c['nist_plane_psd_power_w']*1e9 for c in r['channels']], width, color=COLORS[0], label='NIST plane wave')
        axes[2].bar(ch+width/2, [c['finite_coherent_psd_power_w']*1e9 for c in r['channels']], width, color=COLORS[1], label='Finite covariance')
        axes[2].set(xticks=ch, xticklabels=['CH1', 'CH2', 'CH3'], ylabel='Collected power (nW)', title='Same absolute power scale')
        axes[2].legend(fontsize=8)
        save(fig, 'case04_haze_reference')
    # Replot original and low-RMS hybrid images on the SAME absolute scale.
    original = ROOT/'build/field_coverage_20261016/result'
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0), constrained_layout=True)
    peak = 0
    datasets = []
    for path in (original, LOW):
        with np.load(path/'image_fields.npz') as f:
            datasets.append({k: f[k] for k in ('ch1_image_axis_um', 'ch2_image_axis_um', 'ch1_image_with_w_um2', 'ch2_image_with_w_um2')})
            peak = max(peak, *(float(f['ch%d_image_with_w_um2'%i].max()) for i in (1, 2)))
    for row, (data, rms) in enumerate(zip(datasets, (2., .2))):
        for col in range(2):
            i = col+1; axis = data['ch%d_image_axis_um'%i]
            art = axes[row, col].pcolormesh(axis, axis, data['ch%d_image_with_w_um2'%i]*1e12,
                cmap='magma', shading='auto', vmin=0, vmax=peak*1e12)
            axes[row, col].set(title='RMS %g nm | CH%d'%(rms, i), xlabel='Image x (um)', ylabel='Image y (um)',
                               xlim=(-35,35), ylim=(-35,35), aspect='equal')
    fig.colorbar(art, ax=axes.ravel().tolist(), label='pW/um^2', shrink=.7)
    save(fig, 'case04_haze_images')
    return r


def camera_case():
    from detector_model import DetectorParameters
    from defect_detection import calculate_detection
    from gray_acquisition import acquire_gray
    source = read(LOW/'imaging_report.json')
    p = DetectorParameters(**read(SMALL/'imaging_report.json')['settings']['detector'])
    rows = []; gray = []; archives = {}
    with np.load(LOW/'image_fields.npz') as saved:
        for exposure in (.001, .005, .01):
            parameters = replace(p, exposure_ms=exposure)
            channels = []
            for i, ch in enumerate(source['channels']):
                prefix = 'ch%d_'%(i+1)
                def optical(name):
                    return dict(image_axis_um=saved[prefix+'image_axis_um'], irradiance_w_um2=saved[prefix+'image_'+name+'_w_um2'], power_w=ch[name]['image_power_w'])
                d = calculate_detection(optical('with'), optical('without'), 266, parameters, np.random.default_rng(i), i)
                d.update(acquire_gray(d, parameters, 1000, 20261017, i))
                channels.append(dict(channel=i+1, gray_snr=d['gray_snr'], gray_valid=d['gray_valid'],
                    clipped_values=d['gray_full_well_values']+d['gray_adc_rail_values']))
                if i == 0:
                    gray.append(d['gray_first_sample_dn'])
                prefix_out = 'us%g_ch%d_'%(exposure*1000, i+1)
                for name in ('gray_first_sample_dn', 'gray_first_reference_dn', 'gray_roi_series_dn'):
                    archives[prefix_out+name] = d[name]
            rows.append(dict(exposure_us=exposure*1000, channels=channels))
    np.savez_compressed(OUT/'camera_gray_examples.npz', **archives)
    with (OUT/'camera_snr.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['exposure_us', 'channel', 'temporal_roi_gray_snr', 'gray_valid', 'clipped_values'])
        for row in rows:
            for ch in row['channels']:
                writer.writerow([row['exposure_us'], ch['channel'], ch['gray_snr'], ch['gray_valid'], ch['clipped_values']])
    fig, axes = plt.subplots(1, 4, figsize=(11.4, 3.7), constrained_layout=True)
    ceiling = max(float(a.max()) for a in gray)
    for ax, a, row in zip(axes[:3], gray, rows):
        art = ax.imshow(a, origin='lower', cmap='gray', vmin=60, vmax=ceiling, extent=(-50,50,-50,50))
        metric = row['channels'][0]
        ax.set(title='%g us | CH1\nSNR %s'%(row['exposure_us'], '%.2f'%metric['gray_snr'] if metric['gray_valid'] else 'INVALID'), xlabel='Sensor x (um)', ylabel='Sensor y (um)')
    fig.colorbar(art, ax=list(axes[:3]), label='Raw DN (shared window)', shrink=.65)
    for i in range(3):
        valid = [r for r in rows if r['channels'][i]['gray_valid']]
        axes[3].plot([r['exposure_us'] for r in valid], [r['channels'][i]['gray_snr'] for r in valid],
                      marker='o', color=COLORS[i], label='CH%d'%(i+1))
    axes[3].set(xlabel='Exposure (us)', ylabel='Temporal ROI gray SNR', title='1000 frame pairs\nTemporal ROI SNR')
    axes[3].legend(); axes[3].grid(alpha=.2)
    save(fig, 'case05_camera')
    return dict(source=str(LOW.relative_to(ROOT)), detector=asdict(p), pairs=1000, seed=20261017, results=rows,
        scope='New camera acquisition simulation on saved optical irradiance; no TDI, no new FEM solves. Same static haze in sample/reference, no spatial-background mismatch. Gain/offset assumptions are not camera calibration.')


def pattern_case():
    from fullwave_prototype.defects import footprint_nm
    from matplotlib.patches import Polygon
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.6), constrained_layout=True)
    for row, path in enumerate((SMALL, LARGE)):
        report = read(path/'imaging_report.json'); p = report['scattering_parameters']
        outline = footprint_nm(p)/1000
        limit = p['defect_length_nm']/1000*.65
        axes[row, 0].add_patch(Polygon(outline, facecolor='#cc7197', edgecolor='#873056'))
        axes[row, 0].set(title='Cross | L %g um, W %g um'%(p['defect_length_nm']/1000, p['defect_width_nm']/1000),
            xlim=(-limit,limit), ylim=(-limit,limit), aspect='equal', xlabel='Object x (um)', ylabel='Object y (um)')
        with np.load(path/'image_fields.npz') as data:
            axis = data['ch3_image_axis_um']; val = data['ch3_image_with_w_um2']
            lim = 12 if row == 0 else 32
            art = axes[row, 1].pcolormesh(axis, axis, val, shading='auto', cmap='magma', norm=PowerNorm(.5, vmin=0, vmax=float(val.max())))
            axes[row, 1].set(xlim=(-lim,lim), ylim=(-lim,lim), aspect='equal', title='CH3 optical | 20x, NA 0.5', xlabel='Image x (um)', ylabel='Image y (um)')
            fig.colorbar(art, ax=axes[row, 1], label='W/um^2 (sqrt colors)', shrink=.75)
            raw = data['ch3_detector_gray_first_sample_dn']
            art = axes[row, 2].imshow(raw, origin='lower', cmap='gray', vmin=60, vmax=max(65, float(raw.max())), extent=(-50,50,-50,50))
            axes[row, 2].set(title='One simulated raw frame | 10 us', xlabel='Sensor x (um)', ylabel='Sensor y (um)')
            fig.colorbar(art, ax=axes[row, 2], label='DN (separate display windows)', shrink=.75)
    save(fig, 'case06_patterns')
    return read(ROOT/'build/cross_2um_20261017/mesh_sensitivity.json')


def main():
    ASSETS.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titleweight': 'semibold', 'axes.spines.top': False, 'axes.spines.right': False})
    results = dict(pipeline=workflow_case(), sweep=sweep_and_apertures(), haze=haze_case(), camera=camera_case(), patterns=pattern_case())
    sources = [LOW/'imaging_report.json', LOW/'local_scattering/report.json',
        ROOT/'build/field_coverage_20261016/result/imaging_report.json',
        ROOT/'build/channel_symmetry_20261016/nist_haze_report.json',
        ROOT/'build/channel_symmetry_20261016/ensemble_report.json',
        SMALL/'imaging_report.json', LARGE/'imaging_report.json',
        ROOT/'build/cross_2um_20261017/mesh_sensitivity.json']
    results['provenance'] = dict(showcase_date='2026-10-17',
        inputs_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        generated_now=['18 SCATMECH particle diameters', '3 aperture configurations', '3 camera exposures, 1000 pairs each'],
        archived_solves='FEM, hybrid and statistical PSD examples use preserved actual solver results, not new FEM solves.')
    (OUT/'case_manifest.json').write_text(json.dumps(results, ensure_ascii=False, indent=2, default=lambda x: np.asarray(x).tolist()), encoding='utf-8')
    print('Showcase data and figures complete: '+str(OUT), flush=True)


if __name__ == '__main__':
    main()
