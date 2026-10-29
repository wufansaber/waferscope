"""S4 finite trench motif + sphere: continuous volume radiation and coherent images.

The motif repeats with the S4 supercell; the extracted one-cell current radiates
in a flat-wafer background. This is an isolation approximation, not AFMM/PML.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import numpy as np

from rcwa_prototype.geometry import Scene
from s4_prototype.finite_patch import FinitePatch, volume_quadrature
from s4_prototype.solver import solve
from grating_scene.pipeline import write_json, digest
from grating_scene.runtime import default_s4
from fullwave_prototype.parameters import Parameters
from fullwave_prototype.radiation import radiation_amplitude
from fullwave_prototype.imaging import ImagingSettings, calculate_images
from rcwa_prototype.angular_distribution import hemisphere_grid


def worker(args):
    scene = Scene(lx_um=args.cell_x, periods_y=args.cell_periods_y,
        harmonics_x=args.mx, harmonics_y=args.my, slices=args.slices, gap_um=args.gap_nm/1000,
        grating_method='fft')
    scene.validate()
    patch = FinitePatch(length_um=args.length, trench_count=args.trenches)
    patch.validate(scene)
    settings = ImagingSettings(numerical_apertures=(args.na,)*3, tilt_degrees=(45.,-45.,0.),
        magnification=args.magnification, pupil_samples=args.pupil_samples,
        padding=4, irradiance_w_um2=args.irradiance)
    settings.validate()
    p = Parameters(wavelength_nm=1000*scene.wavelength_um, diameter_nm=1000*scene.diameter_um,
        gap_nm=1000*scene.gap_um, substrate_n=scene.substrate_n, substrate_k=scene.substrate_k,
        particle_n=scene.particle_n, particle_k=scene.particle_k, incidence_deg=scene.incidence_deg)
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use an empty output folder; previous results cannot be overwritten')
    output.mkdir(parents=True, exist_ok=True)
    points, weights, is_particle = volume_quadrature(scene, patch, lateral=args.source_xy,
        vertical=args.source_z)
    write_json(output/'config.json', dict(scene=asdict(scene), finite_patch=asdict(patch),
        imaging=asdict(settings), source_xy=args.source_xy, source_z=args.source_z,
        polar_samples=args.polar_samples, azimuthal_samples=args.azimuthal_samples))
    np.savez_compressed(output/'source_quadrature.npz', points_um=points,
        weights_um3=weights, is_particle=is_particle)
    arrays = hemisphere_grid(args.polar_samples, args.azimuthal_samples)
    rays = arrays['directions'].reshape(-1,3)
    data, native, sources = {}, {}, {}
    for label, particle in (('reference',False),('sample',True)):
        print('FINITE_S4 '+label+' native solve and source sampling', flush=True)
        data[label], native[label] = solve(scene,output/label,particle=particle,
            near_points=points.tolist(),executable=args.executable,timeout=args.timeout,
            stripe_rule=False,finite_patch=patch)
        contrast = np.where(is_particle, scene.eps_particle-1 if particle else 0,
                            1-scene.eps_substrate)
        sources[label] = data[label]['near_electric']*(weights*contrast)[:,None]
        np.savez_compressed(output/(label+'_sources.npz'),points_um=points,
            polarization_weights=sources[label],is_particle=is_particle)
        print('FINITE_S4 '+label+' continuous upper hemisphere', flush=True)
        field = radiation_amplitude(points,sources[label],rays,p)
        arrays[label+'_far_sp_um'] = field.reshape(arrays['directions'].shape[:-1]+(2,))
        arrays[label+'_dsc_um2_sr'] = np.sum(abs(field)**2,axis=-1).reshape(arrays['directions'].shape[:-1])
    delta = arrays['sample_far_sp_um']-arrays['reference_far_sp_um']
    arrays['delta_far_sp_um'] = delta
    arrays['delta_dsc_um2_sr'] = np.sum(abs(delta)**2,axis=-1)
    arrays['difference_dsc_um2_sr'] = arrays['sample_dsc_um2_sr']-arrays['reference_dsc_um2_sr']
    arrays['interference_dsc_um2_sr'] = 2*np.real(np.sum(delta*arrays['reference_far_sp_um'].conj(),axis=-1))
    # Native discrete modes are retained and used only for an extraction diagnostic.
    valid = data['sample']['propagating']
    nr = data['sample']['directions'][valid]
    nf = radiation_amplitude(points,sources['sample']-sources['reference'],nr,p)
    rho = np.hypot(nr[:,0],nr[:,1])
    s = np.column_stack((-nr[:,1],nr[:,0],np.zeros(len(nr))))/np.maximum(rho[:,None],1e-15)
    s[rho<1e-12] = (-1,0,0)
    electric = s*nf[:,:1]+np.cross(nr,s)*nf[:,1:]
    rebuilt = 2j*np.pi/(scene.lx_um*scene.ly_um*p.k0*nr[:,2,None])*electric
    direct = (data['sample']['electric_reflected_z0']-data['reference']['electric_reflected_z0'])[valid]
    norm = np.sum(nr[:,2,None]*abs(direct)**2)
    closure = float(np.sqrt(np.sum(nr[:,2,None]*abs(rebuilt-direct)**2)/norm)) if norm>0 else None
    np.savez_compressed(output/'modal_extraction_diagnostic.npz',directions=nr,
        rebuilt_electric=rebuilt,direct_electric=direct)
    print('FINITE_S4 native difference-mode / volume extraction relative L2: '+str(closure),flush=True)
    extraction_passed = closure is not None and np.isfinite(closure) and closure <= args.extraction_tolerance
    write_json(output/'extraction_diagnostic.json',dict(
        native_delta_vs_volume_relative_complex_l2=closure,
        tolerance=args.extraction_tolerance,passed=bool(extraction_passed),
        scope='Same-cell native difference modes versus volume reconstruction; not isolation certification.'))
    if not extraction_passed:
        raise RuntimeError('Native S4 / volume extraction consistency failed; refine source quadrature '
                           'and Fourier representation before publishing image results')
    np.savez_compressed(output/'continuous_far_field.npz',**arrays)
    def provider(directions,state):
        return radiation_amplitude(points,sources['sample' if state=='with' else 'reference'],directions,p)
    print('FINITE_S4 coherent objective pupils and image planes',flush=True)
    images, imaging = calculate_images({},p,settings,retain_complex_fields=True,
                                        amplitude_provider=provider)
    np.savez_compressed(output/'image_fields.npz',**images)
    write_json(output/'imaging_report.json',imaging)
    polar, azimuth = np.meshgrid(arrays['theta_deg'],arrays['phi_deg'],indexing='ij')
    table = np.column_stack((polar.ravel(),azimuth.ravel(),rays,
        arrays['solid_angle_weights'].ravel(),
        arrays['reference_far_sp_um'].real.reshape(-1,2), arrays['reference_far_sp_um'].imag.reshape(-1,2),
        arrays['sample_far_sp_um'].real.reshape(-1,2), arrays['sample_far_sp_um'].imag.reshape(-1,2),
        *[arrays[key+'_dsc_um2_sr'].ravel() for key in ('reference','sample','delta','difference','interference')]))
    np.savetxt(output/'continuous_far_field.csv',table,delimiter=',',comments='',
        header='theta_deg,phi_deg,nx,ny,nz,solid_angle_sr,ref_Fs_real_um,ref_Fp_real_um,ref_Fs_imag_um,ref_Fp_imag_um,sample_Fs_real_um,sample_Fp_real_um,sample_Fs_imag_um,sample_Fp_imag_um,ref_DSC_um2_sr,sample_DSC_um2_sr,delta_DSC_um2_sr,intensity_difference_um2_sr,interference_um2_sr')
    plot_results(output,scene,patch,arrays,images,imaging)
    report = dict(format='waferscope-s4-finite-mixed',version=1,
        status='EXPERIMENTAL_FINITE_MOTIF_PERIODIC_CURRENT_ISOLATION_UNCHECKED',
        backend='S4 native',scene=asdict(scene),finite_patch=asdict(patch),native=native,
        imaging=imaging,source_quadrature_points=len(points),
        source_pit_volume_um3=float(weights[~is_particle].sum()),
        source_sliced_sphere_volume_um3=float(weights[is_particle].sum()),
        native_delta_vs_volume_relative_complex_l2=closure,
        near_far_extraction_passed=bool(extraction_passed),
        near_far_extraction_tolerance=args.extraction_tolerance,
        angular_integrals_um2={key:float(np.sum(arrays['solid_angle_weights']*arrays[key+'_dsc_um2_sr']))
                               for key in ('reference','sample','delta','difference','interference')},
        decomposition_max_absolute_error=float(np.max(abs(arrays['difference_dsc_um2_sr']
            -arrays['delta_dsc_um2_sr']-arrays['interference_dsc_um2_sr']))),
        limitations=['Finite trench motif and sphere repeat with the S4 supercell; isolation is not certified.',
            'Native S4 orders are discrete; continuous radiation is reconstructed from one-cell pit AND particle currents.',
            'Finite trench end scattering is included in the motif; no analytic invariant-stripe shortcut is used.',
            'Reference is the same finite trench motif without the particle, solved with a common phase origin.',
            'Plane-wave illumination, scattered field relative to a flat wafer; the flat-wafer specular delta beam is excluded.',
            'Classical legacy material constants, sliced sphere and finite Fourier basis require convergence checks.',
            'Images are ideal aplanatic dark-field objective images; no camera noise or real lens prescription.'])
    report['files_sha256']={str(path.relative_to(output)).replace('\\','/'):digest(path)
        for path in output.rglob('*') if path.is_file()}
    write_json(output/'finite_mixed_report.json',report)
    print(json.dumps(dict(status=report['status'],native_extraction_relative_l2=closure,
        image_power_w=[row['with']['image_power_w'] for row in imaging['channels']]),indent=2),flush=True)


def plot_results(output,scene,patch,far,images,imaging):
    from matplotlib.figure import Figure
    from matplotlib.colors import Normalize,TwoSlopeNorm
    from grating_scene.plotting import positive_norm
    fig=Figure(figsize=(12,8))
    values=[far[key+'_dsc_um2_sr'] for key in ('reference','sample','delta','difference')]
    shared=positive_norm(values[:2],True)
    for i,(title,value) in enumerate(zip(('Finite grating reference','Finite grating + sphere',
                                         'Changed-field DSC','Sample - reference intensity'),values),1):
        ax=fig.add_subplot(2,2,i)
        if i==4:
            limit=max(float(abs(value).max()),1e-30)
            norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
        else:
            norm=shared if i<3 else positive_norm([value],True)
        im=ax.pcolormesh(far['phi_deg'],far['theta_deg'],value,shading='nearest',
                         norm=norm,cmap='RdBu_r' if i==4 else 'viridis')
        ax.set(title=title,xlabel='Azimuth phi (deg)',ylabel='Polar theta (deg)',xlim=(0,360),ylim=(0,90))
        fig.colorbar(im,ax=ax,label='um^2 / sr')
    fig.suptitle('S4 finite-motif current radiation | cell {:.3g} x {:.3g} um | isolation unchecked'.format(scene.lx_um,scene.ly_um))
    fig.tight_layout(); fig.savefig(output/'continuous_far_field.png',dpi=160)
    fig=Figure(figsize=(15,10))
    for i in range(1,4):
        axis=images['ch{}_image_axis_um'.format(i)]
        chosen=abs(axis)<=18*imaging['settings']['magnification']/20
        xy=axis[chosen]
        keys=('without','with','difference','change')
        vs=[images['ch{}_image_{}_w_um2'.format(i,key)][np.ix_(chosen,chosen)] for key in keys]
        shared=Normalize(0,max(float(vs[0].max()),float(vs[1].max()),1e-30))
        for j,(key,value) in enumerate(zip(keys,vs),1):
            ax=fig.add_subplot(3,4,(i-1)*4+j)
            if key=='difference':
                limit=max(float(abs(value).max()),1e-30)
                norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
            else:
                norm=shared if j<3 else Normalize(0,max(float(value.max()),1e-30))
            im=ax.pcolormesh(xy,xy,value,norm=norm,cmap='RdBu_r' if j==3 else 'magma',shading='auto')
            labels={'without':'Grating reference','with':'Grating + sphere','difference':'Sample - reference','change':'Changed-field intensity'}
            ax.set(title='CH{} {}'.format(i,labels[key]),xlabel="Image x' (um)",ylabel="Image y' (um)",aspect='equal')
            fig.colorbar(im,ax=ax,label='W / um^2',shrink=.75)
    fig.suptitle('Coherent dark-field images | NA {}, tilt +45/-45/0 deg, M={} | full complex field interference retained'.format(
        imaging['settings']['numerical_apertures'][0],imaging['settings']['magnification']),fontsize=11)
    fig.tight_layout(); fig.savefig(output/'mixed_images.png',dpi=150)
    fig=Figure(figsize=(8,6)); ax=fig.add_subplot(111,projection='3d')
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    length=patch.length_um/2; w=(scene.pitch_um-scene.ridge_width_um)/2
    faces=[]
    for center in patch.centers(scene):
        y0,y1=center-w,center+w; z=-scene.depth_um
        faces.extend([[(a,y0,0),(b,y0,0),(b,y0,z),(a,y0,z)] for a,b in [(-length,length)]])
        faces.append([(-length,y1,0),(length,y1,0),(length,y1,z),(-length,y1,z)])
        faces.append([(-length,y0,z),(length,y0,z),(length,y1,z),(-length,y1,z)])
        for x in (-length,length): faces.append([(x,y0,0),(x,y1,0),(x,y1,z),(x,y0,z)])
    ax.add_collection3d(Poly3DCollection(faces,facecolors='#547e9a',edgecolors='#25394a',alpha=.7))
    phi,theta=np.meshgrid(np.linspace(0,2*np.pi,33),np.linspace(0,np.pi,25))
    r=scene.diameter_um/2
    ax.plot_surface(r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),
                    scene.gap_um+r+r*np.cos(theta),color='#c96845',alpha=.9)
    ax.set(xlim=(-scene.lx_um/2,scene.lx_um/2),ylim=(-scene.ly_um/2,scene.ly_um/2),
        zlim=(-scene.depth_um,scene.cap_top_um*1.5),xlabel='x (um)',ylabel='y (um)',zlabel='z (um)',
        title='Finite trench motif + sphere in S4 cell; outside is flat wafer')
    ax.set_box_aspect((scene.lx_um,scene.ly_um,.5)); fig.tight_layout(); fig.savefig(output/'geometry.png',dpi=160)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--executable',type=Path,default=default_s4())
    parser.add_argument('--cell-x',type=float,default=1.5)
    parser.add_argument('--cell-periods-y',type=int,default=4)
    parser.add_argument('--mx',type=int,default=6)
    parser.add_argument('--my',type=int,default=8)
    parser.add_argument('--slices',type=int,default=8)
    parser.add_argument('--length',type=float,default=.8)
    parser.add_argument('--trenches',type=int,default=2)
    parser.add_argument('--gap-nm',type=float,default=0.)
    parser.add_argument('--source-xy',type=int,default=10)
    parser.add_argument('--source-z',type=int,default=8)
    parser.add_argument('--polar-samples',type=int,default=24)
    parser.add_argument('--azimuthal-samples',type=int,default=72)
    parser.add_argument('--na',type=float,default=.6)
    parser.add_argument('--magnification',type=float,default=20.)
    parser.add_argument('--irradiance',type=float,default=1e-6)
    parser.add_argument('--pupil-samples',type=int,default=64,choices=(32,64,96,128))
    parser.add_argument('--timeout',type=float,default=1800.)
    parser.add_argument('--max-memory-gib',type=float,default=8.)
    parser.add_argument('--extraction-tolerance',type=float,default=.02)
    parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if not np.isfinite(args.extraction_tolerance) or not 0<args.extraction_tolerance<1:
        parser.error('Extraction tolerance must be a finite relative number in (0,1)')
    if args.worker:
        worker(args)
    else:
        from fullwave_prototype.guarded_run import run_guarded
        from rcwa_prototype.geometry import matrix_memory_estimate
        import psutil
        scene=Scene(lx_um=args.cell_x,periods_y=args.cell_periods_y,
            harmonics_x=args.mx,harmonics_y=args.my,slices=args.slices,
            gap_um=args.gap_nm/1000,grating_method='fft')
        scene.validate()
        FinitePatch(length_um=args.length,trench_count=args.trenches).validate(scene)
        estimate=matrix_memory_estimate(scene.requested_modes,scene.slices+1)
        available=psutil.virtual_memory().available/2**30
        blocked=estimate['planning_gib']>min(args.max_memory_gib,available-3.) or available<3.25
        args.output.parent.mkdir(parents=True,exist_ok=True)
        write_json(args.output.with_suffix('.admission.json'),dict(planning=estimate,
            available_gib=available,system_reserve_gib=3.,maximum_gib=args.max_memory_gib,
            status='blocked_memory_preflight' if blocked else 'admitted'))
        if blocked:
            parser.error('S4 planning or system reserve prevents startup; see admission JSON')
        result=run_guarded([sys.executable,'-m','tools.s4_finite_mixed','--worker']+sys.argv[1:],
            args.output.with_suffix('.resources.json'),maximum=args.max_memory_gib,
            reserve=3.,timeout=args.timeout)
        if result['status']!='completed': raise SystemExit(1)


if __name__=='__main__': main()
