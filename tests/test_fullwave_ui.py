"""Qt bridge checks without importing NGSolve or starting a heavy solve."""
import importlib.util
from dataclasses import asdict
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
HAS_QT = importlib.util.find_spec('PyQt5') is not None
if HAS_QT:
    from PyQt5.QtCore import QProcess
    from PyQt5.QtGui import QPixmap
    from PyQt5.QtWidgets import QApplication
    from fullwave_ui.controls import SimulationControls
    from fullwave_ui.results import ResultsTabs
    from fullwave_ui.widget import FullwaveWidget
    from fullwave_prototype.imaging_run import build_parser, scattering_settings
    from fullwave_prototype.imaging import ImagingSettings
    from fullwave_prototype.imaging_options import imaging_settings
    from scatmech_ui.app import ScatmechStandaloneWindow
    from scatmech_ui.native_far_field import calculate_far_field
    from scatmech_ui.parameters import SimulationParameters
    from scatmech_ui.plotting import PlotCanvas, PARTICLE_INTENSITY_CMAP
    from scatmech_ui.far_field_plotting import plot_2d_far_field


@unittest.skipUnless(HAS_QT, 'Qt environment required')
class FullwaveUiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_measured_psd_controls_use_file_rms_and_disable_analytic_inputs(self):
        import numpy as np
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'psd.npz'
            frequency = np.arange(-4, 5)*2.
            x, y = np.meshgrid(frequency, frequency)
            density = 2e-8*np.exp(-(.2*x)**2-(.4*y)**2)
            np.savez(source, fx_per_um=frequency, fy_per_um=frequency, psd_um4=density)
            controls = SimulationControls()
            controls.fields['roughness-model'].setCurrentText('measured_psd')
            controls.measurement_path.setText(str(source))
            arguments = controls.arguments(particle_enabled=False)
            _, roughness = scattering_settings(build_parser().parse_args(arguments), particle_enabled=False)
            self.assertAlmostEqual(roughness.rms_nm, 1000*np.sqrt(density.sum()*4))
            self.assertFalse(controls.fields['rms-nm'].isEnabled())
            self.assertFalse(controls.fields['correlation-length-um'].isEnabled())
            self.assertIn('--measurement-path', arguments)
            controls.close()

    def test_convergence_dialog_preserves_enabled_channels_and_scan_levels(self):
        from fullwave_ui.convergence_dialog import ConvergenceDialog
        from fullwave_prototype.haze_convergence import default_plan
        controls = SimulationControls()
        args = build_parser().parse_args(controls.arguments(particle_enabled=False))
        parameters, roughness = scattering_settings(args, particle_enabled=False)
        dialog = ConvergenceDialog(default_plan(parameters, roughness, 4, (.4, 0, .5), (45, -45, 0)))
        dialog.levels['patch_width_um'].clear()
        dialog.accept()
        self.assertEqual(dialog.result(), dialog.Accepted)
        self.assertNotIn('patch_width_um', dialog.plan['scans'])
        self.assertEqual(dialog.plan['numerical_apertures'], [.4, 0, .5])
        dialog.close()
        controls.close()

    def test_convergence_results_disable_missing_paired_figure(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'haze_convergence_report.json').write_text(
                json.dumps(dict(status='complete')), encoding='utf-8')
            pixmap = QPixmap(64, 64)
            pixmap.fill()
            self.assertTrue(pixmap.save(str(folder/'haze_convergence.png')))
            error_path = folder/'haze_paired_errors.png'
            self.assertTrue(pixmap.save(str(error_path)))
            results = ResultsTabs()
            results.load_convergence(folder)
            index = results.convergence_tabs.indexOf(results.paired_error_view)
            self.assertTrue(results.convergence_tabs.isTabEnabled(index))
            self.assertIsNotNone(results.paired_error_view.item)
            error_path.unlink()
            results.load_convergence(folder)
            self.assertFalse(results.convergence_tabs.isTabEnabled(index))
            results.close()

    def test_finite_illumination_page_is_loaded_and_old_results_clear_it(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'illumination_report.json').write_text(json.dumps(
                dict(model='gaussian', incident_power_fraction_inside_rough_patch=.02)), encoding='utf-8')
            pixmap = QPixmap(64, 64)
            pixmap.fill()
            pixmap.save(str(folder/'illumination.png'))
            results = ResultsTabs()
            results.load_illumination(folder)
            self.assertIsNotNone(results.illumination_view.item)
            self.assertIn('gaussian', results.illumination_details.toPlainText())
            (folder/'illumination.png').unlink()
            results.load_illumination(folder)
            self.assertIsNone(results.illumination_view.item)
            self.assertFalse(results.isTabEnabled(results.indexOf(results.illumination_tabs)))
            results.close()

    def test_finite_haze_launch_does_not_pass_a_plane_wave_nist_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            engine = root/'fem.exe'
            engine.touch()
            widget = FullwaveWidget()
            widget.controls.fields['illumination-fit'].setCurrentText('manual')
            widget.interpreter.setText(str(engine))
            widget.ensemble_count.setValue(2)
            with patch('fullwave_ui.widget.ROOT', root), \
                 patch('fullwave_ui.widget.results_root', return_value=root), \
                 patch('fullwave_prototype.haze_crosscheck.statistical_reference') as reference, \
                 patch.object(widget.process, 'start') as start:
                widget.start_ensemble()
            start.assert_called_once()
            reference.assert_not_called()
            arguments = start.call_args.args[1]
            self.assertIn('gaussian', arguments)
            self.assertNotIn('--reference', arguments)
            widget.close()

    def test_hybrid_launch_keeps_physical_beam_and_uses_guarded_separate_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            engine = root/'fem.exe'
            engine.touch()
            widget = FullwaveWidget()
            widget.interpreter.setText(str(engine))
            widget.controls.fields['na2'].setValue(0)
            with patch('fullwave_ui.hybrid_job.application_root', return_value=root), \
                 patch('fullwave_ui.hybrid_job.results_root', return_value=root), \
                 patch.object(widget.process, 'start') as start:
                widget.start_hybrid_job()
            start.assert_called_once()
            arguments = start.call_args.args[1]
            self.assertIn('fullwave_prototype.guarded_run', arguments)
            self.assertIn('hybrid_run', arguments)
            self.assertNotIn('--source', arguments)
            self.assertEqual(widget.job_kind, 'hybrid')
            inputs = json.loads((widget.output/'hybrid_request_input.json').read_text(encoding='utf-8'))
            self.assertEqual(inputs['parameters']['beam_profile'], 'gaussian')
            self.assertGreater(inputs['parameters']['beam_length_um'], 12)
            self.assertEqual(inputs['imaging_settings']['illumination_fit'], 'camera')
            self.assertEqual(inputs['roughness']['width_um'], .2)
            self.assertAlmostEqual(inputs['roughness']['taper_fraction'], .125)
            self.assertEqual(inputs['imaging_settings']['numerical_apertures'][1], 0)
            self.assertFalse(widget.hybrid_reimage.isEnabled())
            widget.close()

    def test_hybrid_result_is_not_routed_to_local_fem_imaging(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'hybrid_report.json').write_text('{}', encoding='utf-8')
            widget = FullwaveWidget()
            with patch.object(widget, 'load_hybrid_result', return_value=True) as load:
                self.assertTrue(widget.load_result(folder))
            load.assert_called_once_with(folder)
            widget.close()

    def test_hybrid_loader_rejects_other_models_before_loading_figures(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            (folder/'hybrid_report.json').write_text(json.dumps(dict(model='local_fem')), encoding='utf-8')
            widget = FullwaveWidget()
            with patch('fullwave_ui.widget.QMessageBox.warning') as warning, \
                 patch.object(widget.results, 'load') as load:
                self.assertFalse(widget.load_hybrid_result(folder))
            load.assert_not_called()
            warning.assert_called_once()
            widget.close()

    def test_copy_fem_haze_inputs_to_scatmech(self):
        window = ScatmechStandaloneWindow()
        self.assertEqual(window.simulation_widget.roughness_group.psd.currentData(), 'gaussian')
        with tempfile.TemporaryDirectory() as directory:
            result = Path(directory)
            source = result/'scattering'
            source.mkdir()
            physical = dict(wavelength_nm=266, substrate_n=5.4, substrate_k=.3,
                            incidence_deg=75, polarization='p', beam_waist_um=0)
            roughness = dict(model='gaussian', rms_nm=1.5,
                             correlation_length_um=.06, width_um=.24)
            report = dict(parameters=physical, roughness=roughness, angular_step_deg=5)
            (source/'report.json').write_text(json.dumps(report), encoding='utf-8')
            settings = ImagingSettings(numerical_apertures=(.4, .5, .5),
                tilt_degrees=(30., -45., 0.),
                irradiance_w_um2=1/(4*math.cos(math.radians(75))))
            (result/'imaging_report.json').write_text(json.dumps(dict(settings=asdict(settings))),
                                                      encoding='utf-8')
            window.fullwave_widget.loaded_source = source
            window.fullwave_widget.loaded_result = result
            window.fullwave_widget.controls.fields['wavelength-nm'].setValue(355)
            self.assertFalse(window.simulation_widget.roughness_group.match_fem_button.isHidden())
            window.simulation_widget.roughness_group.match_fem_button.click()
            target = window.simulation_widget
            params = target.parameters()
            self.assertEqual((params.wavelength_nm, params.wafer_n, params.wafer_k),
                             (266, 5.4, .3))
            self.assertEqual((params.roughness_psd, params.haze_only,
                              params.roughness_rms_nm, params.roughness_correlation_length_um),
                             ('gaussian', True, 1.5, .06))
            self.assertAlmostEqual(params.illuminated_area_um2, .24**2)
            self.assertEqual(target.objective_collection_parameters().channel_numerical_apertures[0], .4)
            self.assertEqual(target.objective_collection_parameters().tilt_degrees[0], 30)
            self.assertAlmostEqual(target.imaging_group.irradiance.value(),
                                   1000/(4*math.cos(math.radians(75))), places=4)
            report['roughness'].update(model='abc', psd_exponent=3.5)
            (source/'report.json').write_text(json.dumps(report), encoding='utf-8')
            window.simulation_widget.roughness_group.match_fem_button.click()
            abc_params = window.simulation_widget.parameters()
            self.assertEqual(abc_params.roughness_psd, 'abc')
            self.assertAlmostEqual(abc_params.roughness_exponent, 3.5)
            self.assertAlmostEqual(abc_params.roughness_correlation_length_um, .06)
            report['roughness']['model'] = 'legacy'
            (source/'report.json').write_text(json.dumps(report), encoding='utf-8')
            with patch('scatmech_ui.app.QMessageBox.warning') as warning:
                window.simulation_widget.roughness_group.match_fem_button.click()
            self.assertIn('legacy', warning.call_args.args[2])
        window.close()

    def test_haze_only_plot_uses_absolute_dsc_and_correct_label(self):
        result = calculate_far_field(SimulationParameters(
            wavelength_nm=266, roughness_enabled=True, haze_only=True,
            roughness_psd='gaussian', angular_step_deg=15))
        canvas = PlotCanvas()
        plot_2d_far_field(canvas, result, PARTICLE_INTENSITY_CMAP)
        axis = canvas.figure.axes[0]
        self.assertIn('Haze only', axis.get_title())
        self.assertNotIn('Particle', axis.get_title())
        self.assertAlmostEqual(axis.collections[0].get_clim()[1], result['intensity'].max())

    def test_haze_ensemble_launch_and_result_tab(self):
        import numpy as np
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            engine = root/'fem.exe'
            engine.touch()
            widget = FullwaveWidget()
            self.assertEqual(widget.ensemble_count.value(), 16)
            widget.interpreter.setText(str(engine))
            widget.ensemble_count.setValue(2)
            widget.controls.fields['beam-profile'].setCurrentText('plane')
            with patch('fullwave_ui.widget.ROOT', root), \
                 patch('fullwave_ui.widget.results_root', return_value=root), \
                 patch('fullwave_prototype.haze_crosscheck.statistical_reference',
                       return_value=np.ones((19, 72))), \
                 patch('fullwave_prototype.haze_crosscheck.finite_patch_reference',
                       return_value=np.full((19, 72), .8)), \
                 patch.object(widget.process, 'start') as start:
                widget.start_ensemble()
            start.assert_called_once()
            arguments = start.call_args.args[1]
            self.assertIn('haze_ensemble', arguments)
            self.assertIn('--reference', arguments)
            self.assertEqual(widget.job_kind, 'ensemble')
            with np.load(widget.output/'nist_psd_reference.npz') as saved:
                self.assertEqual(saved['dsc'].shape, (19, 72))
                self.assertEqual(saved['finite_patch_dsc'].shape, (19, 72))
                self.assertEqual(str(saved['psd_model']), 'gaussian')
            widget.set_running(False)
            report = dict(status='complete', seeds=[7, 8], warnings=['Preview only.'],
                          nist_reference=True, finite_patch_reference=True,
                          metrics=dict(samples=2, comparison_area_um2=.04,
                                       patch_geometric_area_um2=.04,
                                       patch_effective_area_um2=.03,
                                       fem_mean_hemisphere_dsc_um2=1e-5,
                                       fem_hemisphere_standard_error_um2=1e-6,
                                       nist_hemisphere_dsc_um2=2e-5,
                                       fem_to_nist_hemisphere_ratio=.5,
                                       angular_shape_cosine=.9,
                                       finite_patch_hemisphere_dsc_um2=1.8e-5,
                                       fem_to_finite_patch_hemisphere_ratio=.56,
                                       finite_patch_angular_shape_cosine=.95,
                                       finite_patch_to_nist_hemisphere_ratio=.9,
                                       finite_patch_to_nist_angular_shape_cosine=.99,
                                       input_stationary_rms_nm=2,
                                       realized_patch_rms_mean_nm=1.7))
            (widget.output/'haze_ensemble_report.json').write_text(
                json.dumps(report), encoding='utf-8')
            image = QPixmap(80, 60)
            image.fill()
            for name in ('haze_ensemble.png', 'haze_ensemble_fem.png',
                         'haze_ensemble_nist.png', 'haze_ensemble_finite.png',
                         'haze_ensemble_stderr.png'):
                self.assertTrue(image.save(str(widget.output/name)))
            self.assertTrue(widget.load_ensemble_result(widget.output))
            self.assertTrue(widget.results.isTabEnabled(
                widget.results.indexOf(widget.results.ensemble_tabs)))
            self.assertIn('FEM / NIST', widget.results.ensemble_metrics.toPlainText())
            self.assertIn('FEM / finite patch', widget.results.ensemble_metrics.toPlainText())
            widget.close()

    def test_abc_controls_keep_b_and_exponent_separate_from_gaussian_l(self):
        controls = SimulationControls()
        controls.fields['roughness-model'].setCurrentText('abc')
        controls.fields['correlation-length-um'].setValue(.2)
        controls.fields['psd-exponent'].setValue(4.5)
        arguments = controls.arguments(particle_enabled=False)
        params, roughness = scattering_settings(build_parser().parse_args(arguments),
                                                 particle_enabled=False)
        self.assertEqual(roughness.model, 'abc')
        self.assertEqual((roughness.correlation_length_um, roughness.psd_exponent), (.2, 4.5))
        self.assertEqual(params.beam_waist_um, 0.)
        self.assertTrue(controls.fields['psd-exponent'].isEnabled())
        controls.fields['roughness-model'].setCurrentText('gaussian')
        self.assertFalse(controls.fields['psd-exponent'].isEnabled())
        controls.close()

    def test_default_illumination_and_saved_irradiance(self):
        controls = SimulationControls()
        args = build_parser().parse_args(controls.arguments())
        settings = imaging_settings(args, incidence_deg=75)
        self.assertEqual(settings.illumination_power_w, 1)
        self.assertGreater(settings.illumination_length_um, 12)
        self.assertGreater(settings.illumination_width_um, 17)
        from fullwave_prototype.finite_illumination import power_area
        p, _ = scattering_settings(args)
        self.assertEqual(settings.illumination_model, 'gaussian')
        self.assertAlmostEqual(settings.irradiance_w_um2, 1/power_area(p))
        self.assertIn('{:.6g}'.format(settings.irradiance_w_um2), controls.irradiance_label.text())
        controls.fields['illumination-power-w'].setValue(.5)
        self.assertAlmostEqual(imaging_settings(build_parser().parse_args(controls.arguments()),
                                               incidence_deg=75).irradiance_w_um2,
                               settings.irradiance_w_um2/2)
        controls.load_imaging_settings(asdict(ImagingSettings()))
        restored = imaging_settings(build_parser().parse_args(controls.arguments()),
                                    incidence_deg=75)
        self.assertAlmostEqual(restored.irradiance_w_um2, 1e-6)

    def test_controls_forward_materials_channels_and_polarization(self):
        controls = SimulationControls()
        controls.preset.setCurrentIndex(1)
        controls.fields['accuracy'].setCurrentText('balanced')
        controls.fields['order'].setCurrentText('0')
        controls.fields['direct-solver'].setCurrentText('sparsecholesky')
        controls.fields['pml-inner-fraction'].setValue(.5)
        controls.fields['particle-k'].setValue(.12)
        controls.fields['wavelength-nm'].setValue(266)
        controls.fields['polarization'].setCurrentText('s')
        controls.fields['na2'].setValue(0)
        controls.fields['tilt1'].setValue(30)
        controls.fields['threads'].setValue(1)
        args = build_parser().parse_args(controls.arguments())
        p, _ = scattering_settings(args)
        self.assertEqual(p.threads, 1)
        self.assertEqual((args.accuracy, args.pupil_samples, p.order), ('balanced', 128, 3))
        self.assertEqual((p.pml_strength, p.pml_mesh_scale), (1., 1.25))
        self.assertEqual((p.pml_order, p.pml_inner_fraction), (0, .5))
        self.assertEqual((p.particle_k, p.wavelength_nm, p.polarization), (.12, 266, 's'))
        self.assertEqual(args.na, [.5, 0, .5])
        self.assertEqual(args.tilts, [30, -45, 0])
        controls.fields['pml-order'].setCurrentText('2')
        controls.fields['pml-inner-fraction'].setValue(.25)
        controls.fields['local-mesh-scale'].setValue(.9)
        tuned, rough = scattering_settings(build_parser().parse_args(controls.arguments()))
        self.assertEqual(tuned.pml_order, 2)
        self.assertEqual(tuned.pml_inner_fraction, .25)
        self.assertAlmostEqual(tuned.particle_mesh_um, .0135)
        self.assertAlmostEqual(rough.mesh_um, .0225)
        controls.fields['accuracy'].setCurrentText('preview')
        controls.fields['mesh-scale'].setValue(1.5)
        controls.fields['local-mesh-scale'].setValue(1)
        controls.fields['particle-mesh-um'].setValue(.0015)
        controls.fields['quadrature-order'].setCurrentText('7')
        reproduction = build_parser().parse_args(controls.arguments())
        selected, surface = scattering_settings(reproduction)
        self.assertEqual((selected.order, selected.direct_solver, selected.pml_strength),
                         (2, 'sparsecholesky', 2.))
        self.assertAlmostEqual(selected.particle_mesh_um, .0015)
        self.assertAlmostEqual(selected.air_mesh_um, .12)
        self.assertAlmostEqual(surface.mesh_um, .0375)
        self.assertEqual(reproduction.quadrature_order, 7)
        self.assertNotIn('--particle-mesh-um', controls.arguments(imaging_only=True))
        self.assertNotIn('--quadrature-order', controls.arguments(imaging_only=True))
        controls.fields['gap-nm'].setValue(0)
        with self.assertRaises(ValueError):
            controls.arguments()

    def test_gaussian_roughness_controls_reach_fem(self):
        controls = SimulationControls()
        controls.fields['roughness-model'].setCurrentText('gaussian')
        self.assertEqual(controls.fields['rough-grid-points'].value(), 9)
        self.assertEqual(controls.fields['accuracy'].currentText(), 'preview')
        self.assertEqual(controls.fields['order'].currentText(), '1')
        controls.fields['correlation-length-um'].setValue(.06)
        controls.fields['rough-width-um'].setValue(.24)
        _, surface = scattering_settings(build_parser().parse_args(controls.arguments()))
        self.assertEqual((surface.model, surface.grid_points, surface.width_um,
                          surface.correlation_length_um), ('gaussian', 9, .24, .06))
        self.assertEqual(build_parser().parse_args(controls.arguments()).order, 1)
        controls.fields['rough-grid-points'].setValue(5)
        with self.assertRaisesRegex(ValueError, 'grid intervals'):
            controls.arguments()

    def test_load_is_atomic_and_keeps_result_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            report = dict(channels=[{'with': dict(image_power_w=1e-12*i)}
                                    for i in (1, 2, 3)])
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            tabs = ResultsTabs()
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            for _, name in tabs.FIGURES:
                self.assertTrue(pixmap.save(str(path/name)))
            self.assertEqual(tabs.load(path), [1, 2, 3])
            old = tabs.parameters.toPlainText()
            (path/tabs.FIGURES[2][1]).unlink()
            with self.assertRaises(ValueError):
                tabs.load(path)
            self.assertEqual(tabs.parameters.toPlainText(), old)
            self.assertIn(str(path.resolve()), old)
            report['channels'][0]['with']['image_power_w'] = float('nan')
            (path/'imaging_report.json').write_text(json.dumps(report), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'invalid image powers'):
                tabs.load(path)

    def test_10_nm_controls_reproduce_saved_preview_configuration(self):
        controls = SimulationControls()
        self.assertEqual(controls.preset.currentText(), '10 nm PSL preview')
        args = build_parser().parse_args(controls.arguments())
        p, roughness = scattering_settings(args)
        from fullwave_prototype.numerical_settings import numerical_parameters
        expected = numerical_parameters('preview', mesh_scale=1.5, diameter_nm=10., gap_nm=10.,
            particle_mesh_um=.0015, threads=2, direct_solver='umfpack', order=1,
            pml_inner_fraction=0)
        from dataclasses import replace
        expected = replace(expected, beam_profile='gaussian', beam_length_um=2., beam_width_um=2., beam_order=31)
        from fullwave_prototype.field_coverage import fit_camera_beam
        from detector_model import DetectorParameters
        expected, _ = fit_camera_beam(expected, DetectorParameters(pixels_x=16, pixels_y=16),
                                     20, (.5, .5, .5), (45, -45, 0))
        self.assertEqual(p, expected)
        self.assertEqual((roughness.model, roughness.grid_points), ('gaussian', 9))
        self.assertAlmostEqual(roughness.mesh_um, .0375)
        self.assertEqual((roughness.rms_nm, roughness.seed), (2, 7))
        self.assertEqual(args.quadrature_order, 7)
        self.assertEqual((args.pupil_samples, args.magnification, args.na, args.tilts),
                         (128, 20, [.5, .5, .5], [45, -45, 0]))
        controls.fields['diameter-nm'].setValue(25)
        controls.fields['detector'].setChecked(True)
        self.assertEqual(controls.preset.currentText(), 'Custom')
        controls.preset_button.click()
        self.assertEqual(controls.preset.currentText(), '10 nm PSL preview')
        self.assertEqual(controls.fields['diameter-nm'].value(), 10)
        self.assertFalse(controls.fields['detector'].isChecked())

    def test_runtime_preset_switches_with_parameter_set(self):
        widget = FullwaveWidget()
        self.assertEqual((widget.memory.value(), widget.timeout.value()), (10, 30))
        widget.controls.preset.setCurrentIndex(1)
        self.assertEqual((widget.controls.fields['diameter-nm'].value(), widget.timeout.value()), (100, 30))
        widget.controls.preset.setCurrentIndex(0)
        self.assertEqual((widget.controls.fields['diameter-nm'].value(), widget.timeout.value()), (10, 30))
        widget.memory.setValue(6)
        widget.timeout.setValue(20)
        widget.controls.fields['diameter-nm'].setValue(20)
        self.assertEqual((widget.memory.value(), widget.timeout.value()), (6, 20))
        widget.close()

    def test_reimage_forwards_only_imaging_controls(self):
        controls = SimulationControls()
        controls.fields['beam-waist-um'].setValue(.001)
        controls.fields['coma-nm'].setValue(25)
        controls.fields['detector'].setChecked(True)
        arguments = controls.arguments(imaging_only=True)
        self.assertNotIn('--beam-waist-um', arguments)
        self.assertNotIn('--diameter-nm', arguments)
        self.assertIn('--coma-nm', arguments)
        self.assertIn('--detector', arguments)
        self.assertIn('--illumination-power-w', arguments)
        self.assertIn('--illumination-length-um', arguments)
        self.assertIn('--illumination-width-um', arguments)

    def test_near_plane_tabs_validate_and_clear_legacy_results(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            report = dict(channels=[{'with': dict(image_power_w=1e-12)}]*3,
                          near_field_planes=[dict(plane=p) for p in ('XZ', 'XY', 'YZ')])
            report_path = path/'imaging_report.json'
            report_path.write_text(json.dumps(report), encoding='utf-8')
            tabs = ResultsTabs()
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            for name in [n for _, n in tabs.FIGURES]+['near_field_'+p+'.png' for p in ('xz', 'xy', 'yz')]:
                self.assertTrue(pixmap.save(str(path/name)))
            tabs.load(path)
            self.assertEqual(tabs.tabText(2), 'Three image planes')
            for view in tabs.near_views.values():
                self.assertIsNotNone(view.item)
                self.assertTrue(tabs.near_tabs.isTabEnabled(tabs.near_tabs.indexOf(view)))
            old = tabs.parameters.toPlainText()
            (path/'near_field_xy.png').unlink()
            with self.assertRaisesRegex(ValueError, 'near-field figure: xy'):
                tabs.load(path)
            self.assertEqual(tabs.parameters.toPlainText(), old)
            report.pop('near_field_planes')
            report_path.write_text(json.dumps(report), encoding='utf-8')
            tabs.load(path)
            for view in tabs.near_views.values():
                self.assertIsNone(view.item)
                self.assertFalse(tabs.near_tabs.isTabEnabled(tabs.near_tabs.indexOf(view)))

    def test_component_comparison_tab_is_optional_and_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'imaging_report.json').write_text(json.dumps(dict(
                channels=[{'with': {'image_power_w': 1e-12}}]*3)), encoding='utf-8')
            tabs = ResultsTabs()
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            for _, name in tabs.FIGURES:
                self.assertTrue(pixmap.save(str(path/name)))
            tabs.load(path)
            self.assertFalse(tabs.isTabEnabled(tabs.indexOf(tabs.comparison_tabs)))
            run = path/'component_comparison'/'test'
            run.mkdir(parents=True)
            figures = dict(far='far.png', far_log='far_log.png',
                           image='image.png', image_log='image_log.png')
            for name in figures.values():
                self.assertTrue(pixmap.save(str(run/name)))
            (run/'comparison.json').write_text(json.dumps({'status': 'complete'}), encoding='utf-8')
            (run/'totals.txt').write_text('Particle 1 pW; haze 2 pW', encoding='utf-8')
            (path/'component_comparison_manifest.json').write_text(json.dumps(dict(
                directory='component_comparison/test', report='comparison.json',
                totals='totals.txt', figures=figures)), encoding='utf-8')
            tabs.load(path)
            self.assertTrue(tabs.isTabEnabled(tabs.indexOf(tabs.comparison_tabs)))
            self.assertIn('haze 2 pW', tabs.comparison_totals.toPlainText())
            old = tabs.comparison_totals.toPlainText()
            (run/'image_log.png').unlink()
            with self.assertRaisesRegex(ValueError, 'component comparison figure'):
                tabs.load(path)
            self.assertEqual(tabs.comparison_totals.toPlainText(), old)

    def test_process_command_and_cooperative_cancel(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root/'python.exe'
            python.touch()
            widget = FullwaveWidget()
            widget.controls.fields['illumination-fit'].setCurrentText('manual')
            widget.interpreter.setText(str(python))
            with patch('fullwave_ui.widget.ROOT', root), patch.object(widget.process, 'start') as start:
                widget.start_job()
                interpreter, command = start.call_args[0]
                self.assertEqual(interpreter, str(python))
                self.assertIn('fullwave_prototype.guarded_run', command)
                self.assertIn('--cancel-file', command)
                self.assertIn('--particle-k', command)
                self.assertFalse(widget.controls.isEnabled())
                with patch.object(widget.process, 'state', return_value=QProcess.Running):
                    widget.cancel_job()
                self.assertTrue(widget.cancel_file.exists())
                self.assertFalse(widget.cancel.isEnabled())
                widget.process_error(QProcess.FailedToStart)
                self.assertTrue(widget.controls.isEnabled())

    def test_component_comparison_uses_guarded_existing_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root/'python.exe'
            python.touch()
            result = root/'result'
            source = result/'scattering'
            source.mkdir(parents=True)
            (source/'report.json').write_text(json.dumps(dict(
                source_format_version=1, parameters=dict(incidence_deg=75))), encoding='utf-8')
            (source/'paired_fields.npz').touch()
            (result/'imaging_report.json').write_text(json.dumps(dict(
                channels=[{'with': {'image_power_w': 1e-12}}]*3,
                settings=asdict(ImagingSettings(numerical_apertures=(0, 0, .96),
                                                tilt_degrees=(45, -45, 0))))), encoding='utf-8')
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            widget = FullwaveWidget()
            for _, name in widget.results.FIGURES:
                self.assertTrue(pixmap.save(str(result/name)))
            self.assertTrue(widget.load_result(result))
            self.assertTrue(widget.compare_button.isEnabled())
            self.assertEqual([widget.controls.fields['na'+str(i)].value() for i in (1, 2, 3)],
                             [0, 0, .96])
            widget.controls.fields['na2'].setValue(.4125)
            widget.controls.fields['tilt2'].setValue(-30)
            widget.interpreter.setText(str(python))
            with patch('fullwave_ui.widget.ROOT', root), patch.object(widget.process, 'start') as start:
                widget.start_comparison()
                command = start.call_args[0][1]
                self.assertIn('component_compare', command)
                self.assertEqual(command[command.index('--result')+1], str(result.resolve()))
                settings_path = Path(command[command.index('--imaging-settings-json')+1])
                requested = json.loads(settings_path.read_text(encoding='utf-8'))
                self.assertEqual(requested['numerical_apertures'], [0, .4125, .96])
                self.assertEqual(requested['tilt_degrees'], [45, -30, 0])
                self.assertTrue(widget.active_report.name.startswith('comparison_resources_'))
                self.assertFalse(widget.compare_button.isEnabled())
                widget.process_error(QProcess.FailedToStart)
            widget.close()

    def test_window_keeps_scatmech_and_fem_separate(self):
        from scatmech_ui.app import ScatmechStandaloneWindow
        window = ScatmechStandaloneWindow()
        self.assertIs(window.tabs.widget(0), window.simulation_widget)
        self.assertIs(window.tabs.widget(1), window.fullwave_widget)
        self.assertTrue(window.fullwave_widget.wait_for_shutdown())
        window.close()

    def test_whole_beam_job_is_distinct_and_resource_guarded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root/'python.exe'
            python.touch()
            widget = FullwaveWidget()
            widget.controls.fields['detector'].setChecked(True)
            widget.interpreter.setText(str(python))
            widget.controls.fields['na1'].setValue(0)
            widget.whole_coverage.setValue(.995)
            with patch('fullwave_ui.whole_haze_job.application_root', return_value=root), \
                    patch.object(widget.process, 'start') as start:
                widget.start_whole_haze()
                command = start.call_args[0][1]
                self.assertIn('whole_haze', command)
                self.assertIn('--cancel-file', command)
                self.assertEqual(widget.job_kind, 'whole_haze')
                inputs = json.loads(Path(command[command.index('--input')+1]).read_text(encoding='utf-8'))
                self.assertEqual(inputs['coverage'], .995)
                self.assertEqual(inputs['imaging_settings']['numerical_apertures'][0], 0)
                self.assertFalse(inputs['imaging_settings']['detector']['enabled'])
                self.assertEqual(inputs['imaging_settings']['detector']['pixels_x'], 16)
                self.assertFalse(widget.whole_button.isEnabled())
            widget.process_error(QProcess.FailedToStart)
            widget.close()

    def test_whole_beam_result_tab_does_not_accept_a_fem_report(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            tabs = ResultsTabs()
            report = dict(status='complete', model='whole_beam_first_order_surface')
            (path/'whole_haze_report.json').write_text(json.dumps(report), encoding='utf-8')
            pixmap = QPixmap(80, 60)
            pixmap.fill()
            for name in tabs.whole_haze_views:
                self.assertTrue(pixmap.save(str(path/name)))
            tabs.load_whole_haze(path)
            old = tabs.whole_haze_details.toPlainText()
            report['model'] = 'FEM'
            (path/'whole_haze_report.json').write_text(json.dumps(report), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'first-order surface'):
                tabs.load_whole_haze(path)
            self.assertEqual(tabs.whole_haze_details.toPlainText(), old)


if __name__ == '__main__':
    unittest.main()
