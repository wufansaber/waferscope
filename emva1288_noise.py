"""EMVA 1288 Linear core equations; not a camera certification implementation."""
import numpy as np

STANDARD_URL = 'https://www.emva.org/wp-content/uploads/EMVA1288Linear_4.0Release.pdf'


def noise_metrics(signal, background, parameters, nominal_signal=None):
    """Conditional temporal noise and uncorrected spatial+temporal SNR predictions."""
    gain = parameters.gain_e_dn
    temporal_e2 = signal+background+parameters.read_noise_e**2+gain**2/12
    nominal_signal = signal if nominal_signal is None else nominal_signal
    spatial_e2 = (parameters.dsnu_e**2+(parameters.prnu_fraction*nominal_signal)**2
                  +(parameters.dark_current_std_e_s*parameters.exposure_ms*1e-3)**2)
    total_e2 = nominal_signal+background+parameters.read_noise_e**2+gain**2/12+spatial_e2
    temporal_snr = np.divide(signal, np.sqrt(temporal_e2), out=np.zeros_like(signal), where=temporal_e2 > 0)
    total_snr = np.divide(nominal_signal, np.sqrt(total_e2), out=np.zeros_like(signal), where=total_e2 > 0)
    return dict(temporal_variance_dn2=temporal_e2/gain**2,
                digital_ideal_snr=temporal_snr, emva_total_snr=total_snr,
                digital_ideal_roi_snr=float(signal.sum()/np.sqrt(temporal_e2.sum())),
                emva_total_roi_snr=float(nominal_signal.sum()/np.sqrt(total_e2.sum())),
                emva_system_gain_dn_e=1/gain, emva_quantization_variance_dn2=1/12,
                noise_definition='EMVA 1288 Release 4.0 Linear: temporal Eqs. 13-15; '
                    'DSNU/PRNU and total SNR Eqs. 66-69. Independent-pixel spatial model. '
                    'Total SNR uses nominal average QE, before pixel calibration. '
                    'DCNU is an optional extension, not part of the standard. SNR is pre-clipping; '
                    'quantization variance is an approximation for adequately resolved noise.',
                quantization_noise_resolved=bool(parameters.read_noise_e/gain >= .4),
                noise_standard_url=STANDARD_URL)


def pair_statistics(first, second):
    """Release 4 Eq. 18 removes the inter-frame mean offset before variance."""
    a, b = np.asarray(first, float), np.asarray(second, float)
    if a.shape != b.shape or a.size < 2 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Two equal finite image arrays required')
    difference = a-b
    return dict(mean_dn=float((a.mean()+b.mean())/2),
                temporal_variance_dn2=float(np.mean((difference-difference.mean())**2)/2),
                spatial_covariance_dn2=float(np.mean((a-a.mean())*(b-b.mean()))))
