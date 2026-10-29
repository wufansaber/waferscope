"""Explicit physical, observation and ideal-relay settings for periodic scenes."""
from dataclasses import asdict, dataclass, field

import numpy as np

from rcwa_prototype.geometry import Scene


@dataclass(frozen=True)
class SceneConfig:
    scene: Scene = field(default_factory=lambda: Scene(
        harmonics_x=3, harmonics_y=16, slices=4, gap_um=0., grating_method='analytic_li'))
    particle: bool = True
    polarization: str = 'p'
    numerical_apertures: tuple = (.6, .6, .6)
    tilts_deg: tuple = (45., -45., 0.)
    magnification: float = 20.
    transmission: float = 1.
    analyzer: str = 'none'
    defocus_um: float = 0.
    illumination_power_w: float = 1.
    normalization_length_um: float = 2.
    normalization_width_um: float = 2.
    near_samples: int = 25
    near_enabled: bool = True
    near_xy_z_um: float = .02
    near_air_height_um: float = .15
    image_samples: int = 128

    def validate(self):
        self.scene.validate()
        s = self.scene
        if s.harmonics_x is None or s.grating_method != 'analytic_li':
            raise ValueError('This scene requires an explicit rectangular basis and analytic stripe rule')
        if self.particle and min(s.lx_um, s.ly_um) <= s.diameter_um:
            raise ValueError('The particle must not touch periodic copies; enlarge the cell')
        if not isinstance(self.near_enabled, bool) or not isinstance(self.particle, bool) or self.polarization not in ('s', 'p'):
            raise ValueError('Boolean particle flag and s/p illumination required')
        if len(self.numerical_apertures) != 3 or len(self.tilts_deg) != 3:
            raise ValueError('Exactly three independent collection channels required')
        values = (*self.numerical_apertures, *self.tilts_deg, self.magnification,
            self.transmission, self.defocus_um, self.illumination_power_w,
            self.normalization_length_um, self.normalization_width_um,
            self.near_xy_z_um, self.near_air_height_um)
        if not np.isfinite(values).all():
            raise ValueError('All optical and observation settings must be finite')
        if any(not 0 <= n < 1 for n in self.numerical_apertures):
            raise ValueError('Each air NA must be in [0,1); zero disables the channel')
        if any(abs(t) >= 90 for t in self.tilts_deg):
            raise ValueError('Collection tilts must be strictly between -90 and 90 degrees')
        if self.magnification < 1 or not 0 <= self.transmission <= 1:
            raise ValueError('Magnification >=1 and power transmission in [0,1] required')
        if min(self.illumination_power_w, self.normalization_length_um,
               self.normalization_width_um, self.near_air_height_um) <= 0:
            raise ValueError('Positive illumination normalization and observation height required')
        if self.analyzer not in ('none', 'x', 'y', '45', '135'):
            raise ValueError('Unknown linear analyzer')
        for name, lower, upper in (('near_samples', 5, 129), ('image_samples', 16, 1024)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
                raise ValueError(name+' is outside the supported integer range')
        # Intensities have twice the field bandwidth. Avoid wraparound in pixel integration.
        if self.image_samples < 4*max(s.harmonics_x, s.harmonics_y)+1:
            raise ValueError('Image samples must be >= 4*max(Mx,My)+1 to resolve intensity bandwidth')

    @property
    def irradiance_w_um2(self):
        return self.illumination_power_w/(self.normalization_length_um*self.normalization_width_um
            * np.cos(np.deg2rad(self.scene.incidence_deg)))

    @property
    def incident_cell_power_w(self):
        return self.illumination_power_w*self.scene.lx_um*self.scene.ly_um/(
            self.normalization_length_um*self.normalization_width_um)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        data = dict(data)
        data['scene'] = Scene(**data['scene'])
        for name in ('numerical_apertures', 'tilts_deg'):
            if name in data:
                data[name] = tuple(data[name])
        result = cls(**data)
        result.validate()
        return result
