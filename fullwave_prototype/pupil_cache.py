"""Content-addressed complex pupils; optical/camera changes do not rerun radiation."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
from zipfile import BadZipFile
import numpy as np


class PupilCache:
    def __init__(self, directory, data, parameters):
        self.directory = Path(directory)
        digest = hashlib.sha256()
        for key in ('source_points_with', 'source_polarization_weights_with',
                    'source_points_without', 'source_polarization_weights_without'):
            array = np.ascontiguousarray(data[key])
            digest.update(str((key, array.shape, array.dtype.str)).encode('ascii'))
            digest.update(array.tobytes())
        digest.update(repr((parameters.wavelength_nm, parameters.eps_substrate)).encode('ascii'))
        self.source_digest = digest.hexdigest()

    def path(self, na, tilt, samples):
        payload = json.dumps([2, self.source_digest, float(na), float(tilt), int(samples)])
        return self.directory/(hashlib.sha256(payload.encode('ascii')).hexdigest()+'.npz')

    def load(self, na, tilt, samples):
        path = self.path(na, tilt, samples)
        try:
            with np.load(path, allow_pickle=False) as saved:
                values = {name: saved[name] for name in ('with', 'without')}
                if any(a.shape != (samples, samples, 2) or not np.isfinite(a).all() for a in values.values()):
                    return None
                return values
        except (OSError, ValueError, KeyError, EOFError, BadZipFile):
            return None

    def save(self, na, tilt, samples, values):
        path = self.path(na, tilt, samples)
        temporary = None
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=self.directory, suffix='.npz', delete=False) as stream:
                temporary = Path(stream.name)
                np.savez_compressed(stream, **values)
            os.replace(str(temporary), str(path))
        except OSError as exc:
            print('Pupil cache unavailable; result remains valid: {}'.format(exc), flush=True)
        finally:
            if temporary is not None and temporary.exists():
                try:
                    temporary.unlink()
                except OSError:
                    pass
