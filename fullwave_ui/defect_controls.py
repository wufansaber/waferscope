"""Finite defect controls and a wafer-coordinate footprint preview."""
import math

from PyQt5.QtCore import QPointF, Qt, QSignalBlocker
from PyQt5.QtGui import QColor, QPainter, QPen, QPolygonF
from PyQt5.QtWidgets import QLabel, QWidget

from fullwave_prototype.defects import DEFAULTS, SHAPES, footprint_nm, metadata


class DefectPreview(QWidget):
    def __init__(self, controls):
        super().__init__()
        self.controls = controls
        self.setMinimumHeight(175)
        self.setMaximumHeight(200)
        self.setToolTip('Wafer XY view. line_x is along X; line_y is along Y before Z rotation. '
                        'Illumination propagates from -Y toward +Y.')

    def paintEvent(self, event):
        f = self.controls.fields
        p = {key: f[key.replace('_', '-')].currentText() if key == 'defect_shape'
             else f[key.replace('_', '-')].value() for key in DEFAULTS}
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        width, height = self.width(), self.height()
        center = QPointF(width*.55, height*.52)
        span = f['diameter-nm'].value() if p['defect_shape'] == 'sphere' else p['defect_length_nm']*1.42
        scale = min(width*.7, height*.66)/max(span, 1)
        painter.setPen(QPen(QColor('#167d9a'), 1))
        painter.drawLine(QPointF(15, height-30), QPointF(15, 40))
        painter.drawLine(QPointF(15, 40), QPointF(11, 48))
        painter.drawLine(QPointF(15, 40), QPointF(19, 48))
        painter.drawText(2, height-10, 'k_in')
        painter.drawText(width-25, height-10, '+X')
        painter.drawText(width-25, 20, '+Y')
        painter.setPen(QPen(QColor('#812d59'), 1.5))
        painter.setBrush(QColor('#d78ba6'))
        if p['defect_shape'] == 'sphere':
            radius = f['diameter-nm'].value()*scale/2
            painter.drawEllipse(center, radius, radius)
        else:
            polygon = QPolygonF([QPointF(center.x()+x*scale, center.y()-y*scale)
                                 for x, y in footprint_nm(p)])
            painter.drawPolygon(polygon)
        painter.end()


class DefectControls:
    def __init__(self, controls):
        self.controls = controls
        self._loading = False
        controls.choice('defect-shape', 'Defect shape', SHAPES)
        for name, label, value in (
                ('defect-length-nm', 'Line / arm length (nm)', 200),
                ('defect-width-nm', 'Line / arm width (nm)', 50),
                ('defect-height-nm', 'Defect height (nm)', 50)):
            controls.number(name, label, value, 1, 5000, 2)
        controls.number('defect-rotation-deg', 'Rotation about +Z (deg)', 0, -180, 180, 1)
        controls.fields['defect-shape'].setToolTip('Raised finite solid, not an etched groove. '
            'Cross: two equal orthogonal arms, overlap counted once. FEM only; not a spherical BV model.')
        controls.fields['defect-length-nm'].setToolTip('Total end-to-end length of each cross arm or finite line')
        controls.fields['domain-radius-um'].setToolTip('Automatically enlarged if the defect reaches the physical '
            'domain boundary, with 20% radial margin. Mesh sizes and FEM order are unchanged.')
        self.preview = DefectPreview(controls)
        controls.form.addRow('Defect XY footprint', self.preview)
        self.info = QLabel()
        self.info.setWordWrap(True)
        controls.form.addRow('Geometry', self.info)
        self._move_near_diameter()
        controls.fields['defect-shape'].currentTextChanged.connect(self.update)
        for name in ('diameter-nm', 'gap-nm', 'domain-radius-um')+tuple(
                key.replace('_', '-') for key in DEFAULTS if key != 'defect_shape'):
            controls.fields[name].valueChanged.connect(self.update)
        self.update()

    def _move_near_diameter(self):
        form, fields = self.controls.form, self.controls.fields
        target = form.getWidgetPosition(fields['diameter-nm'])[0]+1
        widgets = [fields[key.replace('_', '-')] for key in DEFAULTS]+[self.preview, self.info]
        for widget in widgets:
            row = form.takeRow(form.getWidgetPosition(widget)[0])
            form.insertRow(target, row.labelItem.widget(), row.fieldItem.widget())
            target += 1

    def update(self, *_):
        f = self.controls.fields
        sphere = f['defect-shape'].currentText() == 'sphere'
        f['diameter-nm'].setEnabled(sphere)
        for key in DEFAULTS:
            if key != 'defect_shape':
                f[key.replace('_', '-')].setEnabled(not sphere)
        p = dict(diameter_nm=f['diameter-nm'].value(), gap_nm=f['gap-nm'].value())
        p.update({key: f[key.replace('_', '-')].currentText() if key == 'defect_shape'
                  else f[key.replace('_', '-')].value() for key in DEFAULTS})
        info = metadata(p)
        error = (not sphere and p['defect_width_nm'] > p['defect_length_nm'])
        if not error and not self._loading and not getattr(self.controls, '_applying_preset', False):
            self.fit_domain(info['required_radius_um'])
        error |= f['domain-radius-um'].value() <= info['required_radius_um']
        placement = ' | bottom z = {:g} nm'.format(p['gap_nm'])
        if sphere:
            placement += ' | center z = {:g} nm'.format(p['gap_nm']+p['diameter_nm']/2)
        self.info.setText('V = {:.4g} um^3 | domain R > {:.4g} um{}{}'.format(
            info['expected_volume_um3'], info['required_radius_um'],
            placement,
            ' | Invalid dimensions / domain' if error else ''))
        self.info.setStyleSheet('color: #b02c30' if error else '')
        self.preview.update()

    def fit_domain(self, minimum):
        radius = self.controls.fields['domain-radius-um']
        if radius.value() > minimum:
            return
        step = 10.**(-radius.decimals())
        target = math.ceil(max(minimum*1.2, minimum+step)/step)*step
        # Keep this a domain correction, never a hidden mesh/order adjustment.
        with QSignalBlocker(radius):
            radius.setValue(min(target, radius.maximum()))

    def load(self, saved):
        self._loading = True
        try:
            for name, default in DEFAULTS.items():
                field = self.controls.fields[name.replace('_', '-')]
                value = saved.get(name, default)
                field.setCurrentText(value) if name == 'defect_shape' else field.setValue(value)
            self.update()
        finally:
            self._loading = False
