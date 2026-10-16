"""Compact numeric comparison for three independent paired camera receivers."""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QAbstractItemView, QHeaderView, QLabel, QTableWidget,
                            QTableWidgetItem, QVBoxLayout, QWidget)


class CameraMetrics(QWidget):
    ROWS = (
        ('Gray-frame SNR status', 'gray_status'),
        ('Gray-frame ROI SNR', 'gray_snr'),
        ('Gray-frame SNR (dB)', 'gray_snr_db'),
        ('Acquired frame pairs', 'gray_frames'),
        ('Gray ROI contrast (DN)', 'gray_contrast_dn'),
        ('Gray ROI noise (DN RMS)', 'gray_noise_dn'),
        ('Status', 'detection_status'),
        ('Defect contrast / ROI (e-)', 'defect_contrast_roi_e'),
        ('Haze / ROI (e-)', 'haze_roi_e'),
        ('Detection noise / ROI (e- RMS)', 'detection_noise_roi_e'),
        ('Predicted ROI detection SNR', 'detection_roi_snr'),
        ('Predicted known-template SNR', 'detection_template_snr'),
        ('Sample photon variance (e-^2)', 'sample_shot_roi_variance_e2'),
        ('Reference haze variance (e-^2)', 'reference_haze_roi_variance_e2'),
        ('Camera variance (e-^2)', 'camera_roi_variance_e2'),
        ('Extra haze mismatch variance (e-^2)', 'haze_mismatch_roi_variance_e2'),
        ('Invalid pixels', 'detection_invalid_pixels'),
    )

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout(self)
        self.mode = QLabel()
        self.mode.setWordWrap(True)
        layout.addWidget(self.mode)
        self.table = QTableWidget(len(self.ROWS), 3)
        self.table.setHorizontalHeaderLabels(('CH1', 'CH2', 'CH3'))
        self.table.setVerticalHeaderLabels([label for label, _ in self.ROWS])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        layout.addWidget(self.table)

    def clear(self):
        self.mode.clear()
        self.table.clearContents()

    def load(self, channels, parameters):
        self.mode.setText('Same-surface reference | {} frames | residual haze RMS {:.3g}%'.format(
            parameters.get('reference_frames', 1), 100*parameters.get('haze_mismatch_fraction', 0)))
        notes = [c.get('detector', {}).get('camera_profile_note') for c in channels]
        if any(notes):
            self.mode.setToolTip(next(n for n in notes if n))
        for col, channel in enumerate(channels):
            detector = channel['detector']
            for row, (_, key) in enumerate(self.ROWS):
                value = detector.get(key)
                text = 'Invalid' if value is None else '{:.5g}'.format(value) if isinstance(value, (int, float)) else str(value)
                if key.startswith('gray_') and key not in detector:
                    text = 'Not acquired'
                if not channel.get('enabled', True) and key in ('gray_status', 'gray_snr', 'gray_snr_db'):
                    text = 'NA=0 (off)'
                if key == 'detection_status' and not channel.get('enabled', True):
                    text = 'NA=0 (off)'
                elif key == 'detection_status':
                    text = 'Valid' if detector.get('detection_roi_valid') else 'Invalid / clipped'
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignCenter)
                item.setToolTip(detector.get('gray_scope', '') if key.startswith('gray_') else detector.get('detection_scope', ''))
                if key.startswith('gray_') and key in detector and not detector.get('gray_valid'):
                    item.setForeground(Qt.red)
                if not detector.get('detection_roi_valid') and key in ('detection_status', 'detection_roi_snr', 'detection_template_snr'):
                    item.setForeground(Qt.red)
                self.table.setItem(row, col, item)
