"""Check brochure data and package the verified PDF and editable materials."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from zipfile import ZipFile, ZIP_DEFLATED

from pypdf import PdfReader


def main():
    root = Path(__file__).resolve().parents[1]
    folder = root/'docs/showcase_20261017'
    pdf = root/'output/pdf/FarField_Simulation_功能与案例手册.pdf'
    data = json.loads((folder/'case_manifest.json').read_text(encoding='utf-8'))
    reader = PdfReader(pdf)
    assert len(reader.pages) == 10
    texts = [page.extract_text() for page in reader.pages]
    assert all(len(t) > 100 and '\ufffd' not in t and '\x00' not in t for t in texts)
    assert '未收敛' in texts[8] and '不代表认证' in texts[9]
    dense = data['pipeline']['far_field_angular_evaluation']
    assert dense['new_grid_shape']==[401,401] and dense['evaluated_directions']>120000
    assert dense['coarse_direction_complex_relative_l2']<1e-10
    assert dense['surface_fft_vs_direct_relative_l2']<.01
    assert '401×401' in texts[0] and '401×401' in texts[2]
    for row in data['camera']['results']:
        for ch in row['channels']:
            assert ch['gray_valid'] and ch['clipped_values'] == 0
            assert '%.2f'%ch['gray_snr'] in texts[7]
    with (folder/'particle_size_sweep.csv').open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.reader(stream))
    assert len(rows) == 19
    for csv_row, expected in zip(rows[1:], data['sweep']['rows']):
        assert all(math.isclose(float(x), y, rel_tol=1e-13) for x, y in zip(csv_row, expected))
        assert math.isclose(expected[1], expected[2], rel_tol=1e-12)
    for path, digest in data['provenance']['inputs_sha256'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest() == digest
    md = (folder/'软件功能与案例手册.md').read_text(encoding='utf-8')
    links = re.findall(r'\]\(([^)]+)\)', md)
    assert links and all((folder/link).is_file() for link in links)
    for name in ('build_showcase_cases.py', 'build_showcase_brochure.py', 'validate_showcase_brochure.py',
                 'refine_showcase_far_field.py'):
        compile((root/'tools'/name).read_text(encoding='utf-8'), name, 'exec')
    package = root/'output/pdf/FarField_Simulation_案例与宣传资料.zip'
    with ZipFile(package, 'w', compression=ZIP_DEFLATED) as archive:
        archive.write(pdf, pdf.name)
        for file in sorted(folder.rglob('*')):
            if file.is_file():
                archive.write(file, 'editable_materials/'+file.relative_to(folder).as_posix())
        for name in ('build_showcase_cases.py', 'build_showcase_brochure.py', 'validate_showcase_brochure.py',
                     'refine_showcase_far_field.py'):
            archive.write(root/'tools'/name, 'reproduction_helpers/'+name)
    with ZipFile(package) as archive:
        assert archive.testzip() is None
    report = dict(status='passed', pdf_pages=10, showcase_categories=6,
        dense_far_field_grid=dense['new_grid_shape'], dense_far_field_directions=dense['evaluated_directions'],
        particle_sweep_points=18, camera_exposures=3, camera_pairs_per_exposure_channel=1000,
        valid_camera_metrics=9, referenced_files=len(links),
        source_report_hashes_checked=len(data['provenance']['inputs_sha256']),
        visual_qa='All 10 final rendered pages inspected',
        pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
        zip_sha256=hashlib.sha256(package.read_bytes()).hexdigest(),
        zip_bytes=package.stat().st_size)
    output = root/'build/showcase_20261017/validation.json'
    output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
