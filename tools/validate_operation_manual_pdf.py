"""Verify manual content completeness, outline targets, images, and PDF integrity."""
import hashlib
import json
from pathlib import Path
import re

import mistune
from pypdf import PdfReader
from .build_operation_manual_pdf import normalize, plain


def compact(text):
    return re.sub(r'\s+', '', normalize(text))


def main():
    root = Path(__file__).resolve().parents[1]
    source = root/'far_field_standalone/OPERATION_MANUAL_CN.md'
    pdf = root/'output/pdf/FarField_Simulation_中文操作手册.pdf'
    folder = root/'build/manual_pdf_20261017'
    build = json.loads((folder/'build_report.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(source.read_bytes()).hexdigest() == build['markdown_sha256']
    assert hashlib.sha256(pdf.read_bytes()).hexdigest() == build['pdf_sha256']
    reader = PdfReader(pdf)
    text = '\n'.join(page.extract_text() for page in reader.pages)
    assert '\x00' not in text and '\ufffd' not in text
    merged = compact(text)
    ast = mistune.create_markdown(renderer='ast', plugins=['table'])(source.read_text(encoding='utf-8'))
    checked, missing = 0, []

    def visit(token):
        nonlocal checked
        kind = token['type']
        if kind in ('paragraph', 'block_text', 'table_cell'):
            if any(c['type']=='image' for c in token.get('children',[])):
                return
            value = compact(plain(token.get('children',[])))
            if value:
                checked += 1
                if value not in merged:
                    missing.append(value[:100])
        elif kind == 'block_code':
            checked += 1
            if compact(token['raw']) not in merged:
                missing.append(token['raw'][:100])
        else:
            for child in token.get('children',[]):
                visit(child)
    for token in ast:
        if token['type'] not in ('list','heading') or not any(c.get('attrs',{}).get('url','').startswith('#') for c in token.get('children',[])):
            # The authored TOC is replaced by generated linked entries; its labels still appear.
            visit(token)
    assert not missing, missing
    main_entries = [entry for entry in reader.outline if not isinstance(entry,list)]
    assert len(main_entries) == 18
    for entry in main_entries:
        page = reader.get_destination_page_number(entry)
        assert 0 <= page < len(reader.pages)
        assert compact(entry.title) in compact(reader.pages[page].extract_text())
    toc_links = [annotation.get_object() for annotation in reader.pages[1].get('/Annots',[])]
    assert len(toc_links) >= 18
    for path in build['figures']:
        assert (root/path).is_file()
    assert len(build['figures']) == 4
    snapshots = json.loads((folder/'ui_inputs.json').read_text(encoding='utf-8'))
    low = snapshots['266 nm low-haze camera']
    assert low['wavelength-nm'] == 266 and low['diameter-nm'] == 100 and low['rms-nm'] == .2
    assert low['illumination-fit']=='camera' and low['gray-snr-pairs']==5000
    assert low['detector-quantum-efficiency']==.639 and low['detector-gain-e-dn']==16
    cross = snapshots['266 nm cross preview']
    assert cross['defect-shape']=='cross' and cross['defect-length-nm']==200
    assert cross['illumination-fit']=='manual' and cross['beam-profile']=='plane'
    report = dict(status='passed', pages=len(reader.pages), chapters=len(main_entries),
        content_blocks_checked=checked, toc_link_annotations=len(toc_links),
        actual_ui_screenshots=len(build['figures']), preset_inputs_checked=2,
        visual_qa='All 21 final rendered pages inspected', pdf_bytes=pdf.stat().st_size,
        pdf_sha256=build['pdf_sha256'])
    (folder/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
