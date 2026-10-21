"""Render the canonical Markdown operation manual as a linked, paginated PDF."""
import hashlib
from html import escape
import json
from pathlib import Path
import re

import mistune
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (CondPageBreak, Image, KeepTogether, LongTable,
    PageBreak, Paragraph, Preformatted, SimpleDocTemplate, Spacer, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'far_field_standalone/OPERATION_MANUAL_CN.md'
OUTPUT = ROOT/'output/pdf/WaferScope_中文操作手册.pdf'
EVIDENCE = ROOT/'build/manual_pdf_20261017'
WIDTH = A4[0]-96


def normalize(text):
    text = text.replace('⁻', '-').replace('−', '-').replace('–', '-').replace('—', '-')
    for old, new in zip('⁰¹²³⁴⁵⁶⁷⁸⁹', ('^0','^1','^2','^3','^4','^5','^6','^7','^8','^9')):
        text = text.replace(old, new)
    return text


def plain(children):
    return ''.join(token.get('raw', plain(token.get('children', []))) for token in children)


def slug(text):
    text = re.sub(r'[^\w\s-]', '', text.lower())
    return re.sub(r'\s+', '-', text.strip())


def inline(children, anchors):
    result = []
    for token in children:
        kind = token['type']
        if kind == 'text':
            result.append(escape(normalize(token['raw'])))
        elif kind == 'codespan':
            result.append('<font color="#225c62">'+escape(normalize(token['raw']))+'</font>')
        elif kind == 'strong':
            result.append('<b>'+inline(token['children'], anchors)+'</b>')
        elif kind == 'emphasis':
            result.append(inline(token['children'], anchors))
        elif kind == 'link':
            target = token['attrs']['url']
            label = inline(token['children'], anchors)
            if target.startswith('#') and target[1:] in anchors:
                result.append('<link href="#'+anchors[target[1:]]+'" color="#007f86">'+label+'</link>')
            elif target.startswith(('https://', 'http://')):
                result.append('<link href="'+escape(target, quote=True)+'" color="#007f86">'+label+'</link>')
            else:
                result.append(label)
        elif kind in ('softbreak', 'linebreak'):
            result.append('<br/>')
        else:
            raise ValueError('Unsupported inline token: '+kind)
    return ''.join(result)


class ManualDocument(SimpleDocTemplate):
    def afterFlowable(self, flowable):
        if getattr(flowable, 'manual_heading', None):
            level, text, key = flowable.manual_heading
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(text, key, level=level, closed=level > 0)
            if level == 0:
                self.notify('TOCEntry', (0, text, self.page, key))


def main():
    pdfmetrics.registerFont(TTFont('CN', 'C:/Windows/Fonts/msyh.ttc'))
    pdfmetrics.registerFont(TTFont('CNB', 'C:/Windows/Fonts/msyhbd.ttc'))
    pdfmetrics.registerFontFamily('CN', normal='CN', bold='CNB', italic='CN', boldItalic='CNB')
    body = ParagraphStyle('ManualBody', fontName='CN', fontSize=9.7, leading=15.1,
        textColor=colors.HexColor('#202529'), spaceAfter=7, wordWrap='CJK',
        allowWidows=0, allowOrphans=0)
    heading = ParagraphStyle('ManualChapter', parent=body, fontName='CNB', fontSize=14,
        leading=21, spaceBefore=16, spaceAfter=10, keepWithNext=True, textColor=colors.black)
    subheading = ParagraphStyle('ManualSubheading', parent=heading, fontSize=11,
        leading=17, spaceBefore=12, spaceAfter=6)
    caption = ParagraphStyle('ManualCaption', parent=body, fontSize=8.2, leading=12.5,
        textColor=colors.HexColor('#525b61'), spaceAfter=12)
    listing = ParagraphStyle('ManualList', parent=body, leftIndent=16, firstLineIndent=-16, spaceAfter=5)
    cell = ParagraphStyle('ManualCell', parent=body, fontSize=8.5, leading=12.5, spaceAfter=0)
    code = ParagraphStyle('ManualCode', parent=body, fontSize=8, leading=12,
        backColor=colors.HexColor('#f1f5f5'), borderPadding=8, spaceBefore=4, spaceAfter=11)
    title = ParagraphStyle('ManualTitle', parent=heading, fontSize=27, leading=38, spaceAfter=20)
    contents = ParagraphStyle('ManualContents', parent=body, leading=22, fontSize=10,
        leftIndent=0, firstLineIndent=0, rightIndent=25, spaceBefore=3)
    source_text = SOURCE.read_text(encoding='utf-8')
    ast = mistune.create_markdown(renderer='ast', plugins=['table'])(source_text)
    headings = [t for t in ast if t['type'] == 'heading' and t['attrs']['level'] in (2,3)
                and plain(t['children']) != '目录']
    anchors = {slug(plain(t['children'])): 'section-%d'%i for i, t in enumerate(headings)}
    heading_keys = {plain(t['children']): 'section-%d'%i for i, t in enumerate(headings)}
    story = [Spacer(1,55), Paragraph('WaferScope<br/>中文操作手册', title),
             Paragraph('模型入口  参数设置  成像与相机  数据保存', subheading), Spacer(1,20)]
    toc = TableOfContents(); toc.levelStyles = [contents]; toc.dotsMinLevel = 0
    toc_pending = False
    images = []

    def make_table(token):
        head = token['children'][0]['children']
        rows = [head]+[row['children'] for row in token['children'][1]['children']]
        values = [[Paragraph(inline(c['children'], anchors), cell) for c in row] for row in rows]
        n = len(head)
        ratios = ([.35,.65] if n == 2 else [.30,.27,.43] if n == 3 else [1/n]*n)
        table = LongTable(values, colWidths=[WIDTH*r for r in ratios], repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e3ecee')),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white, colors.HexColor('#f7f8f8')]),
            ('GRID',(0,0),(-1,-1),.45,colors.HexColor('#d9d9d9')),
            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ('LEFTPADDING',(0,0),(-1,-1),7), ('RIGHTPADDING',(0,0),(-1,-1),7),
            ('TOPPADDING',(0,0),(-1,-1),6), ('BOTTOMPADDING',(0,0),(-1,-1),6)]))
        return [table, Spacer(1,10)]

    i = 0
    while i < len(ast):
        token = ast[i]; i += 1
        kind = token['type']
        if kind == 'blank_line':
            continue
        if kind == 'heading':
            level = token['attrs']['level']; text = plain(token['children'])
            if level == 1:
                continue
            if text == '目录':
                story += [Spacer(1,15), Paragraph('阅读提示', subheading),
                    Paragraph('先按第 2 节选择计算入口，再按第 17 节运行常用任务。'
                              '结果解释与模型边界请结合第 10、14、18 节。', body),
                    PageBreak(), Paragraph('目录', heading), toc, PageBreak()]
                toc_pending = True
                continue
            if level == 2:
                number = int(text.split('.')[0])
                story.append(CondPageBreak(310 if number in (4,5,8,11,15,16,17,18) else 210))
            p = Paragraph(inline(token['children'], anchors), heading if level == 2 else subheading)
            following = next((t for t in ast[i:] if t['type'] != 'blank_line'), None)
            if following and following['type'] == 'table':
                # Keep the heading with initial rows, without moving a whole long table.
                story.append(CondPageBreak(140))
                p.keepWithNext = False
            p.manual_heading = (0 if level == 2 else 1, normalize(text), heading_keys[text])
            story.append(p)
        elif kind == 'paragraph':
            children = token['children']
            if len(children) == 1 and children[0]['type'] == 'image':
                path = SOURCE.parent/children[0]['attrs']['url']
                with PILImage.open(path) as im:
                    pw, ph = im.size
                scale = min(WIDTH/pw, (300 if 'controls' in path.name else 235)/ph)
                figure = Image(str(path), width=pw*scale, height=ph*scale)
                figure.hAlign = 'CENTER'
                group = [figure, Spacer(1,6)]
                j = i
                while j < len(ast) and ast[j]['type'] == 'blank_line':
                    j += 1
                if j < len(ast) and ast[j]['type'] == 'paragraph' and plain(ast[j]['children']).startswith('图 '):
                    group.append(Paragraph(inline(ast[j]['children'], anchors), caption))
                    i = j+1
                story.append(KeepTogether(group)); images.append(str(path.relative_to(ROOT)))
            else:
                if plain(children).lstrip().startswith('|'):
                    raise ValueError('Malformed Markdown table; escape literal cell separators')
                story.append(Paragraph(inline(children, anchors), body))
        elif kind == 'list':
            if toc_pending:
                toc_pending = False
                continue
            ordered = token.get('attrs',{}).get('ordered',False)
            for index, item in enumerate(token['children'],1):
                for block in item['children']:
                    if block['type'] not in ('block_text','paragraph'):
                        raise ValueError('Unsupported list content: '+block['type'])
                    prefix = str(index)+'. ' if ordered else '- '
                    story.append(Paragraph(prefix+inline(block['children'], anchors), listing))
        elif kind == 'table':
            story += make_table(token)
        elif kind == 'block_code':
            # Preserve literal syntax; line wrapping is visual only in this PDF.
            lines = []
            for line in normalize(token['raw']).splitlines():
                while len(line) > 85:
                    lines.append(line[:85]); line = '    '+line[85:]
                lines.append(line)
            story.append(Preformatted('\n'.join(lines), code))
        else:
            raise ValueError('Unsupported block: '+kind)

    def frame(c, doc):
        c.saveState(); c.setFont('CN',8); c.setFillColor(colors.HexColor('#59646a'))
        if doc.page > 1:
            c.drawString(48,A4[1]-30,'WaferScope 中文操作手册')
        c.drawString(48,25,'2026-10-21  |  WaferScope 源码与分文件构建')
        c.drawRightString(A4[0]-48,25,'第 %d 页'%doc.page)
        c.restoreState()
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    doc = ManualDocument(str(OUTPUT), pagesize=A4, rightMargin=48,leftMargin=48,
        topMargin=48,bottomMargin=46, title='WaferScope 中文操作手册',
        author='WaferScope', subject='操作步骤 参数说明 模型边界与常见问题')
    doc.multiBuild(story, onFirstPage=frame, onLaterPages=frame)
    EVIDENCE.mkdir(parents=True,exist_ok=True)
    evidence = dict(source=str(SOURCE.relative_to(ROOT)),
        markdown_sha256=hashlib.sha256(source_text.encode('utf-8')).hexdigest(),
        pdf_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        chapters=sum(t['attrs']['level']==2 for t in headings),
        heading_count=len(headings), figures=images,
        relative_document_links='Kept as labels in PDF; navigable local links remain in Markdown.',
        image_scope='Current Qt UI, saved solver/camera results; no new electromagnetic solves.')
    (EVIDENCE/'build_report.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2),encoding='utf-8')
    print(str(OUTPUT)); print('Main chapters:',evidence['chapters'])


if __name__ == '__main__':
    main()
