from __future__ import annotations

import csv
import math
import os
import subprocess
from pathlib import Path

from PIL import Image
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_ALIGN_VERTICAL, WD_CELL_VERTICAL_ALIGNMENT, WD_ROW_HEIGHT_RULE, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph


ROOT = Path(__file__).resolve().parents[1]
OUT_TABLES = ROOT / "tables"
OUT_FIGURES = ROOT / "figures"
OUT_MANUSCRIPT = ROOT / "manuscript"
SOURCE_TABLES = ROOT / "outputs" / "tables"


def read_csv(name: str) -> list[dict[str, str]]:
    with (SOURCE_TABLES / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def set_cell_text(cell, text: str, *, bold: bool = False, size: float = 8.2,
                  align=WD_ALIGN_PARAGRAPH.LEFT, italic: bool = False):
    cell.text = ""
    p = cell.paragraphs[0]
    p.alignment = align
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.line_spacing = 1.0
    r = p.add_run(str(text))
    r.font.name = "Arial"
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic
    r.font.color.rgb = RGBColor(0, 0, 0)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def set_cell_width(cell, width_inches: float):
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    tcW = tcPr.first_child_found_in("w:tcW")
    if tcW is None:
        tcW = OxmlElement("w:tcW")
        tcPr.append(tcW)
    tcW.set(qn("w:w"), str(int(width_inches * 1440)))
    tcW.set(qn("w:type"), "dxa")


def set_table_borders(table, *, header_bottom: bool = True, final_bottom: bool = True):
    # Submission-style horizontal rules only; no cell shading and no vertical borders.
    for row_idx, row in enumerate(table.rows):
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            shd = tcPr.find(qn("w:shd"))
            if shd is not None:
                tcPr.remove(shd)
            borders = tcPr.find(qn("w:tcBorders"))
            if borders is None:
                borders = OxmlElement("w:tcBorders")
                tcPr.append(borders)
            for edge in ("left", "right", "insideH", "insideV"):
                el = borders.find(qn(f"w:{edge}"))
                if el is None:
                    el = OxmlElement(f"w:{edge}")
                    borders.append(el)
                el.set(qn("w:val"), "nil")
            if row_idx == 0:
                for edge in ("top", "bottom"):
                    el = borders.find(qn(f"w:{edge}"))
                    if el is None:
                        el = OxmlElement(f"w:{edge}")
                        borders.append(el)
                    el.set(qn("w:val"), "single")
                    el.set(qn("w:sz"), "8" if edge == "top" else "6")
                    el.set(qn("w:color"), "000000")
            if final_bottom and row_idx == len(table.rows) - 1:
                el = borders.find(qn("w:bottom"))
                if el is None:
                    el = OxmlElement("w:bottom")
                    borders.append(el)
                el.set(qn("w:val"), "single")
                el.set(qn("w:sz"), "8")
                el.set(qn("w:color"), "000000")


def prevent_row_split(row):
    trPr = row._tr.get_or_add_trPr()
    cant_split = OxmlElement("w:cantSplit")
    trPr.append(cant_split)


def repeat_header(row):
    trPr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    trPr.append(tbl_header)


def add_page_number(section):
    footer = section.footer
    p = footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    r = p.add_run("Page ")
    r.font.name = "Arial"
    r.font.size = Pt(8)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    p._p.append(fld)


def setup_doc(*, landscape_page: bool = True) -> Document:
    doc = Document()
    section = doc.sections[0]
    if landscape_page:
        section.orientation = WD_ORIENT.LANDSCAPE
        section.page_width = Inches(11)
        section.page_height = Inches(8.5)
        section.top_margin = Inches(0.55)
        section.bottom_margin = Inches(0.55)
        section.left_margin = Inches(0.60)
        section.right_margin = Inches(0.60)
    else:
        section.page_width = Inches(8.5)
        section.page_height = Inches(11)
        section.top_margin = Inches(0.75)
        section.bottom_margin = Inches(0.75)
        section.left_margin = Inches(0.85)
        section.right_margin = Inches(0.85)
    section.header_distance = Inches(0.25)
    section.footer_distance = Inches(0.3)
    add_page_number(section)

    normal = doc.styles["Normal"]
    normal.font.name = "Arial"
    normal.font.size = Pt(9)
    normal.paragraph_format.space_after = Pt(3)
    normal.paragraph_format.line_spacing = 1.0
    for style_name in ("Title", "Heading 1", "Heading 2"):
        style = doc.styles[style_name]
        style.font.name = "Arial"
        style.font.color.rgb = RGBColor(0, 0, 0)
    doc.styles["Title"].font.size = Pt(14)
    doc.styles["Title"].font.bold = True
    doc.styles["Heading 1"].font.size = Pt(11)
    doc.styles["Heading 1"].font.bold = True
    doc.styles["Heading 2"].font.size = Pt(9)
    doc.styles["Heading 2"].font.bold = True
    return doc


def add_title(doc: Document, title: str):
    p = doc.add_paragraph(style="Title")
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_after = Pt(5)
    p.add_run(title)


def add_note(doc: Document, text: str):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.line_spacing = 1.0
    r = p.add_run(text)
    r.font.name = "Arial"
    r.font.size = Pt(7.5)


def add_table(doc: Document, headers: list[str], rows: list[list[str]], widths: list[float],
              *, font_size: float = 8.2, alignments=None):
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    table.autofit = False
    table.allow_autofit = False
    hdr = table.rows[0]
    repeat_header(hdr)
    prevent_row_split(hdr)
    for idx, (cell, header, width) in enumerate(zip(hdr.cells, headers, widths)):
        set_cell_width(cell, width)
        set_cell_text(cell, header, bold=True, size=font_size,
                      align=(alignments[idx] if alignments else WD_ALIGN_PARAGRAPH.LEFT))
    for values in rows:
        row = table.add_row()
        prevent_row_split(row)
        for idx, (cell, value, width) in enumerate(zip(row.cells, values, widths)):
            set_cell_width(cell, width)
            set_cell_text(cell, value, size=font_size,
                          align=(alignments[idx] if alignments else WD_ALIGN_PARAGRAPH.LEFT))
    set_table_borders(table)
    return table


def add_section_row(table, text: str, cols: int, *, size: float = 8.2):
    row = table.add_row()
    prevent_row_split(row)
    merged = row.cells[0]
    for idx in range(1, cols):
        merged = merged.merge(row.cells[idx])
    set_cell_text(row.cells[0], text, bold=True, size=size)
    # A single thin separator above section headings, still without shading.
    tcPr = row.cells[0]._tc.get_or_add_tcPr()
    borders = tcPr.find(qn("w:tcBorders"))
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        tcPr.append(borders)
    top = OxmlElement("w:top")
    top.set(qn("w:val"), "single")
    top.set(qn("w:sz"), "4")
    top.set(qn("w:color"), "777777")
    borders.append(top)


def format_p(value: str) -> str:
    x = float(value)
    if x < 0.001:
        exponent = math.floor(math.log10(x))
        mantissa = x / (10 ** exponent)
        return f"{mantissa:.2f} x 10^{exponent}"
    return f"{x:.3f}"


def build_table1():
    data = read_csv("29_profile_baseline_characteristics.csv")
    # Avoid duplicating the identical metastatic-status/cancer-phenotype rows.
    data = [r for r in data if r["variable"] != "Metastatic status"]
    doc = setup_doc(landscape_page=True)
    add_title(doc, "Table 1. Baseline characteristics of the MIMIC-IV discovery cohort by inflammatory-metabolic stress profile")
    headers = ["Characteristic", "Overall\n(n = 4,380)", "Profile 1: lower stress\n(n = 3,197)", "Profile 2: higher stress\n(n = 1,183)"]
    widths = [3.85, 1.72, 1.83, 1.83]
    table = doc.add_table(rows=1, cols=4)
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    table.autofit = False
    table.allow_autofit = False
    for idx, (cell, header, width) in enumerate(zip(table.rows[0].cells, headers, widths)):
        set_cell_width(cell, width)
        set_cell_text(cell, header, bold=True, size=8.0,
                      align=WD_ALIGN_PARAGRAPH.LEFT if idx == 0 else WD_ALIGN_PARAGRAPH.CENTER)
    repeat_header(table.rows[0])
    prevent_row_split(table.rows[0])
    current_section = None
    for r in data:
        if r["section"] != current_section:
            current_section = r["section"]
            add_section_row(table, current_section, 4, size=8.0)
        label = r["variable"]
        if r["level"]:
            label = f"{label}, {r['level']}"
        if r["unit"]:
            label = f"{label}, {r['unit']}"
        row = table.add_row()
        prevent_row_split(row)
        vals = [label, r["overall"], r["profile_1_lower_stress"], r["profile_2_higher_stress"]]
        for idx, (cell, value, width) in enumerate(zip(row.cells, vals, widths)):
            set_cell_width(cell, width)
            set_cell_text(cell, value, size=7.8,
                          align=WD_ALIGN_PARAGRAPH.LEFT if idx == 0 else WD_ALIGN_PARAGRAPH.CENTER)
    set_table_borders(table)
    add_note(doc, "Data are median (IQR) for continuous variables and n (%) for categorical variables. Descriptive medians and quartiles were averaged across the 10 imputed data sets; they are not Rubin-rule estimates. Profile membership is the outcome-blind consensus assignment across the 10 aligned solutions; 26 ties were resolved using the higher mean posterior probability.")
    add_note(doc, "Laboratory values are the first-24-hour measurements used in the prespecified analysis. CRRT was unavailable in the frozen analysis data and was not reconstructed. Abbreviations: CRRT, continuous renal replacement therapy; IQR, interquartile range; SAPS-II, Simplified Acute Physiology Score II; SOFA, Sequential Organ Failure Assessment; WBC, white blood cell count.")
    path = OUT_TABLES / "Table1_baseline_characteristics.docx"
    doc.save(path)
    return path


def build_table2():
    selection = read_csv("27_lpa_profile_selection.csv")
    profiles = read_csv("28_lpa_profile_characteristics.csv")
    one_per_k = {}
    for r in selection:
        one_per_k.setdefault(int(r["profiles"]), r)
    rows_a = []
    for k in range(2, 7):
        r = one_per_k[k]
        rows_a.append([
            str(k), f"{float(r['mean_aic']):,.2f}", f"{float(r['mean_bic']):,.2f}",
            f"{float(r['mean_entropy']):.3f}", f"{float(r['minimum_profile_percent_across_candidate']):.2f}",
            f"{float(r['max_profile_proportion_range_pp']):.2f}",
            "Yes" if r["structure_stable"].upper() == "TRUE" else "No",
            "Yes" if r["admissible_for_later_selection"].upper() == "TRUE" else "No",
            "Yes" if r["recommended"].upper() == "TRUE" else "No",
        ])
    doc = setup_doc(landscape_page=True)
    add_title(doc, "Table 2. Latent profile model selection and standardized characteristics of the selected solution")
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(3)
    rr = p.add_run("Panel A. Candidate model comparison across 10 imputed data sets")
    rr.bold = True
    rr.font.name = "Arial"
    rr.font.size = Pt(9)
    add_table(doc,
              ["Profiles", "Mean AIC", "Mean BIC", "Mean entropy", "Minimum profile, %", "Maximum proportion range, pp", "Stable", "Admissible", "Selected"],
              rows_a, [0.58, 1.05, 1.08, 0.86, 1.12, 1.28, 0.70, 0.78, 0.72], font_size=7.4,
              alignments=[WD_ALIGN_PARAGRAPH.CENTER] * 9)

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(5)
    p.paragraph_format.space_after = Pt(3)
    rr = p.add_run("Panel B. Mean standardized variable values in the selected two-profile solution")
    rr.bold = True
    rr.font.name = "Arial"
    rr.font.size = Pt(9)
    vars_order = ["WBC", "platelet", "lactate", "glucose", "creatinine", "bicarbonate"]
    rows_b = []
    for var in vars_order:
        r1 = next(r for r in profiles if r["profile"] == "1" and r["variable"] == var)
        r2 = next(r for r in profiles if r["profile"] == "2" and r["variable"] == var)
        rows_b.append([var if var == "WBC" else var.capitalize(),
                       f"{float(r1['standardized_mean']):.3f}", f"{float(r2['standardized_mean']):.3f}"])
    add_table(doc, ["Variable", "Profile 1: lower stress", "Profile 2: higher stress"], rows_b,
              [2.60, 2.35, 2.35], font_size=8.0,
              alignments=[WD_ALIGN_PARAGRAPH.LEFT, WD_ALIGN_PARAGRAPH.CENTER, WD_ALIGN_PARAGRAPH.CENTER])
    add_note(doc, "AIC and BIC values are means across the 10 imputed data sets. The minimum-profile criterion required every profile to include at least 5% of the cohort. The maximum proportion range is the largest percentage-point range for an aligned profile across imputations. The two-profile solution was selected because it was the only candidate meeting the minimum-size, cross-imputation stability, and cross-domain interpretability criteria; selection was outcome-blind.")
    add_note(doc, "Mean profile proportions were 72.01% (range, 71.28%-72.58%) for Profile 1 and 27.99% (range, 27.42%-28.72%) for Profile 2. Lactate and creatinine were log1p transformed before within-imputation z-score standardization. Abbreviations: AIC, Akaike information criterion; BIC, Bayesian information criterion; pp, percentage points; WBC, white blood cell count.")
    path = OUT_TABLES / "Table2_lpa_selection_profile_definition.docx"
    doc.save(path)
    return path


def build_table3():
    survival = read_csv("30_profile_survival_results.csv")
    rmst = read_csv("31_profile_rmst.csv")
    primary = [r for r in survival if r["analysis_method"].startswith("Flexible parametric")]
    doc = setup_doc(landscape_page=True)
    add_title(doc, "Table 3. Association of inflammatory-metabolic stress profile with 28-day mortality in MIMIC-IV")
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(3)
    rr = p.add_run("Panel A. Flexible parametric survival models")
    rr.bold = True
    rr.font.name = "Arial"
    rr.font.size = Pt(9)
    rows_a = []
    for r in primary:
        label = r["model"] if "Sensitivity" not in r["model"] else "Sensitivity model"
        rows_a.append([
            label, r["adjustment"], f"{int(float(r['n'])):,}", f"{int(float(r['events'])):,}",
            f"{float(r['estimate']):.2f}",
            f"{float(r['confidence_interval_lower']):.2f}-{float(r['confidence_interval_upper']):.2f}",
            format_p(r["p_value"]),
        ])
    add_table(doc, ["Model", "Adjustment", "N", "Deaths", "HR", "95% CI", "P value"], rows_a,
              [1.20, 3.25, 0.72, 0.72, 0.72, 1.15, 1.28], font_size=7.8,
              alignments=[WD_ALIGN_PARAGRAPH.LEFT, WD_ALIGN_PARAGRAPH.LEFT] + [WD_ALIGN_PARAGRAPH.CENTER] * 5)

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(5)
    p.paragraph_format.space_after = Pt(3)
    rr = p.add_run("Panel B. Unadjusted restricted mean survival time through day 28")
    rr.bold = True
    rr.font.name = "Arial"
    rr.font.size = Pt(9)
    rows_b = []
    for r in rmst:
        label = r["estimand"].replace("RMST: ", "").replace("RMST difference: ", "Difference: ")
        rows_b.append([
            label, f"{int(float(r['n'])):,}", f"{int(float(r['events'])):,}",
            f"{float(r['estimate_days']):.2f}",
            f"{float(r['confidence_interval_lower_days']):.2f} to {float(r['confidence_interval_upper_days']):.2f}",
            "-" if not r["p_value"] else format_p(r["p_value"]),
        ])
    add_table(doc, ["Estimand", "N", "Deaths", "RMST, days", "95% CI, days", "P value"], rows_b,
              [3.15, 0.75, 0.75, 1.10, 1.40, 1.25], font_size=7.8,
              alignments=[WD_ALIGN_PARAGRAPH.LEFT] + [WD_ALIGN_PARAGRAPH.CENTER] * 5)
    add_note(doc, "Profile 1 (lower inflammatory-metabolic stress) is the reference; Profile 2 (higher inflammatory-metabolic stress) is the exposure. Flexible parametric models used a proportional-hazards parameterization with a 4-df baseline spline. The six profile-defining laboratory variables were not included as adjustment covariates. The SOFA and SAPS-II specifications were fitted separately.")
    add_note(doc, "Proportional-hazards tests were significant in Cox sensitivity models; hazard ratios should therefore be interpreted as time-averaged associations over 28 days. The RMST difference is Profile 2 minus Profile 1. Complete-case results are reported separately. Abbreviations: CI, confidence interval; HR, hazard ratio; RMST, restricted mean survival time; SAPS-II, Simplified Acute Physiology Score II; SOFA, Sequential Organ Failure Assessment.")
    path = OUT_TABLES / "Table3_mimic_mortality_rmst.docx"
    doc.save(path)
    return path


def build_table4():
    perf = read_csv("35_profile_transfer_model_performance.csv")
    assign = read_csv("36_eicu_profile_assignment_qc.csv")
    outcome = read_csv("38_eicu_profile_outcome_results.csv")
    selected_metrics = [
        "AUROC", "Brier score", "Classification agreement", "Cohen kappa",
        "Sensitivity for Profile 2", "Specificity for Profile 1",
        "Calibration intercept", "Calibration slope",
    ]
    by_metric = {r["metric"]: r for r in perf}
    rows_a = []
    for metric in selected_metrics:
        r = by_metric[metric]
        rows_a.append([metric, f"{float(r['estimate']):.3f}",
                       f"{float(r['lower_bound']):.3f}-{float(r['upper_bound']):.3f}", r["bound_type"]])
    rows_b = []
    for pnum in ("Profile 1", "Profile 2"):
        r = next(r for r in assign if r["section"] == "profile_distribution_summary"
                 and r["metric"] == "profile_percent_across_imputations" and r["profile"] == pnum)
        label = "Profile 1: lower stress" if pnum == "Profile 1" else "Profile 2: higher stress"
        rows_b.append([label, f"{float(r['percent']):.2f}", f"{float(r['minimum']):.2f}-{float(r['maximum']):.2f}"])
    unanimous = next(r for r in assign if r["metric"] == "unanimous_assignment_n")
    rows_b.append(["Unanimous assignment across 10 imputations", f"{float(unanimous['percent']):.2f}",
                   f"{int(float(unanimous['value'])):,}/{int(float(unanimous['denominator_n'])):,}"])

    primary = [r for r in outcome if r["analysis_set"] == "MICE primary"]
    rows_c = []
    for r in primary:
        rows_c.append([r["model"], r["adjustment"], f"{int(float(r['n'])):,}", f"{int(float(r['events'])):,}",
                       f"{float(r['estimate']):.2f}",
                       f"{float(r['confidence_interval_lower']):.2f}-{float(r['confidence_interval_upper']):.2f}",
                       format_p(r["p_value"])])
    doc = setup_doc(landscape_page=True)
    add_title(doc, "Table 4. MIMIC transfer-classifier performance and external validation in eICU")

    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(2)
    rr = p.add_run("Panel A. Internal performance of the frozen MIMIC transfer classifier"); rr.bold = True; rr.font.name = "Arial"; rr.font.size = Pt(8.8)
    add_table(doc, ["Metric", "Estimate", "Interval", "Interval type"], rows_a,
              [2.55, 1.05, 1.35, 4.25], font_size=7.2,
              alignments=[WD_ALIGN_PARAGRAPH.LEFT, WD_ALIGN_PARAGRAPH.CENTER, WD_ALIGN_PARAGRAPH.CENTER, WD_ALIGN_PARAGRAPH.LEFT])

    p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(3); p.paragraph_format.space_after = Pt(2)
    rr = p.add_run("Panel B. eICU profile assignment across 10 imputed data sets (N = 1,407)"); rr.bold = True; rr.font.name = "Arial"; rr.font.size = Pt(8.8)
    add_table(doc, ["Assignment metric", "Mean, %", "Range or count"], rows_b,
              [4.45, 1.35, 1.80], font_size=7.4,
              alignments=[WD_ALIGN_PARAGRAPH.LEFT, WD_ALIGN_PARAGRAPH.CENTER, WD_ALIGN_PARAGRAPH.CENTER])

    p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(3); p.paragraph_format.space_after = Pt(2)
    rr = p.add_run("Panel C. eICU hospital mortality association in the MICE primary analysis"); rr.bold = True; rr.font.name = "Arial"; rr.font.size = Pt(8.8)
    add_table(doc, ["Model", "Adjustment", "N", "Deaths", "OR", "95% CI", "P value"], rows_c,
              [1.05, 2.65, 0.68, 0.70, 0.65, 1.10, 1.20], font_size=7.4,
              alignments=[WD_ALIGN_PARAGRAPH.LEFT, WD_ALIGN_PARAGRAPH.LEFT] + [WD_ALIGN_PARAGRAPH.CENTER] * 5)
    add_note(doc, "The logistic transfer classifier was trained only in MIMIC using the six frozen profile variables and was applied in eICU without retraining, recalibration, or threshold modification (threshold = 0.50). Internal intervals in Panel A are 95% bootstrap intervals except for the calibration intercept and slope, which use 95% Wald confidence intervals.")
    add_note(doc, "Profile 1 is the reference and Profile 2 is the exposure in Panel C. Eleven eICU records with unknown hospital status were excluded (analysis N = 1,396). Hospital mortality is not harmonized with MIMIC 28-day all-cause mortality. Reliable post-discharge vital status and death time were unavailable; a 28-day all-cause mortality model was therefore not estimated. Complete-case results are reported separately. Abbreviations: AUROC, area under the receiver operating characteristic curve; CI, confidence interval; MICE, multiple imputation by chained equations; OR, odds ratio.")
    path = OUT_TABLES / "Table4_transfer_eicu_validation.docx"
    doc.save(path)
    return path


def build_legends():
    doc = setup_doc(landscape_page=False)
    add_title(doc, "Figure Legends")
    legends = [
        ("Figure 1. Study flow and outcome-blind analytic workflow in the MIMIC-IV discovery and eICU external validation cohorts.",
         "Panel A shows the MIMIC-IV cohort sequence from 94,458 adult ICU admissions to 41,295 Sepsis-3 ICU admissions, 6,691 cancer-plus-sepsis ICU admissions, 5,904 admissions with an ICU length of stay of at least 24 hours, and the final solid-tumor discovery cohort of 4,380 after excluding hematologic malignancy (n = 1,394) and mixed malignancy (n = 130). Multiple imputation was performed in 10 data sets, candidate latent profile models containing two to six profiles were compared without outcome information, and a two-profile solution was frozen (Profile 1, lower inflammatory-metabolic stress, n = 3,197; Profile 2, higher inflammatory-metabolic stress, n = 1,183). The 28-day mortality analysis included all 4,380 admissions; the complete-case sensitivity cohort included 2,649 admissions. Panel B shows the eICU sequence from 200,234 adult ICU admissions to 26,856 admissions meeting the diagnostic sepsis screening proxy, 1,736 solid-tumor-plus-sepsis candidates, and 1,407 admissions with an ICU length of stay of at least 24 hours. After eICU-specific multiple imputation, the frozen MIMIC logistic transfer classifier was applied without retraining or threshold modification. Hospital status was known for 1,396 admissions; 11 with unknown status were excluded. The complete-case outcome sensitivity cohort included 737 admissions. The eICU sepsis screen was a database-specific diagnostic proxy and was not fully harmonized to MIMIC Sepsis-3. Profile discovery and transfer-rule construction were outcome-blind."),
        ("Figure 2. Standardized inflammatory-metabolic patterns of the selected two-profile solution in MIMIC-IV.",
         "Lines show the mean standardized values of white blood cell count, platelet count, lactate, glucose, creatinine, and bicarbonate for the aligned profiles across 10 imputed data sets. Lactate and creatinine were log1p transformed before within-imputation z-score standardization. Profile 1 (lower inflammatory-metabolic stress; mean proportion, 72.01%; range, 71.28%-72.58%) was characterized by lower white blood cell count, lactate, glucose, and creatinine, higher bicarbonate, and platelet values near the cohort mean. Profile 2 (higher inflammatory-metabolic stress; mean proportion, 27.99%; range, 27.42%-28.72%) showed the reciprocal pattern. Labels describe clinically derived states and do not imply mechanistic biological subtypes. Profiles were selected and labeled without mortality outcomes."),
        ("Figure 3. Twenty-eight-day survival by inflammatory-metabolic stress profile in the MIMIC-IV discovery cohort.",
         "Kaplan-Meier curves compare Profile 1 (lower inflammatory-metabolic stress; n = 3,197) and Profile 2 (higher inflammatory-metabolic stress; n = 1,183) for 28-day mortality. Numbers at risk and confidence bands are shown as plotted. Profile 1 is the reference. Adjusted associations are reported in Table 3 using flexible parametric survival models; the six laboratory variables defining the profiles were not included as covariates. Proportional-hazards tests were significant in Cox sensitivity models, so reported hazard ratios represent time-averaged associations over 28 days. The observational comparison does not establish causality."),
        ("Figure 4. External validation of the frozen MIMIC inflammatory-metabolic stress profile classifier in eICU.",
         "The frozen 10-component MIMIC logistic transfer classifier was applied to the eICU solid-tumor-plus-sepsis cohort after eICU-specific multiple imputation, without retraining, recalibration, or threshold modification (classification threshold = 0.50). The figure displays the eICU profile distribution and odds ratios for hospital mortality from logistic regression models with Profile 1 as the reference and Profile 2 as the exposure. The primary MICE analysis included 1,396 admissions with known hospital status and adjusted sequentially for age, sex, and cancer phenotype; the complete-case sensitivity analysis included 737 admissions. Profile-defining laboratory variables were not adjustment covariates. The eICU endpoint is hospital mortality and differs from the MIMIC 28-day all-cause mortality endpoint. Post-discharge vital status and reliable death time were unavailable in eICU; no 28-day all-cause mortality model was estimated."),
    ]
    for idx, (heading, body) in enumerate(legends):
        p = doc.add_paragraph()
        p.paragraph_format.keep_with_next = True
        p.paragraph_format.space_before = Pt(7 if idx else 2)
        p.paragraph_format.space_after = Pt(3)
        r = p.add_run(heading)
        r.font.name = "Arial"
        r.font.size = Pt(10.5)
        r.bold = True
        p2 = doc.add_paragraph()
        p2.paragraph_format.space_after = Pt(4)
        p2.paragraph_format.line_spacing = 1.05
        r2 = p2.add_run(body)
        r2.font.name = "Arial"
        r2.font.size = Pt(9.5)
    path = OUT_MANUSCRIPT / "figure_legends.docx"
    doc.save(path)
    return path


def draw_box(c, x, y, w, h, text, *, fill, border, font_size=8.7, bold=False):
    c.setFillColor(fill)
    c.setStrokeColor(border)
    c.setLineWidth(1.0)
    c.roundRect(x, y, w, h, 5, fill=1, stroke=1)
    style = ParagraphStyle(
        "box", fontName="Helvetica-Bold" if bold else "Helvetica", fontSize=font_size,
        leading=font_size + 2, textColor=colors.HexColor("#17202A"), alignment=TA_CENTER,
        spaceBefore=0, spaceAfter=0,
    )
    para = Paragraph(text, style)
    pw, ph = para.wrap(w - 12, h - 8)
    para.drawOn(c, x + (w - pw) / 2, y + (h - ph) / 2)


def draw_arrow(c, x, y_top, y_bottom, color):
    c.setStrokeColor(color)
    c.setFillColor(color)
    c.setLineWidth(1.25)
    c.line(x, y_top, x, y_bottom + 5)
    c.line(x, y_bottom + 5, x - 3.5, y_bottom + 10)
    c.line(x, y_bottom + 5, x + 3.5, y_bottom + 10)


def build_flowchart():
    OUT_FIGURES.mkdir(parents=True, exist_ok=True)
    pdf_path = OUT_FIGURES / "Figure1_flowchart.pdf"
    page_w, page_h = landscape(letter)
    c = canvas.Canvas(str(pdf_path), pagesize=(page_w, page_h), pageCompression=1)
    c.setTitle("Figure 1. Study flow and outcome-blind analytic workflow")
    # Author metadata is intentionally omitted from generated artifacts;
    # code authorship is recorded in CITATION.cff.
    # until the authors and copyright holder approve the public release.
    c.setFillColor(colors.HexColor("#17202A"))
    c.setFont("Helvetica-Bold", 15)
    c.drawCentredString(page_w / 2, page_h - 25, "Study flow and outcome-blind analytic workflow")

    left_x, right_x = 30, 411
    panel_w = 351
    panel_y, panel_h = 63, 505
    panel_border = colors.HexColor("#5D6D7E")
    for x, title, subtitle, accent in [
        (left_x, "A  MIMIC-IV discovery cohort", "Profile discovery and 28-day outcome analysis", colors.HexColor("#2471A3")),
        (right_x, "B  eICU external validation cohort", "Frozen classifier transport and hospital outcome analysis", colors.HexColor("#148F77")),
    ]:
        c.setFillColor(colors.white)
        c.setStrokeColor(panel_border)
        c.setLineWidth(1)
        c.roundRect(x, panel_y, panel_w, panel_h, 7, fill=1, stroke=1)
        c.setFillColor(accent)
        c.roundRect(x, panel_y + panel_h - 48, panel_w, 48, 7, fill=1, stroke=0)
        c.setFillColor(colors.white)
        c.setFont("Helvetica-Bold", 11.5)
        c.drawString(x + 12, panel_y + panel_h - 20, title)
        c.setFont("Helvetica", 8.4)
        c.drawString(x + 12, panel_y + panel_h - 35, subtitle)

    mimic_boxes = [
        "Adult ICU admissions<br/><b>n = 94,458</b>",
        "Sepsis-3 ICU admissions<br/><b>n = 41,295</b>",
        "Cancer + sepsis ICU admissions<br/><b>n = 6,691</b>",
        "ICU length of stay >=24 h<br/><b>n = 5,904</b>",
        "Solid-tumor discovery cohort<br/><b>n = 4,380</b>",
        "MICE (m = 10); LPA models with 2-6 profiles<br/><b>Outcome-blind selection</b>",
        "Frozen two-profile solution<br/><b>P1 n = 3,197 | P2 n = 1,183</b>",
        "28-day mortality analysis: <b>n = 4,380</b><br/>Complete-case sensitivity: <b>n = 2,649</b>",
    ]
    eicu_boxes = [
        "Adult ICU admissions<br/><b>n = 200,234</b>",
        "Diagnostic sepsis screening proxy<br/><b>n = 26,856</b>",
        "Solid tumor + sepsis candidates<br/><b>n = 1,736</b>",
        "ICU length of stay >=24 h<br/><b>n = 1,407</b>",
        "eICU-specific MICE (m = 10)<br/><b>No mortality in imputation</b>",
        "Frozen MIMIC logistic classifier<br/><b>No retraining or threshold change</b>",
        "Hospital status known: <b>n = 1,396</b><br/>Unknown status excluded: <b>n = 11</b>",
        "Hospital mortality analysis: <b>n = 1,396</b><br/>Complete-case sensitivity: <b>n = 737</b>",
    ]
    box_w, box_h, gap = 293, 42, 10
    start_y = panel_y + panel_h - 63 - box_h
    fill_left = colors.HexColor("#EBF5FB")
    fill_right = colors.HexColor("#E8F8F5")
    for x, boxes, fill, border in [
        (left_x + 29, mimic_boxes, fill_left, colors.HexColor("#2471A3")),
        (right_x + 29, eicu_boxes, fill_right, colors.HexColor("#148F77")),
    ]:
        y = start_y
        for idx, text in enumerate(boxes):
            draw_box(c, x, y, box_w, box_h, text, fill=fill, border=border,
                     font_size=8.1 if idx >= 5 else 8.5)
            if idx < len(boxes) - 1:
                draw_arrow(c, x + box_w / 2, y, y - gap, border)
            y -= box_h + gap

    # Prespecified exclusions and harmonization caveats are printed as linked notes, not inferred counts.
    c.setFillColor(colors.HexColor("#566573"))
    c.setFont("Helvetica", 7.4)
    c.drawString(left_x + 8, 48, "Solid-tumor restriction excluded hematologic malignancy (n = 1,394) and mixed malignancy (n = 130).")
    c.drawString(right_x + 8, 48, "The eICU diagnostic sepsis proxy was database-specific and not fully harmonized to MIMIC Sepsis-3.")
    c.setFont("Helvetica-Oblique", 7.3)
    c.drawCentredString(page_w / 2, 28, "Profile discovery, labeling, and transfer-rule construction were performed without mortality outcomes.")
    c.showPage()
    c.save()

    tiff_base = OUT_FIGURES / "Figure1_flowchart"
    pdftoppm = os.environ.get("PDFTOPPM_BIN", "pdftoppm")
    subprocess.run([pdftoppm, "-f", "1", "-singlefile", "-r", "300", "-tiff", str(pdf_path), str(tiff_base)],
                   check=True)
    tiff_path = OUT_FIGURES / "Figure1_flowchart.tif"
    requested_path = OUT_FIGURES / "Figure1_flowchart.tiff"
    if tiff_path.exists():
        tiff_path.replace(requested_path)
    with Image.open(requested_path) as im:
        if im.mode != "RGB" or im.info.get("dpi") != (300.0, 300.0):
            rgb = im.convert("RGB")
            rgb.save(requested_path, format="TIFF", compression="tiff_lzw", dpi=(300, 300))
    return pdf_path, requested_path


def verify_no_shading(docx_paths: list[Path]):
    import zipfile
    for path in docx_paths:
        with zipfile.ZipFile(path, "r") as zf:
            xml = zf.read("word/document.xml").decode("utf-8")
        if "<w:shd" in xml:
            raise RuntimeError(f"Cell shading found in {path}")


def main():
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUT_FIGURES.mkdir(parents=True, exist_ok=True)
    OUT_MANUSCRIPT.mkdir(parents=True, exist_ok=True)
    flow = build_flowchart()
    docx_paths = [build_legends(), build_table1(), build_table2(), build_table3(), build_table4()]
    verify_no_shading(docx_paths)
    for path in [*flow, *docx_paths]:
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
