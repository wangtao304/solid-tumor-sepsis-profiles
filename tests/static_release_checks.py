#!/usr/bin/env python3
"""Static, non-clinical checks for the v1.0.0 public source release."""

from __future__ import annotations

import ast
import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

FORBIDDEN_SUFFIXES = {
    ".rds", ".rdata", ".rda", ".parquet", ".feather", ".sas7bdat",
    ".dta", ".sav", ".xlsx", ".xls", ".doc", ".docx", ".pdf",
    ".tif", ".tiff", ".png", ".jpg", ".jpeg", ".ttf", ".otf",
    ".woff", ".woff2", ".log",
}

FORBIDDEN_PARTS = {".git", "__pycache__", ".pytest_cache"}

REQUIRED_FILES = {
    "README.md",
    "THIRD_PARTY_NOTICES.md",
    "LICENSE",
    "CITATION.cff",
    "config/config.example.yml",
    "docs/analysis_workflow.md",
    "docs/code_to_results_map.csv",
    "docs/source_release_crosswalk.csv",
    "environment/DEPENDENCIES.md",
    "environment/session_info.txt",
    "tests/run_lpa_em_synthetic_self_test.R",
    "tests/test_probability_ensemble_synthetic.R",
}

TEXT_SUFFIXES = {".r", ".sql", ".py", ".md", ".txt", ".csv", ".yml", ".yaml", ".cff", ""}

PATTERNS = {
    "windows_absolute_path": re.compile(r"(?i)\b[A-Z]:[\\/](?!TO-BE-CONFIRMED)"),
    "credential_uri": re.compile(r"(?i)(postgres(?:ql)?|mysql|mssql)://[^\s/:]+:[^\s/@]+@"),
    "private_key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "aws_access_key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
    "github_token": re.compile(r"\bgh[pousr]_[A-Za-z0-9]{20,}\b"),
    "generic_secret_assignment": re.compile(
        r"(?im)^\s*(?:password|passwd|pwd|token|api[_-]?key|secret)\s*[:=]\s*[\"']?(?!\$\{|TO.BE.CONFIRMED|<)[^\s\"']{8,}"
    ),
}

EXPECTED_SQL = {
    "R/11_run_multiple_imputation.R": "sql/11_mice_input.sql",
    "R/24_run_inflammatory_metabolic_imputation.R": "sql/24_im_new_mice_input.sql",
    "R/30_analyze_profile_survival.R": "sql/30_profile_survival_input.sql",
    "R/33_execute_stage9b.R": "sql/30_profile_survival_input.sql",
    "R/36_prepare_eicu_profile_assignment.R": "sql/36_eicu_profile_assignment_input.sql",
    "R/38_validate_eicu_profile_outcomes.R": "sql/38_eicu_profile_outcome_input.sql",
    "R/41_first_eligible_stay_sensitivity.R": "sql/41_first_eligible_stay_sensitivity_input.sql",
    "R/revision/R1B_cluster_robust_analysis.R": "sql/41_first_eligible_stay_sensitivity_input.sql",
    "R/revision/R1B_time_varying_analysis_ready.R": "sql/revision/R1B_landmark_boundary_input.sql",
    "R/revision/R1C_time_varying_profile_analysis.R": "sql/revision/R1B_landmark_boundary_input.sql",
}


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def main() -> int:
    errors: list[str] = []
    files = [p for p in ROOT.rglob("*") if p.is_file()]

    for path in files:
        rp = rel(path)
        if any(part.lower() in FORBIDDEN_PARTS for part in path.parts):
            errors.append(f"forbidden path component: {rp}")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"forbidden file type: {rp}")

    for required in sorted(REQUIRED_FILES):
        if not (ROOT / required).is_file():
            errors.append(f"missing required file: {required}")

    for path in files:
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8-sig")
        except UnicodeDecodeError:
            errors.append(f"non-UTF-8 text candidate: {rel(path)}")
            continue
        for name, pattern in PATTERNS.items():
            if pattern.search(text):
                errors.append(f"{name}: {rel(path)}")
        if path.suffix.lower() == ".py":
            try:
                ast.parse(text, filename=rel(path))
            except SyntaxError as exc:
                errors.append(f"python syntax: {rel(path)}:{exc.lineno}: {exc.msg}")

    for source, sql in EXPECTED_SQL.items():
        source_path = ROOT / source
        sql_path = ROOT / sql
        if not source_path.is_file():
            errors.append(f"missing source dependency owner: {source}")
        if not sql_path.is_file():
            errors.append(f"missing SQL dependency: {sql}")

    for csv_path in [ROOT / "docs/code_to_results_map.csv", ROOT / "docs/source_release_crosswalk.csv"]:
        if csv_path.is_file():
            try:
                with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
                    rows = list(csv.reader(handle))
                if len(rows) < 2 or len(rows[0]) < 3:
                    errors.append(f"incomplete CSV: {rel(csv_path)}")
            except csv.Error as exc:
                errors.append(f"CSV parse: {rel(csv_path)}: {exc}")

    if errors:
        print("STATIC RELEASE CHECK: FAIL")
        for error in errors:
            print(f"- {error}")
        return 1

    print(f"STATIC RELEASE CHECK: PASS ({len(files)} files inspected)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

