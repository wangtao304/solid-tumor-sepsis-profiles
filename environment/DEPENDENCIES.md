# Environment and dependencies

## Recorded analytical environment

| Component | Recorded version | Role | Evidence status |
|---|---:|---|---|
| R | 4.6.0, Windows UCRT x86-64, locale C | Statistical runtime | Frozen in R1A environment audit |
| mice | 3.19.0 | Multiple imputation | Frozen in R1A environment audit |
| survival | 3.8-6 | Survival objects and Kaplan-Meier calculations | Package DESCRIPTION; session summary reports 3.8.6 |
| rstpm2 | 1.7.1 | Flexible parametric survival models | Frozen; R1B recorded binary built under R 4.6.1 |
| PostgreSQL client/server | Not frozen | Controlled SQL streaming | Version limitation |
| MIMIC-IV | v3.0 | Discovery database | Study definition |
| eICU-CRD | v2.0 | Transport-evaluation database | Study definition |

The exact transitive R package versions were not frozen and no `renv.lock` exists. The release does not claim an exact binary environment reconstruction.

## Python figure and table tooling

`python/stage12a_build_figures_tables.py` imports Pillow, python-docx, and ReportLab and may call Poppler `pdftoppm`. Static inspection found no `pypdf` import or call, so `pypdf` is not a runtime dependency of this release. The exact versions used in the original Stage 12A run were not preserved. No package binary, Poppler executable, or font file is included.

The versions visible on the release-audit host are not evidence of the historical Stage 12A environment and are therefore not presented as the original versions. There is no `requirements.txt` or Python lockfile that can support an exact historical reconstruction.

## Database and upstream-schema dependencies

The R scripts invoke an externally supplied PostgreSQL `psql` client. The SQL assumes the MIMIC-IV and eICU schemas configured by the authorized user. It also refers to `mimiciv_derived` relations produced outside this package. The exact MIMIC Code commit used to build those relations was not recorded.

## Operating-system-specific behavior

The recorded analyses ran on Windows. The R1C source uses Windows `certutil` for historical SHA-256 gates. Users on another operating system must implement an equivalent integrity check only after documenting and approving the change; a failed gate must not be ignored.

The Stage 12A Python source reads `PDFTOPPM_BIN` or searches `PATH`. Word document rendering and final manuscript assembly were performed outside this package and are not a release-package dependency claim.

## Installation boundary

This repository does not install or upgrade software. Users must provision compatible runtimes and packages separately. See `THIRD_PARTY_NOTICES.md` for external license information and `LICENSE` for the original project-code license.

