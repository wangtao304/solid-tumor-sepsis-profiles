# Third-Party Notices

This source release calls external software but does not redistribute third-party package source code, package binaries, fonts, database contents, or clinical data. No third-party source code, binary, or font file was identified as vendored in the audited release. Any later distribution that adds binaries or vendored source requires a new license review.

## R runtime and packages

| Component | Recorded version | Role | License reported by package metadata | Upstream |
|---|---:|---|---|---|
| R | 4.6.0 | Statistical runtime and base/recommended packages | R licensing; GPL version 2 or 3 framework | https://www.r-project.org/ |
| mice | 3.19.0 | Multiple imputation | GPL (>= 2) | https://github.com/amices/mice |
| survival | 3.8-6 | Survival objects and Kaplan-Meier calculations | LGPL (>= 2) | https://github.com/therneau/survival |
| rstpm2 | 1.7.1 | Flexible parametric survival models | GPL-2 or GPL-3 | https://github.com/mclements/rstpm2 |

The project invokes these packages through their public R interfaces. Their source and binaries are not included. Their licenses do not constitute a license grant for this project's own source code.

## Python figure and table tooling

The release scripts import Pillow, python-docx, and ReportLab and call the external Poppler `pdftoppm` utility. The exact versions used for the original Stage 12A run were not preserved. No Python package, Poppler binary, or font file is included in this release.

License types observed during the release audit were MIT-CMU (Pillow), MIT (python-docx), and BSD-style (ReportLab). These observations must not be represented as a historical environment lock. The exact Poppler release and license notice must be confirmed if Poppler is redistributed; this package only calls an externally supplied executable.

An earlier internal dependency narrative listed `pypdf`. Static inspection of the released Python source found no `pypdf` import or call, so `pypdf` is not treated as a runtime dependency of this release.

## PostgreSQL

The analysis scripts use an externally supplied PostgreSQL command-line client. PostgreSQL is copyright the PostgreSQL Global Development Group, with portions copyright the Regents of the University of California, and is distributed under the PostgreSQL License. No PostgreSQL binary is included. The exact client and server versions used by the original analysis were not preserved.

## Clinical databases and upstream schemas

MIMIC-IV v3.1 and eICU-CRD v2.0 data are not included. Access requires the user's own authorization and compliance with applicable PhysioNet, database, and institutional terms.

The SQL refers to externally materialized `mimiciv_derived` relations. Their upstream SQL definitions are not bundled. A bounded local review did not recover the historical MIMIC Code source snapshot or build record, so all 11 consumed relations are documented as `HISTORICAL_VERSION_UNRESOLVED`. MIMIC-IV v3.1 is not treated as a MIMIC Code version, and no current upstream revision is represented as the historical build. If future material adds or adapts upstream definitions, the exact upstream copyright and license notices must be retained.

## Project license boundary

The MIT License in `LICENSE` applies only to original project code and documentation for which Wang Tao and Li Yanan hold the rights to grant permission. It does not relicense R, Python packages, PostgreSQL, Poppler, MIMIC Code, MIMIC-IV, eICU-CRD, other external material, or fonts. Those components remain governed by their respective licenses and terms.

