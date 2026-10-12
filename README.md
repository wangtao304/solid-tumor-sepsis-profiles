# MIMIC cancer sepsis inflammatory-metabolic stress profiles

## Release status

Version 1.0.1 is the current audited source-code release supporting the revised manuscript submitted to *Frontiers in Cellular and Infection Microbiology* (manuscript ID 1999181). The repository is `https://github.com/wangtao304/solid-tumor-sepsis-profiles`, and the fixed release is tagged `v1.0.1`. The previously released v1.0.0 remains preserved.

This documentation-only release corrects the stated clinical dataset version to MIMIC-IV v3.1, as confirmed by the authors. No analytical source code was changed. The previously released v1.0.0 is preserved. The historical build versions of the 11 upstream `mimiciv_derived` relations remain unresolved.

The release contains the reviewed SQL, R, and Python source used to support the revised manuscript and Supplementary Material. It does not contain clinical records, patient-level exports, completed imputations, profile assignments, fitted model objects, database credentials, or frozen manuscript results.

## Study scope

The code supports:

- MIMIC-IV discovery-cohort construction and aggregate quality control;
- six-variable missing-data preparation and multiple imputation;
- outcome-blind latent profile analysis and model selection;
- MIMIC clinical description and mortality analyses;
- MIMIC-only profile-transfer classifier construction;
- eICU-CRD transport evaluation;
- repeated-stay, first-eligible-stay, patient-clustered, landmark, time-varying, and death-date precision sensitivity analyses;
- aggregate table and figure generation.

The six profile-defining variables are white blood cell count, platelet count, lactate, glucose, creatinine, and bicarbonate. No script in this package may use mortality to select or label a profile.

## Data access and protected inputs

The study used MIMIC-IV v3.1 and eICU Collaborative Research Database v2.0. Neither database is redistributed. Users must obtain their own credentialed access, comply with the applicable PhysioNet and institutional terms, and run database code only inside an approved environment.

The workflow also depends on locally materialized `mimiciv_derived` relations. A bounded local provenance review did not recover the historical MIMIC Code checkout, immutable commit or tag, materialization date, or build manifest. The historical build version for all 11 consumed relations is therefore `HISTORICAL_VERSION_UNRESOLVED`: `bg`, `blood_differential`, `charlson`, `chemistry`, `complete_blood_count`, `first_day_sofa`, `icustay_detail`, `sapsii`, `sepsis3`, `vasoactive_agent`, and `ventilation`.

MIMIC-IV v3.1 is the clinical database release and is not evidence of any MIMIC Code tag or commit, including MIMIC Code v3.0.0. No current upstream `HEAD` is presented as the historical build. The upstream definition SQL is not redistributed, historical local modifications cannot be assessed, and byte-for-byte reconstruction of the derived schema is not claimed. See `docs/upstream_schema_provenance.md`.

Several R scripts create or read protected `.rds` objects under `data/processed/` or `data/interim/`. Those objects can contain completed patient-level data, profile assignments, or fitted parameters and are deliberately excluded. Authorized users can regenerate them locally by following `docs/analysis_workflow.md`. They must not commit, upload, or share them.

The SQL used for MICE, outcome linkage, repeated-stay selection, and landmark analyses may stream identifiers or timestamps into the corresponding R process. This is permitted only inside the authorized environment. The release scripts write aggregate outputs only; users must not redirect those streams to patient-level files.

## Software environment

The recorded R analysis environment was:

- R 4.6.0, Windows UCRT x86-64, locale `C`;
- `mice` 3.19.0;
- `survival` 3.8-6 (reported as 3.8.6 in the session summary);
- `rstpm2` 1.7.1.

The R1B audit recorded R 4.6.0 with an `rstpm2` binary built under R 4.6.1. The exact PostgreSQL client/server versions and the original Python/Poppler versions used for Stage 12A were not frozen. There is no `renv.lock`, Python lockfile, or complete transitive-dependency snapshot. See `environment/DEPENDENCIES.md` and `environment/session_info.txt`.

## Configuration

Run all scripts from the package root. The scripts read environment variables directly; `config/config.example.yml` is documentation, not an automatically parsed configuration file. Database variables have empty defaults and scripts fail closed when required values are absent.

Required variables depend on the stage and include:

- `MIMIC_PSQL_BIN`, `MIMIC_PGPASS_FILE`, `MIMIC_DB_HOST`, `MIMIC_DB_PORT`, `MIMIC_DB_NAME`, `MIMIC_DB_USER`;
- `EICU_DB_HOST`, `EICU_DB_PORT`, `EICU_DB_NAME`, `EICU_DB_USER`;
- `PDFTOPPM_BIN` for optional PDF rendering when `pdftoppm` is not on `PATH`.

Reviewer-revision scripts have explicit execution gates. Their names and safety purpose are documented in `config/config.example.yml` and `docs/analysis_workflow.md`.

Do not store a password, token, connection URI, or real local path in tracked files. Use process environment variables and an approved password file outside the repository.

## Run order

The authoritative stage-by-stage order is in `docs/analysis_workflow.md`. The high-level sequence is:

1. Construct and audit the database-resident MIMIC cohort.
2. Create the final solid-tumor six-variable MICE input and protected MICE object.
3. Fit 2-6 profiles and apply the outcome-blind selection rules.
4. Produce aggregate profile descriptions and MIMIC outcome analyses.
5. Build the MIMIC-only 10-component logistic transfer ensemble.
6. Prepare eICU completed datasets, apply the frozen ensemble, and perform the transport evaluation.
7. Run repeated-stay and first-stay checks.
8. Run the executed R1B patient-clustered analysis and R1C landmark/time-varying sensitivity workflow when explicitly authorized.
9. Render aggregate figures and tables.

`docs/code_to_results_map.csv` maps each source file to manuscript tables, figures, supplementary items, and protected intermediates.

### Stage-to-output summary

| Code stage | Principal manuscript or supplementary outputs |
|---|---|
| Cohort construction and missingness | Figure 1 and Supplementary Table S8 |
| Final MIMIC MICE | Supplementary Table S2 |
| LPA fitting and outcome-blind selection | Table 2, Figure 2, and Supplementary Table S1 |
| Profile description | Table 1 |
| MIMIC outcomes and clustered inference | Table 3, Figure 3, and Supplementary Tables S3, S4, and S9 |
| Profile transfer and eICU transport evaluation | Table 4, Figure 4, and Supplementary Tables S5 and S6 |
| First-eligible-stay sensitivity | Supplementary Table S7 |
| Landmark, time-varying, and date-precision analyses | Supplementary Table S10 and Supplementary Figures S1 and S2 |

The complete source-to-result mapping and stage-level limitations are in `docs/code_to_results_map.csv`.

The internal release manifest uses three provenance categories:

- `IDENTICAL`: the release file is byte-identical to the declared reviewed source;
- `SANITIZED_OR_RELOCATED`: the release file is actual analysis source with only the documented configuration, path, integrity-fingerprint, or metadata changes; it is not byte-identical to the source;
- `GENERATED_FOR_RELEASE`: documentation, configuration examples, or static safety tests created for this source package rather than executed study-analysis code.

## Important implementation distinctions

- `R/33_execute_stage9b.R` is the authoritative executed MIMIC outcome implementation. `R/30_analyze_profile_survival.R` is retained as a supporting historical implementation referenced by the final Supplementary Methods.
- `R/revision/R1B_cluster_robust_analysis.R` supplies the patient-clustered confidence intervals and P values used in the R1F manuscript. The point estimates and RMST remain linked to the frozen primary analysis.
- `R/revision/R1B_time_varying_analysis_ready.R` performs the 24-hour landmark population and death-date boundary preflight. It does not fit the R1C model.
- `R/revision/R1C_time_varying_profile_analysis.R` fits the frozen landmark time-varying and date-precision sensitivity models. It requires the R1B aggregate preflight and clustered results.
- `python/stage12a_build_figures_tables.py` is the actual Stage 12A assembly source for Figure 1 and the original table/figure layout. It predates the R1B clustered-interval revision and must not be used alone to recreate the final R1F Table 3. The final table requires the R1B clustered results.
- The R1C script retains integrity checks for the protected frozen objects and executed R1B aggregates. A legitimately rebuilt object may not be byte-identical to the historical object. Do not bypass a failed hash gate silently; document and authorize any provenance reconciliation.

## Reproducibility limits

This package does not claim one-command complete reproduction. The main limits are:

- credentialed clinical databases and a secure execution environment are required;
- the protected fitted/imputed objects are regenerated locally and are not supplied;
- the exact upstream `mimiciv_derived` commit was not recorded;
- the original Python, PostgreSQL, Poppler, and transitive package versions were not completely frozen;
- the final Word manuscript assembly is outside this source package;
- aggregate result CSVs are not bundled, so figure/table renderers run only after authorized upstream stages create them;
- the R1C integrity gates reflect the historical executed files and require explicit provenance review for rebuilt objects.

## Validation performed on this release

The release review performs only tests that do not access clinical data:

- parsing every included R source file;
- parsing every included Python source file;
- running the custom EM synthetic self-test;
- running the transfer-probability ensemble synthetic test;
- checking declared local source and SQL dependencies;
- scanning file names and file contents for protected data, fitted objects, identifiers with literal values, credentials, private paths, logs, manuscripts, reviews, fonts, and Git history;
- checking the final archive contents against the public-release manifest.

Clinical SQL, MICE, LPA, transfer training, and outcome models are not executed during release validation. The executed non-clinical test results are recorded in `docs/validation_results.md`. This release has not undergone a new real-data end-to-end reproduction.

## Contents

- `R/`: final analysis sources and executed reviewer-revision sources.
- `sql/`: cohort, variable, linkage, audit, and landmark SQL.
- `python/`: aggregate figure/table assembly source.
- `tests/`: wholly synthetic, non-clinical tests.
- `docs/`: workflow, method, provenance, and code-to-result mapping.
- `environment/`: recorded environment evidence and dependency limitations.
- `config/`: placeholder-only configuration documentation.
- `data/` and `outputs/`: intentionally empty protected/local working locations.

## Authorship, ownership, and license

Copyright (c) 2026 Wang Tao and Li Yanan. Code contributors: Wang Tao and Li Yanan. Original project-authored code is released under the MIT License; see `LICENSE`.

The MIT License does not cover MIMIC-IV or eICU-CRD data, externally maintained MIMIC Code definitions, third-party software or source outside the authors' ownership, or fonts. External components remain governed by their own terms and licenses; see `THIRD_PARTY_NOTICES.md`.

`CITATION.cff` identifies the software contributors only. It does not alter or replace the manuscript author list.

## Public availability

The source is published at `https://github.com/wangtao304/solid-tumor-sepsis-profiles`. Version 1.0.1 is the current fixed release supporting the revised manuscript; v1.0.0 remains preserved as the original release. The release does not include patient-level data, completed imputations, trained or fitted model objects, or protected database-derived artifacts.

