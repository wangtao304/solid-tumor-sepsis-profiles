# Analysis workflow

## Scope and execution boundary

This document maps the actual source chain supporting the R1F manuscript. It is not an instruction to run clinical analyses outside an approved environment. Database access, patient-level in-memory linkage, and protected object generation require the user's own authorization and credentials. The v1.0.1 public source release contains no input data or fitted object.

Run scripts from the package root. Do not write patient-level SQL output to disk. Aggregate exports belong under `outputs/`; protected objects belong under ignored `data/processed/` or `data/interim/` paths.

## A. Cohort and early feasibility history

The following SQL files preserve the implemented cohort and variable-audit history:

1. `sql/01_cohort_feasibility.sql`
2. `sql/02_variable_coverage.sql`
3. `sql/03_cancer_phenotype_distribution.sql`
4. `sql/04_eicu_validation_feasibility.sql`
5. `sql/05_variable_detailed_coverage.sql`
6. `sql/06_build_immune_metabolic_lca_base.sql`
7. `sql/10_derive_imputation_auxiliary_variables.sql`
8. `sql/10_imputation_auxiliary_variable_qc.sql`

`sql/06_build_immune_metabolic_lca_base.sql` creates the protected all-cancer database-resident base. It depends on externally materialized `mimiciv_derived` relations. The source-referenced relation list and unresolved upstream commit are recorded in `docs/upstream_schema_provenance.md`. `sql/10_derive_imputation_auxiliary_variables.sql` adds the approved first-24-hour auxiliary variables.

The historical all-cancer MICE and solid-tumor restriction branch is retained because it supports the documented cohort-development timeline and Figure 1:

1. `sql/11_mice_input.sql`
2. `R/11_run_multiple_imputation.R`
3. `R/19_solid_tumor_feasibility_audit.R`

`R/19_solid_tumor_feasibility_audit.R` also reads aggregate outputs 03 and 04. The historical all-cancer MICE object is protected and excluded.

## B. Final six-variable MIMIC preparation

1. `sql/22_im_new_variable_missingness.sql` creates the aggregate pre-imputation variable report.
2. `sql/23_im_complete_case_summary.sql` creates the aggregate six-variable complete-case report.
3. `sql/24_im_new_mice_input.sql` streams the approved analysis variables and auxiliary predictors to R inside the authorized environment.
4. `R/24_run_inflammatory_metabolic_imputation.R` performs the frozen MICE procedure and writes the protected MICE object plus aggregate QC table 24.

MIMIC settings are `m=10`, `maxit=20`, predictive mean matching, five donors, and seed 20260721. Mortality outcomes, identifiers, length-of-stay variables, differential counts, NLR, and derived ANC/ALC are excluded from the imputation model as implemented.

Protected output:

`data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds`

This object is deliberately absent from the package and must be rebuilt by an authorized user.

## C. LPA and outcome-blind model selection

1. `R/25_run_inflammatory_metabolic_lpa.R` reads the protected MICE object, applies the frozen transformations and z standardization, fits 2-6 diagonal Gaussian profiles separately in the 10 completed datasets, and writes protected fit records plus aggregate tables 25-26.
2. `R/27_select_inflammatory_metabolic_lpa.R` aligns labels, applies the frozen stability/admissibility rules, selects the candidate without outcomes, and writes tables 27-28 plus the profile-pattern figure.
3. `R/29_describe_lpa_profiles.R` combines the protected objects and aggregate selected-profile characteristics to write aggregate baseline Table 29.

The custom expectation-maximization implementation and synthetic test are documented in `docs/CUSTOM_EM_METHOD.md` and `tests/run_lpa_em_synthetic_self_test.R`.

Protected output:

`data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds`

The object contains fitted records and assignments and is excluded.

## D. MIMIC outcome analysis

The authoritative executed path is:

1. `sql/30_profile_survival_input.sql`
2. `R/33_execute_stage9b.R`

The SQL streams the required outcome and covariate linkage to memory. The R script reads the frozen MICE and LPA objects, validates one row per stay, fits the prespecified flexible parametric models, and writes aggregate tables 30-33 plus the Kaplan-Meier figure. `R/30_analyze_profile_survival.R` is retained as a supporting historical implementation because the final Supplementary Methods trace to it, but it is not the authoritative executor for the frozen Stage 9B results.

The final R1F patient-clustered inference is added by the R1B path in section H. Do not recreate the R1F Table 3 from the historical ordinary-variance output alone.

## E. Transfer classifier

`R/35_build_profile_transfer_model.R` reads the frozen MIMIC MICE and LPA objects. It fits one logistic component in each of the 10 completed MIMIC datasets using the frozen consensus profile label, stores component-specific transformations and scaling constants, and evaluates label fidelity. It averages component predicted probabilities, not coefficients.

Protected output:

`data/processed/inflammatory_metabolic_profile_transfer_logistic_mi10_seed20260723.rds`

The transfer object is excluded. It is regenerated locally by the source code; public distribution of the object itself requires a separate author and data-use decision.

`tests/test_probability_ensemble_synthetic.R` verifies probability averaging and the fixed threshold with wholly synthetic values.

## F. eICU preparation and transport evaluation

1. `sql/34_eicu_validation_feasibility.sql` writes aggregate feasibility information.
2. `sql/36_eicu_external_validation_preflight.sql` audits cohort/outcome feasibility without fitting a model.
3. `sql/36_eicu_profile_assignment_input.sql` streams the approved eICU variables to R.
4. `R/36_prepare_eicu_profile_assignment.R` performs eICU-specific MICE and applies every frozen MIMIC component to every completed eICU dataset. It averages the 10 component probabilities and applies the unchanged threshold of 0.5. It writes aggregate tables 36-37 and a protected local object.
5. `sql/38_eicu_profile_outcome_input.sql` streams hospital outcome linkage.
6. `R/38_validate_eicu_profile_outcomes.R` performs the frozen hospital-mortality logistic analyses and complete-case sensitivity, writing aggregate tables 38-39.
7. `R/38_render_eicu_external_validation.R` renders the external-evaluation figure from aggregate outputs.

eICU MICE settings are `m=10`, `maxit=20`, predictive mean matching, five donors, and seed 20260724. The eICU cohort uses a database-specific diagnostic sepsis screening proxy and is not a full replication of the MIMIC-IV Sepsis-3 operationalization. Hospital mortality is not equivalent to 28-day all-cause mortality.

Protected output:

`data/processed/eicu_inflammatory_metabolic_profile_transfer_mi10_seed20260724.rds`

It is excluded.

## G. Repeated-stay and first-stay sensitivity

1. `sql/40_repeated_icu_stay_cohort_gate.sql` reconciles the frozen discovery cohort.
2. `sql/40_repeated_icu_stay_audit.sql` writes aggregate repeated-stay counts.
3. `sql/41_first_eligible_stay_qc.sql` performs aggregate first-stay QC.
4. `sql/41_first_eligible_stay_sensitivity_input.sql` streams only the first eligible ICU stay per subject to R.
5. `R/41_first_eligible_stay_sensitivity.R` applies the frozen profile assignment and outcome models and writes aggregate Table 41.

The source SQL may return identifiers to the R process to enforce one stay per subject. Those rows must remain in memory inside the secure environment and must not be exported.

## H. R1B patient-clustered inference

`R/revision/R1B_cluster_robust_analysis.R` is the executed R1B source. It uses the original 4,380-stay cohort, unchanged profile assignments, unchanged point-estimate model, and `subject_id`-clustered sandwich covariance. It changes standard errors, confidence intervals, and P values only. No frailty is introduced.

Prerequisites:

- the protected MIMIC MICE and LPA objects;
- `sql/41_first_eligible_stay_sensitivity_input.sql` for the full linkage fields;
- the frozen `outputs/tables/30_profile_survival_results.csv`.

The script is gated by `R1B_EXECUTE_CLUSTER_ROBUST=YES`. It writes aggregate R1B results and diagnostics under `outputs/revision/`. These aggregate files are not bundled.

## I. R1B landmark preflight

`sql/revision/R1B_landmark_boundary_input.sql` provides the fields used to verify the 24-hour landmark population and death-date boundaries. `R/revision/R1B_time_varying_analysis_ready.R` performs the preflight and writes:

- `outputs/revision/R1B_landmark_population_audit.csv`;
- `outputs/revision/R1B_death_date_boundary_audit.csv`.

It does not fit the time-varying model. The script is gated by `R1B_EXECUTE_LANDMARK_AUDIT=YES`; the separate historical placeholder gate `R1B_EXECUTE_TIME_VARYING_MODEL` remains closed in that source.

## J. R1C landmark, time-varying, and date-precision analyses

`R/revision/R1C_time_varying_profile_analysis.R` reads the protected MIMIC objects, the R1B preflight aggregates, and the R1B cluster result. It uses the frozen 24-hour landmark definition, four-degree-of-freedom baseline spline, two-degree-of-freedom profile time interaction, proportional effects for the remaining covariates, patient-clustered covariance, and the prespecified earliest/midpoint/latest death-time scenarios. It writes aggregate R1C tables under `outputs/revision/`.

The script is gated by `R1C_EXECUTE_TIME_VARYING=YES`. It also verifies historical SHA-256 values. Because rebuilt stochastic/protected objects can be scientifically legitimate yet byte-different, a hash mismatch must trigger a documented provenance review. Do not remove or silently bypass the check.

`R/revision/R1C_render_figures_from_aggregates.R` reads only aggregate R1C curve and point files and renders Supplementary Figures S1-S2.

## K. Figure and table assembly

`python/stage12a_build_figures_tables.py` is the actual Stage 12A assembly source. It reads aggregate CSV outputs 27, 28, 29, 30, 31, 35, 36, and 38. It creates Figure 1, other aggregate figures, and editable Word tables in the original Stage 12A layout.

The script predates the R1B patient-clustered update. It must not be used alone to recreate the final R1F Table 3. The R1F table combines the frozen point estimates/RMST with R1B clustered confidence intervals and P values. Final Word document integration is not included in this source package.

## L. Safe validation

The following commands do not access clinical data when run from the package root:

```text
Rscript -e "for (f in list.files('R', pattern='[.]R$', recursive=TRUE, full.names=TRUE)) parse(file=f); for (f in list.files('tests', pattern='[.]R$', recursive=TRUE, full.names=TRUE)) parse(file=f)"
Rscript tests/run_lpa_em_synthetic_self_test.R
Rscript tests/test_probability_ensemble_synthetic.R
python tests/static_release_checks.py
```

Do not execute any other analysis script merely to test the release package. The release audit never connects to a database or runs a real-data model.

