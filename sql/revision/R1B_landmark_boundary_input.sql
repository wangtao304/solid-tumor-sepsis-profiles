-- Reviewer Revision R1B: read-only 24-hour landmark/death-boundary audit input.
-- Identifiers and dates/timestamps are streamed only to controlled R memory for
-- deterministic linkage and boundary classification. No patient-level output
-- is written by the R1B workflow.

BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '5min';

WITH cohort AS (
  SELECT
    b.subject_id,
    b.stay_id,
    b.age,
    b.sex,
    b.cancer_phenotype,
    b.charlson_comorbidity_index,
    b.sofa,
    b.sapsii,
    b.mortality_28d,
    i.icu_intime,
    i.icu_outtime,
    a.dischtime,
    a.deathtime,
    p.dod
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base AS b
  JOIN mimiciv_derived.icustay_detail AS i
    ON i.subject_id = b.subject_id
   AND i.stay_id = b.stay_id
  JOIN mimiciv_hosp.admissions AS a
    ON a.subject_id = i.subject_id
   AND a.hadm_id = i.hadm_id
  JOIN mimiciv_hosp.patients AS p
    ON p.subject_id = i.subject_id
  WHERE b.cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
)
SELECT
  subject_id,
  stay_id,
  age,
  sex,
  cancer_phenotype,
  charlson_comorbidity_index,
  sofa,
  sapsii,
  mortality_28d AS event_28d,
  icu_intime,
  icu_outtime,
  dischtime,
  deathtime,
  dod AS death_date
FROM cohort
ORDER BY stay_id;

COMMIT;
