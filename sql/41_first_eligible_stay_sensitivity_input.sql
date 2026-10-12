-- Stage 12L-2 controlled first-eligible-stay sensitivity input.
-- Identifiers and ICU admission time are streamed to R memory only for
-- deterministic linkage and first-stay selection. They are never exported.
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
  icu_intime,
  age,
  sex,
  cancer_phenotype,
  charlson_comorbidity_index,
  sofa,
  sapsii,
  mortality_28d AS event_28d,
  CASE
    WHEN mortality_28d = 0 THEN 28.0
    WHEN deathtime IS NOT NULL
         AND deathtime >= icu_intime
         AND deathtime <= icu_intime + INTERVAL '28 day'
      THEN LEAST(
        28.0,
        GREATEST(
          1.0 / 1440.0,
          EXTRACT(EPOCH FROM (deathtime - icu_intime)) / 86400.0
        )
      )
    WHEN dod IS NOT NULL
      THEN LEAST(
        28.0,
        GREATEST(
          0.5,
          (dod - icu_intime::date)::double precision + 0.5
        )
      )
    ELSE NULL
  END AS time_28d_days
FROM cohort
ORDER BY stay_id;

COMMIT;
