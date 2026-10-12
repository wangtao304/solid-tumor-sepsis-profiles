-- Stage 12L-2 independent aggregate-only QC for first eligible ICU stays.
-- No identifier or patient-level row is returned.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '5min';

WITH eligible AS (
  SELECT
    b.subject_id,
    b.stay_id,
    b.mortality_28d,
    i.icu_intime
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base AS b
  JOIN mimiciv_derived.icustay_detail AS i
    ON i.subject_id = b.subject_id
   AND i.stay_id = b.stay_id
  WHERE b.cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
),
ranked AS (
  SELECT
    subject_id,
    mortality_28d,
    row_number() OVER (
      PARTITION BY subject_id
      ORDER BY icu_intime, stay_id
    ) AS eligible_stay_rank,
    min(icu_intime) OVER (PARTITION BY subject_id) AS minimum_icu_intime,
    icu_intime
  FROM eligible
),
first_stays AS (
  SELECT subject_id, mortality_28d
  FROM ranked
  WHERE eligible_stay_rank = 1
),
earliest_ties AS (
  SELECT subject_id
  FROM ranked
  WHERE icu_intime = minimum_icu_intime
  GROUP BY subject_id
  HAVING count(*) > 1
)
SELECT
  (SELECT count(*) FROM eligible) AS full_eligible_stays,
  (SELECT count(DISTINCT subject_id) FROM eligible) AS full_unique_subjects,
  (SELECT count(*) FROM first_stays) AS first_stay_rows,
  (
    SELECT count(*) - count(DISTINCT subject_id)
    FROM first_stays
  ) AS duplicate_subject_rows,
  (SELECT sum(mortality_28d) FROM first_stays) AS first_stay_deaths_28d,
  (SELECT count(*) FROM earliest_ties) AS subjects_with_tied_earliest_time;

COMMIT;
