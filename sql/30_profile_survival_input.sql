-- Stage 9B controlled survival input.
-- Patient/stay identifiers are streamed to R memory only for deterministic
-- linkage validation and are never written to an output table or file.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '5min';

WITH cohort AS (
  SELECT
    b.stay_id,
    b.wbc,
    b.platelet_count,
    b.lactate,
    b.glucose,
    b.creatinine,
    b.bicarbonate,
    b.age,
    b.sex,
    b.cancer_phenotype,
    b.charlson_comorbidity_index,
    b.sofa,
    b.sapsii,
    b.vasopressor_0_24,
    b.mechanical_ventilation_0_24,
    b.mortality_28d,
    b.icu_mortality,
    b.hospital_mortality,
    i.icu_intime,
    i.icu_outtime,
    a.dischtime,
    a.deathtime,
    p.dod
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base b
  JOIN mimiciv_derived.icustay_detail i USING (stay_id)
  JOIN mimiciv_hosp.admissions a
    ON a.subject_id = i.subject_id
   AND a.hadm_id = i.hadm_id
  JOIN mimiciv_hosp.patients p
    ON p.subject_id = i.subject_id
  WHERE b.cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
)
SELECT
  stay_id,
  age,
  sex,
  cancer_phenotype,
  charlson_comorbidity_index,
  sofa,
  sapsii,
  vasopressor_0_24,
  mechanical_ventilation_0_24,
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
  END AS time_28d_days,
  icu_mortality AS event_icu,
  CASE
    WHEN icu_mortality = 1
      THEN GREATEST(
        1.0 / 1440.0,
        EXTRACT(EPOCH FROM (deathtime - icu_intime)) / 86400.0
      )
    ELSE EXTRACT(EPOCH FROM (icu_outtime - icu_intime)) / 86400.0
  END AS time_icu_days,
  hospital_mortality AS event_hospital,
  CASE
    WHEN hospital_mortality = 1
      THEN GREATEST(
        1.0 / 1440.0,
        EXTRACT(EPOCH FROM (deathtime - icu_intime)) / 86400.0
      )
    ELSE EXTRACT(EPOCH FROM (dischtime - icu_intime)) / 86400.0
  END AS time_hospital_days,
  (
    wbc IS NOT NULL
    AND platelet_count IS NOT NULL
    AND lactate IS NOT NULL
    AND glucose IS NOT NULL
    AND creatinine IS NOT NULL
    AND bicarbonate IS NOT NULL
  ) AS complete_case_primary
FROM cohort
ORDER BY stay_id;

COMMIT;
