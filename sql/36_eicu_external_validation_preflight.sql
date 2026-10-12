-- Stage 10C execution preflight. Aggregate output only; no identifiers returned.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '10min';

WITH adult AS (
  SELECT
    patientunitstayid,
    gender,
    CASE
      WHEN age = '> 89' THEN 90
      WHEN age ~ '^[0-9]+$' THEN age::integer
      ELSE NULL
    END AS age,
    unitdischargeoffset,
    hospitaldischargeoffset,
    hospitaldischargestatus
  FROM eicu_crd.patient
  WHERE CASE
    WHEN age = '> 89' THEN 90
    WHEN age ~ '^[0-9]+$' THEN age::integer
    ELSE NULL
  END >= 18
),
diagnosis_flags AS (
  SELECT
    patientunitstayid,
    bool_or(
      coalesce(icd9code, '') ~ '(^|[, ]+)(038|995\.9[12]|785\.52)'
      OR lower(coalesce(diagnosisstring, '')) LIKE '%sepsis%'
    ) AS sepsis_proxy,
    bool_or(
      upper(replace(split_part(coalesce(icd9code, ''), ',', 1), '.', ''))
        ~ '^(1[4-9][0-9])'
    ) AS solid_icd,
    bool_or(
      upper(replace(split_part(coalesce(icd9code, ''), ',', 1), '.', ''))
        ~ '^(19[6-9])'
    ) AS metastatic_icd,
    bool_or(
      upper(replace(split_part(coalesce(icd9code, ''), ',', 1), '.', ''))
        ~ '^(20[0-8])'
    ) AS hematologic_icd
  FROM eicu_crd.diagnosis
  GROUP BY patientunitstayid
),
apache AS (
  SELECT
    patientunitstayid,
    bool_or(coalesce(metastaticcancer, 0) = 1) AS metastatic_solid,
    bool_or(coalesce(lymphoma, 0) = 1 OR coalesce(leukemia, 0) = 1)
      AS hematologic_apache
  FROM eicu_crd.apachepredvar
  GROUP BY patientunitstayid
),
target AS (
  SELECT
    a.*,
    CASE
      WHEN d.metastatic_icd OR coalesce(p.metastatic_solid, false)
        THEN 'metastatic solid'
      ELSE 'nonmetastatic solid'
    END AS cancer_phenotype
  FROM adult a
  JOIN diagnosis_flags d USING (patientunitstayid)
  LEFT JOIN apache p USING (patientunitstayid)
  WHERE d.sepsis_proxy
    AND (d.solid_icd OR coalesce(p.metastatic_solid, false))
    AND NOT (d.hematologic_icd OR coalesce(p.hematologic_apache, false))
    AND a.unitdischargeoffset >= 1440
),
labs AS (
  SELECT
    t.patientunitstayid,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY l.labresult)
      FILTER (WHERE lower(l.labname) = 'wbc x 1000') AS wbc,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY l.labresult)
      FILTER (WHERE lower(l.labname) = 'platelets x 1000') AS platelet,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY l.labresult)
      FILTER (WHERE lower(l.labname) = 'lactate') AS lactate,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY l.labresult)
      FILTER (WHERE lower(l.labname) IN ('glucose', 'bedside glucose')) AS glucose,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY l.labresult)
      FILTER (WHERE lower(l.labname) = 'creatinine') AS creatinine,
    percentile_cont(0.5) WITHIN GROUP (ORDER BY l.labresult)
      FILTER (WHERE lower(l.labname) IN ('bicarbonate', 'hco3')) AS bicarbonate
  FROM target t
  LEFT JOIN eicu_crd.lab l
    ON l.patientunitstayid = t.patientunitstayid
   AND l.labresultoffset >= 0
   AND l.labresultoffset < 1440
   AND l.labresult IS NOT NULL
  GROUP BY t.patientunitstayid
),
joined AS (
  SELECT t.*, l.*
  FROM target t
  JOIN labs l USING (patientunitstayid)
)
SELECT
  count(*) AS target_n,
  count(*) FILTER (
    WHERE wbc IS NOT NULL AND platelet IS NOT NULL AND lactate IS NOT NULL
      AND glucose IS NOT NULL AND creatinine IS NOT NULL AND bicarbonate IS NOT NULL
  ) AS six_variable_complete_n,
  count(*) FILTER (WHERE age IS NOT NULL) AS age_available_n,
  count(*) FILTER (WHERE gender IN ('Male', 'Female')) AS binary_sex_available_n,
  count(*) FILTER (WHERE cancer_phenotype IS NOT NULL) AS phenotype_available_n,
  count(*) FILTER (WHERE hospitaldischargestatus IS NOT NULL) AS hospital_status_available_n,
  count(*) FILTER (
    WHERE hospitaldischargestatus = 'Expired'
      AND hospitaldischargeoffset > 0
      AND hospitaldischargeoffset <= 40320
  ) AS observed_in_hospital_death_by_28d_n,
  count(*) FILTER (
    WHERE hospitaldischargestatus = 'Alive'
      AND hospitaldischargeoffset > 0
      AND hospitaldischargeoffset < 40320
  ) AS discharged_alive_before_28d_unknown_vital_status_n,
  count(*) FILTER (WHERE hospitaldischargeoffset >= 40320) AS observed_in_hospital_through_28d_n
FROM joined;

COMMIT;
