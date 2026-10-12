-- Stage 10A: aggregate-only eICU external-validation feasibility audit.
-- Read only. Returns no patientunitstayid and does not use mortality outcomes.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '10min';

WITH adult AS (
  SELECT patientunitstayid, unitdischargeoffset
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
        ~ '^(20[0-8])'
    ) AS hematologic_icd
  FROM eicu_crd.diagnosis
  GROUP BY patientunitstayid
),
apache_cancer AS (
  SELECT
    patientunitstayid,
    bool_or(coalesce(metastaticcancer, 0) = 1) AS metastatic_solid,
    bool_or(coalesce(lymphoma, 0) = 1 OR coalesce(leukemia, 0) = 1)
      AS hematologic_apache
  FROM eicu_crd.apachepredvar
  GROUP BY patientunitstayid
),
sepsis AS (
  SELECT a.patientunitstayid, a.unitdischargeoffset
  FROM adult a
  JOIN diagnosis_flags d USING (patientunitstayid)
  WHERE d.sepsis_proxy
),
solid_tumor AS (
  SELECT s.patientunitstayid, s.unitdischargeoffset
  FROM sepsis s
  JOIN diagnosis_flags d USING (patientunitstayid)
  LEFT JOIN apache_cancer c USING (patientunitstayid)
  WHERE (d.solid_icd OR coalesce(c.metastatic_solid, false))
    AND NOT (d.hematologic_icd OR coalesce(c.hematologic_apache, false))
),
target AS (
  SELECT patientunitstayid
  FROM solid_tumor
  WHERE unitdischargeoffset >= 1440
),
lab_availability AS (
  SELECT
    t.patientunitstayid,
    bool_or(lower(l.labname) = 'wbc x 1000' AND l.labresult IS NOT NULL) AS wbc,
    bool_or(lower(l.labname) = 'platelets x 1000' AND l.labresult IS NOT NULL) AS platelet,
    bool_or(lower(l.labname) = 'lactate' AND l.labresult IS NOT NULL) AS lactate,
    bool_or(
      lower(l.labname) IN ('glucose', 'bedside glucose')
      AND l.labresult IS NOT NULL
    ) AS glucose,
    bool_or(lower(l.labname) = 'creatinine' AND l.labresult IS NOT NULL) AS creatinine,
    bool_or(
      lower(l.labname) IN ('bicarbonate', 'hco3')
      AND l.labresult IS NOT NULL
    ) AS bicarbonate
  FROM target t
  LEFT JOIN eicu_crd.lab l
    ON l.patientunitstayid = t.patientunitstayid
   AND l.labresultoffset >= 0
   AND l.labresultoffset < 1440
  GROUP BY t.patientunitstayid
),
coverage AS (
  SELECT
    v.variable,
    count(*)::bigint AS denominator_n,
    count(*) FILTER (WHERE coalesce(v.available, false))::bigint AS available_n,
    count(*) FILTER (WHERE NOT coalesce(v.available, false))::bigint AS missing_n
  FROM lab_availability a
  CROSS JOIN LATERAL (
    VALUES
      ('WBC', a.wbc),
      ('platelet', a.platelet),
      ('lactate', a.lactate),
      ('glucose', a.glucose),
      ('creatinine', a.creatinine),
      ('bicarbonate', a.bicarbonate)
  ) v(variable, available)
  GROUP BY v.variable
),
rows AS (
  SELECT
    1 AS display_order,
    'cohort_attrition'::text AS section,
    'Adult ICU admissions'::text AS metric,
    count(*)::bigint AS cohort_n,
    NULL::bigint AS denominator_n,
    NULL::bigint AS available_n,
    NULL::bigint AS missing_n,
    NULL::numeric AS missing_percent,
    'Age >=18 years; unit of analysis is ICU stay'::text AS definition_note
  FROM adult

  UNION ALL
  SELECT
    2, 'cohort_attrition', 'Sepsis-3 screening cohort', count(*),
    NULL, NULL, NULL, NULL,
    'eICU has no installed Sepsis-3 derived table; diagnostic ICD-9/text sepsis proxy used for feasibility only'
  FROM sepsis

  UNION ALL
  SELECT
    3, 'cohort_attrition', 'Solid tumor + sepsis screening cohort', count(*),
    NULL, NULL, NULL, NULL,
    'Solid ICD-9 140-199 or APACHE metastatic cancer; hematologic or mixed malignancy excluded'
  FROM solid_tumor

  UNION ALL
  SELECT
    4, 'cohort_attrition', 'Solid tumor + sepsis with ICU LOS >=24h', count(*),
    NULL, NULL, NULL, NULL,
    'unitdischargeoffset >=1440 minutes'
  FROM target

  UNION ALL
  SELECT
    10 + array_position(
      ARRAY['WBC', 'platelet', 'lactate', 'glucose', 'creatinine', 'bicarbonate'],
      variable
    ),
    'variable_coverage', variable, NULL, denominator_n, available_n, missing_n,
    round(100.0 * missing_n / NULLIF(denominator_n, 0), 2),
    'At least one nonmissing harmonized laboratory value during ICU hours 0-24'
  FROM coverage
)
SELECT
  section,
  metric,
  cohort_n,
  denominator_n,
  available_n,
  missing_n,
  missing_percent,
  definition_note
FROM rows
ORDER BY display_order;

COMMIT;
