-- Stage 2 aggregate-only eICU audit. Creates no tables and returns no stay identifiers.
WITH adult AS (
  SELECT patientunitstayid, unitdischargeoffset FROM eicu_crd.patient
  WHERE CASE WHEN age='> 89' THEN 90 WHEN age ~ '^[0-9]+$' THEN age::int ELSE NULL END >=18
),
dx AS (
  SELECT patientunitstayid,
    bool_or(upper(replace(split_part(coalesce(icd9code,''),',',1),'.','')) ~ '^(1[4-9][0-9]|20[0-8])') AS cancer_icd,
    bool_or(coalesce(icd9code,'') ~ '(^|[, ]+)(038|995\.9[12]|785\.52)' OR lower(coalesce(diagnosisstring,'')) LIKE '%sepsis%') AS sepsis
  FROM eicu_crd.diagnosis GROUP BY patientunitstayid
),
apache_cancer AS (
  SELECT patientunitstayid, bool_or(coalesce(lymphoma,0)=1 OR coalesce(leukemia,0)=1 OR coalesce(metastaticcancer,0)=1) AS cancer
  FROM eicu_crd.apachepredvar GROUP BY patientunitstayid
),
candidates AS (
  SELECT a.patientunitstayid, a.unitdischargeoffset FROM adult a JOIN dx d USING(patientunitstayid)
  LEFT JOIN apache_cancer c USING(patientunitstayid)
  WHERE d.sepsis AND (d.cancer_icd OR coalesce(c.cancer,false))
),
target AS (SELECT * FROM candidates WHERE unitdischargeoffset>=1440),
lab_availability AS (
  SELECT t.patientunitstayid,
    bool_or(lower(l.labname)='wbc x 1000' AND l.labresult IS NOT NULL) AS wbc,
    bool_or(lower(l.labname)='-polys' AND l.labresult IS NOT NULL) AS neutrophil,
    bool_or(lower(l.labname)='-lymphs' AND l.labresult IS NOT NULL) AS lymphocyte,
    bool_or(lower(l.labname)='platelets x 1000' AND l.labresult IS NOT NULL) AS platelet,
    bool_or(lower(l.labname)='lactate' AND l.labresult IS NOT NULL) AS lactate,
    bool_or(lower(l.labname) IN ('glucose','bedside glucose') AND l.labresult IS NOT NULL) AS glucose,
    bool_or(lower(l.labname)='creatinine' AND l.labresult IS NOT NULL) AS creatinine,
    bool_or(lower(l.labname) IN ('bicarbonate','hco3') AND l.labresult IS NOT NULL) AS bicarbonate
  FROM target t LEFT JOIN eicu_crd.lab l ON l.patientunitstayid=t.patientunitstayid AND l.labresultoffset>=0 AND l.labresultoffset<1440
  GROUP BY t.patientunitstayid
),
coverage AS (
  SELECT a.patientunitstayid,v.variable,v.available FROM lab_availability a
  CROSS JOIN LATERAL (VALUES ('WBC',a.wbc),('neutrophil',a.neutrophil),('lymphocyte',a.lymphocyte),('platelet',a.platelet),('lactate',a.lactate),('glucose',a.glucose),('creatinine',a.creatinine),('bicarbonate',a.bicarbonate)) v(variable,available)
),
rows AS (
  SELECT 1 ord,'count' row_type,'Adult ICU admissions' metric,COUNT(*)::bigint value,NULL::bigint denominator_n,NULL::bigint missing_n,NULL::numeric missing_percent FROM adult
  UNION ALL SELECT 2,'count','Cancer + sepsis candidates',COUNT(*),NULL,NULL,NULL FROM candidates
  UNION ALL SELECT 3,'count','Cancer + sepsis with ICU LOS >=24h',COUNT(*),NULL,NULL,NULL FROM target
  UNION ALL
  SELECT 10+array_position(ARRAY['WBC','neutrophil','lymphocyte','platelet','lactate','glucose','creatinine','bicarbonate'],variable),
    'variable_coverage',variable,COUNT(*) FILTER(WHERE available),COUNT(*),COUNT(*) FILTER(WHERE NOT coalesce(available,false)),
    round(100.0*COUNT(*) FILTER(WHERE NOT coalesce(available,false))/NULLIF(COUNT(*),0),2)
  FROM coverage GROUP BY variable
)
SELECT row_type,metric,value,denominator_n,missing_n,missing_percent FROM rows ORDER BY ord;
