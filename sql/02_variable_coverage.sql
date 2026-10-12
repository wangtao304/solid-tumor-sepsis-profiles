-- Stage 2 aggregate-only audit. Availability is at least one value in [ICU intime, +24h).
WITH cancer_by_admission AS (
    SELECT hadm_id FROM mimiciv_hosp.diagnoses_icd
    WHERE (icd_version = 9 AND upper(replace(icd_code, '.', '')) ~ '^(1[4-9][0-9]|20[0-8])')
       OR (icd_version = 10 AND upper(replace(icd_code, '.', '')) ~ '^C(0[0-9]|[1-6][0-9]|7[0-9AB]|8[0-9]|9[0-6])')
    GROUP BY hadm_id
),
target AS (
    SELECT i.stay_id, i.hadm_id, i.icu_intime
    FROM mimiciv_derived.icustay_detail i
    JOIN (SELECT DISTINCT stay_id FROM mimiciv_derived.sepsis3 WHERE sepsis3 = true) s USING (stay_id)
    JOIN cancer_by_admission c USING (hadm_id)
    WHERE i.admission_age >= 18 AND i.los_icu >= 1
),
cbc AS (
    SELECT t.stay_id, bool_or(x.wbc IS NOT NULL) AS wbc, bool_or(x.platelet IS NOT NULL) AS platelet
    FROM target t JOIN mimiciv_derived.complete_blood_count x ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
diff AS (
    SELECT t.stay_id, bool_or(x.neutrophils_abs IS NOT NULL OR x.neutrophils IS NOT NULL) AS neutrophil,
           bool_or(x.lymphocytes_abs IS NOT NULL OR x.lymphocytes IS NOT NULL) AS lymphocyte
    FROM target t JOIN mimiciv_derived.blood_differential x ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
blood_gas AS (
    SELECT t.stay_id, bool_or(x.lactate IS NOT NULL) AS lactate
    FROM target t JOIN mimiciv_derived.bg x ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
chem AS (
    SELECT t.stay_id, bool_or(x.glucose IS NOT NULL) AS glucose, bool_or(x.creatinine IS NOT NULL) AS creatinine,
           bool_or(x.bicarbonate IS NOT NULL) AS bicarbonate
    FROM target t JOIN mimiciv_derived.chemistry x ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
availability AS (
    SELECT t.stay_id, coalesce(c.wbc,false) AS wbc, coalesce(d.neutrophil,false) AS neutrophil,
           coalesce(d.lymphocyte,false) AS lymphocyte, coalesce(c.platelet,false) AS platelet,
           coalesce(b.lactate,false) AS lactate, coalesce(h.glucose,false) AS glucose,
           coalesce(h.creatinine,false) AS creatinine, coalesce(h.bicarbonate,false) AS bicarbonate,
           (s.sofa IS NOT NULL) AS sofa, (p.sapsii IS NOT NULL) AS sapsii
    FROM target t LEFT JOIN cbc c USING(stay_id) LEFT JOIN diff d USING(stay_id)
    LEFT JOIN blood_gas b USING(stay_id) LEFT JOIN chem h USING(stay_id)
    LEFT JOIN mimiciv_derived.first_day_sofa s USING(stay_id) LEFT JOIN mimiciv_derived.sapsii p USING(stay_id)
),
long AS (
    SELECT a.stay_id, v.variable, v.available FROM availability a
    CROSS JOIN LATERAL (VALUES ('WBC',a.wbc),('neutrophil',a.neutrophil),('lymphocyte',a.lymphocyte),('platelet',a.platelet),('lactate',a.lactate),('glucose',a.glucose),('creatinine',a.creatinine),('bicarbonate',a.bicarbonate),('SOFA',a.sofa),('SAPS-II',a.sapsii)) v(variable,available)
)
SELECT variable, COUNT(*) FILTER (WHERE available)::bigint AS available_n,
       COUNT(*) FILTER (WHERE NOT available)::bigint AS missing_n,
       round(100.0 * COUNT(*) FILTER (WHERE NOT available) / NULLIF(COUNT(*),0), 2) AS missing_percent
FROM long GROUP BY variable
ORDER BY array_position(ARRAY['WBC','neutrophil','lymphocyte','platelet','lactate','glucose','creatinine','bicarbonate','SOFA','SAPS-II'], variable);
