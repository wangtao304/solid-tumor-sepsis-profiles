-- Stage 3 aggregate-only variable audit.
-- No table is created and no patient/stay identifier is returned.
-- Continuous summaries use each stay's median within [ICU intime, +24h),
-- followed by cohort-level median and Q1-Q3 across available stays.
WITH cancer_by_admission AS (
    SELECT hadm_id
    FROM mimiciv_hosp.diagnoses_icd
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
    SELECT t.stay_id,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.wbc) FILTER (WHERE x.wbc IS NOT NULL) AS wbc,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.platelet) FILTER (WHERE x.platelet IS NOT NULL) AS platelet
    FROM target t
    LEFT JOIN mimiciv_derived.complete_blood_count x
      ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
differential AS (
    SELECT t.stay_id,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.neutrophils) FILTER (WHERE x.neutrophils IS NOT NULL) AS neutrophil,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.lymphocytes) FILTER (WHERE x.lymphocytes IS NOT NULL) AS lymphocyte,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.neutrophils/NULLIF(x.lymphocytes,0))
             FILTER (WHERE x.neutrophils IS NOT NULL AND x.lymphocytes>0) AS nlr
    FROM target t
    LEFT JOIN mimiciv_derived.blood_differential x
      ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
blood_gas AS (
    SELECT t.stay_id,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.lactate) FILTER (WHERE x.lactate IS NOT NULL) AS lactate
    FROM target t
    LEFT JOIN mimiciv_derived.bg x
      ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
chemistry AS (
    SELECT t.stay_id,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.glucose) FILTER (WHERE x.glucose IS NOT NULL) AS glucose,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.creatinine) FILTER (WHERE x.creatinine IS NOT NULL) AS creatinine,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY x.bicarbonate) FILTER (WHERE x.bicarbonate IS NOT NULL) AS bicarbonate
    FROM target t
    LEFT JOIN mimiciv_derived.chemistry x
      ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
scores AS (
    SELECT t.stay_id, max(s.sofa)::double precision AS sofa, max(p.sapsii)::double precision AS sapsii
    FROM target t
    LEFT JOIN mimiciv_derived.first_day_sofa s USING (stay_id)
    LEFT JOIN mimiciv_derived.sapsii p USING (stay_id)
    GROUP BY t.stay_id
),
stay_values AS (
    SELECT t.stay_id, c.wbc, d.neutrophil, d.lymphocyte, c.platelet, d.nlr,
           b.lactate, h.glucose, h.creatinine, h.bicarbonate, s.sofa, s.sapsii
    FROM target t
    LEFT JOIN cbc c USING (stay_id)
    LEFT JOIN differential d USING (stay_id)
    LEFT JOIN blood_gas b USING (stay_id)
    LEFT JOIN chemistry h USING (stay_id)
    LEFT JOIN scores s USING (stay_id)
),
long AS (
    SELECT v.variable, v.value
    FROM stay_values x
    CROSS JOIN LATERAL (VALUES
      ('WBC',x.wbc),('neutrophil',x.neutrophil),('lymphocyte',x.lymphocyte),('platelet',x.platelet),
      ('NLR',x.nlr),('lactate',x.lactate),('glucose',x.glucose),('creatinine',x.creatinine),
      ('bicarbonate',x.bicarbonate),('SOFA',x.sofa),('SAPS-II',x.sapsii)
    ) v(variable,value)
),
summary AS (
    SELECT variable, count(value)::bigint AS available_n,
           round(100.0*(count(*)-count(value))/NULLIF(count(*),0),2) AS missing_percent,
           percentile_cont(0.5) WITHIN GROUP (ORDER BY value) FILTER (WHERE value IS NOT NULL) AS median,
           percentile_cont(0.25) WITHIN GROUP (ORDER BY value) FILTER (WHERE value IS NOT NULL) AS q1,
           percentile_cont(0.75) WITHIN GROUP (ORDER BY value) FILTER (WHERE value IS NOT NULL) AS q3
    FROM long GROUP BY variable
)
SELECT variable, available_n, missing_percent, round(median::numeric,2) AS median,
       concat(round(q1::numeric,2),'-',round(q3::numeric,2)) AS "IQR"
FROM summary
ORDER BY array_position(ARRAY['WBC','neutrophil','lymphocyte','platelet','NLR','lactate','glucose','creatinine','bicarbonate','SOFA','SAPS-II'],variable);
