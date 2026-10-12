-- Stage 4A: create the database-resident pre-LCA base.
-- Patient/stay identifiers remain inside the approved database.
BEGIN;

CREATE SCHEMA IF NOT EXISTS mimic_cancer_sepsis;

CREATE TABLE mimic_cancer_sepsis.immune_metabolic_lca_base AS
WITH cancer_flags AS (
    SELECT hadm_id,
      bool_or((icd_version=9 AND upper(replace(icd_code,'.','')) ~ '^20[0-8]') OR
              (icd_version=10 AND upper(replace(icd_code,'.','')) ~ '^C(8[1-9]|9[0-6])')) AS hematologic,
      bool_or((icd_version=9 AND upper(replace(icd_code,'.','')) ~ '^19[6-9]') OR
              (icd_version=10 AND upper(replace(icd_code,'.','')) ~ '^C(7[7-9B]|80)')) AS metastatic_solid,
      bool_or((icd_version=9 AND upper(replace(icd_code,'.','')) ~ '^1(4[0-9]|5[0-9]|6[0-9]|7[0-9]|8[0-9]|9[0-5])') OR
              (icd_version=10 AND upper(replace(icd_code,'.','')) ~ '^C(0[0-9]|[1-6][0-9]|7[0-6A])')) AS nonmetastatic_solid
    FROM mimiciv_hosp.diagnoses_icd
    GROUP BY hadm_id
),
target AS (
    SELECT i.subject_id, i.hadm_id, i.stay_id, i.icu_intime, i.icu_outtime,
           i.admission_age AS age, i.gender AS sex,
           CASE
             WHEN c.hematologic AND (c.metastatic_solid OR c.nonmetastatic_solid) THEN 'mixed'
             WHEN c.hematologic THEN 'hematologic'
             WHEN c.metastatic_solid THEN 'metastatic solid'
             ELSE 'nonmetastatic solid'
           END AS cancer_phenotype
    FROM mimiciv_derived.icustay_detail i
    JOIN (SELECT DISTINCT stay_id FROM mimiciv_derived.sepsis3 WHERE sepsis3=true) s USING (stay_id)
    JOIN cancer_flags c USING (hadm_id)
    WHERE i.admission_age>=18 AND i.los_icu>=1
      AND (c.hematologic OR c.metastatic_solid OR c.nonmetastatic_solid)
),
cbc AS (
    SELECT t.stay_id,
      percentile_cont(0.5) WITHIN GROUP (ORDER BY x.platelet) FILTER (WHERE x.platelet IS NOT NULL) AS platelet_count,
      percentile_cont(0.5) WITHIN GROUP (ORDER BY x.wbc) FILTER (WHERE x.wbc IS NOT NULL) AS wbc
    FROM target t
    LEFT JOIN mimiciv_derived.complete_blood_count x
      ON x.hadm_id=t.hadm_id AND x.charttime>=t.icu_intime AND x.charttime<t.icu_intime+interval '24 hour'
    GROUP BY t.stay_id
),
differential AS (
    SELECT t.stay_id,
      percentile_cont(0.5) WITHIN GROUP (ORDER BY x.neutrophils_abs) FILTER (WHERE x.neutrophils_abs IS NOT NULL) AS neutrophil_abs,
      percentile_cont(0.5) WITHIN GROUP (ORDER BY x.lymphocytes_abs) FILTER (WHERE x.lymphocytes_abs IS NOT NULL) AS lymphocyte_abs,
      percentile_cont(0.5) WITHIN GROUP (ORDER BY x.neutrophils_abs/NULLIF(x.lymphocytes_abs,0))
        FILTER (WHERE x.neutrophils_abs IS NOT NULL AND x.lymphocytes_abs>0) AS nlr
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
    SELECT t.stay_id, max(s.sofa) AS sofa, max(p.sapsii) AS sapsii
    FROM target t
    LEFT JOIN mimiciv_derived.first_day_sofa s USING (stay_id)
    LEFT JOIN mimiciv_derived.sapsii p USING (stay_id)
    GROUP BY t.stay_id
)
SELECT t.subject_id, t.stay_id,
       d.neutrophil_abs, d.lymphocyte_abs, c.platelet_count,
       b.lactate, h.glucose, h.creatinine, h.bicarbonate,
       c.wbc, d.nlr, s.sofa, s.sapsii,
       t.age, t.sex, t.cancer_phenotype,
       ch.charlson_comorbidity_index,
       CASE
         WHEN a.deathtime IS NOT NULL AND a.deathtime>=t.icu_intime AND a.deathtime<=t.icu_intime+interval '28 day' THEN 1
         WHEN a.deathtime IS NULL AND p.dod IS NOT NULL AND p.dod>=t.icu_intime::date AND p.dod<=(t.icu_intime+interval '28 day')::date THEN 1
         ELSE 0
       END::smallint AS mortality_28d,
       CASE WHEN a.deathtime IS NOT NULL AND a.deathtime>=t.icu_intime AND a.deathtime<=t.icu_outtime THEN 1 ELSE 0 END::smallint AS icu_mortality,
       a.hospital_expire_flag::smallint AS hospital_mortality
FROM target t
LEFT JOIN cbc c USING (stay_id)
LEFT JOIN differential d USING (stay_id)
LEFT JOIN blood_gas b USING (stay_id)
LEFT JOIN chemistry h USING (stay_id)
LEFT JOIN scores s USING (stay_id)
LEFT JOIN mimiciv_derived.charlson ch USING (subject_id,hadm_id)
LEFT JOIN mimiciv_hosp.admissions a USING (subject_id,hadm_id)
LEFT JOIN mimiciv_hosp.patients p USING (subject_id);

ALTER TABLE mimic_cancer_sepsis.immune_metabolic_lca_base
  ADD CONSTRAINT immune_metabolic_lca_base_pk PRIMARY KEY (stay_id);

COMMENT ON TABLE mimic_cancer_sepsis.immune_metabolic_lca_base IS
  'Stage 4A database-resident pre-LCA base; patient-level export prohibited.';

COMMIT;
