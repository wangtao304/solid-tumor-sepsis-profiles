-- Stage 12L-1: aggregate-only cohort reproduction gate.
-- No patient or stay identifier is returned.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '5min';

WITH cancer_flags AS (
    SELECT
        hadm_id,
        bool_or(
            (icd_version = 9 AND upper(replace(icd_code, '.', '')) ~ '^20[0-8]') OR
            (icd_version = 10 AND upper(replace(icd_code, '.', '')) ~ '^C(8[1-9]|9[0-6])')
        ) AS hematologic,
        bool_or(
            (icd_version = 9 AND upper(replace(icd_code, '.', '')) ~ '^19[6-9]') OR
            (icd_version = 10 AND upper(replace(icd_code, '.', '')) ~ '^C(7[7-9B]|80)')
        ) AS metastatic_solid,
        bool_or(
            (icd_version = 9 AND upper(replace(icd_code, '.', '')) ~ '^1(4[0-9]|5[0-9]|6[0-9]|7[0-9]|8[0-9]|9[0-5])') OR
            (icd_version = 10 AND upper(replace(icd_code, '.', '')) ~ '^C(0[0-9]|[1-6][0-9]|7[0-6A])')
        ) AS nonmetastatic_solid
    FROM mimiciv_hosp.diagnoses_icd
    GROUP BY hadm_id
),
raw_target AS (
    SELECT
        i.subject_id,
        i.stay_id,
        CASE
            WHEN c.hematologic AND (c.metastatic_solid OR c.nonmetastatic_solid) THEN 'mixed'
            WHEN c.hematologic THEN 'hematologic'
            WHEN c.metastatic_solid THEN 'metastatic solid'
            ELSE 'nonmetastatic solid'
        END AS cancer_phenotype
    FROM mimiciv_derived.icustay_detail AS i
    JOIN (
        SELECT DISTINCT stay_id
        FROM mimiciv_derived.sepsis3
        WHERE sepsis3 = true
    ) AS s USING (stay_id)
    JOIN cancer_flags AS c USING (hadm_id)
    WHERE i.admission_age >= 18
      AND i.los_icu >= 1
      AND (c.hematologic OR c.metastatic_solid OR c.nonmetastatic_solid)
),
raw_solid AS (
    SELECT subject_id, stay_id, cancer_phenotype
    FROM raw_target
    WHERE cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
),
frozen_solid AS (
    SELECT subject_id, stay_id, cancer_phenotype
    FROM mimic_cancer_sepsis.immune_metabolic_lca_base
    WHERE cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
)
SELECT
    (SELECT count(*) FROM raw_solid) AS raw_rule_eligible_icu_stays,
    (SELECT count(DISTINCT stay_id) FROM raw_solid) AS raw_rule_unique_stay_ids,
    (SELECT count(*) FROM frozen_solid) AS frozen_base_eligible_icu_stays,
    (SELECT count(DISTINCT stay_id) FROM frozen_solid) AS frozen_base_unique_stay_ids,
    (
        SELECT count(*)
        FROM (
            SELECT subject_id, stay_id, cancer_phenotype FROM raw_solid
            EXCEPT
            SELECT subject_id, stay_id, cancer_phenotype FROM frozen_solid
        ) AS raw_only
    ) AS raw_rows_absent_from_frozen_base,
    (
        SELECT count(*)
        FROM (
            SELECT subject_id, stay_id, cancer_phenotype FROM frozen_solid
            EXCEPT
            SELECT subject_id, stay_id, cancer_phenotype FROM raw_solid
        ) AS frozen_only
    ) AS frozen_rows_absent_from_raw_rule;

COMMIT;
