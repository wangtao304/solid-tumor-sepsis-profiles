-- Mutually exclusive cancer categories in the final Stage 2 target cohort.
WITH cancer_flags AS (
    SELECT hadm_id,
      bool_or((icd_version=9 AND upper(replace(icd_code,'.','')) ~ '^20[0-8]') OR (icd_version=10 AND upper(replace(icd_code,'.','')) ~ '^C(8[1-9]|9[0-6])')) AS hematologic,
      bool_or((icd_version=9 AND upper(replace(icd_code,'.','')) ~ '^19[6-9]') OR (icd_version=10 AND upper(replace(icd_code,'.','')) ~ '^C(7[7-9B]|80)')) AS metastatic_solid,
      bool_or((icd_version=9 AND upper(replace(icd_code,'.','')) ~ '^1(4[0-9]|5[0-9]|6[0-9]|7[0-9]|8[0-9]|9[0-5])') OR (icd_version=10 AND upper(replace(icd_code,'.','')) ~ '^C(0[0-9]|[1-6][0-9]|7[0-6A])')) AS nonmetastatic_solid
    FROM mimiciv_hosp.diagnoses_icd GROUP BY hadm_id
),
target AS (
    SELECT i.stay_id, c.*
    FROM mimiciv_derived.icustay_detail i
    JOIN (SELECT DISTINCT stay_id FROM mimiciv_derived.sepsis3 WHERE sepsis3=true) s USING (stay_id)
    JOIN cancer_flags c USING (hadm_id)
    WHERE i.admission_age>=18 AND i.los_icu>=1 AND (c.hematologic OR c.metastatic_solid OR c.nonmetastatic_solid)
),
classified AS (
    SELECT CASE WHEN hematologic AND (metastatic_solid OR nonmetastatic_solid) THEN 'mixed'
                WHEN hematologic THEN 'hematologic'
                WHEN metastatic_solid THEN 'metastatic solid'
                ELSE 'nonmetastatic solid' END AS phenotype
    FROM target
)
SELECT phenotype, COUNT(*)::bigint AS n, round(100.0*COUNT(*)/sum(COUNT(*)) OVER (),2) AS percent
FROM classified GROUP BY phenotype
ORDER BY array_position(ARRAY['hematologic','metastatic solid','nonmetastatic solid','mixed'],phenotype);
