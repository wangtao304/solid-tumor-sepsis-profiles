-- Stage 2 aggregate-only audit. Creates no tables and returns no stay identifiers.
WITH adult_icu AS (
    SELECT stay_id, hadm_id, icu_intime, los_icu
    FROM mimiciv_derived.icustay_detail
    WHERE admission_age >= 18
),
cancer_by_admission AS (
    SELECT hadm_id
    FROM mimiciv_hosp.diagnoses_icd
    WHERE (icd_version = 9 AND upper(replace(icd_code, '.', '')) ~ '^(1[4-9][0-9]|20[0-8])')
       OR (icd_version = 10 AND upper(replace(icd_code, '.', '')) ~ '^C(0[0-9]|[1-6][0-9]|7[0-9AB]|8[0-9]|9[0-6])')
    GROUP BY hadm_id
),
sepsis AS (
    SELECT DISTINCT stay_id
    FROM mimiciv_derived.sepsis3
    WHERE sepsis3 = true
),
counts AS (
    SELECT 1 AS ord, 'Adult ICU admissions' AS stage, COUNT(*)::bigint AS n FROM adult_icu
    UNION ALL
    SELECT 2, 'Sepsis-3 ICU admissions', COUNT(*) FROM adult_icu a JOIN sepsis s USING (stay_id)
    UNION ALL
    SELECT 3, 'Cancer + sepsis ICU admissions', COUNT(*) FROM adult_icu a JOIN sepsis s USING (stay_id) JOIN cancer_by_admission c USING (hadm_id)
    UNION ALL
    SELECT 4, 'Cancer + sepsis with ICU LOS >=24h', COUNT(*) FROM adult_icu a JOIN sepsis s USING (stay_id) JOIN cancer_by_admission c USING (hadm_id) WHERE a.los_icu >= 1
),
with_denominators AS (
    SELECT *, lag(n) OVER (ORDER BY ord) AS previous_n, max(n) FILTER (WHERE ord = 1) OVER () AS adult_n FROM counts
)
SELECT stage, n,
       round(100.0 * n / NULLIF(previous_n, 0), 2) AS retained_percent_previous,
       round(100.0 * n / NULLIF(adult_n, 0), 2) AS retained_percent_adult
FROM with_denominators ORDER BY ord;
