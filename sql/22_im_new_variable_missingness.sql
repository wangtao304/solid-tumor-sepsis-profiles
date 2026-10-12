-- Stage 7A aggregate-only missingness and distribution QC.
-- The solid-tumor cohort is filtered virtually from the existing database base;
-- no table is created and no identifier or patient-level row is returned.
WITH target AS (
    SELECT wbc,
           platelet_count,
           lactate,
           glucose,
           creatinine,
           bicarbonate,
           sofa,
           sapsii
    FROM mimic_cancer_sepsis.immune_metabolic_lca_base
    WHERE cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
),
long_values AS (
    SELECT v.variable,
           v.domain,
           v.role,
           v.sort_order,
           v.value
    FROM target t
    CROSS JOIN LATERAL (VALUES
        ('WBC', 'inflammatory', 'primary', 1, t.wbc::double precision),
        ('platelet', 'inflammatory', 'primary', 2, t.platelet_count::double precision),
        ('lactate', 'metabolic', 'primary', 3, t.lactate::double precision),
        ('glucose', 'metabolic', 'primary', 4, t.glucose::double precision),
        ('creatinine', 'metabolic', 'primary', 5, t.creatinine::double precision),
        ('bicarbonate', 'metabolic', 'primary', 6, t.bicarbonate::double precision),
        ('SOFA', 'severity', 'sensitivity', 7, t.sofa::double precision),
        ('SAPS-II', 'severity', 'sensitivity', 8, t.sapsii::double precision)
    ) v(variable, domain, role, sort_order, value)
),
summary AS (
    SELECT variable,
           domain,
           role,
           sort_order,
           count(*)::bigint AS cohort_n,
           count(value)::bigint AS available_n,
           (count(*) - count(value))::bigint AS missing_n,
           round(100.0 * (count(*) - count(value)) / NULLIF(count(*), 0), 2) AS missing_percent,
           percentile_cont(0.25) WITHIN GROUP (ORDER BY value)
             FILTER (WHERE value IS NOT NULL) AS p25,
           percentile_cont(0.50) WITHIN GROUP (ORDER BY value)
             FILTER (WHERE value IS NOT NULL) AS median,
           percentile_cont(0.75) WITHIN GROUP (ORDER BY value)
             FILTER (WHERE value IS NOT NULL) AS p75
    FROM long_values
    GROUP BY variable, domain, role, sort_order
)
SELECT variable,
       domain,
       role,
       cohort_n,
       available_n,
       missing_n,
       missing_percent,
       round(p25::numeric, 2) AS p25,
       round(median::numeric, 2) AS median,
       round(p75::numeric, 2) AS p75,
       round((p75 - p25)::numeric, 2) AS iqr
FROM summary
ORDER BY sort_order;
