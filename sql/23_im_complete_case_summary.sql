-- Stage 7A aggregate-only joint-completeness QC.
-- SOFA and SAPS-II extensions are evaluated separately by design.
-- No table is created and no identifier or patient-level row is returned.
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
specification_flags AS (
    SELECT v.specification,
           v.required_variables,
           v.sort_order,
           v.is_complete
    FROM target t
    CROSS JOIN LATERAL (VALUES
        (
            'inflammatory_domain',
            'WBC|platelet',
            1,
            t.wbc IS NOT NULL AND t.platelet_count IS NOT NULL
        ),
        (
            'metabolic_domain',
            'lactate|glucose|creatinine|bicarbonate',
            2,
            t.lactate IS NOT NULL
              AND t.glucose IS NOT NULL
              AND t.creatinine IS NOT NULL
              AND t.bicarbonate IS NOT NULL
        ),
        (
            'primary_all_6',
            'WBC|platelet|lactate|glucose|creatinine|bicarbonate',
            3,
            t.wbc IS NOT NULL
              AND t.platelet_count IS NOT NULL
              AND t.lactate IS NOT NULL
              AND t.glucose IS NOT NULL
              AND t.creatinine IS NOT NULL
              AND t.bicarbonate IS NOT NULL
        ),
        (
            'primary_plus_SOFA',
            'WBC|platelet|lactate|glucose|creatinine|bicarbonate|SOFA',
            4,
            t.wbc IS NOT NULL
              AND t.platelet_count IS NOT NULL
              AND t.lactate IS NOT NULL
              AND t.glucose IS NOT NULL
              AND t.creatinine IS NOT NULL
              AND t.bicarbonate IS NOT NULL
              AND t.sofa IS NOT NULL
        ),
        (
            'primary_plus_SAPS-II',
            'WBC|platelet|lactate|glucose|creatinine|bicarbonate|SAPS-II',
            5,
            t.wbc IS NOT NULL
              AND t.platelet_count IS NOT NULL
              AND t.lactate IS NOT NULL
              AND t.glucose IS NOT NULL
              AND t.creatinine IS NOT NULL
              AND t.bicarbonate IS NOT NULL
              AND t.sapsii IS NOT NULL
        )
    ) v(specification, required_variables, sort_order, is_complete)
),
summary AS (
    SELECT specification,
           required_variables,
           sort_order,
           count(*)::bigint AS cohort_n,
           count(*) FILTER (WHERE is_complete)::bigint AS complete_case_n
    FROM specification_flags
    GROUP BY specification, required_variables, sort_order
)
SELECT specification,
       required_variables,
       cohort_n,
       complete_case_n,
       (cohort_n - complete_case_n)::bigint AS incomplete_case_n,
       round(100.0 * complete_case_n / NULLIF(cohort_n, 0), 2) AS complete_case_percent
FROM summary
ORDER BY sort_order;
