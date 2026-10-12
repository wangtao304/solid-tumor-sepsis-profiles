-- Stage 12L-1: repeated ICU stay aggregate audit after cohort gate passes.
-- The query returns aggregate metrics only; no identifier is selected.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '5min';

WITH eligible AS (
    SELECT subject_id
    FROM mimic_cancer_sepsis.immune_metabolic_lca_base
    WHERE cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
),
per_subject AS (
    SELECT subject_id, count(*)::bigint AS eligible_stays
    FROM eligible
    GROUP BY subject_id
),
summary AS (
    SELECT
        (SELECT count(*)::bigint FROM eligible) AS total_stays,
        count(*)::bigint AS unique_subjects,
        count(*) FILTER (WHERE eligible_stays = 1)::bigint AS subjects_one,
        count(*) FILTER (WHERE eligible_stays > 1)::bigint AS subjects_repeat,
        coalesce(sum(eligible_stays) FILTER (WHERE eligible_stays > 1), 0)::bigint AS stays_repeat,
        max(eligible_stays)::bigint AS max_stays,
        count(*) FILTER (WHERE eligible_stays = 2)::bigint AS subjects_two,
        count(*) FILTER (WHERE eligible_stays = 3)::bigint AS subjects_three,
        count(*) FILTER (WHERE eligible_stays >= 4)::bigint AS subjects_four_plus
    FROM per_subject
)
SELECT metric, value, unit, denominator
FROM summary
CROSS JOIN LATERAL (
    VALUES
        ('total_eligible_icu_stays', total_stays::numeric, 'count', NULL::numeric),
        ('unique_subjects', unique_subjects::numeric, 'count', NULL::numeric),
        ('subjects_with_exactly_1_eligible_stay', subjects_one::numeric, 'count', unique_subjects::numeric),
        ('subjects_with_more_than_1_eligible_stay', subjects_repeat::numeric, 'count', unique_subjects::numeric),
        ('percentage_of_subjects_with_more_than_1_eligible_stay', round(100.0 * subjects_repeat / NULLIF(unique_subjects, 0), 6), 'percent', unique_subjects::numeric),
        ('number_of_stays_from_repeat_contributors', stays_repeat::numeric, 'count', total_stays::numeric),
        ('percentage_of_stays_from_repeat_contributors', round(100.0 * stays_repeat / NULLIF(total_stays, 0), 6), 'percent', total_stays::numeric),
        ('excess_repeated_stays', (total_stays - unique_subjects)::numeric, 'count', total_stays::numeric),
        ('percentage_of_excess_repeated_stays', round(100.0 * (total_stays - unique_subjects) / NULLIF(total_stays, 0), 6), 'percent', total_stays::numeric),
        ('maximum_number_of_eligible_stays_per_subject', max_stays::numeric, 'count', NULL::numeric),
        ('subjects_with_1_eligible_stay', subjects_one::numeric, 'count', unique_subjects::numeric),
        ('subjects_with_2_eligible_stays', subjects_two::numeric, 'count', unique_subjects::numeric),
        ('subjects_with_3_eligible_stays', subjects_three::numeric, 'count', unique_subjects::numeric),
        ('subjects_with_4_or_more_eligible_stays', subjects_four_plus::numeric, 'count', unique_subjects::numeric)
) AS metrics(metric, value, unit, denominator);

COMMIT;
