-- Stage 4C aggregate-only QC for the two MICE auxiliary variables.
WITH long_qc AS (
  SELECT
    v.variable,
    v.variable_order,
    v.value
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base b
  CROSS JOIN LATERAL (
    VALUES
      ('vasopressor_0_24'::text, 1, b.vasopressor_0_24),
      ('mechanical_ventilation_0_24'::text, 2, b.mechanical_ventilation_0_24)
  ) v(variable, variable_order, value)
)
SELECT
  variable,
  count(value) AS available_n,
  round(100.0 * count(*) FILTER (WHERE value IS NULL) / count(*), 2) AS missing_percent,
  count(*) FILTER (WHERE value IS TRUE) AS positive_n,
  round(
    100.0 * count(*) FILTER (WHERE value IS TRUE) / nullif(count(value), 0),
    2
  ) AS positive_percent
FROM long_qc
GROUP BY variable, variable_order
ORDER BY variable_order;
