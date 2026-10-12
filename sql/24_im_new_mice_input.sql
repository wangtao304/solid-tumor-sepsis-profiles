-- Stage 7B controlled MICE input.
-- Only the six primary variables and eight frozen auxiliaries are streamed.
-- stay_id is used only for deterministic ordering and is not selected.
BEGIN TRANSACTION READ ONLY;
SET LOCAL statement_timeout = '5min';

SELECT
  wbc,
  platelet_count,
  lactate,
  glucose,
  creatinine,
  bicarbonate,
  age,
  sex,
  cancer_phenotype,
  charlson_comorbidity_index,
  sofa,
  sapsii,
  vasopressor_0_24,
  mechanical_ventilation_0_24
FROM mimic_cancer_sepsis.immune_metabolic_lca_base
WHERE cancer_phenotype IN ('metastatic solid', 'nonmetastatic solid')
ORDER BY stay_id;

COMMIT;
