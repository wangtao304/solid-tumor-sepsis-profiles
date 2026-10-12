-- Stage 5 MICE input query.
-- stay_id is used only to make row order deterministic; it is not selected or transferred.
SELECT
  neutrophil_abs,
  lymphocyte_abs,
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
ORDER BY stay_id;
