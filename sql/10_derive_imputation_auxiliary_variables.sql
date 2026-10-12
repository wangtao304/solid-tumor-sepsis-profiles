-- Stage 4C: database-only derivation of prespecified MICE auxiliaries.
-- No patient-level rows are emitted by this script.
BEGIN;

DO $$
DECLARE
  base_n bigint;
  distinct_stay_n bigint;
  existing_aux_n integer;
BEGIN
  SELECT count(*), count(DISTINCT stay_id)
  INTO base_n, distinct_stay_n
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base;

  IF base_n <> 5904 OR distinct_stay_n <> 5904 THEN
    RAISE EXCEPTION 'Stage 4C precondition failed: expected 5904 unique stays, found % rows and % unique stays',
      base_n, distinct_stay_n;
  END IF;

  SELECT count(*)
  INTO existing_aux_n
  FROM information_schema.columns
  WHERE table_schema = 'mimic_cancer_sepsis'
    AND table_name = 'immune_metabolic_lca_base'
    AND column_name IN ('vasopressor_0_24', 'mechanical_ventilation_0_24');

  IF existing_aux_n <> 0 THEN
    RAISE EXCEPTION 'Stage 4C precondition failed: one or both auxiliary columns already exist';
  END IF;
END $$;

ALTER TABLE mimic_cancer_sepsis.immune_metabolic_lca_base
  ADD COLUMN vasopressor_0_24 boolean,
  ADD COLUMN mechanical_ventilation_0_24 boolean;

WITH target_window AS MATERIALIZED (
  SELECT
    b.stay_id,
    i.icu_intime,
    i.icu_intime + interval '24 hours' AS window_end
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base b
  JOIN mimiciv_derived.icustay_detail i USING (stay_id)
),
vasopressor_positive AS (
  SELECT DISTINCT t.stay_id
  FROM target_window t
  JOIN mimiciv_derived.vasoactive_agent v
    ON v.stay_id = t.stay_id
   AND v.starttime < t.window_end
   AND v.endtime > t.icu_intime
  WHERE coalesce(v.dopamine, 0) > 0
     OR coalesce(v.epinephrine, 0) > 0
     OR coalesce(v.norepinephrine, 0) > 0
     OR coalesce(v.phenylephrine, 0) > 0
     OR coalesce(v.vasopressin, 0) > 0
),
mechanical_ventilation_positive AS (
  SELECT DISTINCT t.stay_id
  FROM target_window t
  JOIN mimiciv_derived.ventilation mv
    ON mv.stay_id = t.stay_id
   AND mv.starttime < t.window_end
   AND mv.endtime > t.icu_intime
  WHERE mv.ventilation_status = 'InvasiveVent'
),
derived_auxiliaries AS (
  SELECT
    t.stay_id,
    (v.stay_id IS NOT NULL) AS vasopressor_0_24,
    (mv.stay_id IS NOT NULL) AS mechanical_ventilation_0_24
  FROM target_window t
  LEFT JOIN vasopressor_positive v USING (stay_id)
  LEFT JOIN mechanical_ventilation_positive mv USING (stay_id)
)
UPDATE mimic_cancer_sepsis.immune_metabolic_lca_base b
SET
  vasopressor_0_24 = d.vasopressor_0_24,
  mechanical_ventilation_0_24 = d.mechanical_ventilation_0_24
FROM derived_auxiliaries d
WHERE b.stay_id = d.stay_id;

ALTER TABLE mimic_cancer_sepsis.immune_metabolic_lca_base
  ALTER COLUMN vasopressor_0_24 SET NOT NULL,
  ALTER COLUMN mechanical_ventilation_0_24 SET NOT NULL;

COMMENT ON COLUMN mimic_cancer_sepsis.immune_metabolic_lca_base.vasopressor_0_24 IS
  'Stage 4C: positive-dose dopamine, epinephrine, norepinephrine, phenylephrine, or vasopressin interval overlapping ICU [0,24h); database-only MICE auxiliary.';

COMMENT ON COLUMN mimic_cancer_sepsis.immune_metabolic_lca_base.mechanical_ventilation_0_24 IS
  'Stage 4C: InvasiveVent interval overlapping ICU [0,24h); database-only MICE auxiliary.';

DO $$
DECLARE
  invalid_n bigint;
BEGIN
  SELECT count(*)
  INTO invalid_n
  FROM mimic_cancer_sepsis.immune_metabolic_lca_base
  WHERE vasopressor_0_24 IS NULL
     OR mechanical_ventilation_0_24 IS NULL;

  IF invalid_n <> 0 THEN
    RAISE EXCEPTION 'Stage 4C validation failed: % rows have a missing auxiliary value', invalid_n;
  END IF;
END $$;

COMMIT;
