# v1.0.0 public-release validation results

## Scope

These checks were run on the finalized source-only public tree without connecting to MIMIC-IV, eICU-CRD, or any other database. No clinical SQL, MICE, LPA, transfer training, outcome model, or manuscript-result calculation was executed.

## Executed checks

| Check | Result | Evidence |
|---|---|---|
| R syntax parsing | PASS | All 19 included R files under `R/` and `tests/` parsed with R 4.6.0 during finalization on 2026-10-12 |
| Python syntax parsing | PASS | Every included Python file parsed through Python `ast` in the final static release check |
| Custom EM synthetic self-test | PASS | `STAGE8A_LPA_EM_SELF_TEST_OK` |
| Transfer probability-ensemble synthetic test | PASS | `SYNTHETIC_PROBABILITY_ENSEMBLE_TEST_OK` |
| Local SQL/source dependency presence | PASS | All declared local SQL dependencies used by R scripts were present |
| Required release-document presence | PASS | README, LICENSE, workflow, mapping, environment, notices, citation, configuration, and tests were present |
| Forbidden file-type scan | PASS | No clinical-data, fitted-object, manuscript, figure, font, log, or Git-history file type was found in the 59-file public tree |
| Content safety scan | PASS | No private absolute path, credential-bearing URI, private key, common token pattern, or non-placeholder secret assignment was found |
| Code-to-results CSV import | PASS | Both public CSV mappings imported successfully with the bundled spreadsheet parser |

R startup emitted warnings that `C.UTF-8` locale components could not be set. This matches the recorded environment limitation; the effective recorded locale was `C`. The parse and synthetic tests still passed. The final static release check reported `STATIC RELEASE CHECK: PASS (59 files inspected)`.

## Checks deliberately not performed

- No database connection or SQL execution.
- No clinical-data access or patient-level extract.
- No read of any original or regenerated patient-derived RDS object.
- No MICE, LPA, profile selection, transfer-classifier fitting, survival model, logistic model, bootstrap, or landmark model execution.
- No recomputation or numerical comparison of frozen manuscript results.
- No package installation or upgrade.
- No remote repository creation, push, upload, or public release.
- No assertion that the exact historical Python/PostgreSQL/Poppler/transitive R environment has been reconstructed.

## Interpretation

The checks establish source syntax, synthetic behavior for the two custom components, local dependency presence, and packaging/privacy safeguards. They do not independently reproduce or validate the clinical findings.

