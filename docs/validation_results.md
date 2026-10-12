# v1.0.1 public-release validation results

## Scope

These checks were run on the finalized source-only v1.0.1 public tree without connecting to MIMIC-IV, eICU-CRD, or any other database. No clinical SQL, MICE, LPA, transfer training, outcome model, or manuscript-result calculation was executed. The v1.0.1 changes are documentation and release-metadata corrections only; the R, SQL, and Python source files remain byte-identical to v1.0.0.

## Validation evidence

The v1.0.1 packaging/privacy static check was rerun on a clean 59-file staging tree. R syntax parsing and the two wholly synthetic component tests were not rerun for this documentation-only correction; the PASS results below are retained from v1.0.0 because every included R, SQL, and Python source file is byte-identical to that preserved release.

| Check | Result | Evidence |
|---|---|---|
| R syntax parsing | PASS (v1.0.0 evidence retained) | All 19 included R files under `R/` and `tests/` parsed with R 4.6.0 during v1.0.0 finalization on 2026-10-12; source hashes are unchanged in v1.0.1 |
| Python syntax parsing | PASS | Every included Python file parsed through Python `ast` in the rerun v1.0.1 static release check |
| Custom EM synthetic self-test | PASS (v1.0.0 evidence retained) | `STAGE8A_LPA_EM_SELF_TEST_OK`; source hashes are unchanged in v1.0.1 |
| Transfer probability-ensemble synthetic test | PASS (v1.0.0 evidence retained) | `SYNTHETIC_PROBABILITY_ENSEMBLE_TEST_OK`; source hashes are unchanged in v1.0.1 |
| Local SQL/source dependency presence | PASS | All declared local SQL dependencies used by R scripts were present in the clean v1.0.1 staging tree |
| Required release-document presence | PASS | README, LICENSE, workflow, mapping, environment, notices, citation, configuration, and tests were present |
| Forbidden file-type scan | PASS | No clinical-data, fitted-object, manuscript, figure, font, log, or Git-history file type was found in the clean 59-file v1.0.1 staging tree |
| Content safety scan | PASS | No private absolute path, credential-bearing URI, private key, common token pattern, or non-placeholder secret assignment was found |
| Code-to-results CSV import | PASS | Both public CSV mappings imported successfully with the bundled spreadsheet parser |

During v1.0.0 finalization, R startup emitted warnings that `C.UTF-8` locale components could not be set. This matches the recorded environment limitation; the effective recorded locale was `C`. The v1.0.1 clean-staging static release check reported `STATIC RELEASE CHECK: PASS (59 files inspected)`.

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

