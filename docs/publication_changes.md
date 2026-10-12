# Publication-finalization changes

## v1.0.1 documentation correction

This documentation-only release corrects the stated clinical dataset version to MIMIC-IV v3.1, as confirmed by the authors. No analytical source code was changed. The previously released v1.0.0 is preserved. The historical build versions of the 11 upstream `mimiciv_derived` relations remain unresolved.

The correction is limited to release documentation, citation/version metadata, the public configuration example, and descriptive source-to-result metadata. It does not modify SQL, R, or Python source, any statistical method, or any scientific result. MIMIC-IV v3.1 is the clinical database version and is not presented as a MIMIC Code tag or commit.

## Scope

The v1.0.0 public tree was created from the audited 58-file `Code_Release_v1` candidate. Research SQL, R, Python analysis logic, configuration defaults already sanitized in the candidate, and the frozen code-to-result mapping were not changed during publication finalization.

## Public-release changes

| Public file | Candidate source | Publication change |
|---|---|---|
| `LICENSE` | New, author-authorized release metadata | Added the standard MIT License with `Copyright (c) 2026 Wang Tao and Li Yanan` |
| `README.md` | `Code_Release_v1/README.md` | Replaced candidate/placeholder language with confirmed copyright, contributors, MIT boundary, repository/version metadata, public data-access requirements, output map, validation limits, and the final upstream-version disclosure |
| `CITATION.cff` | `Code_Release_v1/CITATION.cff` | Replaced placeholders with the confirmed software contributors, v1.0.0, release date, license, and repository/release URLs; no DOI, ORCID, acceptance claim, or manuscript author list was added |
| `THIRD_PARTY_NOTICES.md` | `Code_Release_v1/THIRD_PARTY_NOTICES.md` | Clarified the MIT boundary and non-vendoring finding; corrected the earlier `pypdf` narrative because the released Python source does not import or call it |
| `environment/DEPENDENCIES.md` | `Code_Release_v1/environment/DEPENDENCIES.md` | Corrected the Python dependency narrative and replaced the internal license-review reference with public `LICENSE`/notice references |
| `docs/upstream_schema_provenance.md` | `Code_Release_v1/docs/upstream_schema_provenance.md` | Replaced unresolved author-review placeholders with the final bounded negative finding; all 11 relations are explicitly `HISTORICAL_VERSION_UNRESOLVED` |
| `docs/analysis_workflow.md` | `Code_Release_v1/docs/analysis_workflow.md` | Changed candidate wording to v1.0.0 public-release wording only |
| `data/README.md` | `Code_Release_v1/data/README.md` | Changed candidate wording to public-release wording only |
| `tests/static_release_checks.py` | `Code_Release_v1/tests/static_release_checks.py` | Updated the documentation requirement from internal `LICENSE_REVIEW.md` to public `LICENSE`; no clinical or statistical analysis behavior was added |
| `python/stage12a_build_figures_tables.py` | `Code_Release_v1/python/stage12a_build_figures_tables.py` | Replaced a superseded candidate-status comment with a `CITATION.cff` pointer; executable logic is unchanged |
| `docs/validation_results.md` | `Code_Release_v1/docs/validation_results.md` | Updated after final non-clinical validation to describe the actual public release and final file count |
| `LICENSE_REVIEW.md` | `Code_Release_v1/LICENSE_REVIEW.md` | Excluded from the public tree because it is an internal pre-authorization audit containing superseded placeholders; the original remains unchanged in the project |

No manuscript, frozen scientific result, patient-derived object, patient-level record, credential, raw log, review correspondence, font, or original Git history is included.
