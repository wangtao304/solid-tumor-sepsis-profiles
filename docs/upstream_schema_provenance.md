# Upstream derived-schema provenance

## Status

The released SQL was run against MIMIC-IV v3.1 and an externally materialized `mimiciv_derived` schema. A bounded local review of the available project records, referenced build materials, and version metadata did not recover the historical MIMIC Code checkout, an immutable commit or tag, the materialization date, or a module-level build manifest.

The historical build version for every consumed relation is therefore `HISTORICAL_VERSION_UNRESOLVED`. MIMIC-IV v3.1 is a database release and is not evidence of any MIMIC Code tag or commit, including MIMIC Code v3.0.0. No current or later upstream revision is presented as the historically executed source. Historical local modifications cannot be assessed, and byte-for-byte or one-command reconstruction of the historical derived schema is not claimed.

## Relations referenced by the released source

| Relation | Historical build status |
|---|---|
| `mimiciv_derived.bg` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.blood_differential` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.charlson` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.chemistry` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.complete_blood_count` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.first_day_sofa` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.icustay_detail` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.sapsii` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.sepsis3` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.vasoactive_agent` | `HISTORICAL_VERSION_UNRESOLVED` |
| `mimiciv_derived.ventilation` | `HISTORICAL_VERSION_UNRESOLVED` |

The list was obtained by static inspection of the included R and SQL source. It is a dependency inventory, not a substitute for the upstream SQL definitions. The upstream definition SQL is not redistributed. Users who construct a compatible derived schema must independently select, verify, and document the upstream revision they use; that revision must not be described as the historical study build unless supported by new direct evidence.

