# Protected data boundary

This directory is intentionally empty in the public source release. Authorized users may create `data/processed/` and `data/interim/` locally when running the workflow inside an approved clinical-data environment. Those directories may contain patient-level imputation, profile, and transfer objects and must never be committed, uploaded, or shared.

