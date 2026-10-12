#!/usr/bin/env Rscript

# Runs only the synthetic recovery check embedded in the custom LPA script.
# No clinical data are read when STAGE8A_SELF_TEST_ONLY=1.
Sys.setenv(STAGE8A_SELF_TEST_ONLY = "1")
source(file.path("R", "25_run_inflammatory_metabolic_lpa.R"), chdir = FALSE)
