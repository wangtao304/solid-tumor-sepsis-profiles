#!/usr/bin/env Rscript

# Reviewer Revision R1C: prespecified 24-hour-landmark time-varying analysis.
# Patient identifiers and clinical timestamps exist only in protected process
# memory. Every written result is aggregate. No fitted model object is saved.

options(stringsAsFactors = FALSE, warn = 1)

if (!identical(Sys.getenv("R1C_EXECUTE_TIME_VARYING"), "YES")) {
  cat("R1C execution gate is closed; no query or model was run.\n")
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages({
  library(mice)
  library(survival)
  library(rstpm2)
})

if (packageVersion("rstpm2") != package_version("1.7.1")) {
  stop("R1C was frozen for rstpm2 1.7.1; execution refused for another version.")
}

out_dir <- "outputs/revision"
mice_path <- "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
lpa_path <- "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
sql_path <- "sql/revision/R1B_landmark_boundary_input.sql"
r1b_population_path <- "outputs/revision/R1B_landmark_population_audit.csv"
r1b_boundary_path <- "outputs/revision/R1B_death_date_boundary_audit.csv"
r1b_script_path <- "R/revision/R1B_time_varying_analysis_ready.R"
r1b_cluster_path <- "outputs/revision/R1B_cluster_robust_results.csv"

output_paths <- file.path(out_dir, c(
  "R1C_time_varying_hr_points.csv",
  "R1C_time_varying_hr_curve.csv",
  "R1C_same_cohort_model_comparison.csv",
  "R1C_time_interaction_tests.csv",
  "R1C_date_precision_sensitivity.csv",
  "R1C_model_diagnostics.csv",
  "R1C_landmark_exclusion_breakdown.csv",
  "R1C_time_varying_hr.pdf",
  "R1C_time_varying_hr.tiff",
  "R1C_date_precision_sensitivity.pdf",
  "R1C_date_precision_sensitivity.tiff"
))
if (any(file.exists(output_paths)) &&
    !identical(Sys.getenv("R1C_REPLACE_OWN_OUTPUTS"), "YES")) {
  stop("An R1C output already exists. Set R1C_REPLACE_OWN_OUTPUTS=YES only to replace R1C-owned outputs.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_hashes <- c(
  "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds" =
    "770F7AEA660C2352F697BA37272E9A79971B231858281D3ADB1597B4FB9EE249",
  "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds" =
    "5AF61A8DC687BA6480B62B79A20F5945ECB3217FC611EFBC6D3A550FF8A4D366",
  "R/revision/R1B_time_varying_analysis_ready.R" =
    "9DECF00C0DADB5762D9BD87383DD5A9506AD13FBAD9EB014091279DA7E1C10EE",
  "sql/revision/R1B_landmark_boundary_input.sql" =
    "4B44AF0EA22C45059942FA926069A540FB33160072FD3B904C0C4AAB59FAA8C5",
  "outputs/revision/R1B_landmark_population_audit.csv" =
    "7465D0BA20B14292185D84101409F277DFA5159C4DB67F4BD3F538F5F2A11279",
  "outputs/revision/R1B_death_date_boundary_audit.csv" =
    "AAD8AD2EAE57708903908B1DE988C8353EC51E02E0B21260A5E854CEB20A73A3",
  "outputs/revision/R1B_cluster_robust_results.csv" =
    "D5F8FDE0BC6907AA876867F1DAF8ACD1E2E1D1117627F9B9426BB8484E389202"
)
actual_hashes <- unname(tools::md5sum(names(expected_hashes)))
# tools::md5sum is intentionally not used for the gate. PowerShell generated
# the SHA-256 freeze record; R verifies SHA-256 through certutil below.
sha256_file <- function(path) {
  lines <- system2("certutil", c("-hashfile", path, "SHA256"), stdout = TRUE, stderr = TRUE)
  hit <- grep("^[0-9A-Fa-f ]{64,}$", trimws(lines), value = TRUE)
  if (length(hit) != 1L) stop("Could not obtain SHA-256 for ", path)
  toupper(gsub(" ", "", trimws(hit)))
}
sha_gate <- vapply(names(expected_hashes), sha256_file, character(1))
if (!identical(unname(sha_gate), unname(expected_hashes))) {
  stop("A frozen R1C input SHA-256 does not match the pre-fit specification.")
}
rm(actual_hashes, sha_gate)

profile_variables <- c(
  "wbc", "platelet_count", "lactate", "glucose", "creatinine", "bicarbonate"
)
mi <- readRDS(mice_path)
lpa <- readRDS(lpa_path)
if (
  !inherits(mi, "mids") || mi$m != 10L || nrow(mi$data) != 4380L ||
  !identical(lpa$configuration$variables, profile_variables) ||
  !identical(lpa$outcome_variables_used, FALSE) ||
  !identical(lpa$survival_information_used, FALSE)
) stop("Frozen MICE/LPA safeguard check failed.")

# Reconstruct the frozen two-profile consensus; no discovery or relabeling.
votes <- matrix(NA_integer_, 4380L, 10L)
posterior <- matrix(0, 4380L, 2L)
for (i in seq_len(10L)) {
  key <- sprintf("imp%02d_k2", i)
  record <- lpa$fit_records[[key]]
  alignment <- lpa$alignments[[key]]
  if (is.null(record) || is.null(alignment) || length(alignment) != 2L) {
    stop("Malformed frozen LPA record: ", key)
  }
  votes[, i] <- match(record$model$classification, alignment)
  posterior <- posterior + record$model$posterior[, alignment, drop = FALSE] / 10
}
v1 <- rowSums(votes == 1L)
v2 <- rowSums(votes == 2L)
profile <- ifelse(v1 > v2, 1L, ifelse(v2 > v1, 2L, max.col(posterior)))
if (!identical(as.integer(tabulate(profile, 2L)), c(3197L, 1183L))) {
  stop("Frozen profile consensus could not be reproduced.")
}

# Read-only patient/time linkage. No identifier or timestamp is written.
psql_bin <- Sys.getenv("MIMIC_PSQL_BIN", unset = "")
pgpass_file <- Sys.getenv(
  "MIMIC_PGPASS_FILE",
  unset = ""
)
db_host <- Sys.getenv("MIMIC_DB_HOST", "")
db_port <- Sys.getenv("MIMIC_DB_PORT", "")
db_user <- Sys.getenv("MIMIC_DB_USER", "")
db_name <- Sys.getenv("MIMIC_DB_NAME", "")
if (!file.exists(psql_bin) || !file.exists(pgpass_file) ||
    any(!nzchar(c(db_host, db_port, db_user, db_name)))) {
  stop("Set all MIMIC database environment variables to authorized local values.")
}
Sys.setenv(PGPASSFILE = pgpass_file)
psql_args <- c(
  "-h", db_host,
  "-p", db_port,
  "-U", db_user,
  "-d", db_name,
  "-X", "-w", "-v", "ON_ERROR_STOP=1", "-q", "--csv", "-f", sql_path
)
query_text <- system2(psql_bin, psql_args, stdout = TRUE, stderr = "")
if (!is.null(attr(query_text, "status")) && attr(query_text, "status") != 0L) {
  stop("Read-only landmark-boundary query failed.")
}
db <- utils::read.csv(
  text = paste(query_text, collapse = "\n"), check.names = FALSE, na.strings = ""
)
rm(query_text)

numeric_equal <- function(x, y, tolerance = 1e-8) {
  x <- as.numeric(x); y <- as.numeric(y)
  identical(is.na(x), is.na(y)) &&
    all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}
local <- mi$data
link_checks <- c(
  age = numeric_equal(db$age, local$age),
  sex = identical(as.character(db$sex), as.character(local$sex)),
  cancer = identical(as.character(db$cancer_phenotype), as.character(local$cancer_phenotype)),
  charlson = numeric_equal(db$charlson_comorbidity_index, local$charlson_comorbidity_index),
  sofa = numeric_equal(db$sofa, local$sofa),
  sapsii = numeric_equal(db$sapsii, local$sapsii)
)
db$event_28d <- as.integer(db$event_28d)
if (
  nrow(db) != 4380L || length(unique(db$stay_id)) != 4380L ||
  length(unique(db$subject_id)) != 3772L || anyDuplicated(db$stay_id) ||
  anyNA(db[c("subject_id", "stay_id", "event_28d", "icu_intime")]) ||
  !all(link_checks) || sum(db$event_28d) != 1537L ||
  !all(db$event_28d %in% 0:1)
) stop("R1C landmark linkage/cohort gate failed.")

parse_time <- function(x) {
  x[x == ""] <- NA_character_
  as.POSIXct(x, tz = "UTC", tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%OS"))
}
icu_in <- parse_time(db$icu_intime)
icu_out <- parse_time(db$icu_outtime)
discharge <- parse_time(db$dischtime)
death_exact <- parse_time(db$deathtime)
death_date <- as.Date(db$death_date)
L <- icu_in + 24 * 3600
H <- icu_in + 28 * 24 * 3600
date_lower <- as.POSIXct(death_date, tz = "UTC")
date_upper <- date_lower + 24 * 3600

has_exact <- !is.na(death_exact)
has_date <- !is.na(death_date)
exact_in_followup <- has_exact & death_exact >= icu_in & death_exact <= H
exact_at_or_before_L <- db$event_28d == 1L & exact_in_followup & death_exact <= L
exact_after_L <- db$event_28d == 1L & exact_in_followup & death_exact > L
exact_date_conflict <- db$event_28d == 1L & has_exact & has_date &
  as.Date(death_exact, tz = "UTC") != death_date
invalid_exact_for_event <- db$event_28d == 1L & has_exact & !exact_in_followup
missing_event_time <- db$event_28d == 1L & !has_exact & !has_date
death_evidence_with_zero_flag <- db$event_28d == 0L & (
  exact_in_followup | (has_date & date_upper > icu_in & date_lower <= H)
)
icu_duration_conflict <- !is.na(icu_out) & icu_out < L
other_time_conflict <- exact_date_conflict | invalid_exact_for_event |
  missing_event_time | death_evidence_with_zero_flag | icu_duration_conflict

date_only <- db$event_28d == 1L & !has_exact & has_date & !other_time_conflict
date_before_or_at_L <- date_only & date_upper <= L
date_spans_L <- date_only & date_lower <= L & date_upper > L
post_L_observation <- (!is.na(icu_out) & icu_out > L) |
  (!is.na(discharge) & discharge > L)
date_L_resolved_by_observation <- date_spans_L & post_L_observation
landmark_uncertain <- date_spans_L & !post_L_observation
date_spans_H <- date_only & date_lower <= H & date_upper > H
date_definite_event_by_H <- date_only & date_upper <= H
date_definite_after_H <- date_only & date_lower > H
exclude_pre_or_at_L <- (exact_at_or_before_L | date_before_or_at_L) & !other_time_conflict
retain <- !exclude_pre_or_at_L & !landmark_uncertain & !other_time_conflict

effective_lower <- date_lower
effective_lower[date_L_resolved_by_observation] <- L[date_L_resolved_by_observation]
midpoint <- effective_lower +
  as.numeric(difftime(date_upper, effective_lower, units = "secs")) / 2
latest_legal <- date_upper - 1
icu_shortfall_minutes <- as.numeric(difftime(L, icu_out, units = "mins"))

base_covariates <- data.frame(
  subject_id = db$subject_id[retain],
  profile2 = as.integer(profile[retain] == 2L),
  age10 = (as.numeric(local$age[retain]) - mean(as.numeric(local$age[retain]))) / 10,
  sex = factor(local$sex[retain]),
  charlson = as.numeric(local$charlson_comorbidity_index[retain]) -
    mean(as.numeric(local$charlson_comorbidity_index[retain])),
  cancer_phenotype = factor(local$cancer_phenotype[retain]),
  sofa = as.numeric(local$sofa[retain]) - mean(as.numeric(local$sofa[retain])),
  check.names = FALSE
)
reference_means <- c(
  age_years = mean(as.numeric(local$age[retain])),
  charlson = mean(as.numeric(local$charlson_comorbidity_index[retain])),
  sofa = mean(as.numeric(local$sofa[retain]))
)

make_scenario <- function(scenario) {
  event <- integer(nrow(db))
  end_time <- H
  event[exact_after_L & retain] <- 1L
  end_time[exact_after_L & retain] <- death_exact[exact_after_L & retain]
  if (scenario == "midpoint") {
    idx <- date_only & retain & midpoint <= H
    event[idx] <- 1L
    end_time[idx] <- midpoint[idx]
  } else if (scenario == "earliest") {
    idx <- date_only & retain & effective_lower <= H
    event[idx] <- 1L
    end_time[idx] <- effective_lower[idx]
  } else if (scenario == "latest") {
    idx <- date_only & retain & date_upper <= H
    event[idx] <- 1L
    end_time[idx] <- latest_legal[idx]
  } else stop("Unknown date-time scenario.")
  time <- as.numeric(difftime(end_time[retain], L[retain], units = "days"))
  ans <- base_covariates
  ans$time <- time
  ans$event <- event[retain]
  if (anyNA(ans) || any(ans$time <= 0) || any(ans$time > 27 + 1e-10)) {
    stop("Scenario ", scenario, " produced missing or invalid post-landmark time.")
  }
  ans
}

scenario_data <- setNames(lapply(c("midpoint", "earliest", "latest"), make_scenario),
                          c("midpoint", "earliest", "latest"))

scenario_gate <- vapply(scenario_data, nrow, integer(1))
subject_gate <- vapply(scenario_data, function(x) length(unique(x$subject_id)), integer(1))
event_gate <- vapply(scenario_data, function(x) sum(x$event), integer(1))
profile1_gate <- vapply(scenario_data, function(x) sum(x$profile2 == 0L), integer(1))
profile2_gate <- vapply(scenario_data, function(x) sum(x$profile2 == 1L), integer(1))
if (
  !identical(unname(scenario_gate), rep(4340L, 3L)) ||
  !identical(unname(subject_gate), rep(3742L, 3L)) ||
  !identical(unname(event_gate), c(1500L, 1505L, 1490L)) ||
  !identical(unname(profile1_gate), rep(3181L, 3L)) ||
  !identical(unname(profile2_gate), rep(1159L, 3L)) ||
  !identical(scenario_data$midpoint[c("subject_id", "profile2", "age10", "sex", "charlson", "cancer_phenotype", "sofa")],
             scenario_data$earliest[c("subject_id", "profile2", "age10", "sex", "charlson", "cancer_phenotype", "sofa")]) ||
  !identical(scenario_data$midpoint[c("subject_id", "profile2", "age10", "sex", "charlson", "cancer_phenotype", "sofa")],
             scenario_data$latest[c("subject_id", "profile2", "age10", "sex", "charlson", "cancer_phenotype", "sofa")])
) stop("The three frozen scenario populations or event counts did not reproduce R1B.")

# Reproduce and compare the aggregate R1B audits before fitting.
metric <- function(name, value, unit, status, note) {
  data.frame(metric = name, value = value, unit = unit, status = status,
             note = note, stringsAsFactors = FALSE)
}
population <- do.call(rbind, list(
  metric("original_candidate_stays", nrow(db), "stays", "PASS", "Frozen solid-tumor MIMIC cohort."),
  metric("original_unique_subjects", length(unique(db$subject_id)), "subjects", "PASS", "Identifiers used only in protected memory."),
  metric("original_profile_1_stays", sum(profile == 1L), "stays", "PASS", "Frozen consensus assignment."),
  metric("original_profile_2_stays", sum(profile == 2L), "stays", "PASS", "Frozen consensus assignment."),
  metric("original_28_day_death_records", sum(db$event_28d), "records", "PASS", "Frozen outcome flag before landmark reclassification."),
  metric("exact_timestamp_death_records", sum(db$event_28d == 1L & exact_in_followup & !other_time_conflict), "records", "AUDITED", "Valid exact timestamp has priority."),
  metric("date_only_death_records", sum(date_only), "records", "AUDITED", "Death interval is [death date 00:00, next date 00:00)."),
  metric("administratively_censored_records", sum(db$event_28d == 0L & !other_time_conflict), "records", "AUDITED", "No frozen 28-day death flag or conflicting death evidence."),
  metric("death_at_or_before_landmark", sum(exclude_pre_or_at_L), "records", "EXCLUDED", "No post-landmark risk time; exact and definite date-only records."),
  metric("follow_up_termination_at_or_before_landmark", 0L, "records", "AUDITED", "Primary all-cause follow-up is administratively censored at ICU day 28, not ICU discharge."),
  metric("landmark_vital_status_uncertain", sum(landmark_uncertain), "records", "EXCLUDED", "Date-only interval spans L and no credible post-L observation resolves survival."),
  metric("other_time_conflict", sum(other_time_conflict), "records", "EXCLUDED", "Exact/date disagreement, invalid exact event time, missing event time, zero-flag death evidence, or ICU-duration conflict."),
  metric("landmark_analysis_stays", sum(retain), "stays", "READY", "Post-L risk set under frozen exclusion rules."),
  metric("landmark_analysis_subjects", length(unique(db$subject_id[retain])), "subjects", "READY", "Unique subjects in post-L risk set."),
  metric("landmark_profile_1_stays", sum(retain & profile == 1L), "stays", "READY", "Frozen Profile 1."),
  metric("landmark_profile_2_stays", sum(retain & profile == 2L), "stays", "READY", "Frozen Profile 2."),
  metric("landmark_midpoint_working_events", sum(scenario_data$midpoint$event), "events", "WORKING_CONVENTION", "Exact timestamps plus midpoint convention for eligible date-only intervals."),
  metric("landmark_earliest_scenario_events", sum(scenario_data$earliest$event), "events", "SENSITIVITY_SCENARIO", "Earliest legal time for date-only intervals."),
  metric("landmark_latest_scenario_events", sum(scenario_data$latest$event), "events", "SENSITIVITY_SCENARIO", "Latest legal time for date-only intervals."),
  metric("maximum_landmark_follow_up_days", max(scenario_data$midpoint$time), "days", "PASS", "H remains ICU admission +28 days; maximum post-L follow-up is 27 days."),
  metric("minimum_landmark_follow_up_days", min(scenario_data$midpoint$time), "days", "PASS", "No zero or negative post-landmark risk time."),
  metric("covariate_row_link_checks_passed", sum(link_checks), "checks", "PASS", paste(names(link_checks), collapse = "; "))
))
boundary <- do.call(rbind, list(
  metric("exact_timestamp_valid_events", sum(db$event_28d == 1L & exact_in_followup & !other_time_conflict), "records", "AUDITED", "Exact timestamp within ICU admission through H."),
  metric("date_only_valid_records", sum(date_only), "records", "AUDITED", "No usable exact timestamp and no time conflict."),
  metric("exact_event_at_or_before_L", sum(exact_at_or_before_L & !other_time_conflict), "records", "EXCLUDED", "Exact death at or before ICU admission +24 h."),
  metric("date_only_definitely_at_or_before_L", sum(date_before_or_at_L), "records", "EXCLUDED", "Entire date interval ends by L."),
  metric("date_only_spans_L_resolved_by_post_L_observation", sum(date_L_resolved_by_observation), "records", "INCLUDED", "ICU or hospital timestamp documents observation after L; post-L legal interval is used."),
  metric("date_only_spans_L_unresolved", sum(landmark_uncertain), "records", "LANDMARK_VITAL_STATUS_UNCERTAIN", "Excluded from landmark main analysis."),
  metric("date_only_definite_event_by_H", sum(date_definite_event_by_H & retain), "records", "EVENT_IN_ALL_SCENARIOS", "Entire eligible date interval ends by H."),
  metric("date_only_spans_H", sum(date_spans_H & retain), "records", "DAY28_EVENT_STATUS_UNCERTAIN", "Event status differs between earliest and latest legal-time scenarios."),
  metric("date_only_definitely_after_H", sum(date_definite_after_H & retain), "records", "CENSORED_IN_ALL_SCENARIOS", "Entire date interval begins after H."),
  metric("exact_date_record_conflict", sum(exact_date_conflict), "records", "OTHER_TIME_CONFLICT", "Exact timestamp date and date-only record disagree."),
  metric("invalid_exact_time_for_positive_flag", sum(invalid_exact_for_event), "records", "OTHER_TIME_CONFLICT", "Exact timestamp lies outside ICU admission through H for a positive frozen event flag."),
  metric("positive_flag_missing_any_death_time", sum(missing_event_time), "records", "OTHER_TIME_CONFLICT", "Positive frozen event flag without exact or date-only time."),
  metric("death_evidence_with_zero_flag", sum(death_evidence_with_zero_flag), "records", "OTHER_TIME_CONFLICT", "Death evidence overlaps the 28-day window despite a zero frozen event flag."),
  metric("icu_outtime_before_L", sum(icu_duration_conflict), "records", "OTHER_TIME_CONFLICT", "Contradicts the >=24 h cohort criterion."),
  metric("icu_outtime_shortfall_at_most_1_minute", sum(icu_duration_conflict & icu_shortfall_minutes <= 1, na.rm = TRUE), "records", "OTHER_TIME_CONFLICT_DETAIL", "Aggregate timing detail; no stay identifier exported."),
  metric("icu_outtime_shortfall_over_1_to_60_minutes", sum(icu_duration_conflict & icu_shortfall_minutes > 1 & icu_shortfall_minutes <= 60, na.rm = TRUE), "records", "OTHER_TIME_CONFLICT_DETAIL", "Aggregate timing detail; no stay identifier exported."),
  metric("icu_outtime_shortfall_over_60_minutes", sum(icu_duration_conflict & icu_shortfall_minutes > 60, na.rm = TRUE), "records", "OTHER_TIME_CONFLICT_DETAIL", "Aggregate timing detail; no stay identifier exported."),
  metric("midpoint_minus_earliest_event_count", sum(scenario_data$midpoint$event) - sum(scenario_data$earliest$event), "events", "SCENARIO_DIFFERENCE", "Midpoint working convention minus earliest legal-time scenario."),
  metric("latest_minus_midpoint_event_count", sum(scenario_data$latest$event) - sum(scenario_data$midpoint$event), "events", "SCENARIO_DIFFERENCE", "Latest legal-time scenario minus midpoint working convention."),
  metric("earliest_minus_latest_event_count", sum(scenario_data$earliest$event) - sum(scenario_data$latest$event), "events", "SCENARIO_RANGE", "Maximum aggregate event-count difference induced by date precision."),
  metric("raw_timestamp_or_identifier_fields_exported", 0L, "fields", "PASS", "All identifiers and raw date/time fields were discarded before writing aggregate files.")
))
saved_population <- utils::read.csv(r1b_population_path, check.names = FALSE)
saved_boundary <- utils::read.csv(r1b_boundary_path, check.names = FALSE)
population_comparison <- all.equal(saved_population, population, check.attributes = FALSE)
boundary_comparison <- all.equal(saved_boundary, boundary, check.attributes = FALSE)
if (!isTRUE(population_comparison) || !isTRUE(boundary_comparison)) {
  cat("Population aggregate comparison: ", paste(population_comparison, collapse = " | "), "\n")
  cat("Boundary aggregate comparison: ", paste(boundary_comparison, collapse = " | "), "\n")
  stop("Recomputed R1C landmark audit differs from the frozen R1B aggregate audit.")
}

exclusion_breakdown <- do.call(rbind, list(
  metric("original_eligible_stays", 4380L, "stays", "SOURCE_COHORT", "Frozen solid-tumor discovery cohort with database LOS >=24 h criterion."),
  metric("exact_death_at_or_before_24h_landmark", sum(exact_at_or_before_L & !other_time_conflict), "stays", "EXCLUDED", "Exact deaths at or before L have no positive post-landmark risk time."),
  metric("date_only_death_definitely_at_or_before_landmark", sum(date_before_or_at_L), "stays", "EXCLUDED", "Date interval ends by L."),
  metric("landmark_vital_status_uncertain", sum(landmark_uncertain), "stays", "EXCLUDED", "No unresolved date-only interval spanned L in this cohort."),
  metric("other_time_conflict_total", sum(other_time_conflict), "stays", "EXCLUDED", "Mutually exclusive from landmark-death exclusions."),
  metric("exact_date_disagreement", sum(exact_date_conflict), "stays", "OTHER_TIME_CONFLICT", "Exact timestamp and death-date record disagree."),
  metric("invalid_exact_time_for_positive_event", sum(invalid_exact_for_event), "stays", "OTHER_TIME_CONFLICT", "Exact timestamp outside the eligible follow-up for a positive event flag."),
  metric("positive_event_missing_death_time", sum(missing_event_time), "stays", "OTHER_TIME_CONFLICT", "Positive event flag without exact or date-only time."),
  metric("death_evidence_with_zero_event_flag", sum(death_evidence_with_zero_flag), "stays", "OTHER_TIME_CONFLICT", "Death evidence overlaps the horizon but the frozen event flag is zero."),
  metric("icu_outtime_before_24h_landmark", sum(icu_duration_conflict), "stays", "OTHER_TIME_CONFLICT", "Timestamp-derived outtime is earlier than L despite database LOS eligibility; retained source values are not overwritten."),
  metric("icu_outtime_shortfall_at_most_1_minute", sum(icu_duration_conflict & icu_shortfall_minutes <= 1, na.rm = TRUE), "stays", "CONFLICT_DETAIL", "Near-boundary timestamp or duration-rounding discordance."),
  metric("icu_outtime_shortfall_over_1_to_60_minutes", sum(icu_duration_conflict & icu_shortfall_minutes > 1 & icu_shortfall_minutes <= 60, na.rm = TRUE), "stays", "CONFLICT_DETAIL", "Sub-hour discordance between timestamp-derived duration and database cohort eligibility."),
  metric("icu_outtime_shortfall_over_60_minutes", sum(icu_duration_conflict & icu_shortfall_minutes > 60, na.rm = TRUE), "stays", "CONFLICT_DETAIL", "No such record was observed."),
  metric("landmark_analysis_stays", sum(retain), "stays", "INCLUDED", "Eligible after conservative exclusions."),
  metric("day28_date_interval_boundary_uncertain", sum(date_spans_H & retain), "stays", "SCENARIO_SENSITIVE_NOT_EXCLUDED", "Changes event status between earliest and latest date-precision scenarios; not a cohort exclusion.")
))
if (4380L != sum(retain) + sum(exclude_pre_or_at_L) + sum(landmark_uncertain) + sum(other_time_conflict)) {
  stop("Landmark inclusion/exclusion partition does not reconcile to 4,380.")
}

# Patient-level/time objects no longer needed after protected scenario frames.
rm(db, local, mi, lpa, votes, posterior, v1, v2, profile, icu_in, icu_out,
   discharge, death_exact, death_date, date_lower, date_upper, effective_lower,
   midpoint, latest_legal, L, H, icu_shortfall_minutes)

control_settings <- list(
  maxit = 10000L, use.gr = TRUE,
  reltol.search = 1e-8, reltol.final = 1e-10, reltol.outer = 1e-8
)
full_formula <- survival::Surv(time, event) ~
  profile2 + age10 + sex + charlson + cancer_phenotype + sofa
display_icu_day <- c(2, 7, 14, 21, 28)
prediction_time <- c(1, 6, 13, 20, 27)
curve_time <- sort(unique(c(seq(0.25, 27, by = 0.25), prediction_time)))
z975 <- stats::qnorm(0.975)

make_reference <- function(data, time) {
  data.frame(
    time = time,
    profile2 = 0L,
    age10 = 0,
    sex = factor("F", levels = levels(data$sex)),
    charlson = 0,
    cancer_phenotype = factor("metastatic solid", levels = levels(data$cancer_phenotype)),
    sofa = 0,
    subject_id = data$subject_id[1L],
    check.names = FALSE
  )
}

fit_one <- function(data, scenario, time_varying) {
  warning_text <- character()
  fit <- withCallingHandlers(
    rstpm2::stpm2(
      full_formula, data = data, df = 4L, scale = "hazard",
      tvc = if (time_varying) list(profile2 = 2L) else NULL,
      robust = TRUE, cluster = "subject_id", frailty = FALSE,
      control = control_settings
    ),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  b <- stats::coef(fit)
  V <- stats::vcov(fit)
  conv <- fit@details$convergence
  gradient <- tryCatch(fit@args$gradnegll(b), error = function(e) rep(NA_real_, length(b)))
  Hnum <- tryCatch(
    stats::optimHess(b, fn = fit@minuslogl, gr = fit@args$gradnegll),
    error = function(e) matrix(NA_real_, length(b), length(b))
  )
  h_eigen <- tryCatch(eigen((Hnum + t(Hnum)) / 2, symmetric = TRUE, only.values = TRUE)$values,
                      error = function(e) NA_real_)
  v_eigen <- tryCatch(eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values,
                      error = function(e) NA_real_)
  grid_ref <- make_reference(data, curve_time)
  h0 <- as.numeric(predict(fit, newdata = grid_ref, type = "hazard"))
  s0 <- as.numeric(predict(fit, newdata = grid_ref, type = "surv"))
  grid_exp <- grid_ref
  grid_exp$profile2 <- 1L
  h1 <- as.numeric(predict(fit, newdata = grid_exp, type = "hazard"))
  s1 <- as.numeric(predict(fit, newdata = grid_exp, type = "surv"))
  hr_pred <- predict(
    fit, newdata = grid_ref, type = "hr", var = "profile2",
    exposed = function(x) { x$profile2 <- 1L; x },
    se.fit = TRUE, level = 0.95
  )
  if (!all(c("Estimate", "lower", "upper") %in% names(hr_pred))) {
    stop("Unexpected rstpm2 HR prediction columns for ", scenario)
  }
  hr <- data.frame(
    scenario = scenario,
    days_after_landmark = curve_time,
    ICU_day = curve_time + 1,
    N = nrow(data),
    subjects = length(unique(data$subject_id)),
    events = sum(data$event),
    risk_stays = vapply(curve_time, function(t) sum(data$time >= t), integer(1)),
    risk_subjects = vapply(curve_time, function(t) length(unique(data$subject_id[data$time >= t])), integer(1)),
    HR = as.numeric(hr_pred$Estimate),
    CI_lower = as.numeric(hr_pred$lower),
    CI_upper = as.numeric(hr_pred$upper),
    stringsAsFactors = FALSE
  )
  hr$SE_logHR <- (log(hr$CI_upper) - log(hr$CI_lower)) / (2 * z975)
  finite_hr <- all(is.finite(as.matrix(hr[c("HR", "SE_logHR", "CI_lower", "CI_upper")]))) &&
    all(hr$HR > 0) && all(hr$CI_lower > 0) && all(hr$CI_upper >= hr$CI_lower)
  hazard_ok <- all(is.finite(c(h0, h1))) && all(c(h0, h1) >= -1e-12)
  surv_ok <- all(is.finite(c(s0, s1))) && all(c(s0, s1) >= -1e-10) &&
    all(c(s0, s1) <= 1 + 1e-10) && all(diff(s0) <= 1e-8) && all(diff(s1) <= 1e-8)
  grad_max <- if (all(is.finite(gradient))) max(abs(gradient)) else NA_real_
  h_min <- if (all(is.finite(h_eigen))) min(h_eigen) else NA_real_
  h_max <- if (all(is.finite(h_eigen))) max(h_eigen) else NA_real_
  v_min <- if (all(is.finite(v_eigen))) min(v_eigen) else NA_real_
  v_max <- if (all(is.finite(v_eigen))) max(v_eigen) else NA_real_
  pass <- identical(as.integer(conv), 0L) && all(is.finite(b)) && all(is.finite(V)) &&
    is.finite(grad_max) && is.finite(h_min) && h_min > 0 &&
    is.finite(v_min) && v_min > 0 && hazard_ok && surv_ok && finite_hr
  diagnostics <- data.frame(
    scenario = scenario,
    model = if (time_varying) "SOFA-adjusted Profile time-varying" else "SOFA-adjusted Profile-constant comparison",
    N = nrow(data),
    subjects = length(unique(data$subject_id)),
    events = sum(data$event),
    profile_1_n = sum(data$profile2 == 0L),
    profile_2_n = sum(data$profile2 == 1L),
    optimizer = as.character(fit@optimizer),
    convergence_code = as.integer(conv),
    max_abs_gradient = grad_max,
    gradient_threshold = NA_real_,
    gradient_acceptance_rule = "finite and reported; R1B specified no numeric raw-gradient cutoff",
    hessian_min_eigenvalue = h_min,
    hessian_max_eigenvalue = h_max,
    robust_vcov_min_eigenvalue = v_min,
    robust_vcov_max_eigenvalue = v_max,
    coefficients_finite = all(is.finite(b)),
    robust_covariance_finite = all(is.finite(V)),
    hazard_nonnegative = hazard_ok,
    survival_in_0_1_and_nonincreasing = surv_ok,
    HR_and_pointwise_CI_finite = finite_hr,
    min_risk_stays_at_fixed_points = min(vapply(prediction_time, function(t) sum(data$time >= t), integer(1))),
    min_risk_subjects_at_fixed_points = min(vapply(prediction_time, function(t) length(unique(data$subject_id[data$time >= t])), integer(1))),
    warnings = if (length(warning_text)) paste(unique(warning_text), collapse = " | ") else "none",
    overall_numeric_status = if (pass) "PASS" else "FAIL",
    baseline_spline_df = 4L,
    profile_time_interaction_df = if (time_varying) 2L else 0L,
    variance = "subject-clustered sandwich covariance",
    cluster_unit = "patient (subject_id in protected memory)",
    frailty = FALSE,
    stringsAsFactors = FALSE
  )
  if (!pass) {
    utils::write.csv(diagnostics,
                     file.path(out_dir, paste0("R1C_FAILED_", scenario, "_diagnostics.csv")),
                     row.names = FALSE, na = "")
    stop("Frozen model failed numeric acceptance for ", scenario,
         if (!time_varying) " constant comparison" else " time-varying model")
  }
  list(fit = fit, hr = hr, diagnostics = diagnostics)
}

# Formal execution: midpoint time-varying, same-risk-set Profile-constant
# comparison, then earliest/latest date-precision time-varying scenarios.
fit_mid <- fit_one(scenario_data$midpoint, "midpoint", TRUE)
fit_const <- fit_one(scenario_data$midpoint, "midpoint", FALSE)
fit_early <- fit_one(scenario_data$earliest, "earliest", TRUE)
fit_late <- fit_one(scenario_data$latest, "latest", TRUE)
fits_tv <- list(midpoint = fit_mid, earliest = fit_early, latest = fit_late)

interaction_test <- function(item, scenario) {
  fit <- item$fit
  b <- stats::coef(fit)
  idx <- which(grepl("profile2", names(b), fixed = TRUE) & names(b) != "profile2")
  if (length(idx) != 2L) {
    stop("Expected exactly two Profile time-interaction coefficients for ", scenario,
         "; found: ", paste(names(b)[idx], collapse = ", "))
  }
  V <- stats::vcov(fit)[idx, idx, drop = FALSE]
  stat <- as.numeric(crossprod(b[idx], solve(V, b[idx])))
  data.frame(
    scenario = scenario,
    N = item$diagnostics$N,
    subjects = item$diagnostics$subjects,
    events = item$diagnostics$events,
    test = "patient-clustered robust joint Wald test of Profile time interaction",
    coefficients_tested = paste(names(b)[idx], collapse = " | "),
    chi_square = stat,
    df = length(idx),
    P_value = stats::pchisq(stat, df = length(idx), lower.tail = FALSE),
    covariance = "subject-clustered sandwich covariance",
    ordinary_LRT_used = FALSE,
    stringsAsFactors = FALSE
  )
}
time_tests <- do.call(rbind, Map(interaction_test, fits_tv, names(fits_tv)))

all_curve <- do.call(rbind, lapply(fits_tv, `[[`, "hr"))
all_curve$inference <- "conditional instantaneous HR with patient-clustered robust pointwise 95% CI"
all_curve$prediction_reference <- sprintf(
  "Profile 2 vs 1; age %.4f y; Charlson %.4f; SOFA %.4f; sex F; metastatic solid",
  reference_means["age_years"], reference_means["charlson"], reference_means["sofa"]
)
point_rows <- all_curve[all_curve$days_after_landmark %in% prediction_time, , drop = FALSE]
point_rows <- point_rows[order(match(point_rows$scenario, c("midpoint", "earliest", "latest")),
                               point_rows$days_after_landmark), ]

constant_b <- stats::coef(fit_const$fit)
constant_V <- stats::vcov(fit_const$fit)
constant_logHR <- unname(constant_b["profile2"])
constant_se <- sqrt(constant_V["profile2", "profile2"])
comparison <- rbind(
  data.frame(
    scenario = "midpoint", model = "SOFA-adjusted Profile time-varying",
    N = 4340L, subjects = 3742L, events = 1500L,
    baseline_spline_df = 4L, profile_time_interaction_df = 2L,
    profile_HR = NA_real_, SE_logHR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
    P_value = time_tests$P_value[time_tests$scenario == "midpoint"],
    P_value_meaning = "robust joint Wald test of Profile time interaction",
    log_likelihood = -as.numeric(fit_mid$fit@min),
    AIC = 2 * length(stats::coef(fit_mid$fit)) + 2 * as.numeric(fit_mid$fit@min),
    BIC = log(4340) * length(stats::coef(fit_mid$fit)) + 2 * as.numeric(fit_mid$fit@min),
    variance = "patient-clustered robust", nested_test_with_primary_4380 = FALSE,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = "midpoint", model = "SOFA-adjusted Profile-constant association comparison",
    N = 4340L, subjects = 3742L, events = 1500L,
    baseline_spline_df = 4L, profile_time_interaction_df = 0L,
    profile_HR = exp(constant_logHR), SE_logHR = constant_se,
    CI_lower = exp(constant_logHR - z975 * constant_se),
    CI_upper = exp(constant_logHR + z975 * constant_se),
    P_value = 2 * stats::pnorm(-abs(constant_logHR / constant_se)),
    P_value_meaning = "patient-clustered robust Wald test of constant Profile coefficient",
    log_likelihood = -as.numeric(fit_const$fit@min),
    AIC = 2 * length(stats::coef(fit_const$fit)) + 2 * as.numeric(fit_const$fit@min),
    BIC = log(4340) * length(stats::coef(fit_const$fit)) + 2 * as.numeric(fit_const$fit@min),
    variance = "patient-clustered robust", nested_test_with_primary_4380 = FALSE,
    stringsAsFactors = FALSE
  )
)

date_sensitivity <- reshape(
  point_rows[c("ICU_day", "days_after_landmark", "scenario", "HR", "CI_lower", "CI_upper")],
  idvar = c("ICU_day", "days_after_landmark"), timevar = "scenario", direction = "wide"
)
date_sensitivity <- date_sensitivity[order(date_sensitivity$days_after_landmark), ]
date_sensitivity$minimum_HR_across_scenarios <- apply(
  date_sensitivity[c("HR.midpoint", "HR.earliest", "HR.latest")], 1, min
)
date_sensitivity$maximum_HR_across_scenarios <- apply(
  date_sensitivity[c("HR.midpoint", "HR.earliest", "HR.latest")], 1, max
)
date_sensitivity$absolute_HR_range <-
  date_sensitivity$maximum_HR_across_scenarios - date_sensitivity$minimum_HR_across_scenarios
date_sensitivity$relative_HR_range_vs_midpoint_percent <-
  100 * date_sensitivity$absolute_HR_range / date_sensitivity$HR.midpoint
date_sensitivity$direction_consistent <- apply(
  date_sensitivity[c("HR.midpoint", "HR.earliest", "HR.latest")], 1,
  function(x) if (all(x > 1)) "all HR >1" else if (all(x < 1)) "all HR <1" else "crosses 1 across scenarios"
)
date_sensitivity$interpretation <-
  "Date-precision scenarios; not confidence bounds and not a new confidence interval"

model_diagnostics <- do.call(rbind, list(
  fit_mid$diagnostics, fit_const$diagnostics, fit_early$diagnostics, fit_late$diagnostics
))

# Restrict exported schemas to aggregate fields.
export_objects <- list(
  points = point_rows,
  curve = all_curve,
  comparison = comparison,
  tests = time_tests,
  date_sensitivity = date_sensitivity,
  diagnostics = model_diagnostics,
  exclusions = exclusion_breakdown
)
for (nm in names(export_objects)) {
  bad <- grepl("^(subject_id|stay_id|hadm_id|patient_id|death_date|timestamp)$",
               names(export_objects[[nm]]), ignore.case = TRUE)
  if (any(bad)) stop("Patient-level field detected in aggregate export schema: ", nm)
}

utils::write.csv(point_rows, output_paths[1], row.names = FALSE, na = "")
utils::write.csv(all_curve, output_paths[2], row.names = FALSE, na = "")
utils::write.csv(comparison, output_paths[3], row.names = FALSE, na = "")
utils::write.csv(time_tests, output_paths[4], row.names = FALSE, na = "")
utils::write.csv(date_sensitivity, output_paths[5], row.names = FALSE, na = "")
utils::write.csv(model_diagnostics, output_paths[6], row.names = FALSE, na = "")
utils::write.csv(exclusion_breakdown, output_paths[7], row.names = FALSE, na = "")

draw_main <- function() {
  d <- all_curve[all_curve$scenario == "midpoint", ]
  fixed <- point_rows[point_rows$scenario == "midpoint", ]
  yr <- range(c(d$CI_lower, d$CI_upper, 1), finite = TRUE)
  plot(NA, xlim = c(1, 28), ylim = yr, log = "y", xaxt = "n",
       xlab = "ICU day", ylab = "Hazard ratio for Profile 2 vs Profile 1",
       main = "Time-varying profile association after the 24-hour landmark",
       las = 1, bty = "l")
  axis(1, at = c(2, 7, 14, 21, 28), labels = c("2", "7", "14", "21", "28"))
  polygon(c(d$ICU_day, rev(d$ICU_day)), c(d$CI_lower, rev(d$CI_upper)),
          col = grDevices::adjustcolor("#2A6FBB", alpha.f = 0.18), border = NA)
  lines(d$ICU_day, d$HR, col = "#1F5A94", lwd = 2.2)
  points(fixed$ICU_day, fixed$HR, pch = 21, bg = "white", col = "#1F5A94", cex = 1.0)
  abline(h = 1, col = "#555555", lty = 2)
  abline(v = 1, col = "#333333", lty = 3)
  text(1.15, exp(log(yr[1]) + 0.92 * diff(log(yr))), "24-h\nlandmark",
       pos = 4, cex = 0.72, col = "#333333")
  legend("topright", legend = c("Conditional HR(t)", "95% pointwise CI", "No association"),
         col = c("#1F5A94", grDevices::adjustcolor("#2A6FBB", alpha.f = 0.35), "#555555"),
         lwd = c(2.2, 8, 1), lty = c(1, 1, 2), bty = "n", cex = 0.86)
  mtext("Midpoint date convention; SOFA-adjusted; patient-clustered robust covariance", side = 3, line = 0.25, cex = 0.78)
}

draw_date <- function() {
  cols <- c(midpoint = "#1F5A94", earliest = "#D97706", latest = "#5B4B8A")
  ltys <- c(midpoint = 1, earliest = 2, latest = 3)
  yr <- range(c(all_curve$HR, 1), finite = TRUE)
  plot(NA, xlim = c(1, 28), ylim = yr, log = "y", xaxt = "n",
       xlab = "ICU day", ylab = "Hazard ratio for Profile 2 vs Profile 1",
       main = "Date-precision sensitivity of the time-varying association",
       las = 1, bty = "l")
  axis(1, at = c(2, 7, 14, 21, 28), labels = c("2", "7", "14", "21", "28"))
  for (s in c("midpoint", "earliest", "latest")) {
    d <- all_curve[all_curve$scenario == s, ]
    lines(d$ICU_day, d$HR, col = cols[s], lty = ltys[s], lwd = 2)
  }
  abline(h = 1, col = "#555555", lty = 2)
  abline(v = 1, col = "#333333", lty = 3)
  text(1.15, exp(log(yr[1]) + 0.92 * diff(log(yr))), "24-h\nlandmark",
       pos = 4, cex = 0.72, col = "#333333")
  legend("topright", legend = c("Midpoint", "Earliest legitimate time", "Latest legitimate time"),
         col = cols, lty = ltys, lwd = 2, bty = "n", cex = 0.84)
  mtext("Scenario curves are date-precision sensitivity analyses, not confidence bounds", side = 3, line = 0.25, cex = 0.78)
}

grDevices::pdf(output_paths[8], width = 7.2, height = 5.4, useDingbats = FALSE,
               family = "sans", paper = "special")
draw_main()
grDevices::dev.off()
grDevices::tiff(output_paths[9], width = 7.2, height = 5.4, units = "in",
                res = 300, compression = "lzw", bg = "white")
draw_main()
grDevices::dev.off()
grDevices::pdf(output_paths[10], width = 7.2, height = 5.4, useDingbats = FALSE,
               family = "sans", paper = "special")
draw_date()
grDevices::dev.off()
grDevices::tiff(output_paths[11], width = 7.2, height = 5.4, units = "in",
                res = 300, compression = "lzw", bg = "white")
draw_date()
grDevices::dev.off()

if (!all(file.exists(output_paths)) || any(file.info(output_paths)$size <= 0)) {
  stop("One or more R1C aggregate/figure outputs were not created correctly.")
}

# Remove protected frames and fitted objects before normal exit.
rm(scenario_data, base_covariates, fit_mid, fit_const, fit_early, fit_late,
   fits_tv, retain, exact_after_L, exact_at_or_before_L, date_only,
   other_time_conflict, exclude_pre_or_at_L, landmark_uncertain,
   exact_date_conflict, invalid_exact_for_event, missing_event_time,
   death_evidence_with_zero_flag, icu_duration_conflict)
cat("R1C_TIME_VARYING_ANALYSIS_COMPLETE aggregate outputs and figures written.\n")
