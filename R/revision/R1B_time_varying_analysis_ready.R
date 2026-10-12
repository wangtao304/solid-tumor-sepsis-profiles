#!/usr/bin/env Rscript

# Reviewer Revision R1B: freeze and audit the 24-hour landmark population.
# This script deliberately DOES NOT fit a time-varying or other outcome model.
# Identifiers/timestamps exist only in memory; outputs contain aggregates only.

options(stringsAsFactors = FALSE, warn = 1)

if (!identical(Sys.getenv("R1B_EXECUTE_LANDMARK_AUDIT"), "YES")) {
  cat("R1B landmark audit execution gate is closed; no audit or model executed.\n")
  quit(save = "no", status = 0L)
}

suppressPackageStartupMessages(library(mice))

mice_path <- "data/processed/inflammatory_metabolic_mice_solid_m10_maxit20_seed20260721.rds"
lpa_path <- "data/processed/inflammatory_metabolic_lpa_solid_2to6_mi10_seed20260722.rds"
sql_path <- "sql/revision/R1B_landmark_boundary_input.sql"
population_output <- "outputs/revision/R1B_landmark_population_audit.csv"
boundary_output <- "outputs/revision/R1B_death_date_boundary_audit.csv"
future_model_output <- "outputs/revision/R1B_time_varying_profile_results.csv"
dir.create(dirname(population_output), recursive = TRUE, showWarnings = FALSE)
model_execution_requested <- identical(
  Sys.getenv("R1B_EXECUTE_TIME_VARYING_MODEL"), "SEPARATELY_AUTHORIZED"
)

if (!all(file.exists(c(mice_path, lpa_path, sql_path)))) {
  stop("A required frozen input or reviewed SQL source is missing.")
}
audit_outputs_exist <- file.exists(population_output) && file.exists(boundary_output)
if (xor(file.exists(population_output), file.exists(boundary_output))) {
  stop("Only one landmark audit output exists; resolve the incomplete state first.")
}
if (audit_outputs_exist && !model_execution_requested) {
  stop("Refusing to overwrite existing R1B landmark audit outputs.")
}
if (file.exists(future_model_output)) {
  stop("Refusing to overwrite an existing time-varying model output.")
}

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

# Read-only patient/time linkage. Raw fields are discarded before any output.
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
args <- c(
  "-h", db_host,
  "-p", db_port,
  "-U", db_user,
  "-d", db_name,
  "-X", "-w", "-v", "ON_ERROR_STOP=1", "-q", "--csv", "-f", sql_path
)
txt <- system2(psql_bin, args, stdout = TRUE, stderr = "")
if (!is.null(attr(txt, "status")) && attr(txt, "status") != 0L) {
  stop("Read-only R1B landmark-boundary query failed.")
}
db <- utils::read.csv(text = paste(txt, collapse = "\n"), check.names = FALSE, na.strings = "")
rm(txt)

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
) stop("R1B landmark linkage/cohort gate failed.")

# All MIMIC timestamps are interpreted on the common database clock. Assigning
# tz='UTC' here prevents host-local daylight-saving conversion; it does not
# convert or overwrite a database timestamp.
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
  exact_in_followup |
  (has_date & date_upper > icu_in & date_lower <= H)
)
icu_duration_conflict <- !is.na(icu_out) & icu_out < L

other_time_conflict <- exact_date_conflict | invalid_exact_for_event |
  missing_event_time | death_evidence_with_zero_flag | icu_duration_conflict

# Use date-only timing only when the 28-day event is positive, no exact
# timestamp is available, and the record has passed the conflict screen.
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

# Retain the requested precision metadata in protected memory for deterministic
# future preparation. These patient-level fields are removed before CSV output.
db$death_interval_lower <- date_lower
db$death_interval_upper <- date_upper
db$time_precision <- ifelse(
  db$event_28d == 0L, "administrative_censoring",
  ifelse(exact_in_followup, "exact_timestamp",
         ifelse(date_only, "date_only_interval", "time_conflict"))
)
db$landmark_boundary_flag <- date_spans_L
db$day28_boundary_flag <- date_spans_H

# Landmark membership. A record at or before L contributes no post-landmark
# risk time. Make the reported exclusion categories mutually exclusive: a row
# with an independent time conflict is counted under that conflict, not twice.
exclude_pre_or_at_L <- (exact_at_or_before_L | date_before_or_at_L) &
  !other_time_conflict
retain <- !exclude_pre_or_at_L & !landmark_uncertain & !other_time_conflict

icu_shortfall_minutes <- as.numeric(difftime(L, icu_out, units = "mins"))

# Working midpoint for date-only deaths after landmark. For an L-spanning date
# interval resolved by credible post-L observation, the legal lower bound is L.
effective_lower <- date_lower
effective_lower[date_L_resolved_by_observation] <- L[date_L_resolved_by_observation]
midpoint <- effective_lower + as.numeric(difftime(date_upper, effective_lower, units = "secs")) / 2

working_event <- integer(nrow(db))
working_event[exact_after_L & retain] <- 1L
working_event[date_only & retain & midpoint <= H] <- 1L

earliest_event <- integer(nrow(db))
earliest_event[exact_after_L & retain] <- 1L
earliest_event[date_only & retain & effective_lower <= H] <- 1L

latest_event <- integer(nrow(db))
latest_event[exact_after_L & retain] <- 1L
latest_event[date_only & retain & date_upper <= H] <- 1L

working_death_time <- death_exact
working_death_time[date_only] <- midpoint[date_only]
analysis_end <- H
analysis_end[working_event == 1L] <- working_death_time[working_event == 1L]
follow_up_days <- as.numeric(difftime(analysis_end, L, units = "days"))

if (
  anyNA(profile) || anyNA(retain) || anyNA(working_event) ||
  any(follow_up_days[retain] <= 0) || any(follow_up_days[retain] > 27 + 1e-10)
) stop("The frozen landmark preparation produced an invalid risk time.")

metric <- function(name, value, unit, status, note) {
  data.frame(metric = name, value = value, unit = unit, status = status, note = note,
             stringsAsFactors = FALSE)
}

population <- do.call(rbind, list(
  metric("original_candidate_stays", nrow(db), "stays", ifelse(nrow(db) == 4380L, "PASS", "FAIL"), "Frozen solid-tumor MIMIC cohort."),
  metric("original_unique_subjects", length(unique(db$subject_id)), "subjects", ifelse(length(unique(db$subject_id)) == 3772L, "PASS", "FAIL"), "Identifiers used only in protected memory."),
  metric("original_profile_1_stays", sum(profile == 1L), "stays", ifelse(sum(profile == 1L) == 3197L, "PASS", "FAIL"), "Frozen consensus assignment."),
  metric("original_profile_2_stays", sum(profile == 2L), "stays", ifelse(sum(profile == 2L) == 1183L, "PASS", "FAIL"), "Frozen consensus assignment."),
  metric("original_28_day_death_records", sum(db$event_28d), "records", ifelse(sum(db$event_28d) == 1537L, "PASS", "FAIL"), "Frozen outcome flag before landmark reclassification."),
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
  metric("landmark_midpoint_working_events", sum(working_event[retain]), "events", "WORKING_CONVENTION", "Exact timestamps plus midpoint convention for eligible date-only intervals."),
  metric("landmark_earliest_scenario_events", sum(earliest_event[retain]), "events", "SENSITIVITY_SCENARIO", "Earliest legal time for date-only intervals."),
  metric("landmark_latest_scenario_events", sum(latest_event[retain]), "events", "SENSITIVITY_SCENARIO", "Latest legal time for date-only intervals."),
  metric("maximum_landmark_follow_up_days", max(follow_up_days[retain]), "days", ifelse(max(follow_up_days[retain]) <= 27 + 1e-10, "PASS", "FAIL"), "H remains ICU admission +28 days; maximum post-L follow-up is 27 days."),
  metric("minimum_landmark_follow_up_days", min(follow_up_days[retain]), "days", ifelse(min(follow_up_days[retain]) > 0, "PASS", "FAIL"), "No zero or negative post-landmark risk time."),
  metric("covariate_row_link_checks_passed", sum(link_checks), "checks", ifelse(all(link_checks), "PASS", "FAIL"), paste(names(link_checks), collapse = "; "))
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
  metric("midpoint_minus_earliest_event_count", sum(working_event[retain]) - sum(earliest_event[retain]), "events", "SCENARIO_DIFFERENCE", "Midpoint working convention minus earliest legal-time scenario."),
  metric("latest_minus_midpoint_event_count", sum(latest_event[retain]) - sum(working_event[retain]), "events", "SCENARIO_DIFFERENCE", "Latest legal-time scenario minus midpoint working convention."),
  metric("earliest_minus_latest_event_count", sum(earliest_event[retain]) - sum(latest_event[retain]), "events", "SCENARIO_RANGE", "Maximum aggregate event-count difference induced by date precision."),
  metric("raw_timestamp_or_identifier_fields_exported", 0L, "fields", "PASS", "All identifiers and raw date/time fields were discarded before writing aggregate files.")
))

if (any(grepl("subject_id|stay_id|hadm_id|death_date|timestamp", names(population), ignore.case = TRUE)) ||
    any(grepl("subject_id|stay_id|hadm_id|death_date|timestamp", names(boundary), ignore.case = TRUE))) {
  stop("Patient-level field detected in aggregate output schema.")
}

# Model-ready protected-memory frame. Its construction is part of R1B data
# preparation; it is never saved. Model execution remains under a second,
# separately authorized gate below.
analysis_data <- data.frame(
  subject_id = db$subject_id[retain],
  profile2 = as.integer(profile[retain] == 2L),
  age10 = (as.numeric(local$age[retain]) - mean(as.numeric(local$age[retain]))) / 10,
  sex = factor(local$sex[retain]),
  charlson = as.numeric(local$charlson_comorbidity_index[retain]) -
    mean(as.numeric(local$charlson_comorbidity_index[retain])),
  cancer_phenotype = factor(local$cancer_phenotype[retain]),
  sofa = as.numeric(local$sofa[retain]) - mean(as.numeric(local$sofa[retain])),
  sapsii = as.numeric(local$sapsii[retain]) - mean(as.numeric(local$sapsii[retain])),
  time = follow_up_days[retain],
  event = working_event[retain],
  check.names = FALSE
)
if (
  nrow(analysis_data) != 4340L || anyNA(analysis_data) ||
  !all(analysis_data$profile2 %in% 0:1) ||
  any(analysis_data$time <= 0) || any(analysis_data$time > 27)
) stop("Model-ready protected-memory frame failed the frozen QC gates.")

# Remove protected row-level objects before writing aggregate-only CSVs.
rm(db, local, icu_in, icu_out, discharge, death_exact, death_date, date_lower,
   date_upper, L, H, working_death_time, effective_lower, midpoint, profile,
   votes, posterior, mi, lpa, icu_shortfall_minutes)
if (!audit_outputs_exist) {
  utils::write.csv(population, population_output, row.names = FALSE, na = "")
  utils::write.csv(boundary, boundary_output, row.names = FALSE, na = "")
} else {
  saved_population <- utils::read.csv(population_output, check.names = FALSE)
  saved_boundary <- utils::read.csv(boundary_output, check.names = FALSE)
  if (!isTRUE(all.equal(saved_population, population, check.attributes = FALSE)) ||
      !isTRUE(all.equal(saved_boundary, boundary, check.attributes = FALSE))) {
    stop("Recomputed landmark audit does not match the frozen aggregate outputs.")
  }
}

if (!model_execution_requested) {
  cat("R1B_LANDMARK_AUDIT_COMPLETE aggregate outputs written; time-varying outcome model not fitted.\n")
  quit(save = "no", status = 0L)
}

# The code below is prepared for a separately authorized future stage. It was
# not executed in R1B. Settings are fixed and must not be tuned to the result.
suppressPackageStartupMessages({
  library(survival)
  library(rstpm2)
})
if (packageVersion("rstpm2") != package_version("1.7.1")) {
  stop("The future model was reviewed for rstpm2 1.7.1; another version requires re-review.")
}

spec <- data.frame(
  model = c(
    "Model 1", "Model 2", "Model 3", "Model 4",
    "Sensitivity: SAPS-II instead of SOFA"
  ),
  rhs = c(
    "profile2",
    "profile2 + age10 + sex",
    "profile2 + age10 + sex + charlson + cancer_phenotype",
    "profile2 + age10 + sex + charlson + cancer_phenotype + sofa",
    "profile2 + age10 + sex + charlson + cancer_phenotype + sapsii"
  ),
  stringsAsFactors = FALSE
)
control_settings <- list(
  maxit = 10000L, use.gr = TRUE,
  reltol.search = 1e-8, reltol.final = 1e-10, reltol.outer = 1e-8
)
display_icu_day <- c(2, 7, 14, 21, 28)
prediction_time <- c(1, 6, 13, 20, 27)
reference_row <- data.frame(
  time = prediction_time,
  profile2 = 0L,
  age10 = 0,
  sex = factor(levels(analysis_data$sex)[1L], levels = levels(analysis_data$sex)),
  charlson = 0,
  cancer_phenotype = factor(
    levels(analysis_data$cancer_phenotype)[1L],
    levels = levels(analysis_data$cancer_phenotype)
  ),
  sofa = 0,
  sapsii = 0,
  subject_id = analysis_data$subject_id[1L],
  check.names = FALSE
)

future_rows <- vector("list", nrow(spec))
for (i in seq_len(nrow(spec))) {
  formula_i <- stats::as.formula(paste("survival::Surv(time, event) ~", spec$rhs[i]))
  fit <- rstpm2::stpm2(
    formula_i, data = analysis_data, df = 4L, scale = "hazard",
    tvc = list(profile2 = 2L),
    robust = TRUE, cluster = "subject_id", frailty = FALSE,
    control = control_settings
  )
  if (
    fit@details$convergence != 0L || any(!is.finite(stats::coef(fit))) ||
    any(!is.finite(stats::vcov(fit)))
  ) stop("A future time-varying model failed the frozen convergence criteria.")
  pred <- stats::predict(
    fit, newdata = reference_row, type = "hr", var = "profile2",
    exposed = function(data) { data$profile2 <- 1L; data },
    se.fit = TRUE, level = 0.95
  )
  if (!all(c("Estimate", "lower", "upper") %in% names(pred)) ||
      any(!is.finite(as.matrix(pred[c("Estimate", "lower", "upper")]))) ) {
    stop("Unexpected or nonfinite future time-varying HR prediction.")
  }
  future_rows[[i]] <- data.frame(
    model = spec$model[i],
    time_origin = "ICU admission +24-hour landmark",
    icu_day = display_icu_day,
    days_after_landmark = prediction_time,
    N_stays = nrow(analysis_data),
    N_subjects = length(unique(analysis_data$subject_id)),
    events_midpoint_convention = sum(analysis_data$event),
    HR = pred$Estimate,
    CI_lower = pred$lower,
    CI_upper = pred$upper,
    baseline_spline_df = 4L,
    profile_time_interaction_df = 2L,
    variance_type = "subject-clustered robust",
    date_only_event_time_rule = "midpoint; earliest/latest scenarios required separately",
    stringsAsFactors = FALSE
  )
}
future_result <- do.call(rbind, future_rows)
if (any(grepl("subject_id|stay_id|hadm_id", names(future_result), ignore.case = TRUE))) {
  stop("Identifier field detected in future aggregate model output.")
}
utils::write.csv(future_result, future_model_output, row.names = FALSE, na = "")
cat("R1B_TIME_VARYING_MODEL_COMPLETE aggregate output written.\n")
