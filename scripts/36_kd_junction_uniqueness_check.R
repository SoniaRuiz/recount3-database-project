## ============================================================================
## kd_junction_uniqueness_check_fast.R
##
## Question: how many junctions detected in the knockdown (KD) database are NOT
## found in any other IntroVerse database?  If almost none, a line such as
## "also detected in N non-KD projects" cannot work as a safeguard on the card,
## because large projects (TCGA, GTEx, ...) already contain nearly everything.
##
## WHY THIS IS FAST
## Each database stores one row per DISTINCT junction in two small master
## tables, so the junction lists can be compared directly:
##   novel  - novel donor / novel acceptor junctions   (novel_coordinates, novel_type)
##   other  - unannotated / novel-combination junctions (ref_coordinates, type)
## local_splicing_activity (one row per junction x project x cluster, no index
## on the coordinate columns) is NEVER read, so there is no multi-GB scan.
## Annotated introns are left out on purpose: they come from the shared
## reference annotation, not from detection, so they cannot be "unique to KD".
##
## WHAT THIS DOES NOT REPORT (these need local_splicing_activity, which is the
## slow table - ask for a restricted follow-up on just the unique junctions):
##   - how many individuals outside the KD data carry a junction
##   - whether a junction is KD-only versus present in the matched control arm
## "KD data" here means the whole KD database (KD + control clusters together).
##
## Junction keys are coordinate strings ("chr1:100-200:+"); integer junction
## ids are database-specific and cannot be compared across databases.
## ============================================================================
run = FALSE

if (run) {
  
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
  })

  ## ---- Parameters (edit these) ----------------------------------------------
  database_folder <- "/home/sg2173/PROJECTS/srRNAseq/introverse_v2/database"   # same as config.R
  ensembl_version <- "115"

  kd_db <- "ENCODE_SR_1read_shRNA.sqlite"                # the KD data set

  ## Databases used as "elsewhere". The NMD KD database (SRP083135) is left out
  ## on purpose: it is itself a knockdown data set, so seeing a junction there
  ## says nothing about it occurring outside a KD. Add it back if you disagree.
  other_dbs <- c("TCGA_1read_subsampleFALSE.sqlite",
                "GTEX_1read_subsampleFALSE.sqlite",
                "SRP058181_1read_subsampleFALSE.sqlite",
                "SRP100948_1read_subsampleFALSE.sqlite")

  out_rds <- "kd_junction_uniqueness_fast.rds"
  ## -----------------------------------------------------------------------------

  stopifnot(length(other_dbs) <= 9)      # SQLite allows 10 attached DBs by default

  ## Same path rule as database_connect() in utils.R:
  ##   <database_folder>/<file name without .sqlite>/<ensembl_version>/<file name>
  db_path <- function(db_file) {
    file.path(database_folder, sub("\\.sqlite$", "", db_file), ensembl_version, db_file)
  }

  step <- function(msg) message(format(Sys.time(), "%H:%M:%S"), "  ", msg)
  t_start <- Sys.time()

  ## SQLite silently CREATES an empty database when given a path that does not
  ## exist, so check every file first.
  all_paths <- vapply(c(kd_db, other_dbs), db_path, character(1))
  not_found <- all_paths[!file.exists(all_paths)]
  if (length(not_found) > 0) {
    stop("These database files do not exist (check database_folder / ensembl_version):\n  ",
        paste(not_found, collapse = "\n  "))
  }
  empty <- all_paths[file.size(all_paths) == 0]
  if (length(empty) > 0) {
    stop("These database files are empty (0 bytes) - delete them and fix database_folder:\n  ",
        paste(empty, collapse = "\n  "))
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path(kd_db))
  aliases <- paste0("o", seq_along(other_dbs))
  for (i in seq_along(other_dbs)) {
    DBI::dbExecute(con, sprintf("ATTACH DATABASE '%s' AS %s", db_path(other_dbs[i]), aliases[i]))
  }

  ## ---- Sanity check: the two master tables and the columns used exist ---------
  required <- list(novel = c("novel_coordinates", "novel_type"),
                  other = c("ref_coordinates", "type"))
  for (s in c("main", aliases)) {
    db_name <- if (s == "main") kd_db else other_dbs[match(s, aliases)]
    tables  <- DBI::dbGetQuery(con, sprintf("SELECT name FROM %s.sqlite_master WHERE type = 'table'", s))$name
    for (tb in names(required)) {
      if (!(tb %in% tables)) {
        stop(sprintf("%s has no table '%s'. Tables found: %s", db_name, tb, paste(tables, collapse = ", ")))
      }
      cols <- DBI::dbGetQuery(con, sprintf("PRAGMA %s.table_info('%s')", s, tb))$name
      miss <- setdiff(required[[tb]], cols)
      if (length(miss) > 0) {
        stop(sprintf("%s: table '%s' lacks column(s): %s. Columns found: %s",
                    db_name, tb, paste(miss, collapse = ", "), paste(cols, collapse = ", ")))
      }
    }
  }

  ## ---- KD junction list --------------------------------------------------------
  step("KD database: reading its junction list (novel + other tables)...")
  invisible(DBI::dbExecute(con, '
    CREATE TEMP TABLE kd_j AS
    SELECT novel_coordinates AS coordinates, novel_type AS type FROM main.novel
    UNION ALL
    SELECT ref_coordinates, type FROM main."other"'))
  invisible(DBI::dbExecute(con, "CREATE INDEX kd_j_idx ON kd_j(coordinates)"))
  n_kd <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM kd_j")$n
  step(sprintf("KD database: %s non-annotated junctions.", format(n_kd, big.mark = ",")))

  ## ---- Which KD junctions appear in each other database ------------------------
  ## Each other database's own junction tables are scanned once, probing the KD
  ## list's index - this does not depend on how (or whether) the other
  ## databases' tables are indexed. CROSS JOIN pins that loop order.
  short <- sub("_.*$", "", sub("\\.sqlite$", "", other_dbs))
  for (i in seq_along(other_dbs)) {
    step(sprintf("Other database %d/%d (%s)...", i, length(other_dbs), other_dbs[i]))
    invisible(DBI::dbExecute(con, sprintf('
      CREATE TEMP TABLE hit_%1$d AS
      SELECT j.coordinates
      FROM %2$s.novel n CROSS JOIN kd_j j ON j.coordinates = n.novel_coordinates
      UNION
      SELECT j.coordinates
      FROM %2$s."other" o CROSS JOIN kd_j j ON j.coordinates = o.ref_coordinates', i, aliases[i])))
    invisible(DBI::dbExecute(con, sprintf("CREATE INDEX hit_%1$d_idx ON hit_%1$d(coordinates)", i)))
  }

  step("Combining results...")
  idx <- seq_along(other_dbs)
  final <- DBI::dbGetQuery(con, sprintf("
    SELECT j.coordinates, j.type, %s
    FROM kd_j j
    %s",
    paste0("(h", idx, ".coordinates IS NOT NULL) AS in_", short, collapse = ", "),
    paste0("LEFT JOIN hit_", idx, " h", idx, " ON h", idx, ".coordinates = j.coordinates",
          collapse = "\n  ")))
  DBI::dbDisconnect(con)

  in_cols <- paste0("in_", short)
  final$n_other_dbs <- rowSums(final[, in_cols, drop = FALSE])

  ## ---- Summary -----------------------------------------------------------------
  k <- length(other_dbs)
  summarise_unique <- function(d, label) {
    row <- data.frame(type = label, n_junctions = nrow(d),
                      unique_to_KD_data = sum(d$n_other_dbs == 0),
                      pct_unique = round(100 * mean(d$n_other_dbs == 0), 3))
    for (m in 1:k) row[[paste0("in_", m, "_other_dbs")]] <- sum(d$n_other_dbs == m)
    row
  }
  summary_tbl <- do.call(rbind, c(
    lapply(split(final, final$type), function(d) summarise_unique(d, d$type[1])),
    list(summarise_unique(final, "ALL"))))
  rownames(summary_tbl) <- NULL

  per_db <- data.frame(database = short,
                      n_KD_junctions_found = colSums(final[, in_cols, drop = FALSE]),
                      pct_of_KD_junctions = round(100 * colMeans(final[, in_cols, drop = FALSE]), 2))
  rownames(per_db) <- NULL

  cat("\nKD database: ", kd_db, "  |  other databases: ", paste(other_dbs, collapse = ", "), "\n", sep = "")
  cat("\n=== KD junctions by number of OTHER databases that also contain them ===\n")
  print(summary_tbl)
  cat("\n=== Share of KD junctions found in each other database ===\n")
  print(per_db)

  step(sprintf("Finished in %.1f min.", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  saveRDS(list(summary = summary_tbl, per_database = per_db, per_junction = final), out_rds)
  cat("Per-junction table and summaries saved to", out_rds, "\n")




}



## ============================================================================
## kd_unique_junctions_arm_check.R
##
## Follow-up to kd_junction_uniqueness_check_fast.R.
##
## That script found junctions present in the KD database and in NO other
## IntroVerse database. But "KD database" pools the KD ('case') and control
## clusters, so those junctions are not necessarily knockdown products: a
## junction seen in both arms is a cell-line / depth effect. This script
## splits them, per RBP project (each RBP project holds its own KD and
## matched control cluster), into:
##   KD only       - detected in the KD cluster, absent from the control cluster
##   control only  - detected in the control cluster, absent from the KD cluster
##   both          - detected in both clusters of that project
##
## HOW IT STAYS CHEAP
## Only the junctions found in no other database are looked up, in ONE
## filtered pass over local_splicing_activity (no coordinate index exists on
## that table, so one sequential read of the file is unavoidable). The pass
## keeps just the matching rows - there is no large sort or GROUP BY over the
## whole table, which is what made the first approach slow.
##
## LIMITATION - read this before interpreting "KD only"
## "Absent from the control cluster" can simply mean the control was sequenced
## too shallowly to see the junction. This script does NOT separate those two
## explanations; it reports how well the KD-arm detection is supported
## (individuals, reads) as a first indication only.
## ============================================================================
if (run) {

  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
  })

  ## ---- Parameters (edit these) ----------------------------------------------
  database_folder <- "/home/sg2173/PROJECTS/srRNAseq/introverse_v2/database"   # same as config.R
  ensembl_version <- "115"
  kd_db           <- "ENCODE_SR_1read_shRNA.sqlite"

  ## Output of kd_junction_uniqueness_check_fast.R
  fast_rds <- "kd_junction_uniqueness_fast.rds"
  out_rds  <- "kd_unique_junctions_arm_check.rds"
  ## -----------------------------------------------------------------------------

  db_path <- function(db_file) {
    file.path(database_folder, sub("\\.sqlite$", "", db_file), ensembl_version, db_file)
  }
  step <- function(msg) message(format(Sys.time(), "%H:%M:%S"), "  ", msg)
  t_start <- Sys.time()

  ## ---- Inputs ------------------------------------------------------------------
  if (!file.exists(fast_rds)) stop("Cannot find ", fast_rds, " - run kd_junction_uniqueness_check_fast.R first (or fix the path).")
  if (!file.exists(db_path(kd_db))) stop("Database file does not exist: ", db_path(kd_db))
  if (file.size(db_path(kd_db)) == 0) stop("Database file is empty (0 bytes): ", db_path(kd_db))

  per_junction <- readRDS(fast_rds)$per_junction
  stopifnot(all(c("coordinates", "type", "n_other_dbs") %in% names(per_junction)))
  u <- unique(per_junction[per_junction$n_other_dbs == 0, c("coordinates", "type")])
  step(sprintf("%s junctions found in no other database (from %s).", format(nrow(u), big.mark = ","), fast_rds))
  if (nrow(u) == 0) stop("No junctions to check.")

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path(kd_db))

  required <- c("ref_coordinates", "novel_coordinates", "project", "cluster",
                "ref_n_individuals", "ref_sum_counts", "novel_n_individuals", "novel_sum_counts")
  tables <- DBI::dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type = 'table'")$name
  if (!("local_splicing_activity" %in% tables)) {
    stop(kd_db, " has no table 'local_splicing_activity'. Tables found: ", paste(tables, collapse = ", "))
  }
  cols <- DBI::dbGetQuery(con, "PRAGMA table_info('local_splicing_activity')")$name
  if (length(setdiff(required, cols)) > 0) {
    stop("local_splicing_activity lacks column(s): ", paste(setdiff(required, cols), collapse = ", "),
        ". Columns found: ", paste(cols, collapse = ", "))
  }

  ## ---- The one pass ------------------------------------------------------------
  DBI::dbWriteTable(con, "u", u, temporary = TRUE, overwrite = TRUE)
  invisible(DBI::dbExecute(con, "CREATE INDEX u_idx ON u(coordinates)"))

  step("Scanning local_splicing_activity ONCE for those junctions (the slow step - reads the whole file)...")
  ## A row can carry the junction as its reference junction (unannotated /
  ## novel-combination junctions) or as its novel junction (novel donor/acceptor),
  ## so either coordinate column is tested. Only matching rows are kept.
  invisible(DBI::dbExecute(con, "
    CREATE TEMP TABLE hits_raw AS
    SELECT ref_coordinates, novel_coordinates, project, cluster,
          ref_n_individuals, ref_sum_counts, novel_n_individuals, novel_sum_counts
    FROM main.local_splicing_activity
    WHERE ref_coordinates   IN (SELECT coordinates FROM u)
      OR novel_coordinates IN (SELECT coordinates FROM u)"))
  n_raw <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM hits_raw")$n
  step(sprintf("Scan done: %s matching rows kept.", format(n_raw, big.mark = ",")))

  ## Stack into one row per (junction, project, cluster): the junction's own
  ## individuals/reads come from the ref_* columns when it matched as the
  ## reference junction, from the novel_* columns when it matched as the novel one.
  invisible(DBI::dbExecute(con, "
    CREATE TEMP TABLE jx_rows AS
    SELECT ref_coordinates AS coordinates, project, cluster,
          ref_n_individuals AS n_ind, ref_sum_counts AS n_counts
    FROM hits_raw WHERE ref_coordinates IN (SELECT coordinates FROM u)
    UNION ALL
    SELECT novel_coordinates, project, cluster, novel_n_individuals, novel_sum_counts
    FROM hits_raw WHERE novel_coordinates IS NOT NULL
                    AND novel_coordinates IN (SELECT coordinates FROM u)"))

  cluster_labels <- DBI::dbGetQuery(con, "SELECT cluster, COUNT(*) AS n FROM jx_rows GROUP BY cluster")
  if (!all(cluster_labels$cluster %in% c("case", "control"))) {
    warning("Cluster labels other than 'case'/'control' were found and will be ignored: ",
            paste(setdiff(cluster_labels$cluster, c("case", "control")), collapse = ", "))
  }

  pairs <- DBI::dbGetQuery(con, "
    SELECT coordinates, project,
          SUM(CASE WHEN cluster = 'case'    THEN 1 ELSE 0 END) AS rows_case,
          SUM(CASE WHEN cluster = 'control' THEN 1 ELSE 0 END) AS rows_control,
          MAX(CASE WHEN cluster = 'case'    THEN n_ind    END)  AS ind_case,
          MAX(CASE WHEN cluster = 'control' THEN n_ind    END)  AS ind_control,
          MAX(CASE WHEN cluster = 'case'    THEN n_counts END)  AS reads_case,
          MAX(CASE WHEN cluster = 'control' THEN n_counts END)  AS reads_control
    FROM jx_rows
    GROUP BY coordinates, project")
  DBI::dbDisconnect(con)

  ## ---- Classify ------------------------------------------------------------------
  pairs <- merge(pairs, u, by = "coordinates", all.x = TRUE)
  pairs$status <- ifelse(pairs$rows_case > 0 & pairs$rows_control > 0, "both",
                  ifelse(pairs$rows_case > 0, "KD only",
                  ifelse(pairs$rows_control > 0, "control only", NA)))
  pairs <- pairs[!is.na(pairs$status), ]

  ## Junction level (across all projects) - mutually exclusive classes
  flags <- rowsum(cbind(any_case    = as.numeric(pairs$rows_case > 0),
                        any_control = as.numeric(pairs$rows_control > 0),
                        kd_only_in_a_project = as.numeric(pairs$status == "KD only")),
                  group = pairs$coordinates)
  jx <- data.frame(coordinates = rownames(flags), flags > 0, row.names = NULL)
  jx <- merge(u, jx, by = "coordinates", all.x = TRUE)
  jx$found_in_lsa <- !is.na(jx$any_case)
  for (v in c("any_case", "any_control", "kd_only_in_a_project")) jx[[v]][is.na(jx[[v]])] <- FALSE

  summ_row <- function(d, label) {
    data.frame(type = label,
              n_junctions            = nrow(d),
              not_found_in_lsa       = sum(!d$found_in_lsa),
              kd_only_everywhere     = sum(d$any_case & !d$any_control),
              control_only_everywhere = sum(!d$any_case & d$any_control),
              both_arms_somewhere    = sum(d$any_case & d$any_control),
              kd_only_in_a_project   = sum(d$kd_only_in_a_project))
  }
  table1 <- do.call(rbind, c(lapply(split(jx, jx$type), function(d) summ_row(d, d$type[1])),
                            list(summ_row(jx, "ALL"))))
  rownames(table1) <- NULL

  ## (junction, project) level - what each RBP project's matched control says
  st <- c("KD only", "both", "control only")
  table2 <- do.call(rbind, lapply(split(pairs, pairs$type), function(d) {
    tb <- table(factor(d$status, levels = st))
    data.frame(type = d$type[1], `KD only` = tb[["KD only"]], both = tb[["both"]],
              `control only` = tb[["control only"]], check.names = FALSE)
  }))
  tb_all <- table(factor(pairs$status, levels = st))
  table2 <- rbind(table2, data.frame(type = "ALL", `KD only` = tb_all[["KD only"]], both = tb_all[["both"]],
                                    `control only` = tb_all[["control only"]], check.names = FALSE))
  rownames(table2) <- NULL

  ## Per RBP project
  table3 <- do.call(rbind, lapply(split(pairs, pairs$project), function(d)
    data.frame(project = d$project[1],
              `KD only` = sum(d$status == "KD only"), both = sum(d$status == "both"),
              `control only` = sum(d$status == "control only"), check.names = FALSE)))
  table3 <- table3[order(-table3[["KD only"]]), ]
  rownames(table3) <- NULL

  ## How well supported are the KD-only detections? (first indication only)
  kd_only <- pairs[pairs$status == "KD only", ]
  support_row <- function(d, label) {
    data.frame(type = label, n_kd_only_pairs = nrow(d),
              pct_single_individual = round(100 * mean(d$ind_case == 1, na.rm = TRUE), 1),
              median_individuals    = median(d$ind_case, na.rm = TRUE),
              median_reads          = median(d$reads_case, na.rm = TRUE))
  }
  table4 <- if (nrow(kd_only) > 0) {
    rbind(do.call(rbind, lapply(split(kd_only, kd_only$type), function(d) support_row(d, d$type[1]))),
          support_row(kd_only, "ALL"))
  } else NULL
  if (!is.null(table4)) rownames(table4) <- NULL

  ## ---- Print ------------------------------------------------------------------------
  cat("\nJunctions found in no other database, split by arm. 'Project' = one RBP's KD + matched control.\n")
  cat("\n=== 1. Per junction, across all projects ===\n",
      "kd_only_everywhere: in a KD cluster, in NO control cluster anywhere\n",
      "control_only_everywhere: in a control cluster, in no KD cluster anywhere\n",
      "both_arms_somewhere: in at least one KD and at least one control cluster (possibly different projects)\n",
      "kd_only_in_a_project: has >=1 project where the KD cluster has it and that project's own control does not\n", sep = "")
  print(table1)
  cat("\n=== 2. Per (junction, RBP project) pair ===\n")
  print(table2)
  cat("\n=== 3. Per RBP project (top 25 by 'KD only') ===\n")
  print(head(table3, 25))
  cat("\n=== 4. Support for the 'KD only' pairs (KD-cluster individuals / reads for that junction) ===\n")
  print(table4)

  step(sprintf("Finished in %.1f min.", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  saveRDS(list(per_junction = jx, per_pair = pairs, by_type = table1, by_pair_status = table2,
              by_project = table3, kd_only_support = table4), out_rds)
  cat("Per-junction and per-pair tables and summaries saved to", out_rds, "\n")
}


## ============================================================================
## kd_exclusive_junctions_ranked.R
##
## Follow-up to kd_unique_junctions_arm_check.R. Lists, per KD RBP project, the
## junctions found ONLY in that project's KD cluster, ranked so the best
## supported come first. It reads the .rds saved by the previous script, so
## the ranking itself needs no database access and runs in seconds. (An
## optional gene annotation step does read a few small indexed tables of the KD
## database - see `annotate_genes`.)
##
## Every junction handled here is already known (from the previous two scripts)
## to be absent from every other IntroVerse database.
##
## TWO DEFINITIONS - both are written out, flagged, so you can use either:
##   exclusive (strict): the junction is detected in ONE RBP project only, in
##       its KD cluster only - not in that project's control, and not in any
##       other RBP project (KD or control).
##   KD-only pair (looser): the junction is in a project's KD cluster but not in
##       that same project's control. It may ALSO occur in other RBP projects.
##
## RANKING (within each project): more samples first, then more reads.
##   samples = individuals carrying the junction in that KD cluster (ind_case)
##   reads   = total reads supporting it in that KD cluster (reads_case)
## Samples go first because detection in several samples is stronger evidence
## than one deeply covered sample; change `rank_samples_first` to rank by reads.
## ============================================================================

if (run) {
    


suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

## ---- Parameters (edit these) ----------------------------------------------
arm_rds <- "kd_unique_junctions_arm_check.rds"      # output of kd_unique_junctions_arm_check.R

rank_samples_first <- TRUE     # TRUE: samples then reads; FALSE: reads then samples

## Project(s) left out of the second background comparison (section 2). NULL =
## automatically the project with the most exclusive junctions, because a single
## outlier project can dominate the pooled counts.
exclude_projects <- NULL

## Optional: add gene name, parent intron and distance from the KD database's
## small master tables (novel / intron / other / gene). Skipped with a message
## if anything unexpected is found - the ranking does not depend on it.
annotate_genes  <- TRUE
database_folder <- "/home/sg2173/PROJECTS/srRNAseq/introverse_v2/database"   # same as config.R
ensembl_version <- "115"
kd_db           <- "ENCODE_SR_1read_shRNA.sqlite"

out_exclusive_csv <- "kd_exclusive_junctions_ranked.csv"
out_kd_only_csv   <- "kd_only_pairs_ranked.csv"
out_rds           <- "kd_exclusive_junctions_ranked.rds"
## -----------------------------------------------------------------------------

step <- function(msg) message(format(Sys.time(), "%H:%M:%S"), "  ", msg)

if (!file.exists(arm_rds)) stop("Cannot find ", arm_rds, " - run kd_unique_junctions_arm_check.R first (or fix the path).")
arm <- readRDS(arm_rds)
pairs <- arm$per_pair
need <- c("coordinates", "project", "type", "status", "ind_case", "reads_case")
if (!all(need %in% names(pairs))) stop("per_pair lacks column(s): ", paste(setdiff(need, names(pairs)), collapse = ", "))

step(sprintf("%s (junction, project) pairs read.", format(nrow(pairs), big.mark = ",")))

## ---- Which junctions are in how many projects --------------------------------
n_projects <- table(pairs$coordinates)                  # projects where it is detected, either arm
kd_only <- pairs[pairs$status == "KD only", c("coordinates", "project", "type", "ind_case", "reads_case")]
kd_only$n_projects_detected  <- as.integer(n_projects)[match(kd_only$coordinates, names(n_projects))]
kd_only$exclusive_to_project <- kd_only$n_projects_detected == 1
## Background set: all control-only (junction, project) pairs, i.e. the junction
## is in a project's control cluster but not in that same project's KD cluster.
## Such a junction cannot have been produced by a knockdown, so how many of them
## reach each support level is what chance alone gives. It is compared with ALL
## KD-only pairs, NOT with the "exclusive" subset: control samples appear to be
## shared across RBP projects, so the same control-only junction recurs in many
## projects' control clusters, and a "found in one project only" rule would strip
## the background (but barely touch KD, whose samples are project-specific).
ctrl_pairs <- pairs[pairs$status == "control only", c("coordinates", "project", "type", "ind_control", "reads_control")]
ctrl_pairs$n_projects_detected <- as.integer(n_projects)[match(ctrl_pairs$coordinates, names(n_projects))]

step(sprintf("%s KD-only pairs; %s of them are exclusive to a single project.",
             format(nrow(kd_only), big.mark = ","), format(sum(kd_only$exclusive_to_project), big.mark = ",")))
if (nrow(kd_only) == 0) stop("No KD-only pairs in ", arm_rds)

## ---- Ranking ------------------------------------------------------------------
## NAs (no count available) sort last.
rank_in_project <- function(d) {
  a <- if (rank_samples_first) -d$ind_case else -d$reads_case
  b <- if (rank_samples_first) -d$reads_case else -d$ind_case
  d <- d[order(d$project, a, b, d$coordinates, na.last = TRUE), ]
  d$rank_in_project <- ave(seq_len(nrow(d)), d$project, FUN = seq_along)
  d
}
kd_only   <- rank_in_project(kd_only)
exclusive <- kd_only[kd_only$exclusive_to_project, ]
exclusive <- rank_in_project(exclusive)                  # re-rank within the strict set
exclusive$global_rank <- order(order(-(if (rank_samples_first) exclusive$ind_case else exclusive$reads_case),
                                     -(if (rank_samples_first) exclusive$reads_case else exclusive$ind_case),
                                     exclusive$project, exclusive$coordinates, na.last = TRUE))

## ---- Optional gene annotation ----------------------------------------------------
add_annotation_cols <- function(d, ann) {
  m <- match(d$coordinates, ann$coordinates)
  d$gene_name     <- ann$gene_name[m]
  d$parent_intron <- ann$parent_intron[m]
  d$distance      <- ann$distance[m]
  d
}
ann <- NULL
if (annotate_genes) {
  step("Annotating junctions (gene, parent intron, distance) from the KD database's master tables...")
  ann <- tryCatch({
    path <- file.path(database_folder, sub("\\.sqlite$", "", kd_db), ensembl_version, kd_db)
    if (!file.exists(path) || file.size(path) == 0) stop("database file not found/empty: ", path)
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbWriteTable(con, "ex", data.frame(coordinates = unique(kd_only$coordinates)), temporary = TRUE, overwrite = TRUE)
    invisible(DBI::dbExecute(con, "CREATE INDEX ex_idx ON ex(coordinates)"))
    ## A coordinate lives in 'novel' (novel donor/acceptor) OR in 'other', never both,
    ## so both lookups are run for every junction instead of hard-coding type labels.
    nov <- DBI::dbGetQuery(con, "
      SELECT e.coordinates, g.gene_name, i.ref_coordinates AS parent_intron, n.distance AS distance
      FROM ex e
      JOIN main.novel n        ON n.novel_coordinates = e.coordinates
      LEFT JOIN main.intron i  ON i.ref_junID = n.ref_junID
      LEFT JOIN main.gene g    ON g.gene_id = i.gene_id")
    oth <- DBI::dbGetQuery(con, '
      SELECT e.coordinates, g.gene_name, NULL AS parent_intron, NULL AS distance
      FROM ex e
      JOIN main."other" o      ON o.ref_coordinates = e.coordinates
      LEFT JOIN main.gene g    ON g.gene_id = o.gene_id')
    oth$parent_intron <- as.character(oth$parent_intron); oth$distance <- as.numeric(oth$distance)
    all_ann <- rbind(nov, oth)
    ## a junction can map to several genes: collapse to one row per coordinate
    data.frame(
      coordinates   = sort(unique(all_ann$coordinates)),
      gene_name     = tapply(all_ann$gene_name, all_ann$coordinates,
                             function(x) paste(sort(unique(stats::na.omit(x))), collapse = ";"))[sort(unique(all_ann$coordinates))],
      parent_intron = tapply(all_ann$parent_intron, all_ann$coordinates, function(x) x[1])[sort(unique(all_ann$coordinates))],
      distance      = tapply(all_ann$distance, all_ann$coordinates, function(x) x[1])[sort(unique(all_ann$coordinates))],
      row.names = NULL, stringsAsFactors = FALSE)
  }, error = function(e) { message("  Annotation skipped: ", conditionMessage(e)); NULL })
  if (!is.null(ann)) {
    step(sprintf("Annotated %s of %s junctions.", format(sum(unique(kd_only$coordinates) %in% ann$coordinates), big.mark = ","),
                 format(length(unique(kd_only$coordinates)), big.mark = ",")))
  }
}
if (is.null(ann)) ann <- data.frame(coordinates = character(0), gene_name = character(0),
                                    parent_intron = character(0), distance = numeric(0))
kd_only   <- add_annotation_cols(kd_only, ann)
exclusive <- add_annotation_cols(exclusive, ann)

col_order <- c("project", "rank_in_project", "coordinates", "type", "ind_case", "reads_case",
               "n_projects_detected", "gene_name", "parent_intron", "distance")
exclusive <- exclusive[, c(col_order[1:2], "global_rank", col_order[3:length(col_order)])]
exclusive <- exclusive[order(exclusive$project, exclusive$rank_in_project), ]
kd_only   <- kd_only[order(kd_only$project, kd_only$rank_in_project), c(col_order, "exclusive_to_project")]
names(exclusive)[names(exclusive) == "ind_case"]   <- "samples_KD"
names(exclusive)[names(exclusive) == "reads_case"] <- "reads_KD"
names(kd_only)[names(kd_only) == "ind_case"]       <- "samples_KD"
names(kd_only)[names(kd_only) == "reads_case"]     <- "reads_KD"

## ---- Summaries -----------------------------------------------------------------------
## 1. Support grids: how many junctions reach each minimum samples x reads level.
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
max_s <- safe_max(c(exclusive$samples_KD, kd_only$samples_KD, ctrl_pairs$ind_control))
if (is.na(max_s)) max_s <- 1
s_thr <- seq_len(min(max_s, 5)); r_thr <- c(1, 2, 5, 10, 20, 50)
support_grid <- function(samples, reads) {
  g <- sapply(r_thr, function(r) sapply(s_thr, function(s) sum(samples >= s & reads >= r, na.rm = TRUE)))
  matrix(g, nrow = length(s_thr), dimnames = list(paste0(">= ", s_thr, " sample(s)"), paste0(">= ", r_thr, " reads")))
}
## KD / background in each cell. Inf = background is 0 while KD is not; NA = both 0.
ratio_grid <- function(kd, bg) {
  r <- round(kd / bg, 2)
  r[bg == 0 & kd > 0]  <- Inf
  r[bg == 0 & kd == 0] <- NA
  r
}
grid <- support_grid(exclusive$samples_KD, exclusive$reads_KD)          # exclusive KD junctions (for prioritising)

## Background comparison at (junction, project) pair level: all KD-only pairs vs
## all control-only pairs, with the same thresholds.
grid_kdp    <- support_grid(kd_only$samples_KD, kd_only$reads_KD)
grid_ctp    <- support_grid(ctrl_pairs$ind_control, ctrl_pairs$reads_control)
grid_pratio <- ratio_grid(grid_kdp, grid_ctp)

## Same again without the outlier project(s)
dropped <- if (is.null(exclude_projects)) names(which.max(table(exclusive$project))) else exclude_projects
kd2 <- kd_only[!(kd_only$project %in% dropped), ]
ct2 <- ctrl_pairs[!(ctrl_pairs$project %in% dropped), ]
grid_kdp2    <- support_grid(kd2$samples_KD, kd2$reads_KD)
grid_ctp2    <- support_grid(ct2$ind_control, ct2$reads_control)
grid_pratio2 <- ratio_grid(grid_kdp2, grid_ctp2)

## 2. Per project
ctrl_only <- table(factor(pairs$project[pairs$status == "control only"], levels = unique(pairs$project)))
kd_pairs  <- table(factor(kd_only$project, levels = unique(pairs$project)))
by_project <- do.call(rbind, lapply(split(exclusive, exclusive$project), function(d) {
  data.frame(project = d$project[1],
             n_exclusive        = nrow(d),
             n_ge2_samples      = sum(d$samples_KD >= 2, na.rm = TRUE),
             max_samples        = suppressWarnings(max(d$samples_KD, na.rm = TRUE)),
             max_reads          = suppressWarnings(max(d$reads_KD, na.rm = TRUE)),
             median_reads       = median(d$reads_KD, na.rm = TRUE))
}))
by_project$kd_only_pairs      <- as.integer(kd_pairs[by_project$project])
by_project$control_only_pairs <- as.integer(ctrl_only[by_project$project])
by_project$kd_to_control_only <- ifelse(by_project$control_only_pairs > 0,
                                        round(by_project$kd_only_pairs / by_project$control_only_pairs, 2), NA)
## Stricter tier, same project, same arm size: pairs seen in >= 2 samples
n_ge2 <- function(x, proj) { v <- tapply(x >= 2, proj, sum, na.rm = TRUE); v[is.na(v)] <- 0; v }
kd_ge2 <- n_ge2(kd_only$samples_KD, kd_only$project)
ct_ge2 <- n_ge2(ctrl_pairs$ind_control, ctrl_pairs$project)
by_project$kd_pairs_ge2_samples      <- as.integer(ifelse(is.na(kd_ge2[by_project$project]), 0, kd_ge2[by_project$project]))
by_project$control_pairs_ge2_samples <- as.integer(ifelse(is.na(ct_ge2[by_project$project]), 0, ct_ge2[by_project$project]))
by_project$ratio_ge2_samples <- ifelse(by_project$control_pairs_ge2_samples > 0,
                                       round(by_project$kd_pairs_ge2_samples / by_project$control_pairs_ge2_samples, 2),
                                       ifelse(by_project$kd_pairs_ge2_samples > 0, Inf, NA))
by_project <- by_project[order(-by_project$n_ge2_samples, -by_project$n_exclusive), ]
rownames(by_project) <- NULL

## ---- Print -------------------------------------------------------------------------------
cat("\n'Exclusive' = found in ONE RBP project only, in its KD cluster only (and in no other IntroVerse database).\n")
cat(sprintf("Exclusive junctions: %s   |   all KD-only pairs: %s\n",
            format(nrow(exclusive), big.mark = ","), format(nrow(kd_only), big.mark = ",")))

cat("\n=== 1. Exclusive junctions meeting a minimum number of samples AND reads (in the KD cluster) ===\n")
print(grid)

cat("\n=== 2. Background check: do KD-only pairs exceed what the CONTROL arm gives by chance? ===\n",
    "Unit = (junction, RBP project) pair. KD-only: in the project's KD cluster, not in its control. Control-only:\n",
    "the mirror image - it cannot be produced by a knockdown, so it is the level chance alone reaches. ratio =\n",
    "KD-only / control-only at each tier: ~1 means no more than background, clearly above 1 means more than\n",
    "chance. Inf = background is 0; NA = both are 0. Assumes the arms have similar numbers of samples and depth.\n", sep = "")
share1 <- function(d) round(100 * mean(d$n_projects_detected == 1, na.rm = TRUE), 1)
cat(sprintf("Diagnostic - share of pairs whose junction occurs in only ONE project: KD-only %s%%, control-only %s%%.\n",
            share1(kd_only), share1(ctrl_pairs)),
    "(A much lower share for control-only fits control samples being shared across projects; this is why the\n",
    "background is NOT restricted to single-project junctions.)\n", sep = "")
cat(sprintf("\n-- All projects (KD-only pairs: %s, control-only pairs: %s) --\n",
            format(nrow(kd_only), big.mark = ","), format(nrow(ctrl_pairs), big.mark = ",")))
cat("\nKD-only pairs:\n");      print(grid_kdp)
cat("\nControl-only pairs (background):\n"); print(grid_ctp)
cat("\nRatio KD-only / control-only:\n"); print(grid_pratio)
cat(sprintf("\n-- Without %s (KD-only pairs: %s, control-only pairs: %s) --\n",
            paste(dropped, collapse = ", "), format(nrow(kd2), big.mark = ","), format(nrow(ct2), big.mark = ",")))
cat("\nKD-only pairs:\n");      print(grid_kdp2)
cat("\nControl-only pairs (background):\n"); print(grid_ctp2)
cat("\nRatio KD-only / control-only:\n"); print(grid_pratio2)

cat("\n=== 3. Per RBP project (ordered by number of exclusive junctions with >= 2 samples, then total) ===\n",
    "kd_to_control_only = KD-only pairs / control-only pairs in that project (each project compared with its own\n",
    "control). ratio_ge2_samples is the same ratio restricted to pairs seen in >= 2 samples - the tier where chance\n",
    "should rarely reach. A ratio near 1 means KD-only counts are what chance alone would give.\n", sep = "")
print(head(by_project, 25), row.names = FALSE)

cat("\n=== 4. Top 30 exclusive junctions overall ===\n")
print(head(exclusive[order(exclusive$global_rank), c("global_rank", "project", "coordinates", "type", "samples_KD", "reads_KD", "gene_name")], 30), row.names = FALSE)

write.csv(exclusive, out_exclusive_csv, row.names = FALSE)
write.csv(kd_only,   out_kd_only_csv,   row.names = FALSE)
saveRDS(list(exclusive = exclusive, kd_only_pairs = kd_only, control_only_pairs = ctrl_pairs,
             grid = grid, grid_kd_pairs = grid_kdp, grid_control_pairs = grid_ctp, grid_pair_ratio = grid_pratio,
             excluded_projects = dropped, grid_kd_pairs_without = grid_kdp2, grid_control_pairs_without = grid_ctp2,
             grid_pair_ratio_without = grid_pratio2, by_project = by_project), out_rds)
cat(sprintf("\nWritten: %s (exclusive, ranked within project), %s (all KD-only pairs, flagged), %s\n",
            out_exclusive_csv, out_kd_only_csv, out_rds))
}

## ============================================================================
## kd_arm_depth_check.R
##
## Reads ONLY the small `metadata` table of the ENCODE KD database (instant).
## Answers three questions that decide how to read the KD-only / control-only
## comparison:
##   1. How many samples (and cell lines) does each project's KD and control arm
##      really have?
##   2. Are control samples SHARED between RBP projects? (This explains why
##      control-only junctions recur across projects - the earlier numbers
##      suggested it, this checks it directly.)
##   3. Does the KD arm have more sequencing depth than its control, and does
##      that explain which projects have many more KD-only than control-only
##      junctions (e.g. AQR)?
##
## On (3): for junctions that are only sporadically detected, the chance of
## seeing one in an arm grows roughly in proportion to the arm's total reads.
## So, to first order, KD-only / control-only pairs is expected to be about
## KD depth / control depth by chance alone. `enrichment_over_depth` divides the
## observed pair ratio by that depth ratio. It is a ROUGH first-order guide, not
## a statistical test.
## ============================================================================

if (run) {
  
suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

## ---- Parameters (edit these) ----------------------------------------------
## Any build of the ENCODE KD database will do for the metadata (the samples are
## the same); point this at the one you want.
metadata_db <- "/home/sg2173/PROJECTS/srRNAseq/introverse_v2/database/ENCODE_SR_1read_shRNA/115/ENCODE_SR_1read_shRNA.sqlite"

## Output of kd_exclusive_junctions_ranked.R; if present, its per-project table
## is joined in (section 3). Set to NULL to skip.
ranked_rds <- "kd_exclusive_junctions_ranked.rds"

out_csv <- "kd_arm_depth_summary.csv"
## -----------------------------------------------------------------------------

if (!file.exists(metadata_db)) stop("Database file does not exist: ", metadata_db)
if (file.size(metadata_db) == 0) stop("Database file is empty (0 bytes): ", metadata_db)

con <- DBI::dbConnect(RSQLite::SQLite(), metadata_db)
cols <- DBI::dbGetQuery(con, "PRAGMA table_info(metadata)")$name
need <- c("SRA_project", "cluster", "cell_line", "sample_id", "read_depth")
if (!all(need %in% cols)) {
  stop("metadata lacks column(s): ", paste(setdiff(need, cols), collapse = ", "),
       ". Columns found: ", paste(cols, collapse = ", "))
}
md <- DBI::dbGetQuery(con, "SELECT SRA_project AS project, cluster, cell_line, sample_id, read_depth FROM metadata")
DBI::dbDisconnect(con)

md$read_depth <- suppressWarnings(as.numeric(md$read_depth))
if (!all(md$cluster %in% c("case", "control"))) {
  warning("cluster labels other than case/control present and ignored: ",
          paste(setdiff(unique(md$cluster), c("case", "control")), collapse = ", "))
  md <- md[md$cluster %in% c("case", "control"), ]
}

## ---- 1. Arm sizes and depth per project ----------------------------------------
arm <- function(d) data.frame(samples = length(unique(d$sample_id)),
                              cell_lines = length(unique(d$cell_line)),
                              depth_M = round(sum(d$read_depth[!duplicated(d$sample_id)], na.rm = TRUE) / 1e6, 1))
by_arm <- do.call(rbind, lapply(split(md, list(md$project, md$cluster), drop = TRUE), function(d)
  cbind(project = d$project[1], cluster = d$cluster[1], arm(d))))
rownames(by_arm) <- NULL
wide <- merge(by_arm[by_arm$cluster == "case", -2], by_arm[by_arm$cluster == "control", -2],
              by = "project", suffixes = c("_KD", "_control"), all = TRUE)
wide$depth_ratio_KD_to_control <- round(wide$depth_M_KD / wide$depth_M_control, 2)

cat("\n=== 1. Arm sizes and total sequencing depth (millions of reads) per RBP project ===\n")
cat(sprintf("Projects: %d. KD arm samples - min %s, median %s, max %s. Control arm samples - min %s, median %s, max %s.\n",
            nrow(wide),
            min(wide$samples_KD, na.rm = TRUE), median(wide$samples_KD, na.rm = TRUE), max(wide$samples_KD, na.rm = TRUE),
            min(wide$samples_control, na.rm = TRUE), median(wide$samples_control, na.rm = TRUE), max(wide$samples_control, na.rm = TRUE)))
unequal <- wide[which(wide$samples_KD != wide$samples_control | is.na(wide$samples_KD) | is.na(wide$samples_control)), ]
cat(sprintf("Projects whose two arms do NOT have the same number of samples: %d\n", nrow(unequal)))
if (nrow(unequal) > 0) print(head(unequal, 20), row.names = FALSE)
cat("\nDepth ratio KD/control across projects (quantiles):\n")
print(round(quantile(wide$depth_ratio_KD_to_control, c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE), 2))

## ---- 2. Are control samples shared between projects? ----------------------------------
ctl <- md[md$cluster == "control", ]
n_proj_per_ctl <- tapply(ctl$project, ctl$sample_id, function(x) length(unique(x)))
kd  <- md[md$cluster == "case", ]
n_proj_per_kd  <- tapply(kd$project, kd$sample_id, function(x) length(unique(x)))
cat("\n=== 2. Are samples shared between RBP projects? ===\n")
cat(sprintf("Control samples: %d distinct; used by more than one project: %d (%.1f%%); most projects sharing one control sample: %d (median %s).\n",
            length(n_proj_per_ctl), sum(n_proj_per_ctl > 1), 100 * mean(n_proj_per_ctl > 1),
            max(n_proj_per_ctl), median(n_proj_per_ctl)))
cat(sprintf("KD samples:      %d distinct; used by more than one project: %d (%.1f%%).\n",
            length(n_proj_per_kd), sum(n_proj_per_kd > 1), 100 * mean(n_proj_per_kd > 1)))

## ---- 3. Relate depth to the KD-only / control-only ratios ------------------------------------
joined <- NULL
if (!is.null(ranked_rds) && file.exists(ranked_rds)) {
  bp <- readRDS(ranked_rds)$by_project
  keep <- intersect(c("project", "kd_only_pairs", "control_only_pairs", "kd_to_control_only",
                      "kd_pairs_ge2_samples", "control_pairs_ge2_samples", "ratio_ge2_samples"), names(bp))
  joined <- merge(bp[, keep], wide[, c("project", "samples_KD", "samples_control", "depth_M_KD", "depth_M_control", "depth_ratio_KD_to_control")],
                  by = "project")
  joined$enrichment_over_depth <- round(joined$kd_to_control_only / joined$depth_ratio_KD_to_control, 2)
  joined <- joined[order(-joined$enrichment_over_depth), ]
  rownames(joined) <- NULL
  cat("\n=== 3. Observed KD-only/control-only pair ratio versus the depth ratio, per project ===\n",
      "enrichment_over_depth = kd_to_control_only / depth_ratio_KD_to_control. About 1 = what depth alone predicts\n",
      "(first-order guide only); well above 1 = more KD-only junctions than depth explains.\n", sep = "")
  print(head(joined[, c("project", "kd_only_pairs", "control_only_pairs", "kd_to_control_only",
                        "depth_ratio_KD_to_control", "enrichment_over_depth")], 25), row.names = FALSE)
  ok <- is.finite(joined$kd_to_control_only) & is.finite(joined$depth_ratio_KD_to_control) &
        joined$kd_to_control_only > 0 & joined$depth_ratio_KD_to_control > 0
  if (sum(ok) >= 5) {
    cat(sprintf("\nSpearman correlation, log pair ratio vs log depth ratio, %d projects: %.2f\n", sum(ok),
                suppressWarnings(cor(joined$kd_to_control_only[ok], joined$depth_ratio_KD_to_control[ok], method = "spearman"))))
  }
} else if (!is.null(ranked_rds)) {
  cat("\n(Section 3 skipped: ", ranked_rds, " not found - run kd_exclusive_junctions_ranked.R first, or fix the path.)\n", sep = "")
}

write.csv(if (!is.null(joined)) joined else wide, out_csv, row.names = FALSE)
cat("\nWritten:", out_csv, "\n")
}