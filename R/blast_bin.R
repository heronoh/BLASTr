# Which BLAST+ BLASTr runs. In this order:
# 1. the conda environment of BLASTr (`blastr-blast-env`, made by
#    install_dependencies() with BLAST+ `blast_version_pinned`);
# 2. the BLAST+ of the system (on the PATH), if it is `blast_version_min` or
#    newer;
# 3. otherwise an error that says to run install_dependencies().
# Searches never install or update BLAST+: only install_dependencies() does.

# BLAST+ version installed by install_dependencies(): the one BLASTr is tested
# with (see NEWS.md).
blast_version_pinned <- "2.17.0"

# Lowest BLAST+ version BLASTr runs with: older ones lack options it uses
# (BLAST+ 2.12 does not accept `-mt_mode 2`).
blast_version_min <- "2.16.0"

# Versions already asked to each program in this session, and the notices
# already given.
blastr_state <- new.env(parent = emptyenv())
blastr_state$versions <- list()
blastr_state$notified <- character(0L)

#' @title Is BLAST+ ready for BLASTr?
#'
#' @description
#' Checks, without installing or downloading anything, whether BLASTr has a
#' `blastn` it runs with (BLAST+ 2.16.0 or newer): the one of its conda
#' environment (installed by [install_dependencies()], BLAST+ 2.17.0) or, when
#' that environment has none, the one of the system. Useful to run code only
#' where BLAST+ is ready, as in the examples.
#'
#' @param env_name Name of the conda environment.
#'
#' @returns `TRUE` or `FALSE`.
#'
#' @examples
#' has_blast()
#' @family setup functions
#' @export
has_blast <- function(env_name = "blastr-blast-env") {
  !is.null(tryCatch(resolve_blast_bin("blastn", env_name), error = function(e) NULL))
}

# Directory of a condathis environment, as condathis::get_env_dir() gives it,
# but without creating the condathis directory when it does not exist.
conda_env_dir <- function(env_name) {
  if (!nzchar(Sys.getenv("XDG_DATA_HOME")) && identical(Sys.info()[["sysname"]], "Darwin")) {
    withr::local_envvar(XDG_DATA_HOME = fs::path_home(".local", "share"))
  }
  fs::path(tools::R_user_dir("condathis", which = "data"), "envs", env_name)
}

# "blastn: 2.17.0+ ..." -> "2.17.0"
parse_blast_version <- function(x) {
  v <- regmatches(x, regexpr("[0-9]+\\.[0-9]+\\.[0-9]+(?=\\+)", x, perl = TRUE))
  if (length(v) == 0L) NA_character_ else v[[1L]]
}

# Version of a BLAST+ program ("2.17.0"), or NA when it does not answer;
# asked once per file and modification time.
program_version <- function(bin) {
  key <- paste(bin, file.mtime(bin))
  v <- blastr_state$versions[[key]]
  if (!is.null(v)) {
    return(v)
  }
  out <- tryCatch(
    processx::run(
      bin, "-version",
      error_on_status = FALSE, timeout = 30,
      env = c("current", BLAST_USAGE_REPORT = "false")
    )$stdout,
    error = function(e) ""
  )
  v <- parse_blast_version(out)
  blastr_state$versions[[key]] <- v
  v
}

blast_version_ok <- function(v) {
  !is.na(v) && utils::compareVersion(v, blast_version_min) >= 0L
}

# Version of `blastn` in a conda environment, or NA when it is not there.
env_blast_version <- function(env_name = "blastr-blast-env") {
  bin <- fs::path(conda_env_dir(env_name), "bin", "blastn")
  if (!fs::file_exists(bin)) NA_character_ else program_version(bin)
}

# The BLAST+ program BLASTr runs: list(path, source = "conda" or "system",
# version), or an error saying what was found.
resolve_blast_bin <- function(cmd, env_name = "blastr-blast-env") {
  exe <- if (identical(.Platform$OS.type, "windows")) paste0(cmd, ".exe") else cmd
  env_bin <- fs::path(conda_env_dir(env_name), "bin", exe)
  env_v <- if (fs::file_exists(env_bin)) program_version(env_bin) else NA_character_
  if (blast_version_ok(env_v)) {
    return(list(path = as.character(env_bin), source = "conda", version = env_v))
  }
  sys_bin <- unname(Sys.which(cmd))
  sys_v <- if (nzchar(sys_bin)) program_version(sys_bin) else NA_character_
  if (blast_version_ok(sys_v)) {
    return(list(path = normalizePath(sys_bin, mustWork = FALSE), source = "system", version = sys_v))
  }
  state <- function(found, v) {
    if (!found) "not found" else if (is.na(v)) "found, version unknown" else paste("version", v)
  }
  cli::cli_abort(
    c(
      `x` = "BLAST+ {blast_version_min} or newer ({.code {cmd}}) was not found.",
      `i` = "Conda environment {.val {env_name}}: {state(fs::file_exists(env_bin), env_v)}.",
      `i` = "System: {state(nzchar(sys_bin), sys_v)}.",
      `>` = "Install BLAST+ {blast_version_pinned} with {.code install_dependencies()}."
    ),
    class = "blastr_error_blast_not_found"
  )
}

blast_bin_path <- function(cmd, env_name = "blastr-blast-env") {
  resolve_blast_bin(cmd, env_name)$path
}

# Run a BLAST+ program (blastn, makeblastdb, blastdbcmd, ...) with the BLAST+
# chosen above; returns the result of processx::run().
run_blast_program <- function(cmd, args, env_name = "blastr-blast-env",
                              verbose = "silent", error = c("cancel", "continue")) {
  error <- rlang::arg_match(error)
  bin <- blast_bin_path(cmd, env_name)
  if (verbose %in% c("cmd", "full")) {
    cli::cli_inform(c(`i` = "{.code {cmd} {paste(args, collapse = ' ')}}"))
  }
  res <- processx::run(
    bin, as.character(args),
    error_on_status = FALSE,
    echo = verbose %in% c("output", "full"),
    env = c("current", BLAST_USAGE_REPORT = "false")
  )
  if (identical(error, "cancel") && !identical(res$status, 0L)) {
    err <- utils::tail(strsplit(res$stderr, "\n", fixed = TRUE)[[1L]], 5L)
    err <- gsub("([{}])", "\\1\\1", err)
    cli::cli_abort(
      c(`x` = "{.code {cmd}} failed (status {res$status}).", stats::setNames(err, rep("!", length(err)))),
      class = "blastr_error_blast_program"
    )
  }
  res
}

# Install BLAST+ `blast_version_pinned` in the conda environment: create it,
# update it in place when it has another version, or recreate it (`force`).
install_blast_env <- function(env_name = "blastr-blast-env", verbose = "silent", force = FALSE) {
  pkg <- paste0("bioconda::blast==", blast_version_pinned)
  channels <- c("conda-forge", "bioconda")
  if (isTRUE(force) || isFALSE(condathis::env_exists(env_name))) {
    condathis::create_env(
      packages = pkg, channels = channels, env_name = env_name,
      verbose = verbose, overwrite = isTRUE(force)
    )
  } else {
    current <- env_blast_version(env_name)
    if (!identical(current, blast_version_pinned)) {
      cli::cli_inform(c(
        `i` = "Updating BLAST+ in the conda environment {.val {env_name}} from {if (is.na(current)) 'an unknown version' else current} to {blast_version_pinned}."
      ))
      condathis::install_packages(packages = pkg, env_name = env_name, channels = channels, verbose = verbose)
    }
  }
  blastr_state$versions <- list()
  installed <- env_blast_version(env_name)
  if (!identical(installed, blast_version_pinned)) {
    cli::cli_warn(c(
      `!` = "BLAST+ in {.val {env_name}} is {if (is.na(installed)) 'not there' else installed}, not {blast_version_pinned}.",
      `i` = "Recreate it with {.code install_dependencies(force = TRUE)}."
    ))
  }
  invisible(TRUE)
}
