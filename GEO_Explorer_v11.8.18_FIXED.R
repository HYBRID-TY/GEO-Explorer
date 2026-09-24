############################################################
# v11.8.18 conservative correctness fix (2026-09-24):
# - strengthens raw-count versus normalized/log-scale detection;
# - uses source-name evidence without trusting integer-like values alone;
# - keeps missing normalized/log expression as missing (count merges still use zero);
# - validates sample/metadata alignment before differential analysis;
# - fixes the legacy single-sample merger's uninitialized `failed` counter.
# Original UI, outputs, and user-selectable analysis methods are preserved.
############################################################

############################################################
# GEO Explorer - BEST UNIVERSAL VERSION v9.1 BULK ONLY
# v14 FULL UI preserved + manual DEG + metadata explorer
# v7.4: Fixed vectorized processed-file detection (& not &&); reads processed XLSX before RAW.tar
# Bulk-only patch v6.3:
# - Excludes single-cell/10X matrix.mtx files from bulk parsing.
# - Detects mixed bulk + scRNA supplementary data and keeps only bulk-compatible files.
#
# Fix in v6.4:
# - Entrez conversion is now performed inside the ultra-fast reader.
# - Numeric gene_id columns like 54757 / 714 / 79887 are converted
#   immediately to Gene Symbols before DEG.
# - Message reports Gene ID type and conversion count.
############################################################

############################################################
# GEO Gene Expression Explorer - Auto Differential Analysis
# Upgraded version: auto-detect raw counts / FPKM-TPM / log expression
# raw count -> DESeq2 or limma-voom; normalized/log/microarray -> limma
############################################################

# =========================
# 0. Packages
# =========================
cran_pkgs <- c(
  "shiny", "GEOquery", "data.table", "dplyr", "tibble",
  "ggplot2", "ggpubr", "R.utils", "readxl", "Biobase", "pheatmap", "colourpicker", "RColorBrewer", "svglite"
)

bioc_pkgs <- c("DESeq2", "limma", "edgeR", "org.Hs.eg.db", "AnnotationDbi")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(shiny)
library(GEOquery)
library(data.table)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(R.utils)
library(readxl)
library(Biobase)
library(DESeq2)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(colourpicker)
library(RColorBrewer)
library(svglite)

select <- dplyr::select

cat("RUNNING VERSION = v11.8.18 conservative expression-type, NA, and sample-alignment fixes\n")
# NOTE: v8.3 keeps v7.9 DEG/plotting modules, but adds method-aware data source selection.
# Auto/Force DESeq2 prefer raw count matrices; Force limma re-scans online/local sources and prefers TPM/FPKM/normalized matrices.
# v8.4: show the actual DEG source file and synchronize dashboard/group counts with the samples actually used by DEG.


# =========================
# 0.5. Cache settings
# =========================
# Parsed cache avoids repeatedly scanning and testing supplementary files.
# It stores the successfully parsed expression matrix + metadata per GSE.
GEO_ROOT <- "GEO_downloads"
PARSED_CACHE_ROOT <- file.path(GEO_ROOT, "_parsed_cache")
dir.create(GEO_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(PARSED_CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)


# =========================
# 0.6. Robust download settings
# =========================
# Large GEO supplementary matrices can be hundreds of MB. On Windows/RStudio,
# download.file may fall back to deprecated wininet and emit warnings that can
# break retry logic. This downloader avoids wininet, uses binary mode, keeps a
# long timeout, removes partial files, and verifies .gz integrity when possible.
options(timeout = max(3600, getOption("timeout", 60)))

get_remote_file_size <- function(url) {
  out <- tryCatch({
    h <- utils::capture.output(utils::download.file(url, tempfile(), method = "libcurl", mode = "wb", quiet = TRUE, headers = TRUE))
    m <- regmatches(h, regexpr("(?i)content-length:\\s*[0-9]+", h, perl = TRUE))
    m <- unlist(m)
    if (length(m) == 0) return(NA_real_)
    as.numeric(gsub("[^0-9]", "", m[length(m)]))
  }, error = function(e) NA_real_, warning = function(w) NA_real_)
  out
}

is_gz_file_ok <- function(f) {
  if (!file.exists(f) || file.info(f)$size <= 0) return(FALSE)
  if (!grepl("\\.gz$", f, ignore.case = TRUE)) return(TRUE)
  ok <- tryCatch({
    con <- gzfile(f, open = "rb")
    on.exit(close(con), add = TRUE)
    readBin(con, what = "raw", n = 1)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  ok
}

safe_download_geo_file <- function(url, dest, max_retries = 3) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  # Avoid the deprecated/fragile Windows wininet method. libcurl is bundled with
  # most modern R builds and is much better for large HTTPS GEO files.
  method_order <- c()
  if (isTRUE(capabilities("libcurl"))) method_order <- c(method_order, "libcurl")
  method_order <- c(method_order, "auto")
  method_order <- unique(method_order)

  # If a previous timeout left a partial or corrupt file, remove it before retrying.
  if (file.exists(dest) && !is_gz_file_ok(dest)) {
    message("发现残缺/损坏文件，删除后重新下载: ", basename(dest))
    unlink(dest, force = TRUE)
  }

  for (attempt in seq_len(max_retries)) {
    if (file.exists(dest) && is_gz_file_ok(dest)) {
      message("本地已存在，跳过下载: ", basename(dest), " | ", round(file.info(dest)$size / 1024 / 1024, 1), " MB")
      return(TRUE)
    }

    if (file.exists(dest)) unlink(dest, force = TRUE)

    for (method in method_order) {
      message("下载文件 attempt ", attempt, "/", max_retries, " [", method, "]: ", basename(dest))

      ok <- tryCatch({
        suppressWarnings(utils::download.file(
          url,
          destfile = dest,
          mode = "wb",
          quiet = FALSE,
          method = method
        ))
        TRUE
      }, error = function(e) {
        message("下载失败 [", method, "]: ", conditionMessage(e))
        FALSE
      })

      local_ok <- file.exists(dest) && file.info(dest)$size > 0 && is_gz_file_ok(dest)

      if (ok && local_ok) {
        message("下载成功: ", basename(dest), " | ", round(file.info(dest)$size / 1024 / 1024, 1), " MB")
        return(TRUE)
      }

      if (file.exists(dest)) {
        message("下载未完成或文件不可读，删除残缺文件: ", basename(dest), " | ", round(file.info(dest)$size / 1024 / 1024, 1), " MB")
        unlink(dest, force = TRUE)
      }
    }

    Sys.sleep(min(8 * attempt, 20))
  }

  stop("download failed after retries: ", url,
       "\n建议手动浏览器下载该文件到: ", normalizePath(dirname(dest), winslash = "/", mustWork = FALSE),
       "\n文件名必须保持为: ", basename(dest))
}

cache_file_for_gse <- function(gse_id) {
  file.path(PARSED_CACHE_ROOT, paste0(toupper(trimws(gse_id)), "_v11_8_8_count_merge_rewrite_parsed_expr_meta.rds"))
}

# =========================
# 1. Basic helper functions
# =========================
parse_gene_assignment_symbol <- function(x) {
  # Robust Affymetrix/GEO gene_assignment parser.
  # Common formats include:
  #   "NM_001... // TP53 // tumor protein p53 // ..."
  #   "7892501 // DDX11L1 // DEAD/H-box helicase..."
  #   "TP53 /// TP53P1"
  # The old parser always took the 2nd field. That is often correct, but can
  # fail on some GPLs. This version scores all fields and keeps the most
  # gene-symbol-like token.
  x <- as.character(x)
  x[x == "" | x == "---" | is.na(x)] <- NA_character_

  get_one <- function(s) {
    if (is.na(s) || trimws(s) == "") return(NA_character_)

    # Split on common GEO delimiters.
    parts <- unlist(strsplit(s, "\\s*///\\s*|\\s*//\\s*|;|,"))
    parts <- trimws(parts)
    parts <- parts[!is.na(parts) & parts != "" & parts != "---" & parts != "NA"]

    if (length(parts) == 0) return(NA_character_)

    # Remove descriptions in parentheses and obvious non-symbol tokens.
    parts2 <- gsub("\\s*\\(.*?\\)\\s*", "", parts)
    parts2 <- trimws(parts2)

    # Strong HGNC-like symbols: start with a letter, short, no spaces.
    strong <- parts2[
      grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", parts2) &
        !grepl("^NM_|^NR_|^XM_|^XR_|^ENSG|^ENST|^[0-9]+$", parts2, ignore.case = TRUE) &
        !grepl("hypothetical|predicted|protein|transcript|cluster|orf ", parts2, ignore.case = TRUE)
    ]

    if (length(strong) > 0) return(strong[1])

    # Fallback: if GEO's 2nd field exists and is not numeric/accession-like.
    if (length(parts2) >= 2 &&
        grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", parts2[2]) &&
        !grepl("^NM_|^NR_|^XM_|^XR_|^ENSG|^ENST|^[0-9]+$", parts2[2], ignore.case = TRUE)) {
      return(parts2[2])
    }

    # Last fallback: first cleaned token.
    parts2[1]
  }

  symbol <- vapply(x, get_one, character(1))
  symbol <- clean_symbol_value(symbol)
  symbol
}

clean_group_value <- function(x) {
  x <- as.character(x)
  x <- gsub("^.*?:", "", x)
  trimws(x)
}

clean_symbol_value <- function(x) {
  x <- as.character(x)
  x <- gsub("\\..*$", "", x)
  x <- gsub("///.*$", "", x)
  x <- gsub("//.*$", "", x)
  x <- gsub(";.*$", "", x)
  x <- gsub(",.*$", "", x)
  x <- trimws(x)
  x[x == "" | x == "---" | x == "NA" | is.na(x)] <- NA
  x
}

ensembl_to_symbol <- function(ids) {
  ids2 <- clean_symbol_value(ids)
  is_ens <- grepl("^ENSG", ids2, ignore.case = TRUE)
  out <- ids2
  if (any(is_ens, na.rm = TRUE)) {
    mapped <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = unique(ids2[is_ens]),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    ))
    out[is_ens] <- as.character(mapped[ids2[is_ens]])
  }
  clean_symbol_value(out)
}

is_numeric_like <- function(x, min_prop = 0.5) {
  y <- suppressWarnings(as.numeric(x))
  mean(!is.na(y)) >= min_prop
}

# =========================
# Safe column/index helpers
# =========================
safe_char_cols <- function(x, all_cols = NULL) {
  # GEO supplementary files sometimes make nested/list-like candidate columns
  # after Excel/header repair. Base R data-frame subsetting cannot use list
  # indices, so always flatten to a plain character vector first.
  x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]
  x <- trimws(x)
  if (!is.null(all_cols)) {
    all_cols <- as.character(all_cols)
    x <- intersect(x, all_cols)
  }
  unique(x)
}

safe_one_col <- function(x, all_cols = NULL, fallback = NA_character_) {
  x <- safe_char_cols(x, all_cols = all_cols)
  if (length(x) == 0) return(fallback)
  x[1]
}

# v11.7.8: Encoding-safe string helpers.
# Some GEO supplementary Excel/text files contain Latin-1/Windows encoded text or
# binary workbook bytes. Base string functions such as grepl/tolower/make.names
# can throw "invalid multibyte string". These helpers convert to UTF-8 and drop
# undecodable bytes rather than stopping the whole GEO load.
safe_utf8 <- function(x) {
  x <- as.character(x)
  y <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = ""))
  y[is.na(y)] <- ""
  y
}

safe_make_names <- function(x, unique = TRUE) {
  make.names(safe_utf8(x), unique = unique)
}

safe_read_preview_lines <- function(f, n = 30) {
  if (grepl("\\.xlsx$|\\.xls$|\\.zip$|\\.tar$|\\.tar\\.gz$|\\.tgz$", basename(f), ignore.case = TRUE)) {
    return(character())
  }
  out <- tryCatch(readLines(f, n = n, warn = FALSE, encoding = "UTF-8"), error = function(e) NULL)
  if (is.null(out)) {
    out <- tryCatch(readLines(f, n = n, warn = FALSE, encoding = "latin1"), error = function(e) character())
  }
  safe_utf8(out)
}


# =========================
# Stage diagnosis / structured logging helpers
# =========================
.short_dim_msg <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.data.frame(x) || is.matrix(x)) return(paste0(nrow(x), " rows x ", ncol(x), " cols"))
  if (is.list(x) && !is.null(x$raw)) return(paste0("raw=", nrow(x$raw), " genes x ", max(0, ncol(x$raw) - 1), " samples"))
  paste0(class(x)[1], " length=", length(x))
}

safe_stage <- function(stage, expr, file = NULL, fatal = TRUE) {
  label <- if (!is.null(file)) paste0(stage, " | ", basename(file)) else stage
  message("
========== ", label, " ==========")
  t0 <- Sys.time()
  out <- tryCatch({
    val <- eval.parent(substitute(expr))
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    message("✅ Stage OK: ", label, " | ", .short_dim_msg(val), " | ", elapsed, " sec")
    val
  }, error = function(e) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    msg <- paste0(
      "❌ Stage FAILED: ", label, "
",
      "Elapsed: ", elapsed, " sec
",
      "Error: ", conditionMessage(e)
    )
    message(msg)
    if (isTRUE(fatal)) stop(msg, call. = FALSE)
    return(NULL)
  })
  out
}

safe_parser_try <- function(parser_name, fun, file) {
  safe_stage(
    stage = paste0("Parser: ", parser_name),
    file = file,
    expr = {
      res <- fun(file)
      if (is.null(res)) {
        message("Parser returned NULL: ", parser_name, " | ", basename(file))
      } else {
        message("Parser selected: ", parser_name, " | ", basename(file))
        if (!is.null(res$raw)) message("Parsed matrix: ", nrow(res$raw), " genes x ", ncol(res$raw) - 1, " samples")
        if (!is.null(res$symbol_col)) message("Symbol column: ", res$symbol_col)
        if (!is.null(res$gene_id_type)) message("Gene ID type: ", res$gene_id_type)
      }
      res
    },
    fatal = FALSE
  )
}

# v11.8.0 generic parser-reject logger.
# This is deliberately not GSE-specific. Every parser can now explain why it
# rejected a candidate file instead of silently returning NULL.
log_parser_reject <- function(parser, file, reason, detail = "") {
  msg <- paste0(
    "❌ Parser reject | ", parser,
    " | ", basename(file),
    " | ", reason,
    ifelse(is.null(detail) || detail == "", "", paste0(" | ", detail))
  )
  message(msg)
  invisible(NULL)
}

excel_sheet_preview_msg <- function(f, max_rows = 8, max_cols = 8) {
  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) character())
  if (length(sheets) == 0) return("Excel sheet preview unavailable")
  out <- character()
  for (sh in head(sheets, 5)) {
    dat <- tryCatch(as.data.frame(readxl::read_excel(f, sheet = sh, col_names = FALSE, n_max = max_rows), stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(dat)) {
      out <- c(out, paste0("Sheet ", sh, ": read failed"))
    } else {
      nr <- nrow(dat); nc <- ncol(dat)
      small <- dat[seq_len(min(nr, max_rows)), seq_len(min(nc, max_cols)), drop = FALSE]
      small[] <- lapply(small, safe_utf8)
      out <- c(out, paste0("Sheet ", sh, ": ", nr, " preview rows x ", nc, " cols; first row = ", paste(as.character(unlist(small[1,], use.names = FALSE)), collapse = " | ")))
    }
  }
  paste(out, collapse = "\n")
}

# =========================
# 2. Data type detection
# =========================
detect_expression_type <- function(mat_or_df, source_hint = NULL) {
  x <- as.matrix(mat_or_df)
  suppressWarnings(storage.mode(x) <- "numeric")
  vals <- as.numeric(x)
  vals <- vals[is.finite(vals) & !is.na(vals)]
  if (length(vals) < 100) {
    return(list(
      type = "unknown",
      method = "limma",
      transform_for_limma = "log2(x + 1)",
      collapse_method = "mean",
      reason = "numeric values too few"
    ))
  }
  vals_sample <- if (length(vals) > 200000) sample(vals, 200000) else vals
  nonneg <- mean(vals_sample >= 0) > 0.995
  int_prop <- mean(abs(vals_sample - round(vals_sample)) < 1e-6)
  zero_prop <- mean(vals_sample == 0)
  q <- as.numeric(quantile(vals_sample, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  names(q) <- c("min", "q25", "median", "q75", "q99", "max")

  hint <- paste(c(source_hint, attr(mat_or_df, "source_hint")), collapse = " ")
  hint <- tolower(hint)
  normalized_hint <- grepl(
    "tpm|fpkm|rpkm|rpm|cpm|normalized|normalised|norm[._ -]?count|abundance|expression[._ -]?value",
    hint, perl = TRUE
  )
  raw_hint <- grepl(
    "raw[._ -]?count|read[._ -]?count|gene[._ -]?count|count[._ -]?matrix|counts[._ -]?matrix|featurecounts|htseq|\\bcounts?\\b",
    hint, perl = TRUE
  ) && !normalized_hint
  log_hint <- grepl("log2|\\brlog\\b|\\bvst\\b|\\brma\\b|\\bmas5\\b|microarray|quantile", hint, perl = TRUE)

  col_sums <- suppressWarnings(colSums(x, na.rm = TRUE))
  col_sums <- col_sums[is.finite(col_sums) & col_sums > 0]
  lib_cv <- if (length(col_sums) >= 3 && mean(col_sums) > 0) stats::sd(col_sums) / mean(col_sums) else NA_real_

  # Integer-like values are not sufficient evidence for raw counts: rounded TPM,
  # CPM and normalized counts can look integer-like too.  Auto-DESeq2 therefore
  # requires either an explicit count source name or a very strong count signature.
  strong_numeric_count <- nonneg && int_prop >= 0.995 && zero_prop >= 0.20 &&
    q["q99"] >= 100 && q["max"] >= 1000 &&
    !is.na(lib_cv) && lib_cv >= 0.10

  if (!normalized_hint && !log_hint &&
      ((raw_hint && nonneg && int_prop >= 0.95 && q["max"] >= 20) || strong_numeric_count)) {
    return(list(
      type = "raw_count",
      method = "DESeq2",
      transform_for_limma = "none",
      collapse_method = "sum",
      reason = paste0(
        "conservative raw-count evidence: source_hint=", ifelse(raw_hint, "count-like", "none"),
        "; integers=", round(int_prop * 100, 1), "%; zeros=", round(zero_prop * 100, 1),
        "%; q99=", round(q["q99"], 2), "; max=", round(q["max"], 2),
        ifelse(is.na(lib_cv), "", paste0("; library-size CV=", round(lib_cv, 3)))
      )
    ))
  }

  # Positive RNA abundance is often linear even when its range is low.  A sizeable
  # zero fraction or a long right tail separates low-range TPM/FPKM/CPM from most
  # already-log microarray matrices.  Explicit normalized names take precedence.
  upper_tail_ratio <- q["q99"] / max(q["q75"], 0.25)
  middle_tail_ratio <- q["q75"] / max(q["median"], 0.25)
  looks_linear_positive <- nonneg && (
    normalized_hint ||
      q["q99"] > 25 || q["max"] > 150 ||
      zero_prop >= 0.05 || upper_tail_ratio >= 3.5 || middle_tail_ratio >= 2.5
  )
  if (looks_linear_positive && !log_hint) {
    return(list(
      type = "normalized_expression_FPKM_TPM_or_similar",
      method = "limma",
      transform_for_limma = "log2(x + 1)",
      collapse_method = "mean",
      reason = paste0(
        "positive linear normalized/submitted expression; source_hint=",
        ifelse(normalized_hint, "normalized-like", "none"),
        "; zeros=", round(zero_prop * 100, 1), "%; q99=", round(q["q99"], 2),
        "; max=", round(q["max"], 2), "; upper-tail ratio=", round(upper_tail_ratio, 2),
        "; limma will use log2(x+1)"
      )
    ))
  }

  # Negative values, an explicit log/microarray hint, or a compressed non-zero
  # distribution are the only automatic already-log cases.
  if (q["min"] < 0 || log_hint ||
      (nonneg && zero_prop < 0.05 && q["q99"] <= 25 && q["max"] <= 150 &&
       upper_tail_ratio < 3.5 && middle_tail_ratio < 2.5)) {
    return(list(
      type = "log_or_microarray_expression",
      method = "limma",
      transform_for_limma = "none",
      collapse_method = "mean",
      reason = paste0(
        "values look already log-scale/microarray-like; source_hint=",
        ifelse(log_hint, "log/microarray-like", "none"),
        "; zeros=", round(zero_prop * 100, 1), "; q99=", round(q["q99"], 2),
        "; max=", round(q["max"], 2)
      )
    ))
  }

  # FPKM/TPM/RPM/normalized positive expression.
  return(list(
    type = "normalized_expression_FPKM_TPM_or_similar",
    method = "limma",
    transform_for_limma = "log2(x + 1)",
    collapse_method = "mean",
    reason = paste0(
      "ambiguous positive matrix is handled conservatively as linear normalized expression; ",
      "integers=", round(int_prop * 100, 1), "%; zeros=", round(zero_prop * 100, 1),
      "%; q99=", round(q["q99"], 2), "; source hint did not prove raw counts"
    )
  ))
}




# v11.4 helper: whenever an expression matrix is loaded from old parsed cache,
# re-run expression-type detection using the current conservative thresholds.
# This prevents stale cached objects from keeping an old incorrect label such as
# log_or_microarray_expression for large positive submitted matrices.
refresh_expression_detect_current <- function(obj) {
  if (is.null(obj) || is.null(obj$raw)) return(obj)
  raw <- as.data.frame(obj$raw)
  if (ncol(raw) < 3) return(obj)
  det_new <- tryCatch({
    detect_expression_type(
      raw[, -1, drop = FALSE],
      source_hint = paste(obj$expr_file, attr(obj$raw, "source_hint"), collapse = " ")
    )
  }, error = function(e) NULL)
  if (is.null(det_new)) return(obj)

  old_type <- tryCatch({
    if (!is.null(obj$expr_detect$type)) obj$expr_detect$type else attr(obj$raw, "expr_detect")$type
  }, error = function(e) NA_character_)

  obj$expr_detect <- det_new
  attr(obj$raw, "expr_detect") <- det_new
  attr(obj$raw, "collapse_method") <- det_new$collapse_method

  if (is.null(obj$msg)) obj$msg <- ""
  if (!identical(as.character(old_type), as.character(det_new$type))) {
    obj$msg <- paste0(
      obj$msg,
      "\nExpression type was re-detected by current v11.4 rules: ",
      old_type, " -> ", det_new$type,
      "\nCurrent detection reason: ", det_new$reason,
      "\nCurrent recommended method: ", det_new$method,
      "\nCurrent limma transform: ", det_new$transform_for_limma
    )
  } else {
    obj$msg <- paste0(
      obj$msg,
      "\nExpression type rechecked by current v11.7 rules: ", det_new$type,
      " | ", det_new$reason,
      " | limma transform: ", det_new$transform_for_limma
    )
  }
  obj
}

# =========================
# Entrez ID -> Gene Symbol helper
# =========================
looks_like_entrez_id_vector <- function(x) {
  x <- as.character(x)
  x <- sub("\\.0$", "", x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^[0-9]+$", x)) > 0.70
}

entrez_to_symbol <- function(ids) {
  ids0 <- as.character(ids)
  ids_clean <- sub("\\.0$", "", ids0)

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
  }

  mapped <- tryCatch({
    suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(ids_clean[!is.na(ids_clean) & ids_clean != ""]),
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    ))
  }, error = function(e) NULL)

  if (is.null(mapped)) return(ids0)

  out <- ids_clean
  hit <- mapped[ids_clean]
  ok <- !is.na(hit) & hit != ""
  out[ok] <- as.character(hit[ok])
  out
}

convert_gene_ids_to_symbols_auto <- function(symbols) {
  symbols0 <- as.character(symbols)
  symbols0 <- clean_symbol_value(symbols0)

  if (looks_like_entrez_id_vector(symbols0)) {
    converted <- entrez_to_symbol(symbols0)
    return(list(
      symbols = converted,
      id_type = "ENTREZID",
      reason = "numeric gene_id detected; converted ENTREZID -> SYMBOL",
      converted_n = sum(converted != sub("\\.0$", "", symbols0), na.rm = TRUE)
    ))
  }

  if (mean(grepl("^ENSG", symbols0, ignore.case = TRUE), na.rm = TRUE) > 0.50) {
    converted <- ensembl_to_symbol(symbols0)
    return(list(
      symbols = converted,
      id_type = "ENSEMBL",
      reason = "ENSEMBL detected; converted ENSEMBL -> SYMBOL",
      converted_n = sum(converted != symbols0, na.rm = TRUE)
    ))
  }

  list(
    symbols = symbols0,
    id_type = "SYMBOL_or_other",
    reason = "no ID conversion needed",
    converted_n = 0
  )
}


collapse_by_symbol_auto <- function(raw, forced_method = NULL, source_hint = NULL) {
  raw <- as.data.frame(raw)
  colnames(raw)[1] <- "Symbol"

  raw$Symbol <- clean_symbol_value(raw$Symbol)

  id_conv <- convert_gene_ids_to_symbols_auto(raw$Symbol)
  raw$Symbol <- clean_symbol_value(id_conv$symbols)

  raw <- raw[!is.na(raw$Symbol) & raw$Symbol != "" & raw$Symbol != "NA", , drop = FALSE]

  for (cc in colnames(raw)[-1]) {
    raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
  }

  det <- detect_expression_type(raw[-1], source_hint = source_hint)
  collapse_method <- if (!is.null(forced_method)) forced_method else det$collapse_method

  if (collapse_method == "sum") {
    raw2 <- raw %>%
      dplyr::group_by(Symbol) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")
  } else {
    raw2 <- raw %>%
      dplyr::group_by(Symbol) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  }

  raw2 <- as.data.frame(raw2)

  attr(raw2, "expr_detect") <- det
  attr(raw2, "collapse_method") <- collapse_method
  attr(raw2, "gene_id_type") <- id_conv$id_type
  attr(raw2, "gene_id_conversion_reason") <- id_conv$reason
  attr(raw2, "gene_id_converted_n") <- id_conv$converted_n
  attr(raw2, "source_hint") <- source_hint

  raw2
}

# =========================
# 3. Series matrix parser
# =========================
get_symbol_from_gpl <- function(gpl_id, expr_ids) {
  gpl <- getGEO(gpl_id)
  gpl_tab <- Table(gpl)
  if (is.null(gpl_tab) || nrow(gpl_tab) == 0) return(NULL)
  id_candidates <- colnames(gpl_tab)[
    grepl("^ID$|ID_REF|probe|Probe|transcript|Transcript", colnames(gpl_tab), ignore.case = TRUE)
  ]
  symbol_candidates <- colnames(gpl_tab)[
    grepl("symbol|gene.symbol|gene symbol|GENE_SYMBOL|Symbol|Gene.Symbol|gene_name|Gene.Name",
          colnames(gpl_tab), ignore.case = TRUE)
  ]
  if (length(id_candidates) == 0 || length(symbol_candidates) == 0) return(NULL)
  ids <- as.character(gpl_tab[[id_candidates[1]]])
  symbols <- clean_symbol_value(gpl_tab[[symbol_candidates[1]]])
  symbols[match(expr_ids, ids)]
}

get_series_matrix_expr <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
  eset <- gse[[1]]
  expr <- exprs(eset)
  meta <- pData(eset)
  feature <- fData(eset)
  expr_ids <- rownames(expr)

  gene_symbol <- NULL
  symbol_col <- NA_character_

  # ---------- helper: is parsed symbol acceptable? ----------
  symbol_quality <- function(sym, expr_n) {
    sym <- clean_symbol_value(sym)
    valid <- !is.na(sym) & sym != "" & sym != "NA"
    n_valid <- sum(valid)
    n_unique <- length(unique(sym[valid]))
    numeric_prop <- if (n_valid > 0) mean(grepl("^[0-9]+$", sym[valid])) else 1
    accession_prop <- if (n_valid > 0) {
      mean(grepl("^NM_|^NR_|^XM_|^XR_|^ENSG|^ENST|^AFFX|^ILMN|_at$|_s_at$|_x_at$",
                 sym[valid], ignore.case = TRUE))
    } else 1
    symbol_like_prop <- if (n_valid > 0) {
      mean(grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", sym[valid]) &
             !grepl("^NM_|^NR_|^XM_|^XR_|^ENSG|^ENST|^AFFX|^ILMN|^[0-9]+$",
                    sym[valid], ignore.case = TRUE))
    } else 0
    duplicate_n <- n_valid - n_unique
    list(
      n_valid = n_valid,
      n_unique = n_unique,
      numeric_prop = numeric_prop,
      accession_prop = accession_prop,
      symbol_like_prop = symbol_like_prop,
      duplicate_n = duplicate_n,
      ok = n_valid >= max(100, expr_n * 0.30) &&
        symbol_like_prop >= 0.30 &&
        numeric_prop < 0.50 &&
        accession_prop < 0.70
    )
  }

  # ---------- 1) Try feature gene_assignment ----------
  if (ncol(feature) > 0 && "gene_assignment" %in% colnames(feature)) {
    gene_symbol_try <- parse_gene_assignment_symbol(feature$gene_assignment)
    q <- symbol_quality(gene_symbol_try, nrow(expr))
    message("gene_assignment解析质量: valid=", q$n_valid,
            ", unique=", q$n_unique,
            ", duplicate=", q$duplicate_n,
            ", symbol_like=", round(q$symbol_like_prop, 3),
            ", numeric=", round(q$numeric_prop, 3),
            ", accession=", round(q$accession_prop, 3))
    if (isTRUE(q$ok)) {
      gene_symbol <- gene_symbol_try
      symbol_col <- "gene_assignment"
    }
  }

  # ---------- 2) Try explicit symbol-like feature columns ----------
  if ((is.null(gene_symbol) || length(gene_symbol) != nrow(expr) || all(is.na(gene_symbol))) && ncol(feature) > 0) {
    symbol_candidates <- colnames(feature)[
      grepl("gene.symbol|gene symbol|GENE_SYMBOL|^Symbol$|symbol|gene_name|Gene.Name|GENE_NAME|external_gene_name",
            colnames(feature), ignore.case = TRUE)
    ]

    if (length(symbol_candidates) > 0) {
      for (cc in symbol_candidates) {
        gene_symbol_try <- clean_symbol_value(feature[[cc]])
        q <- symbol_quality(gene_symbol_try, nrow(expr))
        message("Feature symbol候选列检查: ", cc,
                " | valid=", q$n_valid,
                ", unique=", q$n_unique,
                ", duplicate=", q$duplicate_n,
                ", symbol_like=", round(q$symbol_like_prop, 3),
                ", numeric=", round(q$numeric_prop, 3),
                ", accession=", round(q$accession_prop, 3))
        if (isTRUE(q$ok)) {
          gene_symbol <- gene_symbol_try
          symbol_col <- cc
          break
        }
      }
    }
  }

  # ---------- 3) Try GPL annotation ----------
  if (is.null(gene_symbol) || length(gene_symbol) != nrow(expr) || all(is.na(gene_symbol))) {
    gpl_id <- annotation(eset)
    if (!is.null(gpl_id) && gpl_id != "") {
      message("Series Matrix没有可靠symbol，尝试从GPL平台注释解析：", gpl_id)
      gene_symbol_try <- get_symbol_from_gpl(gpl_id, expr_ids)
      q <- symbol_quality(gene_symbol_try, nrow(expr))
      message("GPL symbol解析质量: valid=", q$n_valid,
              ", unique=", q$n_unique,
              ", duplicate=", q$duplicate_n,
              ", symbol_like=", round(q$symbol_like_prop, 3),
              ", numeric=", round(q$numeric_prop, 3),
              ", accession=", round(q$accession_prop, 3))
      if (isTRUE(q$ok)) {
        gene_symbol <- gene_symbol_try
        symbol_col <- paste0("GPL annotation: ", gpl_id)
      }
    }
  }

  # ---------- 4) Ensembl rownames fallback ----------
  if (is.null(gene_symbol) || length(gene_symbol) != nrow(expr) || all(is.na(gene_symbol))) {
    if (any(grepl("^ENSG", expr_ids))) {
      gene_symbol <- ensembl_to_symbol(expr_ids)
      symbol_col <- "rownames Ensembl ID -> SYMBOL"
    }
  }

  if (is.null(gene_symbol) || length(gene_symbol) != nrow(expr) || all(is.na(gene_symbol))) {
    stop(paste0(
      "Series Matrix/GPL无法解析可靠gene symbol。\n",
      "当前fData列名：\n", paste(colnames(feature), collapse = ", "),
      "\n表达矩阵行名示例：\n", paste(head(expr_ids, 20), collapse = ", "),
      "\n建议：检查GPL注释列，或使用probe-level/no-merge模式。"
    ))
  }

  # ---------- Build and collapse expression ----------
  expr_df <- as.data.frame(expr)
  expr_df$Symbol <- gene_symbol
  expr_df <- expr_df[, c("Symbol", setdiff(colnames(expr_df), "Symbol"))]

  before_n <- nrow(expr_df)
  before_unique <- length(unique(clean_symbol_value(expr_df$Symbol)))
  before_dup <- sum(duplicated(clean_symbol_value(expr_df$Symbol)))

  expr_df <- collapse_by_symbol_auto(expr_df, source_hint = paste("Series Matrix", gse_id, symbol_col))
  det <- attr(expr_df, "expr_detect")

  # ---------- Console diagnostics ----------
  cat("\n========== SYMBOL CHECK ==========\n")
  cat("GSE:", gse_id, "\n")
  cat("Symbol source:", symbol_col, "\n")
  cat("Rows before collapse:", before_n, "\n")
  cat("Unique symbols before collapse:", before_unique, "\n")
  cat("Duplicated symbols before collapse:", before_dup, "\n")
  cat("Rows after collapse:", nrow(expr_df), "\n")
  cat("Columns:", ncol(expr_df), "\n")
  cat("First 20 Symbols after collapse:\n")
  print(head(expr_df$Symbol, 20))
  cat("Duplicated Symbols after collapse:", sum(duplicated(expr_df$Symbol)), "\n")
  cat("Unique Symbols after collapse:", length(unique(expr_df$Symbol)), "\n")
  cat("Expression detect:", det$type, "\n")
  cat("Collapse method:", attr(expr_df, "collapse_method"), "\n")
  cat("==================================\n")

  # Warn if still probe-like / almost no collapse occurred.
  warning_msg <- ""
  if (nrow(expr_df) > 45000 || nrow(expr_df) > before_n * 0.90) {
    warning_msg <- paste0(
      "\n⚠️ SYMBOL WARNING:\n",
      "Rows after collapse remains very high (", nrow(expr_df), ").\n",
      "This may indicate that the platform annotation is probe/transcript-cluster level,\n",
      "or Gene Symbol parsing did not reduce probes into gene-level symbols.\n",
      "Please inspect the SYMBOL CHECK output above.\n"
    )
  }

  list(
    raw = expr_df,
    meta = as.data.frame(meta),
    expr_file = "GEO Series Matrix exprs()",
    meta_file = "GEO Series Matrix pData()",
    expr_detect = det,
    msg = paste0(
      "使用 GEO Series Matrix 表达矩阵。\n",
      "Symbol解析方式：", symbol_col, "\n",
      "Rows before collapse: ", before_n, "\n",
      "Rows after collapse: ", nrow(expr_df), "\n",
      "Unique symbols before collapse: ", before_unique, "\n",
      "Duplicated symbols before collapse: ", before_dup, "\n",
      "数据类型自动判断：", det$type, "\n",
      "推荐差异分析方法：", det$method, "\n",
      "重复Symbol合并方法：", attr(expr_df, "collapse_method"), "\n",
      "判断原因：", det$reason,
      warning_msg
    )
  )
}

# =========================
# 4. Supplementary parser
# =========================
read_table_any <- function(f) {
  tryCatch({
    if (grepl("\\.xlsx$|\\.xls$", f, ignore.case = TRUE)) {
      # v7.5: read Excel without assuming the first physical row is the true header.
      # Many GEO supplementary Excel files have section headers (Raw data / Annotation)
      # on row 1 and the real column names on row 2. We repair this later generically.
      as.data.frame(readxl::read_excel(f, col_names = FALSE))
    } else {
      data.table::fread(f, data.table = FALSE, fill = TRUE, check.names = FALSE)
    }
  }, error = function(e) {
    message("读取失败: ", basename(f), " | ", e$message)
    NULL
  })
}

is_generic_excel_section_header_row <- function(x) {
  x <- tolower(paste(as.character(x), collapse = " "))
  grepl("raw data|annotation|gene information|processed data|normalized data|expression data", x)
}

looks_like_expression_header_row <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  txt <- tolower(paste(x, collapse = " "))
  has_gene <- grepl("gene|symbol|ensembl|entrez|id", txt)
  sample_hits <- sum(grepl("^gsm[0-9]+$|control|healthy|normal|case|sepsis|patient|sample|^s[0-9]+$|^c[0-9]+$|^h[0-9]+$", x, ignore.case = TRUE), na.rm = TRUE)
  has_gene && sample_hits >= 2
}

repair_expression_table_header <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  if (is.null(dat) || nrow(dat) < 2 || ncol(dat) < 2) return(dat)

  # Drop fully empty rows/columns first.
  dat <- dat[rowSums(!is.na(dat) & dat != "") > 0, , drop = FALSE]
  dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
  if (nrow(dat) < 2 || ncol(dat) < 2) return(dat)

  # If read_excel(col_names=FALSE) generated .../X names, find the real header row.
  candidate_rows <- seq_len(min(8, nrow(dat)))
  header_idx <- NA_integer_
  for (i in candidate_rows) {
    if (looks_like_expression_header_row(dat[i, ])) {
      header_idx <- i
      break
    }
  }

  if (!is.na(header_idx)) {
    new_names <- as.character(unlist(dat[header_idx, ], use.names = FALSE))
    new_names[is.na(new_names) | trimws(new_names) == ""] <- paste0("X", which(is.na(new_names) | trimws(new_names) == ""))
    colnames(dat) <- make.names(new_names, unique = TRUE)
    dat <- dat[-seq_len(header_idx), , drop = FALSE]
  } else {
    # Older behavior for txt/csv files where the first data row contains sample names.
    first_row <- as.character(dat[1, ])
    if (looks_like_expression_header_row(first_row) ||
        sum(grepl("^GSM|^HC[0-9]+$|^S[0-9]+$|sepsis|control|healthy|sample|patient", first_row, ignore.case = TRUE), na.rm = TRUE) >= 2) {
      new_names <- make.names(ifelse(is.na(first_row) | first_row == "", colnames(dat), first_row), unique = TRUE)
      colnames(dat) <- new_names
      dat <- dat[-1, , drop = FALSE]
    } else {
      colnames(dat) <- make.names(colnames(dat), unique = TRUE)
    }
  }

  dat <- dat[rowSums(!is.na(dat) & dat != "") > 0, , drop = FALSE]
  dat
}

sample_name_from_file <- function(f) {
  bn <- basename(f)
  bn <- sub("\\.gz$", "", bn, ignore.case = TRUE)
  bn <- sub("\\.htseq\\.results$", "", bn, ignore.case = TRUE)
  bn <- sub("\\.txt$|\\.tsv$|\\.csv$|\\.count$|\\.counts$|\\.sf$|\\.out$", "", bn, ignore.case = TRUE)
  bn <- sub("_COUNT$|_counts$|_count$|\\.COUNT$|\\.counts$|\\.count$", "", bn, ignore.case = TRUE)

  # Keep informative full sample name, e.g. GSM7123641_P2 not just GSM7123641.
  # Matching functions know how to match by full/GSM/short aliases.
  bn
}

read_one_single_sample_expr_file <- function(f) {
  dat <- tryCatch({
    fread(f, header = FALSE, data.table = FALSE, fill = TRUE)
  }, error = function(e) {
    message("读取单样本文件失败: ", basename(f), " | ", e$message)
    NULL
  })
  if (is.null(dat) || nrow(dat) < 100 || ncol(dat) < 2) return(NULL)
  dat <- dat[, 1:2, drop = FALSE]
  colnames(dat) <- c("Symbol", "count")
  dat$Symbol <- as.character(dat$Symbol)
  num <- suppressWarnings(as.numeric(dat$count))
  if (is.na(num[1]) && sum(!is.na(num[-1])) >= 100) {
    dat <- dat[-1, , drop = FALSE]
    num <- suppressWarnings(as.numeric(dat$count))
  }
  dat$count <- num
  dat <- dat[
    !is.na(dat$Symbol) & dat$Symbol != "" &
      !grepl("^__", dat$Symbol) &
      !is.na(dat$count),
    , drop = FALSE
  ]
  if (nrow(dat) < 100 || sum(!is.na(dat$count)) < 100) return(NULL)
  dat
}

merge_single_sample_expr_files <- function(files) {
  files <- files[!duplicated(basename(files))]
  count_list <- list()
  failed <- 0L
  for (f in files) {
    dat <- read_one_single_sample_expr_file(f)
    if (is.null(dat)) next
    sample_name <- sample_name_from_file(f)
    dat <- collapse_single_sample_duplicates(dat)
    if (is.null(dat) || nrow(dat) < 100 || anyDuplicated(dat$Symbol) > 0) {
      failed <- failed + 1L
      next
    }

    colnames(dat)[2] <- sample_name
    count_list[[sample_name]] <- dat
  }
  if (length(count_list) < 2) return(NULL)
  raw <- Reduce(function(x, y) full_join(x, y, by = "Symbol"), count_list)
  raw[is.na(raw)] <- 0
  collapse_by_symbol_auto(raw, forced_method = "sum", source_hint = paste(basename(files), collapse = " "))
}


# Fast generic reader for submitted gene x sample expression matrices.
# Handles GEO supplementary files such as:
# - *_AllSampleExpressionSubmitted.tsv.gz
# - *_Processed_data_FPKM.csv
# - *_expression_matrix.tsv.gz
# It is not GSE-specific. It recognizes a common structure:
# first column = gene symbol / gene id, remaining columns = samples.
# It also removes non-expression annotation rows such as SOFAVALUE if present.

# Ultra-fast generic reader for large submitted expression matrices.
# Main purpose:
# - do NOT scan every column repeatedly
# - do NOT call slow symbol conversion unless rownames are Ensembl IDs
# - do NOT fall back to slow parser for files that clearly match submitted matrix format
try_ultra_fast_expression_matrix <- function(f) {
  bn <- basename(f)

  if (!grepl("AllSampleExpressionSubmitted|Processed|processed|mRNA|RNA[_-]?seq|RNAseq|Seq|Genes|GeneList|Excel|FPKM|TPM|RPM|RPKM|CPM|expression.*matrix|Expression.*Matrix|normalized|Normalized|GeneExpression",
             bn, ignore.case = TRUE)) {
    return(NULL)
  }

  message("超快速矩阵读取模式: ", bn)

  dt <- tryCatch({
    data.table::fread(
      f,
      data.table = TRUE,
      fill = TRUE,
      check.names = FALSE,
      showProgress = FALSE
    )
  }, error = function(e) {
    message("超快速矩阵读取失败: ", bn, " | ", e$message)
    NULL
  })

  if (is.null(dt) || nrow(dt) < 100 || ncol(dt) < 3) return(NULL)

  dat_df <- as.data.frame(dt)
  colnames(dat_df) <- make.names(colnames(dat_df), unique = TRUE)

  gene_choice <- choose_gene_symbol_column(dat_df)
  symbol_col_used <- gene_choice$symbol_col
  gene_id_col_used <- gene_choice$id_col

  # Numeric expression/sample columns only, excluding annotation columns.
  all_cols <- colnames(dat_df)
  numeric_cols <- all_cols[sapply(dat_df, is_numeric_like, min_prop = 0.50)]
  numeric_cols <- safe_char_cols(numeric_cols, colnames(dat_df))
  sample_cols <- numeric_cols[!is_annotation_col(numeric_cols)]
  sample_cols <- setdiff(sample_cols, c(symbol_col_used, gene_id_col_used))

  # If sample names are not yet numeric due fread type detection, try all non-annotation columns and convert.
  if (length(sample_cols) < 2) {
    sample_cols <- setdiff(all_cols, c(symbol_col_used, gene_id_col_used))
    sample_cols <- sample_cols[!is_annotation_col(sample_cols)]
  }

  if (length(sample_cols) < 2) return(NULL)

  raw <- dat_df[, c(symbol_col_used, sample_cols), drop = FALSE]
  colnames(raw)[1] <- "Symbol"

  raw$Symbol <- clean_symbol_value(raw$Symbol)

  # Remove metadata/trait rows accidentally included in expression matrix.
  bad_symbol_rows <- is.na(raw$Symbol) |
    raw$Symbol == "" |
    grepl("SOFA|SOFAVALUE|GROUP|CLASS|PHENOTYPE|DISEASE|STATUS|OUTCOME|MORT|AGE|SEX|GENDER|RACE|BATCH|SAMPLE",
          raw$Symbol,
          ignore.case = TRUE)

  raw <- raw[!bad_symbol_rows, , drop = FALSE]

  if (nrow(raw) < 100) return(NULL)

  # Convert gene IDs immediately, before collapse.
  id_conv <- convert_gene_ids_to_symbols_auto(raw$Symbol)
  raw$Symbol <- clean_symbol_value(id_conv$symbols)

  for (cc in colnames(raw)[-1]) {
    raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
  }

  # Remove columns with too few numeric values.
  good_cols <- colnames(raw)[-1][sapply(raw[-1], function(x) mean(!is.na(x)) >= 0.50)]

  if (length(good_cols) < 2) return(NULL)

  raw <- raw[, c("Symbol", good_cols), drop = FALSE]

  keep <- rowSums(!is.na(raw[, good_cols, drop = FALSE])) >= max(2, floor(length(good_cols) * 0.5))
  raw <- raw[keep, , drop = FALSE]

  if (nrow(raw) < 100 || ncol(raw) < 3) return(NULL)

  # Expression type detection on subset before collapsing.
  sub_rows <- if (nrow(raw) > 2000) sample(seq_len(nrow(raw)), 2000) else seq_len(nrow(raw))
  sub_cols <- if (length(good_cols) > 60) sample(good_cols, 60) else good_cols
  det <- detect_expression_type(as.data.frame(raw[sub_rows, sub_cols, drop = FALSE]), source_hint = bn)

  collapse_method <- det$collapse_method

  if (collapse_method == "sum") {
    raw_dt <- data.table::as.data.table(raw)
    raw_dt <- raw_dt[, lapply(.SD, sum, na.rm = TRUE), by = Symbol, .SDcols = good_cols]
  } else {
    raw_dt <- data.table::as.data.table(raw)
    raw_dt <- raw_dt[, lapply(.SD, mean, na.rm = TRUE), by = Symbol, .SDcols = good_cols]
  }

  raw2 <- as.data.frame(raw_dt)

  attr(raw2, "expr_detect") <- det
  attr(raw2, "collapse_method") <- collapse_method
  attr(raw2, "symbol_col_used") <- symbol_col_used
  attr(raw2, "symbol_col_reason") <- gene_choice$reason
  attr(raw2, "gene_id_col_used") <- gene_id_col_used
  attr(raw2, "gene_id_type") <- id_conv$id_type
  attr(raw2, "gene_id_conversion_reason") <- id_conv$reason
  attr(raw2, "gene_id_converted_n") <- id_conv$converted_n

  if (nrow(raw2) > 100 && ncol(raw2) > 3) {
    return(list(
      raw = raw2,
      symbol_col = attr(raw2, "symbol_col_used"),
      symbol_col_reason = attr(raw2, "symbol_col_reason"),
      gene_id_col = attr(raw2, "gene_id_col_used"),
      gene_id_type = attr(raw2, "gene_id_type"),
      gene_id_conversion_reason = attr(raw2, "gene_id_conversion_reason"),
      gene_id_converted_n = attr(raw2, "gene_id_converted_n"),
      expr_detect = det,
      fast_reader = TRUE,
      ultra_fast_reader = TRUE
    ))
  }

  NULL
}


try_fast_submitted_expression_matrix <- function(f) {
  bn <- basename(f)

  # Only apply to likely processed expression matrices.
  likely_matrix <- grepl(
    "AllSampleExpressionSubmitted|Processed|processed|mRNA|RNA[_-]?seq|RNAseq|Seq|Genes|GeneList|Excel|expression|Expression|matrix|Matrix|FPKM|TPM|normalized|Normalized",
    bn,
    ignore.case = TRUE
  )

  if (!likely_matrix) return(NULL)

  message("快速矩阵读取模式: ", bn)

  dat <- tryCatch({
    data.table::fread(f, data.table = FALSE, fill = TRUE, check.names = FALSE)
  }, error = function(e) {
    message("快速矩阵读取失败: ", bn, " | ", e$message)
    NULL
  })

  if (is.null(dat) || nrow(dat) < 100 || ncol(dat) < 3) return(NULL)

  # Remove empty columns.
  dat <- dat[, colSums(!is.na(dat)) > 0, drop = FALSE]
  if (ncol(dat) < 3) return(NULL)

  # Identify symbol column.
  cn <- colnames(dat)
  symbol_candidates <- cn[
    grepl("symbol|gene.symbol|gene_name|gene.name|gene$|Gene$|SYMBOL|external_gene_name|ensembl|id$|ID$|V1$",
          cn, ignore.case = TRUE)
  ]

  symbol_col <- if (length(symbol_candidates) > 0) symbol_candidates[1] else cn[1]

  # First column may have no clean name but values are gene symbols.
  colnames(dat)[colnames(dat) == symbol_col] <- "Symbol"
  symbol_col <- "Symbol"

  # Remove annotation / phenotype rows accidentally included in expression matrix.
  # Example GSE310929: second row SOFAVALUE with values like "sofa: 3".
  dat$Symbol <- as.character(dat$Symbol)
  bad_symbol_rows <- is.na(dat$Symbol) |
    dat$Symbol == "" |
    grepl("SOFA|SOFAVALUE|GROUP|CLASS|PHENOTYPE|DISEASE|STATUS|OUTCOME|MORT|AGE|SEX|GENDER|RACE|BATCH",
          dat$Symbol,
          ignore.case = TRUE)

  # Keep rows that look like real gene names or Ensembl IDs.
  # This prevents metadata rows from slowing numeric conversion.
  dat <- dat[!bad_symbol_rows, , drop = FALSE]

  if (nrow(dat) < 100) return(NULL)

  # Candidate sample columns: columns whose values are mostly numeric after removing bad rows.
  sample_cols <- setdiff(colnames(dat), "Symbol")
  numeric_prop <- sapply(sample_cols, function(cc) {
    x <- suppressWarnings(as.numeric(dat[[as.character(cc)[1]]]))
    mean(!is.na(x))
  })

  sample_cols <- sample_cols[numeric_prop >= 0.50]
  if (length(sample_cols) < 2) return(NULL)

  raw <- dat[, c("Symbol", sample_cols), drop = FALSE]

  # Convert only selected sample columns.
  for (cc in sample_cols) {
    raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
  }

  # Remove genes with too many missing values.
  keep_row <- rowSums(!is.na(raw[, sample_cols, drop = FALSE])) >= max(2, floor(length(sample_cols) * 0.5))
  raw <- raw[keep_row, , drop = FALSE]

  if (nrow(raw) < 100 || ncol(raw) < 3) return(NULL)

  raw <- collapse_by_symbol_auto(raw, source_hint = bn)

  if (nrow(raw) > 100 && ncol(raw) > 3) {
    return(list(
      raw = raw,
      symbol_col = "Symbol",
      expr_detect = attr(raw, "expr_detect"),
      fast_reader = TRUE
    ))
  }

  NULL
}



# =========================
# Robust gene symbol column chooser
# =========================
looks_like_xloc_or_transcript_id <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^XLOC_|^MSTRG\\.|^TCONS_|^CUFF\\.", x, ignore.case = TRUE)) > 0.30
}

looks_like_gene_symbol_vector <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "---"]
  if (length(x) == 0) return(FALSE)
  prop <- mean(grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", x))
  prop > 0.50
}

choose_gene_symbol_column <- function(dat) {
  cn <- colnames(dat)

  strong_symbol_cols <- cn[
    grepl("^gene_name$|^gene.symbol$|^gene_symbol$|^symbol$|^hugo$|^hgnc_symbol$|external_gene_name",
          cn,
          ignore.case = TRUE)
  ]

  if (length(strong_symbol_cols) > 0) {
    for (cc in strong_symbol_cols) {
      if (looks_like_gene_symbol_vector(dat[[as.character(cc)[1]]])) {
        return(list(symbol_col = cc, id_col = cn[1], reason = paste0("preferred symbol column: ", cc)))
      }
    }
  }

  first_col <- cn[1]
  if (looks_like_xloc_or_transcript_id(dat[[first_col]])) {
    candidate_cols <- cn[
      grepl("gene_name|gene.symbol|gene_symbol|symbol|hugo|hgnc|nearest_refseq|refseq",
            cn,
            ignore.case = TRUE)
    ]
    candidate_cols <- setdiff(candidate_cols, first_col)

    for (cc in candidate_cols) {
      if (looks_like_gene_symbol_vector(dat[[as.character(cc)[1]]])) {
        return(list(symbol_col = cc, id_col = first_col, reason = paste0("first column is XLOC-like; using annotation column: ", cc)))
      }
    }
  }

  generic_candidates <- cn[
    grepl("symbol|gene_name|gene.name|gene$|Gene$|SYMBOL|external_gene_name|hugo|hgnc|ensembl|id$|ID$",
          cn,
          ignore.case = TRUE)
  ]

  if (length(generic_candidates) > 0) {
    return(list(symbol_col = generic_candidates[1], id_col = cn[1], reason = paste0("generic candidate: ", generic_candidates[1])))
  }

  non_num <- cn[!sapply(dat, is_numeric_like, min_prop = 0.5)]
  if (length(non_num) > 0) {
    return(list(symbol_col = non_num[1], id_col = cn[1], reason = paste0("first non-numeric column: ", non_num[1])))
  }

  list(symbol_col = cn[1], id_col = cn[1], reason = "fallback first column")
}

prepare_raw_from_table_with_symbol <- function(dat, sample_cols) {
  gene_choice <- choose_gene_symbol_column(dat)

  symbol_col <- safe_one_col(gene_choice$symbol_col, colnames(dat), fallback = colnames(dat)[1])
  id_col <- safe_one_col(gene_choice$id_col, colnames(dat), fallback = colnames(dat)[1])

  sample_cols <- safe_char_cols(sample_cols, colnames(dat))
  sample_cols <- setdiff(sample_cols, c(symbol_col, id_col))
  sample_cols <- sample_cols[!grepl("^gene$|gene_name|orig_id|nearest_refseq|length|gc.content|gc content|symbol|refseq|description",
                                    sample_cols,
                                    ignore.case = TRUE)]

  sample_cols <- safe_char_cols(sample_cols, colnames(dat))
  raw_cols <- safe_char_cols(c(symbol_col, sample_cols), colnames(dat))
  if (length(raw_cols) < 2) return(NULL)
  raw <- dat[, raw_cols, drop = FALSE]
  colnames(raw)[1] <- "Symbol"

  raw$Symbol <- clean_symbol_value(raw$Symbol)

  id_conv <- convert_gene_ids_to_symbols_auto(raw$Symbol)
  raw$Symbol <- clean_symbol_value(id_conv$symbols)
  attr(raw, "gene_id_type") <- id_conv$id_type
  attr(raw, "gene_id_conversion_reason") <- id_conv$reason

  if (looks_like_xloc_or_transcript_id(raw$Symbol)) {
    attr(raw, "gene_id_warning") <- paste0(
      "Detected XLOC/MSTRG/TCONS-like IDs in selected symbol column: ", symbol_col,
      ". Gene symbol annotation may be unavailable."
    )
  }

  for (cc in colnames(raw)[-1]) {
    raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
  }

  attr(raw, "symbol_col_used") <- symbol_col
  attr(raw, "gene_id_col_used") <- id_col
  attr(raw, "symbol_col_reason") <- gene_choice$reason

  raw
}



# =========================
# Robust gene symbol column chooser v5.1
# =========================
looks_like_xloc_or_transcript_id <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^XLOC_|^MSTRG\\.|^TCONS_|^CUFF\\.", x, ignore.case = TRUE)) > 0.30
}

looks_like_gene_symbol_vector <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "---"]
  if (length(x) == 0) return(FALSE)
  mean(grepl("^[A-Za-z][A-Za-z0-9_.-]{1,30}$", x)) > 0.50
}

choose_gene_symbol_column <- function(dat) {
  cn <- colnames(dat)

  strong_symbol_cols <- cn[
    grepl("^gene_name$|^gene.symbol$|^gene_symbol$|^symbol$|^hugo$|^hgnc_symbol$|external_gene_name",
          cn, ignore.case = TRUE)
  ]

  if (length(strong_symbol_cols) > 0) {
    for (cc in strong_symbol_cols) {
      if (looks_like_gene_symbol_vector(dat[[as.character(cc)[1]]])) {
        return(list(symbol_col = cc, id_col = cn[1], reason = paste0("preferred symbol column: ", cc)))
      }
    }
  }

  first_col <- cn[1]
  if (looks_like_xloc_or_transcript_id(dat[[first_col]])) {
    candidate_cols <- setdiff(cn[
      grepl("gene_name|gene.symbol|gene_symbol|symbol|hugo|hgnc|nearest_refseq|refseq",
            cn, ignore.case = TRUE)
    ], first_col)

    for (cc in candidate_cols) {
      if (looks_like_gene_symbol_vector(dat[[as.character(cc)[1]]])) {
        return(list(symbol_col = cc, id_col = first_col, reason = paste0("first column is XLOC-like; using annotation column: ", cc)))
      }
    }
  }

  generic_candidates <- cn[
    grepl("symbol|gene_name|gene.name|gene$|Gene$|SYMBOL|external_gene_name|hugo|hgnc|ensembl|id$|ID$",
          cn, ignore.case = TRUE)
  ]

  if (length(generic_candidates) > 0) {
    return(list(symbol_col = generic_candidates[1], id_col = cn[1], reason = paste0("generic candidate: ", generic_candidates[1])))
  }

  non_num <- cn[!sapply(dat, is_numeric_like, min_prop = 0.5)]
  if (length(non_num) > 0) {
    return(list(symbol_col = non_num[1], id_col = cn[1], reason = paste0("first non-numeric column: ", non_num[1])))
  }

  list(symbol_col = cn[1], id_col = cn[1], reason = "fallback first column")
}

is_annotation_col <- function(cn) {
  grepl(
    paste0(
      "^gene$|gene_name|gene.symbol|gene_symbol|symbol|orig_id|",
      "entrez|ensembl|hgnc|hugo|nearest_refseq|refseq|accession|",
      "length|width|gc\\.?content|gc content|description|annotation|",
      "transcript|tx|chromosome|chrom|^chr$|chr\\.|start|end|strand|",
      "biotype|gene.type|gene_biotype|gene.information|gene_information"
    ),
    cn,
    ignore.case = TRUE
  )
}

prepare_raw_from_table_with_symbol <- function(dat, sample_cols = NULL) {
  gene_choice <- choose_gene_symbol_column(dat)

  if (is.null(sample_cols)) {
    sample_cols <- colnames(dat)[sapply(dat, is_numeric_like, min_prop = 0.50)]
  }

  symbol_col <- safe_one_col(gene_choice$symbol_col, colnames(dat), fallback = colnames(dat)[1])
  id_col <- safe_one_col(gene_choice$id_col, colnames(dat), fallback = colnames(dat)[1])

  sample_cols <- safe_char_cols(sample_cols, colnames(dat))
  sample_cols <- setdiff(sample_cols, c(symbol_col, id_col))
  sample_cols <- sample_cols[!is_annotation_col(sample_cols)]

  sample_cols <- safe_char_cols(sample_cols, colnames(dat))
  raw_cols <- safe_char_cols(c(symbol_col, sample_cols), colnames(dat))
  if (length(raw_cols) < 2) return(NULL)
  raw <- dat[, raw_cols, drop = FALSE]
  colnames(raw)[1] <- "Symbol"

  raw$Symbol <- clean_symbol_value(raw$Symbol)

  id_conv2 <- convert_gene_ids_to_symbols_auto(raw$Symbol)
  raw$Symbol <- clean_symbol_value(id_conv2$symbols)
  attr(raw, "gene_id_type") <- id_conv2$id_type
  attr(raw, "gene_id_conversion_reason") <- id_conv2$reason

  for (cc in colnames(raw)[-1]) {
    raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
  }

  attr(raw, "symbol_col_used") <- symbol_col
  attr(raw, "gene_id_col_used") <- id_col
  attr(raw, "symbol_col_reason") <- gene_choice$reason

  if (looks_like_xloc_or_transcript_id(raw$Symbol)) {
    attr(raw, "gene_id_warning") <- paste0(
      "Detected XLOC/MSTRG/TCONS-like IDs in selected symbol column: ", symbol_col,
      ". Gene symbol annotation may be unavailable."
    )
  }

  raw
}


try_ready_matrix_file <- function(f) {
  dat <- read_table_any(f)

  if (is.null(dat) || nrow(dat) < 10 || ncol(dat) < 3) return(NULL)

  dat <- repair_expression_table_header(dat)
  if (is.null(dat) || nrow(dat) < 10 || ncol(dat) < 3) return(NULL)

  dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
  if (ncol(dat) < 3) return(NULL)
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)

  gene_choice <- choose_gene_symbol_column(dat)

  # v7.5 generic block logic:
  # GEO Excel files often contain multiple blocks side-by-side, e.g.
  # [Gene symbol + sample expression columns] [Annotation columns] [Gene information].
  # We do not hard-code a GSE. We keep numeric columns that look like sample/expression
  # columns and remove annotation columns by name and by position after common annotation starts.
  cn <- colnames(dat)
  gene_idx <- match(gene_choice$symbol_col, cn)
  if (is.na(gene_idx)) gene_idx <- 1

  annotation_start <- which(grepl(
    "^annotation$|annotation|entrez|chromosome|chrom|^chr$|start|end|width|strand|transcript|gene.information|gene_information|description|refseq|accession",
    cn,
    ignore.case = TRUE
  ))
  annotation_start <- annotation_start[annotation_start > gene_idx]

  candidate_range <- seq_along(cn)
  if (length(annotation_start) > 0) {
    candidate_range <- seq(gene_idx + 1, min(annotation_start) - 1)
  }
  candidate_range <- candidate_range[candidate_range >= 1 & candidate_range <= length(cn)]

  numeric_cols <- cn[sapply(dat, is_numeric_like, min_prop = 0.50)]
  numeric_cols <- safe_char_cols(numeric_cols, colnames(dat))
  numeric_cols <- numeric_cols[!is_annotation_col(numeric_cols)]

  sample_name_cols <- grep(
    "^GSM|control|healthy|normal|case|sepsis|patient|sample|^S[0-9]+|^C[0-9]+|^H[0-9]+|^K[0-9]+|^d[0-9]+|^hf[0-9]+",
    cn,
    value = TRUE,
    ignore.case = TRUE
  )
  sample_name_cols <- safe_char_cols(sample_name_cols, colnames(dat))
  sample_name_cols <- sample_name_cols[!is_annotation_col(sample_name_cols)]

  # Prefer the contiguous expression block before the first annotation block.
  block_cols <- cn[candidate_range]
  block_cols <- block_cols[block_cols != gene_choice$symbol_col]
  block_cols <- block_cols[!is_annotation_col(block_cols)]
  block_cols <- safe_char_cols(block_cols, colnames(dat))
  block_cols <- block_cols[sapply(block_cols, function(cc) is_numeric_like(dat[[as.character(cc)[1]]], min_prop = 0.50))]

  sample_cols <- unique(c(block_cols, sample_name_cols, numeric_cols))
  sample_cols <- safe_char_cols(sample_cols, colnames(dat))
  sample_cols <- setdiff(sample_cols, c(gene_choice$symbol_col, gene_choice$id_col))
  sample_cols <- sample_cols[!is_annotation_col(sample_cols)]

  # Final content sanity check: expression sample columns should contain many numeric values,
  # and should not be genomic-coordinate-like columns.
  if (length(sample_cols) > 0) {
    sample_cols <- sample_cols[sapply(sample_cols, function(cc) {
      vals <- suppressWarnings(as.numeric(dat[[as.character(cc)[1]]]))
      finite <- vals[is.finite(vals)]
      if (length(finite) < max(10, nrow(dat) * 0.30)) return(FALSE)
      # Annotation coordinate columns are often huge positive integers with almost no zeros.
      int_prop <- mean(abs(finite - round(finite)) < 1e-6, na.rm = TRUE)
      zero_prop <- mean(finite == 0, na.rm = TRUE)
      coord_like <- median(finite, na.rm = TRUE) > 1000 && int_prop > 0.95 && zero_prop < 0.01
      !coord_like
    })]
  }

  if (length(sample_cols) < 2) return(NULL)

  raw <- prepare_raw_from_table_with_symbol(dat, sample_cols = sample_cols)
  if (is.null(raw) || ncol(raw) < 3) return(NULL)

  keep <- c(TRUE, sapply(raw[-1], function(x) sum(!is.na(x)) >= max(10, floor(nrow(raw) * 0.30))))
  raw <- raw[, keep, drop = FALSE]

  if (ncol(raw) < 3) return(NULL)

  raw <- collapse_by_symbol_auto(raw, source_hint = basename(f))

  attr(raw, "symbol_col_used") <- gene_choice$symbol_col
  attr(raw, "gene_id_col_used") <- gene_choice$id_col
  attr(raw, "symbol_col_reason") <- gene_choice$reason

  if (nrow(raw) > 50 && ncol(raw) > 3) {
    return(list(
      raw = raw,
      symbol_col = attr(raw, "symbol_col_used"),
      symbol_col_reason = attr(raw, "symbol_col_reason"),
      gene_id_col = attr(raw, "gene_id_col_used"),
      gene_id_warning = attr(raw, "gene_id_warning"),
      expr_detect = attr(raw, "expr_detect"),
      fast_reader = FALSE,
      ready_matrix_reader = TRUE
    ))
  }

  NULL
}


# =========================
# Robust RAW/supplementary file discovery helper
# =========================


collapse_single_sample_duplicates <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) return(dat)
  if (!all(c("Symbol", "count") %in% colnames(dat))) return(dat)

  # IMPORTANT v6.7:
  # Do NOT convert Ensembl/Entrez IDs to gene symbols inside every single-sample file.
  # Doing mapIds() 80-200 times is slow and floods the console with
  # "select() returned 1:many mapping" messages.
  # Here we only clean IDs and collapse exact duplicate IDs within each file.
  # The final merged matrix is converted to gene symbols ONCE by collapse_by_symbol_auto().
  dat$Symbol <- clean_symbol_value(dat$Symbol)
  dat$count <- suppressWarnings(as.numeric(dat$count))

  dat <- dat[
    !is.na(dat$Symbol) & dat$Symbol != "" &
      !is.na(dat$count),
    ,
    drop = FALSE
  ]

  if (nrow(dat) == 0) return(dat)

  det <- detect_expression_type(data.frame(x = dat$count))
  if (!is.null(det$type) && det$type == "raw_count") {
    dat <- dat %>%
      dplyr::group_by(Symbol) %>%
      dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  } else {
    dat <- dat %>%
      dplyr::group_by(Symbol) %>%
      dplyr::summarise(count = mean(count, na.rm = TRUE), .groups = "drop")
  }

  as.data.frame(dat)
}


deduplicate_expression_files_by_sample <- function(files) {
  if (length(files) == 0) return(files)

  files <- unique(files)

  get_key <- function(f) {
    bn <- basename(f)
    bn <- sub("\\.gz$", "", bn, ignore.case = TRUE)
    bn <- sub("\\.htseq\\.results$", "", bn, ignore.case = TRUE)
    bn <- sub("\\.txt$|\\.tsv$|\\.csv$|\\.count$|\\.counts$|\\.sf$|\\.out$", "", bn, ignore.case = TRUE)
    bn
  }

  keys <- vapply(files, get_key, character(1))

  # Prefer uncompressed files over gz files.
  is_gz <- grepl("\\.gz$", files, ignore.case = TRUE)
  ord <- order(keys, is_gz)  # FALSE before TRUE
  files2 <- files[ord]
  keys2 <- keys[ord]

  files2[!duplicated(keys2)]
}


safe_uncompress_gz_files <- function(files, max_files = 500) {
  gz_files <- files[grepl("\\.gz$", files, ignore.case = TRUE)]
  if (length(gz_files) == 0) return(invisible(NULL))

  gz_files <- head(gz_files, max_files)

  for (gf in gz_files) {
    out_f <- sub("\\.gz$", "", gf, ignore.case = TRUE)
    if (!file.exists(out_f)) {
      try(R.utils::gunzip(gf, destname = out_f, remove = FALSE, overwrite = TRUE), silent = TRUE)
    }
  }

  invisible(NULL)
}


# Detect 10X / single-cell supplementary structure early.
# This app is designed for sample-level bulk RNA-seq / microarray DEG.
# 10X matrix.mtx + barcodes.tsv + features.tsv/genes.tsv is cell-level data,
# so it should not enter the bulk parser.

is_scRNA_like_filename <- function(x) {
  bn <- basename(x)
  grepl(
    paste0(
      "matrix.*\\.mtx(\\.gz)?$|",
      "\\.mtx(\\.gz)?$|",
      "barcodes.*\\.tsv(\\.gz)?$|",
      "features.*\\.tsv(\\.gz)?$|",
      "genes.*\\.tsv(\\.gz)?$|",
      "hashing.*\\.xlsx$|hashing.*\\.xls$|hashing.*\\.csv$|hashing.*\\.tsv$|",
      "cellranger|filtered_feature_bc_matrix|raw_feature_bc_matrix|",
      "ADT|HTO|CITE|cell[_-]?barcode|cell[_-]?ranger"
    ),
    bn,
    ignore.case = TRUE
  )
}

is_10x_single_cell_dir <- function(supp_dir) {
  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (length(files) == 0) return(FALSE)

  bn <- basename(files)

  # Accept GEO names such as GSM5243248_matrix_Patient1.mtx,
  # not only the strict 10X name matrix.mtx.
  has_mtx <- any(grepl("matrix.*\\.mtx(\\.gz)?$|\\.mtx(\\.gz)?$", bn, ignore.case = TRUE))
  has_barcode <- any(grepl("barcodes.*\\.tsv(\\.gz)?$", bn, ignore.case = TRUE))
  has_feature <- any(grepl("features.*\\.tsv(\\.gz)?$|genes.*\\.tsv(\\.gz)?$", bn, ignore.case = TRUE))

  # 10X structure: sparse matrix + cell barcodes + features/genes.
  if (has_mtx && has_barcode && has_feature) return(TRUE)

  # Strong single-cell supplementary pattern even if .mtx was skipped/removed:
  # many barcode/features files or hashing files.
  n_barcode <- sum(grepl("barcodes.*\\.tsv(\\.gz)?$", bn, ignore.case = TRUE))
  n_feature <- sum(grepl("features.*\\.tsv(\\.gz)?$|genes.*\\.tsv(\\.gz)?$", bn, ignore.case = TRUE))
  n_hashing <- sum(grepl("hashing", bn, ignore.case = TRUE))

  if ((n_barcode >= 2 && n_feature >= 2) || n_hashing >= 2) return(TRUE)

  FALSE
}

report_scRNA_files <- function(supp_dir, gse_id = "") {
  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (length(files) == 0) return(invisible(character()))

  sc_files <- basename(files[is_scRNA_like_filename(files)])
  sc_files <- unique(sc_files)

  if (length(sc_files) > 0) {
    message(
      "检测到单细胞/10X相关文件，将在bulk分析中跳过：",
      length(sc_files), " 个。示例：",
      paste(head(sc_files, 10), collapse = ", ")
    )
  }

  invisible(sc_files)
}

stop_if_only_scRNA_no_bulk_candidates <- function(supp_dir, candidate_files, gse_id = "") {
  if (length(candidate_files) == 0 && is_10x_single_cell_dir(supp_dir)) {
    files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    sc_files <- basename(files[is_scRNA_like_filename(files)])
    sc_files <- head(unique(sc_files), 20)

    stop(paste0(
      "检测到10X/单细胞RNA-seq数据：", gse_id, "\n",
      "但没有找到可用于bulk RNA-seq / Microarray分析的表达矩阵或单样本bulk表达/count文件。\n",
      "已跳过matrix.mtx、barcodes.tsv、features.tsv/genes.tsv、hashing_info等单细胞文件。\n",
      "单细胞相关文件示例：\n",
      paste(sc_files, collapse = "\n"), "\n",
      "如果要分析单细胞数据，建议使用Seurat/Scanpy单独流程。"
    ))
  }
  invisible(FALSE)
}

list_expression_candidate_files <- function(supp_dir) {
  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  report_scRNA_files(supp_dir)

  # Some tar files contain nested tar/zip/gz files. Unpack one more layer when obvious.
  nested_archives <- files[grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$|\\.zip$", files, ignore.case = TRUE)]
  if (length(nested_archives) > 0) {
    for (af in nested_archives) {
      if (grepl("\\.zip$", af, ignore.case = TRUE)) {
        try(utils::unzip(af, exdir = dirname(af)), silent = TRUE)
      } else {
        try(utils::untar(af, exdir = dirname(af)), silent = TRUE)
      }
    }
  }

  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)

  # Try to make plain text versions for common gz text files.
  safe_uncompress_gz_files(files)

  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  report_scRNA_files(supp_dir)

  candidate_files <- files[
    grepl(
      "\\.txt$|\\.txt\\.gz$|\\.tsv$|\\.tsv\\.gz$|\\.csv$|\\.csv\\.gz$|\\.xls$|\\.xlsx$|\\.count$|\\.counts$|\\.out$|\\.sf$|htseq\\.results$|htseq\\.results\\.gz$|featureCounts$|gene_counts$|quant\\.sf$",
      files,
      ignore.case = TRUE
    )
  ]

  # Mixed GSE support: skip only obvious 10X/scRNA files, keep bulk-compatible files in the same series.
  candidate_files <- candidate_files[!is_scRNA_like_filename(candidate_files)]

  deduplicate_expression_files_by_sample(candidate_files)
}

debug_supplementary_files_message <- function(supp_dir, max_show = 80) {
  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)

  if (length(files) == 0) {
    return("解压后没有发现任何文件。")
  }

  info <- data.frame(
    file = basename(files),
    path = dirname(files),
    size_MB = round(file.info(files)$size / 1024 / 1024, 3),
    stringsAsFactors = FALSE
  )

  info <- info[order(-info$size_MB), , drop = FALSE]
  show <- head(info, max_show)

  paste0(
    "解压后发现文件数：", length(files), "\n",
    "最大/最可能相关文件预览：\n",
    paste0(show$file, " | ", show$size_MB, " MB", collapse = "\n")
  )
}



# =========================
# Universal bulk single-sample file validator / cohort selector
# =========================
# Purpose:
# - Do NOT hard-code any GSE ID or tissue words such as Whole-blood/PBMC.
# - Skip obvious single-cell support files by filename.
# - For remaining GSM*.tsv/count files, validate by CONTENT STRUCTURE.
# - If one GEO contains multiple sample-level cohorts/modalities, use GEO metadata
#   to keep the largest homogeneous cohort instead of blindly merging all GSM files.

extract_gsm_from_filename <- function(f) {
  m <- regmatches(basename(f), regexpr("GSM[0-9]+", basename(f), ignore.case = TRUE))
  if (length(m) == 0 || is.na(m) || m == "") return(NA_character_)
  toupper(m)
}

normalize_text_signature <- function(x, max_words = 4) {
  x <- tolower(as.character(x))
  x <- gsub("<[^>]+>", " ", x)
  x <- gsub("\\([^)]*\\)", " ", x)
  x <- gsub("gsm[0-9]+", " ", x)
  x <- gsub("[0-9]+", " ", x)
  x <- gsub("[_./:;,+-]+", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (is.na(x) || x == "") return("unknown")
  w <- unlist(strsplit(x, "\\s+"))
  w <- w[w != ""]
  if (length(w) == 0) return("unknown")
  paste(head(w, max_words), collapse = " ")
}

sample_metadata_signature_for_file <- function(f, series_meta) {
  gsm <- extract_gsm_from_filename(f)
  if (is.na(gsm) || is.null(series_meta) || nrow(series_meta) == 0) return("unknown")

  meta <- as.data.frame(series_meta)
  cn <- colnames(meta)
  sample_col <- NULL
  for (cc in c("geo_accession", "sample", "Sample", "accession", "gsm")) {
    hit <- cn[tolower(cn) == tolower(cc)]
    if (length(hit) > 0) { sample_col <- hit[1]; break }
  }
  if (is.null(sample_col)) return("unknown")

  idx <- which(toupper(as.character(meta[[sample_col]])) == gsm)
  if (length(idx) == 0) return("unknown")
  row <- meta[idx[1], , drop = FALSE]

  # Prefer generic sample description fields if present. This is not tissue/GSE specific;
  # it groups samples by the dataset's own metadata wording.
  text_cols_priority <- c("source_name_ch1", "title", "description", "characteristics_ch1")
  text_cols <- cn[tolower(cn) %in% tolower(text_cols_priority)]
  if (length(text_cols) == 0) {
    text_cols <- cn[sapply(row, function(z) is.character(z) || is.factor(z))]
  }
  if (length(text_cols) == 0) return("unknown")

  txt <- paste(as.character(row[, text_cols, drop = TRUE]), collapse = " ")
  normalize_text_signature(txt, max_words = 4)
}

filter_single_files_by_metadata_cohort <- function(files, series_meta, min_cohort_fraction = 0.30) {
  if (length(files) < 4 || is.null(series_meta) || nrow(as.data.frame(series_meta)) == 0) return(files)

  sig <- vapply(files, sample_metadata_signature_for_file, character(1), series_meta = series_meta)
  tab <- sort(table(sig), decreasing = TRUE)
  if (length(tab) <= 1) return(files)

  best_sig <- names(tab)[1]
  best_n <- as.integer(tab[1])
  total_n <- length(files)

  # Only filter when a clearly dominant cohort exists. This avoids damaging normal datasets
  # where metadata signatures are heterogeneous but all files are same modality.
  if (best_sig != "unknown" && best_n >= 2 && best_n / total_n >= min_cohort_fraction) {
    kept <- files[sig == best_sig]
    skipped <- files[sig != best_sig]
    message("检测到多个GSM单样本文件cohort，按metadata自动保留最大同类cohort：", best_sig,
            "；保留 ", length(kept), " 个，跳过 ", length(skipped), " 个。")
    message("被跳过cohort概览：", paste(head(names(tab)[names(tab) != best_sig], 10), collapse = ", "))
    return(kept)
  }

  files
}

looks_like_bulk_gene_vector <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != "" & !grepl("^__", x)]
  if (length(x) < 30) return(FALSE)

  # Accept common gene identifiers/symbols across bulk RNA-seq files.
  # Do NOT require a specific GSE/sample title.
  prop_gene <- mean(
    grepl("^ENSG[0-9]+", x, ignore.case = TRUE) |                 # Ensembl gene
      grepl("^ENS[A-Z]*G[0-9]+", x, ignore.case = TRUE) |          # other Ensembl-like gene IDs
      grepl("^[0-9]+$", x) |                                      # Entrez ID
      grepl("^[A-Za-z][A-Za-z0-9_.-]{1,45}$", x)                  # HUGO-like symbols, IG/TCR symbols
  )

  prop_gene >= 0.40
}

is_probably_bulk_expression_pair <- function(symbols, values, min_rows = 100) {
  symbols <- as.character(symbols)
  symbols <- trimws(symbols)
  values_num <- suppressWarnings(as.numeric(values))

  keep <- !is.na(symbols) & symbols != "" & !grepl("^__", symbols) & !is.na(values_num)
  n_keep <- sum(keep)
  if (n_keep < min_rows) return(FALSE)

  symbols2 <- symbols[keep]
  values2 <- values_num[keep]

  if (!looks_like_bulk_gene_vector(symbols2)) return(FALSE)
  if (length(unique(symbols2)) < max(30, min_rows * 0.40)) return(FALSE)
  if (mean(is.finite(values2)) < 0.95) return(FALSE)

  # Skip obvious annotation/binary metadata columns, but do not reject sparse raw counts.
  u <- length(unique(values2))
  if (u <= 2 && length(values2) > 1000 && mean(values2 == 0, na.rm = TRUE) < 0.98) return(FALSE)

  TRUE
}

score_bulk_expression_pair <- function(symbols, values) {
  symbols <- as.character(symbols)
  symbols <- trimws(symbols)
  values_num <- suppressWarnings(as.numeric(values))
  keep <- !is.na(symbols) & symbols != "" & !grepl("^__", symbols) & !is.na(values_num)
  n_keep <- sum(keep)
  if (n_keep < 30) return(-Inf)

  symbols2 <- symbols[keep]
  values2 <- values_num[keep]

  gene_prop <- mean(
    grepl("^ENSG[0-9]+", symbols2, ignore.case = TRUE) |
      grepl("^ENS[A-Z]*G[0-9]+", symbols2, ignore.case = TRUE) |
      grepl("^[0-9]+$", symbols2) |
      grepl("^[A-Za-z][A-Za-z0-9_.-]{1,45}$", symbols2)
  )
  uniq_prop <- length(unique(symbols2)) / max(1, length(symbols2))
  finite_prop <- mean(is.finite(values2))
  nonneg_bonus <- ifelse(mean(values2 >= 0, na.rm = TRUE) > 0.95, 0.2, 0)

  # Penalize columns that look like genomic coordinates/length more than expression.
  vals_sample <- if (length(values2) > 10000) sample(values2, 10000) else values2
  int_prop <- mean(abs(vals_sample - round(vals_sample)) < 1e-6, na.rm = TRUE)
  zero_prop <- mean(vals_sample == 0, na.rm = TRUE)
  coord_penalty <- ifelse(median(vals_sample, na.rm = TRUE) > 1000 && zero_prop < 0.01 && int_prop > 0.95, 0.4, 0)

  log10(n_keep + 1) + gene_prop * 3 + uniq_prop + finite_prop + nonneg_bonus - coord_penalty
}

read_table_for_single_sample_expression <- function(f) {
  # Try normal fread first. If a featureCounts file has leading comment lines,
  # fread usually handles them poorly; retry with comment.char="#".
  dat <- tryCatch({
    data.table::fread(f, data.table = FALSE, fill = TRUE, check.names = FALSE)
  }, error = function(e) NULL)

  if (is.null(dat) || ncol(dat) < 2) {
    dat <- tryCatch({
      data.table::fread(f, data.table = FALSE, fill = TRUE, check.names = FALSE, comment.char = "#")
    }, error = function(e) NULL)
  }

  # Headerless fallback.
  if (is.null(dat) || nrow(dat) < 30 || ncol(dat) < 2) {
    dat <- tryCatch({
      data.table::fread(f, header = FALSE, data.table = FALSE, fill = TRUE, check.names = FALSE)
    }, error = function(e) NULL)
  }

  dat
}

choose_best_gene_value_columns <- function(dat) {
  if (is.null(dat) || nrow(dat) < 30 || ncol(dat) < 2) return(NULL)

  colnames(dat) <- make.names(colnames(dat), unique = TRUE)
  cn <- colnames(dat)

  # Candidate gene columns: gene-like names first, then non-numeric columns, then first column.
  gene_named <- cn[grepl("gene|symbol|ensembl|feature|id$|ID$|Geneid|Name", cn, ignore.case = TRUE)]
  non_num <- cn[!sapply(dat, is_numeric_like, min_prop = 0.50)]
  gene_candidates <- unique(c(gene_named, non_num, cn[1]))

  # Candidate value columns: numeric columns. Prefer expression/count-ish names.
  numeric_cols <- cn[sapply(dat, is_numeric_like, min_prop = 0.30)]
  numeric_cols <- setdiff(numeric_cols, gene_candidates[1])
  if (length(numeric_cols) == 0) return(NULL)

  priority_value <- numeric_cols[grepl(
    "count|counts|expected_count|TPM|FPKM|RPKM|CPM|NumReads|expression|expr|abundance|read",
    numeric_cols,
    ignore.case = TRUE
  )]
  value_candidates <- unique(c(priority_value, numeric_cols))

  best <- list(score = -Inf, gene_col = NA_character_, value_col = NA_character_)

  for (gc in gene_candidates) {
    for (vc in setdiff(value_candidates, gc)) {
      sc <- score_bulk_expression_pair(dat[[gc]], dat[[vc]])
      # small bonus for likely value column names
      if (grepl("count|counts|expected_count|TPM|FPKM|RPKM|CPM|NumReads|expression|expr|abundance|read", vc, ignore.case = TRUE)) sc <- sc + 0.5
      if (sc > best$score) {
        best <- list(score = sc, gene_col = gc, value_col = vc)
      }
    }
  }

  if (!is.finite(best$score)) return(NULL)
  best
}


# =========================
# Supplementary expression sanity checks
# =========================
# A true per-sample bulk RNA-seq expression/count file is rarely ~1 KB.
# Tiny GSM*.txt files in GEO RAW.tar are often notes, links, QC stubs, or SRA metadata,
# not gene-level expression/count tables. This is a structural filter, not GSE-specific.
is_too_tiny_for_expression_file <- function(f, min_size_bytes = 5000) {
  sz <- suppressWarnings(file.info(f)$size)
  if (length(sz) == 0 || is.na(sz)) return(FALSE)
  sz < min_size_bytes
}

preview_text_file_lines <- function(f, n = 5) {
  out <- tryCatch(readLines(f, n = n, warn = FALSE), error = function(e) character())
  out <- gsub("[\r\n\t]+", " ", out)
  out <- out[nchar(out) > 0]
  paste(head(out, n), collapse = " | ")
}


# =========================
# Bottom-level failure diagnosis helpers
# =========================
# These functions do not decide by GSE ID, tissue name, or disease name.
# They classify failed supplementary files by file structure/content so the app can
# explain why a GEO cannot be parsed as a bulk gene-level expression dataset.
classify_non_expression_file <- function(f) {
  bn <- safe_utf8(basename(f))
  bn_low <- tolower(bn)
  sz <- suppressWarnings(file.info(f)$size)
  if (length(sz) == 0 || is.na(sz)) sz <- NA_real_

  if (grepl("\\.xlsx$|\\.xls$", bn_low, ignore.case = TRUE)) {
    return("Excel workbook not accepted by expression parsers")
  }
  if (grepl("\\.zip$|\\.tar$|\\.tar\\.gz$|\\.tgz$", bn_low, ignore.case = TRUE)) {
    return("archive file")
  }
  if (grepl("matrix\\.mtx|barcodes|features|genes\\.tsv|hashing|hto|adt|cell[_-]?barcode|barcode[_-]?cell", bn_low)) {
    return("single-cell support file")
  }
  if (grepl("fastq|fq\\.gz|\\.sra$|\\.bam$|\\.sam$", bn_low)) {
    return("raw sequencing file")
  }
  if (grepl("multiqc|fastqc|qc|quality|summary|report|log", bn_low)) {
    return("QC/report file")
  }

  lines <- safe_read_preview_lines(f, n = 30)
  txt_low <- tolower(paste(lines, collapse = " "))

  if (grepl("left reads|right reads|input *:|mapped reads|uniquely mapped|overall alignment rate|alignment rate|concordantly exactly|star|hisat|bowtie|subread", txt_low)) {
    return("alignment/mapping statistics file")
  }
  if (grepl("fastqc|per base sequence quality|sequence duplication levels|adapter content|overrepresented sequences", txt_low)) {
    return("QC/report file")
  }
  if (grepl("run accession|experiment accession|sample accession|sra|srx|srr|ena|biosample", txt_low)) {
    return("SRA metadata/stub file")
  }
  if (!is.na(sz) && sz < 5000) {
    return("tiny non-expression text/stub")
  }
  "unrecognized non-expression file"
}

diagnose_failed_supplementary_files <- function(files, max_show = 25) {
  if (length(files) == 0) {
    return("未发现可检查的supplementary候选文件。")
  }

  files <- unique(files)
  classes <- vapply(files, classify_non_expression_file, character(1))
  tab <- sort(table(classes), decreasing = TRUE)

  info <- data.frame(
    file = safe_utf8(basename(files)),
    class = safe_utf8(classes),
    size_KB = round(suppressWarnings(file.info(files)$size) / 1024, 2),
    preview = vapply(files, function(ff) paste(safe_read_preview_lines(ff, n = 2), collapse = " | "), character(1)),
    stringsAsFactors = FALSE
  )
  info <- info[order(info$class, info$size_KB), , drop = FALSE]

  show <- head(info, max_show)
  paste0(
    "自动诊断结果：\n",
    paste0(names(tab), ": ", as.integer(tab), " 个", collapse = "\n"),
    "\n\n代表性文件预览：\n",
    paste0(show$file, " | ", show$class, " | ", show$size_KB, " KB | ", show$preview, collapse = "\n"),
    if (nrow(info) > max_show) paste0("\n... 共 ", nrow(info), " 个候选文件，仅显示前 ", max_show, " 个") else "",
    "\n\n未检测到可直接用于DESeq2/limma的gene-level count/FPKM/TPM表达矩阵。",
    "\n如果主要是alignment/mapping statistics或SRA metadata/stub，说明GEO没有提供processed表达矩阵，需要从SRA FASTQ/BAM重新定量，或换有processed matrix的GEO。"
  )
}

summarize_rejected_candidate_files <- function(files, max_show = 30) {
  if (length(files) == 0) return("无候选文件。")
  info <- data.frame(
    file = safe_utf8(basename(files)),
    size_KB = round(suppressWarnings(file.info(files)$size) / 1024, 2),
    preview = vapply(files, function(ff) paste(safe_read_preview_lines(ff, n = 2), collapse = " | "), character(1)),
    stringsAsFactors = FALSE
  )
  info <- info[order(info$size_KB), , drop = FALSE]
  show <- head(info, max_show)
  paste0(
    paste0(show$file, " | ", show$size_KB, " KB | ", show$preview, collapse = "\n"),
    if (nrow(info) > max_show) paste0("\n... 共 ", nrow(info), " 个候选文件") else ""
  )
}

# Wider parser for single-sample gene-count/expression files.
# Bottom-level principle:
#   - Never use GSE-specific words such as PBMC/whole blood to decide.
#   - Decide by file structure: one gene-like column + one numeric expression/count column.
#   - Skip cell/barcode/feature/hash metadata by filename and by content.
read_one_single_sample_expr_file_robust <- function(f, min_rows = 100) {
  bn <- basename(f)

  # Skip files that are structurally too small to be gene-level expression/count tables.
  # This prevents wasting time on GEO stub/metadata txt files.
  if (is_too_tiny_for_expression_file(f)) {
    return(NULL)
  }

  # Explicitly skip single-cell support/metadata files. This is structural, not GSE-specific.
  if (grepl("matrix\\.mtx|barcodes|features|genes\\.tsv|hashing|hto|adt|cell[_-]?barcode|barcode[_-]?cell",
            bn,
            ignore.case = TRUE)) {
    return(NULL)
  }

  dat <- read_table_for_single_sample_expression(f)
  if (is.null(dat) || nrow(dat) < 30 || ncol(dat) < 2) return(NULL)

  # Remove completely empty columns.
  dat <- dat[, colSums(!is.na(dat)) > 0, drop = FALSE]
  if (nrow(dat) < 30 || ncol(dat) < 2) return(NULL)

  # Fast path for true 2-column files, including headerless htseq-like outputs.
  if (ncol(dat) == 2) {
    dat2 <- dat[, 1:2, drop = FALSE]
    colnames(dat2) <- c("Symbol", "count")
    dat2$Symbol <- as.character(dat2$Symbol)
    dat2$count <- suppressWarnings(as.numeric(dat2$count))
    dat2 <- dat2[
      !is.na(dat2$Symbol) & dat2$Symbol != "" &
        !grepl("^__|^N_", dat2$Symbol) &
        !is.na(dat2$count),
      ,
      drop = FALSE
    ]

    if (is_probably_bulk_expression_pair(dat2$Symbol, dat2$count, min_rows = min_rows)) {
      return(collapse_single_sample_duplicates(dat2))
    }
    return(NULL)
  }

  best <- choose_best_gene_value_columns(dat)
  if (is.null(best)) return(NULL)

  out <- dat[, c(best$gene_col, best$value_col), drop = FALSE]
  colnames(out) <- c("Symbol", "count")
  out$Symbol <- as.character(out$Symbol)
  out$count <- suppressWarnings(as.numeric(out$count))

  out <- out[
    !is.na(out$Symbol) & out$Symbol != "" &
      !grepl("^__|^N_", out$Symbol) &
      !is.na(out$count),
    ,
    drop = FALSE
  ]

  if (!is_probably_bulk_expression_pair(out$Symbol, out$count, min_rows = min_rows)) {
    return(NULL)
  }

  attr(out, "gene_col_used") <- best$gene_col
  attr(out, "value_col_used") <- best$value_col
  collapse_single_sample_duplicates(out)
}

merge_single_sample_expr_files_robust <- function(files, max_files = 5000) {
  files <- deduplicate_expression_files_by_sample(files)
  files <- files[!duplicated(normalizePath(files, winslash = "/", mustWork = FALSE))]
  tiny_files <- files[vapply(files, is_too_tiny_for_expression_file, logical(1))]
  if (length(tiny_files) > 0) {
    tiny_classes <- sort(table(vapply(tiny_files, classify_non_expression_file, character(1))), decreasing = TRUE)
    message("跳过过小、通常不可能是gene-level表达矩阵的文件：", length(tiny_files), " 个")
    message("过小文件类别：", paste0(names(tiny_classes), "=", as.integer(tiny_classes), collapse = "; "))
    message("过小文件示例：", paste(head(basename(tiny_files), 10), collapse = ", "))
  }
  files <- files[!vapply(files, is_too_tiny_for_expression_file, logical(1))]
  files <- head(files, max_files)

  count_list <- list()
  failed <- 0

  for (f in files) {
    dat <- read_one_single_sample_expr_file_robust(f)

    if (is.null(dat)) {
      failed <- failed + 1
      next
    }

    sample_name <- sample_name_from_file(f)
    sample_name <- make.names(sample_name, unique = FALSE)

    # Avoid duplicate sample names after make.names()
    if (sample_name %in% names(count_list)) {
      sample_name <- make.unique(c(names(count_list), sample_name))[length(count_list) + 1]
    }

    colnames(dat)[2] <- sample_name
    count_list[[sample_name]] <- dat
  }

  message("单样本文件识别成功：", length(count_list), "；失败/跳过：", failed, "；每个文件内重复Symbol已先聚合")
  if (length(count_list) > 0) {
    message("成功识别样本示例：", paste(head(names(count_list), 10), collapse = ", "))
  }

  if (length(count_list) < 2) return(NULL)

  message("正在合并单样本矩阵：", length(count_list), " 个样本")

  # Each single-sample table has unique Symbol now; merge is one-to-one/one-to-zero.
  raw <- Reduce(function(x, y) merge(x, y, by = "Symbol", all = TRUE), count_list)
  raw[is.na(raw)] <- 0

  message("单样本矩阵合并完成：", nrow(raw), " genes x ", ncol(raw) - 1, " samples")

  collapse_by_symbol_auto(raw, forced_method = "sum", source_hint = paste(basename(files), collapse = " "))
}

# Detect 10X-style matrix directories.
try_read_10x_matrix <- function(supp_dir) {
  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)

  mtx <- files[grepl("matrix\\.mtx(\\.gz)?$", basename(files), ignore.case = TRUE)]
  features <- files[grepl("features\\.tsv(\\.gz)?$|genes\\.tsv(\\.gz)?$", basename(files), ignore.case = TRUE)]
  barcodes <- files[grepl("barcodes\\.tsv(\\.gz)?$", basename(files), ignore.case = TRUE)]

  if (length(mtx) == 0 || length(features) == 0 || length(barcodes) == 0) return(NULL)

  # Keep this lightweight and optional.
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    install.packages("Matrix", repos = "https://cloud.r-project.org")
  }

  mtx_f <- mtx[1]
  feat_f <- features[1]
  bc_f <- barcodes[1]

  mat <- tryCatch(Matrix::readMM(mtx_f), error = function(e) NULL)
  if (is.null(mat)) return(NULL)

  feat <- tryCatch(data.table::fread(feat_f, header = FALSE, data.table = FALSE), error = function(e) NULL)
  bc <- tryCatch(data.table::fread(bc_f, header = FALSE, data.table = FALSE), error = function(e) NULL)

  if (is.null(feat) || is.null(bc)) return(NULL)

  genes <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
  samples <- bc[[1]]

  if (length(genes) != nrow(mat) || length(samples) != ncol(mat)) return(NULL)

  rownames(mat) <- make.names(as.character(genes), unique = TRUE)
  colnames(mat) <- make.names(as.character(samples), unique = TRUE)

  dense <- as.data.frame(as.matrix(mat))
  dense$Symbol <- rownames(dense)
  dense <- dense[, c("Symbol", setdiff(colnames(dense), "Symbol"))]

  raw <- collapse_by_symbol_auto(dense, forced_method = "sum", source_hint = paste("10X raw counts", basename(mtx_f)))

  if (nrow(raw) > 100 && ncol(raw) > 3) {
    return(list(
      raw = raw,
      symbol_col = "10X features.tsv",
      expr_detect = attr(raw, "expr_detect"),
      fast_reader = TRUE,
      ultra_fast_reader = FALSE,
      reader_note = "10X matrix.mtx/features/barcodes"
    ))
  }

  NULL
}



# =========================
# Processed-first supplementary download helper
# =========================
get_geo_supplementary_urls <- function(gse_id) {
  prefix <- sub("[0-9]{3}$", "nnn", gse_id)
  base_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", prefix, "/", gse_id, "/suppl/")

  html <- tryCatch(readLines(base_url, warn = FALSE), error = function(e) NULL)
  if (is.null(html)) return(data.frame(file = character(), url = character(), stringsAsFactors = FALSE))

  href <- regmatches(html, gregexpr('href="[^"]+"', html))
  href <- unlist(href)
  href <- gsub('href="|"$', "", href)
  href <- href[!grepl("^/|^\\?|^\\.\\.", href)]
  href <- href[href != ""]

  files <- basename(URLdecode(href))
  urls <- paste0(base_url, href)

  data.frame(file = files, url = urls, stringsAsFactors = FALSE)
}

# =========================
# v7.7: NCBI-generated RNA-seq matrix discovery
# =========================
# Newer GEO pages may provide standardized NCBI-generated matrices separately from
# submitter supplementary files, e.g.:
#   GSE*_raw_counts_GRCh38*_NCBI.tsv.gz
#   GSE*_norm_counts_FPKM*_NCBI.tsv.gz
#   GSE*_norm_counts_TPM*_NCBI.tsv.gz
# These are usually more standard for bulk DEG than submitter XLSX tables or RAW.tar.
# This is a generic GEO Download-page parser, not GSE-specific.
extract_filename_from_geo_href <- function(href) {
  href_dec <- URLdecode(href)
  m <- regmatches(href_dec, regexpr("(?<=file=)[^&]+", href_dec, perl = TRUE))
  if (length(m) > 0 && !is.na(m) && nzchar(m)) return(basename(m))
  basename(href_dec)
}

make_absolute_geo_url <- function(href, base = "https://www.ncbi.nlm.nih.gov") {
  href <- gsub("&amp;", "&", href)
  href <- trimws(href)

  # Some NCBI GEO download-page links are returned as ftp://ftp.ncbi.nlm.nih.gov/...
  # The old code only recognized http/https, so it incorrectly produced:
  # https://www.ncbi.nlm.nih.gov/ftp://ftp.ncbi.nlm.nih.gov/...
  # Convert NCBI ftp links to the equivalent HTTPS FTP-host URL.
  href <- sub("^/+ftp://", "ftp://", href, ignore.case = TRUE)
  if (grepl("^ftp://", href, ignore.case = TRUE)) {
    return(sub("^ftp://", "https://", href, ignore.case = TRUE))
  }

  if (grepl("^https?://", href, ignore.case = TRUE)) return(href)
  if (grepl("^//", href)) return(paste0("https:", href))
  if (grepl("^/", href)) return(paste0(base, href))
  paste0(base, "/", href)
}

is_ncbi_generated_matrix_filename <- function(x) {
  x <- basename(x)
  has_table_ext <- grepl("\\.tsv(\\.gz)?$|\\.txt(\\.gz)?$|\\.csv(\\.gz)?$", x, ignore.case = TRUE)
  has_ncbi_matrix <- grepl(
    "raw[_-]?counts?|norm[_-]?counts?|normalized[_-]?counts?|FPKM|TPM|RPKM|CPM|NCBI",
    x,
    ignore.case = TRUE
  )
  is_annotation_only <- grepl("annot|annotation|gene[_-]?annotation", x, ignore.case = TRUE)
  has_table_ext & has_ncbi_matrix & !is_annotation_only
}

ncbi_generated_priority <- function(file) {
  f <- basename(file)
  if (grepl("raw[_-]?counts?", f, ignore.case = TRUE)) return(1L)   # best for DESeq2
  if (grepl("TPM", f, ignore.case = TRUE)) return(2L)
  if (grepl("FPKM", f, ignore.case = TRUE)) return(3L)
  if (grepl("norm|normalized", f, ignore.case = TRUE)) return(4L)
  9L
}

# v7.8: smart processed-first priority.
# Principle: any organized gene x sample matrix should be tried before RAW archives.
# Among processed matrices, raw count matrices are best for DESeq2, then TPM/FPKM/normalized
# matrices for limma. This is file-structure based, not GSE-specific.
smart_processed_priority <- function(file, source = "") {
  f <- basename(file)
  src <- tolower(as.character(source))

  # NCBI-generated matrices are usually standardized gene x sample matrices.
  src_bonus <- ifelse(grepl("ncbi", src, ignore.case = TRUE), 0L, 10L)

  if (grepl("raw[_-]?counts?|read[_-]?counts?|count[_-]?matrix|counts[_-]?matrix|gene[_-]?counts", f, ignore.case = TRUE)) {
    return(src_bonus + 1L)
  }
  if (grepl("TPM", f, ignore.case = TRUE)) return(src_bonus + 2L)
  if (grepl("FPKM", f, ignore.case = TRUE)) return(src_bonus + 3L)
  if (grepl("RPKM|CPM|norm[_-]?counts?|normalized|normalised", f, ignore.case = TRUE)) return(src_bonus + 4L)
  if (grepl("expression|expr|matrix|processed|mRNA|RNA[_-]?seq|RNAseq|transcriptome|all.*sample|submitted", f, ignore.case = TRUE)) {
    return(src_bonus + 5L)
  }
  src_bonus + 9L
}

# v10.6: method-aware data source priority.
# Principle:
#   Auto / Force DESeq2: prefer raw count matrices -> DESeq2.
#   Force limma: re-scan available online/local files and prefer TPM/FPKM/normalized matrices.
#   If no normalized matrix exists, limma can still fall back to log2(CPM + 1) from raw counts.
method_aware_source_priority <- function(file, source = "", method_preference = "auto") {
  f <- basename(file)
  src <- tolower(as.character(source))
  method_preference <- ifelse(is.null(method_preference) || is.na(method_preference), "auto", method_preference)

  is_submitter <- grepl("submitter|local", src, ignore.case = TRUE)
  is_ncbi <- grepl("ncbi", src, ignore.case = TRUE)

  is_raw_count <- grepl("raw[_-]?counts?|read[_-]?counts?|count[_-]?matrix|counts[_-]?matrix|gene[_-]?counts|featureCounts|htseq", f, ignore.case = TRUE) &&
    !grepl("norm|normalized|normalised|FPKM|TPM|RPKM|RPM|CPM", f, ignore.case = TRUE)
  is_tpm <- grepl("TPM", f, ignore.case = TRUE)
  is_fpkm <- grepl("FPKM", f, ignore.case = TRUE)
  is_rpm <- grepl("RPM", f, ignore.case = TRUE)
  is_rpkm_cpm <- grepl("RPKM|CPM", f, ignore.case = TRUE)
  is_norm <- grepl("RPKM|RPM|CPM|norm[_-]?counts?|normalized|normalised", f, ignore.case = TRUE)
  is_processed_expr <- grepl("expression|expr|matrix|processed|mRNA|RNA[_-]?seq|RNAseq|transcriptome|all.*sample|submitted", f, ignore.case = TRUE)

  # v10.6 IMPORTANT:
  # For limma, prefer submitter-provided processed matrices (RPM/TPM/FPKM/normalized)
  # over NCBI-generated TPM/FPKM matrices when they cover more selected samples.
  # Final selection is coverage-aware at DEG time, not GSE-specific.
  if (method_preference == "limma") {
    source_bonus <- ifelse(is_submitter, 0L, ifelse(is_ncbi, 20L, 10L))
    if (is_rpm) return(source_bonus + 1L)
    if (is_tpm) return(source_bonus + 2L)
    if (is_fpkm) return(source_bonus + 3L)
    if (is_rpkm_cpm) return(source_bonus + 4L)
    if (is_norm) return(source_bonus + 5L)
    if (is_processed_expr) return(source_bonus + 6L)
    if (is_raw_count) return(source_bonus + 80L)  # only last-resort fallback for limma
    return(source_bonus + 99L)
  }

  if (method_preference == "DESeq2") {
    source_bonus <- ifelse(is_ncbi, 0L, 10L)
    if (is_raw_count) return(source_bonus + 1L)
    return(source_bonus + 99L)
  }

  # auto: keep old count-first behavior for broad bulk RNA-seq compatibility.
  # The user can Force limma to select author-provided RPM/normalized matrices.
  source_bonus <- ifelse(is_ncbi, 0L, 10L)
  if (is_raw_count) return(source_bonus + 1L)
  if (is_tpm) return(source_bonus + 2L)
  if (is_fpkm) return(source_bonus + 3L)
  if (is_rpm) return(source_bonus + 4L)
  if (is_norm) return(source_bonus + 5L)
  if (is_processed_expr) return(source_bonus + 6L)
  source_bonus + 90L
}


# Generic parser for submitter normalized Excel matrices.
# These files may use GEO accessions, sample titles, or study-specific sample IDs.
# We keep the true expression column names; choose_sample_mapping() later matches them to metadata.
try_submitter_rpm_excel_matrix <- function(f) {
  bn <- safe_utf8(basename(f))
  if (!grepl("\\.xlsx$|\\.xls$", bn, ignore.case = TRUE)) return(NULL)
  if (!grepl("RPM|RPKM|FPKM|TPM|expression|processed|normalized|normalised|count", bn, ignore.case = TRUE)) return(NULL)

  message("submitter normalized Excel候选读取模式: ", bn)

  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) "Sheet1")
  sheets <- safe_utf8(sheets)
  best_msg <- character()

  parse_one_excel_sheet <- function(sheet_name) {
    dat0 <- tryCatch({
      as.data.frame(readxl::read_excel(f, sheet = sheet_name, col_names = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE)
    }, error = function(e) {
      best_msg <<- c(best_msg, paste0("sheet ", sheet_name, " read failed: ", e$message))
      NULL
    })
    if (is.null(dat0) || nrow(dat0) < 20 || ncol(dat0) < 3) return(NULL)

    # Try several possible header rows. GEO Excel files often have section headers
    # above the real matrix header. This is generic and not GSE-specific.
    header_candidates <- seq_len(min(12, nrow(dat0) - 1))

    for (header_idx in header_candidates) {
      hdr <- safe_utf8(unlist(dat0[header_idx, ], use.names = FALSE))
      hdr[is.na(hdr) | trimws(hdr) == ""] <- paste0("V", which(is.na(hdr) | trimws(hdr) == ""))
      hdr <- safe_make_names(hdr, unique = TRUE)

      dat <- dat0[-seq_len(header_idx), , drop = FALSE]
      colnames(dat) <- hdr
      dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
      if (nrow(dat) < 50 || ncol(dat) < 3) next

      colnames(dat) <- safe_make_names(colnames(dat), unique = TRUE)
      cn <- colnames(dat)

      gene_choice <- choose_gene_symbol_column(dat)
      symbol_col <- safe_one_col(gene_choice$symbol_col, colnames(dat))
      if (is.na(symbol_col) || symbol_col == "") next

      # Prefer true symbol/name columns if present; otherwise keep GeneID/Ensembl/Entrez
      # and let convert_gene_ids_to_symbols_auto() handle conversion later.
      symbol_candidates <- cn[grepl("^symbol$|gene.?symbol|gene.?name|external_gene_name|hgnc", cn, ignore.case = TRUE)]
      if (length(symbol_candidates) > 0) {
        for (cc in safe_char_cols(symbol_candidates, colnames(dat))) {
          vv <- safe_utf8(dat[[as.character(cc)[1]]])
          if (looks_like_gene_symbol_vector(vv)) { symbol_col <- cc; break }
        }
      }

      numeric_cols <- cn[sapply(dat, is_numeric_like, min_prop = 0.35)]
      sample_cols <- setdiff(numeric_cols, symbol_col)
      sample_cols <- sample_cols[!is_annotation_col(sample_cols)]
      sample_cols <- safe_char_cols(sample_cols, colnames(dat))

      # Avoid genomic coordinate/length columns by content.
      if (length(sample_cols) > 0) {
        sample_cols <- sample_cols[sapply(sample_cols, function(cc) {
          vals <- suppressWarnings(as.numeric(dat[[as.character(cc)[1]]]))
          finite <- vals[is.finite(vals)]
          if (length(finite) < max(20, nrow(dat) * 0.20)) return(FALSE)
          int_prop <- mean(abs(finite - round(finite)) < 1e-6, na.rm = TRUE)
          zero_prop <- mean(finite == 0, na.rm = TRUE)
          coord_like <- median(finite, na.rm = TRUE) > 1000 && int_prop > 0.95 && zero_prop < 0.01
          !coord_like
        })]
      }

      # v11.7.8: do NOT require 20 samples. Many GEO validation cohorts have 4-12 samples.
      if (length(sample_cols) < 2) {
        best_msg <<- c(best_msg, paste0("sheet ", sheet_name, " header row ", header_idx, ": sample columns <2; n=", length(sample_cols)))
        next
      }

      raw_cols <- safe_char_cols(c(symbol_col, sample_cols), colnames(dat))
      if (length(raw_cols) < 3) next

      raw <- dat[, raw_cols, drop = FALSE]
      colnames(raw)[1] <- "Symbol"
      raw$Symbol <- clean_symbol_value(safe_utf8(raw$Symbol))
      raw <- raw[!is.na(raw$Symbol) & raw$Symbol != "" & raw$Symbol != "NA", , drop = FALSE]
      if (nrow(raw) < 50) next

      for (cc in colnames(raw)[-1]) {
        raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
      }
      good_cols <- colnames(raw)[-1][sapply(raw[-1], function(x) mean(!is.na(x)) >= 0.35)]
      good_cols <- safe_char_cols(good_cols, colnames(raw))
      if (length(good_cols) < 2) next
      raw <- raw[, c("Symbol", good_cols), drop = FALSE]
      # Let the universal ID/type detector decide raw_count vs normalized/log.
      raw2 <- collapse_by_symbol_auto(raw, source_hint = paste(basename(f), sheet_name))
      if (nrow(raw2) < 50 || ncol(raw2) < 3) next

      det <- attr(raw2, "expr_detect")
      attr(raw2, "symbol_col_used") <- symbol_col
      attr(raw2, "symbol_col_reason") <- paste0("generic multi-sheet Excel parser; sheet=", sheet_name, "; header row=", header_idx)

      return(list(
        raw = raw2,
        symbol_col = symbol_col,
        symbol_col_reason = attr(raw2, "symbol_col_reason"),
        gene_id_col = symbol_col,
        gene_id_type = attr(raw2, "gene_id_type"),
        gene_id_conversion_reason = attr(raw2, "gene_id_conversion_reason"),
        gene_id_converted_n = attr(raw2, "gene_id_converted_n"),
        expr_detect = det,
        fast_reader = FALSE,
        submitter_rpm_excel_reader = TRUE
      ))
    }
    NULL
  }

  for (sh in sheets) {
    res <- parse_one_excel_sheet(sh)
    if (!is.null(res)) return(res)
  }

  if (length(best_msg) > 0) {
    message("Excel解析未接受该文件。最近原因：", paste(tail(best_msg, 5), collapse = " | "))
  } else {
    message("Excel解析未接受该文件。")
  }
  NULL
}


# =========================
# v11.7.5 robust plain TXT/TSV gene x sample matrix reader
# =========================
# Fixes GEO supplementary files such as GSE296830_messengerrna.txt where
# Series Matrix/GPL has no gene symbol and the true expression matrix is a
# plain submitted txt table. The parser is generic: it does not hard-code GSE,
# and it protects all column subscripts from becoming list-type objects.
try_plain_gene_sample_text_matrix <- function(f) {
  bn <- basename(f)

  if (!grepl("messenger|mRNA|rna|expression|expr|matrix|genes|gene|normalized|FPKM|TPM|RPM|RPKM|CPM|count",
             bn, ignore.case = TRUE)) {
    return(NULL)
  }

  message("稳健TXT/TSV大矩阵读取模式: ", bn)

  dat <- tryCatch({
    data.table::fread(
      f,
      data.table = FALSE,
      fill = TRUE,
      check.names = FALSE,
      showProgress = FALSE
    )
  }, error = function(e) {
    message("稳健TXT/TSV读取失败: ", bn, " | ", e$message)
    NULL
  })

  if (is.null(dat) || nrow(dat) < 100 || ncol(dat) < 3) return(NULL)

  # Remove empty columns and normalize column names.
  dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
  if (nrow(dat) < 100 || ncol(dat) < 3) return(NULL)
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)

  # Some submitted TXT files are read with bad header. Try header repair once.
  dat <- repair_expression_table_header(dat)
  if (is.null(dat) || nrow(dat) < 100 || ncol(dat) < 3) return(NULL)
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)

  gene_choice <- choose_gene_symbol_column(dat)
  symbol_col <- as.character(unlist(gene_choice$symbol_col, use.names = FALSE))[1]
  id_col <- as.character(unlist(gene_choice$id_col, use.names = FALSE))[1]
  if (is.na(symbol_col) || !(symbol_col %in% colnames(dat))) return(NULL)

  cn <- colnames(dat)

  # Strong sample-name candidates by column names.
  sample_name_cols <- grep(
    "^GSM[0-9]+|^SRR[0-9]+|^ERR[0-9]+|^SAMN[0-9]+|^sample|^patient|^control|^healthy|^normal|^case|^sepsis|^S[0-9]+|^C[0-9]+|^H[0-9]+|^P[0-9]+",
    cn,
    value = TRUE,
    ignore.case = TRUE
  )

  # Numeric/expression-like columns by content.
  numeric_cols <- cn[sapply(cn, function(cc) is_numeric_like(dat[[as.character(cc)[1]]], min_prop = 0.50))]

  sample_cols <- unique(c(sample_name_cols, numeric_cols))
  sample_cols <- safe_char_cols(sample_cols, colnames(dat))
  sample_cols <- setdiff(sample_cols, c(symbol_col, id_col))
  sample_cols <- sample_cols[!is_annotation_col(sample_cols)]

  # Content sanity check: keep columns with enough numeric entries and remove coordinate-like annotation columns.
  if (length(sample_cols) > 0) {
    sample_cols <- sample_cols[sapply(sample_cols, function(cc) {
      cc <- as.character(cc)[1]
      vals <- suppressWarnings(as.numeric(dat[[cc]]))
      finite <- vals[is.finite(vals)]
      if (length(finite) < max(50, nrow(dat) * 0.30)) return(FALSE)
      int_prop <- mean(abs(finite - round(finite)) < 1e-6, na.rm = TRUE)
      zero_prop <- mean(finite == 0, na.rm = TRUE)
      coord_like <- median(finite, na.rm = TRUE) > 1000 && int_prop > 0.95 && zero_prop < 0.01
      !coord_like
    })]
  }

  if (length(sample_cols) < 2) return(NULL)

  raw_cols <- safe_char_cols(c(symbol_col, sample_cols), colnames(dat))
  if (length(raw_cols) < 3) return(NULL)

  raw <- dat[, raw_cols, drop = FALSE]
  colnames(raw)[1] <- "Symbol"
  raw$Symbol <- clean_symbol_value(raw$Symbol)

  bad_symbol_rows <- is.na(raw$Symbol) |
    raw$Symbol == "" |
    raw$Symbol == "NA" |
    grepl("SOFA|SOFAVALUE|GROUP|CLASS|PHENOTYPE|DISEASE|STATUS|OUTCOME|MORT|AGE|SEX|GENDER|RACE|BATCH|SAMPLE|DESCRIPTION|ANNOTATION",
          raw$Symbol,
          ignore.case = TRUE)
  raw <- raw[!bad_symbol_rows, , drop = FALSE]
  if (nrow(raw) < 100) return(NULL)

  for (cc in colnames(raw)[-1]) {
    cc <- as.character(cc)[1]
    raw[[cc]] <- suppressWarnings(as.numeric(raw[[cc]]))
  }

  good_cols <- colnames(raw)[-1][sapply(raw[-1], function(x) mean(!is.na(x)) >= 0.50)]
  good_cols <- safe_char_cols(good_cols, colnames(raw))
  if (length(good_cols) < 2) return(NULL)

  raw <- raw[, c("Symbol", good_cols), drop = FALSE]
  keep <- rowSums(!is.na(raw[, good_cols, drop = FALSE])) >= max(2, floor(length(good_cols) * 0.50))
  raw <- raw[keep, , drop = FALSE]
  if (nrow(raw) < 100 || ncol(raw) < 3) return(NULL)
  id_conv <- convert_gene_ids_to_symbols_auto(raw$Symbol)
  raw$Symbol <- clean_symbol_value(id_conv$symbols)
  raw <- raw[!is.na(raw$Symbol) & raw$Symbol != "" & raw$Symbol != "NA", , drop = FALSE]

  # Detect data type on a small subset for speed.
  sub_rows <- if (nrow(raw) > 3000) sample(seq_len(nrow(raw)), 3000) else seq_len(nrow(raw))
  sub_cols <- if (length(good_cols) > 80) sample(good_cols, 80) else good_cols
  det <- detect_expression_type(as.data.frame(raw[sub_rows, sub_cols, drop = FALSE]), source_hint = bn)
  collapse_method <- det$collapse_method

  raw_dt <- data.table::as.data.table(raw)
  if (collapse_method == "sum") {
    raw_dt <- raw_dt[, lapply(.SD, sum, na.rm = TRUE), by = Symbol, .SDcols = good_cols]
  } else {
    raw_dt <- raw_dt[, lapply(.SD, mean, na.rm = TRUE), by = Symbol, .SDcols = good_cols]
  }
  raw2 <- as.data.frame(raw_dt)

  attr(raw2, "expr_detect") <- det
  attr(raw2, "collapse_method") <- collapse_method
  attr(raw2, "symbol_col_used") <- symbol_col
  attr(raw2, "symbol_col_reason") <- paste0("robust plain TXT/TSV gene x sample matrix parser; ", gene_choice$reason)
  attr(raw2, "gene_id_col_used") <- id_col
  attr(raw2, "gene_id_type") <- id_conv$id_type
  attr(raw2, "gene_id_conversion_reason") <- id_conv$reason
  attr(raw2, "gene_id_converted_n") <- id_conv$converted_n

  if (nrow(raw2) > 100 && ncol(raw2) > 3) {
    return(list(
      raw = raw2,
      symbol_col = attr(raw2, "symbol_col_used"),
      symbol_col_reason = attr(raw2, "symbol_col_reason"),
      gene_id_col = attr(raw2, "gene_id_col_used"),
      gene_id_type = attr(raw2, "gene_id_type"),
      gene_id_conversion_reason = attr(raw2, "gene_id_conversion_reason"),
      gene_id_converted_n = attr(raw2, "gene_id_converted_n"),
      expr_detect = det,
      fast_reader = TRUE,
      plain_text_matrix_reader = TRUE
    ))
  }

  NULL
}


# =========================
# v11.8.0 generic Excel fallbacks
# =========================
# 1) gene x sample workbook: scans multiple sheets and possible header rows.
# 2) sample x gene workbook: transposes matrices where samples are rows and genes are columns.
# 3) multi-sheet single-sample workbook: merges sheets where each sheet has gene/value columns.
# These are structural parsers only; no GSE ID, file name, or disease label is hard-coded.

try_excel_gene_x_sample_matrix_generic <- function(f) {
  bn <- safe_utf8(basename(f))
  if (!grepl("\\.xlsx$|\\.xls$", bn, ignore.case = TRUE)) return(NULL)

  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) character())
  sheets <- safe_utf8(sheets)
  if (length(sheets) == 0) {
    log_parser_reject("excel_gene_x_sample_generic", f, "no readable sheets")
    return(NULL)
  }

  reject_reasons <- character()

  for (sh in sheets) {
    dat0 <- tryCatch({
      as.data.frame(readxl::read_excel(f, sheet = sh, col_names = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE)
    }, error = function(e) {
      reject_reasons <<- c(reject_reasons, paste0("sheet=", sh, ": read failed: ", conditionMessage(e)))
      NULL
    })
    if (is.null(dat0) || nrow(dat0) < 5 || ncol(dat0) < 3) next

    dat0 <- dat0[rowSums(!is.na(dat0) & dat0 != "") > 0, , drop = FALSE]
    dat0 <- dat0[, colSums(!is.na(dat0) & dat0 != "") > 0, drop = FALSE]
    if (nrow(dat0) < 5 || ncol(dat0) < 3) next

    for (header_idx in seq_len(min(25, nrow(dat0) - 1))) {
      hdr_raw <- safe_utf8(unlist(dat0[header_idx, ], use.names = FALSE))
      empty <- is.na(hdr_raw) | trimws(hdr_raw) == ""
      hdr_raw[empty] <- paste0("V", which(empty))
      hdr <- safe_make_names(hdr_raw, unique = TRUE)

      dat <- dat0[-seq_len(header_idx), , drop = FALSE]
      colnames(dat) <- hdr
      dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
      if (nrow(dat) < 20 || ncol(dat) < 3) next

      gene_choice <- choose_gene_symbol_column(dat)
      symbol_col <- safe_one_col(gene_choice$symbol_col, colnames(dat))
      id_col <- safe_one_col(gene_choice$id_col, colnames(dat))
      if (is.na(symbol_col) || symbol_col == "") next

      cn <- colnames(dat)
      numeric_cols <- cn[sapply(dat, is_numeric_like, min_prop = 0.30)]
      numeric_cols <- safe_char_cols(numeric_cols, cn)
      sample_cols <- setdiff(numeric_cols, c(symbol_col, id_col))
      sample_cols <- sample_cols[!is_annotation_col(sample_cols)]

      # keep small cohorts. Two sample columns is enough for preview and many GEO validation datasets.
      if (length(sample_cols) < 2) {
        reject_reasons <- c(reject_reasons, paste0("sheet=", sh, ", header=", header_idx, ": sample_cols<2 (", length(sample_cols), ")"))
        next
      }

      # Avoid coordinate/length/stat columns.
      sample_cols <- sample_cols[sapply(sample_cols, function(cc) {
        vals <- suppressWarnings(as.numeric(dat[[as.character(cc)[1]]]))
        finite <- vals[is.finite(vals)]
        if (length(finite) < max(10, nrow(dat) * 0.20)) return(FALSE)
        int_prop <- mean(abs(finite - round(finite)) < 1e-6, na.rm = TRUE)
        zero_prop <- mean(finite == 0, na.rm = TRUE)
        coord_like <- median(finite, na.rm = TRUE) > 1000 && int_prop > 0.95 && zero_prop < 0.01
        !coord_like
      })]
      if (length(sample_cols) < 2) next

      raw <- dat[, safe_char_cols(c(symbol_col, sample_cols), cn), drop = FALSE]
      colnames(raw)[1] <- "Symbol"
      raw$Symbol <- clean_symbol_value(safe_utf8(raw$Symbol))
      raw <- raw[!is.na(raw$Symbol) & raw$Symbol != "" & raw$Symbol != "NA", , drop = FALSE]
      if (nrow(raw) < 50) next

      for (cc in colnames(raw)[-1]) raw[[as.character(cc)[1]]] <- suppressWarnings(as.numeric(raw[[as.character(cc)[1]]]))
      good_cols <- colnames(raw)[-1][sapply(raw[-1], function(x) mean(!is.na(x)) >= 0.30)]
      good_cols <- safe_char_cols(good_cols, colnames(raw))
      if (length(good_cols) < 2) next
      raw <- raw[, c("Symbol", good_cols), drop = FALSE]
      raw2 <- collapse_by_symbol_auto(raw, source_hint = paste(basename(f), sh))
      if (nrow(raw2) < 50 || ncol(raw2) < 3) next
      det <- attr(raw2, "expr_detect")
      attr(raw2, "symbol_col_used") <- symbol_col
      attr(raw2, "symbol_col_reason") <- paste0("generic Excel gene x sample parser; sheet=", sh, "; header row=", header_idx)
      return(list(
        raw = raw2,
        symbol_col = symbol_col,
        symbol_col_reason = attr(raw2, "symbol_col_reason"),
        gene_id_col = id_col,
        gene_id_type = attr(raw2, "gene_id_type"),
        gene_id_conversion_reason = attr(raw2, "gene_id_conversion_reason"),
        gene_id_converted_n = attr(raw2, "gene_id_converted_n"),
        expr_detect = det,
        fast_reader = FALSE,
        excel_gene_x_sample_reader = TRUE
      ))
    }
  }

  log_parser_reject("excel_gene_x_sample_generic", f, "no gene x sample matrix accepted", paste(tail(reject_reasons, 6), collapse = " | "))
  NULL
}

try_excel_sample_x_gene_matrix_generic <- function(f) {
  bn <- safe_utf8(basename(f))
  if (!grepl("\\.xlsx$|\\.xls$", bn, ignore.case = TRUE)) return(NULL)

  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) character())
  sheets <- safe_utf8(sheets)
  if (length(sheets) == 0) return(NULL)
  reject_reasons <- character()

  for (sh in sheets) {
    dat0 <- tryCatch(as.data.frame(readxl::read_excel(f, sheet = sh, col_names = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(dat0) || nrow(dat0) < 2 || ncol(dat0) < 100) next
    dat0 <- dat0[rowSums(!is.na(dat0) & dat0 != "") > 0, , drop = FALSE]
    dat0 <- dat0[, colSums(!is.na(dat0) & dat0 != "") > 0, drop = FALSE]
    if (nrow(dat0) < 2 || ncol(dat0) < 100) next

    for (header_idx in seq_len(min(20, nrow(dat0) - 1))) {
      hdr_raw <- safe_utf8(unlist(dat0[header_idx, ], use.names = FALSE))
      empty <- is.na(hdr_raw) | trimws(hdr_raw) == ""
      hdr_raw[empty] <- paste0("V", which(empty))
      hdr <- safe_make_names(hdr_raw, unique = TRUE)
      dat <- dat0[-seq_len(header_idx), , drop = FALSE]
      colnames(dat) <- hdr
      dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
      if (nrow(dat) < 2 || ncol(dat) < 100) next

      cn <- colnames(dat)
      # sample ID/group descriptor column: non-numeric column with mostly unique row labels.
      non_num_cols <- cn[!sapply(dat, is_numeric_like, min_prop = 0.50)]
      sample_col <- NA_character_
      for (cc in non_num_cols) {
        vals <- safe_utf8(dat[[as.character(cc)[1]]])
        vals <- vals[!is.na(vals) & vals != ""]
        if (length(vals) >= 2 && length(unique(vals)) >= max(2, floor(length(vals) * 0.5))) { sample_col <- cc; break }
      }
      if (is.na(sample_col)) sample_col <- cn[1]

      numeric_cols <- cn[sapply(dat, is_numeric_like, min_prop = 0.50)]
      numeric_cols <- safe_char_cols(numeric_cols, cn)
      numeric2 <- numeric_cols[!is_annotation_col(numeric_cols)]
      # For sample x gene matrices, column names themselves should look like gene IDs/symbols.
      gene_cols <- numeric2[sapply(numeric2, function(z) looks_like_bulk_gene_vector(z))]
      if (length(gene_cols) < 100) {
        reject_reasons <- c(reject_reasons, paste0("sheet=", sh, ", header=", header_idx, ": gene_cols<100 (", length(gene_cols), ")"))
        next
      }

      sample_names <- safe_utf8(dat[[as.character(sample_col)[1]]])
      sample_names[is.na(sample_names) | trimws(sample_names) == ""] <- paste0("Sample", seq_along(sample_names))[is.na(sample_names) | trimws(sample_names) == ""]
      sample_names <- make.unique(safe_make_names(sample_names, unique = FALSE))

      mat <- as.matrix(dat[, gene_cols, drop = FALSE])
      storage.mode(mat) <- "numeric"
      if (mean(!is.na(mat)) < 0.30) next
      raw <- as.data.frame(t(mat), stringsAsFactors = FALSE)
      colnames(raw) <- sample_names
      raw <- cbind(Symbol = gene_cols, raw, stringsAsFactors = FALSE)
      raw2 <- collapse_by_symbol_auto(raw, source_hint = paste(basename(f), sh))
      if (nrow(raw2) < 50 || ncol(raw2) < 3) next
      det <- attr(raw2, "expr_detect")
      attr(raw2, "symbol_col_used") <- "column names"
      attr(raw2, "symbol_col_reason") <- paste0("generic Excel sample x gene transpose parser; sheet=", sh, "; header row=", header_idx)
      return(list(
        raw = raw2,
        symbol_col = "column names",
        symbol_col_reason = attr(raw2, "symbol_col_reason"),
        gene_id_col = "column names",
        gene_id_type = attr(raw2, "gene_id_type"),
        gene_id_conversion_reason = attr(raw2, "gene_id_conversion_reason"),
        gene_id_converted_n = attr(raw2, "gene_id_converted_n"),
        expr_detect = det,
        fast_reader = FALSE,
        excel_sample_x_gene_reader = TRUE
      ))
    }
  }

  log_parser_reject("excel_sample_x_gene_generic", f, "no sample x gene matrix accepted", paste(tail(reject_reasons, 6), collapse = " | "))
  NULL
}

try_excel_multisheet_single_sample_matrix_generic <- function(f) {
  bn <- safe_utf8(basename(f))
  if (!grepl("\\.xlsx$|\\.xls$", bn, ignore.case = TRUE)) return(NULL)

  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) character())
  sheets <- safe_utf8(sheets)
  if (length(sheets) < 2) return(NULL)

  count_list <- list()
  reject_reasons <- character()

  for (sh in sheets) {
    dat0 <- tryCatch(as.data.frame(readxl::read_excel(f, sheet = sh, col_names = FALSE, .name_repair = "minimal"), stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(dat0) || nrow(dat0) < 30 || ncol(dat0) < 2) next

    accepted <- FALSE
    for (header_idx in 0:min(20, nrow(dat0) - 1)) {
      dat <- dat0
      if (header_idx > 0) {
        hdr <- safe_utf8(unlist(dat0[header_idx, ], use.names = FALSE))
        empty <- is.na(hdr) | trimws(hdr) == ""
        hdr[empty] <- paste0("V", which(empty))
        colnames(dat) <- safe_make_names(hdr, unique = TRUE)
        dat <- dat[-seq_len(header_idx), , drop = FALSE]
      } else {
        colnames(dat) <- paste0("V", seq_len(ncol(dat)))
      }
      dat <- dat[rowSums(!is.na(dat) & dat != "") > 0, , drop = FALSE]
      dat <- dat[, colSums(!is.na(dat) & dat != "") > 0, drop = FALSE]
      if (nrow(dat) < 30 || ncol(dat) < 2) next

      best <- choose_best_gene_value_columns(dat)
      if (is.null(best) || is.na(best$gene_col) || is.na(best$value_col)) next
      gene <- safe_utf8(dat[[as.character(best$gene_col)[1]]])
      val <- suppressWarnings(as.numeric(dat[[as.character(best$value_col)[1]]]))
      if (!is_probably_bulk_expression_pair(gene, val, min_rows = 50)) {
        reject_reasons <- c(reject_reasons, paste0("sheet=", sh, ", header=", header_idx, ": pair failed structure check"))
        next
      }
      one <- data.frame(Symbol = clean_symbol_value(gene), count = val, stringsAsFactors = FALSE)
      one <- one[!is.na(one$Symbol) & one$Symbol != "" & !is.na(one$count), , drop = FALSE]
      one <- collapse_single_sample_duplicates(one)
      if (is.null(one) || nrow(one) < 50) next
      sample_name <- safe_make_names(sh, unique = FALSE)[1]
      colnames(one)[2] <- sample_name
      count_list[[sample_name]] <- one
      accepted <- TRUE
      break
    }
    if (!accepted) reject_reasons <- c(reject_reasons, paste0("sheet=", sh, ": no gene/value pair accepted"))
  }

  if (length(count_list) < 2) {
    log_parser_reject("excel_multisheet_single_sample_generic", f, "accepted sheets <2", paste(tail(reject_reasons, 6), collapse = " | "))
    return(NULL)
  }

  raw <- Reduce(function(x, y) full_join(x, y, by = "Symbol"), count_list)
  raw2 <- collapse_by_symbol_auto(raw, source_hint = paste(basename(f), names(count_list), collapse = " "))
  if (nrow(raw2) < 50 || ncol(raw2) < 3) return(NULL)
  det <- attr(raw2, "expr_detect")
  attr(raw2, "symbol_col_used") <- "per-sheet gene column"
  attr(raw2, "symbol_col_reason") <- paste0("generic Excel multi-sheet single-sample parser; sheets merged=", length(count_list))
  list(
    raw = raw2,
    symbol_col = attr(raw2, "symbol_col_used"),
    symbol_col_reason = attr(raw2, "symbol_col_reason"),
    gene_id_col = "per-sheet gene column",
    gene_id_type = attr(raw2, "gene_id_type"),
    gene_id_conversion_reason = attr(raw2, "gene_id_conversion_reason"),
    gene_id_converted_n = attr(raw2, "gene_id_converted_n"),
    expr_detect = det,
    fast_reader = FALSE,
    excel_multisheet_single_sample_reader = TRUE
  )
}


# v11.8.7 hotfix: obvious per-sample count files should not be parsed as gene x sample matrices.
# Files such as GSM7123684_H1_COUNT.txt are one sample per file; trying matrix parsers
# can create empty sample_cols and trigger "undefined columns selected".
is_obvious_single_sample_expression_file <- function(f) {
  bn <- basename(f)
  has_table_ext <- grepl("\\.txt(\\.gz)?$|\\.tsv(\\.gz)?$|\\.csv(\\.gz)?$|\\.count(s)?$|htseq\\.results$|quant\\.sf$",
                         bn, ignore.case = TRUE)
  has_single_sample_prefix <- grepl("^GSM[0-9]+", bn, ignore.case = TRUE)
  has_count_keyword <- grepl("_COUNT|counts?|htseq|featureCounts|gene_counts|expected_count|abundance|quant\\.sf|salmon|kallisto",
                             bn, ignore.case = TRUE)
  has_table_ext && has_single_sample_prefix && has_count_keyword
}

parse_expression_candidate_file <- function(f) {
  if (is_raw_archive_filename(f) || is_scRNA_like_filename(f)) return(NULL)

  # Hotfix: per-sample COUNT/TXT files must be handled later by
  # merge_single_sample_expr_files_robust(), not by matrix parsers.
  if (is_obvious_single_sample_expression_file(f)) {
    message("跳过gene×sample矩阵parser，留给单样本合并器处理: ", basename(f))
    return(NULL)
  }

  message("尝试按候选表达矩阵读取: ", basename(f))

  # Each parser is isolated. A failure in one parser should not hide the next parser.
  # The console will now report exactly which parser failed and why.
  parsers <- list(
    excel_gene_x_sample_generic = try_excel_gene_x_sample_matrix_generic,
    excel_sample_x_gene_generic = try_excel_sample_x_gene_matrix_generic,
    excel_multisheet_single_sample_generic = try_excel_multisheet_single_sample_matrix_generic,
    submitter_rpm_excel = try_submitter_rpm_excel_matrix,
    plain_gene_sample_text = try_plain_gene_sample_text_matrix,
    ultra_fast_matrix = try_ultra_fast_expression_matrix,
    fast_submitted_matrix = try_fast_submitted_expression_matrix,
    ready_matrix = try_ready_matrix_file
  )

  for (nm in names(parsers)) {
    res <- safe_parser_try(nm, parsers[[nm]], f)
    if (!is.null(res)) return(res)
  }

  message("No parser accepted file: ", basename(f))
  NULL
}

load_method_preferred_expression_source <- function(gse_id, method_preference = "auto") {
  supp_dir <- file.path(GEO_ROOT, gse_id)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  # Re-scan online GEO data sources for the requested method. This is critical for Force limma:
  # the app may have loaded raw_counts for Auto/DESeq2, but a TPM/FPKM matrix may also exist online.
  try(download_geo_supplementary_processed_first(
    gse_id,
    supp_dir,
    allow_raw_fallback = FALSE,
    method_preference = method_preference
  ), silent = TRUE)

  files <- list_expression_candidate_files(supp_dir)
  files <- files[!is_scRNA_like_filename(files)]
  files <- files[!is_raw_archive_filename(files)]
  files <- files[is_processed_matrix_filename(files) | is_ncbi_generated_matrix_filename(files) |
                   grepl("processed|expression|expr|matrix|mRNA|RNA[_-]?seq|RNAseq|Seq|Genes|GeneList|Excel|FPKM|TPM|RPM|RPKM|CPM|raw[_-]?counts?|norm[_-]?counts?|NCBI|count_matrix|counts_matrix|normalized|featureCounts|gene_counts|abundance|quant|expected_count",
                         basename(files), ignore.case = TRUE)]
  files <- unique(files)
  if (length(files) == 0) return(NULL)

  src <- ifelse(is_ncbi_generated_matrix_filename(files), "NCBI-generated", "submitter/local-processed")
  pri <- mapply(method_aware_source_priority, files, src, MoreArgs = list(method_preference = method_preference))
  files <- files[order(pri, basename(files))]
  pri <- pri[order(pri, basename(files))]

  message("Method-aware Data Source Priority候选读取顺序 (", method_preference, "):")
  message(paste0("priority=", pri, " | ", basename(files), collapse = "\n"))

  for (f in files) {
    res <- parse_expression_candidate_file(f)
    if (is.null(res)) next
    det <- res$expr_detect
    if (is.null(det)) det <- detect_expression_type(res$raw[, -1, drop = FALSE], source_hint = attr(res$raw, "source_hint"))

    if (method_preference == "limma" && identical(det$type, "raw_count")) {
      message("Force limma优先寻找normalized/TPM/FPKM；跳过raw count候选: ", basename(f))
      next
    }
    if (method_preference == "DESeq2" && !identical(det$type, "raw_count")) {
      message("Force DESeq2只接受raw count；跳过非count候选: ", basename(f))
      next
    }

    return(list(
      raw = res$raw,
      expr_detect = det,
      expr_file = basename(f),
      source_file = f,
      msg = paste0(
        "Method-aware source selected for ", method_preference, ": ", basename(f), "\n",
        "数据类型自动判断：", det$type, "\n",
        "推荐差异分析方法：", det$method, "\n"
      )
    ))
  }

  NULL
}

# v10.8 generic candidate chooser for DEG.
# Principle: never choose a matrix only because its filename looks nice.
# For every compatible candidate matrix, parse it, map its sample columns to metadata,
# apply the user's current filter/merge/group settings, and score by actual usable samples.
# This fixes the real bottom-level issue: a normalized/count matrix may exist but cover fewer
# selected metadata samples than a submitter matrix. The best source is the compatible matrix
# with the highest selected-sample coverage, not a GSE-specific file name.
build_selected_meta_for_candidate <- function(meta, raw_candidate,
                                              sample_match_mode = "auto_robust",
                                              group_col,
                                              filter_col = "None",
                                              filter_value = "All",
                                              enable_merge = FALSE,
                                              mergeA_groups = NULL,
                                              mergeB_groups = NULL,
                                              mergeA_name = "Control",
                                              mergeB_name = "Case") {
  raw_sample_names <- colnames(raw_candidate)[-1]

  # v11.8.13 CRITICAL FIX:
  # Apply the user-selected Filter column/value BEFORE sample mapping and before DEG.
  # Previous versions only used the filter in preview tables, while DEG could still rebuild
  # metadata from the full unfiltered object (e.g. GSE205672 PBMC filter preview 284/125
  # but DEG used all 299/161 samples). Filtering first makes DEG, plots and downstream
  # analysis truly use the selected subset.
  if (!is.null(filter_col) && filter_col != "None" &&
      !is.null(filter_value) && filter_value != "All" &&
      filter_col %in% colnames(meta)) {
    filter_value_clean <- clean_group_value(filter_value)
    keep_filter <- clean_group_value(meta[[filter_col]]) == filter_value_clean
    keep_filter[is.na(keep_filter)] <- FALSE
    meta <- meta[keep_filter, , drop = FALSE]
  }

  mapping <- choose_sample_mapping(
    meta = meta,
    raw_sample_names = raw_sample_names,
    requested_mode = sample_match_mode
  )
  meta_map <- mapping$meta
  sample_col <- mapping$sample_col
  if (is.null(group_col) || !(group_col %in% colnames(meta_map))) {
    return(list(meta2 = NULL, mapping = mapping, error = "group column not found in candidate-mapped metadata"))
  }

  meta2 <- meta_map %>%
    dplyr::select(sample = dplyr::all_of(sample_col), group = dplyr::all_of(group_col)) %>%
    dplyr::mutate(
      sample = if (sample_col == "Expression_Sample_Order") {
        as.character(sample)
      } else {
        standardize_sample_for_match(sample, raw_sample_names,
                                     mode = ifelse(sample_match_mode == "auto_robust", "clean", sample_match_mode))
      },
      group = clean_group_value(group)
    ) %>%
    dplyr::filter(!is.na(sample), !is.na(group), group != "", group != "NA")

  if (!is.null(filter_col) && filter_col != "None" &&
      !is.null(filter_value) && filter_value != "All" &&
      filter_col %in% colnames(meta_map)) {
    filter_df <- data.frame(
      sample = if (sample_col == "Expression_Sample_Order") {
        as.character(meta_map[[sample_col]])
      } else {
        standardize_sample_for_match(meta_map[[sample_col]], raw_sample_names,
                                     mode = ifelse(sample_match_mode == "auto_robust", "clean", sample_match_mode))
      },
      filter_value = clean_group_value(meta_map[[filter_col]]),
      stringsAsFactors = FALSE
    )
    meta2 <- dplyr::left_join(meta2, filter_df, by = "sample") %>%
      dplyr::filter(.data$filter_value == filter_value) %>%
      dplyr::select(sample, group)
  }

  meta2 <- apply_manual_group_merge(
    meta2,
    enable_merge = enable_merge,
    mergeA_groups = mergeA_groups,
    mergeB_groups = mergeB_groups,
    mergeA_name = mergeA_name,
    mergeB_name = mergeB_name
  )

  list(meta2 = meta2, mapping = mapping, error = NULL)
}

score_candidate_for_current_deg <- function(res, file, source, method_preference,
                                            meta, sample_match_mode,
                                            group_col, filter_col, filter_value,
                                            enable_merge, mergeA_groups, mergeB_groups,
                                            mergeA_name, mergeB_name,
                                            groupA_run, groupB_run) {
  if (is.null(res) || is.null(res$raw)) return(NULL)
  det <- res$expr_detect
  if (is.null(det)) det <- detect_expression_type(res$raw[, -1, drop = FALSE], source_hint = attr(res$raw, "source_hint"))

  is_raw <- identical(det$type, "raw_count")
  if (method_preference == "DESeq2" && !is_raw) return(NULL)
  if (method_preference == "limma" && is_raw) {
    # limma on raw counts is a fallback, not preferred when normalized candidates exist.
    raw_penalty <- 100000L
  } else {
    raw_penalty <- 0L
  }

  bm <- build_selected_meta_for_candidate(
    meta = meta,
    raw_candidate = res$raw,
    sample_match_mode = sample_match_mode,
    group_col = group_col,
    filter_col = filter_col,
    filter_value = filter_value,
    enable_merge = enable_merge,
    mergeA_groups = mergeA_groups,
    mergeB_groups = mergeB_groups,
    mergeA_name = mergeA_name,
    mergeB_name = mergeB_name
  )
  if (is.null(bm$meta2) || nrow(bm$meta2) == 0) return(NULL)

  df2 <- bm$meta2 %>%
    dplyr::filter(group %in% c(groupA_run, groupB_run)) %>%
    dplyr::distinct(sample, .keep_all = TRUE)
  if (nrow(df2) == 0 || length(unique(df2$group)) < 2) return(NULL)

  resolved <- resolve_samples_for_deg(df2, colnames(res$raw)[-1])
  kept <- resolved$kept
  dropped <- resolved$dropped
  kept_counts <- table(kept$group)
  both_groups_ok <- all(c(groupA_run, groupB_run) %in% names(kept_counts)) && all(kept_counts[c(groupA_run, groupB_run)] >= 2)
  if (!both_groups_ok) return(NULL)

  base_priority <- method_aware_source_priority(file, source, method_preference)
  coverage_n <- nrow(kept)
  selected_n <- nrow(df2)
  coverage_ratio <- ifelse(selected_n > 0, coverage_n / selected_n, 0)

  # Sort ascending: maximize coverage first, then prefer compatible normalized/count priority.
  sort_score <- (-coverage_n * 1000000) + raw_penalty + base_priority

  list(
    raw = res$raw,
    expr_detect = det,
    expr_file = basename(file),
    source_file = file,
    source = source,
    meta2 = bm$meta2,
    mapping = bm$mapping,
    selected_n = selected_n,
    kept_n = coverage_n,
    dropped_n = nrow(dropped),
    coverage_ratio = coverage_ratio,
    base_priority = base_priority,
    sort_score = sort_score,
    msg = paste0(
      "Coverage-aware source selected for ", method_preference, ": ", basename(file), "\n",
      "数据类型自动判断：", det$type, "\n",
      "Selected metadata samples covered: ", coverage_n, "/", selected_n,
      if (nrow(dropped) > 0) paste0("; dropped=", nrow(dropped)) else "",
      "\nBase file priority: ", base_priority, "\n"
    )
  )
}

load_best_expression_source_for_current_deg <- function(gse_id, method_preference,
                                                        meta, sample_match_mode,
                                                        group_col, filter_col, filter_value,
                                                        enable_merge, mergeA_groups, mergeB_groups,
                                                        mergeA_name, mergeB_name,
                                                        groupA_run, groupB_run) {
  supp_dir <- file.path(GEO_ROOT, gse_id)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  try(download_geo_supplementary_processed_first(
    gse_id,
    supp_dir,
    allow_raw_fallback = FALSE,
    method_preference = method_preference
  ), silent = TRUE)

  files <- list_expression_candidate_files(supp_dir)
  files <- files[!is_scRNA_like_filename(files)]
  files <- files[!is_raw_archive_filename(files)]
  files <- files[is_processed_matrix_filename(files) | is_ncbi_generated_matrix_filename(files) |
                   grepl("processed|expression|expr|matrix|mRNA|RNA[_-]?seq|RNAseq|Seq|Genes|GeneList|Excel|FPKM|TPM|RPM|RPKM|CPM|raw[_-]?counts?|norm[_-]?counts?|NCBI|count_matrix|counts_matrix|normalized|featureCounts|gene_counts|abundance|quant|expected_count",
                         basename(files), ignore.case = TRUE)]
  files <- unique(files)
  if (length(files) == 0) return(NULL)

  src <- ifelse(is_ncbi_generated_matrix_filename(files), "NCBI-generated", "submitter/local-processed")
  pri <- mapply(method_aware_source_priority, files, src, MoreArgs = list(method_preference = method_preference))
  ord <- order(pri, basename(files))
  files <- files[ord]
  src <- src[ord]
  pri <- pri[ord]

  message("Coverage-aware DEG candidate scan (", method_preference, "):")
  message(paste0("base_priority=", pri, " | ", src, " | ", basename(files), collapse = "\n"))

  scored <- list()
  for (i in seq_along(files)) {
    f <- files[[i]]
    res <- parse_expression_candidate_file(f)
    if (is.null(res)) next
    sc <- score_candidate_for_current_deg(
      res = res,
      file = f,
      source = src[[i]],
      method_preference = method_preference,
      meta = meta,
      sample_match_mode = sample_match_mode,
      group_col = group_col,
      filter_col = filter_col,
      filter_value = filter_value,
      enable_merge = enable_merge,
      mergeA_groups = mergeA_groups,
      mergeB_groups = mergeB_groups,
      mergeA_name = mergeA_name,
      mergeB_name = mergeB_name,
      groupA_run = groupA_run,
      groupB_run = groupB_run
    )
    if (!is.null(sc)) {
      message("Candidate usable: ", basename(f), " | kept=", sc$kept_n, "/", sc$selected_n,
              " | dropped=", sc$dropped_n, " | score=", sc$sort_score)
      scored[[length(scored) + 1]] <- sc
    }
  }
  if (length(scored) == 0) return(NULL)
  scored[[which.min(vapply(scored, function(z) z$sort_score, numeric(1)))]]
}



# =========================
# v11.0 Manual expression source pool
# =========================
# Build a generic candidate pool for raw counts, TPM/RPM/FPKM/CPM and submitter matrices.
# This is NOT GSE-specific. Auto keeps the Load-stage source; manual choice switches the
# global active matrix used by DEG and all downstream plots.
short_source_type_label <- function(det) {
  if (is.null(det) || is.null(det$type)) return("unknown")
  if (identical(det$type, "raw_count")) return("raw_count")
  if (identical(det$type, "normalized_expression_FPKM_TPM_or_similar")) return("normalized/RPM/TPM/FPKM")
  if (identical(det$type, "log_or_microarray_expression")) return("log/microarray")
  det$type
}

make_expr_source_label <- function(prefix, expr_file, det, raw) {
  paste0(prefix, " | ", basename(expr_file),
         " | ", short_source_type_label(det),
         " | ", ncol(raw) - 1, " samples",
         " | ", nrow(raw), " genes")
}

build_manual_expression_source_pool <- function(gse_id, loaded_raw, loaded_det, loaded_expr_file) {
  pool <- list()
  pool[["auto"]] <- list(
    id = "auto",
    label = make_expr_source_label("Auto/current loaded", loaded_expr_file, loaded_det, loaded_raw),
    raw = loaded_raw,
    expr_detect = loaded_det,
    expr_file = loaded_expr_file,
    source = "auto/current-loaded"
  )

  supp_dir <- file.path(GEO_ROOT, gse_id)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  # Download processed supplementary matrices first. This keeps the generic behavior:
  # author matrices are available in the source selector, but Auto still follows the
  # normal loaded source unless the user manually switches.
  try(download_geo_supplementary_processed_first(
    gse_id,
    supp_dir,
    allow_raw_fallback = FALSE,
    method_preference = "auto"
  ), silent = TRUE)

  files <- list_expression_candidate_files(supp_dir)
  files <- files[!is_scRNA_like_filename(files)]
  files <- files[!is_raw_archive_filename(files)]
  files <- files[is_processed_matrix_filename(files) | is_ncbi_generated_matrix_filename(files) |
                   grepl("processed|expression|expr|matrix|mRNA|RNA[_-]?seq|RNAseq|Seq|Genes|GeneList|Excel|FPKM|TPM|RPM|RPKM|CPM|raw[_-]?counts?|norm[_-]?counts?|NCBI|count_matrix|counts_matrix|normalized|featureCounts|gene_counts|abundance|quant|expected_count",
                         basename(files), ignore.case = TRUE)]
  files <- unique(files)
  if (length(files) == 0) return(pool)

  # Prefer human-readable order in the dropdown: raw first, submitter normalized next, NCBI normalized next.
  src <- ifelse(is_ncbi_generated_matrix_filename(files), "NCBI-generated", "submitter/local-processed")
  pri_auto <- mapply(method_aware_source_priority, files, src, MoreArgs = list(method_preference = "auto"))
  ord <- order(pri_auto, basename(files))
  files <- files[ord]
  src <- src[ord]

  seen <- character()
  for (i in seq_along(files)) {
    f <- files[[i]]
    res <- tryCatch(parse_expression_candidate_file(f), error = function(e) NULL)
    if (is.null(res) || is.null(res$raw)) next
    det <- res$expr_detect
    if (is.null(det)) det <- detect_expression_type(res$raw[, -1, drop = FALSE], source_hint = paste(f, attr(res$raw, "source_hint")))

    # Deduplicate by basename + sample count + type. Keep both raw and normalized when they differ.
    sig <- paste(basename(f), ncol(res$raw) - 1, nrow(res$raw), short_source_type_label(det), sep = "|")
    if (sig %in% seen) next
    seen <- c(seen, sig)

    id <- paste0("src_", length(pool))
    pool[[id]] <- list(
      id = id,
      label = make_expr_source_label(src[[i]], basename(f), det, res$raw),
      raw = res$raw,
      expr_detect = det,
      expr_file = basename(f),
      source_file = f,
      source = src[[i]]
    )
  }

  pool
}

get_selected_expression_source <- function(pool, source_id) {
  if (is.null(pool) || length(pool) == 0) return(NULL)
  if (is.null(source_id) || source_id == "" || !(source_id %in% names(pool))) source_id <- "auto"
  pool[[source_id]]
}

get_geo_ncbi_generated_urls <- function(gse_id) {
  page_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gse_id)
  html <- tryCatch(readLines(page_url, warn = FALSE), error = function(e) NULL)
  if (is.null(html)) return(data.frame(file = character(), url = character(), source = character(), stringsAsFactors = FALSE))

  href <- regmatches(html, gregexpr('href="[^"]+"', html))
  href <- unlist(href)
  href <- gsub('href="|"$', "", href)
  href <- href[href != ""]

  files <- vapply(href, extract_filename_from_geo_href, character(1))
  urls <- vapply(href, make_absolute_geo_url, character(1))

  out <- data.frame(file = files, url = urls, source = "NCBI-generated", stringsAsFactors = FALSE)
  out <- out[is_ncbi_generated_matrix_filename(out$file), , drop = FALSE]
  if (nrow(out) == 0) return(out)

  out$priority <- vapply(out$file, ncbi_generated_priority, integer(1))
  out <- out[order(out$priority, out$file), , drop = FALSE]
  out <- out[!duplicated(out$file), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# v7.5: robust processed-file detection.
# GEO processed files are not always named "expression_matrix".
# Examples include: *_mRNA_Seq_422Genes_Excel.xlsx, *_RNAseq.xlsx,
# *_genes.xlsx, *_readcount.csv, *_TPM.txt.
is_processed_matrix_filename <- function(x) {
  x <- basename(x)

  has_table_ext <- grepl(
    "\\.txt(\\.gz)?$|\\.tsv(\\.gz)?$|\\.csv(\\.gz)?$|\\.xlsx(\\.gz)?$|\\.xls(\\.gz)?$",
    x,
    ignore.case = TRUE
  )

  has_processed_keyword <- grepl(
    paste0(
      "matrix|expression|expr|processed|normalized|normalised|",
      "counts_matrix|count_matrix|gene_counts|gene.counts|",
      "raw[_-]?counts?|norm[_-]?counts?|normalized[_-]?counts?|NCBI|",
      "all.*sample|submitted|FPKM|TPM|RPKM|RPM|CPM|",
      "readcount|read_count|read.count|counts?\\.txt|counts?\\.tsv|counts?\\.csv|",
      "mRNA|RNA[_-]?seq|RNAseq|transcriptome|",
      "genes?|gene[_-]?list|Excel|xlsx"
    ),
    x,
    ignore.case = TRUE
  )

  is_bad <- grepl(
    "RAW\\.tar|_RAW\\.tar|raw\\.tar|fastq|fq\\.gz|bam$|sam$|sra$|matrix\\.mtx|barcodes|features\\.tsv|genes\\.tsv|hashing|cellranger|HTO|ADT",
    x,
    ignore.case = TRUE
  )

  has_table_ext & has_processed_keyword & !is_bad
}

is_raw_archive_filename <- function(x) {
  grepl("RAW\\.tar|_RAW\\.tar|raw\\.tar|\\.tar$|\\.tar\\.gz$|\\.tgz$|\\.zip$",
        basename(x),
        ignore.case = TRUE)
}

download_geo_supplementary_processed_first <- function(gse_id, supp_dir, allow_raw_fallback = TRUE, method_preference = "auto") {
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  remote <- get_geo_supplementary_urls(gse_id)
  ncbi_generated <- get_geo_ncbi_generated_urls(gse_id)

  if (nrow(remote) == 0 && nrow(ncbi_generated) == 0) {
    message("无法预扫描supplementary/NCBI-generated列表，回退到getGEOSuppFiles。")
    GEOquery::getGEOSuppFiles(gse_id, baseDir = dirname(supp_dir), makeDirectory = TRUE)
    return(invisible(TRUE))
  }

  # Submitter-supplied processed matrices from GEO FTP supplementary directory.
  processed <- if (nrow(remote) > 0) remote[is_processed_matrix_filename(remote$file), , drop = FALSE] else data.frame()
  raw_archives <- if (nrow(remote) > 0) remote[is_raw_archive_filename(remote$file), , drop = FALSE] else data.frame()
  other <- if (nrow(remote) > 0) remote[!remote$file %in% c(processed$file, raw_archives$file), , drop = FALSE] else data.frame()

  # v7.8 smart selector:
  # Download and try ALL organized matrices before any RAW/archive fallback.
  # This includes NCBI-generated raw_counts/TPM/FPKM matrices and submitter processed
  # matrix files. RAW.tar is only used if no processed matrix exists or all fail later.
  processed_all <- data.frame(file = character(), url = character(), source = character(), stringsAsFactors = FALSE)

  if (nrow(ncbi_generated) > 0) {
    processed_all <- rbind(
      processed_all,
      data.frame(file = ncbi_generated$file, url = ncbi_generated$url, source = "NCBI-generated", stringsAsFactors = FALSE)
    )
  }

  if (nrow(processed) > 0) {
    processed_all <- rbind(
      processed_all,
      data.frame(file = processed$file, url = processed$url, source = "submitter-processed", stringsAsFactors = FALSE)
    )
  }

  if (nrow(processed_all) > 0) {
    processed_all$priority <- mapply(method_aware_source_priority, processed_all$file, processed_all$source, MoreArgs = list(method_preference = method_preference))
    processed_all <- processed_all[order(processed_all$priority, processed_all$source, processed_all$file), , drop = FALSE]
    processed_all <- processed_all[!duplicated(processed_all$file), , drop = FALSE]

    message("发现整理好的gene × sample候选矩阵，按method-aware processed-first策略下载，暂不下载/解压RAW.tar：")
    msg <- paste0(processed_all$source, " | priority=", processed_all$priority, " | ", processed_all$file)
    message(paste(msg, collapse = "\n"))

    for (i in seq_len(nrow(processed_all))) {
      dest <- file.path(supp_dir, processed_all$file[i])
      if (!file.exists(dest)) {
        safe_download_geo_file(processed_all$url[i], dest)
      }
    }

    # Also download small metadata-like files if present; skip RAW.tar.
    if (nrow(other) > 0) {
      meta_like <- other[
        grepl("meta|metadata|clinical|phenotype|sample|annotation|design",
              other$file,
              ignore.case = TRUE) &
          grepl("\\.txt(\\.gz)?$|\\.tsv(\\.gz)?$|\\.csv(\\.gz)?$|\\.xlsx(\\.gz)?$|\\.xls(\\.gz)?$",
                other$file,
                ignore.case = TRUE),
        ,
        drop = FALSE
      ]

      if (nrow(meta_like) > 0) {
        message("同时下载metadata/clinical/sample文件：")
        message(paste(meta_like$file, collapse = "\n"))
        for (i in seq_len(nrow(meta_like))) {
          dest <- file.path(supp_dir, meta_like$file[i])
          if (!file.exists(dest)) {
            safe_download_geo_file(meta_like$url[i], dest)
          }
        }
      }
    }

    return(invisible(TRUE))
  }

  if (allow_raw_fallback && nrow(raw_archives) > 0) {
    message("未发现任何整理好的gene × sample矩阵，才下载RAW/archive作为兜底：")
    message(paste(raw_archives$file, collapse = "\n"))

    for (i in seq_len(nrow(raw_archives))) {
      dest <- file.path(supp_dir, raw_archives$file[i])
      if (!file.exists(dest)) {
        safe_download_geo_file(raw_archives$url[i], dest)
      }
    }

    return(invisible(TRUE))
  }

  if (nrow(other) > 0) {
    message("未发现processed matrix或RAW archive，下载其它supplementary文件尝试解析。")
    for (i in seq_len(nrow(other))) {

      # v11.7.2 minimal safe fix:
      # Some GEO pages include footer/help links (e.g., HHS vulnerability disclosure).
      # Do NOT reshape/filter the whole supplementary table here, because some columns
      # can be list-like in edge cases. Only coerce the current row URL/file safely.
      u <- tryCatch(as.character(other$url[[i]])[1], error = function(e) NA_character_)
      f <- tryCatch(as.character(other$file[[i]])[1], error = function(e) NA_character_)

      if (is.na(u) || !nzchar(u)) next
      if (is.na(f) || !nzchar(f)) f <- basename(u)

      # Skip obvious non-data web/footer links only. This keeps author-uploaded files,
      # NCBI-generated files, RAW archives, and oddly named real files intact.
      if (grepl("hhs\\.gov|vulnerability|disclaimer|privacy|accessibility|faq|help|login|mailto:",
                u, ignore.case = TRUE)) {
        message("跳过非GEO网页链接: ", u)
        next
      }

      # Also skip absolute http(s) links that accidentally got treated as relative
      # supplementary filenames and would become .../suppl/https://...
      if (grepl("^https?://", f, ignore.case = TRUE) &&
          grepl("hhs\\.gov|vulnerability|disclaimer|privacy|accessibility|faq|help|login", f, ignore.case = TRUE)) {
        message("跳过非GEO网页文件名: ", f)
        next
      }

      dest <- file.path(supp_dir, basename(f))
      if (!file.exists(dest)) {
        tryCatch({
          safe_download_geo_file(u, dest)
        }, error = function(e) {
          message("跳过下载失败链接: ", u, " | ", conditionMessage(e))
        })
      }
    }
    return(invisible(TRUE))
  }

  invisible(FALSE)
}




# v11.8.6: generic multi-source matrix merge helper.
# Some GEO series split one cohort across multiple supplementary sources, for example
# one processed workbook contains later samples and RAW/TXT contains earlier samples.
# When sources are expression-type compatible and have complementary sample columns,
# merge them by common gene Symbol instead of choosing only one source.
normalize_sample_key_for_merge <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.gz$", "", x, ignore.case = TRUE)
  x <- gsub("^X", "", x)
  x <- toupper(trimws(x))
  x <- gsub("[^A-Z0-9]+", "", x)
  # Use GSM when present so GSM6336976_HC1 and GSM6336976 are recognized as same sample.
  gsm <- regmatches(x, regexpr("GSM[0-9]+", x, ignore.case = TRUE))
  has_gsm <- !is.na(gsm) & nzchar(gsm)
  x[has_gsm] <- toupper(gsm[has_gsm])
  x
}

expression_type_family_for_merge <- function(obj) {
  tp <- tryCatch(obj$expr_detect$type, error = function(e) NA_character_)
  tp <- as.character(tp)
  if (grepl("raw_count", tp, ignore.case = TRUE)) return("raw_count")
  if (grepl("normalized|FPKM|TPM|RPKM|RPM|CPM", tp, ignore.case = TRUE)) return("normalized")
  if (grepl("log|microarray", tp, ignore.case = TRUE)) return("log")
  "unknown"
}

can_merge_expression_sources <- function(obj1, obj2, min_common_genes = 500) {
  if (is.null(obj1) || is.null(obj2) || is.null(obj1$raw) || is.null(obj2$raw)) return(FALSE)
  fam1 <- expression_type_family_for_merge(obj1)
  fam2 <- expression_type_family_for_merge(obj2)
  # Do not mix raw counts with normalized/log matrices. Normalized + log is also avoided.
  if (!identical(fam1, fam2)) return(FALSE)
  if (fam1 %in% c("unknown")) return(FALSE)
  g1 <- as.character(obj1$raw[[1]])
  g2 <- as.character(obj2$raw[[1]])
  length(intersect(g1, g2)) >= min_common_genes
}

merge_expression_source_objects <- function(obj1, obj2, meta_total_n = NA_integer_, note = "") {
  if (!can_merge_expression_sources(obj1, obj2)) return(NULL)

  raw1 <- as.data.frame(obj1$raw, stringsAsFactors = FALSE)
  raw2 <- as.data.frame(obj2$raw, stringsAsFactors = FALSE)
  colnames(raw1)[1] <- "Symbol"
  colnames(raw2)[1] <- "Symbol"

  # Avoid duplicated samples when the same sample appears in both sources.
  s1 <- colnames(raw1)[-1]
  s2 <- colnames(raw2)[-1]
  k1 <- normalize_sample_key_for_merge(s1)
  k2 <- normalize_sample_key_for_merge(s2)
  keep2 <- !(k2 %in% k1)
  raw2 <- raw2[, c(TRUE, keep2), drop = FALSE]

  if (ncol(raw2) < 2) return(NULL)

  common_genes <- intersect(as.character(raw1$Symbol), as.character(raw2$Symbol))
  if (length(common_genes) < 500) return(NULL)

  raw1 <- raw1[raw1$Symbol %in% common_genes, , drop = FALSE]
  raw2 <- raw2[raw2$Symbol %in% common_genes, , drop = FALSE]
  merged <- merge(raw1, raw2, by = "Symbol", all = FALSE)

  # Preserve numeric storage.
  for (cc in colnames(merged)[-1]) {
    merged[[cc]] <- suppressWarnings(as.numeric(merged[[cc]]))
  }
  det <- tryCatch(
    detect_expression_type(
      merged[, -1, drop = FALSE],
      source_hint = paste(obj1$expr_file, obj2$expr_file, collapse = " ")
    ),
    error = function(e) obj1$expr_detect
  )
  attr(merged, "expr_detect") <- det
  attr(merged, "collapse_method") <- if (!is.null(det$collapse_method)) det$collapse_method else attr(obj1$raw, "collapse_method")

  n1 <- ncol(raw1) - 1
  n2_added <- sum(keep2)
  nmerged <- ncol(merged) - 1

  obj1$raw <- merged
  obj1$expr_detect <- det
  obj1$expr_file <- paste0(obj1$expr_file, " + ", obj2$expr_file)
  # Refresh the leading sample/gene/file lines in the inherited message so the UI
  # does not keep showing the pre-merge source dimensions.
  obj1$msg <- gsub("共读取样本数：[^\n]+", paste0("共读取样本数：", nmerged), obj1$msg)
  obj1$msg <- gsub("基因数：[^\n]+", paste0("基因数：", nrow(merged)), obj1$msg)
  obj1$msg <- gsub("Expression file: [^\n]+", paste0("Expression file: ", obj1$expr_file), obj1$msg)
  obj1$msg <- paste0(
    obj1$msg,
    "\n\nMulti-source merge note: complementary expression sources were merged by common gene Symbol.\n",
    "Source 1 samples: ", n1, "; source 2 added samples: ", n2_added,
    "; merged samples: ", nmerged,
    ifelse(!is.na(meta_total_n), paste0(" / metadata samples: ", meta_total_n), ""), ".\n",
    "Merged expression file: ", obj1$expr_file,
    ifelse(nzchar(note), paste0("\n", note), "")
  )
  obj1
}

make_supplementary_result_object <- function(raw, meta, expr_file, meta_file, metadata_source,
                                             metadata_matched_n, metadata_sample_col, expr_detect,
                                             symbol_col = NA, gene_id_col = NA, symbol_col_reason = NA,
                                             reader_mode = "generic parser", extra_msg = "") {
  list(
    raw = as.data.frame(raw),
    meta = meta,
    expr_file = expr_file,
    meta_file = meta_file,
    metadata_source = metadata_source,
    metadata_matched_n = metadata_matched_n,
    metadata_sample_col = metadata_sample_col,
    expr_detect = expr_detect,
    msg = paste0(
      "Series Matrix/GPL无symbol，已自动使用 supplementary gene × sample 矩阵。\n",
      "共读取样本数：", ncol(raw) - 1, "\n",
      "基因数：", nrow(raw), "\n",
      "Expression file: ", expr_file, "\n",
      "Metadata file: ", meta_file, "\n",
      "Metadata source: ", metadata_source, "\n",
      "Metadata matched samples: ", metadata_matched_n, "/", ncol(raw) - 1, "\n",
      "Metadata sample column: ", metadata_sample_col, "\n",
      "Symbol列：", symbol_col, "\n",
      "Gene ID列：", gene_id_col, "\n",
      "Symbol选择原因：", symbol_col_reason, "\n",
      "读取模式：", reader_mode, "\n",
      "数据类型自动判断：", expr_detect$type, "\n",
      "推荐差异分析方法：", expr_detect$method, "\n",
      "重复Symbol合并方法：", attr(raw, "collapse_method"), "\n",
      "Gene ID类型：", ifelse(is.null(attr(raw, "gene_id_type")), NA, attr(raw, "gene_id_type")), "\n",
      "Gene ID转换：", ifelse(is.null(attr(raw, "gene_id_conversion_reason")), NA, attr(raw, "gene_id_conversion_reason")), "\n",
      "Gene ID转换数量：", ifelse(is.null(attr(raw, "gene_id_converted_n")), NA, attr(raw, "gene_id_converted_n")), "\n",
      "判断原因：", expr_detect$reason,
      extra_msg
    )
  )
}

# v11.7.9 helper: when a submitter processed workbook exists but none of the
# expression parsers accept it, fall back to RAW/archive files. This is still
# generic: it does not assume the RAW archive is RNA-seq counts. After unpacking,
# the normal structural validators decide whether there are usable single-sample
# bulk expression/count files. CEL/FASTQ/BAM archives will simply not be accepted
# by the bulk expression parsers.
download_geo_raw_archives_only <- function(gse_id, supp_dir) {
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  remote <- get_geo_supplementary_urls(gse_id)
  if (is.null(remote) || nrow(remote) == 0) {
    message("RAW/archive fallback: 无法预扫描supplementary列表。")
    return(invisible(FALSE))
  }

  raw_archives <- remote[is_raw_archive_filename(remote$file), , drop = FALSE]
  if (nrow(raw_archives) == 0) {
    message("RAW/archive fallback: 未发现RAW/archive文件。")
    return(invisible(FALSE))
  }

  message("processed矩阵解析失败，启动RAW/archive兜底下载：")
  message(paste(raw_archives$file, collapse = "\n"))

  for (i in seq_len(nrow(raw_archives))) {
    u <- tryCatch(as.character(raw_archives$url[[i]])[1], error = function(e) NA_character_)
    f <- tryCatch(as.character(raw_archives$file[[i]])[1], error = function(e) NA_character_)
    if (is.na(u) || !nzchar(u) || is.na(f) || !nzchar(f)) next
    if (grepl("hhs\\.gov|vulnerability|disclaimer|privacy|accessibility|faq|help|login|mailto:",
              u, ignore.case = TRUE)) next
    dest <- file.path(supp_dir, basename(f))
    if (!file.exists(dest)) {
      tryCatch({
        safe_download_geo_file(u, dest)
      }, error = function(e) {
        message("RAW/archive fallback下载失败: ", u, " | ", conditionMessage(e))
      })
    }
  }
  invisible(TRUE)
}

get_supplementary_expr <- function(gse_id, force_raw_fallback = FALSE) {
  supp_dir <- file.path("GEO_downloads", gse_id)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

  message("开始尝试 supplementary expression matrix...")
  message("正在预扫描 supplementary files: ", gse_id)

  # New strategy:
  # 1) Download processed/gene x sample matrix first.
  # 2) Skip huge RAW.tar if processed matrix exists.
  # 3) RAW/archive is only fallback when no processed matrix is available.
  download_geo_supplementary_processed_first(
    gse_id = gse_id,
    supp_dir = supp_dir,
    allow_raw_fallback = TRUE
  )

  if (isTRUE(force_raw_fallback)) {
    download_geo_raw_archives_only(gse_id, supp_dir)
  }

  local_files_now <- list.files(supp_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  local_processed <- local_files_now[is_processed_matrix_filename(local_files_now) | is_ncbi_generated_matrix_filename(local_files_now)]
  if (length(local_processed) > 1) {
    lp_source <- ifelse(is_ncbi_generated_matrix_filename(local_processed), "NCBI-generated", "submitter/local-processed")
    lp_pri <- mapply(smart_processed_priority, local_processed, lp_source)
    local_processed <- local_processed[order(lp_pri, basename(local_processed))]
  }

  archives <- list.files(supp_dir, pattern = "\\.tar$|\\.tar\\.gz$|\\.tgz$|\\.zip$",
                         recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  # v7.3 critical fix:
  # If a processed table/XLSX exists, do NOT unpack RAW.tar from an older run/cache.
  # Otherwise tiny mapping-statistics files from RAW.tar can flood the parser and hide
  # the real processed expression matrix.
  if (length(local_processed) > 0 && !isTRUE(force_raw_fallback)) {
    message("检测到本地processed表达文件，先跳过RAW/archive解压；若processed解析失败，将自动进入RAW/archive兜底：")
    message(paste(basename(local_processed), collapse = "\n"))
    archives <- character(0)
  }

  if (length(archives) > 0) {
    for (af in archives) {
      message("正在解压压缩包: ", basename(af))
      if (grepl("\\.zip$", af, ignore.case = TRUE)) {
        try(utils::unzip(af, exdir = dirname(af)), silent = TRUE)
      } else {
        try(utils::untar(af, exdir = dirname(af)), silent = TRUE)
      }
    }
  }

  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)

  # Mixed bulk + scRNA GSE support: report and skip scRNA/10X files, but keep bulk RNA files.
  report_scRNA_files(supp_dir, gse_id)

  candidate_files <- list_expression_candidate_files(supp_dir)
  candidate_files <- candidate_files[!is_scRNA_like_filename(candidate_files)]
  stop_if_only_scRNA_no_bulk_candidates(supp_dir, candidate_files, gse_id)

  if (length(candidate_files) == 0) {
    stop(paste0(
      "supplementary文件下载并解压了，但没有找到可识别的表达文件。\n\n",
      debug_supplementary_files_message(supp_dir)
    ))
  }

  # v11.8.4: if this call is a RAW/archive fallback, do not let an already-successful
  # but low-coverage processed workbook stop the fallback pass again. The RAW archive has
  # already been downloaded/unpacked above, so here we temporarily remove local processed
  # workbooks/tables from candidate priority and let the structural validators evaluate
  # extracted GSM/TXT/count files. This is generic and is triggered by coverage, not GSE ID.
  if (isTRUE(force_raw_fallback) && length(local_processed) > 0) {
    message("RAW/archive fallback mode: temporarily skipping local processed matrices so extracted RAW/TXT files can be evaluated:")
    message(paste(basename(local_processed), collapse = "
"))
    candidate_files <- setdiff(candidate_files, local_processed)
  }

  # v7.9 safe smart selector: processed gene x sample matrices first; RAW/single-sample files second.
  # This only changes file loading order. DEG/metadata modules below remain unchanged.
  matrix_priority <- candidate_files[
    is_ncbi_generated_matrix_filename(candidate_files) |
      is_processed_matrix_filename(candidate_files) |
      grepl("processed|expression|expr|matrix|mRNA|RNA[_-]?seq|RNAseq|Seq|Genes|GeneList|Excel|FPKM|TPM|RPM|RPKM|CPM|raw[_-]?counts?|norm[_-]?counts?|NCBI|count_matrix|counts_matrix|normalized|featureCounts|gene_counts|abundance|quant|expected_count",
            basename(candidate_files),
            ignore.case = TRUE)
  ]

  matrix_priority <- unique(matrix_priority)
  if (length(matrix_priority) > 1) {
    mp_source <- ifelse(is_ncbi_generated_matrix_filename(matrix_priority), "NCBI-generated", "submitter/local-processed")
    mp_pri <- mapply(method_aware_source_priority, matrix_priority, mp_source, MoreArgs = list(method_preference = "auto"))
    matrix_priority <- matrix_priority[order(mp_pri, basename(matrix_priority))]
  }
  candidate_files <- unique(c(matrix_priority, setdiff(candidate_files, matrix_priority)))

  message("Supplementary candidate files after filtering/order: ", length(candidate_files))
  if (length(candidate_files) > 0) {
    message(paste0(seq_along(head(candidate_files, 20)), ". ", basename(head(candidate_files, 20)), collapse = "
"))
  }

  for (f in candidate_files) {
    # Skip obvious RAW/archive files here; they are handled by unpacking + single-file logic.
    if (is_raw_archive_filename(f)) next

    message("尝试读取 supplementary 整理矩阵文件: ", basename(f))

    res <- tryCatch(parse_expression_candidate_file(f), error = function(e) {
      message("候选矩阵解析失败，跳过该文件继续尝试后续文件: ", basename(f), " | ", conditionMessage(e))
      NULL
    })

    if (!is.null(res)) {
      gse <- getGEO(gse_id, GSEMatrix = TRUE)
      series_meta <- get_all_series_metadata(gse)

      expr_samples <- colnames(res$raw)[-1]
      best_meta <- find_best_supplementary_metadata(
        supp_dir = supp_dir,
        expr_samples = expr_samples,
        series_meta = series_meta
      )

      if (is.null(best_meta) || is.null(best_meta$meta)) {
        meta <- series_meta
        meta_file_used <- "GEO Series Matrix pData()"
        meta_source <- "series_matrix"
        meta_matched_n <- 0
        meta_sample_col <- NA_character_
      } else {
        meta <- best_meta$meta
        meta_file_used <- best_meta$file
        meta_source <- best_meta$source
        meta_matched_n <- best_meta$matched_n
        meta_sample_col <- best_meta$sample_col
      }

      det <- res$expr_detect

      # v11.8.6 source-coverage-aware fallback + complementary multi-source merge:
      # A processed supplementary workbook can be valid but contain only part of the GEO cohort.
      # If expression samples cover much less than the Series Matrix metadata, continue to RAW.tar
      # fallback. If processed and RAW sources contain complementary samples and compatible
      # expression types, merge them by common gene Symbol instead of choosing only one.
      meta_total_n <- tryCatch(nrow(series_meta), error = function(e) NA_integer_)
      expr_n_current <- length(expr_samples)
      low_coverage <- !isTRUE(force_raw_fallback) &&
        !is.na(meta_total_n) && meta_total_n >= 4 && expr_n_current < ceiling(0.80 * meta_total_n)

      if (isTRUE(low_coverage)) {
        message(
          "Processed matrix parsed but sample coverage is low: ",
          expr_n_current, "/", meta_total_n,
          " (<80%). Continuing RAW/archive fallback; compatible complementary sources will be merged."
        )

        processed_obj_current <- make_supplementary_result_object(
          raw = as.data.frame(res$raw),
          meta = meta,
          expr_file = basename(f),
          meta_file = meta_file_used,
          metadata_source = meta_source,
          metadata_matched_n = meta_matched_n,
          metadata_sample_col = meta_sample_col,
          expr_detect = det,
          symbol_col = res$symbol_col,
          gene_id_col = ifelse(is.null(res$gene_id_col), NA, res$gene_id_col),
          symbol_col_reason = ifelse(is.null(res$symbol_col_reason), NA, res$symbol_col_reason),
          reader_mode = ifelse(isTRUE(res$ultra_fast_reader), "ultra-fast submitted expression matrix", ifelse(isTRUE(res$fast_reader), "fast submitted expression matrix", "generic parser"))
        )

        raw_fallback_obj <- tryCatch({
          download_geo_raw_archives_only(gse_id, supp_dir)
          get_supplementary_expr(gse_id, force_raw_fallback = TRUE)
        }, error = function(e) {
          message("RAW/archive fallback after low processed coverage failed: ", conditionMessage(e))
          NULL
        })

        if (!is.null(raw_fallback_obj) && !is.null(raw_fallback_obj$raw)) {
          raw_fallback_n <- ncol(raw_fallback_obj$raw) - 1

          merged_obj <- tryCatch({
            merge_expression_source_objects(
              obj1 = raw_fallback_obj,
              obj2 = processed_obj_current,
              meta_total_n = meta_total_n,
              note = paste0("Processed source: ", basename(f), "; RAW/archive source: ", raw_fallback_obj$expr_file)
            )
          }, error = function(e) {
            message("Multi-source merge failed: ", conditionMessage(e))
            NULL
          })

          if (!is.null(merged_obj) && !is.null(merged_obj$raw)) {
            merged_n <- ncol(merged_obj$raw) - 1
            if (!is.na(merged_n) && merged_n > max(expr_n_current, raw_fallback_n, na.rm = TRUE)) {
              merged_obj$msg <- paste0(
                merged_obj$msg,
                "

Source selection note: processed matrix ", basename(f),
                " had low sample coverage (", expr_n_current, "/", meta_total_n,
                "); RAW/archive fallback provided ", raw_fallback_n,
                " samples; compatible sources were merged to ", merged_n, " samples."
              )
              return(merged_obj)
            }
          }

          if (!is.na(raw_fallback_n) && raw_fallback_n > expr_n_current) {
            raw_fallback_obj$msg <- paste0(
              raw_fallback_obj$msg,
              "

Source selection note: processed matrix ", basename(f),
              " had low sample coverage (", expr_n_current, "/", meta_total_n,
              "); RAW/archive fallback provided more samples (", raw_fallback_n,
              ") and was selected. Multi-source merge was not used because sources were not compatible or not complementary."
            )
            return(raw_fallback_obj)
          } else {
            message(
              "RAW/archive fallback did not improve sample coverage (",
              raw_fallback_n, " vs ", expr_n_current,
              "). Keeping processed matrix: ", basename(f)
            )
          }
        }
      }

      return(list(
        raw = as.data.frame(res$raw),
        meta = meta,
        expr_file = basename(f),
        meta_file = meta_file_used,
        metadata_source = meta_source,
        metadata_matched_n = meta_matched_n,
        metadata_sample_col = meta_sample_col,
        expr_detect = det,
        msg = paste0(
          "Series Matrix/GPL无symbol，已自动使用 supplementary gene × sample 矩阵。\n",
          "共读取样本数：", ncol(res$raw) - 1, "\n",
          "基因数：", nrow(res$raw), "\n",
          "Expression file: ", basename(f), "\n",
          "Metadata file: ", meta_file_used, "\n",
          "Metadata source: ", meta_source, "\n",
          "Metadata matched samples: ", meta_matched_n, "/", length(expr_samples), "\n",
          "Metadata sample column: ", meta_sample_col, "\n",
          "Symbol列：", res$symbol_col, "\n",
          "Gene ID列：", ifelse(is.null(res$gene_id_col), NA, res$gene_id_col), "\n",
          "Symbol选择原因：", ifelse(is.null(res$symbol_col_reason), NA, res$symbol_col_reason), "\n",
          "读取模式：", ifelse(isTRUE(res$ultra_fast_reader), "ultra-fast submitted expression matrix", ifelse(isTRUE(res$fast_reader), "fast submitted expression matrix", "generic parser")), "\n",
          "数据类型自动判断：", det$type, "\n",
          "推荐差异分析方法：", det$method, "\n",
          "重复Symbol合并方法：", attr(res$raw, "collapse_method"), "\n",
          "Gene ID类型：", ifelse(is.null(attr(res$raw, "gene_id_type")), ifelse(is.null(res$gene_id_type), NA, res$gene_id_type), attr(res$raw, "gene_id_type")), "\n",
          "Gene ID转换：", ifelse(is.null(attr(res$raw, "gene_id_conversion_reason")), ifelse(is.null(res$gene_id_conversion_reason), NA, res$gene_id_conversion_reason), attr(res$raw, "gene_id_conversion_reason")), "\n",
          "Gene ID转换数量：", ifelse(is.null(attr(res$raw, "gene_id_converted_n")), ifelse(is.null(res$gene_id_converted_n), NA, res$gene_id_converted_n), attr(res$raw, "gene_id_converted_n")), "\n",
          "判断原因：", det$reason
        )
      ))
    }
  }

  single_files <- candidate_files[
    grepl("^GSM[0-9]+", basename(candidate_files), ignore.case = TRUE)
  ]

  htseq_count_files <- candidate_files[
    grepl("htseq\\.results|HTSeq|COUNT|count|counts|featureCounts|gene_counts|abundance|quant\\.sf|expected_count|salmon|kallisto",
          basename(candidate_files), ignore.case = TRUE)
  ]

  obvious_single_files <- candidate_files[vapply(candidate_files, is_obvious_single_sample_expression_file, logical(1))]

  single_files <- unique(c(obvious_single_files, single_files, htseq_count_files))
  single_files <- single_files[!is_scRNA_like_filename(single_files)]
  single_files <- deduplicate_expression_files_by_sample(single_files)

  # Bottom-logic rule:
  # File discovery should NOT decide which biological cohort is "correct".
  # PBMC / whole blood / tissue / cell line / time point can all be valid bulk expression files.
  # Here we only remove structural non-bulk/scRNA support files and duplicated filenames.
  # Cohort selection belongs to metadata Filter column / Compare column after loading.

  message("候选单样本表达/count文件数（去重后/结构筛选后）：", length(single_files))
  if (length(single_files) > 0) {
    message("候选单样本文件示例：", paste(head(basename(single_files), 10), collapse = ", "))
  }

  # Bulk-only: 10X/scRNA-seq parsing is intentionally disabled.
  tenx_res <- NULL
  if (!is.null(tenx_res)) {
    raw <- tenx_res$raw
    gse <- getGEO(gse_id, GSEMatrix = TRUE)
    series_meta <- get_all_series_metadata(gse)

    expr_samples <- colnames(raw)[-1]
    best_meta <- find_best_supplementary_metadata(
      supp_dir = supp_dir,
      expr_samples = expr_samples,
      series_meta = series_meta
    )

    if (is.null(best_meta) || is.null(best_meta$meta)) {
      meta <- series_meta
      meta_file_used <- "GEO Series Matrix pData()"
      meta_source <- "series_matrix"
      meta_matched_n <- 0
      meta_sample_col <- NA_character_
    } else {
      meta <- best_meta$meta
      meta_file_used <- best_meta$file
      meta_source <- best_meta$source
      meta_matched_n <- best_meta$matched_n
      meta_sample_col <- best_meta$sample_col
    }

    det <- attr(raw, "expr_detect")
    return(list(
      raw = as.data.frame(raw),
      meta = meta,
      expr_file = "10X matrix.mtx/features/barcodes",
      meta_file = meta_file_used,
      metadata_source = meta_source,
      metadata_matched_n = meta_matched_n,
      metadata_sample_col = meta_sample_col,
      expr_detect = det,
      msg = paste0(
        "Series Matrix/GPL无symbol，已自动读取 supplementary 10X matrix.mtx 表达矩阵。\n",
        "共读取样本数：", ncol(raw) - 1, "\n",
        "基因数：", nrow(raw), "\n",
        "Metadata file: ", meta_file_used, "\n",
        "Metadata source: ", meta_source, "\n",
        "Metadata matched samples: ", meta_matched_n, "/", length(expr_samples), "\n",
        "Metadata sample column: ", meta_sample_col, "\n",
        "数据类型自动判断：", det$type, "\n",
        "推荐差异分析方法：", det$method
      )
    ))
  }

  if (length(single_files) >= 2) {
    message("检测到可能的单样本表达/count 文件，开始自动合并...")
    raw <- merge_single_sample_expr_files_robust(single_files)

    if (!is.null(raw)) {
      gse <- getGEO(gse_id, GSEMatrix = TRUE)
      series_meta <- get_all_series_metadata(gse)

      expr_samples <- colnames(raw)[-1]
      best_meta <- find_best_supplementary_metadata(
        supp_dir = supp_dir,
        expr_samples = expr_samples,
        series_meta = series_meta
      )

      if (is.null(best_meta) || is.null(best_meta$meta)) {
        meta <- series_meta
        meta_file_used <- "GEO Series Matrix pData()"
        meta_source <- "series_matrix"
        meta_matched_n <- 0
        meta_sample_col <- NA_character_
      } else {
        meta <- best_meta$meta
        meta_file_used <- best_meta$file
        meta_source <- best_meta$source
        meta_matched_n <- best_meta$matched_n
        meta_sample_col <- best_meta$sample_col
      }

      det <- attr(raw, "expr_detect")
      return(list(
        raw = as.data.frame(raw),
        meta = meta,
        expr_file = paste0(ncol(raw) - 1, " single-sample expression/count files"),
        meta_file = meta_file_used,
        metadata_source = meta_source,
        metadata_matched_n = meta_matched_n,
        metadata_sample_col = meta_sample_col,
        expr_detect = det,
        msg = paste0(
          "Series Matrix/GPL无symbol，已自动合并 supplementary 单样本表达/count 文件。\n",
          "共读取样本数：", ncol(raw) - 1, "\n",
          "基因数：", nrow(raw), "\n",
          "Metadata file: ", meta_file_used, "\n",
          "Metadata source: ", meta_source, "\n",
          "Metadata matched samples: ", meta_matched_n, "/", length(expr_samples), "\n",
          "Metadata sample column: ", meta_sample_col, "\n",
          "数据类型自动判断：", det$type, "\n",
          "推荐差异分析方法：", det$method
        )
      ))
    }
  }

  final_candidates <- list_expression_candidate_files(supp_dir)

  if (!isTRUE(force_raw_fallback)) {
    message("processed/organized supplementary文件未被接受，准备尝试RAW/archive兜底。")
    download_geo_raw_archives_only(gse_id, supp_dir)
    return(get_supplementary_expr(gse_id, force_raw_fallback = TRUE))
  }

  stop(paste0(
    "找到supplementary文件，但未识别到合适的 gene × sample 表达矩阵。
",
    "已经尝试processed矩阵和RAW/archive兜底；如果RAW是CEL/FASTQ/BAM或只包含mapping statistics，则当前bulk快速流程不会直接解析。

",
    diagnose_failed_supplementary_files(final_candidates),
    "

解压后文件概览：
",
    debug_supplementary_files_message(supp_dir)
  ))
}

load_geo_expr_auto <- function(gse_id, use_cache = TRUE, rebuild_cache = FALSE) {
  gse_id <- toupper(trimws(gse_id))
  cf <- cache_file_for_gse(gse_id)

  if (isTRUE(use_cache) && !isTRUE(rebuild_cache) && file.exists(cf)) {
    obj <- readRDS(cf)
    obj <- refresh_expression_detect_current(obj)
    obj$loaded_from_cache <- TRUE
    obj$msg <- paste0(
      obj$msg,
      "\n\n本次从本地解析缓存读取，未重新扫描 supplementary 文件。"
    )
    return(obj)
  }

  obj <- tryCatch({
    safe_stage("Stage 1: GEO Series Matrix + GPL annotation", get_series_matrix_expr(gse_id))
  }, error = function(e1) {
    message("Series Matrix读取失败，进入supplementary兜底。原因：", conditionMessage(e1))
    safe_stage("Stage 2: Supplementary expression matrix", get_supplementary_expr(gse_id))
  })

  obj$loaded_from_cache <- FALSE
  saveRDS(obj, cf)
  obj
}



# =========================
# Multiple Series Matrix metadata helper
# =========================
get_all_series_metadata <- function(gse_obj) {
  if (is.null(gse_obj)) return(data.frame())

  # GEOquery may return a list of ExpressionSet objects if a GSE has multiple platforms.
  if (inherits(gse_obj, "ExpressionSet")) {
    meta <- as.data.frame(Biobase::pData(gse_obj))
    meta$SeriesMatrix_Index <- 1
    meta$SeriesMatrix_Platform <- tryCatch(Biobase::annotation(gse_obj), error = function(e) NA_character_)
    return(meta)
  }

  if (is.list(gse_obj)) {
    metas <- list()

    for (i in seq_along(gse_obj)) {
      eset <- gse_obj[[i]]
      if (!inherits(eset, "ExpressionSet")) next

      md <- as.data.frame(Biobase::pData(eset))
      md$SeriesMatrix_Index <- i
      md$SeriesMatrix_Platform <- tryCatch(Biobase::annotation(eset), error = function(e) NA_character_)
      metas[[i]] <- md
    }

    if (length(metas) == 0) return(data.frame())

    return(dplyr::bind_rows(metas))
  }

  data.frame()
}

choose_best_series_eset_for_expression <- function(gse_obj, expr_samples = NULL) {
  if (inherits(gse_obj, "ExpressionSet")) return(gse_obj)
  if (!is.list(gse_obj) || length(gse_obj) == 0) stop("GEOquery returned no ExpressionSet.")

  if (is.null(expr_samples)) return(gse_obj[[1]])

  scores <- sapply(seq_along(gse_obj), function(i) {
    eset <- gse_obj[[i]]
    if (!inherits(eset, "ExpressionSet")) return(-1)
    md <- as.data.frame(Biobase::pData(eset))
    s <- metadata_match_score(md, expr_samples)
    s$matched_n
  })

  gse_obj[[which.max(scores)]]
}


# =========================
# Supplementary metadata helper
# =========================
normalize_id_for_match <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^.*?:", "", x)
  x <- gsub("^['\"]|['\"]$", "", x)
  x <- make.names(x, unique = FALSE)
  tolower(x)
}

extract_first_gsm <- function(x) {
  x <- as.character(x)
  m <- regexpr("GSM[0-9]+", x, ignore.case = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- !is.na(m) & m > 0
  out[ok] <- toupper(regmatches(x, m)[ok])
  out
}


# =========================
# Sample alias helper
# =========================
make_expr_sample_alias_table <- function(expr_samples) {
  expr_samples <- as.character(expr_samples)

  no_gz <- sub("\\.gz$", "", expr_samples, ignore.case = TRUE)
  no_ext <- no_gz
  no_ext <- sub("\\.htseq\\.results$", "", no_ext, ignore.case = TRUE)
  no_ext <- sub("_COUNT$", "", no_ext, ignore.case = TRUE)
  no_ext <- sub("_counts$", "", no_ext, ignore.case = TRUE)
  no_ext <- sub("_count$", "", no_ext, ignore.case = TRUE)
  no_ext <- sub("\\.txt$|\\.tsv$|\\.csv$|\\.count$|\\.counts$|\\.sf$|\\.out$", "", no_ext, ignore.case = TRUE)

  # v10.3: extract GSM accession even when it is not at the very beginning
  # of the expression column name, e.g. X.GSM1548283, GSM1548283_sample,
  # or sample-GSM1548283. This prevents DEG from dropping samples that Load
  # already matched by cleaned/alias mapping.
  gsm_any <- extract_first_gsm(no_ext)
  gsm <- ifelse(!is.na(gsm_any) & gsm_any != "", gsm_any, no_ext)

  short <- sub("^.*?GSM[0-9]+[_\\.-]?", "", no_ext, ignore.case = TRUE)

  # v11.7.7: generic processed-matrix bridge.
  # Some submitter matrices use internal names such as SZ329mes while GEO metadata
  # only contains SZ329 somewhere in title/source/characteristics.  Create a
  # conservative base alias by removing common assay suffixes and by stripping a
  # trailing letter block after a letter+number sample ID. This is NOT GSE-specific.
  base <- no_ext
  base <- sub("([_\\.-]?(mrna|messengerrna|rna|rnaseq|seq|counts?|fpkm|tpm|rpm|rpkm|mes|mirna|mir))$", "", base, ignore.case = TRUE)
  base2 <- sub("^([A-Za-z]+[0-9]+)[A-Za-z]+$", "\\1", base, ignore.case = TRUE)
  base <- ifelse(nchar(base2) >= 4, base2, base)

  data.frame(
    expr_sample = expr_samples,
    alias_full = no_ext,
    alias_gsm = gsm,
    alias_short = short,
    alias_base = base,
    alias_full_norm = normalize_id_for_match(no_ext),
    alias_gsm_norm = normalize_id_for_match(gsm),
    alias_short_norm = normalize_id_for_match(short),
    alias_base_norm = normalize_id_for_match(base),
    stringsAsFactors = FALSE
  )
}

match_meta_values_to_expr_samples <- function(meta_values, expr_samples) {
  alias <- make_expr_sample_alias_table(expr_samples)

  v <- as.character(meta_values)
  v_norm <- normalize_id_for_match(v)
  v_gsm_any <- extract_first_gsm(v)
  v_gsm <- ifelse(!is.na(v_gsm_any) & v_gsm_any != "", v_gsm_any, NA_character_)
  v_gsm_norm <- normalize_id_for_match(v_gsm)

  out <- rep(NA_character_, length(v))

  # First, match full normalized names.
  m1 <- match(v_norm, alias$alias_full_norm)
  ok1 <- !is.na(m1)
  out[ok1] <- alias$expr_sample[m1[ok1]]

  # Then match extracted GSM accessions from metadata values to extracted GSM
  # accessions from expression columns. This is the critical DEG-stage fix.
  m2 <- match(v_gsm_norm, alias$alias_gsm_norm)
  ok2 <- !is.na(m2) & is.na(out)
  out[ok2] <- alias$expr_sample[m2[ok2]]

  # Also allow original metadata value to equal the expression GSM alias.
  m2b <- match(v_norm, alias$alias_gsm_norm)
  ok2b <- !is.na(m2b) & is.na(out)
  out[ok2b] <- alias$expr_sample[m2b[ok2b]]

  m3 <- match(v_norm, alias$alias_short_norm)
  ok3 <- !is.na(m3) & is.na(out)
  out[ok3] <- alias$expr_sample[m3[ok3]]

  # v11.7.7: allow exact match to conservative base alias, e.g. metadata SZ329
  # matching expression column SZ329mes.
  if ("alias_base_norm" %in% colnames(alias)) {
    m4 <- match(v_norm, alias$alias_base_norm)
    ok4 <- !is.na(m4) & is.na(out)
    out[ok4] <- alias$expr_sample[m4[ok4]]
  }

  out
}



metadata_column_quality <- function(cn, matched_n = 0, n_expr = 1) {
  cn0 <- as.character(cn)
  cn2 <- tolower(cn0)
  coverage <- matched_n / max(1, n_expr)

  bonus <- 0

  # Strong, explicit sample identifier names. These should beat accidental
  # matches from Excel description/header rows when coverage is similar.
  if (grepl("^geo_accession$|^gsm$|^accession$|^sample[._ ]?id$|^sampleid$|^sample[._ ]?name$|^sample$|^run$|^sra[._ ]?run$", cn2)) {
    bonus <- bonus + 250
  } else if (grepl("geo|gsm|accession|sample|specimen|subject|patient|run", cn2)) {
    bonus <- bonus + 80
  }

  # Clean, short names are usually intentional headers.
  if (nchar(cn2) <= 30) bonus <- bonus + 20
  if (nchar(cn2) <= 15) bonus <- bonus + 10

  # Penalize obvious Excel header failures / description rows.
  if (grepl("^x[0-9]+$|^\\.\\.\\.[0-9]+$", cn2)) bonus <- bonus - 120
  if (grepl("^gsm[0-9]+", cn2)) bonus <- bonus - 120  # a sample ID used as a column header usually means skip/header is wrong
  if (nchar(cn2) > 80) bonus <- bonus - 180
  if (grepl("merged\\.metadata|common\\.columns|added\\.columns|description|readme|note|instructions|disease\\.column", cn2)) bonus <- bonus - 220

  # Very low coverage should never win simply because the name looks nice.
  if (coverage < 0.10) bonus <- bonus - 200

  bonus
}

metadata_match_score <- function(meta, expr_samples) {
  if (is.null(meta) || nrow(meta) == 0 || ncol(meta) == 0) {
    return(list(score = 0, column = NA_character_, matched_n = 0,
                rank_score = -Inf, quality_bonus = 0, high_confidence = FALSE))
  }

  cn_all <- colnames(meta)

  # For very wide metadata sheets, test high-priority columns first. If a
  # clear sample ID column gives near-complete coverage, avoid scanning hundreds
  # of extra columns. This is important for large GEO meta files.
  priority <- cn_all[grepl("^geo_accession$|^gsm$|^accession$|^sample[._ ]?id$|^sampleid$|^sample[._ ]?name$|^sample$|^run$|^sra[._ ]?run$",
                          tolower(cn_all))]
  other <- setdiff(cn_all, priority)
  scan_cols <- c(priority, other)

  best <- list(score = 0, column = NA_character_, matched_n = 0,
               rank_score = -Inf, quality_bonus = 0, high_confidence = FALSE)

  for (cn in scan_cols) {
    matched <- match_meta_values_to_expr_samples(meta[[cn]], expr_samples)
    matched_n <- sum(!is.na(matched))
    coverage <- matched_n / max(1, length(expr_samples))
    q_bonus <- metadata_column_quality(cn, matched_n, length(expr_samples))

    # Rank uses coverage first, but a strong semantic column can beat an
    # accidental Excel-description column with 1-2 more matches.
    rank_score <- matched_n + q_bonus

    if (rank_score > best$rank_score) {
      best <- list(
        score = coverage,
        column = cn,
        matched_n = as.integer(matched_n),
        rank_score = rank_score,
        quality_bonus = q_bonus,
        high_confidence = (coverage >= 0.95 && q_bonus >= 50)
      )
    }

    # Early stop: trusted sample-ID column with near-perfect match.
    if (coverage >= 0.98 && q_bonus >= 80) break
  }

  if (best$matched_n <= 0) {
    return(list(score = 0, column = NA_character_, matched_n = 0,
                rank_score = -Inf, quality_bonus = 0, high_confidence = FALSE))
  }

  best
}

# Read metadata in multiple safe ways.
# Excel files from GEO often have a description row before the real header.
# We try several skip values and later choose the version with the best sample match
# and the cleanest sample column name.
read_metadata_variants <- function(f, max_excel_skip = 5) {
  out <- list()

  if (grepl("\\.xlsx$|\\.xls$", f, ignore.case = TRUE)) {
    for (sk in 0:max_excel_skip) {
      md <- tryCatch({
        as.data.frame(readxl::read_excel(f, skip = sk, .name_repair = "unique_quiet"))
      }, error = function(e) NULL)

      if (!is.null(md) && nrow(md) > 0 && ncol(md) > 0) {
        colnames(md) <- make.names(colnames(md), unique = TRUE)
        out[[paste0("excel_skip_", sk)]] <- md
      }
    }
  } else {
    md <- tryCatch({
      data.table::fread(f, data.table = FALSE, fill = TRUE, check.names = FALSE)
    }, error = function(e) NULL)

    if (!is.null(md) && nrow(md) > 0 && ncol(md) > 0) {
      colnames(md) <- make.names(colnames(md), unique = TRUE)
      out[["table_default"]] <- md
    }
  }

  out
}

read_metadata_any <- function(f) {
  vars <- read_metadata_variants(f)
  if (length(vars) == 0) {
    message("读取metadata失败: ", basename(f))
    return(NULL)
  }
  vars[[1]]
}


sample_column_name_bonus <- function(cn) {
  metadata_column_quality(cn, matched_n = 1, n_expr = 1)
}


find_best_supplementary_metadata <- function(supp_dir, expr_samples, series_meta = NULL) {
  files <- list.files(supp_dir, recursive = TRUE, full.names = TRUE)

  meta_files <- files[
    grepl("\\.xlsx$|\\.xls$|\\.csv$|\\.csv\\.gz$|\\.tsv$|\\.tsv\\.gz$|\\.txt$|\\.txt\\.gz$",
          files, ignore.case = TRUE) &
      grepl("metadata|meta|sample|clinical|phenotype|reanalysis|annotation",
            basename(files), ignore.case = TRUE)
  ]

  # Avoid using expression matrix as metadata.
  meta_files <- meta_files[
    !grepl("expression|expr|matrix|FPKM|TPM|count|counts|AllSampleExpression",
           basename(meta_files), ignore.case = TRUE)
  ]

  best <- NULL

  if (!is.null(series_meta)) {
    s <- metadata_match_score(series_meta, expr_samples)
    best <- list(
      meta = as.data.frame(series_meta),
      file = "GEO Series Matrix pData()",
      sample_col = s$column,
      matched_n = s$matched_n,
      score = s$score,
      rank_score = s$rank_score,
      quality_bonus = s$quality_bonus,
      high_confidence = s$high_confidence,
      source = "series_matrix"
    )
  }

  if (length(meta_files) > 0) {
    for (mf in meta_files) {
      message("尝试读取 supplementary metadata: ", basename(mf))

      variants <- read_metadata_variants(mf)
      if (length(variants) == 0) {
        message("metadata读取失败: ", basename(mf))
        next
      }

      for (variant_name in names(variants)) {
        md <- variants[[variant_name]]
        if (is.null(md) || nrow(md) == 0 || ncol(md) == 0) next

        s <- metadata_match_score(md, expr_samples)

        rank_score <- s$rank_score

        message("metadata匹配: ", basename(mf),
                " | variant=", variant_name,
                " | ", s$matched_n, "/", length(expr_samples),
                " | column=", s$column,
                " | semantic_bonus=", s$quality_bonus)

        best_rank_score <- if (is.null(best) || is.null(best$rank_score)) -Inf else best$rank_score

        if (is.null(best) || rank_score > best_rank_score) {
          best <- list(
            meta = as.data.frame(md),
            file = paste0(basename(mf), " [", variant_name, "]"),
            sample_col = s$column,
            matched_n = s$matched_n,
            score = s$score,
            rank_score = rank_score,
            quality_bonus = s$quality_bonus,
            high_confidence = s$high_confidence,
            source = "supplementary_metadata"
          )
        }

        # Early stop for large metadata files: once a semantic sample-ID column
        # matches almost all expression samples, later Excel skip variants are
        # usually header-shift artifacts and only waste time.
        if (isTRUE(s$high_confidence)) {
          message("metadata高可信匹配已找到，提前停止该文件的其他skip尝试: ",
                  basename(mf), " | variant=", variant_name,
                  " | column=", s$column,
                  " | ", s$matched_n, "/", length(expr_samples))
          break
        }
      }
    }
  }

  best
}



# =========================
# Universal group backup helper
# =========================
make_group_backup_summary <- function(meta, max_unique = 80) {
  if (is.null(meta) || nrow(meta) == 0) return(NULL)

  out <- list()

  for (cn in colnames(meta)) {
    x <- clean_group_value(meta[[cn]])
    x <- x[!is.na(x) & x != "" & x != "NA"]

    if (length(x) == 0) next

    tab <- sort(table(x), decreasing = TRUE)
    n_unique <- length(tab)

    important_name <- grepl(
      "group|condition|disease|sepsis|control|shock|sofa|mort|death|survival|outcome|phenotype|subtype|cluster|class|dataset|study|source|cohort|batch|sex|gender|race|age|ards|infection|organ|severity",
      cn,
      ignore.case = TRUE
    )

    if (n_unique >= 2 && (n_unique <= max_unique || important_name)) {
      out[[cn]] <- data.frame(
        column = cn,
        n_unique = n_unique,
        level = names(tab),
        n = as.integer(tab),
        percent = round(as.integer(tab) / sum(tab) * 100, 2),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) return(NULL)
  dplyr::bind_rows(out)
}

make_group_column_overview <- function(meta, max_unique = 120) {
  if (is.null(meta) || nrow(meta) == 0) return(NULL)

  res <- data.frame(
    column = colnames(meta),
    n_unique = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      length(unique(x[!is.na(x) & x != "" & x != "NA"]))
    }),
    missing_n = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      sum(is.na(x) | x == "" | x == "NA")
    }),
    top_values = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      x <- x[!is.na(x) & x != "" & x != "NA"]
      if (length(x) == 0) return("")
      tab <- sort(table(x), decreasing = TRUE)
      paste(paste0(names(head(tab, 10)), " (", as.integer(head(tab, 10)), ")"), collapse = " | ")
    }),
    stringsAsFactors = FALSE
  )

  res %>%
    filter(n_unique >= 2, n_unique <= max_unique) %>%
    arrange(n_unique, column)
}

choose_default_group_column <- function(meta, candidates, overview) {
  # Prefer Auto_Group if it has at least 2 valid groups.
  if ("Auto_Group" %in% colnames(meta)) {
    x <- clean_group_value(meta$Auto_Group)
    x <- x[!is.na(x) & x != "" & x != "NA"]
    if (length(unique(x)) >= 2) return("Auto_Group")
  }

  # Otherwise use highest scoring candidate.
  if (!is.null(candidates) && nrow(candidates) > 0) {
    return(candidates$column[1])
  }

  # Otherwise use first overview grouping column.
  if (!is.null(overview) && nrow(overview) > 0) {
    return(overview$column[1])
  }

  colnames(meta)[1]
}



# =========================
# Universal group discovery / metadata explorer helper
# =========================
score_group_column_v2 <- function(meta, cn) {
  x <- clean_group_value(meta[[cn]])
  x <- x[!is.na(x) & x != "" & x != "NA"]

  n <- length(x)
  n_unique <- length(unique(x))
  if (n == 0 || n_unique < 2) return(-999)

  tab <- sort(table(x), decreasing = TRUE)
  top_prop <- max(tab) / sum(tab)

  cn_low <- tolower(cn)

  bad_name <- grepl(
    "contact|email|phone|address|submission|last_update|relation|supplementary|data_processing|protocol|instrument|library|taxid|description|title|abstract|summary",
    cn_low
  )
  if (bad_name) return(-999)

  score <- 0

  # Best broad clinical grouping columns.
  if (grepl("disease|diagnosis|condition|group|phenotype|subtype|class|category", cn_low)) score <- score + 40

  # Useful binary/clinical columns.
  if (grepl("shock|outcome|mort|death|survival|ards|aki|infection|severity|status", cn_low)) score <- score + 26

  # Useful covariates.
  if (grepl("dataset|study|cohort|source|batch|platform|sex|gender|race|age", cn_low)) score <- score + 12

  # Values with biomedical group names.
  if (any(grepl(
    "sepsis|septic|control|healthy|normal|shock|ards|aki|sirs|covid|trauma|anaphylaxis|cardiogenic|critical|infection|mortality|survivor|non.?survivor|dead|alive",
    x,
    ignore.case = TRUE
  ))) score <- score + 28

  # Number of groups:
  # For an atlas, 3-30 level disease/phenotype columns are often more useful than a binary Auto_Group.
  if (n_unique >= 3 && n_unique <= 30) score <- score + 25
  else if (n_unique == 2) score <- score + 12
  else if (n_unique > 30 && n_unique <= 80) score <- score + 5
  else if (n_unique > 80) score <- score - 20

  # Avoid columns that are almost unique IDs.
  if (n_unique > min(100, n * 0.5)) score <- score - 50

  # Avoid extremely imbalanced columns unless name is clearly important.
  if (top_prop > 0.95 && !grepl("disease|condition|phenotype|subtype|class|group", cn_low)) {
    score <- score - 15
  }

  # Auto_Group is useful but should NOT hide richer original metadata.
  if (cn == "Auto_Group") score <- score - 8

  score
}

make_group_discovery_table <- function(meta, max_unique = 120) {
  if (is.null(meta) || nrow(meta) == 0) return(NULL)

  res <- data.frame(
    column = colnames(meta),
    score_v2 = sapply(colnames(meta), function(cn) score_group_column_v2(meta, cn)),
    n_unique = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      length(unique(x[!is.na(x) & x != "" & x != "NA"]))
    }),
    missing_n = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      sum(is.na(x) | x == "" | x == "NA")
    }),
    top_values = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      x <- x[!is.na(x) & x != "" & x != "NA"]
      if (length(x) == 0) return("")
      tab <- sort(table(x), decreasing = TRUE)
      paste(paste0(names(head(tab, 12)), " (", as.integer(head(tab, 12)), ")"), collapse = " | ")
    }),
    stringsAsFactors = FALSE
  )

  res %>%
    filter(score_v2 > -999, n_unique >= 2, n_unique <= max_unique) %>%
    arrange(desc(score_v2), n_unique, column)
}

make_group_backup_summary <- function(meta, max_unique = 120) {
  if (is.null(meta) || nrow(meta) == 0) return(NULL)

  out <- list()

  for (cn in colnames(meta)) {
    x <- clean_group_value(meta[[cn]])
    x <- x[!is.na(x) & x != "" & x != "NA"]

    if (length(x) == 0) next

    tab <- sort(table(x), decreasing = TRUE)
    n_unique <- length(tab)

    important_name <- grepl(
      "group|condition|disease|diagnosis|sepsis|control|shock|sofa|mort|death|survival|outcome|phenotype|subtype|cluster|class|category|dataset|study|source|cohort|batch|sex|gender|race|age|ards|aki|infection|organ|severity|status",
      cn,
      ignore.case = TRUE
    )

    if (n_unique >= 2 && (n_unique <= max_unique || important_name)) {
      out[[cn]] <- data.frame(
        column = cn,
        n_unique = n_unique,
        level = names(tab),
        n = as.integer(tab),
        percent = round(as.integer(tab) / sum(tab) * 100, 2),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0) return(NULL)
  dplyr::bind_rows(out)
}

make_group_column_overview <- function(meta, max_unique = 120) {
  make_group_discovery_table(meta, max_unique = max_unique)
}

choose_default_group_column <- function(meta, candidates, overview) {
  discovery <- make_group_discovery_table(meta)

  if (!is.null(discovery) && nrow(discovery) > 0) {
    # Prefer rich disease/phenotype/condition columns if available.
    rich <- discovery %>%
      filter(
        column != "Auto_Group",
        grepl("disease|diagnosis|condition|phenotype|subtype|class|category|group",
              column,
              ignore.case = TRUE),
        n_unique >= 3
      )

    if (nrow(rich) > 0) return(rich$column[1])

    # Otherwise use best non-Auto_Group column.
    non_auto <- discovery %>% filter(column != "Auto_Group")
    if (nrow(non_auto) > 0) return(non_auto$column[1])

    # Finally fallback to Auto_Group.
    return(discovery$column[1])
  }

  if (!is.null(candidates) && nrow(candidates) > 0) return(candidates$column[1])
  if (!is.null(overview) && nrow(overview) > 0) return(overview$column[1])
  colnames(meta)[1]
}

pick_default_groups <- function(vals) {
  vals <- vals[!is.na(vals) & vals != "" & vals != "NA"]

  if (length(vals) < 2) return(list(A = vals[1], B = vals[1]))

  control_hits <- vals[grepl("healthy|control|normal|ctrl|con$", vals, ignore.case = TRUE)]
  sepsis_hits <- vals[grepl("^sepsis$|septic|sepsis - shock|sepsis-shock|sepsis - ards|sepsis - aki", vals, ignore.case = TRUE)]

  if (length(control_hits) > 0 && length(sepsis_hits) > 0) {
    # Prefer plain "Control" and plain "Sepsis" when present.
    a <- if ("Control" %in% control_hits) "Control" else control_hits[1]
    b <- if ("Sepsis" %in% sepsis_hits) "Sepsis" else sepsis_hits[1]
    return(list(A = a, B = b))
  }

  list(A = vals[1], B = vals[2])
}

# =========================
# 5. Metadata helper functions
# =========================
score_group_column <- function(x, colname) {
  x <- clean_group_value(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  n <- length(x)
  n_unique <- length(unique(x))
  if (n == 0 || n_unique < 2) return(-999)
  bad_cols <- c(
    "submission", "last_update", "contact", "email", "phone", "address",
    "platform", "series", "supplementary", "relation", "data_processing",
    "extract_protocol", "label_protocol", "hyb_protocol", "scan_protocol",
    "instrument", "library", "taxid"
  )
  if (any(grepl(paste(bad_cols, collapse = "|"), colname, ignore.case = TRUE))) return(-999)
  score <- 0
  if (tolower(colname) == "title") score <- score - 10
  if (grepl("disease|diagnosis|group|condition|phenotype|class|outcome|mort|survival|shock|sepsis|control|healthy|case|status|source_name|characteristics",
            colname, ignore.case = TRUE)) score <- score + 12
  if (any(grepl("sepsis|septic|control|healthy|normal|case|patient|shock|survivor|non.?survivor|death|dead|alive|yes|no|con|cap|no.?cap|HC|IC|PIC|PHC",
                x, ignore.case = TRUE))) score <- score + 10
  if (n_unique >= 2 && n_unique <= 8) score <- score + 8
  else if (n_unique <= 15) score <- score + 3
  else if (n_unique > min(30, n * 0.6)) score <- score - 8
  score
}

get_group_candidates <- function(meta) {
  scores <- sapply(colnames(meta), function(cn) score_group_column(meta[[cn]], cn))
  res <- data.frame(
    column = names(scores),
    score = as.numeric(scores),
    n_unique = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      length(unique(x[!is.na(x) & x != "" & x != "NA"]))
    }),
    example_values = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      x <- unique(x[!is.na(x) & x != "" & x != "NA"])
      paste(head(x, 10), collapse = " | ")
    }),
    stringsAsFactors = FALSE
  )
  res %>% filter(score > -999, n_unique >= 2, n_unique <= 100) %>% arrange(desc(score), n_unique, column)
}

get_filter_candidates <- function(meta) {
  res <- data.frame(
    column = colnames(meta),
    n_unique = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      length(unique(x[!is.na(x) & x != "" & x != "NA"]))
    }),
    example_values = sapply(colnames(meta), function(cn) {
      x <- clean_group_value(meta[[cn]])
      x <- unique(x[!is.na(x) & x != "" & x != "NA"])
      paste(head(x, 10), collapse = " | ")
    }),
    stringsAsFactors = FALSE
  )
  bad_cols <- c(
    "submission", "last_update", "contact", "email", "phone", "address",
    "platform", "series", "supplementary", "relation", "data_processing",
    "protocol", "instrument", "library", "taxid"
  )
  res %>%
    filter(n_unique >= 2, n_unique <= 100,
           !grepl(paste(bad_cols, collapse = "|"), column, ignore.case = TRUE)) %>%
    arrange(n_unique, column)
}

clean_sample_value <- function(x) {
  x <- as.character(x)
  x <- gsub("^.*?:", "", x)
  x <- trimws(x)
  x <- gsub("^['\"]|['\"]$", "", x)
  x
}

sample_match_count <- function(meta_values, raw_sample_names, mode = "strict") {
  v <- as.character(meta_values)
  raw <- as.character(raw_sample_names)

  if (mode == "strict") {
    return(sum(v %in% raw, na.rm = TRUE))
  }

  matched <- match_meta_values_to_expr_samples(v, raw)
  sum(!is.na(matched))
}

standardize_sample_for_match <- function(x, raw_sample_names, mode = "clean") {
  x0 <- as.character(x)
  raw <- as.character(raw_sample_names)

  if (mode == "strict") {
    return(x0)
  }

  matched <- match_meta_values_to_expr_samples(x0, raw)

  out <- x0
  ok <- !is.na(matched)
  out[ok] <- matched[ok]

  out
}

# v10.3: DEG-safe sample resolver.
# The UI/Load stage may display metadata sample IDs (GSM...), while the actual
# expression matrix may use slightly different column names. DEG must use the
# true expression column names, otherwise selected samples can be silently dropped.
resolve_samples_for_deg <- function(df2, expr_sample_names) {
  expr_sample_names <- as.character(expr_sample_names)
  df2 <- as.data.frame(df2, stringsAsFactors = FALSE)
  if (!"sample" %in% colnames(df2)) stop("df2 must contain a sample column")

  df2$metadata_sample <- as.character(df2$sample)
  mapped <- match_meta_values_to_expr_samples(df2$metadata_sample, expr_sample_names)

  direct_ok <- df2$metadata_sample %in% expr_sample_names
  mapped[is.na(mapped) & direct_ok] <- df2$metadata_sample[is.na(mapped) & direct_ok]

  df2$expr_sample <- mapped

  dropped <- df2[is.na(df2$expr_sample) | df2$expr_sample == "", , drop = FALSE]

  kept <- df2[!is.na(df2$expr_sample) & df2$expr_sample != "", , drop = FALSE]
  if (nrow(kept) > 0) {
    # If multiple metadata rows map to the same expression column, keep the first
    # and report the others as dropped/duplicated rather than duplicating columns.
    dup <- duplicated(kept$expr_sample)
    duplicated_rows <- kept[dup, , drop = FALSE]
    kept <- kept[!dup, , drop = FALSE]
    if (nrow(duplicated_rows) > 0) {
      duplicated_rows$drop_reason <- "duplicate metadata rows mapped to the same expression column"
      dropped$drop_reason <- ifelse("drop_reason" %in% colnames(dropped), dropped$drop_reason, "not present in expression matrix")
      dropped <- dplyr::bind_rows(dropped, duplicated_rows)
    }
    kept$sample <- kept$expr_sample
  }

  if (nrow(dropped) > 0 && !"drop_reason" %in% colnames(dropped)) {
    dropped$drop_reason <- "not present in expression matrix after alias/GSM matching"
  }

  list(kept = kept, dropped = dropped)
}


# v11.8.2: explicit expression-column -> metadata-row mapping audit.
# This is the source of truth for debugging sample counts.  It keeps the
# expression matrix columns in order and shows which metadata row/group each
# expression sample actually maps to.  This prevents misleading messages such as
# full-metadata Auto_Group counts being interpreted as DEG-available counts.
build_expression_metadata_map <- function(meta, sample_col, group_col, expr_sample_names) {
  expr_sample_names <- as.character(expr_sample_names)
  if (is.null(meta) || nrow(meta) == 0 || is.null(sample_col) || is.na(sample_col) ||
      !sample_col %in% colnames(meta)) {
    return(data.frame(
      expr_sample = expr_sample_names,
      metadata_row = NA_integer_,
      metadata_sample = NA_character_,
      group = NA_character_,
      geo_accession = NA_character_,
      title = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  metadata_sample <- as.character(meta[[sample_col]])

  if (identical(sample_col, "Expression_Sample_Order")) {
    mapped_expr <- metadata_sample
  } else {
    mapped_expr <- match_meta_values_to_expr_samples(metadata_sample, expr_sample_names)
    direct_ok <- metadata_sample %in% expr_sample_names
    mapped_expr[is.na(mapped_expr) & direct_ok] <- metadata_sample[is.na(mapped_expr) & direct_ok]
  }

  group_value <- if (!is.null(group_col) && !is.na(group_col) && group_col %in% colnames(meta)) {
    clean_group_value(meta[[group_col]])
  } else {
    rep(NA_character_, nrow(meta))
  }

  geo_value <- if ("geo_accession" %in% colnames(meta)) as.character(meta$geo_accession) else rep(NA_character_, nrow(meta))
  title_value <- if ("title" %in% colnames(meta)) as.character(meta$title) else rep(NA_character_, nrow(meta))

  row_map <- data.frame(
    expr_sample = mapped_expr,
    metadata_row = seq_len(nrow(meta)),
    metadata_sample = metadata_sample,
    group = group_value,
    geo_accession = geo_value,
    title = title_value,
    stringsAsFactors = FALSE
  )
  row_map <- row_map[!is.na(row_map$expr_sample) & row_map$expr_sample != "", , drop = FALSE]

  out_list <- lapply(expr_sample_names, function(es) {
    hit <- row_map[row_map$expr_sample == es, , drop = FALSE]
    if (nrow(hit) == 0) {
      return(data.frame(
        expr_sample = es,
        metadata_row = NA_integer_,
        metadata_sample = NA_character_,
        group = NA_character_,
        geo_accession = NA_character_,
        title = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      expr_sample = es,
      metadata_row = paste(hit$metadata_row, collapse = ";"),
      metadata_sample = paste(unique(hit$metadata_sample), collapse = ";"),
      group = paste(unique(hit$group), collapse = ";"),
      geo_accession = paste(unique(hit$geo_accession), collapse = ";"),
      title = paste(unique(hit$title), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(out_list)
  out$mapping_status <- ifelse(is.na(out$metadata_sample) | out$metadata_sample == "", "unmatched", "matched")
  out
}

expression_available_counts_from_map <- function(meta, sample_col, group_col, expr_sample_names) {
  mp <- build_expression_metadata_map(meta, sample_col, group_col, expr_sample_names)
  mp <- mp[mp$mapping_status == "matched" & !is.na(mp$group) & mp$group != "" & mp$group != "NA", , drop = FALSE]
  if (nrow(mp) == 0) {
    return(data.frame(group = character(), n_expression_available = integer(), stringsAsFactors = FALSE))
  }
  # If a duplicated metadata mapping produces combined group strings, keep them visible
  # instead of silently splitting; duplicated mappings should be investigated.
  mp %>% dplyr::count(group, name = "n_expression_available")
}


# v10.4: last-resort order repair for GEO Series Matrix metadata vs supplementary matrix.
# Use only when the full metadata row count equals expression sample count.
# This handles GEO files where expression columns were parsed/renamed but the matrix order
# still follows pData() order. It prevents selected samples from being silently dropped.
repair_deg_samples_by_full_metadata_order <- function(kept, dropped, full_meta, metadata_sample_col, expr_sample_names) {
  kept <- as.data.frame(kept, stringsAsFactors = FALSE)
  dropped <- as.data.frame(dropped, stringsAsFactors = FALSE)
  full_meta <- as.data.frame(full_meta, stringsAsFactors = FALSE)
  expr_sample_names <- as.character(expr_sample_names)

  if (nrow(dropped) == 0) {
    return(list(kept = kept, dropped = dropped, repaired_n = 0, note = "no dropped samples to repair"))
  }
  if (is.null(full_meta) || nrow(full_meta) != length(expr_sample_names)) {
    return(list(kept = kept, dropped = dropped, repaired_n = 0,
                note = "order repair skipped: metadata row count != expression sample count"))
  }
  if (is.null(metadata_sample_col) || is.na(metadata_sample_col) ||
      !(metadata_sample_col %in% colnames(full_meta))) {
    return(list(kept = kept, dropped = dropped, repaired_n = 0,
                note = "order repair skipped: metadata sample column unavailable"))
  }

  order_meta_sample <- as.character(full_meta[[metadata_sample_col]])
  order_map <- data.frame(
    order_meta_sample = order_meta_sample,
    order_meta_norm = normalize_id_for_match(order_meta_sample),
    order_meta_gsm = normalize_id_for_match(extract_first_gsm(order_meta_sample)),
    order_expr_sample = expr_sample_names,
    stringsAsFactors = FALSE
  )

  if (!"metadata_sample" %in% colnames(dropped)) dropped$metadata_sample <- as.character(dropped$sample)
  drop_norm <- normalize_id_for_match(dropped$metadata_sample)
  drop_gsm <- normalize_id_for_match(extract_first_gsm(dropped$metadata_sample))

  m <- match(drop_norm, order_map$order_meta_norm)
  m2 <- match(drop_gsm, order_map$order_meta_gsm)
  m[is.na(m) & !is.na(m2)] <- m2[is.na(m) & !is.na(m2)]

  can_repair <- !is.na(m) & order_map$order_expr_sample[m] %in% expr_sample_names
  if (!any(can_repair)) {
    return(list(kept = kept, dropped = dropped, repaired_n = 0,
                note = "order repair tried but no dropped samples matched full metadata order"))
  }

  repaired <- dropped[can_repair, , drop = FALSE]
  repaired$expr_sample <- order_map$order_expr_sample[m[can_repair]]
  repaired$sample <- repaired$expr_sample
  repaired$drop_reason <- NULL
  repaired$order_repaired <- TRUE

  still_dropped <- dropped[!can_repair, , drop = FALSE]
  if (nrow(still_dropped) > 0) {
    still_dropped$drop_reason <- "not present after alias matching or full-metadata order repair"
  }

  kept2 <- dplyr::bind_rows(kept, repaired)
  dup <- duplicated(kept2$sample)
  if (any(dup)) {
    dup_rows <- kept2[dup, , drop = FALSE]
    dup_rows$drop_reason <- "duplicate after order repair"
    still_dropped <- dplyr::bind_rows(still_dropped, dup_rows)
    kept2 <- kept2[!dup, , drop = FALSE]
  }

  list(
    kept = kept2,
    dropped = still_dropped,
    repaired_n = nrow(repaired),
    note = paste0("order repair used full metadata row order; repaired ", nrow(repaired), " sample(s)")
  )
}

guess_sample_column <- function(meta, raw_sample_names = NULL, match_mode = "clean") {
  cn <- colnames(meta)

  if (!is.null(raw_sample_names)) {
    match_counts <- sapply(cn, function(cc) {
      sample_match_count(meta[[as.character(cc)[1]]], raw_sample_names, mode = match_mode)
    })

    if (length(match_counts) > 0 && max(match_counts, na.rm = TRUE) > 0) {
      return(names(which.max(match_counts)))
    }
  }

  priority <- c(
    "sample", "Sample", "sample_id", "Sample.ID", "sampleid", "sample.name",
    "SampleName", "GSM", "geo_accession", "description", "GEO.Accession", "title"
  )
  sample_col <- priority[priority %in% cn][1]

  if (is.na(sample_col) || length(sample_col) == 0) {
    sample_col <- grep("geo_accession|accession|description|sample|title", cn, value = TRUE, ignore.case = TRUE)[1]
  }

  if (is.na(sample_col) || length(sample_col) == 0) sample_col <- cn[1]

  sample_col
}

infer_auto_group <- function(meta) {
  n <- nrow(meta)
  out <- rep(NA_character_, n)

  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  cn <- colnames(meta)

  # Generic but safer: infer broad groups only from sample-level phenotype-like columns.
  # Do NOT paste the whole metadata row, because relation/summary fields may contain
  # words such as "healthy control" in every row and mislabel COVID/sepsis samples.
  phenotype_cols <- cn[grepl(
    paste0(
      "disease|condition|group|status|phenotype|diagnosis|clinical|",
      "characteristics|source_name|title|sample.type|sample_type|state|class"
    ),
    cn,
    ignore.case = TRUE
  )]

  if (length(phenotype_cols) == 0) {
    phenotype_cols <- cn[grepl("title|source|description|characteristics", cn, ignore.case = TRUE)]
  }

  if (length(phenotype_cols) == 0) return(out)

  txt <- apply(meta[, phenotype_cols, drop = FALSE], 1, function(z) {
    paste(tolower(as.character(z)), collapse = " | ")
  })

  # Specific disease labels first; control is assigned last and only if no disease label matched.
  out[grepl("covid|sars[- ]?cov[- ]?2|coronavirus", txt, ignore.case = TRUE)] <- "COVID-19"
  out[grepl("sepsis|septic|septic shock|septicemia", txt, ignore.case = TRUE)] <- "Sepsis"

  control_hit <- grepl("healthy control|healthy|normal|\\bcontrol\\b|\\bctrl\\b|\\bhc\\b", txt, ignore.case = TRUE)
  out[control_hit & (is.na(out) | out == "")] <- "Healthy"

  out
}

add_auto_group_columns <- function(meta) {
  auto <- infer_auto_group(meta)
  if (sum(!is.na(auto)) >= 2 && length(unique(auto[!is.na(auto)])) >= 2) {
    meta$Auto_Group <- auto
  }
  meta
}



# =========================
# Expression-column metadata fallback
# =========================
# Some GEO supplementary processed matrices use column names such as
# Control1, Control2, Sepsis1... but the GEO Series Matrix metadata still has
# GSM accessions or even a different number of samples. In that case sample-name
# matching should not fail the dataset. Instead, create a minimal metadata table
# directly from expression column names and infer broad groups from those names.
# This is generic: it does not use any GSE ID or tissue-specific keyword.
infer_group_from_expression_sample_names <- function(sample_names) {
  x <- tolower(as.character(sample_names))
  out <- rep(NA_character_, length(x))

  out[grepl('sepsis|septic|infection|infected|case|disease|patient|\bpt\b|\bdis\b', x, ignore.case = TRUE)] <- 'Sepsis'
  out[grepl('healthy|control|ctrl|normal|\bhc\b|\bcon\b|vehicle|sham|baseline', x, ignore.case = TRUE)] <- 'Control'

  # If both a control and a disease word appear, keep disease/case as the safer label.
  out[grepl('sepsis|septic|infection|infected|case|disease|patient', x, ignore.case = TRUE)] <- 'Sepsis'
  out
}

make_expression_sample_metadata <- function(raw_sample_names) {
  raw_sample_names <- as.character(raw_sample_names)
  meta <- data.frame(
    Expression_Sample_Order = raw_sample_names,
    sample = raw_sample_names,
    title = raw_sample_names,
    stringsAsFactors = FALSE
  )
  g <- infer_group_from_expression_sample_names(raw_sample_names)
  if (sum(!is.na(g) & g != '') >= 2 && length(unique(g[!is.na(g) & g != ''])) >= 2) {
    meta$Auto_Group <- g
  }
  meta
}


looks_like_anonymous_expression_sample_names <- function(x) {
  x <- as.character(x)
  x2 <- trimws(x)
  if (length(x2) == 0) return(FALSE)
  prop_anon <- mean(grepl("^(\\.{3}|X|V)?[0-9]+$|^X[0-9]+\\.[0-9]+$", x2, ignore.case = TRUE), na.rm = TRUE)
  prop_no_gsm <- mean(!grepl("GSM[0-9]+", x2, ignore.case = TRUE), na.rm = TRUE)
  prop_anon >= 0.80 && prop_no_gsm >= 0.95
}

has_usable_expression_sample_groups <- function(raw_sample_names) {
  g <- infer_group_from_expression_sample_names(raw_sample_names)
  sum(!is.na(g) & g != '') >= 2 && length(unique(g[!is.na(g) & g != ''])) >= 2
}

metadata_row_bridge_sample_mapping <- function(meta, raw_sample_names, min_match_ratio = 0.80) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  raw_sample_names <- as.character(raw_sample_names)
  n_raw <- length(raw_sample_names)
  if (n_raw == 0 || nrow(meta) == 0) return(NULL)

  alias <- make_expr_sample_alias_table(raw_sample_names)
  alias_long <- unique(data.frame(
    expr_sample = rep(alias$expr_sample, 4),
    alias_norm = c(alias$alias_full_norm, alias$alias_gsm_norm, alias$alias_short_norm, alias$alias_base_norm),
    stringsAsFactors = FALSE
  ))
  alias_long <- alias_long[!is.na(alias_long$alias_norm) & alias_long$alias_norm != "" & nchar(alias_long$alias_norm) >= 4, , drop = FALSE]
  alias_long <- alias_long[!duplicated(alias_long), , drop = FALSE]
  if (nrow(alias_long) == 0) return(NULL)

  # If one alias maps to multiple expression columns, remove it to avoid unsafe matching.
  alias_tab <- table(alias_long$alias_norm)
  alias_long <- alias_long[alias_long$alias_norm %in% names(alias_tab)[alias_tab == 1], , drop = FALSE]
  if (nrow(alias_long) == 0) return(NULL)

  row_text <- apply(meta, 1, function(z) paste(as.character(z), collapse = " | "))
  row_norm <- normalize_id_for_match(row_text)

  assigned <- rep(NA_character_, nrow(meta))
  matched_expr <- character(0)

  for (i in seq_len(nrow(alias_long))) {
    a <- alias_long$alias_norm[i]
    e <- alias_long$expr_sample[i]
    # Exact token/substring match after make.names normalization.
    hit <- grepl(a, row_norm, fixed = TRUE)
    if (sum(hit, na.rm = TRUE) == 1 && !(e %in% matched_expr)) {
      row_i <- which(hit)[1]
      if (is.na(assigned[row_i]) || assigned[row_i] == "") {
        assigned[row_i] <- e
        matched_expr <- c(matched_expr, e)
      }
    }
  }

  matched_n <- length(unique(assigned[!is.na(assigned) & assigned != ""]))
  ratio <- matched_n / max(1, n_raw)
  if (matched_n >= 2 && ratio >= min_match_ratio) {
    meta$Expression_Sample_Order <- assigned
    return(list(
      meta = meta,
      sample_col = "Expression_Sample_Order",
      mode_used = "metadata_row_bridge",
      matched_n = matched_n,
      warning = paste0(
        "No direct GSM/sample-column match, but metadata row-text bridge matched ",
        matched_n, "/", n_raw,
        " expression samples using aliases from all metadata fields. Please verify once."
      )
    ))
  }

  NULL
}

choose_sample_mapping <- function(meta, raw_sample_names, requested_mode = "auto_robust") {
  raw_sample_names <- as.character(raw_sample_names)
  n_raw <- length(raw_sample_names)
  n_meta <- nrow(meta)

  strict_col <- guess_sample_column(meta, raw_sample_names, match_mode = "strict")
  strict_n <- sample_match_count(meta[[strict_col]], raw_sample_names, mode = "strict")

  clean_col <- guess_sample_column(meta, raw_sample_names, match_mode = "clean")
  clean_n <- sample_match_count(meta[[clean_col]], raw_sample_names, mode = "clean")

  # Forced order mode.
  if (requested_mode == "order") {
    if (n_meta == n_raw) {
      meta$Expression_Sample_Order <- raw_sample_names
      return(list(
        meta = meta,
        sample_col = "Expression_Sample_Order",
        mode_used = "order",
        matched_n = n_raw,
        warning = "Forced order matching was used because user selected order mode."
      ))
    } else {
      return(list(
        meta = meta,
        sample_col = clean_col,
        mode_used = "clean_fallback_order_failed",
        matched_n = clean_n,
        warning = paste0("Order matching was requested, but metadata rows (", n_meta,
                         ") != expression samples (", n_raw, "). Used cleaned-name matching instead.")
      ))
    }
  }

  # Strict only.
  if (requested_mode == "strict") {
    return(list(
      meta = meta,
      sample_col = strict_col,
      mode_used = "strict",
      matched_n = strict_n,
      warning = NULL
    ))
  }

  # Clean only.
  if (requested_mode == "clean") {
    return(list(
      meta = meta,
      sample_col = clean_col,
      mode_used = "clean",
      matched_n = clean_n,
      warning = NULL
    ))
  }

  # Auto robust v11.1:
  # Prefer true sample-name/GSM matching.
  # Only use row-order fallback when expression column names are anonymous Excel placeholders
  # (for example ...2, ...3, X2, V2) AND metadata row count exactly equals expression sample count.
  # This fixes processed matrices where Excel headers were not recoverable, while still preventing
  # dangerous order matching for real but discordant GSM columns such as GSE63042 raw counts.
  clean_ratio <- ifelse(n_raw > 0, clean_n / n_raw, 0)

  # v11.7.7: if expression names are internal IDs (e.g. SZ329mes) and GEO metadata
  # has GSM rows plus the internal ID hidden in title/source/characteristics, scan the
  # full metadata row text and map only uniquely matched rows. This avoids unsafe order fallback.
  if (clean_ratio < 0.80) {
    row_bridge <- metadata_row_bridge_sample_mapping(meta, raw_sample_names, min_match_ratio = 0.80)
    if (!is.null(row_bridge)) return(row_bridge)
  }

  if (clean_ratio < 0.80 && n_meta == n_raw && looks_like_anonymous_expression_sample_names(raw_sample_names)) {
    meta$Expression_Sample_Order <- raw_sample_names
    return(list(
      meta = meta,
      sample_col = "Expression_Sample_Order",
      mode_used = "safe_order_anonymous_columns",
      matched_n = n_raw,
      warning = paste0(
        "Cleaned/GSM matching was low (", clean_n, "/", n_raw, "), but expression sample names look anonymous ",
        "(...2/...3/X2/V2) and metadata rows exactly match expression samples. ",
        "Used safe row-order matching for this processed matrix. Please verify group order once."
      )
    ))
  }

  # Final fallback for processed matrices whose expression column names do not match GEO GSM metadata.
  # If sample names themselves contain enough group information, use them as metadata.
  # This prevents wrong dashboards such as 12 expression samples but 28 GEO metadata rows.
  if (clean_ratio < 0.80 && has_usable_expression_sample_groups(raw_sample_names)) {
    expr_meta <- make_expression_sample_metadata(raw_sample_names)
    return(list(
      meta = expr_meta,
      sample_col = "Expression_Sample_Order",
      mode_used = "expression_sample_name_metadata",
      matched_n = n_raw,
      warning = paste0(
        "No reliable match between GEO Series metadata and expression matrix columns (",
        clean_n, "/", n_raw, "). Created metadata from expression column names. ",
        "Please confirm group labels inferred from sample names."
      )
    ))
  }

  return(list(
    meta = meta,
    sample_col = clean_col,
    mode_used = ifelse(clean_ratio >= 0.80, "clean", "clean_low_confidence"),
    matched_n = clean_n,
    warning = if (clean_ratio >= 0.80) {
      paste0("Cleaned/GSM sample-name matching was used. Matched ", clean_n, "/", n_raw,
             " expression samples. Unsafe order fallback is disabled.")
    } else {
      paste0(
        "Low-confidence sample matching: only ", clean_n, "/", n_raw,
        " samples matched. Row-order fallback was NOT used because expression sample names are not anonymous ",
        "or metadata rows (", n_meta, ") do not exactly equal expression samples (", n_raw, ")."
      )
    }
  ))
}


# =========================
# Manual group merge helper
# =========================
get_analysis_group_names <- function(enable_merge,
                                     groupA, groupB,
                                     mergeA_name, mergeB_name) {
  if (isTRUE(enable_merge)) {
    a <- trimws(as.character(mergeA_name))
    b <- trimws(as.character(mergeB_name))
    if (is.na(a) || a == "") a <- "Group_A"
    if (is.na(b) || b == "") b <- "Group_B"
    return(list(A = a, B = b))
  }

  list(A = groupA, B = groupB)
}

apply_manual_group_merge <- function(meta2, enable_merge,
                                     mergeA_groups, mergeB_groups,
                                     mergeA_name = "Control",
                                     mergeB_name = "Case") {
  if (is.null(meta2) || nrow(meta2) == 0) return(meta2)

  out <- meta2
  out$group_original <- out$group

  if (!isTRUE(enable_merge)) {
    out$group_for_analysis <- out$group_original
    out$group <- out$group_for_analysis
    return(out)
  }

  mergeA_groups <- as.character(mergeA_groups)
  mergeB_groups <- as.character(mergeB_groups)
  mergeA_groups <- mergeA_groups[!is.na(mergeA_groups) & mergeA_groups != ""]
  mergeB_groups <- mergeB_groups[!is.na(mergeB_groups) & mergeB_groups != ""]

  if (length(mergeA_groups) == 0 || length(mergeB_groups) == 0) {
    out$group_for_analysis <- NA_character_
    return(out[0, , drop = FALSE])
  }

  # If a group is selected in both A and B, keep it only in A.
  overlap <- intersect(mergeA_groups, mergeB_groups)
  if (length(overlap) > 0) {
    mergeB_groups <- setdiff(mergeB_groups, overlap)
  }

  if (length(mergeB_groups) == 0) {
    out$group_for_analysis <- NA_character_
    return(out[0, , drop = FALSE])
  }

  names <- get_analysis_group_names(TRUE, NULL, NULL, mergeA_name, mergeB_name)
  mergeA_name <- names$A
  mergeB_name <- names$B

  out$group_for_analysis <- NA_character_
  out$group_for_analysis[out$group_original %in% mergeA_groups] <- mergeA_name
  out$group_for_analysis[out$group_original %in% mergeB_groups] <- mergeB_name

  out <- out[!is.na(out$group_for_analysis), , drop = FALSE]
  out$group <- out$group_for_analysis
  out
}


# =========================
# 6. Differential analysis engines
# =========================
prepare_expression_for_limma <- function(raw, samples, det) {
  mat <- raw %>% select(Symbol, all_of(samples)) %>% column_to_rownames("Symbol") %>% as.matrix()
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA

  transform_used <- "none"

  # v11.5 runtime scale guard: candidate labels may come from older cache or
  # ambiguous submitted matrices. Re-check the numeric scale here, immediately
  # before limma, and override stale log/microarray labels when needed.
  vals_guard <- as.numeric(mat)
  vals_guard <- vals_guard[is.finite(vals_guard) & !is.na(vals_guard)]
  q99_guard <- if (length(vals_guard) > 100) as.numeric(stats::quantile(vals_guard, 0.99, na.rm = TRUE)) else NA_real_
  max_guard <- if (length(vals_guard) > 100) max(vals_guard, na.rm = TRUE) else NA_real_
  min_guard <- if (length(vals_guard) > 100) min(vals_guard, na.rm = TRUE) else NA_real_
  large_positive_guard <- is.finite(q99_guard) && is.finite(max_guard) && is.finite(min_guard) &&
    min_guard >= 0 && (q99_guard > 25 || max_guard > 150)
  # v11.7 hard guard: any non-count, large positive submitted matrix is treated
  # as LINEAR normalized expression and must be log2(x+1) transformed before limma.
  # Example: GSE310929_AllSampleExpressionSubmitted.tsv has S100A8 values >100.
  if (!is.null(det$type) && det$type != "raw_count" && large_positive_guard) {
    det$type <- "normalized_expression_FPKM_TPM_or_similar"
    det$transform_for_limma <- "log2(x + 1)"
    det$reason <- paste0("v11.7 hard runtime override: large positive non-count submitted matrix; q99=",
                         round(q99_guard, 2), ", max=", round(max_guard, 2),
                         "; limma will use log2(x+1)")
  }

  # v11.7 final transform rule:
  # raw_count -> log2(CPM + 1) fallback for ordinary limma
  # normalized FPKM/TPM/RPM-like -> log2(x + 1)
  # already log/microarray -> none
  if (!is.null(det$type) && det$type == "raw_count") {
    lib_size <- colSums(pmax(mat, 0), na.rm = TRUE)
    lib_size[lib_size <= 0 | !is.finite(lib_size)] <- median(lib_size[lib_size > 0 & is.finite(lib_size)], na.rm = TRUE)
    cpm <- t(t(pmax(mat, 0)) / lib_size * 1e6)
    mat <- log2(cpm + 1)
    transform_used <- "log2(CPM + 1) fallback from raw count"
  } else if (!is.null(det$type) && det$type == "normalized_expression_FPKM_TPM_or_similar") {
    mat <- log2(pmax(mat, 0) + 1)
    transform_used <- "log2(x + 1) for normalized/RPM/TPM/FPKM matrix"
  } else if (!is.null(det$type) && det$type == "log_or_microarray_expression") {
    transform_used <- "none; input already appears log-scale/microarray-like"
  } else if (!is.null(det$transform_for_limma) && det$transform_for_limma == "log2(x + 1)") {
    mat <- log2(pmax(mat, 0) + 1)
    transform_used <- "log2(x + 1) fallback for unknown positive matrix"
  }

  attr(mat, "limma_transform_used") <- transform_used
  attr(mat, "limma_runtime_type") <- if (!is.null(det$type)) det$type else NA_character_
  attr(mat, "limma_runtime_reason") <- if (!is.null(det$reason)) det$reason else NA_character_
  mat
}

validate_deg_inputs <- function(raw, df2, groupA, groupB) {
  if (is.null(raw) || ncol(raw) < 3 || nrow(raw) < 2) stop("Expression matrix is empty or has fewer than two samples.")
  if (!all(c("sample", "group") %in% colnames(df2))) stop("Metadata must contain sample and group columns.")
  if (anyDuplicated(colnames(raw)[-1])) stop("Expression matrix contains duplicated sample names.")
  if (anyDuplicated(as.character(df2$sample))) stop("Metadata contains duplicated matched sample names.")
  if (anyNA(df2$sample) || anyNA(df2$group) || any(trimws(as.character(df2$sample)) == "") ||
      any(trimws(as.character(df2$group)) == "")) stop("Matched metadata contains missing sample/group values.")
  missing_samples <- setdiff(as.character(df2$sample), colnames(raw)[-1])
  if (length(missing_samples) > 0) {
    stop("Metadata/expression alignment failed; missing expression samples: ", paste(head(missing_samples, 10), collapse = ", "))
  }
  unexpected_groups <- setdiff(unique(as.character(df2$group)), c(groupA, groupB))
  if (length(unexpected_groups) > 0) stop("Unexpected groups remained after DEG filtering: ", paste(unexpected_groups, collapse = ", "))
  group_n <- table(factor(as.character(df2$group), levels = c(groupA, groupB)))
  if (any(as.integer(group_n) < 2)) stop("Each DEG group must have at least two matched expression samples.")
  symbols <- as.character(raw[[1]])
  if (anyNA(symbols) || any(trimws(symbols) == "") || anyDuplicated(symbols)) {
    stop("Gene Symbol column must be non-missing and unique before differential analysis.")
  }
  invisible(TRUE)
}

run_limma_analysis <- function(raw, df2, groupA, groupB, det) {
  validate_deg_inputs(raw, df2, groupA, groupB)
  samples <- df2$sample

  # v11.7: do NOT rely only on cached expr_detect labels here.
  # Rebuild the numeric limma matrix from the selected active expression source,
  # inspect its real scale, and decide the transform immediately before fitting.
  mat0 <- raw %>% select(Symbol, all_of(samples)) %>% column_to_rownames("Symbol") %>% as.matrix()
  suppressWarnings(storage.mode(mat0) <- "numeric")
  mat0[!is.finite(mat0)] <- NA

  vals_guard <- as.numeric(mat0)
  vals_guard <- vals_guard[is.finite(vals_guard) & !is.na(vals_guard)]
  q50_guard <- if (length(vals_guard) > 100) as.numeric(stats::quantile(vals_guard, 0.50, na.rm = TRUE)) else NA_real_
  q90_guard <- if (length(vals_guard) > 100) as.numeric(stats::quantile(vals_guard, 0.90, na.rm = TRUE)) else NA_real_
  q99_guard <- if (length(vals_guard) > 100) as.numeric(stats::quantile(vals_guard, 0.99, na.rm = TRUE)) else NA_real_
  max_guard <- if (length(vals_guard) > 100) max(vals_guard, na.rm = TRUE) else NA_real_
  min_guard <- if (length(vals_guard) > 100) min(vals_guard, na.rm = TRUE) else NA_real_
  int_prop_guard <- if (length(vals_guard) > 100) mean(abs(vals_guard - round(vals_guard)) < 1e-6, na.rm = TRUE) else NA_real_

  # Never infer raw counts again from integer-likeness at fit time.  The central
  # detector already used source evidence plus stricter distribution checks.
  det_type_original <- if (!is.null(det$type)) det$type else "unknown"
  raw_like_guard <- identical(det_type_original, "raw_count")

  # Linear normalized guard: positive non-count matrices with large values are NOT log-scale,
  # even if cached detection says log_or_microarray_expression. Example GSE310929:
  # median ~3, q99 ~37, max >300, S100A8 values >100.
  linear_normalized_guard <- is.finite(min_guard) && is.finite(q99_guard) && is.finite(max_guard) &&
    min_guard >= 0 && !isTRUE(raw_like_guard) && (q99_guard > 25 || max_guard > 150)

  det_reason_original <- if (!is.null(det$reason)) det$reason else "unknown"
  runtime_reason <- paste0(
    "v11.7 runtime scale check: min=", round(min_guard, 3),
    ", median=", round(q50_guard, 3),
    ", q90=", round(q90_guard, 3),
    ", q99=", round(q99_guard, 3),
    ", max=", round(max_guard, 3),
    ", integer_prop=", round(int_prop_guard, 3),
    "; cached_type=", det_type_original
  )

  mat <- mat0
  limma_runtime_type <- det_type_original
  limma_runtime_reason <- runtime_reason
  limma_transform_used <- "none"

  if (raw_like_guard || identical(det_type_original, "raw_count")) {
    lib_size <- colSums(pmax(mat0, 0), na.rm = TRUE)
    ok_lib <- lib_size > 0 & is.finite(lib_size)
    if (!any(ok_lib)) stop("All library sizes are zero/non-finite; cannot run limma on count-like matrix.")
    lib_size[!ok_lib] <- median(lib_size[ok_lib], na.rm = TRUE)
    cpm <- t(t(pmax(mat0, 0)) / lib_size * 1e6)
    mat <- log2(cpm + 1)
    limma_runtime_type <- "raw_count"
    limma_transform_used <- "log2(CPM + 1) fallback from raw count"
  } else if (linear_normalized_guard || identical(det_type_original, "normalized_expression_FPKM_TPM_or_similar")) {
    mat <- log2(pmax(mat0, 0) + 1)
    limma_runtime_type <- "normalized_expression_FPKM_TPM_or_similar"
    limma_transform_used <- "log2(x + 1) for linear normalized/RPM/TPM/FPKM/submitted matrix"
  } else if (identical(det_type_original, "log_or_microarray_expression")) {
    mat <- mat0
    limma_runtime_type <- "log_or_microarray_expression"
    limma_transform_used <- "none; input appears already log-scale/microarray-like"
  } else {
    # conservative fallback for unknown non-negative matrix
    if (is.finite(min_guard) && min_guard >= 0) {
      mat <- log2(pmax(mat0, 0) + 1)
      limma_runtime_type <- "normalized_expression_FPKM_TPM_or_similar"
      limma_transform_used <- "log2(x + 1) fallback for unknown positive matrix"
    } else {
      mat <- mat0
      limma_runtime_type <- "log_or_microarray_expression"
      limma_transform_used <- "none fallback for unknown matrix with negative values"
    }
  }

  run_fit_once <- function(expr_mat_input) {
    expr_mat_all0 <- expr_mat_input[rowSums(is.finite(expr_mat_input)) == ncol(expr_mat_input), , drop = FALSE]
    keep0 <- rowSums(is.finite(expr_mat_all0)) == ncol(expr_mat_all0) & apply(expr_mat_all0, 1, stats::sd, na.rm = TRUE) > 0
    expr_mat0 <- expr_mat_all0[keep0, , drop = FALSE]
    df_fit0 <- df2[match(colnames(expr_mat0), df2$sample), , drop = FALSE]
    group0 <- factor(df_fit0$group, levels = c(groupA, groupB))
    design0 <- model.matrix(~ 0 + group0)
    colnames(design0) <- make.names(levels(group0))
    contrast_name0 <- paste0(make.names(groupB), "-", make.names(groupA))
    cont0 <- limma::makeContrasts(contrasts = contrast_name0, levels = design0)
    fit0 <- limma::lmFit(expr_mat0, design0)
    fit20 <- limma::contrasts.fit(fit0, cont0)
    fit20 <- limma::eBayes(fit20)
    tt0 <- limma::topTable(fit20, number = Inf, sort.by = "P") %>% rownames_to_column("Symbol")
    baseMean0 <- rowMeans(expr_mat0[tt0$Symbol, , drop = FALSE], na.rm = TRUE)
    res0 <- data.frame(
      Symbol = tt0$Symbol,
      baseMean = as.numeric(baseMean0),
      log2FoldChange = tt0$logFC,
      lfcSE = NA_real_,
      stat = tt0$t,
      pvalue = tt0$P.Value,
      padj = tt0$adj.P.Val,
      method = "limma",
      stringsAsFactors = FALSE
    )
    list(res_df = res0, expr_mat = expr_mat0, expr_mat_all = expr_mat_all0, group = group0, df_fit = df_fit0)
  }

  fit_obj <- run_fit_once(mat)
  max_abs_lfc <- suppressWarnings(max(abs(fit_obj$res_df$log2FoldChange), na.rm = TRUE))

  # Final automatic rescue: if limma still returns impossible logFC from a non-count
  # submitted matrix, rerun once using log2(x+1) on the original selected matrix.
  if (is.finite(max_abs_lfc) && max_abs_lfc > 20 && !identical(limma_runtime_type, "raw_count")) {
    mat_rescue <- log2(pmax(mat0, 0) + 1)
    fit_rescue <- run_fit_once(mat_rescue)
    max_abs_rescue <- suppressWarnings(max(abs(fit_rescue$res_df$log2FoldChange), na.rm = TRUE))
    if (is.finite(max_abs_rescue) && max_abs_rescue < max_abs_lfc) {
      fit_obj <- fit_rescue
      limma_runtime_type <- "normalized_expression_FPKM_TPM_or_similar"
      limma_transform_used <- paste0(
        "log2(x + 1) AUTO-RESCUE after abnormal limma logFC; previous max |log2FC|=",
        round(max_abs_lfc, 2), ", rescued max |log2FC|=", round(max_abs_rescue, 2)
      )
      max_abs_lfc <- max_abs_rescue
    } else {
      limma_transform_used <- paste0(
        limma_transform_used,
        "\nWARNING: max |log2FC| = ", round(max_abs_lfc, 2),
        ". Rescue log2(x+1) did not improve enough; verify expression source/mapping."
      )
    }
  }

  if (is.finite(max_abs_lfc) && max_abs_lfc > 20) {
    limma_transform_used <- paste0(
      limma_transform_used,
      "\nWARNING: max |log2FC| = ", round(max_abs_lfc, 2),
      ". This is unusually large; please verify expression scale / source selection."
    )
  }

  expr_mat <- fit_obj$expr_mat
  expr_mat_all <- fit_obj$expr_mat_all
  group <- fit_obj$group
  df_fit <- fit_obj$df_fit

  list(
    method = "limma",
    res_df = fit_obj$res_df,
    plot_expr_mat = expr_mat,
    plot_expr_mat_all = expr_mat_all,
    limma_transform_used = limma_transform_used,
    coldata = data.frame(sample = df_fit$sample, group = group, stringsAsFactors = FALSE),
    message = paste0(
      "limma完成。\n",
      "Contrast: ", groupB, " vs ", groupA, "\n",
      "Samples: ", paste(table(group), collapse = " vs "), "\n",
      "Genes after filtering: ", nrow(expr_mat), "\n",
      "Input判断：", limma_runtime_type, "\n",
      "输入判断原因：", limma_runtime_reason, "\n",
      "limma输入转换：", limma_transform_used
    )
  )
}


run_limma_voom_analysis <- function(raw, df2, groupA, groupB, det) {
  validate_deg_inputs(raw, df2, groupA, groupB)
  samples <- df2$sample

  count_mat <- raw %>% select(Symbol, all_of(samples)) %>% column_to_rownames("Symbol") %>% as.matrix()
  storage.mode(count_mat) <- "numeric"
  count_mat[is.na(count_mat)] <- 0
  count_mat[count_mat < 0] <- 0
  count_mat <- round(count_mat)

  # limma-voom is designed for raw count matrices. Keep genes with a minimal count signal.
  keep <- rowSums(count_mat) > 0
  count_mat <- count_mat[keep, , drop = FALSE]

  df2 <- df2[match(colnames(count_mat), df2$sample), , drop = FALSE]
  group <- factor(df2$group, levels = c(groupA, groupB))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- make.names(levels(group))

  # edgeR is commonly used with voom to compute library-size normalization factors.
  # It is installed automatically here if missing, without changing the rest of the app.
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR", ask = FALSE, update = FALSE)
  }

  dge <- edgeR::DGEList(counts = count_mat, group = group)
  keep2 <- edgeR::filterByExpr(dge, design = design)
  dge <- dge[keep2, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  v <- limma::voom(dge, design, plot = FALSE)

  contrast_name <- paste0(make.names(groupB), "-", make.names(groupA))
  cont <- limma::makeContrasts(contrasts = contrast_name, levels = design)
  fit <- limma::lmFit(v, design)
  fit2 <- limma::contrasts.fit(fit, cont)
  fit2 <- limma::eBayes(fit2)

  tt <- limma::topTable(fit2, number = Inf, sort.by = "P") %>% rownames_to_column("Symbol")
  baseMean <- rowMeans(v$E[tt$Symbol, , drop = FALSE], na.rm = TRUE)

  res_df <- data.frame(
    Symbol = tt$Symbol,
    baseMean = as.numeric(baseMean),
    log2FoldChange = tt$logFC,
    lfcSE = NA_real_,
    stat = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    method = "limma-voom",
    stringsAsFactors = FALSE
  )

  list(
    method = "limma-voom",
    res_df = res_df,
    plot_expr_mat = v$E,
    voom = v,
    coldata = data.frame(sample = df2$sample, group = group, stringsAsFactors = FALSE),
    message = paste0(
      "limma-voom完成。\n",
      "Contrast: ", groupB, " vs ", groupA, "\n",
      "Samples: ", paste(table(group), collapse = " vs "), "\n",
      "Genes after filtering: ", nrow(v$E), "\n",
      "Input判断：", det$type, "\n",
      "limma输入转换：raw count -> TMM normalization -> voom precision weights"
    )
  )
}

run_deseq2_analysis <- function(raw, df2, groupA, groupB, det) {
  validate_deg_inputs(raw, df2, groupA, groupB)
  samples <- df2$sample
  count_mat <- raw %>% select(Symbol, all_of(samples)) %>% column_to_rownames("Symbol") %>% as.matrix()
  storage.mode(count_mat) <- "numeric"
  count_mat[is.na(count_mat)] <- 0
  count_mat <- round(count_mat)
  count_mat[count_mat < 0] <- 0
  count_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]

  df2 <- df2[match(colnames(count_mat), df2$sample), , drop = FALSE]
  coldata <- data.frame(
    sample = df2$sample,
    group = factor(df2$group, levels = c(groupA, groupB))
  )
  rownames(coldata) <- coldata$sample

  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
  # Many GEO count matrices have at least one zero for every gene.
  # Use poscounts size-factor estimation to avoid:
  # "every gene contains at least one zero, cannot compute log geometric means"
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", groupB, groupA))
  res_df <- as.data.frame(res) %>% rownames_to_column("Symbol") %>% arrange(padj)
  res_df$method <- "DESeq2"
  norm_counts <- counts(dds, normalized = TRUE)
  plot_expr_mat <- log2(norm_counts + 1)

  list(
    method = "DESeq2",
    dds = dds,
    res_df = res_df,
    norm_counts = norm_counts,
    plot_expr_mat = plot_expr_mat,
    coldata = coldata,
    message = paste0(
      "DESeq2完成。\n",
      "Contrast: ", groupB, " vs ", groupA, "\n",
      "Samples: ", paste(table(coldata$group), collapse = " vs "), "\n",
      "Genes after filtering: ", nrow(dds), "\n",
      "Input判断：", det$type
    )
  )
}

run_auto_differential <- function(raw, df2, groupA, groupB, user_method = "auto", det = NULL) {
  if (is.null(det)) det <- detect_expression_type(raw[, -1, drop = FALSE], source_hint = attr(raw, "source_hint"))
  if (user_method == "auto") {
    method <- det$method
  } else {
    method <- user_method
  }
  if (method == "DESeq2") {
    if (det$type != "raw_count" && user_method == "DESeq2") {
      warning("你强制使用DESeq2，但输入不像raw count。结果只适合参考。")
    }
    run_deseq2_analysis(raw, df2, groupA, groupB, det)
  } else if (method == "limma_voom") {
    if (det$type != "raw_count") {
      warning("limma-voom适合raw count；当前输入不像raw count。建议改用普通limma。")
    }
    run_limma_voom_analysis(raw, df2, groupA, groupB, det)
  } else {
    run_limma_analysis(raw, df2, groupA, groupB, det)
  }
}


# =========================
# 6.5 Plot helper functions
# =========================
make_heatmap_palette <- function(palette_name = "blue_white_red", n = 101) {
  if (palette_name == "green_black_red") {
    return(colorRampPalette(c("#00A087", "black", "#E64B35"))(n))
  }
  if (palette_name == "purple_white_orange") {
    return(colorRampPalette(c("#7E57C2", "white", "#F39B7F"))(n))
  }
  colorRampPalette(c("#0072B5", "white", "#E64B35"))(n)
}

prepare_heatmap_data <- function(diff_res, expr_mat, coldata,
                                 top_n = 50,
                                 padj_cutoff = 0.05,
                                 logfc_cutoff = 0,
                                 rank_by = "padj") {
  if (is.null(diff_res) || is.null(expr_mat) || is.null(coldata)) return(NULL)

  res <- diff_res %>%
    filter(!is.na(Symbol), Symbol %in% rownames(expr_mat))

  if ("padj" %in% colnames(res)) {
    res <- res %>% filter(is.na(padj) | padj <= padj_cutoff)
  }

  if ("log2FoldChange" %in% colnames(res)) {
    res <- res %>% filter(!is.na(log2FoldChange), abs(log2FoldChange) >= logfc_cutoff)
  }

  if (nrow(res) == 0) return(NULL)

  if (rank_by == "abslogfc") {
    res <- res %>% arrange(desc(abs(log2FoldChange)))
  } else if (rank_by == "pvalue") {
    res <- res %>% arrange(pvalue)
  } else {
    res <- res %>% arrange(padj)
  }

  genes <- unique(head(res$Symbol, top_n))
  genes <- genes[genes %in% rownames(expr_mat)]

  if (length(genes) < 2) return(NULL)

  mat <- expr_mat[genes, coldata$sample, drop = FALSE]

  annotation_col <- data.frame(Group = coldata$group)
  rownames(annotation_col) <- coldata$sample

  list(mat = mat, annotation_col = annotation_col, genes = genes)
}



make_deg_summary <- function(dobj, padj_cutoff = 0.05, logfc_cutoff = 1) {
  # v11.8.15 logic fix:
  # DEG Summary uses padj/FDR only and is no longer affected by the UI |log2FC| cutoff.
  # The |log2FC| slider is kept for Search DEG / Top tables / volcano / heatmap filtering.
  # This matches common DESeq2 reporting and avoids hiding genes such as PTK2
  # when padj < 0.05 but |log2FC| is slightly below the display cutoff.
  if (is.null(dobj) || !is.null(dobj$error) || is.null(dobj$res_df)) {
    return(data.frame(
      Metric = c("Significant DEGs", "Upregulated", "Downregulated", "Top up", "Top down"),
      Value = c(NA, NA, NA, NA, NA),
      stringsAsFactors = FALSE
    ))
  }

  # Convert to ordinary data.frame to avoid Bioconductor Rle/list-column issues.
  res <- as.data.frame(dobj$res_df, stringsAsFactors = FALSE)

  need_cols <- c("Symbol", "log2FoldChange", "padj")
  if (!all(need_cols %in% colnames(res))) {
    return(data.frame(
      Metric = c("Significant DEGs", "Upregulated", "Downregulated", "Top up", "Top down"),
      Value = c(NA, NA, NA, "Missing required columns", "Missing required columns"),
      stringsAsFactors = FALSE
    ))
  }

  res$Symbol <- as.character(res$Symbol)
  res$log2FoldChange <- suppressWarnings(as.numeric(res$log2FoldChange))
  res$padj <- suppressWarnings(as.numeric(res$padj))

  res <- res[
    !is.na(res$log2FoldChange) &
      !is.na(res$padj),
    ,
    drop = FALSE
  ]

  # IMPORTANT: Summary counts are based on adjusted P value only.
  # Do not apply the UI |log2FC| cutoff here.
  sig <- res[
    res$padj < padj_cutoff,
    ,
    drop = FALSE
  ]

  up <- sig[sig$log2FoldChange > 0, , drop = FALSE]
  down <- sig[sig$log2FoldChange < 0, , drop = FALSE]

  top_up <- "None"
  if (nrow(up) > 0) {
    up <- up[order(up$padj, -up$log2FoldChange), , drop = FALSE]
    top_up <- as.character(up$Symbol[1])
  }

  top_down <- "None"
  if (nrow(down) > 0) {
    down <- down[order(down$padj, down$log2FoldChange), , drop = FALSE]
    top_down <- as.character(down$Symbol[1])
  }

  data.frame(
    Metric = c("Significant DEGs", "Upregulated", "Downregulated", "Top up", "Top down"),
    Value = c(as.character(nrow(sig)), as.character(nrow(up)), as.character(nrow(down)), top_up, top_down),
    stringsAsFactors = FALSE
  )
}

get_top_deg_table <- function(dobj, direction = "up", top_n = 20, padj_cutoff = 0.05, logfc_cutoff = 1) {
  if (is.null(dobj) || !is.null(dobj$error) || is.null(dobj$res_df)) return(NULL)

  res <- as.data.frame(dobj$res_df, stringsAsFactors = FALSE)

  need_cols <- c("Symbol", "method", "baseMean", "log2FoldChange", "pvalue", "padj")
  keep_cols <- intersect(need_cols, colnames(res))
  if (!all(c("Symbol", "log2FoldChange", "padj") %in% colnames(res))) return(NULL)

  res$Symbol <- as.character(res$Symbol)
  res$log2FoldChange <- suppressWarnings(as.numeric(res$log2FoldChange))
  res$padj <- suppressWarnings(as.numeric(res$padj))
  if ("pvalue" %in% colnames(res)) res$pvalue <- suppressWarnings(as.numeric(res$pvalue))
  if ("baseMean" %in% colnames(res)) res$baseMean <- suppressWarnings(as.numeric(res$baseMean))

  res <- res[
    !is.na(res$log2FoldChange) &
      !is.na(res$padj) &
      res$padj < padj_cutoff &
      abs(res$log2FoldChange) >= logfc_cutoff,
    ,
    drop = FALSE
  ]

  if (direction == "up") {
    res <- res[res$log2FoldChange > 0, , drop = FALSE]
    res <- res[order(res$padj, -res$log2FoldChange), , drop = FALSE]
  } else {
    res <- res[res$log2FoldChange < 0, , drop = FALSE]
    res <- res[order(res$padj, res$log2FoldChange), , drop = FALSE]
  }

  if (nrow(res) == 0) return(NULL)

  res <- head(res, top_n)
  res$fold_change <- 2^res$log2FoldChange
  res$log2FoldChange <- round(res$log2FoldChange, 3)
  res$fold_change <- round(res$fold_change, 3)
  if ("pvalue" %in% colnames(res)) res$pvalue <- signif(res$pvalue, 3)
  res$padj <- signif(res$padj, 3)
  if ("baseMean" %in% colnames(res)) res$baseMean <- round(res$baseMean, 3)

  show_cols <- c("Symbol", "method", "baseMean", "log2FoldChange", "fold_change", "pvalue", "padj")
  res[, intersect(show_cols, colnames(res)), drop = FALSE]
}

search_deg_table <- function(dobj, query) {
  if (is.null(dobj) || !is.null(dobj$error) || is.null(dobj$res_df)) return(NULL)
  if (is.null(query) || trimws(query) == "") return(NULL)

  res <- as.data.frame(dobj$res_df, stringsAsFactors = FALSE)
  if (!"Symbol" %in% colnames(res)) return(NULL)

  q <- toupper(trimws(query))
  res$Symbol <- as.character(res$Symbol)

  res <- res[grepl(q, toupper(res$Symbol), fixed = TRUE), , drop = FALSE]
  if (nrow(res) == 0) return(NULL)

  if ("padj" %in% colnames(res)) {
    res$padj <- suppressWarnings(as.numeric(res$padj))
    res <- res[order(res$padj), , drop = FALSE]
  }

  if ("log2FoldChange" %in% colnames(res)) {
    res$log2FoldChange <- suppressWarnings(as.numeric(res$log2FoldChange))
    res$fold_change <- 2^res$log2FoldChange
    res$log2FoldChange <- round(res$log2FoldChange, 3)
    res$fold_change <- round(res$fold_change, 3)
  }

  if ("pvalue" %in% colnames(res)) {
    res$pvalue <- signif(suppressWarnings(as.numeric(res$pvalue)), 3)
  }

  if ("padj" %in% colnames(res)) {
    res$padj <- signif(suppressWarnings(as.numeric(res$padj)), 3)
  }

  if ("baseMean" %in% colnames(res)) {
    res$baseMean <- round(suppressWarnings(as.numeric(res$baseMean)), 3)
  }

  show_cols <- c("Symbol", "method", "baseMean", "log2FoldChange", "fold_change", "pvalue", "padj")
  res[, intersect(show_cols, colnames(res)), drop = FALSE]
}


add_logfc_arrow_html <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!"log2FoldChange" %in% colnames(df)) return(df)

  logfc_num <- suppressWarnings(as.numeric(df$log2FoldChange))

  arrow <- ifelse(
    is.na(logfc_num),
    "",
    ifelse(
      logfc_num > 0,
      " <span style='color:#E64B35;font-weight:bold;'>↑</span>",
      ifelse(logfc_num < 0, " <span style='color:#0072B5;font-weight:bold;'>↓</span>", " →")
    )
  )

  df$log2FC_direction <- paste0(
    sprintf("%.3f", logfc_num),
    arrow
  )

  # Put the decorated column beside/inside the useful position.
  show_cols <- colnames(df)
  show_cols <- show_cols[show_cols != "log2FoldChange"]

  if ("fold_change" %in% show_cols) {
    pos <- match("fold_change", show_cols)
    show_cols <- append(show_cols, "log2FC_direction", after = max(0, pos - 1))
  } else {
    show_cols <- c("log2FC_direction", show_cols)
  }

  df[, unique(show_cols), drop = FALSE]
}


prepare_pca_data <- function(dobj, top_var_genes = 1000) {
  if (is.null(dobj) || !is.null(dobj$error) || is.null(dobj$plot_expr_mat) || is.null(dobj$coldata)) return(NULL)

  mat <- dobj$plot_expr_mat
  mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
  if (nrow(mat) < 3 || ncol(mat) < 3) return(NULL)

  vars <- apply(mat, 1, var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]

  genes <- names(sort(vars, decreasing = TRUE))
  genes <- head(genes, min(top_var_genes, length(genes)))

  mat2 <- mat[genes, , drop = FALSE]
  pca <- prcomp(t(mat2), center = TRUE, scale. = TRUE)

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100

  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    left_join(dobj$coldata, by = "sample")

  list(
    pca_df = pca_df,
    var_exp = var_exp,
    top_var_genes = length(genes)
  )
}


make_single_gene_group_stats <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(expression, na.rm = TRUE),
      median = median(expression, na.rm = TRUE),
      sd = sd(expression, na.rm = TRUE),
      sem = sd(expression, na.rm = TRUE) / sqrt(n()),
      min = min(expression, na.rm = TRUE),
      q25 = quantile(expression, 0.25, na.rm = TRUE),
      q75 = quantile(expression, 0.75, na.rm = TRUE),
      max = max(expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean = round(mean, 4),
      median = round(median, 4),
      sd = round(sd, 4),
      sem = round(sem, 4),
      min = round(min, 4),
      q25 = round(q25, 4),
      q75 = round(q75, 4),
      max = round(max, 4)
    )
}

make_single_gene_comparison_stats <- function(df, groupA, groupB, gene, gse_id = NA) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df2 <- df %>%
    filter(group %in% c(groupA, groupB)) %>%
    filter(!is.na(expression), !is.na(group))

  if (length(unique(df2$group)) < 2) return(NULL)

  x <- df2$expression[df2$group == groupA]
  y <- df2$expression[df2$group == groupB]

  if (length(x) < 2 || length(y) < 2) return(NULL)

  wilcox_p <- tryCatch(wilcox.test(y, x, exact = FALSE)$p.value, error = function(e) NA_real_)
  t_p <- tryCatch(t.test(y, x)$p.value, error = function(e) NA_real_)

  mean_A <- mean(x, na.rm = TRUE)
  mean_B <- mean(y, na.rm = TRUE)
  median_A <- median(x, na.rm = TRUE)
  median_B <- median(y, na.rm = TRUE)

  diff_mean <- mean_B - mean_A
  diff_median <- median_B - median_A

  # Since plotting expression is generally log2 scale, this is the log2-scale difference.
  log2FC_plot_scale <- diff_mean
  fold_change_plot_scale <- 2^log2FC_plot_scale

  data.frame(
    GSE = gse_id,
    Gene = gene,
    Group_A = groupA,
    Group_B = groupB,
    N_A = length(x),
    N_B = length(y),
    Mean_A = round(mean_A, 4),
    Mean_B = round(mean_B, 4),
    Median_A = round(median_A, 4),
    Median_B = round(median_B, 4),
    Mean_difference_B_minus_A = round(diff_mean, 4),
    Median_difference_B_minus_A = round(diff_median, 4),
    Log2FC_plot_scale_B_vs_A = round(log2FC_plot_scale, 4),
    Fold_change_plot_scale_B_vs_A = round(fold_change_plot_scale, 4),
    Wilcoxon_P = signif(wilcox_p, 4),
    T_test_P = signif(t_p, 4),
    stringsAsFactors = FALSE
  )
}

make_single_gene_combined_stats <- function(df, groupA, groupB, gene, gse_id = NA) {
  group_stats <- make_single_gene_group_stats(df)
  comp_stats <- make_single_gene_comparison_stats(df, groupA, groupB, gene, gse_id)

  list(group_stats = group_stats, comparison_stats = comp_stats)
}



# =========================
# 6.6 Correlation and enrichment helper functions
# =========================
get_plot_expr_matrix_for_modules <- function(obj, dobj = NULL) {
  if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$plot_expr_mat)) {
    return(dobj$plot_expr_mat)
  }

  if (is.null(obj) || is.null(obj$raw)) return(NULL)

  raw <- obj$raw
  det <- obj$expr_detect

  mat <- raw %>%
    column_to_rownames("Symbol") %>%
    as.matrix()

  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA

  if (is.null(det) || det$type == "raw_count" || det$type == "normalized_expression_FPKM_TPM_or_similar") {
    mat <- log2(pmax(mat, 0) + 1)
  }

  # Plotting libraries often require finite matrices.  Impute only for display,
  # using each gene's observed median rather than interpreting missing expression
  # as biological zero.  Statistical fitting above drops incomplete genes.
  if (anyNA(mat)) {
    for (i in seq_len(nrow(mat))) {
      miss <- is.na(mat[i, ])
      if (!any(miss)) next
      observed <- mat[i, !miss]
      fill <- if (length(observed) > 0) stats::median(observed, na.rm = TRUE) else NA_real_
      mat[i, miss] <- fill
    }
    mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
  }

  mat
}

find_gene_row <- function(mat, gene) {
  if (is.null(mat) || is.null(gene)) return(NA_character_)
  gene <- toupper(trimws(gene))
  hits <- rownames(mat)[toupper(rownames(mat)) == gene]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

make_gene_gene_correlation <- function(mat, meta2, gene_a, gene_b, method = "spearman") {
  if (is.null(mat) || is.null(meta2)) return(NULL)

  ga <- find_gene_row(mat, gene_a)
  gb <- find_gene_row(mat, gene_b)

  if (is.na(ga) || is.na(gb)) return(NULL)

  samples <- intersect(meta2$sample, colnames(mat))
  if (length(samples) < 3) return(NULL)

  df <- data.frame(
    sample = samples,
    Gene_A = as.numeric(mat[ga, samples]),
    Gene_B = as.numeric(mat[gb, samples]),
    stringsAsFactors = FALSE
  ) %>%
    left_join(meta2, by = "sample")

  ct <- tryCatch(
    cor.test(df$Gene_A, df$Gene_B, method = method),
    error = function(e) NULL
  )

  if (is.null(ct)) return(NULL)

  list(
    df = df,
    gene_a = ga,
    gene_b = gb,
    method = method,
    r = unname(ct$estimate),
    p = ct$p.value
  )
}

make_gene_all_correlation <- function(mat, meta2, gene, method = "spearman", top_n = 30) {
  if (is.null(mat) || is.null(meta2)) return(NULL)

  g <- find_gene_row(mat, gene)
  if (is.na(g)) return(NULL)

  samples <- intersect(meta2$sample, colnames(mat))
  if (length(samples) < 3) return(NULL)

  target <- as.numeric(mat[g, samples])

  cors <- apply(mat[, samples, drop = FALSE], 1, function(x) {
    suppressWarnings(cor(target, as.numeric(x), method = method, use = "pairwise.complete.obs"))
  })

  pvals <- apply(mat[, samples, drop = FALSE], 1, function(x) {
    tryCatch(cor.test(target, as.numeric(x), method = method)$p.value, error = function(e) NA_real_)
  })

  res <- data.frame(
    Symbol = names(cors),
    correlation = as.numeric(cors),
    pvalue = as.numeric(pvals),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(correlation), Symbol != g) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    arrange(desc(correlation))

  list(
    target_gene = g,
    full = res,
    top_positive = head(res %>% arrange(desc(correlation)), top_n),
    top_negative = head(res %>% arrange(correlation), top_n)
  )
}

parse_gene_list_input <- function(x) {
  if (is.null(x) || trimws(x) == "") return(character(0))
  x <- gsub("\n", ",", x)
  x <- gsub(";", ",", x)
  genes <- unlist(strsplit(x, ","))
  genes <- toupper(trimws(genes))
  genes <- genes[genes != ""]
  unique(genes)
}

make_gene_family_correlation <- function(mat, meta2, genes, method = "spearman") {
  if (is.null(mat) || is.null(meta2)) return(NULL)

  genes <- parse_gene_list_input(genes)
  if (length(genes) < 2) return(NULL)

  rows <- sapply(genes, function(g) find_gene_row(mat, g))
  rows <- rows[!is.na(rows)]
  rows <- unique(rows)

  if (length(rows) < 2) return(NULL)

  samples <- intersect(meta2$sample, colnames(mat))
  if (length(samples) < 3) return(NULL)

  m <- t(mat[rows, samples, drop = FALSE])
  cm <- cor(m, method = method, use = "pairwise.complete.obs")

  cm
}

get_deg_genes_for_enrichment <- function(dobj, direction = "all", padj_cutoff = 0.05, logfc_cutoff = 1) {
  if (is.null(dobj) || !is.null(dobj$error) || is.null(dobj$res_df)) return(character(0))

  res <- as.data.frame(dobj$res_df, stringsAsFactors = FALSE)

  if (!all(c("Symbol", "log2FoldChange", "padj") %in% colnames(res))) return(character(0))

  res$Symbol <- as.character(res$Symbol)
  res$log2FoldChange <- suppressWarnings(as.numeric(res$log2FoldChange))
  res$padj <- suppressWarnings(as.numeric(res$padj))

  res <- res %>%
    filter(!is.na(Symbol), Symbol != "", !is.na(padj), !is.na(log2FoldChange)) %>%
    filter(padj < padj_cutoff, abs(log2FoldChange) >= logfc_cutoff)

  if (direction == "up") {
    res <- res %>% filter(log2FoldChange > 0)
  } else if (direction == "down") {
    res <- res %>% filter(log2FoldChange < 0)
  }

  unique(res$Symbol)
}

symbol_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != ""]

  if (length(symbols) == 0) return(character(0))

  ids <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )

  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  unique(ids)
}


ensure_enrichment_packages <- function() {
  needed <- c("clusterProfiler", "enrichplot")
  missing <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }

    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }

  invisible(TRUE)
}

run_go_enrichment_safe <- function(gene_symbols, ont = "BP", p_cutoff = 0.05, q_cutoff = 0.2) {
  ensure_enrichment_packages()
  entrez <- symbol_to_entrez(gene_symbols)
  if (length(entrez) < 3) return(NULL)

  tryCatch({
    clusterProfiler::enrichGO(
      gene = entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
  }, error = function(e) NULL)
}

run_kegg_enrichment_safe <- function(gene_symbols, p_cutoff = 0.05, q_cutoff = 0.2) {
  ensure_enrichment_packages()
  entrez <- symbol_to_entrez(gene_symbols)
  if (length(entrez) < 3) return(NULL)

  tryCatch({
    clusterProfiler::enrichKEGG(
      gene = entrez,
      organism = "hsa",
      pvalueCutoff = p_cutoff,
      pAdjustMethod = "BH",
      qvalueCutoff = q_cutoff
    )
  }, error = function(e) NULL)
}

enrich_result_table <- function(x) {
  if (is.null(x)) return(NULL)
  df <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}


# =========================
# 7. Shiny UI
# =========================

# =========================
# v11.8.8 FINAL HOTFIX: rewrite single-sample COUNT reader/merger
# =========================
# Why this block exists:
# GSE228542 provides many per-sample GSM*_COUNT.txt files. Earlier generic
# parsers correctly rejected them as gene x sample matrices, but the downstream
# single-sample merger could still hit "undefined columns selected" inside a
# candidate file. This override uses only guarded column access and never subsets
# data frames by unchecked column names/indices.

safe_df_take_cols <- function(dat, cols) {
  if (is.null(dat) || is.null(cols)) return(NULL)
  cols <- unlist(cols, use.names = FALSE)
  cols <- cols[!is.na(cols)]
  if (is.numeric(cols)) {
    cols <- cols[cols >= 1 & cols <= ncol(dat)]
    if (length(cols) == 0) return(NULL)
    return(dat[, cols, drop = FALSE])
  }
  cols <- as.character(cols)
  cols <- cols[cols %in% colnames(dat)]
  if (length(cols) == 0) return(NULL)
  dat[, cols, drop = FALSE]
}

collapse_single_sample_duplicates <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  if (!all(c("Symbol", "count") %in% colnames(dat))) return(NULL)

  dat$Symbol <- clean_symbol_value(dat$Symbol)
  dat$count <- suppressWarnings(as.numeric(dat$count))
  dat <- dat[!is.na(dat$Symbol) & dat$Symbol != "" & dat$Symbol != "NA" & !is.na(dat$count), , drop = FALSE]
  if (nrow(dat) < 2) return(NULL)

  # Per-sample files are usually raw counts; summing duplicated gene IDs is safest.
  out <- stats::aggregate(count ~ Symbol, data = dat, FUN = function(z) sum(z, na.rm = TRUE))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- c("Symbol", "count")
  out
}

read_one_single_sample_expr_file_robust <- function(f, min_rows = 100) {
  bn <- basename(f)
  if (is_scRNA_like_filename(f)) return(NULL)
  if (is_too_tiny_for_expression_file(f)) return(NULL)

  dat <- tryCatch({
    data.table::fread(
      f,
      header = FALSE,
      data.table = FALSE,
      fill = TRUE,
      check.names = FALSE,
      showProgress = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(dat) || nrow(dat) < 30 || ncol(dat) < 2) return(NULL)
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)

  # Remove fully empty columns with a length-checked logical vector.
  keep_cols <- vapply(dat, function(z) any(!is.na(z) & as.character(z) != ""), logical(1))
  if (length(keep_cols) != ncol(dat)) return(NULL)
  dat <- dat[, keep_cols, drop = FALSE]
  if (nrow(dat) < 30 || ncol(dat) < 2) return(NULL)

  # For obvious COUNT/htseq files, the first column is gene ID and the second
  # numeric-like column is count. This covers GSE228542 GSM*_COUNT.txt.
  obvious_count <- grepl("(_COUNT|_count|counts?|htseq\\.results|featureCounts)", bn, ignore.case = TRUE)

  if (obvious_count) {
    gene_col <- 1
    numeric_props <- vapply(seq_len(ncol(dat)), function(i) {
      mean(!is.na(suppressWarnings(as.numeric(dat[[i]]))))
    }, numeric(1))
    numeric_props[gene_col] <- -Inf
    value_col <- which.max(numeric_props)
    if (length(value_col) == 0 || !is.finite(numeric_props[value_col]) || numeric_props[value_col] < 0.30) return(NULL)
  } else {
    # Generic fallback: first mostly non-numeric column + best numeric column.
    numeric_props <- vapply(seq_len(ncol(dat)), function(i) {
      mean(!is.na(suppressWarnings(as.numeric(dat[[i]]))))
    }, numeric(1))
    non_numeric_cols <- which(numeric_props < 0.50)
    gene_col <- if (length(non_numeric_cols) > 0) non_numeric_cols[1] else 1
    numeric_props[gene_col] <- -Inf
    value_col <- which.max(numeric_props)
    if (length(value_col) == 0 || !is.finite(numeric_props[value_col]) || numeric_props[value_col] < 0.50) return(NULL)
  }

  dat2 <- safe_df_take_cols(dat, c(gene_col, value_col))
  if (is.null(dat2) || ncol(dat2) < 2) return(NULL)
  dat2 <- dat2[, 1:2, drop = FALSE]
  colnames(dat2) <- c("Symbol", "count")
  dat2$Symbol <- as.character(dat2$Symbol)
  dat2$count <- suppressWarnings(as.numeric(dat2$count))

  # Drop header row if present and standard non-gene htseq summary rows.
  dat2 <- dat2[
    !is.na(dat2$Symbol) & dat2$Symbol != "" &
      !grepl("^__|^N_", dat2$Symbol) &
      !is.na(dat2$count),
    , drop = FALSE
  ]

  if (nrow(dat2) < min_rows) return(NULL)
  if (!is_probably_bulk_expression_pair(dat2$Symbol, dat2$count, min_rows = min_rows)) return(NULL)

  collapse_single_sample_duplicates(dat2)
}

merge_single_sample_expr_files_robust <- function(files, max_files = 5000) {
  files <- unique(files)
  files <- deduplicate_expression_files_by_sample(files)
  files <- files[!duplicated(normalizePath(files, winslash = "/", mustWork = FALSE))]
  files <- files[!vapply(files, is_too_tiny_for_expression_file, logical(1))]
  files <- head(files, max_files)

  count_list <- list()
  failed <- 0L

  for (f in files) {
    dat <- tryCatch(read_one_single_sample_expr_file_robust(f), error = function(e) {
      message("单样本文件读取失败并跳过: ", basename(f), " | ", conditionMessage(e))
      NULL
    })
    if (is.null(dat) || nrow(dat) < 2 || ncol(dat) < 2) {
      failed <- failed + 1L
      next
    }

    dat <- dat[, c("Symbol", "count"), drop = FALSE]
    sample_name <- make.names(sample_name_from_file(f), unique = FALSE)
    if (sample_name %in% names(count_list)) {
      sample_name <- make.unique(c(names(count_list), sample_name))[length(count_list) + 1]
    }
    colnames(dat)[2] <- sample_name
    count_list[[sample_name]] <- dat
  }

  message("单样本文件识别成功：", length(count_list), "；失败/跳过：", failed)
  if (length(count_list) > 0) {
    message("成功识别样本示例：", paste(head(names(count_list), 10), collapse = ", "))
  }
  if (length(count_list) < 2) return(NULL)

  message("正在合并单样本矩阵：", length(count_list), " 个样本")
  raw <- Reduce(function(x, y) merge(x, y, by = "Symbol", all = TRUE, sort = FALSE), count_list)
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  sample_cols <- setdiff(colnames(raw), "Symbol")
  for (cc in sample_cols) {
    raw[[cc]] <- suppressWarnings(as.numeric(raw[[cc]]))
    raw[[cc]][is.na(raw[[cc]])] <- 0
  }

  message("单样本矩阵合并完成：", nrow(raw), " genes x ", length(sample_cols), " samples")
  collapse_by_symbol_auto(raw, forced_method = "sum", source_hint = paste(basename(files), collapse = " "))
}


ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .well { padding: 12px; }
      .form-group { margin-bottom: 10px; }
      .irs { margin-top: -6px; }
      h4 { margin-top: 14px; }
      .tab-content { padding-top: 12px; }
      table td span { font-size: 18px; }
      .stat-card-grid {
        display: grid;
        grid-template-columns: repeat(4, minmax(140px, 1fr));
        gap: 12px;
        margin-top: 10px;
      }
      .stat-card {
        border: 1px solid #e5e5e5;
        border-radius: 10px;
        padding: 12px 14px;
        background: #ffffff;
        box-shadow: 0 1px 4px rgba(0,0,0,0.06);
        min-height: 88px;
      }
      .stat-card-title {
        font-size: 13px;
        color: #666;
        margin-bottom: 6px;
      }
      .stat-card-value {
        font-size: 22px;
        font-weight: 700;
        line-height: 1.2;
      }
      .stat-card-subtitle {
        font-size: 12px;
        color: #888;
        margin-top: 4px;
      }
      .comparison-summary-title {
        font-size: 16px;
        margin-bottom: 10px;
        padding: 8px 10px;
        background: #f7f7f7;
        border-radius: 8px;
      }
      @media (max-width: 1100px) {
        .stat-card-grid {
          grid-template-columns: repeat(2, minmax(140px, 1fr));
        }
      }

      .run-status-box {
        border: 1px solid #e5e5e5;
        border-left: 6px solid #4C78A8;
        border-radius: 10px;
        padding: 10px 14px;
        margin: 8px 0 14px 0;
        background: #ffffff;
        box-shadow: 0 1px 6px rgba(0,0,0,0.06);
      }
      .run-status-box.running {
        border-left-color: #F39C12;
        background: #fffaf0;
      }
      .run-status-box.done {
        border-left-color: #2E7D32;
        background: #f4fff7;
      }
      .run-status-box.error {
        border-left-color: #C62828;
        background: #fff5f5;
      }
      .run-status-title {
        font-weight: 700;
        font-size: 15px;
        margin-bottom: 4px;
      }
      .run-status-detail {
        color: #555;
        font-size: 13px;
      }
      .run-status-progress {
        width: 100%;
        height: 7px;
        border-radius: 999px;
        background: #eef1f5;
        overflow: hidden;
        margin-top: 8px;
      }
      .run-status-progress-fill {
        height: 100%;
        border-radius: 999px;
        background: linear-gradient(90deg, #4C78A8, #72B7B2);
        transition: width 0.25s ease;
      }
      .run-status-box.running .run-status-progress-fill {
        background: linear-gradient(90deg, #F39C12, #FFD166);
      }
      .run-status-box.done .run-status-progress-fill {
        background: linear-gradient(90deg, #2E7D32, #81C784);
      }
      .run-status-box.error .run-status-progress-fill {
        background: linear-gradient(90deg, #C62828, #EF5350);
      }

      /* =========================================================
         UI polish layer - visual only, no analysis/server changes
         ========================================================= */

      html, body {
        background: #f5f7fa;
      }

      body {
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI',
                     'Microsoft YaHei', 'PingFang SC', Arial, sans-serif;
        color: #1f2937;
      }

      .container-fluid {
        padding-left: 20px;
        padding-right: 20px;
      }

      /* App title */
      h2 {
        font-size: 25px;
        font-weight: 700;
        letter-spacing: -0.2px;
        color: #172033;
        margin-top: 16px;
        margin-bottom: 14px;
      }

      /* Main tabs */
      .nav-tabs {
        border-bottom: 1px solid #dce2e8;
        margin-bottom: 8px;
      }

      .nav-tabs > li > a {
        color: #5f6b7a;
        font-weight: 600;
        border: none !important;
        border-radius: 8px 8px 0 0;
        padding: 10px 15px;
        transition: all 0.15s ease;
      }

      .nav-tabs > li > a:hover {
        background: #eef3f8;
        color: #234a76;
      }

      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover,
      .nav-tabs > li.active > a:focus {
        color: #1d5f96;
        background: #ffffff;
        border: 1px solid #dce2e8 !important;
        border-bottom-color: #ffffff !important;
        box-shadow: 0 -2px 0 #4c86b8 inset;
      }

      /* Sidebar as one clean control card */
      .col-sm-4 > .well,
      .col-sm-3 > .well {
        border-radius: 10px;
      }

      .well {
        background: #ffffff;
        border: 1px solid #e1e6eb;
        border-radius: 10px;
        box-shadow: 0 2px 9px rgba(31, 41, 55, 0.045);
      }

      .well hr {
        border-top: 1px solid #edf0f3;
      }

      /* Input labels */
      .control-label,
      .radio > label,
      .checkbox > label {
        color: #344054;
      }

      .control-label {
        font-size: 13px;
        font-weight: 650;
        margin-bottom: 6px;
      }

      .form-control,
      .selectize-input {
        border: 1px solid #d7dee5;
        border-radius: 7px;
        box-shadow: none !important;
        min-height: 36px;
        transition: border-color 0.15s ease, box-shadow 0.15s ease;
      }

      .form-control:focus,
      .selectize-input.focus {
        border-color: #6d9fca;
        box-shadow: 0 0 0 3px rgba(76, 120, 168, 0.10) !important;
      }

      /* Buttons */
      .btn {
        border-radius: 7px;
        font-weight: 600;
        transition: all 0.15s ease;
      }

      .btn-default {
        background: #ffffff;
        border-color: #d6dde4;
        color: #344054;
      }

      .btn-default:hover {
        background: #f6f8fa;
        border-color: #bdc8d2;
      }

      .btn-success {
        background: #2f7d64;
        border-color: #2f7d64;
      }

      .btn-success:hover {
        background: #276b56;
        border-color: #276b56;
      }

      /* Headings inside tool panels */
      .tab-pane > h4,
      .well > h4 {
        font-weight: 700;
        color: #223047;
      }

      /* Make slider UI slightly calmer */
      .irs--shiny .irs-bar,
      .irs-bar {
        background: #4f86b6 !important;
        border-top-color: #4f86b6 !important;
        border-bottom-color: #4f86b6 !important;
      }

      .irs--shiny .irs-single,
      .irs-single {
        background: #4f86b6 !important;
      }

      .irs--shiny .irs-handle,
      .irs-slider {
        border-color: #8da8bf !important;
      }

      /* Tables */
      table.table {
        background: #ffffff;
        border-radius: 8px;
        overflow: hidden;
      }

      table.table > thead > tr > th {
        background: #f7f9fb;
        color: #344054;
        font-weight: 700;
        border-bottom: 1px solid #dde3e9;
      }

      /* =========================================================
         Single-gene tab
         ========================================================= */

      .single-gene-control-pane > h4 {
        margin: 2px 0 12px 2px;
        font-size: 19px;
        font-weight: 750;
        color: #172033;
      }

      .single-gene-control-pane > .well {
        padding: 14px 15px 8px 15px;
        background: transparent;
        border: none;
        box-shadow: none;
      }

      /* Each horizontal settings group behaves like a card */
      .single-gene-control-pane > .well > .row {
        background: #ffffff;
        border: 1px solid #e0e6ec;
        border-radius: 10px;
        margin: 0 0 12px 0;
        padding: 13px 10px 7px 10px;
        box-shadow: 0 2px 8px rgba(31, 41, 55, 0.04);
      }

      .single-gene-control-pane > .well > hr {
        display: none;
      }

      .single-gene-control-pane > .well > .row > [class*='col-'] {
        padding-left: 13px;
        padding-right: 13px;
      }

      .single-gene-control-pane > .well > .row > [class*='col-']:first-child {
        border-right: 1px solid #eef1f4;
      }

      .single-gene-control-pane h4 {
        font-size: 15px;
        font-weight: 750;
        color: #315d84;
        margin-top: 2px;
        margin-bottom: 10px;
        padding-bottom: 7px;
        border-bottom: 1px solid #edf1f4;
      }

      .single-gene-control-pane .form-group {
        margin-bottom: 12px;
      }

      .single-gene-control-pane .checkbox,
      .single-gene-control-pane .radio {
        margin-top: 7px;
        margin-bottom: 7px;
      }

      /* Checkbox/radio spacing */
      .single-gene-control-pane input[type='checkbox'],
      .single-gene-control-pane input[type='radio'] {
        margin-top: 3px;
      }

      /* Colour input blocks: more compact */
      .single-gene-control-pane .colourpicker-input-container {
        margin-bottom: 9px;
      }

      /* Sticky preview becomes a proper figure card */
      .single-gene-sticky-preview {
        top: 14px !important;
        border: 1px solid #dce3e9 !important;
        border-radius: 12px !important;
        padding: 13px 13px 15px 13px !important;
        box-shadow: 0 5px 18px rgba(31, 41, 55, 0.085) !important;
        background: #ffffff !important;
      }

      .single-gene-sticky-preview h4 {
        font-size: 17px;
        font-weight: 750;
        color: #223047;
        padding: 0 0 10px 1px;
        margin: 0 0 9px 0 !important;
        border-bottom: 1px solid #edf1f4;
      }

      /* Top quick switches on the single-gene page */
      .tab-pane .single-gene-layout-row {
        margin-top: 4px;
      }

      /* Export and stats cards below the plot */
      .single-gene-layout-row + .well,
      .single-gene-layout-row + .well + .well {
        border-radius: 10px;
        border: 1px solid #e0e6ec;
        box-shadow: 0 2px 8px rgba(31, 41, 55, 0.04);
      }

      /* Download buttons align more like a toolbar */
      .single-gene-layout-row + .well .btn,
      .single-gene-layout-row + .well + .well .btn {
        min-width: 145px;
      }

      /* Small help text */
      .help-block {
        color: #7b8794;
        font-size: 12px;
        line-height: 1.45;
      }

      /* Soft scrollbar for modern browsers */
      ::-webkit-scrollbar {
        width: 10px;
        height: 10px;
      }

      ::-webkit-scrollbar-track {
        background: #f3f5f7;
      }

      ::-webkit-scrollbar-thumb {
        background: #c7d0d9;
        border-radius: 999px;
        border: 2px solid #f3f5f7;
      }

      ::-webkit-scrollbar-thumb:hover {
        background: #aebbc6;
      }

      @media (max-width: 1200px) {
        .single-gene-control-pane > .well > .row > [class*='col-']:first-child {
          border-right: none;
        }
      }

      /* Final plotting-workbench polish */
      .plot-toolbar {
        background: linear-gradient(180deg, #ffffff 0%, #f8fafc 100%);
        border: 1px solid #dde5ec;
        border-radius: 11px;
        padding: 12px 12px 4px 12px;
        margin-bottom: 12px;
        box-shadow: 0 2px 8px rgba(31, 41, 55, 0.045);
      }

      .plot-toolbar .form-group {
        margin-bottom: 8px;
      }

      .plot-toolbar .btn {
        width: 100%;
        min-height: 36px;
        margin-top: 25px;
      }

      .plot-status-bar {
        display: flex;
        flex-wrap: wrap;
        gap: 7px;
        margin: 0 0 10px 0;
      }

      .plot-status-chip {
        display: inline-block;
        padding: 5px 9px;
        border-radius: 999px;
        border: 1px solid #dce5ec;
        background: #f7fafc;
        color: #516173;
        font-size: 11.5px;
        font-weight: 650;
      }

      .single-gene-control-pane .advanced-note {
        color: #7a8794;
        font-size: 11.5px;
        margin-top: -2px;
        margin-bottom: 9px;
      }

      .single-gene-control-pane .btn-primary {
        background: #3f78a8;
        border-color: #3f78a8;
      }

      .single-gene-control-pane .btn-primary:hover {
        background: #356b96;
        border-color: #356b96;
      }
    "))
  ),
  titlePanel("GEO Gene Expression & Differential Analysis Explorer"),
  uiOutput("run_status_ui"),
  sidebarLayout(
    sidebarPanel(
      textInput("gse", "GSE ID", value = "GSE171741"),
      checkboxInput("clean_cache", "加载前自动删除当前GSE下载文件", value = FALSE),
      checkboxInput("use_parsed_cache", "优先读取已解析缓存（推荐，更快）", value = TRUE),
      checkboxInput("rebuild_parsed_cache", "强制重新解析当前GSE", value = FALSE),
      actionButton("run", "1. Load / Reload GSE"),
      hr(),
      textInput("gene", "Gene symbol（换基因不用重新Run）", value = "CPNE2"),
      hr(),
      selectInput("filter_col", "Filter column（先筛选，可选）", choices = NULL),
      selectInput("filter_value", "Filter value", choices = NULL),
      hr(),
      selectInput("group_col", "Compare column（比较分组）", choices = NULL),
      selectInput("groupA", "Seed Group A / Reference seed（初始选择；合并开启后不一定是最终分析组）", choices = NULL),
      selectInput("groupB", "Seed Group B / Case seed（初始选择；合并开启后不一定是最终分析组）", choices = NULL),

      hr(),
      checkboxInput("enable_group_merge", "Enable manual group merge（定义最终分析组）", value = FALSE),
      conditionalPanel(
        condition = "input.enable_group_merge == true",
        textInput("mergeA_name", "Final Group A name（used in DEG/plots）", value = "Control"),
        selectizeInput("mergeA_groups", "Metadata levels included in Final Group A", choices = NULL, multiple = TRUE),
        textInput("mergeB_name", "Final Group B name（used in DEG/plots）", value = "Sepsis"),
        selectizeInput("mergeB_groups", "Metadata levels included in Final Group B", choices = NULL, multiple = TRUE),
        helpText("开启合并后，DEG、单基因图、火山图、热图使用 Final Group A/B；上面的 Seed Group A/B 只用于快速初始化选择，不代表最终分析组。原始 metadata 不会被覆盖。")
      ),

      hr(),
      selectInput(
        "expr_source",
        "Expression source（表达矩阵来源；Auto默认优先raw，可手动切换RPM/TPM/FPKM）",
        choices = c("Auto recommended after Load" = "auto"),
        selected = "auto"
      ),
      helpText("Auto保持原始底层逻辑：优先使用Load阶段自动推荐的数据源；如果某个GSE的raw/RPM样本不一致，可在这里手动切换矩阵。选中的矩阵会同步用于DEG、单基因图、热图、PCA等。"),

      radioButtons(
        "method_mode",
        "Differential analysis method",
        choices = c(
          "Auto: raw count -> DESeq2; normalized/log -> limma" = "auto",
          "Force DESeq2 (only for raw count)" = "DESeq2",
          "Force limma on best matched normalized/RPM matrix" = "limma",
          "Force limma-voom on raw counts (keeps count-matrix samples)" = "limma_voom"
        ),
        selected = "auto"
      ),

      radioButtons(
        "sample_match_mode",
        "Sample matching mode",
        choices = c(
          "Auto strict: cleaned/GSM/sample-title matching only (recommended)" = "auto_robust",
          "Cleaned/GSM sample-name matching only" = "clean",
          "Strict exact sample-name matching only" = "strict",
          "Force metadata order if sample count is identical (dangerous; only if verified)" = "order"
        ),
        selected = "auto_robust"
      ),

      actionButton("run_deg", "2. Run DEG after manual group selection", class = "btn-success"),
      br(), br(),

      hr(),
      helpText("Load只读取数据并自动推荐 Seed Groups；如果启用Manual Merge，请确认 Final Groups 后点击 Run DEG。单基因图显示设置在“单基因图”标签页中。"),
      hr(),
      downloadButton("downloadDEG", "Download all differential results"),
      br(), br(),
      downloadButton("downloadSigDEG", "Download significant DEGs"),
      br(), br(),
      hr(),
      h4("Meta / downstream export"),
      downloadButton("downloadMergedExpressionMatrix", "导出合并表达矩阵 CSV"),
      br(), br(),
      downloadButton("downloadExpressionAvailableMetadata", "导出表达矩阵可匹配分组 CSV"),
      br(), br(),
      downloadButton("downloadPlotScaleExpressionMatrix", "导出log/plot-scale矩阵 CSV"),
      br(), br(),
      downloadButton("downloadMetaReadyDegMatrix", "导出当前Run DEG实际Meta矩阵 CSV（推荐）"),
      br(), br(),
      helpText("推荐做Meta/WGCNA时使用“当前Run DEG实际Meta矩阵”：它直接来自最后一次Run DEG的实际样本、实际分组和plot-scale表达矩阵，不会重新匹配metadata，避免GROUP与表达列错位。"),
      br(),
      helpText("合并表达矩阵 = 当前 Expression source 的 Symbol × samples 矩阵；raw count 保持原始count，TPM/FPKM保持原值。log/plot-scale矩阵用于相关性、热图、PCA和单基因展示。"),
      hr(),
      helpText("图像下载按钮已整合到对应作图模块中。")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Dashboard",
          h4("Dataset Dashboard"),
          fluidRow(
            column(3, wellPanel(h4("Samples"), h3(textOutput("dash_samples", inline = TRUE)))),
            column(3, wellPanel(h4("Genes"), h3(textOutput("dash_genes", inline = TRUE)))),
            column(3, wellPanel(h4("Method"), h3(textOutput("dash_method", inline = TRUE)))),
            column(3, wellPanel(h4("Sample match"), h3(textOutput("dash_match", inline = TRUE))))
          ),
          fluidRow(
            column(3, wellPanel(h4("Final Group A"), h3(textOutput("dash_groupA", inline = TRUE)))),
            column(3, wellPanel(h4("Final Group B"), h3(textOutput("dash_groupB", inline = TRUE)))),
            column(3, wellPanel(h4("Up DEG"), h3(textOutput("dash_up", inline = TRUE)))),
            column(3, wellPanel(h4("Down DEG"), h3(textOutput("dash_down", inline = TRUE))))
          ),
          fluidRow(
            column(6, wellPanel(h4("Seed Group A / B"), textOutput("dash_seed_groups"))),
            column(6, wellPanel(h4("Group sample status"), tableOutput("dash_group_sample_status")))
          ),

          hr(),
          h4("DEG Summary"),
          tableOutput("deg_summary_table"),

          hr(),
          h4("Search DEG"),
          fluidRow(
            column(4, textInput("deg_search", "Search gene symbol / family", value = "CPNE")),
            column(4, sliderInput("deg_table_padj", "Padj cutoff", min = 0.001, max = 1, value = 0.05, step = 0.001)),
            column(4, sliderInput("deg_table_logfc", "|log2FC| cutoff", min = 0, max = 3, value = 1, step = 0.1))
          ),
          tableOutput("deg_search_table"),

          hr(),
          fluidRow(
            column(
              6,
              h4("Top 20 Upregulated"),
              tableOutput("top_up_table")
            ),
            column(
              6,
              h4("Top 20 Downregulated"),
              tableOutput("top_down_table")
            )
          ),

          hr(),
          h4("Selected gene differential result"),
          tableOutput("gene_deg_result"),

          h4("Group summary"),
          tableOutput("summary"),

          h4("Single-gene quick comparison"),
          tableOutput("compare"),

          hr(),
          h4("Message / Log"),
          uiOutput("message_ui")
        ),

        tabPanel(
          "单基因图",

          tags$style(HTML("
            .single-gene-control-pane .shiny-input-container { width: 100%; }
            /* Make the preview column as tall as the control column.
               Without this, CSS sticky is constrained by a short Bootstrap column
               and appears not to follow the page scroll. */
            .single-gene-layout-row {
              display: -webkit-flex;
              display: flex;
              -webkit-align-items: stretch;
              align-items: stretch;
              overflow: visible !important;
            }
            .single-gene-layout-row > [class*='col-'] {
              overflow: visible !important;
            }
            .single-gene-sticky-preview {
              position: -webkit-sticky;
              position: sticky;
              top: 12px;
              z-index: 20;
              align-self: flex-start;
              background: #ffffff;
              border: 1px solid #dddddd;
              border-radius: 6px;
              padding: 10px;
              box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            }
            .single-gene-sticky-preview h4 {
              margin-top: 0;
              margin-bottom: 8px;
            }
          ")),

          fluidRow(
            column(4, checkboxInput("show_all", "Show all groups in single-gene plot", value = TRUE)),
            column(4, checkboxInput("show_plot_p", "Show Wilcoxon P on single-gene plot", value = FALSE))
          ),

          div(
            class = "row single-gene-layout-row",
            column(
              7,
              div(
                class = "single-gene-control-pane",
                h4("单基因图调节面板"),
                div(
                  class = "plot-toolbar",
                  fluidRow(
                    column(
                      4,
                      selectInput(
                        "plot_preset",
                        "快速绘图风格",
                        choices = c(
                          "Prism / Publication" = "prism",
                          "Nature / Clean" = "nature",
                          "Minimal" = "minimal",
                          "当前自定义" = "custom"
                        ),
                        selected = "prism"
                      )
                    ),
                    column(
                      4,
                      actionButton("apply_plot_preset", "应用风格", class = "btn-primary")
                    ),
                    column(
                      4,
                      actionButton("reset_plot_style", "恢复默认")
                    )
                  )
                ),
                wellPanel(
                  fluidRow(
                    column(
                      6,
                      h4("图形类型"),
                      radioButtons(
                        "plot_style",
                        "图形样式",
                        choices = c(
                          "仅箱线图" = "box",
                          "仅小提琴图" = "violin",
                          "小提琴图 + 箱线图" = "violin_box"
                        ),
                        selected = "violin_box"
                      ),
                      checkboxInput("show_legend", "显示图例", value = TRUE),
                      checkboxInput("show_major_grid", "显示主网格线", value = TRUE),
                      checkboxInput("show_minor_grid", "显示次网格线", value = FALSE),
                      checkboxInput("show_plot_title", "显示标题", value = TRUE),
                      checkboxInput("show_plot_subtitle", "显示副标题", value = TRUE)
                    ),
                    column(
                      6,
                      h4("颜色"),
                      colourpicker::colourInput("groupA_color", "A组颜色", value = "#E64B35"),
                      colourpicker::colourInput("groupB_color", "B组颜色", value = "#0072B5"),
                      colourpicker::colourInput("other_color", "其他组颜色", value = "#00A087"),
                      sliderInput("fill_alpha", "箱体/小提琴填充透明度", min = 0, max = 1, value = 0.22, step = 0.05)
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      h4("点和组间距"),
                      sliderInput("point_size", "散点大小", min = 1, max = 8, value = 3, step = 0.2),
                      sliderInput("point_alpha", "散点透明度", min = 0.1, max = 1, value = 0.9, step = 0.05),
                      sliderInput("jitter_width", "散点左右抖动范围", min = 0.01, max = 0.5, value = 0.16, step = 0.01),
                      sliderInput("x_spacing", "两组之间的视觉间距", min = 0.6, max = 2.5, value = 1.0, step = 0.1),
                      selectInput(
                        "point_shape",
                        "散点形状",
                        choices = c("实心圆" = 16, "空心圆" = 1, "方形" = 15, "菱形" = 18, "三角形" = 17),
                        selected = 16
                      ),
                      sliderInput("point_stroke", "散点边线粗细", min = 0, max = 2, value = 0.5, step = 0.1)
                    ),
                    column(
                      6,
                      h4("箱线图/小提琴"),
                      sliderInput("box_line_width", "箱线图框线粗细", min = 0.2, max = 3, value = 1.1, step = 0.1),
                      sliderInput("box_width", "箱线图宽度", min = 0.15, max = 0.9, value = 0.42, step = 0.05),
                      sliderInput("violin_line_width", "小提琴图边框粗细", min = 0.2, max = 3, value = 0.9, step = 0.1),
                      sliderInput("violin_width", "小提琴图宽度", min = 0.2, max = 1.2, value = 0.85, step = 0.05)
                    )
                  ),

                  hr(),

                  fluidRow(
                    column(
                      6,
                      h4("制图比例"),
                      sliderInput("plot_height_px", "图像显示高度(px)", min = 400, max = 1000, value = 650, step = 50),
                      sliderInput("export_width", "导出宽度(inch)", min = 4, max = 16, value = 10, step = 0.5),
                      sliderInput("export_height", "导出高度(inch)", min = 3, max = 12, value = 7, step = 0.5)
                    ),
                    column(
                      6,
                      h4("字体"),
                      sliderInput("base_font_size", "基础字体大小", min = 8, max = 24, value = 13, step = 1),
                      sliderInput("axis_text_angle", "X轴文字角度", min = 0, max = 90, value = 45, step = 5),
                      sliderInput("title_font_size", "标题字体大小", min = 10, max = 30, value = 16, step = 1),
                      sliderInput("subtitle_font_size", "副标题字体大小", min = 7, max = 24, value = 12, step = 1),
                      selectInput(
                        "plot_font_family",
                        "字体",
                        choices = c("系统默认" = "", "Arial" = "Arial", "Helvetica" = "Helvetica",
                                    "Times New Roman" = "Times New Roman", "Georgia" = "Georgia"),
                        selected = ""
                      ),
                      selectInput(
                        "title_hjust",
                        "标题对齐",
                        choices = c("左" = 0, "居中" = 0.5, "右" = 1),
                        selected = 0
                      )
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      h4("边距"),
                      sliderInput("plot_margin_top", "上边距", min = 0, max = 40, value = 10, step = 2),
                      sliderInput("plot_margin_right", "右边距", min = 0, max = 40, value = 10, step = 2),
                      sliderInput("plot_margin_bottom", "下边距", min = 0, max = 40, value = 10, step = 2),
                      sliderInput("plot_margin_left", "左边距", min = 0, max = 40, value = 10, step = 2)
                    ),
                    column(
                      6,
                      h4("坐标轴"),
                      checkboxInput("free_y_expand", "Y轴自动留白", value = TRUE),
                      sliderInput("y_expand_ratio", "Y轴留白比例", min = 0.01, max = 0.30, value = 0.08, step = 0.01),
                      checkboxInput("manual_y_limits", "手动设置Y轴范围", value = FALSE),
                      conditionalPanel(
                        condition = "input.manual_y_limits == true",
                        fluidRow(
                          column(6, numericInput("y_axis_min", "Y最小值", value = 0)),
                          column(6, numericInput("y_axis_max", "Y最大值", value = 10))
                        )
                      ),
                      numericInput("y_decimal_digits", "Y轴数字小数位", value = 1, min = 0, max = 6, step = 1),
                      checkboxInput("hide_both_axis_titles", "隐藏X/Y两个轴标题文字", value = FALSE),
                      checkboxInput("hide_x_title", "隐藏X轴标题", value = FALSE),
                      checkboxInput("hide_y_title", "隐藏Y轴标题", value = FALSE),
                      checkboxInput("hide_x_tick_labels", "隐藏X轴刻度文字（组名）", value = FALSE),
                      checkboxInput("hide_y_tick_labels", "隐藏Y轴刻度文字（表达量数字）", value = FALSE)
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      h4("轴线粗细"),
                      sliderInput("x_axis_line_width", "X轴线粗细", min = 0, max = 3, value = 0.8, step = 0.1),
                      sliderInput("y_axis_line_width", "Y轴线粗细", min = 0, max = 3, value = 0.8, step = 0.1),
                      sliderInput("axis_tick_width", "刻度线粗细", min = 0, max = 3, value = 0.8, step = 0.1),
                      sliderInput("axis_tick_length", "刻度线长度(pt)", min = 0, max = 12, value = 4, step = 0.5)
                    ),
                    column(
                      6,
                      h4("刻度间距"),
                      textInput("y_major_interval", "Y轴主刻度间距（留空=自动）", value = ""),
                      textInput("y_minor_interval", "Y轴次刻度间距（留空=自动/不指定）", value = ""),
                      helpText("例如主刻度填 1 或 0.5；必须为正数。")
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      h4("坐标轴字体"),
                      sliderInput("axis_text_size", "刻度文字大小", min = 6, max = 24, value = 11, step = 1),
                      sliderInput("axis_title_size", "轴标题大小", min = 6, max = 26, value = 13, step = 1),
                      checkboxInput("bold_axis_text", "刻度文字加粗", value = FALSE),
                      checkboxInput("bold_axis_title", "轴标题加粗", value = FALSE)
                    ),
                    column(
                      6,
                      h4("刻度显示"),
                      checkboxInput("show_x_ticks", "显示X轴刻度线", value = TRUE),
                      checkboxInput("show_y_ticks", "显示Y轴刻度线", value = TRUE),
                      checkboxInput("show_x_axis_line", "显示X轴线", value = TRUE),
                      checkboxInput("show_y_axis_line", "显示Y轴线", value = TRUE)
                    )
                  ),

                  fluidRow(
                    column(
                      6,
                      h4("标题与图例"),
                      textInput("custom_plot_title", "自定义标题（留空=自动）", value = ""),
                      textInput("custom_plot_subtitle", "自定义副标题（留空=自动）", value = ""),
                      textInput("custom_x_title", "自定义X轴标题（留空=Group）", value = ""),
                      textInput("custom_y_title", "自定义Y轴标题（留空=自动表达量标题）", value = ""),
                      selectInput(
                        "legend_position",
                        "图例位置",
                        choices = c("右侧" = "right", "左侧" = "left", "顶部" = "top", "底部" = "bottom"),
                        selected = "right"
                      ),
                      sliderInput("legend_text_size", "图例文字大小", min = 6, max = 22, value = 11, step = 1),
                      sliderInput("legend_title_size", "图例标题大小", min = 6, max = 24, value = 12, step = 1)
                    ),
                    column(
                      6,
                      h4("统计标注与均值"),
                      selectInput(
                        "plot_p_label_type",
                        "P值显示格式",
                        choices = c("精确P值" = "p.format", "显著性星号" = "p.signif"),
                        selected = "p.format"
                      ),
                      sliderInput("plot_p_text_size", "P值文字大小", min = 2, max = 10, value = 4.5, step = 0.5),
                      textInput("plot_p_y", "P值高度（留空=自动）", value = ""),
                      checkboxInput("show_mean_marker", "显示各组均值标记", value = FALSE),
                      sliderInput("mean_marker_size", "均值标记大小", min = 1, max = 10, value = 4, step = 0.5),
                      selectInput(
                        "mean_marker_shape",
                        "均值标记形状",
                        choices = c("菱形" = 18, "实心圆" = 16, "方形" = 15, "三角形" = 17),
                        selected = 18
                      ),
                      colourpicker::colourInput("mean_marker_color", "均值标记颜色", value = "#000000")
                    )
                  )
                )
              )
            ),
            column(
              5,
              div(
                class = "single-gene-sticky-preview",
                h4("实时图预览"),
                uiOutput("plot_status_ui"),
                uiOutput("plot_ui")
              )
            )
          ),

          wellPanel(
            h4("图像导出"),
            fluidRow(
              column(3, downloadButton("downloadPlotPNG", "下载单基因图 PNG")),
              column(3, downloadButton("downloadPlotPDF", "下载单基因图 PDF")),
              column(3, downloadButton("downloadPlotSVG", "下载单基因图 SVG"))
            )
          ),

          wellPanel(
            h4("单基因统计信息"),

            fluidRow(
              column(
                12,
                h4("各组描述性统计"),
                div(
                  style = "overflow-x:auto; width:100%;",
                  tableOutput("single_gene_group_stats")
                )
              )
            ),

            hr(),

            fluidRow(
              column(
                12,
                h4("两组比较统计"),
                uiOutput("single_gene_comparison_cards")
              )
            ),

            br(),
            downloadButton("downloadSingleGeneStats", "导出单基因统计 CSV")
          )
        ),

        tabPanel(
          "火山图",
          h4("火山图调节面板"),
          wellPanel(
            fluidRow(
              column(
                3,
                h4("阈值"),
                sliderInput("volcano_logfc_cutoff", "log2FC阈值", min = 0, max = 3, value = 1, step = 0.1),
                sliderInput("volcano_padj_cutoff", "Padj阈值", min = 0.001, max = 0.2, value = 0.05, step = 0.001),
                radioButtons(
                  "volcano_p_type",
                  "Y轴使用",
                  choices = c("Adjusted P value" = "padj", "Raw P value" = "pvalue"),
                  selected = "padj"
                )
              ),
              column(
                3,
                h4("颜色"),
                colourpicker::colourInput("volcano_up_color", "上调颜色", value = "#E64B35"),
                colourpicker::colourInput("volcano_down_color", "下调颜色", value = "#0072B5"),
                colourpicker::colourInput("volcano_ns_color", "不显著颜色", value = "#BDBDBD"),
                colourpicker::colourInput("volcano_target_color", "目标基因颜色", value = "#000000")
              ),
              column(
                3,
                h4("点和标签"),
                sliderInput("volcano_point_size", "点大小", min = 0.3, max = 4, value = 1.2, step = 0.1),
                sliderInput("volcano_point_alpha", "点透明度", min = 0.1, max = 1, value = 0.7, step = 0.05),
                sliderInput("volcano_target_size", "目标基因点大小", min = 1, max = 8, value = 3.5, step = 0.2),
                checkboxInput("volcano_show_target_label", "显示目标基因标签", value = TRUE)
              ),
              column(
                3,
                h4("显示"),
                checkboxInput("volcano_show_cutoff_lines", "显示阈值虚线", value = TRUE),
                checkboxInput("volcano_show_legend", "显示图例", value = TRUE),
                sliderInput("volcano_base_font_size", "字体大小", min = 8, max = 24, value = 13, step = 1),
                sliderInput("volcano_height_px", "显示高度(px)", min = 400, max = 1000, value = 650, step = 50)
              )
            ),
            hr(),
            fluidRow(
              column(3, sliderInput("volcano_export_width", "导出宽度(inch)", min = 4, max = 16, value = 10, step = 0.5)),
              column(3, sliderInput("volcano_export_height", "导出高度(inch)", min = 3, max = 12, value = 7, step = 0.5)),
              column(3, selectInput("volcano_export_dpi", "导出DPI", choices = c(300, 600, 1200), selected = 600)),
              column(3, br(), checkboxInput("volcano_export_label", "导出时显示目标标签", value = TRUE))
            )
          ),
          wellPanel(
            h4("火山图导出"),
            fluidRow(
              column(3, downloadButton("downloadVolcanoPNG", "下载火山图 PNG")),
              column(3, downloadButton("downloadVolcanoPDF", "下载火山图 PDF")),
              column(3, downloadButton("downloadVolcanoSVG", "下载火山图 SVG"))
            ),
            br(),
            fluidRow(
              column(3, downloadButton("downloadUpGenes", "导出上调基因 CSV")),
              column(3, downloadButton("downloadDownGenes", "导出下调基因 CSV")),
              column(3, downloadButton("downloadVolcanoSigGenes", "导出显著基因 CSV"))
            )
          ),
          uiOutput("volcano_ui")
        ),

        tabPanel(
          "聚类热图",
          h4("聚类热图调节面板"),
          wellPanel(
            fluidRow(
              column(
                3,
                h4("基因筛选"),
                numericInput("heatmap_top_n", "Top基因数量", value = 50, min = 5, max = 500, step = 5),
                sliderInput("heatmap_logfc_cutoff", "最小|log2FC|", min = 0, max = 3, value = 0, step = 0.1),
                sliderInput("heatmap_padj_cutoff", "最大Padj", min = 0.001, max = 1, value = 0.05, step = 0.001),
                radioButtons(
                  "heatmap_rank_by",
                  "排序方式",
                  choices = c("Padj最小" = "padj", "|log2FC|最大" = "abslogfc", "P值最小" = "pvalue"),
                  selected = "padj"
                )
              ),
              column(
                3,
                h4("聚类"),
                checkboxInput("heatmap_cluster_rows", "聚类基因", value = TRUE),
                checkboxInput("heatmap_cluster_cols", "聚类样本", value = TRUE),
                checkboxInput("heatmap_show_rownames", "显示基因名", value = TRUE),
                checkboxInput("heatmap_show_colnames", "显示样本名", value = FALSE),
                radioButtons(
                  "heatmap_scale",
                  "数据标准化",
                  choices = c("按基因Z-score" = "row", "不标准化" = "none"),
                  selected = "row"
                )
              ),
              column(
                3,
                h4("颜色"),
                selectInput(
                  "heatmap_palette",
                  "热图配色",
                  choices = c("蓝-白-红" = "blue_white_red", "绿-黑-红" = "green_black_red", "紫-白-橙" = "purple_white_orange"),
                  selected = "blue_white_red"
                ),
                colourpicker::colourInput("heatmap_groupA_color", "A组注释颜色", value = "#E64B35"),
                colourpicker::colourInput("heatmap_groupB_color", "B组注释颜色", value = "#0072B5"),
                sliderInput("heatmap_fontsize", "整体字体大小", min = 5, max = 18, value = 9, step = 1),
                sliderInput("heatmap_fontsize_row", "基因名字体大小", min = 3, max = 16, value = 7, step = 1)
              ),
              column(
                3,
                h4("尺寸"),
                sliderInput("heatmap_height_px", "显示高度(px)", min = 400, max = 1200, value = 700, step = 50),
                sliderInput("heatmap_cellwidth", "单元格宽度", min = 4, max = 30, value = 10, step = 1),
                sliderInput("heatmap_cellheight", "单元格高度", min = 4, max = 30, value = 10, step = 1),
                sliderInput("heatmap_border_width", "网格边框宽度", min = 0, max = 1, value = 0, step = 0.1)
              )
            ),
            hr(),
            fluidRow(
              column(3, sliderInput("heatmap_export_width", "导出宽度(inch)", min = 4, max = 20, value = 10, step = 0.5)),
              column(3, sliderInput("heatmap_export_height", "导出高度(inch)", min = 4, max = 20, value = 9, step = 0.5)),
              column(3, selectInput("heatmap_export_dpi", "导出DPI", choices = c(300, 600, 1200), selected = 600))
            )
          ),
          wellPanel(
            h4("热图导出"),
            fluidRow(
              column(3, downloadButton("downloadHeatmapPNG", "下载热图 PNG")),
              column(3, downloadButton("downloadHeatmapPDF", "下载热图 PDF")),
              column(3, downloadButton("downloadHeatmapSVG", "下载热图 SVG")),
              column(3, downloadButton("downloadHeatmapGenes", "导出热图基因 CSV"))
            )
          ),
          uiOutput("heatmap_ui")
        ),

        tabPanel(
          "PCA",
          h4("PCA调节面板"),
          wellPanel(
            fluidRow(
              column(
                3,
                numericInput("pca_top_var_genes", "用于PCA的高变基因数量", value = 1000, min = 100, max = 10000, step = 100),
                sliderInput("pca_point_size", "点大小", min = 1, max = 8, value = 3.5, step = 0.2),
                sliderInput("pca_point_alpha", "点透明度", min = 0.1, max = 1, value = 0.85, step = 0.05)
              ),
              column(
                3,
                colourpicker::colourInput("pca_groupA_color", "A组颜色", value = "#E64B35"),
                colourpicker::colourInput("pca_groupB_color", "B组颜色", value = "#0072B5"),
                colourpicker::colourInput("pca_other_color", "其他组颜色", value = "#00A087")
              ),
              column(
                3,
                checkboxInput("pca_show_labels", "显示样本名", value = FALSE),
                checkboxInput("pca_show_ellipse", "显示分组椭圆", value = TRUE),
                checkboxInput("pca_show_grid", "显示网格线", value = TRUE)
              ),
              column(
                3,
                sliderInput("pca_base_font_size", "字体大小", min = 8, max = 24, value = 13, step = 1),
                sliderInput("pca_height_px", "显示高度(px)", min = 400, max = 1000, value = 650, step = 50),
                sliderInput("pca_export_width", "导出宽度(inch)", min = 4, max = 16, value = 8, step = 0.5),
                sliderInput("pca_export_height", "导出高度(inch)", min = 3, max = 12, value = 6, step = 0.5)
              )
            )
          ),
          wellPanel(
            h4("PCA导出"),
            fluidRow(
              column(3, downloadButton("downloadPCAPNG", "下载PCA PNG")),
              column(3, downloadButton("downloadPCAPDF", "下载PCA PDF")),
              column(3, downloadButton("downloadPCASVG", "下载PCA SVG"))
            )
          ),
          uiOutput("pca_ui")
        ),

        tabPanel(
          "Correlation",
          h4("相关性分析模块"),
          tabsetPanel(
            tabPanel(
              "Gene vs Gene",
              wellPanel(
                fluidRow(
                  column(3, textInput("corr_gene_a", "Gene A", value = "CPNE2")),
                  column(3, textInput("corr_gene_b", "Gene B", value = "IL6")),
                  column(3, selectInput("corr_method", "相关性方法", choices = c("Spearman" = "spearman", "Pearson" = "pearson"), selected = "spearman")),
                  column(3, checkboxInput("corr_use_selected_groups", "只使用当前A/B两组", value = TRUE))
                ),
                fluidRow(
                  column(3, sliderInput("corr_point_size", "点大小", min = 1, max = 8, value = 3, step = 0.2)),
                  column(3, sliderInput("corr_point_alpha", "点透明度", min = 0.1, max = 1, value = 0.85, step = 0.05)),
                  column(3, checkboxInput("corr_show_lm", "显示线性拟合线", value = TRUE)),
                  column(3, checkboxInput("corr_show_labels", "显示样本名", value = FALSE))
                )
              ),
              h4("相关性统计"),
              tableOutput("corr_gene_gene_stats"),
              plotOutput("corr_gene_gene_plot", height = "600px")
            ),

            tabPanel(
              "Gene vs All Genes",
              wellPanel(
                fluidRow(
                  column(3, textInput("corr_all_gene", "Target gene", value = "CPNE2")),
                  column(3, numericInput("corr_all_top_n", "Top数量", value = 30, min = 5, max = 200, step = 5)),
                  column(3, selectInput("corr_all_method", "相关性方法", choices = c("Spearman" = "spearman", "Pearson" = "pearson"), selected = "spearman")),
                  column(3, checkboxInput("corr_all_current_groups", "只使用当前A/B两组", value = TRUE))
                )
              ),
              fluidRow(
                column(6, h4("Top positive correlated genes"), tableOutput("corr_top_positive")),
                column(6, h4("Top negative correlated genes"), tableOutput("corr_top_negative"))
              ),
              downloadButton("downloadCorrAll", "导出全基因相关性 CSV")
            ),

            tabPanel(
              "Gene Family Matrix",
              wellPanel(
                fluidRow(
                  column(6, textAreaInput("corr_family_genes", "输入基因列表（逗号或换行分隔）", value = "CPNE1, CPNE2, CPNE3, CPNE4, CPNE5, CPNE6, CPNE7, CPNE8, CPNE9", rows = 4)),
                  column(3, selectInput("corr_family_method", "相关性方法", choices = c("Spearman" = "spearman", "Pearson" = "pearson"), selected = "spearman")),
                  column(3, checkboxInput("corr_family_cluster", "聚类基因", value = TRUE))
                )
              ),
              plotOutput("corr_family_heatmap", height = "700px"),
              downloadButton("downloadCorrFamily", "导出相关矩阵 CSV")
            )
          )
        ),

        tabPanel(
          "GO/KEGG",
          h4("功能富集分析模块"),
          wellPanel(
            fluidRow(
              column(3, selectInput("enrich_direction", "基因集", choices = c("All DEGs" = "all", "Up genes" = "up", "Down genes" = "down"), selected = "all")),
              column(3, sliderInput("enrich_padj_cutoff", "DEG Padj阈值", min = 0.001, max = 1, value = 0.05, step = 0.001)),
              column(3, sliderInput("enrich_logfc_cutoff", "DEG |log2FC|阈值", min = 0, max = 3, value = 1, step = 0.1)),
              column(3, numericInput("enrich_show_n", "显示Top terms", value = 15, min = 5, max = 50, step = 5))
            ),
            fluidRow(
              column(3, sliderInput("enrich_p_cutoff", "富集pvalueCutoff", min = 0.001, max = 1, value = 0.05, step = 0.001)),
              column(3, sliderInput("enrich_q_cutoff", "富集qvalueCutoff", min = 0.001, max = 1, value = 0.2, step = 0.001)),
              column(3, selectInput("enrich_go_ont", "GO类型", choices = c("BP", "CC", "MF"), selected = "BP")),
              column(3, actionButton("run_enrichment", "运行GO/KEGG富集"))
            )
          ),

          h4("当前用于富集的基因数"),
          verbatimTextOutput("enrich_gene_count"),

          tabsetPanel(
            tabPanel("GO Dotplot", plotOutput("go_dotplot", height = "650px")),
            tabPanel("GO Table", tableOutput("go_table")),
            tabPanel("KEGG Dotplot", plotOutput("kegg_dotplot", height = "650px")),
            tabPanel("KEGG Table", tableOutput("kegg_table"))
          ),

          br(),
          downloadButton("downloadGO", "导出GO结果 CSV"),
          downloadButton("downloadKEGG", "导出KEGG结果 CSV")
        ),

        tabPanel(
          "Metadata",
          h4("Candidate compare columns"),
          tableOutput("candidate_table"),

          h4("Candidate filter columns"),
          tableOutput("filter_table"),

          hr(),
          h4("Group Explorer: recommended grouping columns"),
          helpText("自动扫描 metadata 中所有可用于分组的列。优先推荐 Disease / condition / phenotype / subtype 等信息量高的列；Auto_Group 只是备选。"),
          tableOutput("group_column_overview_table"),
          downloadButton("downloadGroupColumnOverview", "导出 group column overview"),

          hr(),
          h4("All group levels and sample counts"),
          tableOutput("group_backup_summary_table"),
          downloadButton("downloadGroupBackupSummary", "导出所有分组计数"),

          hr(),
          h4("Metadata preview"),
          tableOutput("meta_preview"),
          downloadButton("downloadFullMetadata", "导出完整metadata")
        )
      )
    )
  )
)

# =========================
# 8. Shiny server
# =========================
server <- function(input, output, session) {
  data_obj <- reactiveVal(NULL)
  diff_obj <- reactiveVal(NULL)
  # v8.9: explicit DEG state flag.
  # TRUE only after user clicks Run DEG and a DEG result is stored.
  # It is reset on Load / group / method / merge changes.
  deg_has_run <- reactiveVal(FALSE)
  # v11.0: parsed expression source candidates for manual raw/RPM/TPM/FPKM switching.
  expr_source_pool <- reactiveVal(NULL)

  run_status <- reactiveVal("Ready")
  run_detail <- reactiveVal("等待加载 GEO 数据")
  run_percent <- reactiveVal(0)
  run_state <- reactiveVal("idle")
  run_start_time <- reactiveVal(NULL)

  set_run_status <- function(status, detail = NULL, percent = NULL, state = "running") {
    run_status(status)
    if (!is.null(detail)) run_detail(detail)
    if (!is.null(percent)) run_percent(max(0, min(100, round(percent))))
    run_state(state)
    try(session$flushReact(), silent = TRUE)
  }

  output$run_status_ui <- renderUI({
    pct <- run_percent()
    state <- run_state()

    # Hide before first load
    if (state == "idle" && pct == 0) return(NULL)

    icon <- if (state == "done") {
      "✅"
    } else if (state == "error") {
      "⚠️"
    } else {
      "⏳"
    }

    rt <- run_start_time()
    elapsed <- if (!is.null(rt)) {
      paste0(round(as.numeric(difftime(Sys.time(), rt, units = "secs")), 1), " sec")
    } else {
      "NA"
    }

    tags$div(
      class = paste("run-status-box", state),
      tags$div(class = "run-status-title", paste0(icon, " ", run_status(), "  |  ", pct, "%")),
      tags$div(class = "run-status-detail", paste0(run_detail(), "  |  耗时: ", elapsed)),
      tags$div(
        class = "run-status-progress",
        tags$div(class = "run-status-progress-fill", style = paste0("width:", pct, "%;"))
      )
    )
  })

  observeEvent(input$run, {
    gse_id <- toupper(trimws(input$gse))
    save_dir <- GEO_ROOT

    # Soft reset: when switching GSE, clear old backend objects and old UI choices.
    # This keeps the app behavior similar to a light restart without closing Shiny.
    data_obj(NULL)
    diff_obj(NULL)
    deg_has_run(FALSE)
    updateSelectInput(session, "filter_col", choices = character(0))
    updateSelectInput(session, "filter_value", choices = character(0))
    updateSelectInput(session, "group_col", choices = character(0))
    updateSelectInput(session, "groupA", choices = character(0))
    updateSelectInput(session, "groupB", choices = character(0))
    updateSelectInput(session, "expr_source", choices = c("Auto recommended after Load" = "auto"), selected = "auto")
    expr_source_pool(NULL)

    run_start_time(Sys.time())
    set_run_status("软重置", paste0("已清空旧结果，准备加载 ", gse_id), 3, "running")

    if (isTRUE(input$clean_cache)) {
      set_run_status("清理缓存", "正在删除当前GSE下载文件和临时文件", 8, "running")
      if (dir.exists(file.path(save_dir, gse_id))) {
        unlink(file.path(save_dir, gse_id), recursive = TRUE, force = TRUE)
      }
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      try({
        tmp_files <- list.files(tempdir(), full.names = TRUE)
        if (length(tmp_files) > 0) unlink(tmp_files, recursive = TRUE, force = TRUE)
      }, silent = TRUE)
    }

    if (isTRUE(input$rebuild_parsed_cache) && file.exists(cache_file_for_gse(gse_id))) {
      unlink(cache_file_for_gse(gse_id), force = TRUE)
    }

    set_run_status("下载/读取 GEO", "正在读取缓存、Series Matrix 或 supplementary expression matrix", 18, "running")

    expr_res <- tryCatch({
      load_geo_expr_auto(
        gse_id,
        use_cache = isTRUE(input$use_parsed_cache),
        rebuild_cache = isTRUE(input$rebuild_parsed_cache)
      )
    }, error = function(e) {
      set_run_status("GEO读取失败", paste0("错误：", e$message), 100, "error")
      data_obj(list(error = paste0("GEO读取失败：", e$message)))
      return(NULL)
    })
    if (is.null(expr_res)) return()

    set_run_status("整理表达矩阵", "正在识别基因名、表达矩阵和数据类型", 45, "running")

    raw <- expr_res$raw
    meta <- expr_res$meta
    meta_original <- as.data.frame(meta)
    colnames(raw)[1] <- "Symbol"
    raw$Symbol <- as.character(raw$Symbol)
    colnames(meta) <- make.names(colnames(meta), unique = TRUE)

    meta <- add_auto_group_columns(meta)

    set_run_status("匹配样本和分组", "正在匹配 expression samples 与 metadata，并自动识别分组", 65, "running")

    match_mode <- ifelse(is.null(input$sample_match_mode), "auto_robust", input$sample_match_mode)

    # If the supplementary metadata matcher already found a strong sample ID column,
    # put that column first so generic sample matching will prefer it.
    if (!is.null(expr_res$metadata_sample_col) &&
        !is.na(expr_res$metadata_sample_col) &&
        expr_res$metadata_sample_col %in% colnames(meta)) {
      priority_col <- expr_res$metadata_sample_col
      meta <- meta[, c(priority_col, setdiff(colnames(meta), priority_col)), drop = FALSE]
    }

    mapping <- choose_sample_mapping(
      meta = meta,
      raw_sample_names = colnames(raw)[-1],
      requested_mode = match_mode
    )

    meta <- mapping$meta
    sample_col <- mapping$sample_col

    candidates <- get_group_candidates(meta)
    filter_candidates <- get_filter_candidates(meta)
    group_discovery_table <- make_group_discovery_table(meta)
    group_column_overview <- group_discovery_table
    group_backup_summary <- make_group_backup_summary(meta)
    if (nrow(candidates) == 0) {
      candidates <- filter_candidates %>% mutate(score = 0) %>% select(column, score, n_unique, example_values)
    }

    if ("Auto_Group" %in% colnames(meta) &&
        length(unique(meta$Auto_Group[!is.na(meta$Auto_Group) & meta$Auto_Group != ""])) >= 2) {
      if (!"Auto_Group" %in% candidates$column) {
        auto_row <- data.frame(
          column = "Auto_Group",
          score = 999,
          n_unique = length(unique(meta$Auto_Group[!is.na(meta$Auto_Group) & meta$Auto_Group != ""])),
          example_values = paste(unique(meta$Auto_Group[!is.na(meta$Auto_Group) & meta$Auto_Group != ""]), collapse = " | "),
          stringsAsFactors = FALSE
        )
        candidates <- bind_rows(auto_row, candidates)
      } else {
        candidates <- candidates %>%
          mutate(score = ifelse(column == "Auto_Group", 999, score)) %>%
          arrange(desc(score), n_unique, column)
      }
    }

    selected_group <- choose_default_group_column(meta, candidates, group_column_overview)

    # Important: do not hide other grouping columns.
    # Auto_Group is first, then high-scoring candidates, then all overview fields.
    group_choices <- unique(c(
      if (!is.null(group_discovery_table)) group_discovery_table$column else character(0),
      candidates$column,
      if ("Auto_Group" %in% colnames(meta)) "Auto_Group" else character(0)
    ))
    group_choices <- group_choices[group_choices %in% colnames(meta)]
    if (length(group_choices) == 0) group_choices <- colnames(meta)

    filter_choices <- unique(c(
      filter_candidates$column,
      if (!is.null(group_column_overview)) group_column_overview$column else character(0)
    ))
    filter_choices <- filter_choices[filter_choices %in% colnames(meta)]

    updateSelectInput(session, "group_col", choices = group_choices, selected = selected_group)
    updateSelectInput(session, "filter_col", choices = c("None", filter_choices), selected = "None")
    updateSelectInput(session, "filter_value", choices = "All", selected = "All")

    set_run_status("扫描表达矩阵来源", "正在建立 raw/RPM/TPM/FPKM 候选表达矩阵列表", 88, "running")
    source_pool <- tryCatch({
      build_manual_expression_source_pool(
        gse_id = gse_id,
        loaded_raw = raw,
        loaded_det = expr_res$expr_detect,
        loaded_expr_file = expr_res$expr_file
      )
    }, error = function(e) {
      message("Expression source pool build failed: ", e$message)
      list(auto = list(
        id = "auto",
        label = make_expr_source_label("Auto/current loaded", expr_res$expr_file, expr_res$expr_detect, raw),
        raw = raw,
        expr_detect = expr_res$expr_detect,
        expr_file = expr_res$expr_file,
        source = "auto/current-loaded"
      ))
    })
    expr_source_pool(source_pool)
    source_choices <- vapply(source_pool, function(z) z$label, character(1))
    names(source_choices) <- source_choices
    # selectInput expects named vector label -> value, so invert here.
    source_choices2 <- stats::setNames(names(source_pool), source_choices)
    updateSelectInput(session, "expr_source", choices = source_choices2, selected = "auto")

    set_run_status("加载完成", "数据已加载并自动推荐 Seed Groups。请确认表达矩阵来源/Final Groups，然后点击 Run DEG。", 95, "running")

    data_obj(list(
      raw = raw,
      meta = meta,
      gse_id = gse_id,
      sample_col = sample_col,
      sample_mapping_mode_used = mapping$mode_used,
      sample_mapping_warning = mapping$warning,
      sample_mapping_matched_n = mapping$matched_n,
      sample_mapping_audit = build_expression_metadata_map(meta, sample_col, selected_group, colnames(raw)[-1]),
      candidates = candidates,
      filter_candidates = filter_candidates,
      group_discovery_table = group_discovery_table,
      group_column_overview = group_column_overview,
      group_backup_summary = group_backup_summary,
      meta_original = meta_original,
      expr_file = expr_res$expr_file,
      meta_file = expr_res$meta_file,
      metadata_source = expr_res$metadata_source,
      metadata_matched_n = expr_res$metadata_matched_n,
      metadata_sample_col = expr_res$metadata_sample_col,
      download_msg = expr_res$msg,
      expr_detect = expr_res$expr_detect,
      loaded_from_cache = ifelse(is.null(expr_res$loaded_from_cache), FALSE, expr_res$loaded_from_cache)
    ))

    set_run_status("加载完成", paste0(gse_id, " 已加载。请手动确认分组后点击 Run DEG。"), 100, "done")
  })

  observeEvent(input$filter_col, {
    diff_obj(NULL)
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return()
    if (is.null(input$filter_col) || input$filter_col == "None") {
      updateSelectInput(session, "filter_value", choices = "All", selected = "All")
      return()
    }
    vals <- unique(clean_group_value(obj$meta[[input$filter_col]]))
    vals <- vals[!is.na(vals) & vals != "" & vals != "NA"]
    vals <- sort(vals)
    updateSelectInput(session, "filter_value", choices = vals, selected = vals[1])
  })

  selected_sample_meta_raw <- reactive({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)
    req(input$group_col)
    meta <- obj$meta
    sample_col <- obj$sample_col
    group_col <- input$group_col
    if (!group_col %in% colnames(meta)) return(NULL)
    match_mode <- ifelse(is.null(input$sample_match_mode), "auto_robust", input$sample_match_mode)

    meta2 <- meta %>%
      select(sample = all_of(sample_col), group = all_of(group_col)) %>%
      mutate(
        sample = if (sample_col == "Expression_Sample_Order") {
          as.character(sample)
        } else {
          standardize_sample_for_match(sample, colnames(obj$raw)[-1], mode = ifelse(match_mode == "auto_robust", "clean", match_mode))
        },
        group = clean_group_value(group)
      ) %>%
      filter(!is.na(sample), !is.na(group), group != "", group != "NA")

    if (!is.null(input$filter_col) && input$filter_col != "None" &&
        !is.null(input$filter_value) && input$filter_value != "All" &&
        input$filter_col %in% colnames(meta)) {
      filter_df <- data.frame(
        sample = if (sample_col == "Expression_Sample_Order") {
          as.character(meta[[sample_col]])
        } else {
          standardize_sample_for_match(meta[[sample_col]], colnames(obj$raw)[-1], mode = ifelse(match_mode == "auto_robust", "clean", match_mode))
        },
        filter_value = clean_group_value(meta[[input$filter_col]]),
        stringsAsFactors = FALSE
      )
      meta2 <- left_join(meta2, filter_df, by = "sample") %>%
        filter(filter_value == input$filter_value) %>%
        select(sample, group)
    }
    meta2
  })

  selected_sample_meta <- reactive({
    meta2 <- selected_sample_meta_raw()
    if (is.null(meta2) || nrow(meta2) == 0) return(meta2)

    apply_manual_group_merge(
      meta2,
      enable_merge = isTRUE(input$enable_group_merge),
      mergeA_groups = input$mergeA_groups,
      mergeB_groups = input$mergeB_groups,
      mergeA_name = input$mergeA_name,
      mergeB_name = input$mergeB_name
    )
  })

  analysis_group_names <- reactive({
    get_analysis_group_names(
      enable_merge = isTRUE(input$enable_group_merge),
      groupA = input$groupA,
      groupB = input$groupB,
      mergeA_name = input$mergeA_name,
      mergeB_name = input$mergeB_name
    )
  })

  # v11.8.2: expression-aware group helper.
  # Metadata selected samples are not necessarily present in the active expression matrix.
  # All default group selection, dashboard preview, and DEG safety checks should use
  # expression-available samples as the source of truth.
  get_expression_available_meta <- function(meta2, raw_mat) {
    if (is.null(meta2) || nrow(meta2) == 0 || is.null(raw_mat) || ncol(raw_mat) < 2) {
      return(data.frame(sample = character(), group = character(), stringsAsFactors = FALSE))
    }
    meta2 <- meta2 %>% dplyr::filter(!is.na(sample), !is.na(group), group != "", group != "NA") %>%
      dplyr::distinct(sample, .keep_all = TRUE)
    if (nrow(meta2) == 0) {
      return(data.frame(sample = character(), group = character(), stringsAsFactors = FALSE))
    }
    resolved <- tryCatch({
      resolve_samples_for_deg(meta2, colnames(raw_mat)[-1])
    }, error = function(e) {
      list(kept = data.frame(sample = character(), group = character(), stringsAsFactors = FALSE),
           dropped = meta2)
    })
    kept <- resolved$kept
    if (is.null(kept) || nrow(kept) == 0) {
      return(data.frame(sample = character(), group = character(), stringsAsFactors = FALSE))
    }
    kept
  }

  expression_group_count_table <- function(meta2, raw_mat) {
    if (is.null(meta2) || nrow(meta2) == 0) {
      return(data.frame(group = character(), n_metadata_selected = integer(),
                        n_expression_available = integer(), stringsAsFactors = FALSE))
    }
    meta_counts <- meta2 %>% dplyr::count(group, name = "n_metadata_selected")
    expr_meta <- get_expression_available_meta(meta2, raw_mat)
    expr_counts <- if (!is.null(expr_meta) && nrow(expr_meta) > 0) {
      expr_meta %>% dplyr::count(group, name = "n_expression_available")
    } else {
      data.frame(group = character(), n_expression_available = integer(), stringsAsFactors = FALSE)
    }
    out <- dplyr::left_join(meta_counts, expr_counts, by = "group")
    out$n_expression_available[is.na(out$n_expression_available)] <- 0L
    out
  }

  pick_default_groups_expression_aware <- function(meta2, raw_mat) {
    tab <- expression_group_count_table(meta2, raw_mat)
    tab <- tab[tab$n_expression_available > 0, , drop = FALSE]
    if (nrow(tab) < 2) {
      groups_all <- sort(unique(meta2$group))
      return(pick_default_groups(groups_all))
    }
    groups_avail <- sort(unique(tab$group))

    control_hits <- groups_avail[grepl("healthy|control|normal|sham|vehicle|baseline", groups_avail, ignore.case = TRUE)]
    disease_hits <- groups_avail[grepl("sepsis|septic|case|patient|disease|shock|ards|aki|sirs|infection|covid", groups_avail, ignore.case = TRUE)]
    disease_hits <- setdiff(disease_hits, control_hits)

    if (length(control_hits) > 0 && length(disease_hits) > 0) {
      # Prefer the largest expression-covered control and disease groups.
      ctrl_tab <- tab[tab$group %in% control_hits, , drop = FALSE]
      dis_tab <- tab[tab$group %in% disease_hits, , drop = FALSE]
      A <- ctrl_tab$group[order(-ctrl_tab$n_expression_available, ctrl_tab$group)][1]
      B <- dis_tab$group[order(-dis_tab$n_expression_available, dis_tab$group)][1]
      return(list(A = A, B = B))
    }

    tab <- tab[order(-tab$n_expression_available, tab$group), , drop = FALSE]
    list(A = tab$group[1], B = tab$group[2])
  }


  observeEvent(selected_sample_meta_raw(), {
    diff_obj(NULL)
    df <- selected_sample_meta_raw()
    obj <- data_obj()
    if (is.null(df) || nrow(df) == 0 || is.null(obj) || !is.null(obj$error)) return()
    groups <- sort(unique(df$group))

    # v11.8.2: choose defaults from groups that are actually present in the expression matrix.
    # Keep all metadata groups in the dropdown, but do not auto-select a group with 0 expression samples.
    defaults <- pick_default_groups_expression_aware(df, obj$raw)
    selected_A <- defaults$A
    selected_B <- defaults$B

    if (identical(selected_A, selected_B) && length(groups) >= 2) {
      available_tab <- expression_group_count_table(df, obj$raw)
      available_groups <- available_tab$group[available_tab$n_expression_available > 0]
      selected_B <- setdiff(available_groups, selected_A)[1]
      if (is.na(selected_B) || length(selected_B) == 0) selected_B <- setdiff(groups, selected_A)[1]
    }

    updateSelectInput(session, "groupA", choices = groups, selected = selected_A)
    updateSelectInput(session, "groupB", choices = groups, selected = selected_B)

    # Manual merge choices are always based on original groups, but default mergeB
    # also prefers expression-covered disease/case groups.
    available_tab <- expression_group_count_table(df, obj$raw)
    available_groups <- available_tab$group[available_tab$n_expression_available > 0]
    b_default <- available_groups[grepl("sepsis|septic|case|patient|disease|shock|ards|aki|sirs|infection|covid", available_groups, ignore.case = TRUE)]
    b_default <- setdiff(b_default, selected_A)
    if (length(b_default) == 0) {
      b_default <- groups[grepl("sepsis|septic|case|patient|disease|shock|ards|aki|sirs|infection|covid", groups, ignore.case = TRUE)]
      b_default <- setdiff(b_default, selected_A)
    }
    if (length(b_default) == 0) b_default <- selected_B

    updateSelectizeInput(session, "mergeA_groups", choices = groups, selected = selected_A, server = TRUE)
    updateSelectizeInput(session, "mergeB_groups", choices = groups, selected = b_default, server = TRUE)

    updateTextInput(session, "mergeA_name", value = ifelse(grepl("healthy|control|normal", selected_A, ignore.case = TRUE), "Control", selected_A))
    updateTextInput(session, "mergeB_name", value = ifelse(any(grepl("sepsis|septic", b_default, ignore.case = TRUE)), "Sepsis", selected_B))
  })

  observeEvent(list(input$groupA, input$groupB, input$method_mode, input$expr_source,
                    input$filter_col, input$filter_value, input$group_col), {
    diff_obj(NULL)
    deg_has_run(FALSE)
    set_run_status("等待手动DEG", "分组、筛选、表达矩阵来源或分析方法已改变。请确认后点击 Run DEG。", 100, "done")
  }, ignoreInit = TRUE)

  # -------------------------
  # Manual differential analysis
  # DEG only runs after user clicks Run DEG.
  # This prevents Load/group UI updates from repeatedly triggering limma/DESeq2.
  # -------------------------

  observeEvent(list(input$enable_group_merge, input$mergeA_groups, input$mergeB_groups, input$mergeA_name, input$mergeB_name), {
    diff_obj(NULL)
    deg_has_run(FALSE)
    set_run_status("等待手动DEG", "手动合并分组设置已改变。请确认后点击 Run DEG。", 100, "done")
  }, ignoreInit = TRUE)


  observeEvent(input$run_deg, {
    # v8.9 audited: this is the ONLY place that calls run_auto_differential().
    # Load, group changes, method changes, and merge changes only clear old DEG results.

    obj <- data_obj()
    meta2 <- selected_sample_meta()

    if (is.null(obj) || !is.null(obj$error)) {
      diff_obj(NULL)
      return()
    }

    if (is.null(meta2) || nrow(meta2) == 0) {
      set_run_status("差异分析未运行", "没有可用样本用于差异分析。请检查样本匹配和分组。", 100, "error")
      diff_obj(list(error = "没有可用样本用于差异分析。请检查样本匹配和分组。"))
      return()
    }

    req(input$groupA, input$groupB, input$method_mode)

    names_run <- analysis_group_names()
    groupA_run <- names_run$A
    groupB_run <- names_run$B

    raw <- obj$raw
    det_for_deg <- obj$expr_detect
    deg_expr_file_used <- obj$expr_file
    deg_source_switched <- FALSE
    deg_source_note <- "Using Auto/current loaded expression matrix."

    # v11.0 manual expression source selector.
    # Auto keeps the original app logic (the Load-stage recommended source, usually raw count when present).
    # Manual choices allow switching to submitter RPM/TPM/FPKM/raw matrices without any GSE-specific rules.
    source_pool <- expr_source_pool()
    selected_src <- get_selected_expression_source(source_pool, input$expr_source)
    if (!is.null(selected_src)) {
      raw <- selected_src$raw
      det_for_deg <- selected_src$expr_detect
      deg_expr_file_used <- selected_src$expr_file
      deg_source_switched <- !identical(ifelse(is.null(input$expr_source), "auto", input$expr_source), "auto")
      deg_source_note <- paste0(
        "Manual/Auto expression source selected: ", selected_src$label, "\n",
        "Data type: ", ifelse(is.null(det_for_deg$type), "unknown", det_for_deg$type), "\n",
        "This selected matrix is used consistently for DEG, single-gene plot, heatmap, PCA and correlation after Run DEG."
      )

      # Rebuild selected metadata against the selected matrix's own sample names.
      # This is essential when raw and RPM/TPM matrices have different sample columns.
      match_mode2 <- ifelse(is.null(input$sample_match_mode), "auto_robust", input$sample_match_mode)
      bm_src <- build_selected_meta_for_candidate(
        meta = obj$meta,
        raw_candidate = raw,
        sample_match_mode = match_mode2,
        group_col = input$group_col,
        filter_col = input$filter_col,
        filter_value = input$filter_value,
        enable_merge = isTRUE(input$enable_group_merge),
        mergeA_groups = input$mergeA_groups,
        mergeB_groups = input$mergeB_groups,
        mergeA_name = input$mergeA_name,
        mergeB_name = input$mergeB_name
      )
      if (!is.null(bm_src$meta2) && nrow(bm_src$meta2) > 0) {
        meta2 <- bm_src$meta2
        message("v11.8.13 filter-aware selected metadata after source rebuild: ", nrow(meta2), " samples")
        if (!is.null(input$filter_col) && input$filter_col != "None" &&
            !is.null(input$filter_value) && input$filter_value != "All") {
          message("v11.8.13 active filter applied to DEG: ", input$filter_col, " = ", input$filter_value)
        }
      }
    }

    # v11.0: automatic coverage-aware switching is disabled by default.
    # Users can manually choose raw/RPM/TPM/FPKM in Expression source.
    if (FALSE && !is.null(input$method_mode) && input$method_mode %in% c("limma", "DESeq2")) {
      set_run_status("选择最佳数据源", paste0("Force ", input$method_mode, "：按当前分组覆盖率评估候选表达矩阵"), 88, "running")
      match_mode2 <- ifelse(is.null(input$sample_match_mode), "auto_robust", input$sample_match_mode)
      alt_src <- tryCatch(load_best_expression_source_for_current_deg(
        gse_id = obj$gse_id,
        method_preference = input$method_mode,
        meta = obj$meta,
        sample_match_mode = match_mode2,
        group_col = input$group_col,
        filter_col = input$filter_col,
        filter_value = input$filter_value,
        enable_merge = isTRUE(input$enable_group_merge),
        mergeA_groups = input$mergeA_groups,
        mergeB_groups = input$mergeB_groups,
        mergeA_name = input$mergeA_name,
        mergeB_name = input$mergeB_name,
        groupA_run = groupA_run,
        groupB_run = groupB_run
      ), error = function(e) {
        message("Coverage-aware source selection failed: ", e$message)
        NULL
      })

      if (!is.null(alt_src) && !is.null(alt_src$raw)) {
        current_score <- tryCatch({
          current_sc <- score_candidate_for_current_deg(
            res = list(raw = raw, expr_detect = det_for_deg),
            file = deg_expr_file_used,
            source = "currently-loaded",
            method_preference = input$method_mode,
            meta = obj$meta,
            sample_match_mode = match_mode2,
            group_col = input$group_col,
            filter_col = input$filter_col,
            filter_value = input$filter_value,
            enable_merge = isTRUE(input$enable_group_merge),
            mergeA_groups = input$mergeA_groups,
            mergeB_groups = input$mergeB_groups,
            mergeA_name = input$mergeA_name,
            mergeB_name = input$mergeB_name,
            groupA_run = groupA_run,
            groupB_run = groupB_run
          )
          if (is.null(current_sc)) Inf else current_sc$sort_score
        }, error = function(e) Inf)

        if (alt_src$sort_score < current_score || input$method_mode %in% c("limma", "DESeq2")) {
          raw <- alt_src$raw
          det_for_deg <- alt_src$expr_detect
          deg_expr_file_used <- alt_src$expr_file
          deg_source_switched <- TRUE
          deg_source_note <- alt_src$msg
          meta2 <- alt_src$meta2
          message("Coverage-aware DEG source switched to: ", alt_src$expr_file)
          message(alt_src$msg)
        }
      } else {
        deg_source_note <- paste0(
          "Coverage-aware source selection did not find a compatible better matrix for ", input$method_mode,
          "; using currently loaded matrix."
        )
        message(deg_source_note)
      }
    }

    df2 <- meta2 %>%
      filter(group %in% c(groupA_run, groupB_run)) %>%
      distinct(sample, .keep_all = TRUE)

    # v8.5 diagnostics: keep the selected metadata samples BEFORE expression matching.
    # This tells the user when selected/merged metadata groups contain samples that are
    # not present in the actual expression matrix used for DEG (common with NCBI TPM/raw matrices).
    selected_before_match <- df2
    selected_before_counts <- table(selected_before_match$group)

    if (length(unique(df2$group)) < 2) {
      set_run_status("差异分析未运行", "差异分析需要两个组。", 100, "error")
      diff_obj(list(error = "差异分析需要两个组。"))
      return()
    }

    # v10.3 FIX: do not use intersect(df2$sample, colnames(raw)) here.
    # That loses samples when metadata IDs and expression column names differ slightly.
    # Instead, resolve every selected metadata sample to the actual expression column
    # using the same alias/GSM matching logic as Load.
    resolved_deg_samples <- resolve_samples_for_deg(df2, colnames(raw)[-1])

    # v10.5: Order repair is disabled by default.
    # We only trust explicit sample-name/GSM matching for DEG.
    order_repair_note <- "disabled in v10.5; using cleaned/GSM sample-name matching only"
    order_repaired_n <- 0

    df2 <- resolved_deg_samples$kept
    dropped_before_match <- resolved_deg_samples$dropped
    dropped_before_counts <- if (nrow(dropped_before_match) > 0) table(dropped_before_match$group) else table(character())
    matched_samples <- df2$sample

    # v10.5: Do not stop when metadata-selected samples are absent from the expression matrix.
    # Instead, continue with the true expression-available intersection and report dropped samples clearly.
    # This avoids both false order matches and endless stop/retry cycles.
    if (!is.null(dropped_before_match) && nrow(dropped_before_match) > 0) {
      message(paste0(
        "v10.5 notice: ", nrow(dropped_before_match),
        " selected metadata samples are not present in the DEG expression matrix; ",
        "continuing with expression-available samples only."
      ))
    }

    if (length(matched_samples) < 4) {
      set_run_status("样本匹配不足", "匹配样本太少，无法运行差异分析。", 100, "error")
      diff_obj(list(error = "匹配样本太少，无法运行差异分析。"))
      return()
    }

    # v11.8.2: use zero-filled expression-available counts for the two requested groups.
    # table()[missing_name] returns NA; this explicitly prevents a group with 0 expression
    # samples from slipping through or producing an unclear NA error.
    group_counts <- table(factor(df2$group, levels = c(groupA_run, groupB_run)))
    names(group_counts) <- c(groupA_run, groupB_run)

    if (any(as.integer(group_counts) < 2)) {
      msg_counts <- paste0(names(group_counts), "=", as.integer(group_counts), collapse = "; ")
      msg <- paste0(
        "当前选择的分组在表达矩阵中样本不足，不能运行DEG。",
        " 每组至少需要2个表达样本。Expression-available counts: ", msg_counts,
        "。请更换 Group A/B，或检查该 supplementary matrix 是否包含这些组。"
      )
      set_run_status("样本数不足", msg, 100, "error")
      diff_obj(list(
        error = msg,
        selected_before_match_n = nrow(selected_before_match),
        selected_before_match_counts = as.list(selected_before_counts),
        dropped_due_to_missing_expression_n = nrow(dropped_before_match),
        dropped_due_to_missing_expression_counts = as.list(dropped_before_counts),
        dropped_due_to_missing_expression = dropped_before_match
      ))
      return()
    }

    df2 <- df2[match(matched_samples, df2$sample), , drop = FALSE]

    set_run_status("差异分析中", paste0("正在分析：", groupB_run, " vs ", groupA_run), 92, "running")

    tryCatch({
      res <- run_auto_differential(
        raw = raw,
        df2 = df2,
        groupA = groupA_run,
        groupB = groupB_run,
        user_method = input$method_mode,
        det = det_for_deg
      )

      res$groupA <- groupA_run
      res$groupB <- groupB_run
      res$coldata <- data.frame(sample = df2$sample, group = df2$group, stringsAsFactors = FALSE)
      res$run_time <- as.character(Sys.time())
      res$filter_col <- input$filter_col
      res$filter_value <- input$filter_value
      res$group_col <- input$group_col
      res$method_mode <- input$method_mode
      res$expr_detect_used <- det_for_deg$type
      res$deg_expr_file <- deg_expr_file_used
      res$deg_source_switched <- deg_source_switched
      res$deg_source_note <- deg_source_note
      res$deg_sample_n <- nrow(df2)
      res$deg_group_counts <- as.list(table(df2$group))
      res$selected_before_match_n <- nrow(selected_before_match)
      res$selected_before_match_counts <- as.list(selected_before_counts)
      res$dropped_due_to_missing_expression_n <- nrow(dropped_before_match)
      res$dropped_due_to_missing_expression_counts <- as.list(dropped_before_counts)
      res$dropped_due_to_missing_expression <- dropped_before_match

      diff_obj(res)
      deg_has_run(TRUE)

      rt <- run_start_time()
      elapsed <- if (!is.null(rt)) round(as.numeric(difftime(Sys.time(), rt, units = "secs")), 1) else NA
      set_run_status(
        "分析完成",
        paste0("已完成 ", res$method, " 分析", if (!is.na(elapsed)) paste0("，耗时 ", elapsed, " 秒") else ""),
        100,
        "done"
      )
    }, error = function(e) {
      set_run_status("差异分析失败", paste0("错误：", e$message), 100, "error")
      diff_obj(list(error = paste0("差异分析运行失败：", e$message)))
      deg_has_run(FALSE)
    })
  })

  gene_expr_df <- reactive({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)
    req(input$gene)
    gene <- toupper(trimws(input$gene))

    # Prefer differential-analysis processed expression matrix when available.
    dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
    if (!is.null(dobj) && is.null(dobj$error) && (!is.null(dobj$plot_expr_mat_all) || !is.null(dobj$plot_expr_mat))) {
      mat <- if (!is.null(dobj$plot_expr_mat_all)) dobj$plot_expr_mat_all else dobj$plot_expr_mat
      gene_row <- rownames(mat)[toupper(rownames(mat)) == gene]
      if (length(gene_row) == 0) return(NULL)
      gene_row <- gene_row[1]
      return(data.frame(sample = colnames(mat), expression = as.numeric(mat[gene_row, ]), stringsAsFactors = FALSE))
    }

    # Before differential analysis: use detected safe plotting scale.
    raw <- obj$raw
    gene_row <- raw %>% filter(toupper(Symbol) == gene)
    if (nrow(gene_row) == 0) return(NULL)
    expr <- gene_row[1, -1] %>% unlist() %>% as.numeric()
    det <- obj$expr_detect
    if (is.null(det) || det$type == "raw_count" || det$type == "normalized_expression_FPKM_TPM_or_similar") {
      expr <- log2(pmax(expr, 0) + 1)
    }
    data.frame(sample = colnames(raw)[-1], expression = expr, stringsAsFactors = FALSE)
  })

  analysis_df <- reactive({
    df_expr <- gene_expr_df()
    if (is.null(df_expr)) return(NULL)

    # v10.9 adaptive active-matrix logic:
    # After Run DEG, all single-gene plots, summaries, and quick comparisons use
    # the exact matrix and sample/group mapping selected by DEG. This keeps raw,
    # TPM/RPM/FPKM, and candidate-switched sources synchronized.
    dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
    if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$coldata)) {
      meta2 <- as.data.frame(dobj$coldata, stringsAsFactors = FALSE)
      if (!all(c("sample", "group") %in% colnames(meta2))) return(NULL)
      meta2$sample <- as.character(meta2$sample)
      meta2$group <- as.character(meta2$group)
      df <- left_join(df_expr, meta2[, c("sample", "group"), drop = FALSE], by = "sample") %>%
        filter(!is.na(expression), !is.na(group), group != "", group != "NA")
      if (nrow(df) == 0) return(NULL)
      return(df)
    }

    # Before Run DEG, use the currently loaded matrix and current UI-selected metadata.
    meta2 <- selected_sample_meta()
    if (is.null(meta2)) return(NULL)
    df <- left_join(df_expr, meta2, by = "sample") %>%
      filter(!is.na(expression), !is.na(group), group != "", group != "NA")
    if (nrow(df) == 0) return(NULL)
    df
  })

  message_text <- reactive({
    obj <- data_obj()
    add <- function(lines, ...) c(lines, paste0(...))
    add_blank <- function(lines) c(lines, "")
    cap <- function(expr) paste(capture.output(expr), collapse = "\n")

    out <- character()

    if (is.null(obj)) {
      return("点击 1. Load / Reload GSE 加载数据。加载后只推荐 Seed Groups；如启用Manual Merge，请确认 Final Groups 后点击 Run DEG。")
    }
    if (!is.null(obj$error)) {
      return(as.character(obj$error))
    }

    out <- add(out, obj$download_msg)
    out <- add(out, "Expression file: ", obj$expr_file)
    out <- add(out, "Metadata file: ", obj$meta_file)
    if (!is.null(obj$metadata_source)) out <- add(out, "Metadata source: ", obj$metadata_source)
    if (!is.null(obj$metadata_matched_n)) out <- add(out, "Metadata matched samples: ", obj$metadata_matched_n, " / ", ncol(obj$raw) - 1)
    if (!is.null(obj$metadata_sample_col)) out <- add(out, "Metadata matched sample column: ", obj$metadata_sample_col)
    out <- add(out, "Sample column used by app: ", obj$sample_col)
    out <- add(out, "Expression samples example: ", paste(head(colnames(obj$raw)[-1], 10), collapse = ", "))
    out <- add(out, "Metadata sample example: ", paste(head(as.character(obj$meta[[obj$sample_col]]), 10), collapse = ", "))
    mm <- ifelse(is.null(input$sample_match_mode), "auto_robust", input$sample_match_mode)
    out <- add(out, "Requested sample matching mode: ", mm)
    out <- add(out, "Actual sample mapping used: ", obj$sample_mapping_mode_used)
    out <- add(out, "Sample-name matches: ", obj$sample_mapping_matched_n, " / ", ncol(obj$raw) - 1)
    if (!is.null(obj$sample_mapping_warning) && !is.na(obj$sample_mapping_warning)) {
      out <- add(out, "Sample mapping warning: ", obj$sample_mapping_warning)
    }
    if (!is.null(input$group_col) && input$group_col %in% colnames(obj$meta)) {
      map_now <- build_expression_metadata_map(obj$meta, obj$sample_col, input$group_col, colnames(obj$raw)[-1])
      out <- add(out, "Expression sample -> metadata row -> selected group mapping preview:")
      show_cols <- c("expr_sample", "metadata_sample", "geo_accession", "title", "group", "mapping_status")
      show_cols <- show_cols[show_cols %in% colnames(map_now)]
      out <- add(out, cap(print(head(map_now[, show_cols, drop = FALSE], 30), row.names = FALSE)))
      map_count <- map_now[map_now$mapping_status == "matched" & !is.na(map_now$group) & map_now$group != "", , drop = FALSE]
      out <- add(out, "Expression-matched counts for selected group column:")
      if (nrow(map_count) > 0) {
        out <- add(out, cap(print(table(map_count$group, useNA = "ifany"))))
      } else {
        out <- add(out, "No matched expression samples for this group column.")
      }
    }
    if ("Auto_Group" %in% colnames(obj$meta) && isTRUE(input$group_col == "Auto_Group")) {
      auto_map <- build_expression_metadata_map(obj$meta, obj$sample_col, "Auto_Group", colnames(obj$raw)[-1])
      auto_map2 <- auto_map[auto_map$mapping_status == "matched" & !is.na(auto_map$group) & auto_map$group != "", , drop = FALSE]
      out <- add(out, "Auto_Group expression-available counts only:")
      if (nrow(auto_map2) > 0) {
        out <- add(out, cap(print(table(auto_map2$group, useNA = "ifany"))))
      } else {
        out <- add(out, "No expression-matched Auto_Group samples.")
      }
      out <- add(out, "Note: Auto_Group is only shown when selected as Compare column. Prefer original metadata columns such as condition/disease/state when available.")
    }

    out <- add(out, "Gene: ", toupper(trimws(input$gene)))
    out <- add(out, "Filter column: ", input$filter_col)
    out <- add(out, "Filter value: ", input$filter_value)
    out <- add(out, "Group column: ", input$group_col)
    out <- add(out, "Selected method mode: ", input$method_mode)

    if (isTRUE(input$enable_group_merge)) {
      effectiveA <- ifelse(!is.null(input$mergeA_name) && nzchar(trimws(input$mergeA_name)), trimws(input$mergeA_name), input$groupA)
      effectiveB <- ifelse(!is.null(input$mergeB_name) && nzchar(trimws(input$mergeB_name)), trimws(input$mergeB_name), input$groupB)
      out <- add(out, "Manual merge enabled: YES")
      out <- add(out, "Seed Group A selector: ", input$groupA)
      out <- add(out, "Seed Group B selector: ", input$groupB)
      out <- add(out, "Final Group A used for DEG/plots: ", effectiveA)
      out <- add(out, "Final Group B used for DEG/plots: ", effectiveB)
      out <- add(out, "Note: Seed Group dropdowns only initialize choices; Dashboard and DEG use Final Group labels when manual merge is enabled.")
    }

    names_preview <- analysis_group_names()
    meta_preview <- selected_sample_meta()
    if (!is.null(meta_preview) && nrow(meta_preview) > 0 && !isTRUE(deg_has_run())) {
      meta_preview2 <- meta_preview %>% filter(group %in% c(names_preview$A, names_preview$B))
      out <- add_blank(out)
      out <- add(out, "--- Pre-DEG metadata selection preview ---")
      out <- add(out, "Run DEG has not been executed after the latest setting change.")
      out <- add(out, "Selected metadata samples before expression matching: ", nrow(meta_preview2))
      out <- add(out, "Selected metadata group counts:")
      out <- add(out, cap(print(table(meta_preview2$group))))
      expr_tab_preview <- expression_group_count_table(meta_preview2, obj$raw)
      out <- add(out, "Expression-available group counts before DEG:")
      out <- add(out, cap(print(expr_tab_preview)))
      zero_groups <- expr_tab_preview$group[expr_tab_preview$n_expression_available == 0]
      if (length(zero_groups) > 0) {
        out <- add(out, "Warning: group(s) with 0 expression samples: ", paste(zero_groups, collapse = ", "))
      }
      out <- add(out, "Note: DEG will only use expression-available samples. Groups with 0 expression samples cannot be analyzed.")
    }

    dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
    if (!is.null(dobj)) {
      out <- add_blank(out)
      out <- add(out, "--- Differential analysis status ---")
      if (!is.null(dobj$error)) {
        out <- add(out, dobj$error)
      } else {
        if (!is.null(dobj$deg_expr_file)) out <- add(out, "DEG / active single-gene expression file actually used: ", dobj$deg_expr_file)
        if (!is.null(dobj$expr_detect_used)) out <- add(out, "DEG / active single-gene input type actually used: ", dobj$expr_detect_used)
        if (!is.null(dobj$deg_source_note)) out <- add(out, "DEG source note: ", dobj$deg_source_note)
        if (!is.null(dobj$selected_before_match_n)) {
          out <- add(out, "Selected metadata samples before expression matching: ", dobj$selected_before_match_n)
          if (!is.null(dobj$selected_before_match_counts)) {
            out <- add(out, "Selected group counts before expression matching:")
            out <- add(out, cap(print(unlist(dobj$selected_before_match_counts))))
          }
        }
        if (!is.null(dobj$coldata)) {
          out <- add(out, "DEG samples actually used: ", nrow(dobj$coldata))
          out <- add(out, "DEG group counts actually used:")
          out <- add(out, cap(print(table(dobj$coldata$group))))
        }
        if (!is.null(dobj$dropped_due_to_missing_expression_n) && dobj$dropped_due_to_missing_expression_n > 0) {
          out <- add(out, "Samples dropped because they are not present in the DEG expression matrix: ", dobj$dropped_due_to_missing_expression_n)
          if (!is.null(dobj$dropped_due_to_missing_expression_counts)) {
            out <- add(out, "Dropped sample counts by group:")
            out <- add(out, cap(print(unlist(dobj$dropped_due_to_missing_expression_counts))))
          }
          if (!is.null(dobj$dropped_due_to_missing_expression) && nrow(dobj$dropped_due_to_missing_expression) > 0) {
            out <- add(out, "Dropped sample IDs (first 30):")
            out <- add(out, cap(print(head(dobj$dropped_due_to_missing_expression[, c("sample", "group"), drop = FALSE], 30))))
          }
        }

        out <- add(out, paste0(dobj$method, "完成。"))
        if (!is.null(dobj$groupB) && !is.null(dobj$groupA)) out <- add(out, "Contrast: ", dobj$groupB, " vs ", dobj$groupA)
        if (!is.null(dobj$coldata)) out <- add(out, "Samples: ", paste(as.integer(table(dobj$coldata$group)), collapse = " vs "))
        if (!is.null(dobj$plot_expr_mat)) out <- add(out, "Genes after filtering: ", nrow(dobj$plot_expr_mat))
        if (!is.null(dobj$expr_detect_used)) out <- add(out, "Input判断：", dobj$expr_detect_used)
        if (identical(dobj$method, "limma")) {
          if (!is.null(dobj$limma_transform_used)) {
            out <- add(out, "limma输入转换：", dobj$limma_transform_used)
          } else {
            out <- add(out, "limma输入转换：根据表达矩阵类型自动决定")
          }
        }
        if (!is.null(dobj$run_time)) out <- add(out, "Run time: ", dobj$run_time)
        out <- add(out, "Run trigger: manual Run DEG button")
      }
    }

    # Final safety: remove only trailing naked UI/group labels. Informative lines with ':' are kept.
    known_group_values <- character()
    try({
      if (!is.null(obj) && is.null(obj$error) && !is.null(input$group_col) && input$group_col %in% colnames(obj$meta)) {
        known_group_values <- unique(trimws(as.character(obj$meta[[input$group_col]])))
      }
    }, silent = TRUE)
    tail_values <- unique(trimws(as.character(c(
      input$groupA, input$groupB,
      input$mergeA_name, input$mergeB_name,
      input$mergeA_groups, input$mergeB_groups,
      known_group_values
    ))))
    tail_values <- tail_values[!is.na(tail_values) & nzchar(tail_values)]
    tail_values_lower <- tolower(tail_values)
    guard <- 0
    while (length(out) > 0 && guard < 100) {
      z <- trimws(tail(out, 1))
      if (nzchar(z) && !grepl(":|：", z) && tolower(z) %in% tail_values_lower) {
        out <- head(out, -1)
        guard <- guard + 1
      } else {
        break
      }
    }

    paste(out, collapse = "\n")
  })

  output$message_ui <- renderUI({
    tags$pre(
      style = "white-space: pre-wrap; word-break: break-word; max-height: 700px; overflow-y: auto; background: #f8f8f8; border: 1px solid #ddd; padding: 10px;",
      message_text()
    )
  })

  output$candidate_table <- renderTable({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)
    obj$candidates %>% mutate(Recommended = ifelse(row_number() == 1, "YES", "")) %>%
      select(Recommended, column, score, n_unique, example_values)
  })

  output$filter_table <- renderTable({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)
    obj$filter_candidates
  })


  output$group_column_overview_table <- renderTable({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)

    if (!is.null(obj$group_discovery_table)) {
      return(obj$group_discovery_table)
    }

    if (!is.null(obj$group_column_overview)) {
      return(obj$group_column_overview)
    }

    make_group_column_overview(obj$meta)
  })

  output$group_backup_summary_table <- renderTable({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)

    df <- if (!is.null(obj$group_backup_summary)) {
      obj$group_backup_summary
    } else {
      make_group_backup_summary(obj$meta)
    }

    if (is.null(df)) return(NULL)
    head(df, 800)
  })

  output$meta_preview <- renderTable({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)
    head(obj$meta, 30)
  })

  output$downloadGroupColumnOverview <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_group_column_overview.csv")
    },
    content = function(file) {
      obj <- data_obj()
      if (is.null(obj) || !is.null(obj$error)) {
        write.csv(data.frame(message = "No metadata"), file, row.names = FALSE)
      } else {
        df <- if (!is.null(obj$group_column_overview)) obj$group_column_overview else make_group_column_overview(obj$meta)
        write.csv(df, file, row.names = FALSE)
      }
    }
  )

  output$downloadGroupBackupSummary <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_all_group_counts.csv")
    },
    content = function(file) {
      obj <- data_obj()
      if (is.null(obj) || !is.null(obj$error)) {
        write.csv(data.frame(message = "No metadata"), file, row.names = FALSE)
      } else {
        df <- if (!is.null(obj$group_backup_summary)) obj$group_backup_summary else make_group_backup_summary(obj$meta)
        write.csv(df, file, row.names = FALSE)
      }
    }
  )

  output$downloadFullMetadata <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_full_metadata_used_by_app.csv")
    },
    content = function(file) {
      obj <- data_obj()
      if (is.null(obj) || !is.null(obj$error)) {
        write.csv(data.frame(message = "No metadata"), file, row.names = FALSE)
      } else {
        write.csv(obj$meta, file, row.names = FALSE)
      }
    }
  )

  # v11.8.11: Export the active merged expression matrix for Meta/WGCNA/GSVA/CIBERSORT.
  # This uses the currently selected Expression source when available; otherwise it exports obj$raw.
  active_expression_source_for_export <- reactive({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error) || is.null(obj$raw)) return(NULL)

    src <- tryCatch({
      get_selected_expression_source(expr_source_pool(), ifelse(is.null(input$expr_source), "auto", input$expr_source))
    }, error = function(e) NULL)

    if (!is.null(src) && !is.null(src$raw)) {
      return(list(
        raw = src$raw,
        expr_detect = src$expr_detect,
        expr_file = src$expr_file,
        source_label = src$label
      ))
    }

    list(
      raw = obj$raw,
      expr_detect = obj$expr_detect,
      expr_file = obj$expr_file,
      source_label = "Auto/current loaded expression matrix"
    )
  })

  output$downloadMergedExpressionMatrix <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_merged_expression_matrix_current_source.csv")
    },
    content = function(file) {
      src <- active_expression_source_for_export()
      if (is.null(src) || is.null(src$raw)) {
        write.csv(data.frame(message = "No expression matrix available. Please Load a GEO dataset first."), file, row.names = FALSE)
      } else {
        raw <- as.data.frame(src$raw, stringsAsFactors = FALSE)
        if (!"Symbol" %in% colnames(raw)) colnames(raw)[1] <- "Symbol"
        write.csv(raw, file, row.names = FALSE)
      }
    }
  )

  output$downloadExpressionAvailableMetadata <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_expression_available_group_metadata.csv")
    },
    content = function(file) {
      # v11.8.14 FIX:
      # If Run DEG has been executed, export the exact metadata used by DEG.
      # This keeps metadata synchronized with DEG/single-gene plot/exported Meta-ready matrix.
      dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
      if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$coldata)) {
        kept <- as.data.frame(dobj$coldata, stringsAsFactors = FALSE)
        if (!"metadata_sample" %in% colnames(kept)) kept$metadata_sample <- kept$sample
        if (!"expr_sample" %in% colnames(kept)) kept$expr_sample <- kept$sample
        write.csv(kept, file, row.names = FALSE)
      } else {
        src <- active_expression_source_for_export()
        meta2 <- selected_sample_meta()
        names_run <- analysis_group_names()

        if (is.null(src) || is.null(src$raw) || is.null(meta2) || nrow(meta2) == 0) {
          write.csv(data.frame(message = "No matched metadata available. Please Load a dataset and select groups first."), file, row.names = FALSE)
        } else {
          meta_sel <- meta2 %>%
            dplyr::filter(group %in% c(names_run$A, names_run$B)) %>%
            dplyr::distinct(sample, .keep_all = TRUE)

          resolved <- tryCatch({
            resolve_samples_for_deg(meta_sel, colnames(src$raw)[-1])
          }, error = function(e) list(kept = data.frame(), dropped = meta_sel))

          kept <- as.data.frame(resolved$kept, stringsAsFactors = FALSE)
          if (nrow(kept) == 0) {
            write.csv(data.frame(message = "No selected metadata samples matched the active expression matrix."), file, row.names = FALSE)
          } else {
            keep_cols <- intersect(c("sample", "group", "metadata_sample", "expr_sample"), colnames(kept))
            extra_cols <- setdiff(colnames(kept), keep_cols)
            kept <- kept[, c(keep_cols, extra_cols), drop = FALSE]
            write.csv(kept, file, row.names = FALSE)
          }
        }
      }
    }
  )

  output$downloadPlotScaleExpressionMatrix <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_log_plot_scale_expression_matrix_current_source.csv")
    },
    content = function(file) {
      # v11.8.14 FIX:
      # Prefer the exact plot-scale matrix stored by the last successful Run DEG.
      # This prevents export from re-matching a different source and breaking group/sample alignment.
      dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
      if (!is.null(dobj) && is.null(dobj$error) &&
          (!is.null(dobj$plot_expr_mat_all) || !is.null(dobj$plot_expr_mat))) {
        mat <- if (!is.null(dobj$plot_expr_mat_all)) dobj$plot_expr_mat_all else dobj$plot_expr_mat
        out <- data.frame(Symbol = rownames(mat), mat, check.names = FALSE)
        write.csv(out, file, row.names = FALSE)
      } else {
        src <- active_expression_source_for_export()
        if (is.null(src) || is.null(src$raw)) {
          write.csv(data.frame(message = "No expression matrix available. Please Load a GEO dataset first."), file, row.names = FALSE)
        } else {
          tmp_obj <- list(raw = src$raw, expr_detect = src$expr_detect)
          mat <- get_plot_expr_matrix_for_modules(tmp_obj)
          if (is.null(mat)) {
            write.csv(data.frame(message = "Could not create plot-scale expression matrix."), file, row.names = FALSE)
          } else {
            out <- data.frame(Symbol = rownames(mat), mat, check.names = FALSE)
            write.csv(out, file, row.names = FALSE)
          }
        }
      }
    }
  )

  # v11.8.14: Meta-ready export from the last actual Run DEG object.
  # Row 1: Symbol,GSM...
  # Row 2: GROUP,HEALTHY/SEPSIS...
  # Row 3+: expression values from dobj$plot_expr_mat / plot_expr_mat_all
  output$downloadMetaReadyDegMatrix <- downloadHandler(
    filename = function() {
      paste0(input$gse, "_META_READY_from_last_RunDEG_plot_scale.csv")
    },
    content = function(file) {
      dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL

      if (is.null(dobj) || !is.null(dobj$error) ||
          is.null(dobj$coldata) ||
          (is.null(dobj$plot_expr_mat_all) && is.null(dobj$plot_expr_mat))) {
        write.csv(
          data.frame(message = "Please click Run DEG successfully first, then export this Meta-ready matrix."),
          file,
          row.names = FALSE
        )
        return()
      }

      mat <- if (!is.null(dobj$plot_expr_mat_all)) dobj$plot_expr_mat_all else dobj$plot_expr_mat
      coldata <- as.data.frame(dobj$coldata, stringsAsFactors = FALSE)

      if (!all(c("sample", "group") %in% colnames(coldata))) {
        write.csv(data.frame(message = "Run DEG object has no sample/group coldata."), file, row.names = FALSE)
        return()
      }

      samples <- as.character(coldata$sample)
      samples <- samples[samples %in% colnames(mat)]
      if (length(samples) < 2) {
        write.csv(data.frame(message = "No Run DEG samples were found in the plot-scale expression matrix."), file, row.names = FALSE)
        return()
      }

      coldata <- coldata[match(samples, as.character(coldata$sample)), , drop = FALSE]
      mat <- mat[, samples, drop = FALSE]

      group_out <- toupper(trimws(as.character(coldata$group)))
      groupA_upper <- toupper(trimws(as.character(dobj$groupA)))
      groupB_upper <- toupper(trimws(as.character(dobj$groupB)))

      group_out[group_out == groupA_upper | group_out %in% c("HEALTHY", "HEALTHY_CONTROL", "HEALTHY CONTROLS", "CONTROL", "CONTROLS", "HC", "HLTY")] <- "HEALTHY"
      group_out[group_out == groupB_upper | group_out %in% c("SEPSIS", "SEPSIS_PATIENT", "SEPSIS PATIENT", "SEPTIC_SHOCK", "SEVERE_SEPSIS", "CASE", "PATIENT")] <- "SEPSIS"

      expr_df <- data.frame(Symbol = rownames(mat), mat, check.names = FALSE)
      group_row <- as.data.frame(as.list(c("GROUP", group_out)), stringsAsFactors = FALSE)
      colnames(group_row) <- colnames(expr_df)

      out <- rbind(group_row, expr_df)
      write.csv(out, file, row.names = FALSE)
    }
  )



  output$gene_deg_result <- renderTable({
    dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
    if (is.null(dobj) || !is.null(dobj$error)) return(NULL)
    gene <- toupper(trimws(input$gene))
    res <- dobj$res_df %>% filter(toupper(Symbol) == gene)
    if (nrow(res) == 0) return(NULL)
    res %>%
      mutate(
        fold_change = 2^log2FoldChange,
        direction = case_when(
          !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ paste0("Up in ", dobj$groupB),
          !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ paste0("Down in ", dobj$groupB),
          TRUE ~ "NS"
        )
      ) %>%
      select(Symbol, method, baseMean, log2FoldChange, fold_change, lfcSE, stat, pvalue, padj, direction)
  }, digits = 4)

  summary_table <- reactive({
    df <- analysis_df()
    meta2 <- selected_sample_meta()
    names_run <- analysis_group_names()

    if (is.null(df) && is.null(meta2)) return(NULL)

    # v8.8: before Run DEG, show both metadata-selected N and expression-available N.
    # This avoids the confusing situation where the Dashboard shows metadata-selected
    # samples, while the single-gene plot can only use samples present in the expression matrix.
    if (is.null(diff_obj())) {
      meta_pre <- meta2 %>% filter(group %in% c(names_run$A, names_run$B))
      meta_counts <- meta_pre %>% dplyr::count(group, name = "n_metadata_selected")

      if (!is.null(df)) {
        expr_stats <- df %>%
          filter(group %in% c(names_run$A, names_run$B)) %>%
          group_by(group) %>%
          summarise(
            n_with_expression = n(),
            mean_expression_preview = round(mean(expression, na.rm = TRUE), 3),
            median_expression_preview = round(median(expression, na.rm = TRUE), 3),
            .groups = "drop"
          )
      } else {
        expr_stats <- data.frame(group = character(), n_with_expression = integer(),
                                 mean_expression_preview = numeric(), median_expression_preview = numeric())
      }

      out <- left_join(meta_counts, expr_stats, by = "group")

      # v9.0 safety fix:
      # During Shiny startup or immediately after changing group/merge controls,
      # meta_counts can temporarily be a 0-row data frame. Assigning a scalar string
      # to a new column of a 0-row data frame triggers:
      #   replacement has 1 row, data has 0
      # Return a correctly shaped empty table instead of crashing the app.
      if (is.null(out) || nrow(out) == 0) {
        return(data.frame(
          group = character(),
          n_metadata_selected = integer(),
          n_with_expression = integer(),
          mean_expression_preview = numeric(),
          median_expression_preview = numeric(),
          note = character(),
          stringsAsFactors = FALSE
        ))
      }

      out$n_with_expression[is.na(out$n_with_expression)] <- 0
      out$note <- rep("Pre-DEG preview; click Run DEG for final DEG sample counts", nrow(out))
      return(out)
    }

    if (is.null(df)) return(NULL)
    df %>% group_by(group) %>% summarise(
      n = n(),
      mean_expression = round(mean(expression, na.rm = TRUE), 3),
      median_expression = round(median(expression, na.rm = TRUE), 3),
      .groups = "drop"
    )
  })
  output$summary <- renderTable(summary_table())

  compare_table <- reactive({
    obj <- data_obj()
    df <- analysis_df()
    meta2 <- selected_sample_meta()
    if (is.null(obj) || !is.null(obj$error)) return(NULL)
    req(input$groupA, input$groupB)
    names_run <- analysis_group_names()

    if (is.null(diff_obj())) {
      if (is.null(meta2)) return(NULL)
      meta_pre <- meta2 %>% filter(group %in% c(names_run$A, names_run$B))
      expr_pre <- if (!is.null(df)) df %>% filter(group %in% c(names_run$A, names_run$B)) else NULL
      return(data.frame(
        GSE = obj$gse_id,
        Gene = toupper(trimws(input$gene)),
        Group_column = input$group_col,
        Group_A = names_run$A,
        Group_B = names_run$B,
        N_A_metadata_selected = sum(meta_pre$group == names_run$A),
        N_B_metadata_selected = sum(meta_pre$group == names_run$B),
        N_A_with_expression = if (!is.null(expr_pre)) sum(expr_pre$group == names_run$A) else NA_integer_,
        N_B_with_expression = if (!is.null(expr_pre)) sum(expr_pre$group == names_run$B) else NA_integer_,
        Note = "Pre-DEG preview; click Run DEG for final DEG comparison",
        stringsAsFactors = FALSE
      ))
    }

    if (is.null(df)) return(NULL)
    df2 <- df %>% filter(group %in% c(names_run$A, names_run$B))
    if (length(unique(df2$group)) < 2) return(NULL)
    test <- wilcox.test(expression ~ group, data = df2)
    data.frame(
      GSE = obj$gse_id,
      Gene = toupper(trimws(input$gene)),
      Filter_column = input$filter_col,
      Filter_value = input$filter_value,
      Group_column = input$group_col,
      Group_A = names_run$A,
      Group_B = names_run$B,
      N_A = sum(df2$group == names_run$A),
      N_B = sum(df2$group == names_run$B),
      Mean_A = round(mean(df2$expression[df2$group == names_run$A], na.rm = TRUE), 3),
      Mean_B = round(mean(df2$expression[df2$group == names_run$B], na.rm = TRUE), 3),
      Median_A = round(median(df2$expression[df2$group == names_run$A], na.rm = TRUE), 3),
      Median_B = round(median(df2$expression[df2$group == names_run$B], na.rm = TRUE), 3),
      Difference_B_vs_A = round(mean(df2$expression[df2$group == names_run$B], na.rm = TRUE) -
                                  mean(df2$expression[df2$group == names_run$A], na.rm = TRUE), 3),
      Wilcoxon_P_for_plot_only = signif(test$p.value, 4)
    )
  })
  output$compare <- renderTable(compare_table())


  output$single_gene_group_stats <- renderTable({
    df <- analysis_df()
    if (is.null(df)) return(NULL)
    if (!input$show_all) {
      names_run <- analysis_group_names()
      df <- df %>% filter(group %in% c(names_run$A, names_run$B))
    }
    make_single_gene_group_stats(df)
  }, digits = 4)

  output$single_gene_comparison_stats <- renderTable({
    obj <- data_obj()
    df <- analysis_df()
    if (is.null(obj) || is.null(df)) return(NULL)
    req(input$groupA, input$groupB)

    make_single_gene_comparison_stats(
      df = df,
      groupA = analysis_group_names()$A,
      groupB = analysis_group_names()$B,
      gene = toupper(trimws(input$gene)),
      gse_id = obj$gse_id
    )
  }, digits = 4)


  output$single_gene_comparison_cards <- renderUI({
    obj <- data_obj()
    df <- analysis_df()
    if (is.null(obj) || is.null(df)) return(NULL)
    req(input$groupA, input$groupB)

    st <- make_single_gene_comparison_stats(
      df = df,
      groupA = analysis_group_names()$A,
      groupB = analysis_group_names()$B,
      gene = toupper(trimws(input$gene)),
      gse_id = obj$gse_id
    )

    if (is.null(st) || nrow(st) == 0) {
      return(tags$div("两组比较统计暂不可用。"))
    }

    logfc <- suppressWarnings(as.numeric(st$Log2FC_plot_scale_B_vs_A[1]))
    fc <- suppressWarnings(as.numeric(st$Fold_change_plot_scale_B_vs_A[1]))

    direction_html <- if (is.na(logfc)) {
      "<span style='color:#777;'>NA</span>"
    } else if (logfc > 0) {
      paste0("<span style='color:#E64B35;font-weight:bold;'>↑ Up in ", st$Group_B[1], "</span>")
    } else if (logfc < 0) {
      paste0("<span style='color:#0072B5;font-weight:bold;'>↓ Down in ", st$Group_B[1], "</span>")
    } else {
      "<span style='color:#777;font-weight:bold;'>No change</span>"
    }

    card <- function(title, value, subtitle = NULL, color = "#333") {
      tags$div(
        class = "stat-card",
        tags$div(class = "stat-card-title", title),
        tags$div(class = "stat-card-value", style = paste0("color:", color, ";"), value),
        if (!is.null(subtitle)) tags$div(class = "stat-card-subtitle", subtitle)
      )
    }

    tagList(
      tags$div(
        class = "comparison-summary-title",
        HTML(paste0(
          "<b>", st$Gene[1], "</b> | ",
          st$Group_B[1], " vs ", st$Group_A[1],
          " &nbsp;&nbsp; ", direction_html
        ))
      ),

      tags$div(
        class = "stat-card-grid",
        card("Wilcoxon P", st$Wilcoxon_P[1]),
        card("T-test P", st$T_test_P[1]),
        card("log2FC", st$Log2FC_plot_scale_B_vs_A[1],
             "plot-scale mean difference",
             ifelse(!is.na(logfc) && logfc > 0, "#E64B35", ifelse(!is.na(logfc) && logfc < 0, "#0072B5", "#333"))),
        card("Fold change", st$Fold_change_plot_scale_B_vs_A[1]),
        card("Mean difference", st$Mean_difference_B_minus_A[1], "B - A"),
        card("Median difference", st$Median_difference_B_minus_A[1], "B - A"),
        card(paste0("Mean ", st$Group_A[1]), st$Mean_A[1], paste0("n = ", st$N_A[1])),
        card(paste0("Mean ", st$Group_B[1]), st$Mean_B[1], paste0("n = ", st$N_B[1]))
      )
    )
  })

  # -------------------------
  # Plotting style presets (UI/plot appearance only)
  # -------------------------
  apply_single_gene_plot_values <- function(session, preset = "prism") {
    if (identical(preset, "nature")) {
      vals <- list(
        base_font_size = 12, axis_text_size = 10, axis_title_size = 11,
        title_font_size = 14, subtitle_font_size = 10,
        x_axis_line_width = 0.7, y_axis_line_width = 0.7,
        axis_tick_width = 0.7, axis_tick_length = 3.5,
        point_size = 2.6, point_alpha = 0.8, jitter_width = 0.13,
        box_line_width = 0.8, violin_line_width = 0.7,
        box_width = 0.36, violin_width = 0.78,
        fill_alpha = 0.16, x_spacing = 1.0,
        axis_text_angle = 0, show_major_grid = FALSE, show_minor_grid = FALSE,
        show_legend = FALSE, show_plot_title = FALSE, show_plot_subtitle = FALSE,
        bold_axis_text = FALSE, bold_axis_title = FALSE,
        legend_position = "right", title_hjust = "0"
      )
      cols <- c("#D55E00", "#0072B2", "#009E73")
    } else if (identical(preset, "minimal")) {
      vals <- list(
        base_font_size = 12, axis_text_size = 10, axis_title_size = 11,
        title_font_size = 14, subtitle_font_size = 10,
        x_axis_line_width = 0.6, y_axis_line_width = 0.6,
        axis_tick_width = 0.6, axis_tick_length = 3,
        point_size = 2.5, point_alpha = 0.72, jitter_width = 0.12,
        box_line_width = 0.7, violin_line_width = 0.6,
        box_width = 0.34, violin_width = 0.74,
        fill_alpha = 0.12, x_spacing = 1.0,
        axis_text_angle = 0, show_major_grid = FALSE, show_minor_grid = FALSE,
        show_legend = FALSE, show_plot_title = FALSE, show_plot_subtitle = FALSE,
        bold_axis_text = FALSE, bold_axis_title = FALSE,
        legend_position = "right", title_hjust = "0.5"
      )
      cols <- c("#666666", "#222222", "#999999")
    } else {
      vals <- list(
        base_font_size = 13, axis_text_size = 11, axis_title_size = 13,
        title_font_size = 16, subtitle_font_size = 12,
        x_axis_line_width = 0.8, y_axis_line_width = 0.8,
        axis_tick_width = 0.8, axis_tick_length = 4,
        point_size = 3, point_alpha = 0.9, jitter_width = 0.16,
        box_line_width = 1.1, violin_line_width = 0.9,
        box_width = 0.42, violin_width = 0.85,
        fill_alpha = 0.22, x_spacing = 1.0,
        axis_text_angle = 45, show_major_grid = TRUE, show_minor_grid = FALSE,
        show_legend = TRUE, show_plot_title = TRUE, show_plot_subtitle = TRUE,
        bold_axis_text = FALSE, bold_axis_title = FALSE,
        legend_position = "right", title_hjust = "0"
      )
      cols <- c("#E64B35", "#0072B5", "#00A087")
    }

    for (nm in names(vals)) {
      val <- vals[[nm]]
      if (nm %in% c("show_major_grid","show_minor_grid","show_legend","show_plot_title",
                    "show_plot_subtitle","bold_axis_text","bold_axis_title")) {
        updateCheckboxInput(session, nm, value = isTRUE(val))
      } else if (nm %in% c("legend_position")) {
        updateSelectInput(session, nm, selected = val)
      } else if (nm %in% c("title_hjust")) {
        updateSelectInput(session, nm, selected = val)
      } else {
        updateSliderInput(session, nm, value = val)
      }
    }

    colourpicker::updateColourInput(session, "groupA_color", value = cols[1])
    colourpicker::updateColourInput(session, "groupB_color", value = cols[2])
    colourpicker::updateColourInput(session, "other_color", value = cols[3])
  }

  observeEvent(input$apply_plot_preset, {
    preset <- ifelse(is.null(input$plot_preset), "prism", input$plot_preset)
    if (!identical(preset, "custom")) {
      apply_single_gene_plot_values(session, preset)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$reset_plot_style, {
    updateSelectInput(session, "plot_preset", selected = "prism")
    apply_single_gene_plot_values(session, "prism")

    updateSelectInput(session, "plot_style", selected = "violin_box")
    updateSelectInput(session, "point_shape", selected = "16")
    updateSliderInput(session, "point_stroke", value = 0.5)
    updateSelectInput(session, "plot_font_family", selected = "")
    updateCheckboxInput(session, "free_y_expand", value = TRUE)
    updateSliderInput(session, "y_expand_ratio", value = 0.08)
    updateCheckboxInput(session, "manual_y_limits", value = FALSE)
    updateNumericInput(session, "y_decimal_digits", value = 1)
    updateCheckboxInput(session, "hide_both_axis_titles", value = FALSE)
    updateCheckboxInput(session, "hide_x_title", value = FALSE)
    updateCheckboxInput(session, "hide_y_title", value = FALSE)
    updateCheckboxInput(session, "hide_x_tick_labels", value = FALSE)
    updateCheckboxInput(session, "hide_y_tick_labels", value = FALSE)
    updateCheckboxInput(session, "show_x_ticks", value = TRUE)
    updateCheckboxInput(session, "show_y_ticks", value = TRUE)
    updateCheckboxInput(session, "show_x_axis_line", value = TRUE)
    updateCheckboxInput(session, "show_y_axis_line", value = TRUE)
    updateTextInput(session, "y_major_interval", value = "")
    updateTextInput(session, "y_minor_interval", value = "")
    updateTextInput(session, "custom_plot_title", value = "")
    updateTextInput(session, "custom_plot_subtitle", value = "")
    updateTextInput(session, "custom_x_title", value = "")
    updateTextInput(session, "custom_y_title", value = "")
    updateSelectInput(session, "plot_p_label_type", selected = "p.format")
    updateSliderInput(session, "plot_p_text_size", value = 4.5)
    updateTextInput(session, "plot_p_y", value = "")
    updateCheckboxInput(session, "show_mean_marker", value = FALSE)
    updateSliderInput(session, "mean_marker_size", value = 4)
    updateSelectInput(session, "mean_marker_shape", selected = "18")
    colourpicker::updateColourInput(session, "mean_marker_color", value = "#000000")
    updateSliderInput(session, "legend_text_size", value = 11)
    updateSliderInput(session, "legend_title_size", value = 12)
    updateSliderInput(session, "plot_margin_top", value = 10)
    updateSliderInput(session, "plot_margin_right", value = 10)
    updateSliderInput(session, "plot_margin_bottom", value = 10)
    updateSliderInput(session, "plot_margin_left", value = 10)
  }, ignoreInit = TRUE)

  plot_reactive <- reactive({
    obj <- data_obj()
    df <- analysis_df()
    if (is.null(obj) || !is.null(obj$error) || is.null(df)) return(NULL)
    req(input$groupA, input$groupB)

    if (!input$show_all) {
      df <- df %>% filter(group %in% c(analysis_group_names()$A, analysis_group_names()$B))
    }

    if (nrow(df) == 0) return(NULL)

    gene <- toupper(trimws(input$gene))
    dobj <- diff_obj()
    method_label <- if (!is.null(dobj) && is.null(dobj$error)) dobj$method else obj$expr_detect$method

    ylab_text <- paste0(gene, " expression used for plotting / ", method_label)

    subtitle_text <- paste0(
      obj$gse_id,
      " | Filter: ", input$filter_col, " = ", input$filter_value,
      " | Group: ", input$group_col
    )

    plot_style <- ifelse(is.null(input$plot_style), "violin_box", input$plot_style)

    groupA_color <- ifelse(is.null(input$groupA_color), "#E64B35", input$groupA_color)
    groupB_color <- ifelse(is.null(input$groupB_color), "#0072B5", input$groupB_color)
    other_color <- ifelse(is.null(input$other_color), "#00A087", input$other_color)

    fill_alpha <- ifelse(is.null(input$fill_alpha), 0.22, input$fill_alpha)
    box_line_width <- ifelse(is.null(input$box_line_width), 1.1, input$box_line_width)
    violin_line_width <- ifelse(is.null(input$violin_line_width), 0.9, input$violin_line_width)
    box_width <- ifelse(is.null(input$box_width), 0.42, input$box_width)
    violin_width <- ifelse(is.null(input$violin_width), 0.85, input$violin_width)
    jitter_width <- ifelse(is.null(input$jitter_width), 0.16, input$jitter_width)
    pt_size <- ifelse(is.null(input$point_size), 3, input$point_size)
    pt_alpha <- ifelse(is.null(input$point_alpha), 0.9, input$point_alpha)
    pt_shape <- suppressWarnings(as.numeric(ifelse(is.null(input$point_shape), 16, input$point_shape)))
    if (!is.finite(pt_shape)) pt_shape <- 16
    pt_stroke <- ifelse(is.null(input$point_stroke), 0.5, input$point_stroke)

    # Keep plot order consistent with the DEG setup:
    # left = Final/Seed Group A (reference), right = Final/Seed Group B (case).
    # If "show all groups" is enabled, other groups are appended after A/B.
    names_run <- analysis_group_names()
    preferred_levels <- c(names_run$A, names_run$B)
    present_groups <- unique(as.character(df$group))
    group_levels <- c(
      preferred_levels[preferred_levels %in% present_groups],
      sort(setdiff(present_groups, preferred_levels))
    )
    groups_now <- group_levels

    # Assign colors:
    # selected/final Group A gets input$groupA_color,
    # selected/final Group B gets input$groupB_color,
    # additional groups get a reusable color series.
    group_cols <- rep(other_color, length(groups_now))
    names(group_cols) <- groups_now

    if (names_run$A %in% names(group_cols)) group_cols[names_run$A] <- groupA_color
    if (names_run$B %in% names(group_cols)) group_cols[names_run$B] <- groupB_color

    extra_groups <- setdiff(groups_now, c(names_run$A, names_run$B))
    if (length(extra_groups) > 0) {
      extra_palette <- c(other_color, "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
      group_cols[extra_groups] <- extra_palette[seq_along(extra_groups)]
    }

    # Adjust visual spacing between groups by using numeric x positions.
    x_spacing <- ifelse(is.null(input$x_spacing), 1.0, input$x_spacing)
    df$group <- factor(df$group, levels = group_levels)
    df$xpos <- as.numeric(df$group) * x_spacing
    x_breaks <- seq_along(group_levels) * x_spacing

    p <- ggplot(df, aes(x = xpos, y = expression, color = group, fill = group))

    # Violin layer: transparent fill + colored border.
    if (plot_style %in% c("violin", "violin_box")) {
      p <- p +
        geom_violin(
          trim = FALSE,
          alpha = fill_alpha,
          linewidth = 0,
          width = violin_width * x_spacing
        ) +
        geom_violin(
          trim = FALSE,
          fill = NA,
          linewidth = violin_line_width,
          width = violin_width * x_spacing
        )
    }

    # Boxplot layer: colored outline + adjustable transparent fill.
    if (plot_style %in% c("box", "violin_box")) {
      p <- p +
        geom_boxplot(
          outlier.shape = NA,
          alpha = fill_alpha,
          linewidth = box_line_width,
          width = box_width * x_spacing
        )
    }

    # Points: adjustable size/transparency/spread.
    p <- p +
      geom_jitter(
        aes(color = group),
        width = jitter_width * x_spacing,
        size = pt_size,
        alpha = pt_alpha,
        shape = pt_shape,
        stroke = pt_stroke
      ) +
      scale_color_manual(values = group_cols) +
      scale_fill_manual(values = group_cols) +
      scale_x_continuous(breaks = x_breaks, labels = group_levels) +
      theme_bw(
        base_size = ifelse(is.null(input$base_font_size), 13, input$base_font_size),
        base_family = ifelse(is.null(input$plot_font_family), "", input$plot_font_family)
      ) +
      theme(
        axis.line.x = if (isTRUE(input$show_x_axis_line)) element_line(
          color = "black",
          linewidth = ifelse(is.null(input$x_axis_line_width), 0.8, input$x_axis_line_width)
        ) else element_blank(),
        axis.line.y = if (isTRUE(input$show_y_axis_line)) element_line(
          color = "black",
          linewidth = ifelse(is.null(input$y_axis_line_width), 0.8, input$y_axis_line_width)
        ) else element_blank(),
        axis.ticks.x = if (isTRUE(input$show_x_ticks)) element_line(
          color = "black",
          linewidth = ifelse(is.null(input$axis_tick_width), 0.8, input$axis_tick_width)
        ) else element_blank(),
        axis.ticks.y = if (isTRUE(input$show_y_ticks)) element_line(
          color = "black",
          linewidth = ifelse(is.null(input$axis_tick_width), 0.8, input$axis_tick_width)
        ) else element_blank(),
        axis.ticks.length = grid::unit(ifelse(is.null(input$axis_tick_length), 4, input$axis_tick_length), "pt"),
        axis.text.x = if (isTRUE(input$hide_x_tick_labels)) element_blank() else element_text(
          angle = ifelse(is.null(input$axis_text_angle), 45, input$axis_text_angle),
          hjust = ifelse(ifelse(is.null(input$axis_text_angle), 45, input$axis_text_angle) == 0, 0.5, 1),
          color = "black",
          size = ifelse(is.null(input$axis_text_size), 11, input$axis_text_size),
          face = ifelse(isTRUE(input$bold_axis_text), "bold", "plain")
        ),
        axis.text.y = if (isTRUE(input$hide_y_tick_labels)) element_blank() else element_text(
          color = "black",
          size = ifelse(is.null(input$axis_text_size), 11, input$axis_text_size),
          face = ifelse(isTRUE(input$bold_axis_text), "bold", "plain")
        ),
        axis.title.x = if (isTRUE(input$hide_both_axis_titles) || isTRUE(input$hide_x_title)) element_blank() else element_text(
          color = "black",
          size = ifelse(is.null(input$axis_title_size), 13, input$axis_title_size),
          face = ifelse(isTRUE(input$bold_axis_title), "bold", "plain")
        ),
        axis.title.y = if (isTRUE(input$hide_both_axis_titles) || isTRUE(input$hide_y_title)) element_blank() else element_text(
          color = "black",
          size = ifelse(is.null(input$axis_title_size), 13, input$axis_title_size),
          face = ifelse(isTRUE(input$bold_axis_title), "bold", "plain")
        ),
        plot.title = element_text(
          face = "bold",
          size = ifelse(is.null(input$title_font_size), 16, input$title_font_size),
          hjust = suppressWarnings(as.numeric(ifelse(is.null(input$title_hjust), 0, input$title_hjust)))
        ),
        plot.subtitle = element_text(
          size = ifelse(is.null(input$subtitle_font_size), 12, input$subtitle_font_size),
          hjust = suppressWarnings(as.numeric(ifelse(is.null(input$title_hjust), 0, input$title_hjust)))
        ),
        legend.position = if (isTRUE(input$show_legend)) ifelse(is.null(input$legend_position), "right", input$legend_position) else "none",
        legend.text = element_text(size = ifelse(is.null(input$legend_text_size), 11, input$legend_text_size)),
        legend.title = element_text(size = ifelse(is.null(input$legend_title_size), 12, input$legend_title_size), face = "bold"),
        panel.grid.major = if (isTRUE(input$show_major_grid)) element_line(linewidth = 0.3, color = "grey85") else element_blank(),
        panel.grid.minor = if (isTRUE(input$show_minor_grid)) element_line(linewidth = 0.2, color = "grey92") else element_blank(),
        plot.margin = margin(
          t = ifelse(is.null(input$plot_margin_top), 10, input$plot_margin_top),
          r = ifelse(is.null(input$plot_margin_right), 10, input$plot_margin_right),
          b = ifelse(is.null(input$plot_margin_bottom), 10, input$plot_margin_bottom),
          l = ifelse(is.null(input$plot_margin_left), 10, input$plot_margin_left)
        )
      ) +
      labs(
        title = if (isTRUE(input$show_plot_title)) {
          if (!is.null(input$custom_plot_title) && nzchar(trimws(input$custom_plot_title))) input$custom_plot_title
          else paste0(gene, " expression among groups")
        } else NULL,
        subtitle = if (isTRUE(input$show_plot_subtitle)) {
          if (!is.null(input$custom_plot_subtitle) && nzchar(trimws(input$custom_plot_subtitle))) input$custom_plot_subtitle
          else subtitle_text
        } else NULL,
        x = if (!is.null(input$custom_x_title) && nzchar(trimws(input$custom_x_title))) input$custom_x_title else "Group",
        y = if (!is.null(input$custom_y_title) && nzchar(trimws(input$custom_y_title))) input$custom_y_title else ylab_text,
        fill = "Group",
        color = "Group"
      )

    # Optional user-defined Y-axis unit spacing. Blank/invalid values preserve ggplot automatic breaks.
    parse_positive_interval <- function(x) {
      if (is.null(x) || length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") return(NULL)
      v <- suppressWarnings(as.numeric(x))
      if (!is.finite(v) || v <= 0) return(NULL)
      v
    }

    y_major_interval <- parse_positive_interval(input$y_major_interval)
    y_minor_interval <- parse_positive_interval(input$y_minor_interval)

    y_scale_args <- list()
    if (isTRUE(input$free_y_expand)) {
      y_scale_args$expand <- expansion(mult = c(0.03, ifelse(is.null(input$y_expand_ratio), 0.08, input$y_expand_ratio)))
    }

    # Optional manual Y limits.
    if (isTRUE(input$manual_y_limits)) {
      ymin_user <- suppressWarnings(as.numeric(input$y_axis_min))
      ymax_user <- suppressWarnings(as.numeric(input$y_axis_max))
      if (is.finite(ymin_user) && is.finite(ymax_user) && ymax_user > ymin_user) {
        y_scale_args$limits <- c(ymin_user, ymax_user)
      }
    }

    # User-controlled decimal places for Y tick labels.
    y_digits <- suppressWarnings(as.integer(ifelse(is.null(input$y_decimal_digits), 1, input$y_decimal_digits)))
    if (!is.finite(y_digits) || y_digits < 0) y_digits <- 1
    y_digits <- min(y_digits, 6)
    y_scale_args$labels <- function(x) formatC(x, format = "f", digits = y_digits)

    if (!is.null(y_major_interval)) {
      y_range <- range(df$expression, na.rm = TRUE)
      if (all(is.finite(y_range))) {
        lower <- floor(y_range[1] / y_major_interval) * y_major_interval
        upper <- ceiling(y_range[2] / y_major_interval) * y_major_interval
        y_scale_args$breaks <- seq(lower, upper, by = y_major_interval)
      }
    }

    if (!is.null(y_minor_interval)) {
      y_range <- range(df$expression, na.rm = TRUE)
      if (all(is.finite(y_range))) {
        lower_minor <- floor(y_range[1] / y_minor_interval) * y_minor_interval
        upper_minor <- ceiling(y_range[2] / y_minor_interval) * y_minor_interval
        y_scale_args$minor_breaks <- seq(lower_minor, upper_minor, by = y_minor_interval)
      }
    }

    if (length(y_scale_args) > 0) {
      p <- p + do.call(scale_y_continuous, y_scale_args)
    }

    # Optional group mean marker.
    if (isTRUE(input$show_mean_marker)) {
      mean_shape <- suppressWarnings(as.numeric(ifelse(is.null(input$mean_marker_shape), 18, input$mean_marker_shape)))
      if (!is.finite(mean_shape)) mean_shape <- 18
      p <- p + stat_summary(
        fun = mean,
        geom = "point",
        shape = mean_shape,
        size = ifelse(is.null(input$mean_marker_size), 4, input$mean_marker_size),
        color = ifelse(is.null(input$mean_marker_color), "#000000", input$mean_marker_color),
        fill = ifelse(is.null(input$mean_marker_color), "#000000", input$mean_marker_color),
        show.legend = FALSE
      )
    }

    if (isTRUE(input$show_plot_p) && !input$show_all && length(unique(df$group)) == 2) {
      p_label <- ifelse(is.null(input$plot_p_label_type), "p.format", input$plot_p_label_type)
      p_size <- ifelse(is.null(input$plot_p_text_size), 4.5, input$plot_p_text_size)

      p_y <- NULL
      if (!is.null(input$plot_p_y) && nzchar(trimws(as.character(input$plot_p_y)))) {
        p_y_try <- suppressWarnings(as.numeric(input$plot_p_y))
        if (is.finite(p_y_try)) p_y <- p_y_try
      }

      if (is.null(p_y)) {
        p <- p + stat_compare_means(method = "wilcox.test", label = p_label, size = p_size)
      } else {
        p <- p + stat_compare_means(method = "wilcox.test", label = p_label, size = p_size, label.y = p_y)
      }
    }

    p
  })
  output$plot_status_ui <- renderUI({
    w <- ifelse(is.null(input$export_width), 10, input$export_width)
    h <- ifelse(is.null(input$export_height), 7, input$export_height)
    style_name <- if (is.null(input$plot_preset)) "custom" else input$plot_preset
    tagList(
      div(
        class = "plot-status-bar",
        span(class = "plot-status-chip", paste0(w, " x ", h, " inch")),
        span(class = "plot-status-chip", "PNG 600 dpi"),
        span(class = "plot-status-chip", paste0("Preset: ", style_name))
      )
    )
  })

  output$plot_ui <- renderUI({
    h <- ifelse(is.null(input$plot_height_px), 650, input$plot_height_px)
    plotOutput("plot", height = paste0(h, "px"))
  })

  output$plot <- renderPlot(plot_reactive())




  # -------------------------
  # Correlation module
  # -------------------------
  corr_meta_selected <- reactive({
    meta2 <- selected_sample_meta()
    if (is.null(meta2)) return(NULL)

    use_ab <- isTRUE(input$corr_use_selected_groups) ||
      isTRUE(input$corr_all_current_groups)

    if (use_ab && !is.null(input$groupA) && !is.null(input$groupB)) {
      meta2 <- meta2 %>% filter(group %in% c(analysis_group_names()$A, analysis_group_names()$B))
    }

    meta2
  })

  corr_expr_mat <- reactive({
    obj <- data_obj()
    dobj <- diff_obj()
    get_plot_expr_matrix_for_modules(obj, dobj)
  })

  corr_gene_gene_obj <- reactive({
    mat <- corr_expr_mat()
    meta2 <- selected_sample_meta()
    if (is.null(meta2)) return(NULL)

    if (isTRUE(input$corr_use_selected_groups) && !is.null(input$groupA) && !is.null(input$groupB)) {
      meta2 <- meta2 %>% filter(group %in% c(analysis_group_names()$A, analysis_group_names()$B))
    }

    make_gene_gene_correlation(
      mat = mat,
      meta2 = meta2,
      gene_a = input$corr_gene_a,
      gene_b = input$corr_gene_b,
      method = input$corr_method
    )
  })

  output$corr_gene_gene_stats <- renderTable({
    obj <- corr_gene_gene_obj()
    if (is.null(obj)) return(NULL)

    data.frame(
      Gene_A = obj$gene_a,
      Gene_B = obj$gene_b,
      Method = obj$method,
      Correlation = round(obj$r, 4),
      P_value = signif(obj$p, 4),
      N = nrow(obj$df),
      stringsAsFactors = FALSE
    )
  })

  output$corr_gene_gene_plot <- renderPlot({
    obj <- corr_gene_gene_obj()
    if (is.null(obj)) return(NULL)

    p <- ggplot(obj$df, aes(x = Gene_A, y = Gene_B, color = group)) +
      geom_point(
        size = ifelse(is.null(input$corr_point_size), 3, input$corr_point_size),
        alpha = ifelse(is.null(input$corr_point_alpha), 0.85, input$corr_point_alpha)
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste0(obj$gene_a, " vs ", obj$gene_b),
        subtitle = paste0(obj$method, " r = ", round(obj$r, 3), ", P = ", signif(obj$p, 3)),
        x = paste0(obj$gene_a, " expression"),
        y = paste0(obj$gene_b, " expression"),
        color = "Group"
      )

    if (isTRUE(input$corr_show_lm)) {
      p <- p + geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed")
    }

    if (isTRUE(input$corr_show_labels)) {
      p <- p + geom_text(aes(label = sample), vjust = -0.8, size = 3, show.legend = FALSE)
    }

    print(p)
  })

  corr_all_obj <- reactive({
    mat <- corr_expr_mat()
    meta2 <- selected_sample_meta()
    if (is.null(meta2)) return(NULL)

    if (isTRUE(input$corr_all_current_groups) && !is.null(input$groupA) && !is.null(input$groupB)) {
      meta2 <- meta2 %>% filter(group %in% c(analysis_group_names()$A, analysis_group_names()$B))
    }

    make_gene_all_correlation(
      mat = mat,
      meta2 = meta2,
      gene = input$corr_all_gene,
      method = input$corr_all_method,
      top_n = ifelse(is.null(input$corr_all_top_n), 30, input$corr_all_top_n)
    )
  })

  output$corr_top_positive <- renderTable({
    obj <- corr_all_obj()
    if (is.null(obj)) return(NULL)
    obj$top_positive %>%
      mutate(correlation = round(correlation, 4), pvalue = signif(pvalue, 3), padj = signif(padj, 3))
  })

  output$corr_top_negative <- renderTable({
    obj <- corr_all_obj()
    if (is.null(obj)) return(NULL)
    obj$top_negative %>%
      mutate(correlation = round(correlation, 4), pvalue = signif(pvalue, 3), padj = signif(padj, 3))
  })

  output$corr_family_heatmap <- renderPlot({
    mat <- corr_expr_mat()
    meta2 <- selected_sample_meta()
    cm <- make_gene_family_correlation(
      mat = mat,
      meta2 = meta2,
      genes = input$corr_family_genes,
      method = input$corr_family_method
    )

    if (is.null(cm)) {
      plot.new()
      text(0.5, 0.5, "相关矩阵不可用：请检查基因名或样本数量。", cex = 1.1)
      return()
    }

    pheatmap::pheatmap(
      cm,
      cluster_rows = isTRUE(input$corr_family_cluster),
      cluster_cols = isTRUE(input$corr_family_cluster),
      display_numbers = TRUE,
      number_format = "%.2f",
      color = colorRampPalette(c("#0072B5", "white", "#E64B35"))(101),
      main = paste0("Gene family correlation matrix (", input$corr_family_method, ")")
    )
  })

  output$downloadCorrAll <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$corr_all_gene)), "_all_gene_correlation.csv"),
    content = function(file) {
      obj <- corr_all_obj()
      if (is.null(obj)) {
        write.csv(data.frame(message = "Correlation result is not available."), file, row.names = FALSE)
      } else {
        write.csv(obj$full, file, row.names = FALSE)
      }
    }
  )

  output$downloadCorrFamily <- downloadHandler(
    filename = function() paste0(input$gse, "_gene_family_correlation_matrix.csv"),
    content = function(file) {
      mat <- corr_expr_mat()
      meta2 <- selected_sample_meta()
      cm <- make_gene_family_correlation(
        mat = mat,
        meta2 = meta2,
        genes = input$corr_family_genes,
        method = input$corr_family_method
      )

      if (is.null(cm)) {
        write.csv(data.frame(message = "Correlation matrix is not available."), file, row.names = FALSE)
      } else {
        write.csv(cm, file)
      }
    }
  )

  # -------------------------
  # Enrichment module
  # -------------------------
  enrich_obj <- reactiveVal(NULL)

  observeEvent(input$run_enrichment, {
    dobj <- diff_obj()

    genes <- get_deg_genes_for_enrichment(
      dobj = dobj,
      direction = input$enrich_direction,
      padj_cutoff = input$enrich_padj_cutoff,
      logfc_cutoff = input$enrich_logfc_cutoff
    )

    if (length(genes) < 3) {
      enrich_obj(list(error = "用于富集分析的基因少于3个。请放宽padj/log2FC阈值。", genes = genes))
      return()
    }

    go_res <- run_go_enrichment_safe(
      gene_symbols = genes,
      ont = input$enrich_go_ont,
      p_cutoff = input$enrich_p_cutoff,
      q_cutoff = input$enrich_q_cutoff
    )

    kegg_res <- run_kegg_enrichment_safe(
      gene_symbols = genes,
      p_cutoff = input$enrich_p_cutoff,
      q_cutoff = input$enrich_q_cutoff
    )

    enrich_obj(list(
      genes = genes,
      go = go_res,
      kegg = kegg_res,
      error = NULL
    ))
  })

  output$enrich_gene_count <- renderText({
    dobj <- diff_obj()
    genes <- get_deg_genes_for_enrichment(
      dobj = dobj,
      direction = input$enrich_direction,
      padj_cutoff = input$enrich_padj_cutoff,
      logfc_cutoff = input$enrich_logfc_cutoff
    )

    paste0("当前阈值下可用于富集的基因数：", length(genes))
  })

  output$go_dotplot <- renderPlot({
    eo <- enrich_obj()
    if (is.null(eo)) {
      plot.new()
      text(0.5, 0.5, "点击“运行GO/KEGG富集”后显示结果。", cex = 1.1)
      return()
    }

    if (!is.null(eo$error)) {
      plot.new()
      text(0.5, 0.5, eo$error, cex = 1.1)
      return()
    }

    if (is.null(eo$go) || nrow(as.data.frame(eo$go)) == 0) {
      plot.new()
      text(0.5, 0.5, "没有显著GO富集结果。", cex = 1.1)
      return()
    }

    print(enrichplot::dotplot(eo$go, showCategory = input$enrich_show_n) +
            ggplot2::ggtitle(paste0("GO ", input$enrich_go_ont, " enrichment")))
  })

  output$kegg_dotplot <- renderPlot({
    eo <- enrich_obj()
    if (is.null(eo)) {
      plot.new()
      text(0.5, 0.5, "点击“运行GO/KEGG富集”后显示结果。", cex = 1.1)
      return()
    }

    if (!is.null(eo$error)) {
      plot.new()
      text(0.5, 0.5, eo$error, cex = 1.1)
      return()
    }

    if (is.null(eo$kegg) || nrow(as.data.frame(eo$kegg)) == 0) {
      plot.new()
      text(0.5, 0.5, "没有显著KEGG富集结果。", cex = 1.1)
      return()
    }

    print(enrichplot::dotplot(eo$kegg, showCategory = input$enrich_show_n) +
            ggplot2::ggtitle("KEGG enrichment"))
  })

  output$go_table <- renderTable({
    eo <- enrich_obj()
    if (is.null(eo) || !is.null(eo$error)) return(NULL)
    df <- enrich_result_table(eo$go)
    if (is.null(df)) return(NULL)
    head(df, 50)
  })

  output$kegg_table <- renderTable({
    eo <- enrich_obj()
    if (is.null(eo) || !is.null(eo$error)) return(NULL)
    df <- enrich_result_table(eo$kegg)
    if (is.null(df)) return(NULL)
    head(df, 50)
  })

  output$downloadGO <- downloadHandler(
    filename = function() paste0(input$gse, "_GO_", input$enrich_go_ont, "_enrichment.csv"),
    content = function(file) {
      eo <- enrich_obj()
      df <- if (is.null(eo) || !is.null(eo$error)) NULL else enrich_result_table(eo$go)
      if (is.null(df)) df <- data.frame(message = "No GO enrichment result.")
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$downloadKEGG <- downloadHandler(
    filename = function() paste0(input$gse, "_KEGG_enrichment.csv"),
    content = function(file) {
      eo <- enrich_obj()
      df <- if (is.null(eo) || !is.null(eo$error)) NULL else enrich_result_table(eo$kegg)
      if (is.null(df)) df <- data.frame(message = "No KEGG enrichment result.")
      write.csv(df, file, row.names = FALSE)
    }
  )


  # -------------------------
  # Dashboard / DEG summary / PCA
  # -------------------------
  output$dash_samples <- renderText({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return("NA")
    ncol(obj$raw) - 1
  })

  output$dash_genes <- renderText({
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return("NA")
    nrow(obj$raw)
  })

  output$dash_method <- renderText({
    dobj <- diff_obj()
    if (is.null(dobj) || !is.null(dobj$error)) return("NA")
    dobj$method
  })

  output$dash_match <- renderText({
    dobj <- diff_obj()
    if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$coldata)) {
      return(paste0(nrow(dobj$coldata), "/", nrow(dobj$coldata), " used in DEG"))
    }
    obj <- data_obj()
    if (is.null(obj) || !is.null(obj$error)) return("NA")
    paste0(obj$sample_mapping_matched_n, "/", ncol(obj$raw) - 1)
  })

  output$dash_seed_groups <- renderText({
    if (isTRUE(input$enable_group_merge)) {
      paste0("Seed A: ", input$groupA, "  |  Seed B: ", input$groupB,
             "  ->  Final A/B: ", analysis_group_names()$A, " / ", analysis_group_names()$B)
    } else {
      paste0("Manual merge OFF. Final groups equal Seed groups: ", input$groupA, " / ", input$groupB)
    }
  })

  output$dash_group_sample_status <- renderTable({
    names_run <- analysis_group_names()
    if (is.null(names_run$A) || is.null(names_run$B)) return(NULL)

    obj <- data_obj()
    df <- selected_sample_meta()
    metadata_A <- if (!is.null(df)) sum(df$group == names_run$A, na.rm = TRUE) else NA_integer_
    metadata_B <- if (!is.null(df)) sum(df$group == names_run$B, na.rm = TRUE) else NA_integer_

    expr_A <- NA_integer_
    expr_B <- NA_integer_
    if (!is.null(obj) && is.null(obj$error) && !is.null(df)) {
      expr_meta <- get_expression_available_meta(df, obj$raw)
      expr_A <- sum(expr_meta$group == names_run$A, na.rm = TRUE)
      expr_B <- sum(expr_meta$group == names_run$B, na.rm = TRUE)
    }

    dobj <- diff_obj()
    used_A <- NA_integer_
    used_B <- NA_integer_
    dropped_A <- NA_integer_
    dropped_B <- NA_integer_

    if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$coldata)) {
      used_A <- sum(dobj$coldata$group == names_run$A, na.rm = TRUE)
      used_B <- sum(dobj$coldata$group == names_run$B, na.rm = TRUE)
      dropped_A <- ifelse(is.na(metadata_A), NA_integer_, metadata_A - used_A)
      dropped_B <- ifelse(is.na(metadata_B), NA_integer_, metadata_B - used_B)
    } else {
      dropped_A <- ifelse(is.na(metadata_A) || is.na(expr_A), NA_integer_, metadata_A - expr_A)
      dropped_B <- ifelse(is.na(metadata_B) || is.na(expr_B), NA_integer_, metadata_B - expr_B)
    }

    data.frame(
      Group = c(names_run$A, names_run$B),
      Metadata_selected = c(metadata_A, metadata_B),
      Expression_available = c(expr_A, expr_B),
      Used_in_DEG = c(used_A, used_B),
      Dropped_by_expression_matrix = c(dropped_A, dropped_B),
      check.names = FALSE
    )
  })

  output$dash_groupA <- renderText({
    names_run <- analysis_group_names()
    if (is.null(names_run$A)) return("NA")
    dobj <- diff_obj()
    if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$coldata)) {
      return(paste0(names_run$A, " (", sum(dobj$coldata$group == names_run$A), ")"))
    }
    df <- selected_sample_meta()
    obj <- data_obj()
    if (is.null(df) || is.null(obj) || !is.null(obj$error)) return("NA")
    expr_meta <- get_expression_available_meta(df, obj$raw)
    paste0(names_run$A, " (expr ", sum(expr_meta$group == names_run$A), " / meta ", sum(df$group == names_run$A), ")")
  })

  output$dash_groupB <- renderText({
    names_run <- analysis_group_names()
    if (is.null(names_run$B)) return("NA")
    dobj <- diff_obj()
    if (!is.null(dobj) && is.null(dobj$error) && !is.null(dobj$coldata)) {
      return(paste0(names_run$B, " (", sum(dobj$coldata$group == names_run$B), ")"))
    }
    df <- selected_sample_meta()
    obj <- data_obj()
    if (is.null(df) || is.null(obj) || !is.null(obj$error)) return("NA")
    expr_meta <- get_expression_available_meta(df, obj$raw)
    paste0(names_run$B, " (expr ", sum(expr_meta$group == names_run$B), " / meta ", sum(df$group == names_run$B), ")")
  })

  output$dash_up <- renderText({
    dobj <- diff_obj()
    if (is.null(dobj) || !is.null(dobj$error)) return("NA")
    tbl <- make_deg_summary(dobj, input$deg_table_padj, input$deg_table_logfc)
    tbl$Value[tbl$Metric == "Upregulated"]
  })

  output$dash_down <- renderText({
    dobj <- diff_obj()
    if (is.null(dobj) || !is.null(dobj$error)) return("NA")
    tbl <- make_deg_summary(dobj, input$deg_table_padj, input$deg_table_logfc)
    tbl$Value[tbl$Metric == "Downregulated"]
  })

  output$deg_summary_table <- renderTable({
    dobj <- diff_obj()
    make_deg_summary(
      dobj,
      padj_cutoff = ifelse(is.null(input$deg_table_padj), 0.05, input$deg_table_padj),
      logfc_cutoff = ifelse(is.null(input$deg_table_logfc), 1, input$deg_table_logfc)
    )
  })

  output$deg_search_table <- renderTable({
    dobj <- diff_obj()
    out <- search_deg_table(dobj, input$deg_search)
    if (is.null(out) || nrow(out) == 0) return(NULL)

    out <- add_logfc_arrow_html(out)

    # Rename for display.
    if ("log2FC_direction" %in% colnames(out)) {
      colnames(out)[colnames(out) == "log2FC_direction"] <- "log2FC"
    }

    out
  }, digits = 4, sanitize.text.function = function(x) x)

  output$top_up_table <- renderTable({
    out <- get_top_deg_table(
      diff_obj(),
      direction = "up",
      top_n = 20,
      padj_cutoff = ifelse(is.null(input$deg_table_padj), 0.05, input$deg_table_padj),
      logfc_cutoff = ifelse(is.null(input$deg_table_logfc), 1, input$deg_table_logfc)
    )

    if (is.null(out) || nrow(out) == 0) return(NULL)
    out <- add_logfc_arrow_html(out)
    if ("log2FC_direction" %in% colnames(out)) {
      colnames(out)[colnames(out) == "log2FC_direction"] <- "log2FC"
    }
    out
  }, digits = 4, sanitize.text.function = function(x) x)

  output$top_down_table <- renderTable({
    out <- get_top_deg_table(
      diff_obj(),
      direction = "down",
      top_n = 20,
      padj_cutoff = ifelse(is.null(input$deg_table_padj), 0.05, input$deg_table_padj),
      logfc_cutoff = ifelse(is.null(input$deg_table_logfc), 1, input$deg_table_logfc)
    )

    if (is.null(out) || nrow(out) == 0) return(NULL)
    out <- add_logfc_arrow_html(out)
    if ("log2FC_direction" %in% colnames(out)) {
      colnames(out)[colnames(out) == "log2FC_direction"] <- "log2FC"
    }
    out
  }, digits = 4, sanitize.text.function = function(x) x)

  build_pca_plot <- function() {
    dobj <- diff_obj()
    pca_obj <- prepare_pca_data(
      dobj,
      top_var_genes = ifelse(is.null(input$pca_top_var_genes), 1000, input$pca_top_var_genes)
    )

    if (is.null(pca_obj)) return(NULL)

    df <- pca_obj$pca_df
    var_exp <- pca_obj$var_exp

    groups_now <- sort(unique(as.character(df$group)))
    group_cols <- rep(ifelse(is.null(input$pca_other_color), "#00A087", input$pca_other_color), length(groups_now))
    names(group_cols) <- groups_now
    if (!is.null(dobj)) {
      if (dobj$groupA %in% names(group_cols)) group_cols[dobj$groupA] <- ifelse(is.null(input$pca_groupA_color), "#E64B35", input$pca_groupA_color)
      if (dobj$groupB %in% names(group_cols)) group_cols[dobj$groupB] <- ifelse(is.null(input$pca_groupB_color), "#0072B5", input$pca_groupB_color)
    }

    p <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
      geom_point(
        size = ifelse(is.null(input$pca_point_size), 3.5, input$pca_point_size),
        alpha = ifelse(is.null(input$pca_point_alpha), 0.85, input$pca_point_alpha)
      ) +
      scale_color_manual(values = group_cols) +
      theme_bw(base_size = ifelse(is.null(input$pca_base_font_size), 13, input$pca_base_font_size)) +
      theme(
        panel.grid.major = if (isTRUE(input$pca_show_grid)) element_line(color = "grey85", linewidth = 0.3) else element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      ) +
      labs(
        title = paste0("PCA based on top variable genes (n=", pca_obj$top_var_genes, ")"),
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
        color = "Group"
      )

    if (isTRUE(input$pca_show_ellipse) && length(unique(df$group)) >= 2) {
      p <- p + stat_ellipse(aes(group = group), linetype = "dashed", linewidth = 0.7, show.legend = FALSE)
    }

    if (isTRUE(input$pca_show_labels)) {
      p <- p + geom_text(aes(label = sample), vjust = -0.8, size = 3, show.legend = FALSE)
    }

    p
  }

  output$pca_ui <- renderUI({
    h <- ifelse(is.null(input$pca_height_px), 650, input$pca_height_px)
    plotOutput("pca_plot", height = paste0(h, "px"))
  })

  output$pca_plot <- renderPlot({
    p <- build_pca_plot()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "PCA数据不足或差异分析尚未完成。", cex = 1.1)
    } else {
      print(p)
    }
  })

  build_volcano_plot <- function() {
    dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
    if (is.null(dobj) || !is.null(dobj$error)) return(NULL)

    gene <- toupper(trimws(input$gene))

    lfc_cut <- ifelse(is.null(input$volcano_logfc_cutoff), 1, input$volcano_logfc_cutoff)
    p_cut <- ifelse(is.null(input$volcano_padj_cutoff), 0.05, input$volcano_padj_cutoff)
    p_type <- ifelse(is.null(input$volcano_p_type), "padj", input$volcano_p_type)

    y_p <- if (p_type == "pvalue" && "pvalue" %in% colnames(dobj$res_df)) {
      dobj$res_df$pvalue
    } else {
      dobj$res_df$padj
    }

    vol <- dobj$res_df %>%
      mutate(
        p_for_plot = y_p,
        p_plot = ifelse(is.na(p_for_plot), 1, pmax(p_for_plot, 1e-300)),
        neglog10p = -log10(p_plot),
        sig = case_when(
          !is.na(p_for_plot) & p_for_plot < p_cut & log2FoldChange >= lfc_cut ~ "Up",
          !is.na(p_for_plot) & p_for_plot < p_cut & log2FoldChange <= -lfc_cut ~ "Down",
          TRUE ~ "NS"
        ),
        target = ifelse(toupper(Symbol) == gene, "Target", "Other")
      )

    volcano_cols <- c(
      Up = ifelse(is.null(input$volcano_up_color), "#E64B35", input$volcano_up_color),
      Down = ifelse(is.null(input$volcano_down_color), "#0072B5", input$volcano_down_color),
      NS = ifelse(is.null(input$volcano_ns_color), "#BDBDBD", input$volcano_ns_color)
    )

    p <- ggplot(vol, aes(x = log2FoldChange, y = neglog10p, color = sig)) +
      geom_point(
        size = ifelse(is.null(input$volcano_point_size), 1.2, input$volcano_point_size),
        alpha = ifelse(is.null(input$volcano_point_alpha), 0.7, input$volcano_point_alpha)
      ) +
      geom_point(
        data = vol %>% filter(target == "Target"),
        aes(x = log2FoldChange, y = neglog10p),
        color = ifelse(is.null(input$volcano_target_color), "#000000", input$volcano_target_color),
        size = ifelse(is.null(input$volcano_target_size), 3.5, input$volcano_target_size)
      ) +
      scale_color_manual(values = volcano_cols) +
      theme_bw(base_size = ifelse(is.null(input$volcano_base_font_size), 13, input$volcano_base_font_size)) +
      theme(
        legend.position = ifelse(isTRUE(input$volcano_show_legend), "right", "none"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = paste0("Volcano plot: ", dobj$groupB, " vs ", dobj$groupA, " (", dobj$method, ")"),
        x = paste0("log2 Fold Change (", dobj$groupB, " / ", dobj$groupA, ")"),
        y = ifelse(p_type == "pvalue", "-log10 P value", "-log10 adjusted P value"),
        color = "Significance"
      )

    if (isTRUE(input$volcano_show_cutoff_lines)) {
      p <- p +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
        geom_hline(yintercept = -log10(p_cut), linetype = "dashed")
    }

    if (isTRUE(input$volcano_show_target_label) || isTRUE(input$volcano_export_label)) {
      p <- p +
        geom_text(
          data = vol %>% filter(target == "Target"),
          aes(label = Symbol),
          color = ifelse(is.null(input$volcano_target_color), "#000000", input$volcano_target_color),
          vjust = -1,
          size = 4
        )
    }

    p
  }

  get_volcano_sig_table <- function(direction = "all") {
    dobj <- diff_obj()
    if (is.null(dobj) || !is.null(dobj$error)) return(data.frame())

    lfc_cut <- ifelse(is.null(input$volcano_logfc_cutoff), 1, input$volcano_logfc_cutoff)
    p_cut <- ifelse(is.null(input$volcano_padj_cutoff), 0.05, input$volcano_padj_cutoff)
    p_type <- ifelse(is.null(input$volcano_p_type), "padj", input$volcano_p_type)

    res <- dobj$res_df
    pval <- if (p_type == "pvalue" && "pvalue" %in% colnames(res)) res$pvalue else res$padj

    out <- res %>%
      mutate(plot_p = pval) %>%
      filter(!is.na(plot_p), plot_p < p_cut, !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_cut)

    if (direction == "up") {
      out <- out %>% filter(log2FoldChange >= lfc_cut)
    } else if (direction == "down") {
      out <- out %>% filter(log2FoldChange <= -lfc_cut)
    }

    out
  }

  get_heatmap_object <- function() {
    dobj <- if (isTRUE(deg_has_run())) diff_obj() else NULL
    if (is.null(dobj) || !is.null(dobj$error)) return(NULL)
    if (is.null(dobj$plot_expr_mat) || is.null(dobj$coldata)) return(NULL)

    prepare_heatmap_data(
      diff_res = dobj$res_df,
      expr_mat = dobj$plot_expr_mat,
      coldata = dobj$coldata,
      top_n = ifelse(is.null(input$heatmap_top_n), 50, input$heatmap_top_n),
      padj_cutoff = ifelse(is.null(input$heatmap_padj_cutoff), 0.05, input$heatmap_padj_cutoff),
      logfc_cutoff = ifelse(is.null(input$heatmap_logfc_cutoff), 0, input$heatmap_logfc_cutoff),
      rank_by = ifelse(is.null(input$heatmap_rank_by), "padj", input$heatmap_rank_by)
    )
  }

  draw_heatmap <- function(silent = FALSE) {
    dobj <- diff_obj()
    hm <- get_heatmap_object()
    if (is.null(hm)) {
      plot.new()
      text(0.5, 0.5, "没有足够基因用于热图。请放宽Padj/log2FC阈值或增加Top基因数量。", cex = 1.1)
      return(NULL)
    }

    group_levels <- unique(as.character(hm$annotation_col$Group))
    ann_colors <- list(Group = setNames(rep("#00A087", length(group_levels)), group_levels))
    if (!is.null(dobj)) {
      if (dobj$groupA %in% group_levels) ann_colors$Group[dobj$groupA] <- ifelse(is.null(input$heatmap_groupA_color), "#E64B35", input$heatmap_groupA_color)
      if (dobj$groupB %in% group_levels) ann_colors$Group[dobj$groupB] <- ifelse(is.null(input$heatmap_groupB_color), "#0072B5", input$heatmap_groupB_color)
    }

    border_col <- ifelse(ifelse(is.null(input$heatmap_border_width), 0, input$heatmap_border_width) <= 0, NA, "grey90")

    pheatmap::pheatmap(
      hm$mat,
      scale = ifelse(is.null(input$heatmap_scale), "row", input$heatmap_scale),
      color = make_heatmap_palette(ifelse(is.null(input$heatmap_palette), "blue_white_red", input$heatmap_palette), 101),
      cluster_rows = isTRUE(input$heatmap_cluster_rows),
      cluster_cols = isTRUE(input$heatmap_cluster_cols),
      show_rownames = isTRUE(input$heatmap_show_rownames),
      show_colnames = isTRUE(input$heatmap_show_colnames),
      annotation_col = hm$annotation_col,
      annotation_colors = ann_colors,
      fontsize = ifelse(is.null(input$heatmap_fontsize), 9, input$heatmap_fontsize),
      fontsize_row = ifelse(is.null(input$heatmap_fontsize_row), 7, input$heatmap_fontsize_row),
      cellwidth = ifelse(is.null(input$heatmap_cellwidth), NA, input$heatmap_cellwidth),
      cellheight = ifelse(is.null(input$heatmap_cellheight), NA, input$heatmap_cellheight),
      border_color = border_col,
      main = if (!is.null(dobj)) paste0("Top DEGs heatmap: ", dobj$groupB, " vs ", dobj$groupA) else "Top DEGs heatmap",
      silent = silent
    )
  }

  output$volcano_ui <- renderUI({
    h <- ifelse(is.null(input$volcano_height_px), 650, input$volcano_height_px)
    plotOutput("volcano", height = paste0(h, "px"))
  })

  output$volcano <- renderPlot({
    p <- build_volcano_plot()
    if (is.null(p)) return(NULL)
    print(p)
  })

  output$heatmap_ui <- renderUI({
    h <- ifelse(is.null(input$heatmap_height_px), 700, input$heatmap_height_px)
    plotOutput("heatmap", height = paste0(h, "px"))
  })

  output$heatmap <- renderPlot({
    draw_heatmap(silent = FALSE)
  })

  output$downloadDEG <- downloadHandler(
    filename = function() paste0(input$gse, "_autoDiff_all_results_", input$groupB, "_vs_", input$groupA, ".csv"),
    content = function(file) {
      dobj <- diff_obj()
      if (is.null(dobj) || !is.null(dobj$error)) {
        write.csv(data.frame(message = "Differential analysis result is not available."), file, row.names = FALSE)
      } else {
        write.csv(dobj$res_df, file, row.names = FALSE)
      }
    }
  )
  output$downloadSigDEG <- downloadHandler(
    filename = function() paste0(input$gse, "_autoDiff_DEGs_padj0.05_log2FC1_", input$groupB, "_vs_", input$groupA, ".csv"),
    content = function(file) {
      dobj <- diff_obj()
      if (is.null(dobj) || !is.null(dobj$error)) {
        write.csv(data.frame(message = "Differential analysis result is not available."), file, row.names = FALSE)
      } else {
        sig <- dobj$res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
        write.csv(sig, file, row.names = FALSE)
      }
    }
  )
  output$downloadPDF <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$gene)), "_group_plot.pdf"),
    content = function(file) ggsave(file, plot = plot_reactive(), width = input$export_width, height = input$export_height, device = "pdf")
  )
  output$downloadPNG <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$gene)), "_group_plot.png"),
    content = function(file) ggsave(file, plot = plot_reactive(), width = input$export_width, height = input$export_height, dpi = 600, device = "png")
  )

  output$downloadPlotPNG <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$gene)), "_group_plot.png"),
    content = function(file) {
      ggsave(
        file,
        plot = plot_reactive(),
        width = input$export_width,
        height = input$export_height,
        dpi = 600,
        device = "png"
      )
    }
  )

  output$downloadPlotPDF <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$gene)), "_group_plot.pdf"),
    content = function(file) {
      ggsave(
        file,
        plot = plot_reactive(),
        width = input$export_width,
        height = input$export_height,
        device = "pdf"
      )
    }
  )

  output$downloadPlotSVG <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$gene)), "_group_plot.svg"),
    content = function(file) {
      svglite::svglite(file, width = input$export_width, height = input$export_height)
      print(plot_reactive())
      dev.off()
    }
  )

  output$downloadSingleGeneStats <- downloadHandler(
    filename = function() paste0(input$gse, "_", toupper(trimws(input$gene)), "_single_gene_stats.csv"),
    content = function(file) {
      obj <- data_obj()
      df <- analysis_df()

      if (is.null(obj) || is.null(df)) {
        write.csv(data.frame(message = "Single-gene statistics are not available."), file, row.names = FALSE)
      } else {
        if (!input$show_all) {
          df <- df %>% filter(group %in% c(analysis_group_names()$A, analysis_group_names()$B))
        }

        stats_obj <- make_single_gene_combined_stats(
          df = df,
          groupA = input$groupA,
          groupB = input$groupB,
          gene = toupper(trimws(input$gene)),
          gse_id = obj$gse_id
        )

        group_stats <- stats_obj$group_stats
        comp_stats <- stats_obj$comparison_stats

        # Write both sections into one CSV-like file.
        con <- file(file, open = "w")
        writeLines("Group descriptive statistics", con)
        close(con)
        write.table(group_stats, file = file, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)

        con <- file(file, open = "a")
        writeLines("", con)
        writeLines("Two-group comparison statistics", con)
        close(con)
        write.table(comp_stats, file = file, sep = ",", row.names = FALSE, col.names = TRUE, append = TRUE)
      }
    }
  )

  output$downloadVolcanoPNG <- downloadHandler(
    filename = function() paste0(input$gse, "_volcano_", input$groupB, "_vs_", input$groupA, ".png"),
    content = function(file) {
      ggsave(
        file,
        plot = build_volcano_plot(),
        width = ifelse(is.null(input$volcano_export_width), 10, input$volcano_export_width),
        height = ifelse(is.null(input$volcano_export_height), 7, input$volcano_export_height),
        dpi = as.numeric(ifelse(is.null(input$volcano_export_dpi), 600, input$volcano_export_dpi)),
        device = "png"
      )
    }
  )

  output$downloadVolcanoPDF <- downloadHandler(
    filename = function() paste0(input$gse, "_volcano_", input$groupB, "_vs_", input$groupA, ".pdf"),
    content = function(file) {
      ggsave(
        file,
        plot = build_volcano_plot(),
        width = ifelse(is.null(input$volcano_export_width), 10, input$volcano_export_width),
        height = ifelse(is.null(input$volcano_export_height), 7, input$volcano_export_height),
        device = "pdf"
      )
    }
  )

  output$downloadVolcanoSVG <- downloadHandler(
    filename = function() paste0(input$gse, "_volcano_", input$groupB, "_vs_", input$groupA, ".svg"),
    content = function(file) {
      svglite::svglite(
        file,
        width = ifelse(is.null(input$volcano_export_width), 10, input$volcano_export_width),
        height = ifelse(is.null(input$volcano_export_height), 7, input$volcano_export_height)
      )
      print(build_volcano_plot())
      dev.off()
    }
  )

  output$downloadUpGenes <- downloadHandler(
    filename = function() paste0(input$gse, "_upregulated_genes_", input$groupB, "_vs_", input$groupA, ".csv"),
    content = function(file) {
      write.csv(get_volcano_sig_table("up"), file, row.names = FALSE)
    }
  )

  output$downloadDownGenes <- downloadHandler(
    filename = function() paste0(input$gse, "_downregulated_genes_", input$groupB, "_vs_", input$groupA, ".csv"),
    content = function(file) {
      write.csv(get_volcano_sig_table("down"), file, row.names = FALSE)
    }
  )

  output$downloadVolcanoSigGenes <- downloadHandler(
    filename = function() paste0(input$gse, "_significant_genes_", input$groupB, "_vs_", input$groupA, ".csv"),
    content = function(file) {
      write.csv(get_volcano_sig_table("all"), file, row.names = FALSE)
    }
  )

  output$downloadHeatmapPNG <- downloadHandler(
    filename = function() paste0(input$gse, "_heatmap_", input$groupB, "_vs_", input$groupA, ".png"),
    content = function(file) {
      png(
        file,
        width = ifelse(is.null(input$heatmap_export_width), 10, input$heatmap_export_width),
        height = ifelse(is.null(input$heatmap_export_height), 9, input$heatmap_export_height),
        units = "in",
        res = as.numeric(ifelse(is.null(input$heatmap_export_dpi), 600, input$heatmap_export_dpi))
      )
      draw_heatmap(silent = FALSE)
      dev.off()
    }
  )

  output$downloadHeatmapPDF <- downloadHandler(
    filename = function() paste0(input$gse, "_heatmap_", input$groupB, "_vs_", input$groupA, ".pdf"),
    content = function(file) {
      pdf(
        file,
        width = ifelse(is.null(input$heatmap_export_width), 10, input$heatmap_export_width),
        height = ifelse(is.null(input$heatmap_export_height), 9, input$heatmap_export_height)
      )
      draw_heatmap(silent = FALSE)
      dev.off()
    }
  )

  output$downloadHeatmapSVG <- downloadHandler(
    filename = function() paste0(input$gse, "_heatmap_", input$groupB, "_vs_", input$groupA, ".svg"),
    content = function(file) {
      svglite::svglite(
        file,
        width = ifelse(is.null(input$heatmap_export_width), 10, input$heatmap_export_width),
        height = ifelse(is.null(input$heatmap_export_height), 9, input$heatmap_export_height)
      )
      draw_heatmap(silent = FALSE)
      dev.off()
    }
  )

  output$downloadHeatmapGenes <- downloadHandler(
    filename = function() paste0(input$gse, "_heatmap_genes_", input$groupB, "_vs_", input$groupA, ".csv"),
    content = function(file) {
      hm <- get_heatmap_object()
      if (is.null(hm)) {
        write.csv(data.frame(message = "No heatmap genes available."), file, row.names = FALSE)
      } else {
        write.csv(data.frame(Symbol = hm$genes), file, row.names = FALSE)
      }
    }
  )
  output$downloadPCAPNG <- downloadHandler(
    filename = function() paste0(input$gse, "_PCA_", input$groupB, "_vs_", input$groupA, ".png"),
    content = function(file) {
      ggsave(
        file,
        plot = build_pca_plot(),
        width = ifelse(is.null(input$pca_export_width), 8, input$pca_export_width),
        height = ifelse(is.null(input$pca_export_height), 6, input$pca_export_height),
        dpi = 600,
        device = "png"
      )
    }
  )

  output$downloadPCAPDF <- downloadHandler(
    filename = function() paste0(input$gse, "_PCA_", input$groupB, "_vs_", input$groupA, ".pdf"),
    content = function(file) {
      ggsave(
        file,
        plot = build_pca_plot(),
        width = ifelse(is.null(input$pca_export_width), 8, input$pca_export_width),
        height = ifelse(is.null(input$pca_export_height), 6, input$pca_export_height),
        device = "pdf"
      )
    }
  )

  output$downloadPCASVG <- downloadHandler(
    filename = function() paste0(input$gse, "_PCA_", input$groupB, "_vs_", input$groupA, ".svg"),
    content = function(file) {
      svglite::svglite(
        file,
        width = ifelse(is.null(input$pca_export_width), 8, input$pca_export_width),
        height = ifelse(is.null(input$pca_export_height), 6, input$pca_export_height)
      )
      print(build_pca_plot())
      dev.off()
    }
  )

}

shinyApp(ui, server)
