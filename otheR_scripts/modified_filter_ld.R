ld_filter_attempt <- function (data, interactive.filter = TRUE, filter.short.ld = "mac", 
          filter.long.ld = NULL, parallel.core = parallel::detectCores() - 
            1, filename = NULL, verbose = TRUE, ...) 
{
  if (!interactive.filter && is.null(filter.short.ld) && is.null(filter.long.ld)) {
    return(data)
  }
  if (interactive.filter) 
    verbose <- TRUE
  if (verbose) {
    cat("################################################################################\n")
    cat("############################## radiator::filter_ld #############################\n")
    cat("################################################################################\n")
  }
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) 
    message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", 
                               round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################# completed filter_ld ##############################\n"), 
          add = TRUE)
  rad.dots <- radiator_dots(func.name = as.list(sys.call())[[1]], 
                            fd = rlang::fn_fmls_names(), args.list = as.list(environment()), 
                            dotslist = rlang::dots_list(..., .homonyms = "error", 
                                                        .check_assign = TRUE), keepers = c("long.ld.missing", 
                                                                                           "ld.method", "ld.figures", "subsample.markers.stats", 
                                                                                           "path.folder", "parameters", "internal"), verbose = FALSE)
  if (missing(data)) 
    rlang::abort("data is missing")
  path.folder <- generate_folder(f = path.folder, rad.folder = "filter_ld", 
                                 internal = internal, file.date = file.date, verbose = verbose)
  write_rad(data = rad.dots, path = path.folder, filename = stringi::stri_join("radiator_filter_ld_args_", 
                                                                               file.date, ".tsv"), tsv = TRUE, write.message = "Function call and arguments stored in: ", 
            internal = internal, verbose = verbose)
  if (interactive.filter) {
    message("\nInteractive mode: on\n")
    message("Step 1. Short distance LD threshold selection")
    message("Step 2. Filtering markers based on short distance LD")
    message("Step 3. Long distance LD pruning selection")
    message("Step 4. Threshold selection")
    message("Step 5. Filtering markers based on long distance LD\n\n")
  }
  if (!is.null(filter.short.ld)) {
    if (filter.short.ld == "maf") {
      message("\n\nPlease update your code to use filter.short.ld = 'mac'\n\n")
      filter.short.ld <- "mac"
    }
    filter.short.ld <- match.arg(filter.short.ld, c("first", 
                                                    "random", "last", "middle", "mac"))
  }
  if (!is.null(ld.method)) {
    ld.method <- match.arg(ld.method, c("composite", "r", 
                                        "r2", "dprime", "corr"))
  }
  if (is.null(filename)) {
    write.ld <- FALSE
    filename <- stringi::stri_join("radiator_", file.date, 
                                   ".ld")
  }
  else {
    write.ld <- TRUE
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_", file.date, 
                                     ".ld")
    }
    else {
      filename <- stringi::stri_join(filename, ".ld")
    }
  }
  filename.gds <- stringi::stri_join(filename, ".gds")
  filename.gds <- file.path(path.folder, filename.gds)
  filename <- file.path(path.folder, filename)
  filename.gds.rad <- stringi::stri_join(filename.gds, ".rad")
  data.type <- radiator::detect_genomic_format(data)
  if (!data.type %in% c("tbl_df", "fst.file", "SeqVarGDSClass", 
                        "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  if (data.type %in% c("tbl_df", "fst.file")) {
    if (is.vector(data)) 
      data <- radiator::tidy_wide(data = data, import.metadata = FALSE)
    data.type <- "tbl_df"
    wl <- bl <- dplyr::select(data, MARKERS, CHROM, LOCUS, 
                              POS) %>% dplyr::distinct(MARKERS, .keep_all = TRUE) %>% 
      dplyr::arrange(LOCUS, MARKERS)
  }
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (!"SeqVarTools" %in% utils::installed.packages()[, 
                                                        "Package"]) {
      rlang::abort("Please install SeqVarTools for this option:\n\n           install.packages(\"BiocManager\")\n           BiocManager::install(\"SeqVarTools\")")
    }
    if (data.type == "gds.file") {
      data <- radiator::read_rad(data, verbose = verbose)
      data.type <- "SeqVarGDSClass"
    }
    wl <- bl <- extract_markers_metadata(data, whitelist = TRUE) %>% 
      dplyr::arrange(LOCUS, MARKERS)
  }
  filters.parameters <- radiator_parameters(generate = TRUE, 
                                            initiate = TRUE, update = FALSE, parameter.obj = parameters, 
                                            data = data, path.folder = path.folder, file.date = file.date, 
                                            internal = internal, verbose = verbose)
  if (interactive.filter || !is.null(filter.short.ld)) {
    if (verbose) 
      message("Minimizing short distance LD...")
    locus.stats <- dplyr::group_by(.data = wl, LOCUS) %>% 
      dplyr::summarise(SNP_N = n()) %>% dplyr::count(SNP_N) %>% 
      readr::write_tsv(x = ., path = file.path(path.folder, 
                                               "short.ld.locus.stats.tsv"), append = FALSE, 
                       col_names = TRUE)
    if (nrow(locus.stats) > 1) {
      range.number.snp.locus <- range(locus.stats$SNP_N, 
                                      na.rm = TRUE)
      if (verbose) 
        message("    The range in the number of SNP/locus is: ", 
                stringi::stri_join(range.number.snp.locus, 
                                   collapse = "-"))
      if (interactive.filter) {
        if (verbose) 
          message("\nStep 1. Short distance LD threshold selection")
        if (verbose) 
          message("the goal is to keep only 1 SNP per read/locus")
      }
      if (interactive.filter) {
        filter.short.ld <- radiator_question(x = "Choose the filter.short.ld threshold\nOptions include:\n1: mac (Not sure ? use mac...)\n2: random\n3: first\n4: middle\n5: last", 
                                             answer.opt = c("1", "2", "3", "4", "5"))
        filter.short.ld <- stringi::stri_replace_all_fixed(str = filter.short.ld, 
                                                           pattern = c("1", "2", "3", "4", "5"), replacement = c("mac", 
                                                                                                                 "random", "first", "middle", "last"), vectorize_all = FALSE)
        filter.short.ld <- match.arg(filter.short.ld, 
                                     c("mac", "random", "first", "middle", "last"))
      }
      if (verbose) 
        message("\nStep 2. Filtering markers based on short distance LD")
      if (verbose) 
        message("filter.short.ld = ", filter.short.ld)
      if (filter.short.ld == "random") {
        wl %<>% dplyr::group_by(LOCUS) %>% dplyr::sample_n(tbl = ., 
                                                           size = 1, replace = FALSE) %>% dplyr::ungroup(.)
      }
      if (filter.short.ld == "first") {
        wl %<>% dplyr::group_by(LOCUS) %>% dplyr::summarise(POS = min(POS)) %>% 
          dplyr::ungroup(.)
      }
      if (filter.short.ld == "last") {
        wl %<>% dplyr::group_by(LOCUS) %>% dplyr::summarise(POS = max(POS)) %>% 
          dplyr::ungroup(.)
      }
      if (filter.short.ld == "middle") {
        snp.locus.prep <- dplyr::group_by(.data = wl, 
                                          LOCUS) %>% dplyr::tally(.) %>% dplyr::ungroup(.)
        pick.middle <- snp.locus.prep %>% dplyr::filter(n > 
                                                          2) %>% dplyr::select(LOCUS)
        if (nrow(pick.middle) == 0) {
          if (verbose) 
            message("IMPORTANT: the data doesn't have more than 3 SNPs per locus")
          if (verbose) 
            message("    First SNP will be selected instead...")
          wl %<>% dplyr::group_by(LOCUS) %>% dplyr::summarise(POS = min(POS)) %>% 
            dplyr::ungroup(.)
        }
        else {
          keep.first <- snp.locus.prep %>% dplyr::filter(n <= 
                                                           2) %>% dplyr::select(LOCUS)
          if (verbose) 
            message("    Number of locus with first SNP selected: ", 
                    nrow(keep.first))
          keep.first.select <- wl %>% dplyr::filter(LOCUS %in% 
                                                      keep.first$LOCUS) %>% dplyr::group_by(LOCUS) %>% 
            dplyr::summarise(POS = min(POS)) %>% dplyr::ungroup(.)
          pick.middle.select <- wl %>% dplyr::filter(LOCUS %in% 
                                                       pick.middle$LOCUS) %>% dplyr::group_by(LOCUS) %>% 
            dplyr::filter(POS != min(POS)) %>% dplyr::filter(POS != 
                                                               max(POS)) %>% dplyr::sample_n(tbl = ., size = 1, 
                                                                                             replace = FALSE)
          if (verbose) 
            message("    Number of locus with random middle SNP selected: ", 
                    nrow(pick.middle))
          wl <- dplyr::bind_rows(keep.first.select, pick.middle.select) %>% 
            dplyr::arrange(LOCUS, POS)
        }
        pick.middle <- snp.locus.prep <- keep.first.select <- pick.middle.select <- keep.first <- NULL
      }
      if (filter.short.ld == "mac") {
        if (data.type == "tbl_df") {
          one.snp <- dplyr::group_by(.data = wl, LOCUS) %>% 
            dplyr::tally(.) %>% dplyr::filter(n == 1) %>% 
            dplyr::left_join(wl, by = "LOCUS") %>% dplyr::select(-n) %>% 
            dplyr::distinct(MARKERS, .keep_all = TRUE)
          if (tibble::has_name(data, "GT_BIN")) {
            more.snp <- dplyr::distinct(wl, MARKERS) %>% 
              dplyr::filter(!MARKERS %in% one.snp$MARKERS)
            n.markers <- nrow(more.snp)
            global_mac <- function(x) {
              mac.data <- dplyr::group_by(x, MARKERS) %>% 
                dplyr::summarise(PP = as.numeric(2 * 
                                                   length(GT_BIN[GT_BIN == 0])), PQ = as.numeric(length(GT_BIN[GT_BIN == 
                                                                                                                 1])), QQ = as.numeric(2 * length(GT_BIN[GT_BIN == 
                                                                                                                                                           2]))) %>% dplyr::mutate(PP = PP + PQ, 
                                                                                                                                                                                   QQ = QQ + PQ, PQ = NULL, MAC_GLOBAL = dplyr::if_else(PP < 
                                                                                                                                                                                                                                          QQ, PP, QQ), PP = NULL, QQ = NULL) %>% 
                dplyr::ungroup(.)
              return(mac.data)
            }
            if (n.markers > 10000) {
              split.vec <- more.snp %>% dplyr::mutate(SPLIT_VEC = split_vec_row(more.snp, 
                                                                                cpu.rounds = ceiling(n.markers/10000), 
                                                                                parallel.core = parallel.core))
              mac.data <- data %>% dplyr::filter(!is.na(GT_BIN)) %>% 
                dplyr::filter(!MARKERS %in% one.snp$MARKERS) %>% 
                dplyr::left_join(split.vec, by = "MARKERS") %>% 
                split(x = ., f = .$SPLIT_VEC) %>% .radiator_parallel_mc(X = ., 
                                                                        FUN = global_mac, mc.cores = parallel.core) %>% 
                dplyr::bind_rows(.)
              more.snp <- split.vec <- NULL
            }
            else {
              mac.data <- global_mac(x = dplyr::filter(data, 
                                                       !is.na(GT_BIN)) %>% dplyr::filter(!MARKERS %in% 
                                                                                           one.snp$MARKERS))
            }
          }
          else {
            mac.data <- data %>% dplyr::filter(GT != 
                                                 "000000") %>% dplyr::filter(!MARKERS %in% 
                                                                               one.snp$MARKERS) %>% dplyr::select(MARKERS, 
                                                                                                                  INDIVIDUALS, GT) %>% dplyr::mutate(A1 = stringi::stri_sub(GT, 
                                                                                                                                                                            1, 3), A2 = stringi::stri_sub(GT, 4, 6)) %>% 
              dplyr::select(MARKERS, INDIVIDUALS, A1, 
                            A2) %>% tidyr::gather(data = ., key = ALLELES, 
                                                  value = GT, -c(MARKERS, INDIVIDUALS)) %>% 
              dplyr::group_by(MARKERS, GT) %>% dplyr::tally(.) %>% 
              dplyr::group_by(MARKERS) %>% dplyr::filter(n == 
                                                           min(n)) %>% dplyr::distinct(MARKERS, .keep_all = TRUE) %>% 
              dplyr::summarise(MAC_GLOBAL = n) %>% dplyr::ungroup(.) %>% 
              dplyr::select(MARKERS, MAC_GLOBAL)
          }
          wl %<>% dplyr::filter(!MARKERS %in% one.snp$MARKERS) %>% 
            dplyr::left_join(mac.data, by = "MARKERS") %>% 
            dplyr::group_by(LOCUS) %>% dplyr::filter(MAC_GLOBAL == 
                                                       max(MAC_GLOBAL)) %>% dplyr::ungroup(.) %>% 
            dplyr::select(-MAC_GLOBAL) %>% dplyr::distinct(LOCUS, 
                                                           .keep_all = TRUE) %>% dplyr::bind_rows(one.snp)
          one.snp <- mac.data <- NULL
        }
        else {
          n.markers <- nrow(wl)
          if (!tibble::has_name(wl, "MAC_GLOBAL")) {
            wl %<>% dplyr::bind_cols(SeqArray::seqAlleleCount(gdsfile = data, 
                                                              ref.allele = NULL, .progress = TRUE, parallel = parallel.core) %>% 
                                       unlist(.) %>% matrix(data = ., nrow = n.markers, 
                                                            ncol = 2, byrow = TRUE, dimnames = list(rownames = wl$MARKERS, 
                                                                                                    colnames = c("REF_COUNT", "ALT_COUNT"))) %>% 
                                       tibble::as_tibble(.)) %>% dplyr::mutate(MAC_GLOBAL = dplyr::if_else(ALT_COUNT < 
                                                                                                             REF_COUNT, ALT_COUNT, REF_COUNT), ALT_COUNT = NULL, 
                                                                               REF_COUNT = NULL)
          }
          wl %<>% dplyr::group_by(LOCUS) %>% dplyr::filter(MAC_GLOBAL == 
                                                             max(MAC_GLOBAL)) %>% dplyr::ungroup(.) %>% 
            dplyr::distinct(LOCUS, .keep_all = TRUE)
        }
      }
      readr::write_tsv(x = wl, path = file.path(path.folder, 
                                                "whitelist.short.ld.tsv"), append = FALSE, col_names = TRUE)
      bl %<>% dplyr::filter(!MARKERS %in% wl$MARKERS) %>% 
        dplyr::mutate(FILTER = "filter.short.ld")
      readr::write_tsv(x = bl, path = file.path(path.folder, 
                                                "blacklist.short.ld.tsv"), append = FALSE, col_names = TRUE)
      if (verbose) 
        message("File written: whitelist.short.ld.tsv")
      if (verbose) 
        message("File written: blacklist.short.ld.tsv")
      if (data.type == "tbl_df") {
        data <- dplyr::filter(data, MARKERS %in% wl$MARKERS)
      }
      else {
        markers.meta <- extract_markers_metadata(gds = data) %>% 
          dplyr::mutate(FILTERS = dplyr::if_else(MARKERS %in% 
                                                   bl$MARKERS, "filter.short.ld", FILTERS))
        update_radiator_gds(gds = data, node.name = "markers.meta", 
                            value = markers.meta, sync = TRUE)
      }
    }
    else {
      if (verbose) 
        message("\nThere is no variation in the number of SNP/locus across the data\n")
    }
    locus.stats <- mac.data <- bl <- wl <- short.ld.fig <- NULL
    filters.parameters <- radiator_parameters(generate = FALSE, 
                                              initiate = FALSE, update = TRUE, parameter.obj = filters.parameters, 
                                              data = data, filter.name = "Filter short ld", param.name = "filter.short.ld", 
                                              values = filter.short.ld, path.folder = path.folder, 
                                              file.date = file.date, internal = internal, verbose = verbose)
    radiator_results_message(rad.message = stringi::stri_join("\nFilter short ld threshold: ", 
                                                              filter.short.ld), filters.parameters, internal, verbose)
  }
  long.ld <- "y"
  if (interactive.filter) {
    long.ld <- radiator_question(x = "\nDo you want to continue filtering using long distance ld  ? (y/n):", 
                                 answer.opt = c("y", "n"))
    if (long.ld == "n") {
      return(data)
    }
  }
  if (interactive.filter || !is.null(filter.long.ld)) {
    ref.genome <- detect_ref_genome(data = data, verbose = FALSE)
    if (interactive.filter) {
      message("\nStep 3. Long distance LD pruning selection")
      if (!ref.genome) {
        message("Using missingness to select SNPs in LD is still under construction with de novo data")
        message("Basic pruning using SNPRelate is used")
        long.ld.missing <- FALSE
      }
      else {
        message("With a reference genome, pruning is done by chromosome/scaffolds")
        message("Pruning method can randomly choose to keep 1 SNP or")
        message("select the SNP based on missing data statistics")
        long.ld.missing <- radiator_question(x = "\nDo you want to use missing data statistics ? (y/n):", 
                                             answer.opt = c("y", "n"))
        if (long.ld.missing == "y") {
          long.ld.missing <- TRUE
        }
        else {
          long.ld.missing <- FALSE
        }
      }
    }
    if (data.type == "tbl_df") {
      biallelic <- radiator::detect_biallelic_markers(data = data)
      if (!biallelic) 
        rlang::abort("Long distance LD: biallelic genotypes required")
      if (verbose) 
        message("\nPreparing the data long LD filtering...")
      data.gds <- radiator::write_gds(data = data, filename = filename, 
                                      verbose = FALSE)
      markers$MARKERS <- markers$MARKERS
      data$MARKERS <- data$MARKERS
      if (verbose) 
        message("SNPRelate GDS file generated: ", filename.gds)
      if (verbose) 
        message("To close the connection use SNPRelate::snpgdsClose(filename)")
    }
    wl <- extract_markers_metadata(data, whitelist = TRUE)
    n.chrom <- length(unique(wl$CHROM))
    if (ref.genome) {
      denovo <- FALSE
      chrom.tick <- dplyr::distinct(wl, CHROM) %>% dplyr::mutate(CHROM_TICK = stringi::stri_join(seq(from = 1, 
                                                                                                     to = n(), by = 1), n(), sep = "/"))
      ld.sample <- dplyr::sample_frac(tbl = chrom.tick, 
                                      size = 0.2) %>% dplyr::select(CHROM) %>% purrr::flatten_chr(.)
      wl %<>% dplyr::left_join(chrom.tick, by = "CHROM") %>% 
        dplyr::mutate(LD_SUBSAMPLE = dplyr::if_else(CHROM %in% 
                                                      ld.sample, TRUE, FALSE))
      chrom.tick <- chrom.tick$CHROM
      n.markers <- nrow(wl)
      wl %<>% dplyr::arrange(CHROM, MARKERS) %>% dplyr::mutate(CHROM = factor(x = CHROM, 
                                                                                 levels = chrom.tick, ordered = TRUE)) %>% dplyr::arrange(CHROM)
    }
    bl <- wl
    ld.sample <- NULL
    if (long.ld.missing) {
      if (verbose) 
        message("\nLong distance LD pruning with missing data")
      wl.bl.ld <- ld_missing(wl = wl, data = data, ld.threshold = if (interactive.filter) {
        seq(0.1, 0.9, by = 0.1)
      }
      else {
        filter.long.ld
      }, denovo = denovo, ld.method = ld.method, ld.figures = ld.figures, 
      parallel.core = parallel.core, verbose = verbose, 
      path.folder = path.folder)
      if (interactive.filter) {
        if (verbose) 
          message("\nStep 4. Threshold selection")
        if (verbose) 
          message("Look at the boxplot, a threshold of 0.2 will blacklist more markers than a threshold of 0.8")
        filter.long.ld <- radiator_question(x = "\nEnter the long LD threshold (filter.long.ld threshold, double/proportion):", 
                                            minmax = c(0, 1))
      }
      if (interactive.filter) 
        message("\nStep 5. Filtering markers based on long distance LD")
      wl.bl.ld <- magrittr::extract2(wl.bl.ld, as.name(filter.long.ld))
      wl <- wl.bl.ld %$% wl
      bl <- wl.bl.ld %$% bl
      if (data.type == "tbl_df") {
        data <- dplyr::filter(data, MARKERS %in% wl$MARKERS)
      }
      else {
        markers.meta <- extract_markers_metadata(gds = data) %>% 
          dplyr::mutate(FILTERS = dplyr::if_else(MARKERS %in% 
                                                   bl, "filter.long.ld", FILTERS))
        update_radiator_gds(gds = data, node.name = "markers.meta", 
                            value = markers.meta, sync = TRUE)
      }
    }
    if (!long.ld.missing) {
      if (verbose) 
        message("\nLong distance LD pruning WITHOUT missing data stats")
      if (is.null(subsample.markers.stats)) {
        subsample.markers.stats <- 0.2
      }
      else {
        if (!is.double(subsample.markers.stats)) 
          subsample.markers.stats <- 0.2
      }
      bp <- ld_boxplot(gds = data, subsample.markers.stats = subsample.markers.stats, 
                       ld.method = ld.method, path.folder = path.folder, 
                       parallel.core = parallel.core, verbose = verbose)
      bp <- NULL
      if (interactive.filter) {
        if (verbose) 
          message("\nStep 4. Threshold selection")
        if (verbose) 
          message("Look at the boxplot, a threshold of 0.2 blacklist more markers than 0.8")
        filter.long.ld <- radiator_question(x = "\nEnter the long LD threshold (filter.long.ld threshold, double/proportion):", 
                                            minmax = c(0, 1))
      }
      if (ld.method == "r2") {
        ld.m <- "r"
      }
      else {
        ld.m <- ld.method
      }
      id <- extract_individuals_metadata(gds = data, ind.field.select = "INDIVIDUALS", 
                                         whitelist = TRUE) %$% INDIVIDUALS
      if (verbose) 
        message("Pruning with SNPRelate...")
      timing <- proc.time()
      wl.variant.id <- SNPRelate::snpgdsLDpruning(gdsobj = data, 
                                                  snp.id = wl$MARKERS, sample.id = id, autosome.only = FALSE, 
                                                  remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, 
                                                  method = ld.m, ld.threshold = filter.long.ld, 
                                                  num.thread = 1L, verbose = TRUE) %>% unlist(.)
      timing <- proc.time() - timing
      if (verbose) 
        message("LD pruning computation time: ", round(timing[[3]]), 
                " sec")
      wl.n <- length(wl.variant.id)
      message("Number of markers whitelised: ", wl.n)
      wl %<>% dplyr::filter(MARKERS %in% wl.variant.id)
      bl %<>% dplyr::setdiff(wl) %>% dplyr::mutate(FILTERS = "filter.long.ld")
      write_rad(data = wl, path = path.folder, filename = "whitelist.long.ld.tsv", 
                tsv = TRUE, internal = FALSE, verbose = verbose)
      write_rad(data = bl, path = path.folder, filename = "blacklist.long.ld.tsv", 
                tsv = TRUE, internal = FALSE, verbose = verbose)
      wl.variant.id <- NULL
      if (data.type == "tbl_df") {
        data <- dplyr::filter(data, MARKERS %in% wl$MARKERS)
      }
      else {
        markers.meta <- extract_markers_metadata(gds = data) %>% 
          dplyr::mutate(FILTERS = dplyr::if_else(MARKERS %in% 
                                                   bl$MARKERS, "filter.long.ld", FILTERS))
        update_radiator_gds(gds = data, node.name = "markers.meta", 
                            value = markers.meta, sync = TRUE)
      }
    }
    filters.parameters <- radiator_parameters(generate = FALSE, 
                                              initiate = FALSE, update = TRUE, parameter.obj = filters.parameters, 
                                              data = data, filter.name = "Filter long ld", param.name = paste0("filter.long.ld / long.ld.missing"), 
                                              values = stringi::stri_join(filter.long.ld, long.ld.missing, 
                                                                          collapse = " / ", ignore_null = FALSE), path.folder = path.folder, 
                                              file.date = file.date, internal = internal, verbose = verbose)
    radiator_results_message(rad.message = stringi::stri_join("\nFilter long ld threshold: ", 
                                                              filter.long.ld), filters.parameters, internal, verbose)
  }
  return(data)
}
