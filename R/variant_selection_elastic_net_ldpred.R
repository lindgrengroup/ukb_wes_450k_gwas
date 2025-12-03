#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Swiss-Army-Knife R Script for Gene-Centric Analysis
# Modes: 
#   1. Elastic Net (Variable Selection)
#   2. LDpred2 (Polygenic Scores / Posterior Effects)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
	  library(optparse)
	    library(data.table)
	    library(bigsnpr)
	      library(bigstatsr)
})

# --- 0. ARGUMENT PARSING ---

option_list <- list(
		      make_option(c("-g", "--gene"), type="character", default=NULL, help="Target Gene Symbol (e.g. HNF1A)", metavar="character"),
		        make_option(c("-r", "--radius"), type="integer", default=2000000, help="Radius in bp (default 2MB)", metavar="number"),
		        make_option(c("-m", "--mode"), type="character", default="elastic_net", help="Analysis mode: 'elastic_net' or 'ldpred2'", metavar="mode"),
			  make_option(c("--test_n"), type="integer", default=NULL, help="Number of samples for test run (subset)", metavar="number"),
			  make_option(c("--seed"), type="integer", default=NULL, help="Random seed for reproducibility", metavar="number"),
			    make_option(c("--maf"), type="double", default=0.001, help="MAF threshold", metavar="number"),
			    make_option(c("--hwe"), type="double", default=1e-10, help="HWE threshold", metavar="number"),
			      make_option(c("--info"), type="double", default=0.8, help="Info score threshold", metavar="number"),
			      make_option(c("--pheno_file"), type="character", help="Path to phenotype file", metavar="path"),
			        make_option(c("--pheno_col"), type="character", help="Phenotype column name", metavar="string"),
			        make_option(c("--gff_path"), type="character", default="/mnt/project/resources/ensembl/Homo_sapiens.GRCh37.87.gff3.gz", help="Path to GFF3 file"),
				  make_option(c("--dnax_dir"), type="character", help="DNAnexus directory containing genotype files"),
				  make_option(c("--out_dir"), type="character", default=".", help="Local output directory (default: current dir)")
				  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene) || is.null(opt$pheno_file) || is.null(opt$dnax_dir)) {
	  print_help(opt_parser)
  stop("Missing required arguments: --gene, --pheno_file, or --dnax_dir", call.=FALSE)
}

# Set Seed if provided
if (!is.null(opt$seed)) {
	  cat(paste0("[INFO] Setting random seed to ", opt$seed, "\n"))
  set.seed(opt$seed)
}

# Create output directory (if not current dir)
if (opt$out_dir != ".") {
	    dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# --- 1. HELPER FUNCTIONS ---
# ------------------------------------------------------------------------------

get_gene_coords_r <- function(gene_symbol, gff_path) {
	  cat(paste0("[INFO] Searching for ", gene_symbol, " in GFF...\n"))
  cmd <- paste0("gunzip -c ", shQuote(gff_path), " | grep -v '^#' | awk '$3 == \"gene\" {print $1, $3, $4, $5, $9}' OFS='\t'")
    
    tryCatch({
	        gene_data <- fread(cmd = cmd, header = FALSE, sep = "\t", 
				                          col.names = c("chr", "feature", "start", "end", "attributes"),
							                         colClasses = c("character", "character", "integer", "integer", "character"))
		  }, error = function(e) stop(paste("Error reading GFF:", e$message)))
    
    regex_pattern <- paste0("Name=", gene_symbol, "(;|$)")
      target_row <- gene_data[grepl(regex_pattern, attributes), ]
      
      if (nrow(target_row) == 0) stop(paste("Gene", gene_symbol, "not found."))
        
        result <- target_row[1, ]
        return(list(chr = result$chr, start = result$start, end = result$end))
}

construct_paths <- function(gene_symbol, radius, maf, hwe, info, gff_path, dnax_dir) {
	  gene_info <- get_gene_coords_r(gene_symbol, gff_path)
  
  maf_label  <- paste0("maf", as.character(maf))
    hwe_label  <- paste0("hwe", as.character(hwe))
    info_label <- paste0("info", as.character(info))
      
      output_prefix <- sprintf("imputed_v3_qced_snps_%s_%s_%s_%s_radius%dbp_chr%s",
			                                  maf_label, hwe_label, info_label, gene_symbol, radius, gene_info$chr)
      
      return(list(
		      pgen_prefix_path = file.path(dnax_dir, output_prefix),
		          bed_prefix_path  = paste0(output_prefix, "_plink1"), # Store locally
		          output_prefix    = output_prefix,
			      chrom            = gene_info$chr
			    ))
}

load_genotype_data <- function(paths, plink2_exe) {
	  pgen_prefix <- paths$pgen_prefix_path
  bed_prefix  <- paths$bed_prefix_path
    
    if (!file.exists(paste0(pgen_prefix, ".pgen"))) stop(paste("PGEN not found:", pgen_prefix))
    
    if (!file.exists(paste0(bed_prefix, ".bed"))) {
	        cat("[INFO] Converting PGEN to BED...\n")
        cmd <- sprintf("%s --pfile %s --make-bed --out %s --allow-extra-chr --memory 14000", 
		                          shQuote(plink2_exe), shQuote(pgen_prefix), shQuote(bed_prefix))
	    system(cmd, ignore.stdout = FALSE)
	  }
      
      rds_path <- paste0(bed_prefix, ".rds")
      if (!file.exists(rds_path)) {
	          cat("[INFO] Reading BED into bigSNP...\n")
          snp_readBed(paste0(bed_prefix, ".bed"), backingfile = sub("\\.rds$", "", rds_path))
	    }
        
        # Attach and Fix Code256 for imputation safety
        obj <- snp_attach(rds_path)
        obj$genotypes$code256 <- bigsnpr::CODE_012
	  
	  # Simple Mean Imputation
	  cat("[INFO] Imputing missing values (mean2)...\n")
	  obj$genotypes <- snp_fastImputeSimple(obj$genotypes, method = "mean2", ncores = nb_cores())
	    
	    return(obj)
}

filter_invariant_variants <- function(G, ind_row) {
	  cat(paste("[INFO] Filtering invariant variants in subset N=", length(ind_row), "...\n"))
  scaler <- big_scale(center = TRUE, scale = TRUE)(G, ind.row = ind_row)
    keep_indices <- which(scaler$scale > 0 & !is.na(scaler$scale))
    cat(paste("       Removed", ncol(G) - length(keep_indices), "invariant variants.\n"))
      return(keep_indices)
}

# ------------------------------------------------------------------------------
# --- 2. CORE ALGORITHM FUNCTIONS ---
# ------------------------------------------------------------------------------

run_elastic_net_workflow <- function(G, map, y, covars, ind_row, ind_col, alpha=0.5) {
	  cat(paste0("[", Sys.time(), "] Starting Elastic Net (Alpha=", alpha, ")\n"))
  
  # Adjust phenotype
  y_adj <- if (!is.null(covars)) residuals(lm(y ~ covars)) else y
    
    # Run
    mod <- big_spLinReg(X = G, y.train = y_adj, ind.train = ind_row, ind.col = ind_col,
			                      alphas = c(alpha), K = 10, warn = FALSE, ncores = nb_cores())
    
    # Extract
    best_mod <- summary(mod, best.only = TRUE)
      coeffs <- best_mod$beta[[1]]
      
      # Map back indices
      effective_ind_col <- attr(mod, "ind.col")
        selected_rel_idx <- which(coeffs != 0)
        
        if (length(selected_rel_idx) == 0) return(NULL)
	  
	  selected_global_idx <- effective_ind_col[selected_rel_idx]
	  selected_betas <- coeffs[selected_rel_idx]
	    
	    var_ids <- if("rsid" %in% names(map)) map$rsid[selected_global_idx] else map$marker.ID[selected_global_idx]
	    
	    return(data.frame(
			          varid = var_ids,
				      chr   = map$chromosome[selected_global_idx],
				      pos   = map$physical.pos[selected_global_idx],
				          beta  = selected_betas
				        ))
}

run_ldpred2_full_workflow <- function(G, map, y, covars, ind_row, ind_col) {
	  cat(paste0("[", Sys.time(), "] Starting LDpred2 Workflow...\n"))
  
  # A. Marginal GWAS
  cat("       Step 1: Marginal GWAS...\n")
    gwas <- big_univLinReg(X = G, y.train = y, covar.train = covars, 
			                            ind.train = ind_row, ind.col = ind_col, ncores = nb_cores())
    
    # Calculate n_eff
    n_eff <- if(length(unique(y)) == 2) {
	        t <- table(y); 4 / (1/min(t) + 1/max(t))
      } else { length(ind_row) }
      
      # Format Summary Stats
      map_sub <- map[ind_col, ]
      df_beta <- data.frame(beta = gwas$estim, beta_se = gwas$std.err, n_eff = n_eff)
        
        # Filter Bad SE
        valid_stats <- which(df_beta$beta_se > 0 & !is.na(df_beta$beta_se))
        df_beta <- df_beta[valid_stats, ]
	  ind_col_clean <- ind_col[valid_stats]
	  map_sub <- map_sub[valid_stats, ]
	    
	    # B. In-Sample LD
	    cat("       Step 2: Calculating In-Sample LD...\n")
	    tmp_file <- tempfile(pattern = "corr")
	      corr <- snp_cor(G, ind.row = ind_row, ind.col = ind_col_clean, 
			                        size = 2000, infos.pos = map_sub$physical.pos, ncores = nb_cores())
	      corr_sfbm <- as_SFBM(corr, tmp_file, compact = TRUE)
	        
	        # C. LDpred2 Auto
	        cat("       Step 3: Running LDpred2-Auto...\n")
	        h2_init <- 0.1 # Scalar init
		  vec_p_init <- seq_log(1e-4, 0.9, length.out = 10)
		  
		  multi_auto <- snp_ldpred2_auto(corr = corr_sfbm, df_beta = df_beta, h2_init = h2_init,
						                                  vec_p_init = vec_p_init, ncores = nb_cores(), 
										                                   allow_jump_sign = FALSE, shrink_corr = 0.95)
		    
		    # D. Post-process
		    range_p <- sapply(multi_auto, function(auto) diff(range(auto$path_p_est)))
		    keep <- (range_p < 0.2)
		      if (sum(keep) == 0) keep <- rep(TRUE, length(multi_auto))
		      
		      final_betas <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
		        
		        # Cleanup
		        file.remove(paste0(tmp_file, ".sbk"))
		        
		        # Return Result
		        var_ids <- if("rsid" %in% names(map_sub)) map_sub$rsid else map_sub$marker.ID
			  
			  return(data.frame(
					        varid = var_ids,
						    chr = map_sub$chromosome,
						    pos = map_sub$physical.pos,
						        marginal_beta = df_beta$beta,
						        marginal_p = predict(gwas, log10=FALSE)[valid_stats],
							    ldpred2_beta = final_betas
							  ))
}

# ------------------------------------------------------------------------------
# --- 3. EXECUTION ---
# ------------------------------------------------------------------------------

# Setup
# plink2 <- ensure_plink2() # Removed: Assumes available in path
plink2 <- "plink2"
paths  <- construct_paths(opt$gene, opt$radius, opt$maf, opt$hwe, opt$info, opt$gff_path, opt$dnax_dir)

# Construct output filename suffix
out_suffix <- if(!is.null(opt$seed)) paste0("_seed", opt$seed) else ""

# Load Data
cat("[INFO] Loading Genotypes...\n")
obj <- load_genotype_data(paths, plink2)
G <- obj$genotypes
map <- obj$map
fam <- obj$fam

# Load Phenotype
cat("[INFO] Loading Phenotypes...\n")
pheno_dt <- fread(opt$pheno_file)
merged <- merge(fam, pheno_dt, by.x = "sample.ID", by.y = "IID", sort=FALSE)
y <- merged[[opt$pheno_col]]

# Handle Covariates (Tranche)
covars <- if ("sequencing_tranche" %in% names(merged)) {
	  model.matrix(~ as.factor(merged$sequencing_tranche) - 1)
} else { NULL }

# Filter Complete Cases
complete_idx <- which(!is.na(y) & (is.null(covars) | complete.cases(covars)))
ind_row <- match(merged$sample.ID[complete_idx], fam$sample.ID)
y <- y[complete_idx]
covars <- if(!is.null(covars)) covars[complete_idx, , drop=FALSE] else NULL

# --- OPTIONAL: SUBSET FOR TEST RUN ---
if (!is.null(opt$test_n)) {
	  n_samples <- min(opt$test_n, length(ind_row))
  cat(paste0("[TEST MODE] Subsetting to first ", n_samples, " samples...\n"))
    
    # Subset Indices
    ind_row <- ind_row[1:n_samples]
    y       <- y[1:n_samples]
      if(!is.null(covars)) covars <- covars[1:n_samples, , drop=FALSE]
}

# Filter Variants (Must happen AFTER subsetting ind_row)
keep_cols <- filter_invariant_variants(G, ind_row)

# Dispatch Mode
if (opt$mode == "elastic_net") {
	  
	  results <- run_elastic_net_workflow(G, map, y, covars, ind_row, keep_cols)
  
  if (!is.null(results)) {
	      fname <- file.path(opt$out_dir, paste0(paths$output_prefix, "_elastic_net", out_suffix, ".csv"))
      fwrite(results, fname)
          
          # Write just the list of IDs for easy extraction
          list_name <- file.path(opt$out_dir, paste0(paths$output_prefix, "_selected_variants", out_suffix, ".txt"))
          fwrite(list(results$varid), list_name, col.names=FALSE)
	      
	      cat(paste("[SUCCESS] Elastic Net results written to", fname, "\n"))
	    } else {
		        cat("[WARN] No variants selected.\n")
	      }
    
} else if (opt$mode == "ldpred2") {
	  
	  results <- run_ldpred2_full_workflow(G, map, y, covars, ind_row, keep_cols)
  
  fname <- file.path(opt$out_dir, paste0(paths$output_prefix, "_ldpred2_results", out_suffix, ".csv"))
    fwrite(results, fname)
    cat(paste("[SUCCESS] LDpred2 results written to", fname, "\n"))
      
} else {
	  stop("Unknown mode. Use 'elastic_net' or 'ldpred2'.")
}

# ------------------------------------------------------------------------------
# --- 4. CLEANUP ---
# ------------------------------------------------------------------------------

cat("[INFO] Cleaning up intermediate files...\n")
extensions <- c(".bed", ".bim", ".fam", ".log", ".bk", ".rds")
files_to_remove <- paste0(paths$bed_prefix_path, extensions)
removed_count <- 0

for (f in files_to_remove) {
	  if (file.exists(f)) {
		      file.remove(f)
    removed_count <- removed_count + 1
      }
}
cat(paste("[SUCCESS] Cleaned up", removed_count, "intermediate files.\n"))
