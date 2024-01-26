


#-------------------------------------------------------------------------------
#   simple_functions
#-------------------------------------------------------------------------------

# ensure_finished_jobs: Make sure all qsub jobs finish running before proceeding
# message2:             Just like message, but with the date and a gap pasted in
# message_header:       Nicely-formatted text to break up your log files into 
#                       readable chunks
# n_eff:                Get the effective sample size of a case-control GWAS
# p_to_z:               Convert a P value into a z-score

# ensure_finished_jobs: Make sure all qsub jobs finish running before proceeding
ensure_finished_jobs <- function(identifier){
  
  external.call <- paste0( "squeue | grep ", identifier, " | wc -l" )
  running.jobs  <- as.numeric( system( external.call, intern=TRUE ) ) 
  total_sleep   <- 0
  while( running.jobs > 0){
    message( "Waited for ", total_sleep, " minutes, there are still ",
             running.jobs, " jobs running")
    Sys.sleep(60)
    total_sleep <- total_sleep + 1
    running.jobs <- as.numeric( system( external.call, intern=TRUE ) )
  }
}

# message2: Just like message, but with the date and a gap pasted in
message2 <- function(...) message( date(), "     ", ...)

# message_header: 
# Nicely-formatted text to break up your log files into readable chunks
message_header <- function(...){
  message( "\n\n#--------------------------------------------------" )
  message( "#   ", ...)
  message( "#--------------------------------------------------" )
}

# n_eff: Get the effective sample size of a case-control GWAS
n_eff <- function( n.cases, n.controls ){
  n.eff <- 4 / ( (1/n.cases) + (1/n.controls) )
  round(n.eff)
}

# p_to_z: Convert a P value into a z-score
p_to_z <- function( p, direction=NULL, limit=.Machine$double.xmin, log.p=FALSE ){
  
  # Set P value lower limit to avoid Inf/-Inf
  if( !is.null( limit ) )  p[ which( p < limit ) ] <- limit
  
  # Get z
  if(log.p){
    z <- -qnorm( p - log(2), log.p=TRUE )
  }else{
    z <- -qnorm(p/2)
  }
  
  # Correct sign, return
  if ( !is.null( direction) )  z <-  z * sign(direction)
  z
}


#-------------------------------------------------------------------------------
#   check_arguments:         Check validity of inputs
#-------------------------------------------------------------------------------

check_arguments <- function( ld_panel         = NULL, 
                             population       = NULL,
                             gw_file          = NULL,
                             chr_bp_col       = NULL,
                             chr_col          = NULL,
                             bp_col           = NULL,
                             a1_col           = NULL,
                             a2_col           = NULL,
                             p_col            = NULL,
                             eaf_col          = NULL,
                             n1_col           = NULL,
                             n0_col           = NULL,
                             n_col            = NULL,
                             n                = NULL,
                             p_threshold      = NULL, 
                             inflation_factor = NULL ){
  
  # ld_panel must be either "hrc" or "g1000"
  if( ld_panel != "hrc" & ld_panel != "g1000" ) stop("ld_panel must be either 'hrc' or 'g1000'")
  
  # population must be one of: "eur", "eas", or "eur_eas"
  if( population != "eur" & population != "eas" & population != "eur_eas" ){
    stop("population must be one of: 'eur', 'eas', or 'eur_eas'")
  }
  
  # Does the GWAS file exist?
  if( !file.exists(gw_file) )  stop("GWAS file does not exist")
  
  # Either chr_bp_col or (chr_col + bp_col) must be specified
  if( is.null(chr_bp_col) ){
    if( is.null(chr_col) | is.null(bp_col) ){
      stop("chr_bp_col is not specified so both chr_col and bp_col must be specified")
    }
  }else{
    if( !is.null(chr_col) | !is.null(bp_col) ){
      stop("chr_bp_col is specified so chr_col and bp_col must not be specified")
    }
  }
  
  # One of the following must be specified:
  # n or n_col or (n1_col + n0_col)
  if( !is.null(n) ){
    if( !is.null(n_col) | !is.null(n1_col) | !is.null(n0_col) ){
      stop("n is specified so n_col, n1_col, and n0_col must not be specified")
    }
  }else if( !is.null(n_col) ){
    if( !is.null(n1_col) | !is.null(n0_col) ){
      stop("n_col is specified so n1_col and n0_col must not be specified")
    }
  }else{
    if( is.null(n1_col) | is.null(n0_col) ){
      stop("Neither n nor n_col are specified so both n1_col and n0_col must be specified")
    }
  }
  
  # Check that GWAS file column names exist
  col.names <- c( chr_bp_col, chr_col, bp_col, a1_col, a2_col, 
                  p_col, eaf_col, n1_col, n0_col, n_col )
  library(data.table)
  gw <- fread( file=gw.file, nrows=100 )
  bad_col_names <- setdiff( col.names, names(gw) )
  if( length(bad_col_names) > 0 ){
    stop("The following specified column names do not exist: ", 
         paste( bad_col_names, collapse=", " ) )
  }
  
  # Position must be a positive integer
  if( !is.null(bp_col) ){
    if( !is.integer( gw[[bp_col]] ) ) stop("Positions must be integers")
    if( any( gw[[bp_col]] ) < 1  )    stop("Positions must be > 0")
  }
  
  # Alleles must be characters
  if( !is.character( gw[[a1_col]] ) ) stop("Effect alleles must be characters")
  if( !is.character( gw[[a2_col]] ) ) stop("Non-effect alleles must be characters")
  
  # P value must be >= 0 and <= 1
  if( !is.numeric( gw[[p_col]] ) ) stop("P values must be numeric")
  if( any(   gw[[p_col]] ) < 0 )   stop("P values must be >= 0")
  if( any(   gw[[p_col]] ) > 1 )   stop("P values must be <= 1")
  
  # Effect allele frequency must be >= 0 and <= 1
  if( eaf_col %in% names(gw) ){
    if( !is.numeric( gw[[eaf_col]] ) ) stop("Effect allele frequencies must be numeric")
    if( any(   gw[[eaf_col]] ) < 0 )   stop("Effect allele frequencies must be >= 0")
    if( any(   gw[[eaf_col]] ) > 1 )   stop("Effect allele frequencies must be <= 1")
  }
  
  # Case counts must be positive integers
  if( !is.null(n1_col) ){
    if( !is.integer( gw[[n1_col]] ) ) stop("Case counts must be integers")
    if( any( gw[[n1_col]] < 1 ) )     stop("Case counts must be > 0")
  }
  
  # Control counts must be positive integers
  if( !is.null(n0_col) ){
    if( !is.integer( gw[[n0_col]] ) ) stop("Control counts must be integers")
    if( any( gw[[n0_col]] < 1 ) )     stop("Control counts must be > 0")
  }
  
  # Effective sample sizes must be positive numbers
  if( !is.null(n_col) ){
    if( !is.numeric( gw[[n_col]] ) ) stop("Effective sample sizes must be numbers")
    if( any( gw[[n_col]] <= 0 ) )    stop("Effective sample sizes must be > 0")
  }
  
  # Effective sample size must be a positive number
  if( !is.null(n) ){
    if( !is.numeric(n) ) stop("Effective sample size must be a number")
    if( n <= 0 )         stop("Effective sample size must be > 0")
  }
  
  # P value threshold must be >= 0 and <= 1
  if( !is.null(p_threshold) ){
    if( !is.numeric(p_threshold) ) stop("P value threshold must be numeric")
    if( p_threshold < 0 )          stop("P value threshold must be >= 0")
    if( p_threshold > 1 )          stop("P value threshold must be <= 1")
  }
  
  # Inflation factor must be a positive number
  if( !is.null(inflation_factor) ){
    if( !is.numeric(inflation_factor) ) stop("Inflation factor must be a number")
    if( inflation_factor <= 0 )         stop("Inflation factor must be > 0")
  }
}


#-------------------------------------------------------------------------------
#   format_gwas:             Write GWAS sumstats in COJO format
#-------------------------------------------------------------------------------

format_gwas <- function( maindir          = "~/cojo",
                         ld_panel         = "hrc",
                         population       = "eur_eas",
                         gw_file          = "/home/heilbron/projects/pops/analyses/pd/meta5_raw.tab.gz",
                         chr_bp_col       = "SNP",
                         chr_col          = NULL,
                         bp_col           = NULL,
                         a1_col           = "A1",
                         a2_col           = "A2",
                         b_col            = NULL,
                         or_col           = NULL,
                         se_col           = NULL,
                         p_col            = "p",
                         eaf_col          = "freq",
                         n1_col           = "N_cases",
                         n0_col           = "N_controls",
                         n_col            = NULL,
                         n                = NULL,
                         p_threshold      = 1, 
                         inflation_factor = 1 ){
  
  # If output files already exist, skip
  # If specified, correct for inflation
  # Write the list of SNP names to file
  
  
  #-----------------------------------------------------------------------------
  #   Input descriptions
  #-----------------------------------------------------------------------------
  
  #   maindir:    Main directory in which to store results
  #   ld_panel:   Which LD reference panel should be used? Options are either: 
  #               "hrc" or "g1000".
  #   gw_file:    GWAS file name
  #   chr_bp_col: Optional. The name of a GWAS column containing chromosome and bp 
  #               information separated by a punctuation character. Must be
  #               specified if chr_col and bp_col are not specified.
  #   chr_col:    Optional. The name of a GWAS column containing chromosome
  #               information. Must be specified if chr_bp_col is not specified.
  #   bp_col:     Optional. The name of a GWAS column containing position (bp)
  #               information. Must be specified if chr_bp_col is not specified.
  #   a1_col:     The name of a GWAS column containing the effect (alt) allele
  #   a2_col:     The name of a GWAS column containing the non-effect (ref) allele
  #   p_col:      The name of a GWAS column containing the P value
  #   eaf_col     Optional. The name of a GWAS column containing the frequency of
  #               the effect (alt) allele. If not specified, all palindromic SNPs
  #               will be removed. Otherwise palindromic SNPs with similar 
  #               frequencies in GWAS and HRC will be preserved.
  #   n1_col:     The name of a GWAS column containing the per-SNP number of 
  #               cases. Must be specified with n0_col, or must specify n_col 
  #               or n.
  #   n0_col:     The name of a GWAS column containing the per-SNP number of 
  #               controls. Must be specified with n1_col, or must specify n_col 
  #               or n.
  #   n_col:      The name of a GWAS column containing the per-SNP effective
  #               sample size. If not specified, must specify both n1_col and 
  #               n0_col, or n.
  #   n:          If per-SNP sample size information is not available, this
  #               specifies the study-wide effective sample size. If not 
  #               specified, must specify both n1_col and n0_col, or n_col.
  #   check_args: Logical. Check whether arguments are valid?
  
  
  #-----------------------------------------------------------------------------
  #   Read in GWAS and HRC, format columns
  #-----------------------------------------------------------------------------
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  library(data.table)
  
  # Create main directory
  message2("Create main directory")
  dir.create( path=maindir, showWarnings=FALSE, recursive=TRUE )
  
  # Read in reference panel SNPs
  rare.or.common.snps <- "rare"
  if( ld_panel == "hrc" ){
    
    # Population: Merged EUR and EAS
    if( population == "eur_eas" ){
      
      # Common SNPs (MAF >= 1%)
      if( rare.or.common.snps == "common" ){
        message2("Read in HRC SNPs with EUR+EAS MAF >= 1%")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eur_eas_snps_maf_ge_0.01.tsv")
        
      # Rare SNPs (MAC >= 10)
      }else if( rare.or.common.snps == "rare" ){
        message2("Read in HRC SNPs with EUR+EAS MAC >= 10")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eur_eas_snps_mac_ge_10.tsv")
      }
    
    # Population: EUR
    }else if( population == "eur" ){
      
      # Common SNPs (MAF >= 1%)
      if( rare.or.common.snps == "common" ){
        message2("Read in HRC SNPs with EUR MAF >= 1%")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eur_snps_maf_ge_0.01.tsv")
        
      # Rare SNPs (MAC >= 10)
      }else if( rare.or.common.snps == "rare" ){
        message2("Read in HRC SNPs with EUR MAC >= 10")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eur_snps_mac_ge_10.tsv")
      }
      
    # Population: EAS
    }else if( population == "eas" ){
      
      # Common SNPs (MAF >= 1%)
      if( rare.or.common.snps == "common" ){
        message2("Read in HRC SNPs with EAS MAF >= 1%")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eas_snps_maf_ge_0.01.tsv")
        
        # Rare SNPs (MAC >= 10)
      }else if( rare.or.common.snps == "rare" ){
        message2("Read in HRC SNPs with EAS MAC >= 10")
        hrc <- fread("/projects/0/prjs0817/projects/pops/data/hrc_eas_snps_mac_ge_10.tsv")
      }
    }
    
  }else if( ld_panel == "g1000" ){
    message2("Read in 1000 Genomes SNPs with EUR MAC >= 10")
    hrc <- fread("/home/heilbron/projects/pops/data/g1000_eur_snps_mac_ge_10.tsv")
  }
  
  # Read in GWAS
  message2("Read in GWAS")
  gw0 <- fread(gw_file)
  
  # Re-name GWAS columns
  message2("Re-name GWAS columns")
  names(gw0)[ names(gw0) == chr_bp_col ] <- "chr_bp"
  names(gw0)[ names(gw0) == chr_col    ] <- "chr"
  names(gw0)[ names(gw0) == bp_col     ] <- "bp"
  names(gw0)[ names(gw0) == a1_col     ] <- "A1"
  names(gw0)[ names(gw0) == a2_col     ] <- "A2"
  names(gw0)[ names(gw0) == b_col      ] <- "b"
  names(gw0)[ names(gw0) == or_col     ] <- "or"
  names(gw0)[ names(gw0) == se_col     ] <- "se"
  names(gw0)[ names(gw0) == p_col      ] <- "p"
  names(gw0)[ names(gw0) == eaf_col    ] <- "freq"
  names(gw0)[ names(gw0) == n1_col     ] <- "n1"
  names(gw0)[ names(gw0) == n0_col     ] <- "n0"
  names(gw0)[ names(gw0) == n_col      ] <- "N"
  
  # If specified, remove SNPs with weak P values to speed up computation
  if( p_threshold < 1 ){
    message2("Remove SNPs with weak P values to speed up computation")
    gw <- gw0[ gw0$p < p_threshold , ]
  }else{
    gw <- gw0
  }
  
  # Make columns for chromosome and position
  if( "chr_bp" %in% names(gw) ){
    message2("Make columns for chromosome and position")
    pattern <- "^chr([[:alnum:]]+)[[:punct:]]([[:digit:]]+)$"
    gw$chr  <- as.integer( sub( pattern     = pattern, 
                                replacement = "\\1", 
                                x=gw$chr_bp ) )
    gw$bp   <- as.integer( sub( pattern     = pattern, 
                                replacement = "\\2", 
                                x=gw$chr_bp ) )
  }
  
  # Create a per-SNP effective sample size column
  #   If this column already exists, do nothing
  #   If an overall study N is provided, use it for all SNPs
  #   Otherwise, compute the effective N from the number of cases and controls
  if( "N" %in% names(gw) ){
    message2("An effective sample size column has been provided and will be used")
  }else if( !is.null(n) ){
    message2("A study-wide effective sample size has been provided and will be ",
             "applied to all SNPs")
    gw$N <- n
  }else if( "n1" %in% names(gw) & "n0" %in% names(gw) ){
    message2("Columns for number of cases and controls have been provided for ",
             "each SNP, computing the effective sample size")
    gw$N <- n_eff( gw$n1, gw$n0 )
  }
  
  # If ORs have been provided, convert to logORs
  if( "or" %in% names(gw) ){
    message2("Convert ORs to logORs")
    gw$b <- log(gw$or)
  }
  
  # Replace SEs to be compatible with betas and P values?
  replace_se <- TRUE
  if(replace_se){
    new_se <- gw$b / p_to_z( p=gw$p, direction=gw$b )
    gw$se <- ifelse( is.na(new_se), gw$se, new_se )
  }
  
  
  #-----------------------------------------------------------------------------
  #   Harmonize GWAS and reference panel alleles
  #-----------------------------------------------------------------------------
  
  # Subset GWAS and reference panel to shared SNPs based on chromosome and position
  message2("Subset GWAS and reference panel to shared SNPs based on chromosome ",
           "and position")
  cpab_gw  <- paste( gw$chr,  gw$bp,  
                     ifelse( gw$A1   < gw$A2,   gw$A1,   gw$A2  ), 
                     ifelse( gw$A1   < gw$A2,   gw$A2,   gw$A1  ),  sep="_" )
  cpab_hrc <- paste( hrc$chr, hrc$bp, 
                     ifelse( hrc$alt < hrc$ref, hrc$alt, hrc$ref ), 
                     ifelse( hrc$alt < hrc$ref, hrc$ref, hrc$alt ), sep="_" )
  cpab_both <- intersect( cpab_hrc, cpab_gw )
  hrc2 <- hrc[ match( cpab_both, cpab_hrc ) , ]
  gw2  <- gw[  match( cpab_both, cpab_gw  ) , ]
  message2( "Of the ", NROW(gw), " GWAS SNPs, ", NROW(gw2), 
            " (", round( 100*NROW(gw2)/NROW(gw), 2 ), 
            "%) were found in the reference panel" )
  message2( "Of the ", NROW(hrc), " reference panel SNPs, ", NROW(hrc2), 
            " (", round( 100*NROW(hrc2)/NROW(hrc), 2 ), 
            "%) were found in the GWAS" )
  
  # Flip SNPs with discordant alleles (flipped in reference v. GWAS)
  message2("Flip SNPs with discordant alleles (flipped in reference v. GWAS)")
  discord <- hrc2$alt != gw2$A1
  original_gwas_A2 <- gw2$A2
  original_gwas_A1 <- gw2$A1
  gw2$A2[discord]   <- original_gwas_A1[discord]
  gw2$A1[discord]   <- original_gwas_A2[discord]
  gw2$b[discord]    <- -1 * gw2$b[discord]
  gw2$freq[discord] <- 1  - gw2$freq[discord]
  
  
  #-----------------------------------------------------------------------------
  #   Wrap up harmonization
  #-----------------------------------------------------------------------------
  
  # Check that CPRA is 100% identical now
  cpra_gw  <- paste( gw2$chr,  gw2$bp,  gw2$A2,   gw2$A1,   sep="_" )
  cpra_hrc <- paste( hrc2$chr, hrc2$bp, hrc2$ref, hrc2$alt, sep="_" )
  if( all( cpra_gw == cpra_hrc ) ){
    message2("CPRA is now identical for all GWAS and reference panel SNPs")
  }else{
    diff_cpra <- head( which( cpra_gw != cpra_hrc ) )
    stop( paste( "Error: not all CPRA are identical for GWAS and reference panel",
                 "SNPs. Here are (up to) the first 6:", diff_cpra, collapse=" " ) )
  }
  
  # Replace GWAS SNP names with reference panel names
  message2("Replace GWAS SNP names with reference panel names")
  gw2$SNP <- hrc2$snp
  
  
  #-----------------------------------------------------------------------------
  #   Format and write outputs
  #-----------------------------------------------------------------------------
  
  # Dump GWAS summary statistics
  # Column names/order: SNP, A1, A2, freq, b, se, p, N, chr, bp
  message2("Dump GWAS summary statistics")
  gw_outfile <- file.path( maindir, "gwas_sumstats.tsv" )
  gw_out <- gw2[ , c( "SNP", "A1", "A2", "freq", "b", "se", "p", "N", "chr", "bp" ) ]
  gw_out <- gw_out[ order( gw_out$chr, gw_out$bp ) , ]
  fwrite( x=gw_out, file=gw_outfile, sep="\t" )
  
  
  #-----------------------------------------------------------------------------
  #   Done
  #-----------------------------------------------------------------------------
}


#-------------------------------------------------------------------------------
#   qc_gwas:                 Write GWAS sumstats in COJO format
#-------------------------------------------------------------------------------

qc_gwas <- function( maindir, ld_panel ){
  
  #-----------------------------------------------------------------------------
  #   Read in GWAS and HRC
  #-----------------------------------------------------------------------------
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  library(data.table)
  
  # Read in reference panel SNPs
  if( ld_panel == "hrc" ){
    rare.or.common.snps <- "rare"
    if( rare.or.common.snps == "common"){
      message2("Read in HRC SNPs with EUR MAF >= 1%")
      hrc <- fread("/home/heilbron/projects/pops/data/hrc_eur_eas_snps_maf_ge_0.01.tsv")
    }else if( rare.or.common.snps == "rare"){
      message2("Read in HRC SNPs with EUR MAC >= 10")
      hrc <- fread("/home/heilbron/projects/pops/data/hrc_eur_eas_snps_mac_ge_10.tsv")
    }else{
      stop("rare.or.common must be 'rare' or 'common'")
    }
  }else if( ld_panel == "g1000" ){
    message2("Read in 1000 Genomes SNPs with EUR MAC >= 10")
    hrc <- fread("/home/heilbron/projects/pops/data/g1000_eur_snps_mac_ge_10.tsv")
  }else{
    stop("ld_panel must be either 'hrc' or 'g1000'")
  }
  
  # Read in GWAS
  message2("Read in GWAS")
  gw_file <- file.path( maindir, "gwas_sumstats.tsv" )
  gw <- fread(gw_file)
  
  
  #-----------------------------------------------------------------------------
  #   Harmonize GWAS and reference panel alleles
  #-----------------------------------------------------------------------------
  
  # Subset GWAS and reference panel to shared SNPs
  message2("Subset GWAS and reference panel to shared SNPs")
  snps_both <- intersect( hrc$snp, gw$SNP )
  hrc2 <- hrc[ match( snps_both, hrc$snp ) , ]
  gw2  <- gw[  match( snps_both, gw$SNP  ) , ]
  message2( "Of the ", NROW(gw), " GWAS SNPs, ", NROW(gw2), 
            " (", round( 100*NROW(gw2)/NROW(gw), 2 ), 
            "%) were found in the reference panel" )
  message2( "Of the ", NROW(hrc), " reference panel SNPs, ", NROW(hrc2), 
            " (", round( 100*NROW(hrc2)/NROW(hrc), 2 ), 
            "%) were found in the GWAS" )
  
  # Remove SNPs with a large AF difference between the GWAS and reference
  big_af_diff  <- 0.1
  big_af_ratio <- 3
  abs_af   <- abs( gw2$freq - hrc2$af )               > big_af_diff
  ratio_af <- exp( abs( log( gw2$freq / hrc2$af ) ) ) > big_af_ratio
  diff_af  <- abs_af | ratio_af
  message2( sum(diff_af), "/", NROW(hrc2), " SNPs (", 
            round( 100 * sum(diff_af) / NROW(hrc2), 2 ), 
            "%) have an absolute AF difference > ", big_af_diff, 
            " or an absolute AF ratio > ", big_af_ratio,
            ": flagging for removal" )
  
  # Report rare compatible SNPs
  common_maf <- 0.4
  comm <- ifelse( hrc2$af <= 0.5, 
                  hrc2$af         > common_maf, 
                  ( 1 - hrc2$af ) > common_maf )
  rare_comp <- !diff_af & !comm
  message2( sum(rare_comp), "/", NROW(hrc2), " SNPs (", 
            round( 100 * sum(rare_comp) / NROW(hrc2), 2 ), 
            "%) have similar AFs and MAF <= ", common_maf, ": no action" )
  
  # Report non-palindromic common compatible SNPs
  pal <- ( hrc2$alt=="A" & hrc2$ref=="T" ) | 
    ( hrc2$alt=="T" & hrc2$ref=="A" ) | 
    ( hrc2$alt=="C" & hrc2$ref=="G" ) | 
    ( hrc2$alt=="G" & hrc2$ref=="C" )
  ncc <- !pal & comm & !diff_af
  message2( sum(ncc), "/", NROW(hrc2), " SNPs (", 
            round( 100 * sum(ncc) / NROW(hrc2), 2 ), 
            "%) have similar AFs, MAF > ", common_maf, 
            ", and are non-palindromic: no action" )
  
  # Prepare to report palindromic common compatible SNPs
  pcc <- pal & comm & !diff_af
  base_msg <- paste0( sum(pcc), "/", NROW(hrc2), " SNPs (", 
                      round( 100 * sum(pcc) / NROW(hrc2), 2 ), 
                      "%) have similar AFs, MAF > ", common_maf, 
                      ", and are palindromic: " )
  
  # If using DENTIST, keep common palindromic SNPs
  # If not using DENTIST, remove common palindromic SNPs
  dentist <- TRUE
  if(dentist){
    message2( base_msg, "no action (DENTIST will be used for further SNP QC)")
    bad_snps <- diff_af
  }else{
    message2( base_msg, "flagging for removal")
    bad_snps <- diff_af | pcc
  }
  
  # Remove flagged SNPs
  message2("Remove flagged SNPs")
  gw3  <- gw2[  !bad_snps , ]
  hrc3 <- hrc2[ !bad_snps , ]
  
  
  #-----------------------------------------------------------------------------
  #   Format and write outputs
  #-----------------------------------------------------------------------------
  
  # Dump GWAS summary statistics
  message2("Dump GWAS summary statistics")
  gw_outfile <- file.path( maindir, "gwas_sumstats.qc.tsv" )
  fwrite( x=gw3, file=gw_outfile, sep="\t" )
  
  
  #-----------------------------------------------------------------------------
  #   Done
  #-----------------------------------------------------------------------------
}


#-------------------------------------------------------------------------------
#   define_loci_clump:       Clump SNPs, add buffer, merge overlapping regions
#-------------------------------------------------------------------------------

define_loci_clump <- function( maindir, population, do.test=FALSE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  suppressMessages( library(IRanges) )
  
  # If output file already exists, skip
  clmp_loci_file <- file.path( maindir, "loci_clumped.tsv" )
  if( file.exists(clmp_loci_file) ){
    message2("Clumped loci file already exists, skipping")
    return()
  }
  
  # Create a directory for per-chromosome clumping outputs
  message2("Create a directory for per-chromosome clumping outputs")
  clump_dir <- file.path( maindir, "clumps" )
  dir.create( clump_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Set static clumping arguments
  message2("Set static clumping arguments")
  plink    <- file.path( "/projects/0/prjs0817/software/plink/plink" )
  gw_file  <- file.path( maindir, "gwas_sumstats.tsv" )
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
  # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
  # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  
  #-----------------------------------------------------------------------------
  #   Loop through chromosomes and clump, then collate results
  #-----------------------------------------------------------------------------
  
  chromosomes <- 1:22
  if(do.test) chromosomes <- 22
  for( i in chromosomes ){
    
    # Set per-chromosome clumping arguments
    bed_file <- file.path( bed_dir, paste0( bed_pref, i, bed_suff ) )
    out_pre  <- file.path( clump_dir, paste0( "chr", i ) )
    outfile  <- paste0( out_pre, ".clumped" )
    
    # If the output file already exists, skip
    if( file.exists(outfile) ){
      message2( "Output file for chromosome ", i, " already exists, skipping" )
      next
    }else{
      message2( "Clumping chromosome: ", i )
    }
    
    # Run clumping
    cmd <- paste( plink, 
                  "--bfile",          bed_file, 
                  "--clump",          gw_file,
                  "--clump-p1",       5e-8,
                  "--clump-field",    "p",
                  "--clump-r2",       0.1,
                  "--clump-kb",       3000,
                  "--out",            out_pre )
    system( command=cmd, intern=TRUE )
  }
  
  # Collate results across chromosomes
  message2("Collate results across chromosomes")
  clump_files <- list.files( path=clump_dir, pattern=".clumped$", full.names=TRUE )
  clmp0 <- lapply( clump_files, fread )
  clmp  <- do.call( rbind, clmp0 )
  clmp$hi.pos <- clmp$lo.pos <- NA
  clmp <- clmp[ order( clmp$CHR, clmp$BP ) , ]
  
  
  #-----------------------------------------------------------------------------
  #   Loop through LD-independent variants and define loci
  #-----------------------------------------------------------------------------
  
  # Read in GWAS sumstats
  message2("Read in GWAS sumstats")
  gw <- fread(gw_file)
  
  # For each LD-independent variant, define locus boundaries
  message2("For each LD-independent variant, define locus boundaries")
  for( i in seq_len( NROW(clmp) ) ){
    
    # Subset
    message2( "Starting SNP ", i, "/", NROW(clmp) )
    sub <- clmp[i,]
    
    # Split LD friends, remove the suffix, add the chosen SNP
    ldf  <- strsplit( x=sub$SP2, split="," )[[1]]
    ldf2 <- sub( pattern="\\([[:digit:]]+\\)$", replacement="", x=ldf )
    ldf3 <- c( sub$SNP, ldf2 )
    ldf4 <- setdiff( x=ldf3, y="NONE" )
    
    # Find position of each SNP in the GWAS sumstats
    pos <- gw$bp[ match( ldf4, gw$SNP ) ]
    
    # Assign lo.pos and hi.pos, including 500kb buffer
    clmp$lo.pos[i] <- min(pos) - 5e5
    clmp$hi.pos[i] <- max(pos) + 5e5
    if( clmp$lo.pos[i] < 1 ) clmp$lo.pos[i] <- 1
  }
  
  
  #-----------------------------------------------------------------------------
  #   Merge overlapping loci and write to file
  #-----------------------------------------------------------------------------
  
  # Merge overlapping loci
  message2("Merge overlapping loci")
  gcols <- c( "SNP", "P", "CHR", "BP", "lo.pos", "hi.pos" )
  clmp2 <- clmp[ , ..gcols ]
  clmp3 <- list()
  for( i in unique(clmp2$CHR) ){
    ir  <- IRanges( start = clmp2$lo.pos[ clmp2$CHR == i ], 
                    end   = clmp2$hi.pos[ clmp2$CHR == i ] )
    ir2 <- reduce(ir)
    clmp3[[i]] <- data.table( chr=i, lo.pos=ir2@start, 
                              hi.pos=ir2@start + ir2@width - 1 )
  }
  clmp4 <- do.call( rbind, clmp3 )
  
  # Write merged loci to file
  message2("Write merged loci to file")
  fwrite( x=clmp4, file=clmp_loci_file, sep="\t" )
}


#-------------------------------------------------------------------------------
#   dentist:                 Run DENTIST and remove SNPs that fail
#-------------------------------------------------------------------------------

dentist <- function( maindir, do.test=FALSE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  suppressMessages( library(parallel) )
  
  # Create a directory to DENTIST results
  message2("Create a directory to DENTIST results")
  dent_dir <- file.path( maindir, "dentist" )
  dir.create( path=dent_dir, showWarnings=FALSE, recursive=TRUE )
  
  # If output file already exists, skip
  dent_outfile <- file.path( maindir, "gwas_sumstats.dentist.tsv" )
  if( file.exists(dent_outfile) ){
    message2("DENTIST-QC'd GWAS sumstats file already exists, skipping")
    return()
  }
  
  # Set global parameters
  dentist_binary <- "/projects/0/prjs0817/software/DENTIST/DENTIST_1.3.0.0"
  n_threads      <- detectCores()
  
  # Read in the GWAS
  message2("Read in the GWAS")
  gw_file <- file.path( maindir, "gwas_sumstats.qc.tsv" )
  gw      <- fread(gw_file)
  
  # Write out a version of the GWAS with no chromosome and position columns
  message2("Write out a version of the GWAS with no chromosome and position columns")
  gcols <- setdiff( names(gw), c( "chr", "bp" ) )
  gw2 <- gw[ , ..gcols ]
  dent_in_ss_file <- file.path( dent_dir, "gwas_sumstats.qc.tsv" )
  fwrite( x=gw2, file=dent_in_ss_file, sep="\t" )
  
  # Read in loci
  message2("Read in loci")
  loci_file <- file.path( maindir, "loci_clumped.tsv" )
  loci <- fread(loci_file)
  
  
  #-----------------------------------------------------------------------------
  #   Loop through loci
  #-----------------------------------------------------------------------------
  
  locus_idx <- seq_along(loci$chr)
  if(do.test) locus_idx <- locus_idx[ loci$chr == 22 ]
  for( i in locus_idx ){
    
    # Extract locus information
    chr        <- loci$chr[i]
    lo         <- loci$lo.pos[i]
    hi         <- loci$hi.pos[i]
    
    # If output file already exists, skip
    out_name <- paste0( "chr", chr, "_", lo, "_", hi )
    out_pre  <- file.path( dent_dir, out_name )
    outfile  <- paste0( out_pre, ".DENTIST.full.txt" )
    if( file.exists(outfile) ){
      message2("DENTIST file already exists for ", out_name, ", skipping")
      next
    }else{
      message2( "Running DENTIST on locus ", i, "/", NROW(loci), ": ",
                out_name )
    }
    
    # Extract local SNP names
    idx  <- gw$chr == chr & gw$bp >= lo & gw$bp <= hi
    snps <- gw$SNP[idx]
    
    # Write local SNPs to file
    locus_name <- paste0( "chr", chr, "_", lo, "_", hi )
    snp_file   <- file.path( dent_dir, paste0( locus_name, ".snplist" ) )
    writeLines( text=snps, con=snp_file )
    
    # Set BED file path
    bed_suff <- paste0( "HRC.r1-1.EGA.GRCh37.chr", chr, 
                        ".impute.plink.combined.EUR.2191.EAS.538" )
    bed_file <- file.path( "/gpfs/work5/0/pgcdac/imputation_references/",
                           "HRC.r1-1_merged_EUR_EAS_panel/", bed_suff )
    
    # Run DENTIST
    # message2( "Running DENTIST on locus ", i, "/", NROW(loci), ": ", outname )
    cmd <- paste( dentist_binary,
                  "--gwas-summary", dent_in_ss_file,
                  "--bfile",        bed_file,
                  "--chrID",        chr,
                  "--extract",      snp_file,
                  "--thread-num",   n_threads,
                  "--out",          out_pre,
                  "--no-missing-genotype" )
    system( command=cmd, intern=FALSE )
  }
  
  # Store SNPs that failed DENTIST
  fail_files <- list.files( path=dent_dir, pattern=".DENTIST.short.txt$", 
                            full.names=TRUE )
  failed0 <- lapply( fail_files, readLines )
  failed  <- unlist(failed0)
  
  # Remove SNPs that fail DENTIST
  gw_post <- gw[ -match( failed, gw$SNP ) , ]
  
  # Write QC'd GWAS to file
  fwrite( x=gw_post, file=dent_outfile, sep="\t" )
}


#-------------------------------------------------------------------------------
#   run_cojo_genome:         Get independent hits across the whole genome
#-------------------------------------------------------------------------------

run_cojo_genome <- function( maindir, population="eur_eas" ){
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # If output file already exists, skip
  jma_outfile <- file.path( maindir, "hits_genome.jma.cojo" )
  if( file.exists(jma_outfile) ){
    message2("Locus-agnostic hits file already exists, skipping")
    return()
  }
  
  # Set static COJO arguments
  message2("Set static COJO arguments")
  gcta_binary    <- file.path( "/projects/0/prjs0817/software/gcta/",
                               "gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" )
  cojo_file      <- file.path( maindir, "gwas_sumstats.tsv" )
  cojo_p         <- 5e-8
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
    # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
    # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  # Create a directory for per-chromosome COJO outputs
  message2("Create a directory for per-chromosome COJO outputs")
  hits_dir <- file.path( maindir, "hits_genome" )
  dir.create( hits_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Loop through chromosomes
  for( i in 1:22 ){
    
    # Set per-chromosome COJO arguments
    bed_file <- file.path( bed_dir, paste0( bed_pref, i, bed_suff ) )
    out_pre  <- file.path( hits_dir, paste0( "chr", i ) )
    outfile  <- paste0( out_pre, ".jma.cojo" )
    
    # If the output file already exists, skip
    if( file.exists(outfile) ){
      message2( "Output file for chromosome: ", i, " already exists, skipping" )
      next
    }
    
    # Run COJO
    message2( "Running COJO on chromosome: ", i )
    cmd <- paste( gcta_binary, 
                  "--bfile",          bed_file, 
                  "--cojo-file",      cojo_file,
                  "--cojo-slct",
                  "--cojo-p",         cojo_p,
                  "--cojo-collinear", 0.9,
                  "--out",            out_pre )
    system(cmd)
  }
  
  # Collate results across chromosomes
  message2("Collate results across chromosomes")
  jma_files <- list.files( path=hits_dir, pattern=".jma.cojo$", full.names=TRUE )
  jma0 <- lapply( X=jma_files, FUN=fread )
  jma  <- do.call( rbind, jma0 )
  jma  <- jma[ order( jma$Chr, jma$bp ) , ]
  
  # Write results to file
  message2("Write results to file")
  fwrite( x=jma, file=jma_outfile, sep="\t" )
}


#-------------------------------------------------------------------------------
#   run_cojo_local:          Get independent hits in each clumping-defined locus
#-------------------------------------------------------------------------------

run_cojo_local <- function( maindir, population, do.test=FALSE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries and sources
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # If output file already exists, skip
  jma_outfile <- file.path( maindir, "hits_locus.jma.cojo" )
  if( file.exists(jma_outfile) ){
    message2("Locus-agnostic hits file already exists, skipping")
    return()
  }
  
  # Set static COJO arguments
  message2("Set static COJO arguments")
  gcta_binary <- file.path( "/projects/0/prjs0817/software/gcta/",
                            "gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" )
  cojo_file   <- file.path( maindir, "gwas_sumstats.tsv" )
  
  # Create a directory for per-locus COJO outputs
  message2("Create a directory for per-locus COJO outputs")
  hits_dir <- file.path( maindir, "hits_locus" )
  dir.create( hits_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
    # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
    # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  # Read in the GWAS
  message2("Read in the GWAS")
  gw <- fread(cojo_file)
  
  # Read in loci
  message2("Read in loci")
  loci_file <- file.path( maindir, "loci_clumped.tsv" )
  loci <- fread(loci_file)
  
  
  #-----------------------------------------------------------------------------
  #   Loop through loci
  #-----------------------------------------------------------------------------
  
  locus_idx <- seq_along(loci$chr)
  if(do.test) locus_idx <- locus_idx[ loci$chr == 22 ]
  for( i in locus_idx ){
    
    # Set per-locus COJO arguments
    chr      <- loci$chr[i]
    lo       <- loci$lo.pos[i]
    hi       <- loci$hi.pos[i]
    bed_file <- file.path( bed_dir, paste0( bed_pref, chr, bed_suff ) )
    outname  <- paste0( "chr", chr, "_", lo, "_", hi )
    out_pre  <- file.path( hits_dir, outname )
    outfile  <- paste0( out_pre, ".jma.cojo" )
    
    # If the output file already exists, skip
    if( file.exists(outfile) ){
      message2( "Output file for: ", outname, " already exists, skipping" )
      next
    }
    
    # Write to file all SNPs in the locus
    snp_file <- sub( pattern=".jma.cojo$", replacement=".snplist", x=outfile )
    idx <- gw$chr == chr & gw$bp >= lo & gw$bp <= hi
    writeLines( text=gw$SNP[idx], con=snp_file )
    
    # Run COJO
    message2( "Running COJO on locus ", i, "/", NROW(loci), ": ", outname )
    cmd <- paste( gcta_binary,
                  "--bfile",     bed_file,
                  "--cojo-file", cojo_file,
                  "--chr",       chr,
                  "--extract",   snp_file,
                  "--cojo-slct",
                  "--out",       out_pre )
    system( command=cmd, intern=TRUE )
  }
  
  # Collate results across chromosomes
  message2("Collate results across chromosomes")
  jma_files <- list.files( path=hits_dir, pattern=".jma.cojo$", full.names=TRUE )
  jma0 <- lapply( X=jma_files, FUN=fread )
  jma  <- do.call( rbind, jma0 )
  jma  <- jma[ order( jma$Chr, jma$bp ) , ]
  
  # Write results to file
  message2("Write results to file")
  fwrite( x=jma, file=jma_outfile, sep="\t" )
}


#-------------------------------------------------------------------------------
#   run_cojo_cluster:        Get independent hits in each clumping-defined locus
#-------------------------------------------------------------------------------

run_cojo_cluster <- function( maindir, population="eur_eas", do.test=FALSE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries and sources
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # If output file already exists, skip
  jma_outfile <- file.path( maindir, "hits_locus.jma.cojo" )
  if( file.exists(jma_outfile) ){
    message2("Locus-agnostic hits file already exists, skipping")
    return()
  }
  
  # Set static COJO arguments
  message2("Set static COJO arguments")
  gcta_binary <- file.path( "/projects/0/prjs0817/software/gcta/",
                            "gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" )
  cojo_file   <- file.path( maindir, "gwas_sumstats.tsv" )
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
    # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
    # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  # Create a directory for per-locus COJO outputs
  message2("Create a directory for per-locus COJO outputs")
  hits_dir <- file.path( maindir, "hits_locus" )
  dir.create( hits_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Read in the GWAS
  message2("Read in the GWAS")
  gw <- fread(cojo_file)
  
  # Read in loci
  message2("Read in loci")
  loci_file <- file.path( maindir, "loci_clumped.tsv" )
  loci <- fread(loci_file)
  
  
  #-----------------------------------------------------------------------------
  #   Loop through loci
  #-----------------------------------------------------------------------------
  
  # Create a job identifier
  job.id <- paste0( "c", sample( x=1:999, size=1 ) )
  
  # Make a log directory
  logdir <- file.path( maindir, "log", "run_cojo" )
  dir.create( path=logdir, showWarnings=FALSE, recursive=TRUE )
  
  # Loop through loci
  locus_idx <- seq_along(loci$chr)
  if(do.test) locus_idx <- locus_idx[ loci$chr == 22 ]
  for( i in locus_idx ){
    
    # Set per-locus COJO arguments
    chr      <- loci$chr[i]
    lo       <- loci$lo.pos[i]
    hi       <- loci$hi.pos[i]
    bed_file <- file.path( bed_dir, paste0( bed_pref, chr, bed_suff ) )
    outname  <- paste0( "chr", chr, "_", lo, "_", hi )
    out_pre  <- file.path( hits_dir, outname )
    outfile  <- paste0( out_pre, ".jma.cojo" )
    
    # If the output file already exists, skip
    if( file.exists(outfile) ){
      message2( "Output file for: ", outname, " already exists, skipping" )
      next
    }
    
    # Write to file all SNPs in the locus
    snp_file <- sub( pattern=".jma.cojo$", replacement=".snplist", x=outfile )
    idx <- gw$chr == chr & gw$bp >= lo & gw$bp <= hi
    writeLines( text=gw$SNP[idx], con=snp_file )
    
    # Create job name and log file name
    jobname <- paste0( job.id, ".", i )
    logfile <- file.path( logdir, paste0( outname, ".log" ) )
    
    # Run COJO
    message2( "Running COJO on locus ", i, "/", NROW(loci), ": ", outname )
    cmd <- paste( gcta_binary,
                  "--bfile",     bed_file,
                  "--cojo-file", cojo_file,
                  "--chr",       chr,
                  "--extract",   snp_file,
                  "--cojo-slct",
                  "--out",       out_pre )
    cmd2 <- paste0( '"', cmd, '"' )
    cmd3 <- paste( "sbatch",
                   "-J", jobname,
                   "-o", logfile,
                   "-e", logfile,
                   "/projects/0/prjs0817/repos/cojo_pipe/b_run_cojo.sh",
                   cmd2 )
    # cmd3 <- paste( "sbatch",
    #                "-J", jobname,
    #                "-o", logfile,
    #                "-e", logfile,
    #                "/projects/0/prjs0817/repos/cojo_pipe/b_run_cojo.sh",
    #                bed_file, cojo_file, chr, snp_file, out_pre )
    system( command=cmd3, intern=FALSE )
  }
  
  # Wait until all jobs are finished before proceeding
  message2("Wait until all jobs are finished before proceeding")
  ensure_finished_jobs(job.id)
  
  # Check that all jobs completed successfully
  message2("Check that all jobs completed successfully")
  jma_files <- list.files( path=hits_dir, pattern=".jma.cojo$", full.names=TRUE )
  if( length(jma_files) != length(locus_idx) ){
    all_outnames <- paste0( "chr", loci$chr, "_", loci$lo.pos, "_", 
                            loci$hi.pos, ".jma.cojo" )
    all_outfiles <- file.path( hits_dir, all_outnames )
    outnames_miss <- all_outnames[ -match( jma_files, all_outfiles ) ]
    stop( "The following jobs failed: ", paste( outnames_miss, collapse=", " ) )
  }
  
  # Collate results across loci
  message2("Collate results across loci")
  jma_files <- list.files( path=hits_dir, pattern=".jma.cojo$", full.names=TRUE )
  jma0 <- lapply( X=jma_files, FUN=fread )
  jma  <- do.call( rbind, jma0 )
  jma  <- jma[ order( jma$Chr, jma$bp ) , ]
  
  # Write results to file
  message2("Write results to file")
  fwrite( x=jma, file=jma_outfile, sep="\t" )
}


#-------------------------------------------------------------------------------
#   rm_weak_hits:            Get independent hits in each clumping-defined locus
#-------------------------------------------------------------------------------

rm_weak_hits <- function(maindir){
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries and sources
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # If output file already exists, skip
  weak_file <- file.path( maindir, "hits_locus_no_weak.jma.cojo" )
  if( file.exists(weak_file) ){
    message2("Weak hits have already been removed, skipping")
    return()
  }
  
  # Read in hits
  message2("Read in hits")
  hits_file <- file.path( maindir, "hits_locus.jma.cojo" )
  hits <- fread(hits_file)
  
  # Remove weak hits
  message2("Remove weak hits")
  hits2 <- hits[ hits$p < 5e-8 , ]
  
  # Write pruned hits to file
  message2("Write pruned hits to file")
  fwrite( x=hits2, file=weak_file, sep="\t" )
}


#-------------------------------------------------------------------------------
#   isolate_signals_local:   Extract conditioned sumstats for each independent hit
#-------------------------------------------------------------------------------

isolate_signals_local <- function( maindir, population="eur_eas", do.test=FALSE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries and sources
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # Set static COJO arguments
  message2("Set static COJO arguments")
  gcta_binary <- file.path( "/projects/0/prjs0817/software/gcta/",
                            "gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" )
  cojo_file   <- file.path( maindir, "gwas_sumstats.tsv" )
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
    # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
    # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  # Create a directory for leave-one-hit-out conditioned sumstats
  message2("Create a directory for leave-one-hit-out conditioned sumstats")
  loho_dir <- file.path( maindir, "loho_conditioned_sumstats" )
  dir.create( loho_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Read in the GWAS
  message2("Read in the GWAS")
  gw <- fread(cojo_file)
  
  # Read in loci
  message2("Read in loci")
  loci_file <- file.path( maindir, "loci_clumped.tsv" )
  loci <- fread(loci_file)
  
  # Read in (marginally significant) COJO hits
  all_hits_file <- file.path( maindir, "hits_locus_no_weak.jma.cojo" )
  all_hits <- fread(all_hits_file)
  
  
  #-----------------------------------------------------------------------------
  #   Loop through loci
  #-----------------------------------------------------------------------------
  
  locus_idx <- seq_along(loci$chr)
  if(do.test) locus_idx <- locus_idx[ loci$chr == 22 ]
  for( i in locus_idx ){
    
    # Subset to local COJO hits
    chr        <- loci$chr[i]
    lo         <- loci$lo.pos[i]
    hi         <- loci$hi.pos[i]
    idx        <- all_hits$Chr == chr & all_hits$bp >= lo & all_hits$bp <= hi
    hits       <- all_hits[ idx , ]
    
    # Find the file containing all SNPs in the locus
    hits_dir   <- file.path( maindir, "hits_locus" )
    locus_name <- paste0( "chr", chr, "_", lo, "_", hi )
    snp_file   <- file.path( hits_dir, paste0( locus_name, ".snplist" ) )
    
    # Set BED file path
    bed_file <- file.path( bed_dir, paste0( bed_pref, chr, bed_suff ) )
    
    # Loop through hits
    for( j in seq_along(hits$SNP) ){
      
      # If the output file already exists, skip
      outname <- paste0( locus_name, "_", j, "_", hits$SNP[j] )
      out_pre <- file.path( loho_dir, outname )
      outfile <- paste0( out_pre, ".cma.cojo" )
      if( file.exists(outfile) ){
        message2( "Output file for locus ", i, "/", NROW(loci), ": ", 
                  locus_name, " already exists, skipping" )
        next
      }else{
        message2( "Starting locus ", i, "/", NROW(loci), ": ", locus_name )
      }
      
      # If there is only one hit, write locus GWAS sumstats
      # If there are multiple hits, sub-loop through hits and run LOHO COJO
      if( NROW(hits) == 1 ){
        
        # Subset GWAS sumstats to the locus, write to file
        message2( "Creating 'conditioned' sumstats for hit ", 
                  j, "/", NROW(hits), ": ", outname )
        idx <- gw$chr == chr & gw$bp >= lo & gw$bp <= hi
        gw2 <- gw[ idx , ]
        fwrite( x=gw2, file=outfile, sep="\t" )
        
      }else{
        
        # Write a file containing all hits in the locus except the focal hit
        loho_snps     <- setdiff( hits$SNP, hits$SNP[j] )
        loho_snp_file <- file.path( loho_dir, paste0( outname, ".txt" ) )
        writeLines( text=loho_snps, con=loho_snp_file )
        
        # Re-run COJO conditioning on the LOHO SNPs
        message2( "Creating conditioned sumstats for hit ", 
                  j, "/", NROW(hits), ": ", outname )
        cmd <- paste( gcta_binary,
                      "--bfile",     bed_file,
                      "--cojo-file", cojo_file,
                      "--chr",       chr,
                      "--extract",   snp_file,
                      "--cojo-cond", loho_snp_file,
                      "--out",       out_pre )
        system( command=cmd, intern=TRUE )
      }
    }
  }
}


#-------------------------------------------------------------------------------
#   isolate_signals_cluster: Extract conditioned sumstats for each independent hit
#-------------------------------------------------------------------------------

isolate_signals_cluster <- function( maindir, population="eur_eas", do.test=FALSE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries and sources
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # Set static COJO arguments
  message2("Set static COJO arguments")
  gcta_binary <- file.path( "/projects/0/prjs0817/software/gcta/",
                            "gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" )
  cojo_file   <- file.path( maindir, "gwas_sumstats.tsv" )
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
    # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
    # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  # Create a directory for leave-one-hit-out conditioned sumstats
  message2("Create a directory for leave-one-hit-out conditioned sumstats")
  loho_dir <- file.path( maindir, "loho_conditioned_sumstats" )
  dir.create( loho_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Read in the GWAS
  message2("Read in the GWAS")
  gw <- fread(cojo_file)
  
  # Read in loci
  message2("Read in loci")
  loci_file <- file.path( maindir, "loci_clumped.tsv" )
  loci <- fread(loci_file)
  
  # Read in (marginally significant) COJO hits
  all_hits_file <- file.path( maindir, "hits_locus_no_weak.jma.cojo" )
  all_hits <- fread(all_hits_file)
  
  
  #-----------------------------------------------------------------------------
  #   Loop through hits
  #-----------------------------------------------------------------------------
  
  # Create a job identifier
  job.id <- paste0( "i", sample( x=1:999, size=1 ) )
  
  # Make a log directory
  logdir <- file.path( maindir, "log", "isolate_signals" )
  dir.create( path=logdir, showWarnings=FALSE, recursive=TRUE )
  
  # Loop through hits
  hits_idx <- seq_along(all_hits$SNP)
  if(do.test) hits_idx <- hits_idx[ all_hits$Chr == 22 ]
  for( i in hits_idx ){
    
    #-----------------------------------------------------------------------------
    #   Check whether the output file already exists
    #-----------------------------------------------------------------------------
    
    # Subset to focal hit
    hit <- all_hits[i,]
    
    # Find the corresponding locus
    locus <- loci[ loci$chr    == hit$Chr & 
                   loci$lo.pos <= hit$bp &
                   loci$hi.pos >= hit$bp , ]
    chr        <- locus$chr
    lo         <- locus$lo.pos
    hi         <- locus$hi.pos
    locus_name <- paste0( "chr", chr, "_", lo, "_", hi )
    
    # Find the hit number
    hits <- all_hits[ all_hits$Chr == chr & 
                      all_hits$bp >= lo & 
                      all_hits$bp <= hi , ]
    hit_num <- which( hits$SNP == hit$SNP )
    
    # If the output file already exists, skip
    outname <- paste0( locus_name, "_", hit_num, "_", hit$SNP )
    out_pre <- file.path( loho_dir, outname )
    outfile <- paste0( out_pre, ".cma.cojo" )
    if( file.exists(outfile) ){
      message2( "Output file for: ", outname, " already exists, skipping" )
      next
    }else{
      message2( "Creating conditioned sumstats for hit ", 
                i, "/", NROW(all_hits), ": ", outname )
    }
    
    
    #-----------------------------------------------------------------------------
    #   Write sumstats or run COJO depending on number of hits in the locus
    #-----------------------------------------------------------------------------
    
    # Find the file containing all SNPs in the locus
    hits_dir   <- file.path( maindir, "hits_locus" )
    snp_file   <- file.path( hits_dir, paste0( locus_name, ".snplist" ) )
    
    # Set BED file path
    bed_file <- file.path( bed_dir, paste0( bed_pref, chr, bed_suff ) )
    
    # If there is only one hit, write locus GWAS sumstats
    # If there are multiple hits, sub-loop through hits and run LOHO COJO
    if( NROW(hits) == 1 ){
      
      # Subset GWAS sumstats to the locus, write to file
      idx <- gw$chr == chr & gw$bp >= lo & gw$bp <= hi
      gw2 <- gw[ idx , ]
      fwrite( x=gw2, file=outfile, sep="\t" )
      
    }else{
      
      # Write a file containing all hits in the locus except the focal hit
      loho_snps     <- setdiff( hits$SNP, hit$SNP )
      loho_snp_file <- file.path( loho_dir, paste0( outname, ".txt" ) )
      writeLines( text=loho_snps, con=loho_snp_file )
      
      # Create job name and log file name
      jobname <- paste0( job.id, ".", i )
      logfile <- file.path( logdir, paste0( outname, ".log" ) )
      
      # Re-run COJO conditioning on the LOHO SNPs
      cmd <- paste( gcta_binary,
                    "--bfile",     bed_file,
                    "--cojo-file", cojo_file,
                    "--chr",       chr,
                    "--extract",   snp_file,
                    "--cojo-cond", loho_snp_file,
                    "--out",       out_pre )
      cmd2 <- paste0( '"', cmd, '"' )
      cmd3 <- paste( "sbatch",
                     "-J", jobname,
                     "-o", logfile,
                     "-e", logfile,
                     "/projects/0/prjs0817/repos/cojo_pipe/c_isolate_signals.sh",
                     cmd2 )
      system( command=cmd3, intern=FALSE )
    }
  }
  
  # Wait until all jobs are finished before proceeding
  message2("Wait until all jobs are finished before proceeding")
  ensure_finished_jobs(job.id)
  
  # Check that all jobs completed successfully
  message2("Check that all jobs completed successfully")
  cma_files <- list.files( path=loho_dir, pattern=".cma.cojo$", full.names=TRUE )
  if( length(cma_files) != length(hits_idx) ){
    pattern <- ".*chr[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_(.*).cma.cojo"
    cma_snps <- sub( pattern=pattern, replacement="\\1", x=cma_files )
    snps_miss <- setdiff( all_hits$SNP, cma_snps )
    stop( "The following hits did not complete successfully: ",
          paste( snps_miss, collapse=", " ) )
  }
}


#-------------------------------------------------------------------------------
#   ld_for_hits:             Compute LD between hits and the rest of the locus
#-------------------------------------------------------------------------------

ld_for_hits <- function( maindir, population ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries and sources
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # Create a directory for LD
  message2("Create a directory for LD")
  ld_dir <- file.path( maindir, "ld" )
  dir.create( ld_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Set static plink arguments
  message2("Set static plink arguments")
  plink    <- file.path( "/projects/0/prjs0817/software/plink/plink" )
  
  # Set reference panel arguments based on population
  message2("Set reference panel arguments based on population")
  bed_pref <- "HRC.r1-1.EGA.GRCh37.chr"
  
  # Merged EUR + EAS
  if( population == "eur_eas" ){
    bed_dir  <- "/gpfs/work5/0/pgcdac/imputation_references/HRC.r1-1_merged_EUR_EAS_panel/"
    bed_suff <- ".impute.plink.combined.EUR.2191.EAS.538"
    
    # Population: EUR
  }else if( population == "eur" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EUR/"
    bed_suff <- ".impute.plink.EUR"
    
    # Population: EAS
  }else if( population == "eas" ){
    bed_dir <- "/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/HRC_reference.r1-1/pop_EAS/"
    bed_suff <- ".impute.plink.EAS"
  }
  
  # Find all conditioned sumstats files
  loho_dir  <- file.path( maindir, "loho_conditioned_sumstats" )
  cma_names <- list.files( path=loho_dir, pattern=".cma.cojo$" )
  cma_files <- file.path( loho_dir, cma_names )
  
  
  #-----------------------------------------------------------------------------
  #   Loop through conditioned sumstats files
  #-----------------------------------------------------------------------------
  
  for( i in seq_along(cma_names) ){
    
    # If output file already exists, skip
    out_name <- sub( pattern=".cma.cojo$", replacement="", x=cma_names[i] )
    out_pre  <- file.path( ld_dir, out_name )
    outfile  <- paste0( out_pre, ".ld" )
    if( file.exists(outfile) ){
      message2("LD file already exists for ", out_name, ", skipping")
      next
    }else{
      message2( "Computing LD for locus ", i, "/", length(cma_names), ": ",
                cma_names[i] )
    }
    
    # Read in conditioned sumstats
    ss <- fread( cma_files[i] )
    
    # Write SNP names to file
    snp_name <- sub( pattern=".cma.cojo$", replacement=".snplist", x=cma_names[i] )
    snp_file <- file.path( ld_dir, snp_name )
    writeLines( text=ss$SNP, con=snp_file )
    
    # Extract hit SNP and chromosome
    pattern <- "chr([[:digit:]]+)_[[:digit:]]+_[[:digit:]]+_[[:digit:]]+_(.*).cma.cojo"
    snp <- sub( pattern=pattern, replacement="\\2", x=cma_names[i] )
    chr <- sub( pattern=pattern, replacement="\\1", x=cma_names[i] )
    
    # Set BED file path
    bed_file <- file.path( bed_dir, paste0( bed_pref, chr, bed_suff ) )
    
    # Run plink to get LD
    cmd <- paste( plink, 
                  "--bfile",          bed_file,
                  "--extract",        snp_file,
                  "--r2",
                  "--ld-snp",         snp,
                  "-ld-window-r2",    0,
                  "--ld-window",      99999,
                  "--ld-window-kb",   99999,
                  "--out",            out_pre )
    system( command=cmd, intern=TRUE )
  }
}


#-------------------------------------------------------------------------------
#   credible_sets:           Compute credible sets for each (conditioned) hit
#-------------------------------------------------------------------------------

credible_sets <- function(maindir){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  suppressMessages( library(coloc) )
  
  # Create a directory to house credible sets
  cs_dir <- file.path( maindir, "cred_sets" )
  dir.create( path=cs_dir, showWarnings=FALSE, recursive=TRUE )
  
  # If output files already exist, skip
  
  # Read in coding and promoter SNPs
  hhc <- fread("/projects/0/prjs0817/projects/cojo/data/high_h2_coding_SNPs.tsv")
  hhp <- fread("/projects/0/prjs0817/projects/cojo/data/high_h2_promoter_SNPs.tsv")
  
  # Find all conditioned sumstats files
  loho_dir  <- file.path( maindir, "loho_conditioned_sumstats" )
  cma_names <- list.files( path=loho_dir, pattern=".cma.cojo$" )
  cma_files <- file.path( loho_dir, cma_names )
  
  
  #-----------------------------------------------------------------------------
  #   Loop through conditioned sumstats files
  #-----------------------------------------------------------------------------
  
  for( i in seq_along(cma_names) ){
    
    # If output file already exists, skip
    out_name <- sub( pattern=".cma.cojo$", replacement="", x=cma_names[i] )
    out_pre  <- file.path( cs_dir, out_name )
    outfile  <- paste0( out_pre, ".cs" )
    if( file.exists(outfile) ){
      message2("Credible set file already exists for ", out_name, ", skipping")
      next
    }else{
      message2( "Computing credible set for locus ", i, "/", 
                length(cma_names), ": ", out_name )
    }
    
    # Read in conditioned sumstats
    ss <- fread( cma_files[i] )
    
    # If sumstats have been conditioned, put into the unconditioned format
    if( "pC" %in% names(ss) ){
      ss$b <- ss$se <- ss$p <- NULL
      names(ss)[ names(ss) == "bC" ]    <- "b"
      names(ss)[ names(ss) == "bC_se" ] <- "se"
      names(ss)[ names(ss) == "pC" ]    <- "p"
      names(ss)[ names(ss) == "n" ]     <- "N"
      names(ss)[ names(ss) == "Chr" ]   <- "chr"
      ss <- ss[ !is.na(ss$b) , ]
    }
    
    # Compute credible set
    dataset <- list( snp     = ss$SNP, 
                     beta    = ss$b, 
                     varbeta = ss$se^2, 
                     N       = median(ss$N),
                     s       = 0.5,
                     type    = "cc" )
    cs <- finemap.abf( dataset=dataset )
    
    # Add PIP to the conditioned sumstats
    ss$pip <- cs$SNP.PP[ -NROW(cs) ]
    
    # Subset to the 95% CS
    ss2 <- ss[ order(-ss$pip) , ]
    ss2$sum <- cumsum(ss2$pip)
    idx <- which( ss2$sum >= 0.95 )
    if( length(idx) == 0 ){
      ss3 <- ss2
    }else{
      ss3 <- head( ss2, idx[1] )
    }
    
    # Annotate CS with whether SNPs are coding or promoter
    ss3$coding   <- hhc$gene[   match( ss3$SNP, hhc$SNP ) ]
    ss3$promoter <- hhp$gene[   match( ss3$SNP, hhp$SNP ) ]
    ss3$ensgid_c <- hhc$ensgid[ match( ss3$SNP, hhc$SNP ) ]
    ss3$ensgid_p <- hhp$ensgid[ match( ss3$SNP, hhp$SNP ) ]
    
    # Write sumstats for credible set SNPs to file
    fwrite( x=ss3, file=outfile, sep="\t" )
  }
}


#-------------------------------------------------------------------------------
#   define_loci_cs:          Assign locus boundaries based on credible sets
#-------------------------------------------------------------------------------

define_loci_cs <- function(maindir){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  
  # If output files already exist, skip
  cs_loci_file <- file.path( maindir, "loci_cs.tsv" )
  if( file.exists(cs_loci_file) ){
    message2("Credible set-based locus file already exists, skipping")
    stop()
  }
  
  # Read in gene locations
  message2("Read in gene locations")
  genes <- fread("/projects/0/prjs0817/projects/pops/data/gene_locations.tsv")
  
  # Find all credible set files
  message2("Find all credible set files")
  cs_dir  <- file.path( maindir, "cred_sets" )
  cs_names <- list.files( path=cs_dir, pattern=".cs$" )
  cs_files <- file.path( cs_dir, cs_names )
  
  
  #-----------------------------------------------------------------------------
  #   Loop through credible set files
  #-----------------------------------------------------------------------------
  
  message2("Loop through credible set files")
  cs_loci0 <- list()
  for( i in seq_along(cs_names) ){
    
    # Read in CS
    cs <- fread( cs_files[i] )
    
    # Extract the hit name
    hit_name <- sub( pattern=".cs$", replacement="", x=cs_names[i])
    
    # Define CS-based locus boundaries
    window <- 3e5
    chr <- unique(cs$chr)
    lo  <- min(cs$bp) - window
    hi  <- max(cs$bp) + window
    p   <- min(cs$p)
    
    # Find the PIP-weighted central position
    bp_mid <- round( sum( cs$bp * cs$pip / sum(cs$pip) ) )
    
    # Find genes affected by coding/promoter SNPs
    c_genes <- paste( unique( cs$coding[   cs$coding   != "" ] ), collapse=";" )
    p_genes <- paste( unique( cs$promoter[ cs$promoter != "" ] ), collapse=";" )
    
    # Find the nearest gene
    genes2 <- genes[ genes$CHR == chr , ]
    genes2$dist_start <- abs( genes2$START - bp_mid )
    genes2$dist_stop  <- abs( genes2$END   - bp_mid )
    genes2$min_dist   <- ifelse( genes2$dist_start < genes2$dist_stop,
                                 genes2$dist_start,  genes2$dist_stop )
    genes2$min_dist   <- ifelse( bp_mid < genes2$END & 
                                 bp_mid > genes2$START,
                                 0, genes2$min_dist )
    min_dist <- min(genes2$min_dist)
    ng_idx <- which( genes2$min_dist == min_dist )
    ng <- paste( genes2$NAME[  ng_idx], collapse="; " )
    
    # Stick all results into the list
    cs_loci0[[i]] <- data.table( hit=hit_name, chr=chr, centre=bp_mid, lo=lo, 
                                 hi=hi, nearest=ng, dist=min_dist, p=p, 
                                 coding=c_genes, promoter=p_genes )
  }
  
  
  #-----------------------------------------------------------------------------
  #   Write credible set-based loci to file
  #-----------------------------------------------------------------------------
  
  message2("Write credible set-based loci to file")
  cs_loci <- do.call( rbind, cs_loci0 )
  cs_loci <- cs_loci[ order( cs_loci$chr, cs_loci$centre ) , ]
  cs_loci$coding[   cs_loci$coding   == "NA" ] <- NA
  cs_loci$promoter[ cs_loci$promoter == "NA" ] <- NA
  fwrite( x=cs_loci, file=cs_loci_file, sep="\t" )
}


#-------------------------------------------------------------------------------
#   lz_plot:                 Create LocusZoom plots
#-------------------------------------------------------------------------------

lz_plot <- function( maindir, cond.or.uncond, merge.loci=TRUE ){
  
  #-----------------------------------------------------------------------------
  #   Get set up
  #-----------------------------------------------------------------------------
  
  # message2: Just like message, but with the date and a gap pasted in
  message2 <- function(...) message( date(), "     ", ...)
  
  # Load libraries
  # source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  message2("Load libraries and sources")
  suppressMessages( library(data.table) )
  suppressMessages( library(locuszoomr) )
  suppressMessages( library(EnsDb.Hsapiens.v75) )
  
  # Create a directory to house LZPs
  message2("Create a directory to house LZPs")
  if(merge.loci){
    dir1 <- "merged"
  }else{
    dir1 <- "separated"
  }
  if( cond.or.uncond == "uncond" ){
    dir2 <- "unconditioned"
  }else if( cond.or.uncond == "cond" ){
    dir2 <- "conditioned"
  }
  lzp_dir <- file.path( maindir, "lzp", dir1, dir2 )
  dir.create( path=lzp_dir, showWarnings=FALSE, recursive=TRUE )
  
  # Read in the GWAS
  if( cond.or.uncond == "uncond" ){
    message2("Read in the GWAS")
    gw_file <- file.path( maindir, "gwas_sumstats.tsv" )
    gw <- fread(gw_file)
  }
  
  # Read in credible set-based loci
  if( !merge.loci ){
    message2("Read in credible set-based loci")
    cs_file <- file.path( maindir, "loci_cs.tsv" )
    cs <- fread(cs_file)
  }
  
  # Read in coding and promoter SNPs
  hhc <- fread("/projects/0/prjs0817/projects/cojo/data/high_h2_coding_SNPs.tsv")
  hhp <- fread("/projects/0/prjs0817/projects/cojo/data/high_h2_promoter_SNPs.tsv")
  
  # Find all conditioned sumstats files
  loho_dir <- file.path( maindir, "loho_conditioned_sumstats" )
  ss_cond_names <- list.files( path=loho_dir, pattern=".cma.cojo$" )
  ss_cond_files <- file.path( loho_dir, ss_cond_names )
  
  
  #-----------------------------------------------------------------------------
  #   Loop through conditioned sumstats files (even if we want unconditioned)
  #-----------------------------------------------------------------------------
  
  for( i in seq_along(ss_cond_files) ){
    
    #-----------------------------------------------------------------------------
    #   Check if output exists, extract locus info from file name
    #-----------------------------------------------------------------------------
    
    # If the output already exists, skip
    lzp_name <- sub( pattern = ".cma.cojo$", replacement = ".jpg", 
                     x       = ss_cond_names[i] )
    lzp_out  <- file.path( lzp_dir, lzp_name )
    if( file.exists(lzp_out) ){
      message2( "LocusZoom plot: ", lzp_name, " already exists, skipping" )
      next
    }else{
      message2( "Making LocusZoom plot ", i, "/", length(ss_cond_files), 
                ": ", lzp_name )
    }
    
    # Extract locus boundaries
    pattern <- "^(chr[[:digit:]]+)_([[:digit:]]+)_([[:digit:]]+)_[[:digit:]]+_(.*).jpg$"
    chr       <- sub( pattern=pattern, replacement="\\1", x=lzp_name )
    loc_start <- sub( pattern=pattern, replacement="\\2", x=lzp_name )
    loc_end   <- sub( pattern=pattern, replacement="\\3", x=lzp_name )
    snp       <- sub( pattern=pattern, replacement="\\4", x=lzp_name )
    
    
    #-----------------------------------------------------------------------------
    #   Get summary statistics for the locus
    #-----------------------------------------------------------------------------
    
    if( cond.or.uncond == "uncond" ){
      
      # Find unconditioned sumstats file
      hits_dir <- file.path( maindir, "hits_locus" )
      ss_name <- paste0( chr, "_", loc_start, "_", loc_end, ".snplist" )
      ss_file <- file.path( hits_dir, ss_name )
      
      # Subset GWAS to SNPs in locus
      snps <- readLines(ss_file)
      ss <- gw[ match( snps, gw$SNP ) , ]
      
    }else if( cond.or.uncond == "cond" ){
      
      # Read in conditioned sumstats file
      # Re-format column names since they differ between single and multi-hit loci
      ss <- fread( ss_cond_files[i] )
      names(ss)[ names(ss) == "Chr" ] <- "chr"
      if( "pC" %in% names(ss) ){
        ss$p <- NULL
        names(ss)[ names(ss) == "pC" ] <- "p"
      }
      
    }
    
    
    #-----------------------------------------------------------------------------
    #   Add LD, coding SNPs, and promoter SNPs to sumstats
    #-----------------------------------------------------------------------------
    
    # Read in LD
    ld_name <- sub( pattern=".jpg$", replacement=".ld", x=lzp_name )
    ld_file <- file.path( maindir, "ld", ld_name )
    ld <- fread(ld_file)
    
    # Add LD to conditioned sumstats
    ss$r2 <- ld$R2[ match( ss$SNP, ld$SNP_B ) ]
    
    # Annotate CS with whether SNPs are coding or promoter
    ss$pch <- 21
    ss$pch[ ss$SNP %in% hhp$SNP ] <- 24
    ss$pch[ ss$SNP %in% hhc$SNP ] <- 25
    
    
    #-----------------------------------------------------------------------------
    #   Plot
    #-----------------------------------------------------------------------------
    
    # If using separated loci, update locus boundaries
    if( !merge.loci ){
      hit_name  <- sub( pattern=".jpg$", replacement="", x=lzp_name )
      loc_start <- cs$lo[ match( hit_name, cs$hit ) ]
      loc_end   <- cs$hi[ match( hit_name, cs$hit ) ]
    }
    
    # Put sumstats into LZP format
    loc <- suppressMessages(
      locus( data  = as.data.frame(ss),
             chrom = "chr", seqname   = ss$chr[1],
             pos   = "bp",  xrange    = c( as.integer(loc_start),
                                           as.integer(loc_end) ),
             p     = "p",   ens_db    = "EnsDb.Hsapiens.v75",
             LD    = "r2",  index_snp = snp ) )
    loc <- suppressMessages( link_recomb( loc, genome="hg19" ) )
    
    # Make LZP
    jpeg( filename=lzp_out, width=480*4, height=480*4, res=75*4 )
    locus_plot( loc, border=TRUE, filter_gene_biotype="protein_coding" )
    dev.off()
  }
}


#-------------------------------------------------------------------------------
#   cojo_wrapper:            Runs the pipeline
#-------------------------------------------------------------------------------

cojo_wrapper <- function(){
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  
  # Make a project area
  dir.create( path=maindir, showWarnings=FALSE, recursive=TRUE )
  
  # Check inputs
  message_header("Check inputs")
  check_arguments()
  
  # Format GWAS
  message_header("Format GWAS")
  format_gwas( maindir, gwas_file, p_threshold, inflation_factor,
               column_names, ld_ref_file )
  
  # Run COJO
  message_header("Run COJO")
  run_cojo(maindir)
  
  # Isolate signals
  message_header("Isolate signals")
  isolate_signals(maindir)
  
  # Compute credible sets
  message_header("Compute credible sets")
  credible_sets(maindir)
  
  # Done
  message_header("Done!")
}


#-------------------------------------------------------------------------------
#   Done
#-------------------------------------------------------------------------------












