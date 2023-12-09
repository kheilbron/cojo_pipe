


#-------------------------------------------------------------------------------
#   Simple functions
#-------------------------------------------------------------------------------

# message2:             Just like message, but with the date and a gap pasted in
# message_header:       Nicely-formatted text to break up your log files into 
#                       readable chunks
# n_eff:                Get the effective sample size of a case-control GWAS

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


#-------------------------------------------------------------------------------
#   check_inputs:     Check validity of inputs
#-------------------------------------------------------------------------------

check_inputs <- function(){
  #
}


#-------------------------------------------------------------------------------
#   format_gwas:      Write GWAS sumstats in COJO format
#-------------------------------------------------------------------------------

format_gwas <- function( maindir          = "~/cojo",
                         ld_panel         = "hrc",
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
  # Format GWAS to have columns: SNP, A1, A2, p, b, se, freq, N
  # Write the list of SNP names to file
  # Write the GWAS summary statistics to file
  
  
  #-------------------------------------------------------------------------------
  #   Input descriptions
  #-------------------------------------------------------------------------------
  
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
  
  
  #-------------------------------------------------------------------------------
  #   Read in GWAS and HRC, format columns, subset to shared SNPs
  #-------------------------------------------------------------------------------
  
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
    gw <- gw0[ gw0$P < p_threshold , ]
  }else{
    gw <- gw0
  }
  
  # Make columns for chromosome and position
  if( "chr_bp" %in% names(gw) ){
    message2("Make columns for chromosome and position")
    gw$chr <- as.integer( sub( pattern     = "^chr([[:alnum:]]+)[[:punct:]]([[:digit:]]+)$", 
                               replacement = "\\1", 
                               x=gw$chr_bp ) )
    gw$bp  <- as.integer( sub( pattern     = "^chr([[:alnum:]]+)[[:punct:]]([[:digit:]]+)$", 
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
    message2("A study-wide effective sample size has been provided and will be applied to all SNPs")
    gw$N <- n
  }else if( "n1" %in% names(gw) & "n0" %in% names(gw) ){
    message2("Columns for number of cases and controls have been provided for each SNP, computing the effective sample size")
    gw$N <- n_eff( gw$n1, gw$n0 )
  }
  
  # If ORs have been provided, convert to logORs
  if( "or" %in% names(gw) ){
    message2("Convert ORs to logORs")
    gw$b <- log(gw$or)
  }
  
  
  #-------------------------------------------------------------------------------
  #   Harmonize GWAS and reference panel alleles: without allele frequency information
  #-------------------------------------------------------------------------------
  
  # Subset GWAS and reference panel to shared SNPs based on chromosome and position
  message2("Subset GWAS and reference panel to shared SNPs based on chromosome and position")
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
            " (", round( 100*NROW(gw2)/NROW(gw), 2 ), "%) were found in the reference panel" )
  message2( "Of the ", NROW(hrc), " reference panel SNPs, ", NROW(hrc2), 
            " (", round( 100*NROW(hrc2)/NROW(hrc), 2 ), "%) were found in the GWAS" )
  
  # Find palindromic SNPs and SNPs with alleles that are flipped in 
  # the reference panel v. GWAS ('discordant')
  message2("Find palindromic SNPs and SNPs with alleles that are flipped in ",
           "the reference panel v. GWAS ('discordant')")
  pal     <- ( hrc2$alt=="A" & hrc2$ref=="T" ) | 
    ( hrc2$alt=="T" & hrc2$ref=="A" ) | 
    ( hrc2$alt=="C" & hrc2$ref=="G" ) | 
    ( hrc2$alt=="G" & hrc2$ref=="C" )
  discord <- hrc2$alt != gw2$A1
  message2( sum(pal),     "/", NROW(hrc2), " (", round( 100 * sum(pal)     / NROW(hrc2), 2 ), "%) SNPs are palindromic" )
  message2( sum(discord), "/", NROW(hrc2), " (", round( 100 * sum(discord) / NROW(hrc2), 2 ), "%) SNPs have discordant alleles" )
  
  # For discordant non-palindromic SNPs: flip GWAS alleles
  message2("For discordant non-palindromic SNPs: flip GWAS alleles")
  disc_nonpal <- discord & !pal
  original_gwas_ref <- gw2$A2
  original_gwas_alt <- gw2$A1
  gw2$A2[disc_nonpal] <- original_gwas_alt[disc_nonpal]
  gw2$A1[disc_nonpal] <- original_gwas_ref[disc_nonpal]
  if( "freq" %in% names(gw2) ){
    gw2$freq[disc_nonpal] <- 1  - gw2$freq[disc_nonpal]
  }
  message2( sum(disc_nonpal), "/", sum(discord), 
            " (", round( 100 * sum(disc_nonpal) / sum(discord), 2 ), 
            "%) discordant SNPs were non-palindromic, flipping alleles" )
  
  
  #-------------------------------------------------------------------------------
  #   Harmonize GWAS and HRC alleles: with allele frequency information
  #-------------------------------------------------------------------------------
  
  # If allele frequency data is available for the GWAS
  if( "freq" %in% names(gw2) ){
    
    # For discordant palindromic SNPs with compatible AFs: flip GWAS alleles
    diff_af_disc    <- abs( ( 1 - gw2$freq ) - hrc2$af ) > 0.2 | hrc2$af > 0.4
    disc_pal_compat <- discord & pal & !diff_af_disc
    gw2$A2[disc_pal_compat]  <- original_gwas_alt[disc_pal_compat]
    gw2$A1[disc_pal_compat]  <- original_gwas_ref[disc_pal_compat]
    gw2$freq[disc_pal_compat] <- 1  - gw2$freq[disc_pal_compat]
    message2( sum(disc_pal_compat), "/", sum(discord), 
              " (", round( 100 * sum(disc_pal_compat) / sum(discord), 2 ), 
              "%) discordant SNPs were palindromic but with AFs that were clearly ",
              "compatible with the reference panel, flipping alleles" )
    
    # Flag discordant palindromic SNPs with incompatible AFs for removal
    disc_pal_incompat <- discord & pal & diff_af_disc
    message2( sum(disc_pal_incompat), "/", sum(discord), 
              " (", round( 100 * sum(disc_pal_incompat) / sum(discord), 2 ), 
              "%) discordant SNPs were palindromic and had AFs that were not ",
              "clearly compatible with the reference panel, flagging for removal" )
    
    # For concordant palindromic SNPs with clearly incompatible AFs: flip GWAS alleles
    diff_af_conc <- abs( gw2$freq - hrc2$af ) > 0.2
    common_af    <- hrc2$af > 0.4
    conc_pal_incompat <- !discord & pal & diff_af_conc & !common_af
    gw2$freq[conc_pal_incompat] <- 1  - gw2$freq[conc_pal_incompat]
    message2( sum(conc_pal_incompat), "/", sum( !discord & pal ), 
              " (", round( 100 * sum(conc_pal_incompat) / sum( !discord & pal ), 2 ), 
              "%) concordant palindromic SNPs had AFs that were clearly ",
              "different from the reference panel, flipped alleles" )
    
    # Report the number of concordant palindromic SNPs with clearly compatible AFs
    conc_pal_compat <- !discord & pal & !diff_af_conc & !common_af
    message2( sum(conc_pal_compat), "/", sum( !discord & pal ), 
              " (", round( 100 * sum(conc_pal_compat) / sum( !discord & pal ), 2 ), 
              "%) concordant palindromic SNPs had AFs that were clearly similar ",
              "to the reference panel, no action" )
    
    # Flag concordant palindromic SNPs with HRC MAF > 40% for removal
    conc_pal_ambig <- !discord & pal & common_af
    message2( sum(conc_pal_ambig), "/", sum( !discord & pal ), 
              " (", round( 100 * sum(conc_pal_ambig) / sum( !discord & pal ), 2 ), 
              "%) concordant palindromic SNPs had reference panel MAF > 40%, ",
              "flagging for removal" )
    
    # Remove flagged SNPs
    message2("Remove flagged SNPs")
    bad_snps <- disc_pal_incompat | conc_pal_ambig
    gw3  <- gw2[  !bad_snps , ]
    hrc3 <- hrc2[ !bad_snps , ]
    
  }else{
    
    # Remove palindromic SNPs
    message2("Remove palindromic SNPs")
    gw3  <- gw2[  !pal , ]
    hrc3 <- hrc2[ !pal , ]
  }
  
  
  #-------------------------------------------------------------------------------
  #   Wrap up harmonization
  #-------------------------------------------------------------------------------
  
  # Report the change in number of SNPs
  message2( "After harmonizing GWAS and reference panel SNPs, ", 
            NROW(gw3), "/", NROW(gw2), " (", round( 100 * NROW(gw3) / NROW(gw2), 2 ), 
            "%) remain" )
  
  # Check that CPRA is 100% identical now
  cpra_gw  <- paste( gw3$chr,  gw3$bp,  gw3$A2,   gw3$A1,   sep="_" )
  cpra_hrc <- paste( hrc3$chr, hrc3$bp, hrc3$ref, hrc3$alt, sep="_" )
  if( all( cpra_gw == cpra_hrc ) ){
    message2("CPRA is now identical for all GWAS and reference panel SNPs")
  }else{
    diff_cpra <- head( which( cpra_gw != cpra_hrc ) )
    stop( paste( "Error: not all CPRA are identical for GWAS and reference panel",
                 "SNPs. Here are (up to) the first 6:", diff_cpra, collapse=" " ) )
  }
  
  # Replace GWAS SNP names with reference panel names
  message2("Replace GWAS SNP names with reference panel names")
  gw3$SNP <- hrc3$snp
  
  
  #-------------------------------------------------------------------------------
  #   Format and write outputs
  #-------------------------------------------------------------------------------
  
  # Dump GWAS summary statistics
  # Column names/order: SNP, A1, A2, p, b, se, freq, N
  message2("Dump GWAS summary statistics")
  gw_outfile <- file.path( maindir, "gwas_sumstats.tsv" )
  gw_out <- gw3[ , c( "SNP", "A1", "A2", "freq", "b", "se", "p", "N" ) ]
  fwrite( x=gw_out, file=gw_outfile, sep="\t" )
  
  
  #-------------------------------------------------------------------------------
  #   Done
  #-------------------------------------------------------------------------------
}


#-------------------------------------------------------------------------------
#   run_cojo:         Get independent hits
#-------------------------------------------------------------------------------

run_cojo <- function(maindir){
  
  # Load libraries and sources
  library(data.table)
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  
  # Create a directory for per-chromosome COJO outputs
  
  # Set static COJO arguments
  gcta_binary    <- ""                        # path to GCTA software
  cojo_file      <- file.path( maindir, "" )  # from previous step
  cojo_p         <- 5e-8                      # p-value cutoff for significance
  max_hits       <- 100                       # prevents runaway computation
  
  # Loop through chromosomes
  for( i in chromosomes ){
    
    # Set per-chromosome COJO arguments
    bed_file       <- ""
    outfile_prefix <- ""
    
    # If the output file already exists, skip
    
    # Run COJO
    message2( "Running COJO on chromosome: ", i )
    cmd <- paste( gcta_binary, 
                  "--bfile",         bed_file, 
                  "--cojo-file",     cojo_file,
                  "--cojo-slct",
                  "--cojo-p",        cojo_p,
                  "--cojo-top-SNPs", max_hits,
                  "--out",           outfile_prefix )
    system(cmd)
  }
  
  # Collate results across chromosomes into a single file
}


#-------------------------------------------------------------------------------
#   isolate_signals:  Extract conditioned sumstats for each independent hit
#-------------------------------------------------------------------------------

isolate_signals <- function(maindir){
  
  # Load libraries and sources
  library(data.table)
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  
  # If output files already exist, skip
  # Read in collated COJO results
  # Group hits into loci
  # Loop through loci
  #   If there is only one hit in the locus, subset the formatted GWAS sumstats
  #   to ±1Mb around the hit
  #   If there are multiple hits in the locus, sub-loop through hits
  #      Write a file containing all hits in the locus except the focal hit
  #      Re-run COJO with the --cojo-cond argument pointing to this file
}


#-------------------------------------------------------------------------------
#   credible_sets:    Compute credible sets for each (conditioned) hit
#-------------------------------------------------------------------------------

credible_sets <- function(maindir){
  
  # Load libraries and sources
  library(data.table)
  library(coloc)
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  
  # If output files already exist, skip
  # Loop through independent hits
  #   Read in conditioned GWAS summary statistics
  #   Use the coloc R package to derive credible sets
  #   Write sumstats for credible set SNPs to file
  dataset <- list( SNP     = data$snp, 
                   beta    = data$b, 
                   varbeta = data$se, 
                   N       = N,
                   s       = s,
                   type    = "cc" )
  out <- finemap.abf( dataset=dataset )
}


#-------------------------------------------------------------------------------
#   cojo_wrapper:     Runs the pipeline
#-------------------------------------------------------------------------------

cojo_wrapper <- function(){
  
  # Load libraries and sources
  source("/projects/0/prjs0817/repos/cojo_pipe/z_cojo_pipe.R")
  
  # Make a project area
  dir.create( path=maindir, showWarnings=FALSE, recursive=TRUE )
  
  # Check inputs
  message_header("Check inputs")
  check_inputs()
  
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












