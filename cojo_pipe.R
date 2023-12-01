


#-------------------------------------------------------------------------------
#   Simple functions
#-------------------------------------------------------------------------------

# message2:             Just like message, but with the date and a gap pasted in
# message_header:       Nicely-formatted text to break up your log files into 
#                       readable chunks

# message2: Just like message, but with the date and a gap pasted in
message2 <- function(...) message(date(), "     ", ...)

# message_header: 
# Nicely-formatted text to break up your log files into readable chunks
message_header <- function(...){
  message( "\n\n#--------------------------------------------------" )
  message( "#   ", ...)
  message( "#--------------------------------------------------" )
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

format_gwas <- function( maindir, gwas_file, p_threshold, inflation_factor,
                         column_names, ld_ref_file ){
  
  # Load libraries and sources
  library(data.table)
  source("") # source this script
  
  # If output files already exist, skip
  # Read in reference panel SNPs
  # Read in GWAS file
  # If specified, remove SNPs with weak P values (to speed up computation)
  # If specified, correct for inflation
  # Compute the effective sample size for each SNP (N)
  # Format GWAS to have columns: SNP, A1, A2, p, b, se, freq, N
  # Harmonize GWAS with reference panel
  # Write the list of SNP names to file
  # Write the GWAS summary statistics to file
}


#-------------------------------------------------------------------------------
#   run_cojo:         Get independent hits
#-------------------------------------------------------------------------------

run_cojo <- function(maindir){
  
  # Load libraries and sources
  library(data.table)
  source("") # source this script
  
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
  source("") # source this script
  
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
  source("") # source this script
  
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
  source("") # source this script
  
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












