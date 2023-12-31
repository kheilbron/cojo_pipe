
#   a_find_EUR_plus_EAS_HRC_SNPs_w_MAC_ge_10.R


# Before firing up R
# module load PLINK/2.00a3.6-GCC-11.3.0

# Make scratch directories to work in
scratch_dir <- file.path( "/scratch-shared", Sys.getenv("USER"), "HRC_EUR_EAS" )
plink_dir   <- file.path( scratch_dir, "plink_freqs" )
rare_dir    <- file.path( scratch_dir, "rare_snp_files" )
common_dir  <- file.path( scratch_dir, "common_snp_files" )
dir.create( path=plink_dir,  showWarnings=FALSE, recursive=TRUE )
dir.create( path=rare_dir,   showWarnings=FALSE, recursive=TRUE )
dir.create( path=common_dir, showWarnings=FALSE, recursive=TRUE )

# Function: make a .bim file including AF and AC
bim_w_af_and_ac <- function(chromosome){
  
  # If output exists, skip
  rare_file   <- file.path( rare_dir,   paste0( chromosome, ".tsv" ) )
  common_file <- file.path( common_dir, paste0( chromosome, ".tsv" ) )
  if( all( file.exists( rare_file, common_file ) ) ){
    message( "Output already exists for chromosome ", chromosome, ", skipping" )
    return()
  }
  
  # Run plink2 to get EAF of HRC SNPs
  bfile     <- paste0( "/gpfs/work5/0/pgcdac/imputation_references/",
                       "HRC.r1-1_merged_EUR_EAS_panel/HRC.r1-1.EGA.GRCh37.chr", 
                       chromosome, ".impute.plink.combined.EUR.2191.EAS.538" )
  plink_out <- file.path( plink_dir, chromosome )
  plink_cmd <- paste( "plink2 --bfile", bfile, "--freq --out", plink_out )
  system(plink_cmd)
  
  # Read in AF file, calculate allele count
  library(data.table)
  af <- fread( paste0( plink_out, ".afreq" ) )
  names(af) <- c( "chr", "snp", "ref", "alt", "af", "n" )
  af$ac <- round( af$af * af$n )
  
  # Read in .bim file
  bim <- fread( paste0( bfile, ".bim" ) )
  names(bim) <- c( "chr", "snp", "cm", "bp", "alt", "ref" )
  
  # Merge
  gcols_bim <- c( "snp", "chr", "bp", "ref", "alt" )
  gcols_af  <- c( "af", "ac" )
  m <- cbind( bim[ , ..gcols_bim ],
              af[  , ..gcols_af  ] )
  
  # Remove MAC < 10
  # Remove MAF < 1%
  rare_df   <- m[ m$ac >= 10 , ]
  common_df <- m[ m$af >= 0.01 , ]
  
  # Write
  fwrite( x=rare_df,   file=rare_file,   sep="\t" )
  fwrite( x=common_df, file=common_file, sep="\t" )
}

# Loop through chromosomes making improved .bim files
for( chromosome in 1:22 ){
  message("\n\n\n")
  message( "Starting chromosome: ", chromosome )
  bim_w_af_and_ac(chromosome)
}

# Prepare to concatenate files
rare_snp_files   <- list.files( path=rare_dir,   full.names=FALSE )
common_snp_files <- list.files( path=common_dir, full.names=FALSE )
cat_rare_files   <- paste( rare_snp_files,   collapse=" " )
cat_common_files <- paste( common_snp_files, collapse=" " )
rare_outfile     <- file.path( scratch_dir, "hrc_eur_eas_snps_mac_ge_10.tsv" )
common_outfile   <- file.path( scratch_dir, "hrc_eur_eas_snps_maf_ge_0.01.tsv" )

# Concatenate rare files
setwd(rare_dir)
system( paste0( "head -n1 ", rare_snp_files[1], " > ", rare_outfile ) )
system( paste0( 'cat ', cat_rare_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', rare_outfile ) )

# Concatenate common files
setwd(common_dir)
system( paste0( "head -n1 ", common_snp_files[1], " > ", common_outfile ) )
system( paste0( 'cat ', cat_common_files, ' | grep -v "snp\tchr\tbp\tref\talt\taf\tac" >> ', common_outfile ) )

# Move back to my home directory
rare_home_dir_file   <- "/home/heilbron/projects/pops/data/hrc_eur_eas_snps_mac_ge_10.tsv"
common_home_dir_file <- "/home/heilbron/projects/pops/data/hrc_eur_eas_snps_maf_ge_0.01.tsv"
file.copy( from=rare_outfile,   to=rare_home_dir_file )
file.copy( from=common_outfile, to=common_home_dir_file )






