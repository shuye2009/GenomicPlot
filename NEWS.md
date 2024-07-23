# CHANGES IN VERSION 1.2.4
   - Modify pseudo count for ratio matrix
   
# CHANGES IN VERSION 1.2.3
   - Add gene id as rownames to heatmap

# CHANGES IN VERSION 1.2.2
   - Fixed a bug in prepare_3parts_genomic_features
   - Add messages to each function 
   
# CHANGES IN VERSION 1.2.1
   - Updated dependency on R >= 4.4.0 and genomation >= 1.36.0

# CHANGES IN VERSION 1.1.10
   - Handle NCBI style of seqlevels (seqname)
   
# CHANGES IN VERSION 1.1.9
   - Fixed a bug in plot_region
   
# CHANGES IN VERSION 1.1.8
## NEW FEATURES
   - Add handling of .gz files for bed and bedGraph format.
   - Add input filter for chromosomes. Only specified chromosomes will be included in visualization and analysis.
   
# CHANGES IN VERSION 1.1.7
## NEW FEATURES
   - Add handle_bedGraph to enable data input of bedGraph format
   
# CHANGES IN VERSION 1.1.6
   - Merged a pull request from Hervé Pagès <notifications@github.com>

# CHANGES IN VERSION 1.1.5
   - Added a function to obtain chromosome info from cached data in the 'circlize' 
   package to produce a Seqinfo object, which is applied to all GRanges and TxDB 
   object
   - Added another function to generate a customized TxDb object from a genome
   annotation (GTF or GFF) file.
   
# CHANGES IN VERSION 1.1.4
   Change chromosome size information source from UCSC web service to cached data in the 'circlize' package, to avoid internet connection issues and web service issues.
   
# CHANGES IN VERSION 1.1.3
   Increase font size for axis labels, align profile with heatmap
   
# CHANGES IN VERSION 1.1.1
   Removed some suggests in DESCRIPTION to reduce installation time and size
## BUG FIXES
   Fixed misalignment between profiles and heatmaps
   
# CHANGES IN VERSION 1.0.0
   Officially released on Bioconductor

# CHANGES IN VERSION 0.99.15

## NEW FEATURES
   None
## SIGNIFICANT CHANGES
   - Provide link to all external and internal data in manual
   - Add examples to all plot functions
## BUG FIXES
   None

# CHANGES IN VERSION 0.99.14
## NEW FEATURES
   None
## SIGNIFICANT CHANGES
   - Make txdb available in /inst/data to avoid repetitive generation of it in examples
   - Store gf5_meta and gf5_genomic in /data so that they can be loaded by calling data()
   - Create /R/data.R for data documentations
## BUG FIXES
   None

# CHANGES IN VERSION 0.99.13
## NEW FEATURES
   None
## SIGNIFICANT CHANGES
   - Add 'title' to draw_combo_plot function arguments.
## BUG FIXES
   - Fixed a bug in plot_peak_annotation, such that txdb$user_genome is treated as a vector of strings rather than a single string
   
# CHANGES IN VERSION 0.99.12
## NEW FEATURES
   - Add PCA plot in plot_bam_correlation.
   - The function plot_5parts_metagene can generate profile plots for both
   5 parts (Promoter, 5'UTR, CDS, 3'UTR, TTS) and 3 parts (Promoter, Gene, TTS) metagene.
## SIGNIFICANT CHANGES
   - The function plot_3parts_metagene is removed.
   - Add setImportParams function to provide default import parameters.
   - Add saveRds option to import parameters to control saving of imported data, the default is FALSE.
   - Add data type and missing file checking for function arguments
## BUG FIXES
   None
