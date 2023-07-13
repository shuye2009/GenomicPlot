# CHANGES IN VERSION 0.99.13

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
