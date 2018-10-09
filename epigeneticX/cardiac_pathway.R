# @author Yisong Zhen
# @since  2018-10-08
# @parent
#  read cardiac_mef2.ChIP.R
#---

library(chipenrich)
macs2.rat.peaks <-  .get_macs2_result_filenames(raw_data_path) %>%
                   map(.get_macs_summit_result) %>% 
                   map(function(x) x + 200)
# find the proper set for rat and annotatin database
supported_genomes()
supported_genesets()

# Run enrichment analysis
enrich_pathway <- chipenrich( macs2.rat.peaks[[1]] %>% data.frame(), 
                              genome   = 'rn6', 
                              genesets = 'biocarta_pathway', 
                              out_name = NULL, 
                              locusdef = 'nearest_tss', qc_plots=FALSE)

# Print the results of the analysis
print(enrich_pathway$results)