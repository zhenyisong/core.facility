# @author Yisong Zhen
# @since  2018-09-26
#---

pkgs <- c('UpSetR','DiffBind','ChIPQC')

qc_result <- ChIPQC("sample_info.csv", "hg19")
counts    <- dba.count(qc_results, summits=250)

distance <- dist(t(coverage))

# > class(peaks)
# [1] "DBA"

dba.plotHeatmap(peaks, maxSites = peak_count, correlations = FALSE)