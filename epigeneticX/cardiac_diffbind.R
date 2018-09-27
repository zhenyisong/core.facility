# @author Yisong Zhen
# @since  2018-09-26
#---

pkgs <- c('UpSetR','DiffBind','ChIPQC')

qc_result <- ChIPQC("sample_info.csv", "hg19")
counts    <- dba.count(qc_results, summits=250)


Load the data:

reads <- readGAlignments(bam)
reads_gr <- granges(reads[[1]])
Obtain average fragment length:

frag_length <- fragmentlength(qc_report)["GSM1598218"]
Extend reads and compute coverage:

reads_ext <- resize(reads_gr, width=frag_length)
cover_ext <- coverage(reads_ext)
distance <- dist(t(coverage))

# > class(peaks)
# [1] "DBA"

dba.plotHeatmap(peaks, maxSites = peak_count, correlations = FALSE)