# @author Yisong Zhen
# @since  2018-09-28
#---

pkgs      <- c('UpSetR','DiffBind','ChIPQC')
samples   <- data.frame( SampleID   = ,
                         Tissue     = ,
                         Factor     = ,
                         Replicate  =,
                         bamReads   = ,
                         Peaks      =  )
qc_result <- ChIPQC('sample_info.csv', annotaiton = 'hg19')
counts    <- dba.count(qc_results, summits = 250)
contrast <- DBA____

# Establish the contrast to compare the two tumor types
dba_peaks <- dba.___(ar_binding, ___=contrast, minMembers=2)

ar_diff <- dba.analyze(ar_binding)

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