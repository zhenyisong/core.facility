
# @refrences
#   Biostatistical Design and Analysis Using R: a practical guide
# raw data was downloaded from here
# http://users.monash.edu.au/~murray/BDAR/
#---

library(ggfortify)
library(multcomp)
setwd('E:\\Downloads\\dataAndScripts\\Chapter10\\Data')

medley      <- read.table('medley.csv', header=T, sep=',')

medley$ZINC <- factor( medley$ZINC, 
                       levels  = c('HIGH', 'MED', 'LOW', 'BACK'), 
                       ordered = F)
medley.aov  <- aov(DIVERSITY~ZINC, medley)

autoplot(medley.aov)

anova(medley.aov)

summary( glht(medley.aov, linfct = mcp(ZINC = 'Tukey')))

