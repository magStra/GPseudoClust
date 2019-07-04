#code to fetch Shalek data: https://github.com/JohnReid/DeLorean/blob/master/inst/scripts/fetch-Shalek-2014.Rmd
library(DeLorean)
library(dplyr)
options(stringsAsFactors = F)
load('../data/ShalekDeLorean.rda')

#lines 9 to 38 were taken from the DeLorean package documentation
#https://github.com/JohnReid/DeLorean/blob/master/inst/scripts/Shalek-DeLorean.Rmd
dl <- de.lorean(
  shalek.A.expr,
  shalek.A.gene.meta,
  shalek.A.cell.meta)
dl$cell.meta <- mutate(dl$cell.meta,
                       precocious=cell %in% c('LPS_1h_S51', 'LPS_1h_S52'))

time.course.cells <- (
  dl$cell.meta
  %>% filter(! is.na(total),
             "" == assay,
             "LPS" == stimulant | "" == stimulant,
             "" == ko,
             FALSE == disrupted,
             total > 1e6,
             "" == replicate))
dl <- filter_cells(dl, cells=time.course.cells$cell)


dl$cell.meta$cell <- factor(
  dl$cell.meta$cell,
  levels=(shalek.A.cell.meta %>% arrange(capture))$cell)


model.name <- getOption("Shalek.model", 'lowrank')
dl <- estimate.hyper(
  dl,
  sigma.tau=getOption("Shalek.sigma.tau", 1),
  length.scale=getOption("Shalek.length.scale", 5),
  model.name=model.name)

geneSel <- read.table("genesShalek.csv",sep=",",header=T)$genesShalek
dl <- filter_genes(dl, genes=geneSel)
dl <- adjust.by.cell.sizes(dl)
ShalekData <- dl$expr

write.table(ShalekData[geneSel,],file="Shalek.csv",sep=",")
