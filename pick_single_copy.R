library(tidyverse)

write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}


setwd("/Volumes/HDD2/")
a2a=read_tsv("GC_2017/CS_cds/blastAHL2AHL.tsv",col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                                             "q_start","q_end","db_start","db_end","evalue","bitscore"))
b2b=read_tsv("GC_2017/CS_cds/blastBHL2BHL.tsv",col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                                             "q_start","q_end","db_start","db_end","evalue","bitscore"))
d2d=read_tsv("GC_2017/CS_cds/blastDHL2DHL.tsv",col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                                             "q_start","q_end","db_start","db_end","evalue","bitscore"))

a_length = a2a %>%
  dplyr::filter(qid==dbid) %>%
  dplyr::select(qid,align_length) %>%
  dplyr::rename(id=qid,cds_length=align_length)

a_perfectly_singlecopy = a2a %>%
  dplyr::count(qid) %>%
  dplyr::filter(n==1)
