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


#setwd("/Volumes/HDD2/GC_2017/")
a2a=read_tsv("CS_cds/blastAHL2AHL.tsv",col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                                             "q_start","q_end","db_start","db_end","evalue","bitscore"))
b2b=read_tsv("CS_cds/blastBHL2BHL.tsv",col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                                             "q_start","q_end","db_start","db_end","evalue","bitscore"))
d2d=read_tsv("CS_cds/blastDHL2DHL.tsv",col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                                             "q_start","q_end","db_start","db_end","evalue","bitscore"))
pick_singlecopy=function(.blast,.genome){
cds_length = .blast %>%
  dplyr::filter(qid==dbid) %>%
  dplyr::select(qid,align_length) %>%
  dplyr::rename(id=qid,cds_length=align_length)
print(paste0(.genome," genome whole CDS number = ",length(cds_length$id),
             ", and sum of length = ",sum(cds_length$cds_length)))
perfectly_singlecopy = .blast %>%
  dplyr::count(qid) %>%
  dplyr::filter(n==1) %>%
  dplyr::left_join(cds_length,by=c(qid="id")) %>%
  dplyr::select(-n)
print(paste0(.genome," genome perfectrly singlecopy CDS = ",length(perfectly_singlecopy$qid),
             ", and sum of length = ",sum(perfectly_singlecopy$cds_length)))
write_df(perfectly_singlecopy,paste0("singlecopy_CS_cds/perfectly_singlecopy_",.genome,".tsv"))
mostly_singlecopy = .blast %>%
  dplyr::filter(match_per >90,align_length > 50) %>%
  dplyr::count(qid) %>%
  dplyr::filter(n==1) %>%
  dplyr::left_join(cds_length,by=c(qid="id")) %>%
  dplyr::select(-n)
print(paste0(.genome," genome perfectrly singlecopy CDS = ",length(mostly_singlecopy$qid),
             ", and sum of length = ",sum(mostly_singlecopy$cds_length)))
write_df(mostly_singlecopy,paste0("singlecopy_CS_cds/mostly_singlecopy_",.genome,".tsv"))
}
pick_singlecopy(a2a,"A")
pick_singlecopy(b2b,"B")
pick_singlecopy(d2d,"D")
