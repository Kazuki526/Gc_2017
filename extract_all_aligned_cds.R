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
#setwd("/Volumes/HDD2/GC_2017)
a_length = read_tsv("singlecopy_CS_cds/mostly_singlecopy_A.tsv") %>%
  dplyr::rename(id=qid)
b_length = read_tsv("singlecopy_CS_cds/mostly_singlecopy_B.tsv") %>%
  dplyr::rename(id=qid)
d_length = read_tsv("singlecopy_CS_cds/mostly_singlecopy_D.tsv") %>%
  dplyr::rename(id=qid)

single_filter = commandArgs(trailingOnly = T)[1]
#single_filter = "perfectly"

blast_a2b = read_tsv(paste0("singlecopy_CS_cds/blast_A2B_",single_filter,".tsv"),
                     col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                   "q_start","q_end","db_start","db_end","evalue","bitscore"))
blast_a2d = read_tsv(paste0("singlecopy_CS_cds/blast_A2D_",single_filter,".tsv"),
                     col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                   "q_start","q_end","db_start","db_end","evalue","bitscore"))
b2a_double = read_tsv(paste0("singlecopy_CS_cds/blast_B2A_",single_filter,".tsv"),
                     col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                   "q_start","q_end","db_start","db_end","evalue","bitscore")) %>%
  dplyr::count(qid) %>% dplyr::filter(n>1) %>%
  dplyr::mutate(focal="double") %>% dplyr::select(-n)
b2d_double = read_tsv(paste0("singlecopy_CS_cds/blast_B2D_",single_filter,".tsv"),
                     col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                   "q_start","q_end","db_start","db_end","evalue","bitscore")) %>%
  dplyr::count(qid) %>% dplyr::filter(n>1) %>%
  dplyr::mutate(focal="double") %>% dplyr::select(-n)
d2b_double = read_tsv(paste0("singlecopy_CS_cds/blast_D2B_",single_filter,".tsv"),
                     col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                   "q_start","q_end","db_start","db_end","evalue","bitscore")) %>%
  dplyr::count(qid) %>% dplyr::filter(n>1) %>%
  dplyr::mutate(focal="double") %>% dplyr::select(-n)
d2a_double = read_tsv(paste0("singlecopy_CS_cds/blast_D2A_",single_filter,".tsv"),
                     col_names = c("qid","dbid","match_per","align_length","mismatch","gap",
                                   "q_start","q_end","db_start","db_end","evalue","bitscore")) %>%
  dplyr::count(qid) %>% dplyr::filter(n>1) %>%
  dplyr::mutate(focal="double") %>% dplyr::select(-n)
b_double = rbind(b2a_double,b2d_double) %>%dplyr::distinct()
d_double = rbind(d2a_double,d2b_double) %>%dplyr::distinct()

#single copyだが2つのCDSと相同性があるものがあるので除去
blast1 = blast_a2b %>%
  dplyr::left_join(b_double %>%rename(dbid=qid)) %>%
  dplyr::filter(is.na(focal)) %>% dplyr::select(-focal) %>%
  dplyr::group_by(qid) %>%
  dplyr::mutate(n=n()) %>% dplyr::ungroup() %>%
  dplyr::filter(n==1) %>%
  dplyr::left_join(a_length %>% dplyr::rename(qid=id,a_length=cds_length)) %>%
  dplyr::left_join(b_length %>% dplyr::rename(dbid=id,b_length=cds_length))

blast2 =blast_a2d %>%
  dplyr::left_join(d_double %>%rename(dbid=qid)) %>%
  dplyr::filter(is.na(focal)) %>% dplyr::select(-focal) %>%
  dplyr::group_by(qid) %>%
  dplyr::mutate(n=n()) %>% dplyr::ungroup() %>%
  dplyr::filter(n==1) %>%
  dplyr::left_join(a_length %>% dplyr::rename(qid=id,a_length=cds_length)) %>%
  dplyr::left_join(d_length %>% dplyr::rename(dbid=id,d_length=cds_length))

matched_list = blast1 %>% dplyr::rename(aid=qid,bid=dbid,b_match=match_per,b_align_length=align_length,
                         a_start_2b=q_start,a_end_2b=q_end,b_start=db_start,b_end=db_end) %>%
  dplyr::select(aid,bid,b_match,b_align_length,a_length,b_length,
                a_start_2b,a_end_2b,b_start,b_end) %>%
  dplyr::full_join(blast2 %>% dplyr::rename(aid=qid,did=dbid,d_match=match_per,d_align_length=align_length,
                         a_start_2d=q_start,a_end_2d=q_end,d_start=db_start,d_end=db_end) %>%
              dplyr::select(aid,did,d_match,d_align_length,a_length,d_length,
                            a_start_2d,a_end_2d,d_start,d_end)) %>% 
  dplyr::filter(!is.na(bid),!is.na(did))

matched_cds = matched_list %>%
  dplyr::filter(a_start_2b < 20,a_start_2d < 20,b_start < 20,d_start < 20,
                a_end_2b > a_length-20,a_end_2d > a_length-20,b_end > b_length-20,d_end > d_length-20)

matched_cds_95 = matched_list %>%
  dplyr::filter(a_length * 0.95 < b_align_length, a_length * 0.95 < d_align_length,
                b_length * 0.95 < b_align_length, d_length * 0.95 < d_align_length)
matched_region = matched_list %>%
  dplyr::mutate(a_start=ifelse(a_start_2b > a_start_2d,a_start_2b,a_start_2d),
                a_end=ifelse(a_end_2b < a_end_2d,a_end_2b,a_end_2d)) %>%
  dplyr::mutate(all_alingn_length=a_end - a_start) 

print(paste0(single_filter," singlecopy cds align result"))
print(paste0("High&Low conf: matches to each other CDS = ",length(matched_region$a_length)," genes"))
print(paste0("               and sum of matched region = ",sum(matched_region$a_length)," bp"))
print(paste0("High&Low conf: CDS length 95% matches to each other CDS = ",length(matched_cds_95$a_length)," genes"))
print(paste0("               and its length = ",sum(matched_cds_95$a_length), " bp"))
print(paste0("High&Low conf: CDS length 95% matches to each other CDS = ",length(matched_cds$a_length)," genes"))
print(paste0("               and its length = ",sum(matched_cds$a_length), " bp"))

matched_cdsh = matched_cds %>%
  dplyr::filter(!str_detect(aid,"LC"),!str_detect(bid,"LC"),!str_detect(did,"LC"))
matched_cdsh_95 = matched_cds_95 %>%
  dplyr::filter(!str_detect(aid,"LC"),!str_detect(bid,"LC"),!str_detect(did,"LC"))
matched_regionh = matched_region %>%
  dplyr::filter(!str_detect(aid,"LC"),!str_detect(bid,"LC"),!str_detect(did,"LC"))

print(paste0(single_filter," singlecopy cds align result"))
print(paste0("only High conf: matches to each other CDS = ",length(matched_regionh$a_length)," genes"))
print(paste0("                and sum of matched region = ",sum(matched_regionh$a_length)," bp"))
print(paste0("only High conf: CDS length 95% matches to each other CDS = ",length(matched_cdsh_95$a_length)," genes"))
print(paste0("                and its length = ",sum(matched_cdsh_95$a_length), " bp"))
print(paste0("only High conf: all site CDS matches to each other CDS = ",length(matched_cdsh$a_length)," genes"))
print(paste0("                and its length = ",sum(matched_cdsh$a_length), " bp"))

print(paste0("make singlecopy_CS_cds/",single_filter,"_singlecopy_all_aligned.tsv (only high conf & all matched)"))
write_df(matched_cdsh,paste0("singlecopy_CS_cds/",single_filter,"_singlecopy_all_aligned.tsv"))

