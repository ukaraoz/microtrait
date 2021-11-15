library(dplyr)
library(ggplot2)
library(RColorBrewer)

base = "/Users/ukaraoz/Work/microtrait/code/microtrait"
# 6. read performance metrics for the microtraithmms
hmm_performance = readr::read_delim("./data-raw/microtrait_hmm-performance.txt",
                                    delim = "\t",
                                    col_names = T,
                                    col_types = readr::cols(`microtrait_hmm-dbxref_kegg` = readr::col_character(),
                                                            `microtrait_hmm-dbxref_kegg_score` = readr::col_double(),
                                                            `microtrait_hmm-dbxref_kegg_fscore` = readr::col_double(),
                                                            `microtrait_hmm-dbxref_kegg_tpr` = readr::col_double(),
                                                            `microtrait_hmm-dbxref_kegg_tnr` = readr::col_double(),
                                                            `microtrait_hmm-dbxref_kegg_fpr` = readr::col_double(),
                                                            `microtrait_hmm-dbxref_kegg_fnr` = readr::col_double(),
                                                            `microtrait_hmm-dbxref_kegg_accuracy` = readr::col_double()
                                    )
                            ) %>%
  dplyr::inner_join(hmms_fromrules, by = c("microtrait_hmm-dbxref_kegg" = "microtrait_hmm-dbxref_kegg"))
p = ggplot(hmms_fromrules, aes(x=`microtrait_hmm-dbxref_kegg_fpr`*100, y=`microtrait_hmm-dbxref_kegg_tpr`*100)) +
  geom_vline(xintercept = 0.1, size = 0.2, colour = "red") + geom_hline(yintercept = 75, size = 0.2, colour = "red") +
  geom_point(size = 1, alpha = 0.2) +
  xlab("1-specificity at Fmax") + ylab("sensitivity at Fmax") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
ggsave(p, width = 8, height = 8, file = file.path(base, "data-raw", "hmms_fromrules_FPRvsTPRatFmax.pdf"))

#a = hmms_fromrules %>% pull(`microtrait_hmm-dbxref_kegg_fscore`)
p = ggplot(hmms_fromrules, aes(x=`microtrait_hmm-dbxref_kegg_fscore`)) +
      stat_ecdf() +
      geom_vline(xintercept = 0.8, size = 0.2, colour = "red") +
      geom_hline(yintercept = 0.0835609, size = 0.2, colour = "red") +
      xlab("F") + ylab("Fraction of HMMs with F-score>F") +
      theme(axis.text=element_text(size=16),
            axis.title=element_text(size=20))
ggsave(p, width = 8, height = 8, file = file.path(base, "data-raw", "hmms_fromrules_Fscoredistribution.pdf"))




