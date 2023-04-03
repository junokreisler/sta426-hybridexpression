# RNAseq of synthetic Arabidopsis halleri x lyrata hybrid expression patterns
STA426 project, report comments

Author: Elina Jansone (junokreisler), student MSc CBB

Supervisor: Stefan Milosavljevic (supermaxiste), PhD

The report was compiled using LaTeX, with in-line code pasted separately. Therefore, figures do not appear right after their code.
The Quarto output would have led to a 20-page report. Therefore, it was decided to split the text into two columns, which reduced the size to 10, as well as forced the figures to be smaller for the sake of being able to fit into the pages properly.

The inline code in the report uses the following files:
* `data_halleri_hybrid.rds`
* `data_lyrata_hybrid.rds`
* `rbh_Ahal_Athal_221215` - reciprocal BLAST hits between genes of *A. halleri and A. thaliana*
* `rbh_Alyr_Athal_221215` - reciprocal BLAST hits between genes of *A. lyrata and A. thaliana*
* `GOslim.txt` - the source file for all *A. thaliana* gene names, required by the TopGO package.

In the analysis, reciprocal BLAST hits were used for "translating" *A. halleri* and *A. lyrata* genes to their closest matches in *A. thaliana*. The "translation" had to be done because currently, the base data of the GSEA analysis is more robust and better described for *A. thaliana* than the researched species.

There is also an option to run the analysis with one-way hit files, albeit yielding less confident results:

* `one_way_Ahal_Athal_221215` - one-way BLAST hits between genes of *A. halleri and A. thaliana*, *halleri* in *thaliana*
* `one_way_Alyr_Athal_221215` - one-way BLAST hits between genes of *A. lyrata and A. thaliana*, *lyrata* in *thaliana*

As there are more one-way hits than reciprocal, it is expected that more of the original significant genes will be retained. While around 60% of the significant genes remained after reciprocal-hit-based "translation" to *A. thaliana* equivalents, the percentage of retained genes rises to around 80% when the one-hit files are used instead. However, reliability of the one-way hit results might be misleading, considering that the genes present in the one-way dataset but not in the reciprocal one likely have a low BLAST score. A low BLAST score would indicate that the target sequence might actually have a different kind of behavior than in *A. thaliana*.

The RDS files were generated from count tables using the following code:

```
folder_halleri <- "./countTables/halleri/"
folder_lyrata <- "./countTables/lyrata/"
folder_hybrid <- "./countTables/hybrid/"

files_halleri <- list.files(folder_halleri, pattern = ".txt$")
files_lyrata <- list.files(folder_lyrata, pattern = ".txt$")
files_hybrid_halleri <- list.files(folder_hybrid, pattern = "[[:graph:]]halleri_[1245].txt$")
files_hybrid_lyrata <- list.files(folder_hybrid, pattern = "[[:graph:]]lyrata_[1245].txt$")

paths_halleri <- paste0(folder_halleri, files_halleri)
paths_lyrata <- paste0(folder_lyrata, files_lyrata)
paths_hybrid_halleri <- paste0(folder_hybrid, files_hybrid_halleri)
paths_hybrid_lyrata <- paste0(folder_hybrid, files_hybrid_lyrata)

# return DGEList geneID with all sample counts given the file path list
concat_samples <- function(paths_list) {
  dataset <- data.table::fread(paths_list[1])
  for (path in paths_list[2:length(paths_list)]) {
  dataset <- cbind(dataset,
                   data.table::fread(path)[,7])
  }
  return (dataset[,-c(2:6)])
}

# make DGELists
data_halleri <- concat_samples(paths_halleri)
data_lyrata <- concat_samples(paths_lyrata)
data_hybrid_halleri <- concat_samples(paths_hybrid_halleri)
data_hybrid_lyrata <- concat_samples(paths_hybrid_lyrata)

data_halleri_hybrid <- cbind(data_halleri,
                             data_hybrid_halleri[,2:5])
colnames(data_halleri_hybrid) <- c("geneID", "Hal_1","Hal_2","Hal_3", "Hal_4",
                                   "Hyb_Hal_1", "Hyb_Hal_2", "Hyb_Hal_3", "Hyb_Hal_4")
DGEList_hal_hyb <- DGEList(counts = as.matrix(data_halleri_hybrid[,2:9]),
                           genes = data_halleri_hybrid[,1],
                           samples = c("Hal_1","Hal_2","Hal_3", "Hal_4",
                                       "Hyb_Hal_1", "Hyb_Hal_2", "Hyb_Hal_3", "Hyb_Hal_4"),
                           group = rep(1:2, each = 4))

data_lyrata_hybrid <- cbind(data_lyrata,
                            data_hybrid_lyrata[,2:5])
colnames(data_lyrata_hybrid) <- c("geneID", "Lyr_1","Lyr_2","Lyr_3", "Lyr_4",
                                   "Hyb_Lyr_1", "Hyb_Lyr_2", "Hyb_Lyr_3", "Hyb_Lyr_4")
DGEList_lyr_hyb <- DGEList(counts = as.matrix(data_lyrata_hybrid[,2:9]),
                           genes = data_lyrata_hybrid[,1],
                           samples = c("Lyr_1","Lyr_2","Lyr_3", "Lyr_4",
                                       "Hyb_Lyr_1", "Hyb_Lyr_2", "Hyb_Lyr_3", "Hyb_Lyr_4"),
                           group = rep(1:2, each = 4))
write_rds(DGEList_hal_hyb, 'data_halleri_hybrid.rds')
write_rds(DGEList_lyr_hyb, 'data_lyrata_hybrid.rds')
```

Volcano plots with all original significant points as native gene names (mind that gene "translation" caused a notable loss in gene number):

```
all_tags_hh <- topTags(
        test_halleri_hybrid, n = 17167)
all_tags_lh <- topTags(
        test_lyrata_hybrid, n = 16965)

ggplot(data = all_tags_hh$table, 
        aes(y = -log10(FDR), x = logFC)) +
        geom_point(size = 0.5) + 
        geom_text(aes(label = geneID), 
        check_overlap = TRUE, size = 1.5, 
        nudge_y = 3) +
        gghighlight(-log10(FDR) > 2 & 
        ((abs(logFC) >= 5 | -log10(FDR) > 50) | 
        logFC < -3.7), 
        unhighlighted_params = aes(colour = 'grey'))

ggplot(data = all_tags_lh$table, 
        aes(y = -log10(FDR), x = logFC)) +
        geom_point(size = 0.5) + 
        geom_text(aes(label = geneID),
        check_overlap = TRUE, size = 1.5, 
        nudge_y = 3) +
        gghighlight(-log10(FDR) > 50 | 
        ((abs(logFC) >= 5 | -log10(FDR) > 50) | 
        logFC < -3.7), 
        unhighlighted_params = aes(colour = 'grey'))
```
