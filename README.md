# Project Proposal

**Aim:** Find RNAs incorrectly labelled as non-coding RNAs (ncRNAs), when they are actually translated. 

**Importance:** Recent research has shown that more regions of the genome are translated than was originally believed. Many RNAs thought to be non-coding have been found to contain small open-reading frames (ORFs) which encode for micropeptides. These mimcropeptides have been found to play important biological roles [1]. Thus, identifying incorrectly labelled ncRNAs is an important step towards understanding the biological roles these micropeptides play.

**Proposed Methodology:** Using several machine- and deep-learning trained classifiers  to distinguish between mRNAs and ncRNAs, find the ncRNAs misclassified by all classifiers. These are the candidate incorrectly classified ncRNAs. Use Ribo-seq data (which isolates only RNA being translated by ribosomes) to match candidate RNAs obtained with RNAs being translated. In case of a positive match, a mis-labelled ncRNA will be found. 

---

# Weekly Progress Log

**Dataset:** In this study [EBI](https://www.ebi.ac.uk/ena/data/view/PRJNA278389) [NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66929) the authors performed both RNA-seq, GRO-seq and Ribo-seq. Control samples from each have been downloaded and are being analyzed for quality. 

## Week 3

**Quality Control** 

GRO-seq samples show better quality when compared to Ribo-seq samples.

## Week 4

**Adapter Trimming**

Done

## Week 5

**Alignment**

| Seq type  | File       |  Alignment tool  | Overall alignment (Unique alignment) |
| --------  | ----       | ---------------- | ----------------- | 
| Gro-seq   | SRR1916554 | Bowtie2          |97.92% (34.34%)    |
| Gro-seq   | SRR1916552 | Bowtie2          | 97.45%  (32.01%)  |
| Ribo-seq  | SRR1916542 | Tophat           |  87.9% (~39%)     |
| Ribo-seq  | SRR1916544 | Tophat           |  71.7% (~38%)     |
| Ribo-seq  | SRR1916545 | Tophat           |  72.5% (~37%)     |
| Ribo-seq  | SRR1916542 | Bowtie2          | 66.30% (21.20%)   |
| Ribo-seq  | SRR1916544 | Bowtie2          | 69.30% (20.98%)   |
| Ribo-seq  | SRR1916545 | Bowtie2          | 70.21% (21.13%)   |

---

## Week 8

Added data from this paper [CPPred: coding potential prediction based on the global description of RNA sequence](https://academic.oup.com/nar/article/47/8/e43/5314020) to classify RNAs as coding or non-coding. Replicate the results from this paper. Use SVM to classify, F1 scores are given below:

| Data      |    Coding RNA    | Non-Coding RNA   |  
| --------  | ----------------------   |------------------------  |
|  Test     |    0.95                  |  0.94                    |     
|  Train    |   0.96                   |  0.94                    |   


Length distribution differs significantly though. Not sure why this is so.

![image](non_coding_coding_RNA_sequences/length_distribution.png)

---

## Week 9

Started working on CNN classifier with sequence embeddings for ncRNA and cRNA.

Model training is going to take time, especially since I still need to tune parameters. 

This figure explains the workflow for the complete model and why I'm using SVM.

![](non_coding_coding_RNA_sequences/project_workflow.jpg)

---
## Week 10

Downloaded gtf file from Ensemble. See this [biostars answer](https://www.biostars.org/p/195317/).

Counted ncRNAs using htseq. 

Finished working on xgboost and generating figures. (See presentation for more details)

---
**References:**

1. Hartford, C.C.R. and Lal, A., 2020. When long non-coding becomes protein-coding. Molecular and Cellular Biology.

---

