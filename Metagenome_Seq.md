
# **MetagenomeSeq**

# **1.Introduction**

Metagenomics is the genetic analysis that aims to determine a microbial
population in a given environment (Sleator, Shortall and Hill, 2008).
The goals of Metagenomics generally include: identifying functional
genes and/or new metabolic pathways, estimating microbial diversity,
understanding population dynamics for an entire community, assembling
the genome of an organism, and identifying useful biomarkers to classify
a type of process that occurred in a specific environment (Schulz *et
al.*, 2020).

Currently, metagenomics has benefited from technological advances in DNA
sequencing and statistical packages that allow its data’s
reproducibility (Calle, 2019). One of the packages used to reproduce
metagenomics data is metagenomeSeq. The main objective of MetagenomeSeq
is to determine the similarity and abundance of microorganisms in two or
more groups of samples; address subsampling effects normalization in the
detection of diseases, and determine correlations of characteristics. In
this portfolio, data from different samples were analyzed and reproduced
using statistical tools provided from the MetagemeSeq package.

# **2.Meterial and Methods**

The data used in this portfolio were obtained from the **Bioconductor**
[Page](http://www.bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf),
where the first step was the installation of the MetagemeSeq package on
Rstudio.

``` r
library(metagenomeSeq)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    ## Loading required package: glmnet

    ## Loading required package: Matrix

    ## Loaded glmnet 4.0-2

    ## Loading required package: RColorBrewer

In MetagenomeSeq Package, the data must be converted into MRexperiment
objects, so it will be easy to normalize data, run statistical tests and
visualize the results.

To to convert data into MRexperiment objects, the **BIOM Format
package** must be installed because it serves as a bridge to get
MRexperiment-class object.

## 2.1.Load BiomFile

``` r
library(biomformat)
biom_file <- system.file("extdata", "min_sparse_otu_table.biom", package = "biomformat")
b <- read_biom(biom_file)
biom2MRexperiment(b)
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 5 features, 6 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData: none
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

Follows an example of recording an mMR experiment object from a BIOM
file.

``` r
data(mouseData)
b <- MRexperiment2biom(mouseData)
```

## 2.2.Loading Count Data

After pre-processing and annotating the sequencing samples,
MetagenomeSeq requires a counting matrix with samples along the columns
and resources along the lines.

In this example, using a pulmonary microbial, the OTU matrix is stored
as a tab-delimited file. LoadMeta is used to loads the Taxa and counts
as a list.

``` r
dataDirectory <- system.file("extdata", package = "metagenomeSeq")
lung = loadMeta(file.path(dataDirectory, "CHK_NAME.otus.count.csv"))
dim(lung$counts)
```

    ## [1] 1000   78

## 2.3.Loading Taxonomy

To load the annotated taxonomy, it is necessary to ensure that the OTUs
and taxa annotations are in the same order as the matrix rows

``` r
taxa = read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"),
stringsAsFactors = FALSE)
```

## 2.4.Loading metadata

This function provide the data as list.

``` r
clin = loadPhenoData(file.path(dataDirectory, "CHK_clinical.csv"),
tran = TRUE)
ord = match(colnames(lung$counts), rownames(clin))
clin = clin[ord, ]
head(clin[1:2, ])
```

    ##                                          SampleType          SiteSampled
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash Bronchoscope.Channel
    ## CHK_6467_E3B11_OW_V1V2                           OW           OralCavity
    ##                                     SmokingStatus
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
    ## CHK_6467_E3B11_OW_V1V2                     Smoker

## 2.5.Creating new MRExperiment

The newMRexperiment function takes as input, a counting matrix,
featureData, and phenoData. Normalization factors and cover depths are
also input options. Biobase can be used to create annotated data frames.

``` r
phenotypeData = AnnotatedDataFrame(clin)
phenotypeData
```

    ## An object of class 'AnnotatedDataFrame'
    ##   rowNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 CHK_6467_E3B11_OW_V1V2
    ##     ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription

It is possible to get annotated features using taxonomic annotation as
shown below.

``` r
OTUdata = AnnotatedDataFrame(taxa)
OTUdata
```

    ## An object of class 'AnnotatedDataFrame'
    ##   rowNames: 1 2 ... 1000 (1000 total)
    ##   varLabels: OTU Taxonomy ... strain (10 total)
    ##   varMetadata: labelDescription

And then, we can view the MRexperiment data, using the option bellow.

``` r
obj = newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)
obj
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 1000 features, 78 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1 2 ... 1000 (1000 total)
    ##   fvarLabels: OTU Taxonomy ... strain (10 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

## 2.6.Data Sets

Two types of data set were used as an example in the MetagenomeSeq
package. One from the Human Lung microbiome and another from the
Humanized gnotobiotic mouse gut.

The human Lung microbiome data were obtained from samples of the
respiratory flora from six healthy individuals, three smokers and three
non-smokers. Swabs were collected from the oral cavity and
bronchoalveolar lavage.The data were presented in MRexperiment object
format.

``` r
data(lungData)
lungData
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 51891 features, 78 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1 2 ... 51891 (51891 total)
    ##   fvarLabels: taxa
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

For Humanized gnotobiotic mouse gut, twelve germ-free adult male mice
were supplemented with a diet rich in low-fat vegetable polysaccharides.
Each mouse was inoculated with fecal material from a healthy adult
human.After inoculation, the mice remained on the same diet for four
weeks. After four weeks, a subset of six went on a high fat and sugar
diet for eight weeks.

Fecal samples of each mouse were submitted to PCR weekly.

``` r
data(mouseData)
mouseData
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 10172 features, 139 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: PM1:20080107 PM1:20080108 ... PM9:20080303 (139 total)
    ##   varLabels: mouseID date ... status (5 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Prevotellaceae:1 Lachnospiraceae:1 ...
    ##     Parabacteroides:956 (10172 total)
    ##   fvarLabels: superkingdom phylum ... OTU (7 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

## 2.7.Useful Comands

Information of Phenotype can be accessed with the pData ad phenoData
methods:

``` r
phenoData(obj)
```

    ## An object of class 'AnnotatedDataFrame'
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription

``` r
head(pData(obj),3)
```

    ##                                          SampleType          SiteSampled
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash Bronchoscope.Channel
    ## CHK_6467_E3B11_OW_V1V2                           OW           OralCavity
    ## CHK_6467_E3B08_OW_V1V2                           OW           OralCavity
    ##                                     SmokingStatus
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
    ## CHK_6467_E3B11_OW_V1V2                     Smoker
    ## CHK_6467_E3B08_OW_V1V2                  NonSmoker

Feature information can be accessed using the fData and featureData
functions:

``` r
featureData(obj)
```

    ## An object of class 'AnnotatedDataFrame'
    ##   featureNames: 1 2 ... 1000 (1000 total)
    ##   varLabels: OTU Taxonomy ... strain (10 total)
    ##   varMetadata: labelDescription

``` r
head(fData(obj)[, -c(2, 10)], 3)
```

    ##   OTU superkingdom         phylum                  class             order
    ## 1   1     Bacteria Proteobacteria  Epsilonproteobacteria Campylobacterales
    ## 2   2         <NA>           <NA>                   <NA>              <NA>
    ## 3   3     Bacteria Actinobacteria Actinobacteria (class)   Actinomycetales
    ##               family         genus                  species
    ## 1 Campylobacteraceae Campylobacter     Campylobacter rectus
    ## 2               <NA>          <NA>                     <NA>
    ## 3   Actinomycetaceae   Actinomyces Actinomyces radicidentis

The MRcounts function can be used to access the raw or normalized counts
matrix

``` r
head(MRcounts(obj[, 1:2]))
```

    ##   CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 CHK_6467_E3B11_OW_V1V2
    ## 1                                   0                      0
    ## 2                                   0                      0
    ## 3                                   0                      0
    ## 4                                   0                      0
    ## 5                                   0                      0
    ## 6                                   0                      0

We can easily subdivide the MRexperimental-class object as follows:

``` r
featuresToKeep = which(rowSums(obj)>=100)
samplesToKeep = which(pData(obj)$SmokingStatus=="Smoker")
obj_smokers = obj[featuresToKeep,samplesToKeep]
obj_smokers
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 1 features, 33 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (33 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 570
    ##   fvarLabels: OTU Taxonomy ... strain (10 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
head(pData(obj_smokers),3)
```

    ##                                          SampleType          SiteSampled
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash Bronchoscope.Channel
    ## CHK_6467_E3B11_OW_V1V2                           OW           OralCavity
    ## CHK_6467_E3B11_BAL_A_V1V2                     BAL.A                 Lung
    ##                                     SmokingStatus
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
    ## CHK_6467_E3B11_OW_V1V2                     Smoker
    ## CHK_6467_E3B11_BAL_A_V1V2                  Smoker

NormFactors function can be used to access normalization scaling
factors:

``` r
head(normFactors(obj))
```

    ##                                     [,1]
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2   NA
    ## CHK_6467_E3B11_OW_V1V2                NA
    ## CHK_6467_E3B08_OW_V1V2                NA
    ## CHK_6467_E3B07_BAL_A_V1V2             NA
    ## CHK_6467_E3B11_BAL_A_V1V2             NA
    ## CHK_6467_E3B09_OP_V1V2                NA

``` r
normFactors(obj) <- rnorm(ncol(obj))
head(normFactors(obj))
```

    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2              CHK_6467_E3B11_OW_V1V2 
    ##                         -2.98463022                         -0.01278957 
    ##              CHK_6467_E3B08_OW_V1V2           CHK_6467_E3B07_BAL_A_V1V2 
    ##                          0.55527449                          1.45615949 
    ##           CHK_6467_E3B11_BAL_A_V1V2              CHK_6467_E3B09_OP_V1V2 
    ##                         -0.99524266                         -0.98216661

Sequence depth can be accessed using libSize method:

``` r
head(libSize(obj))
```

    ##                                     [,1]
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2    0
    ## CHK_6467_E3B11_OW_V1V2                16
    ## CHK_6467_E3B08_OW_V1V2                 1
    ## CHK_6467_E3B07_BAL_A_V1V2              2
    ## CHK_6467_E3B11_BAL_A_V1V2            118
    ## CHK_6467_E3B09_OP_V1V2                 5

``` r
libSize(obj) <- rnorm(ncol(obj))
head(libSize(obj))
```

    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2              CHK_6467_E3B11_OW_V1V2 
    ##                          0.52896745                          0.49708163 
    ##              CHK_6467_E3B08_OW_V1V2           CHK_6467_E3B07_BAL_A_V1V2 
    ##                          1.24062957                         -0.18027182 
    ##           CHK_6467_E3B11_BAL_A_V1V2              CHK_6467_E3B09_OP_V1V2 
    ##                          0.06437852                          0.94046164

Besides, to preserve a threshold of minimum depth or OTU presence, data
can be filtered

``` r
data(mouseData)
filterData(mouseData,present=10,depth=1000)
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 1057 features, 137 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: PM1:20080108 PM1:20080114 ... PM9:20080303 (137 total)
    ##   varLabels: mouseID date ... status (5 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Erysipelotrichaceae:8 Lachnospiraceae:129 ...
    ##     Collinsella:34 (1057 total)
    ##   fvarLabels: superkingdom phylum ... OTU (7 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

The mergeMRexperiments function can be used to merge two objects of the
MRexperiment class

``` r
data(mouseData)
newobj = mergeMRexperiments(mouseData,mouseData)
```

    ## MRexperiment 1 and 2 share sample ids; adding labels to sample ids.

``` r
newobj
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 10172 features, 278 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: PM1:20080107.x PM1:20080108.x ... PM9:20080303.y (278
    ##     total)
    ##   varLabels: mouseID date ... status (5 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Prevotellaceae:1 Lachnospiraceae:1 ...
    ##     Parabacteroides:956 (10172 total)
    ##   fvarLabels: superkingdom phylum ... OTU (7 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

## 2.8.Normalization

Normalization is essential, as samples can have different depths
coverage. cumNorm is a function used to calculate normalization factors.
Alternatively, Wrench can be used.

Calculating Normalization Factors

``` r
data(lungData)
p=cumNormStatFast(lungData)
```

``` r
lungData = cumNorm(lungData,p=p)
```

Calculating normalization factors using Wrench

wrench function is similar to cumNorm; however, it uses the `condition`
as an argument instead of `p`. `Condition` separates samples according
to the phenotypic groups of interest.

``` r
condition = mouseData$diet
mouseData = wrenchNorm(mouseData,condition=condition)
```

    ## Warning: glm.fit: algorithm did not converge

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: glm.fit: algorithm did not converge

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

    ## Warning: Partial NA coefficients for 8430 probe(s)

## 2.9.Exporting data

Normalized count matrices can be exported using the following commands:

``` r
mat = MRcounts(lungData, norm = TRUE, log = TRUE)[1:5, 1:5]
exportMat(mat, file = file.path(dataDirectory, "tmp.tsv"))
```

To save the statistic of exported data (library size,quantile value,
sample scaling factor, number of identified features) `exportStats`
function

``` r
exportStats(lungData[,1:5],file=file.path(dataDirectory,"tmp.tsv"))
```

    ## Default value being used.

``` r
head(read.csv(file=file.path(dataDirectory,"tmp.tsv"),sep="\t"))
```

    ##                               Subject Scaling.factor Quantile.value
    ## 1 CHK_6467_E3B11_BRONCH2_PREWASH_V1V2             67              2
    ## 2              CHK_6467_E3B11_OW_V1V2           2475              1
    ## 3              CHK_6467_E3B08_OW_V1V2           2198              1
    ## 4           CHK_6467_E3B07_BAL_A_V1V2            836              1
    ## 5           CHK_6467_E3B11_BAL_A_V1V2           1008              1
    ##   Number.of.identified.features Library.size
    ## 1                            60          271
    ## 2                          3299         7863
    ## 3                          2994         8360
    ## 4                          1188         5249
    ## 5                          1098         3383

``` r
system(paste("rm",file.path(dataDirectory,"tmp.tsv")))
```

    ## Warning in system(paste("rm", file.path(dataDirectory, "tmp.tsv"))): 'rm' not
    ## found

    ## [1] 127

## 3.Statistical testing

After normalization, we can evaluate subsampling effects on detecting
the differentially abundant characteristic. And for that, we can use
**fitFeatureModel** or **fitzig**. The **MRfulltable**, **MRcoefs**,
**MRtable** are summary tables for outputs.

follows an example of pulmonary microbiome comparison between smokers
and nonsmokers

``` r
data(lungData)
lungData = lungData[,-which(is.na(pData(lungData)$SmokingStatus))]
lungData=filterData(lungData,present=30,depth=1)
lungData <- cumNorm(lungData, p=.5)
pd <- pData(lungData)
mod <- model.matrix(~1+SmokingStatus, data=pd)
lungres1 = fitFeatureModel(lungData,mod)
head(MRcoefs(lungres1))
```

    ##           logFC        se      pvalues   adjPvalues
    ## 3465  -4.824949 0.5697511 0.000000e+00 0.000000e+00
    ## 35827 -4.304266 0.5445548 2.664535e-15 1.079137e-13
    ## 2817   2.320656 0.4324661 8.045793e-08 1.629273e-06
    ## 2735   2.260203 0.4331098 1.803341e-07 2.921412e-06
    ## 5411   1.748296 0.3092461 1.572921e-08 4.246888e-07
    ## 48745 -1.645805 0.3293117 5.801451e-07 7.831959e-06

Using **fitZig** for differential abundance testing

The user must restrict significant resources. In the pulmonary
microbiome analysis, we can remove controls, and characteristics absent
in many samples, and calculate the normalization factors.

``` r
data(lungData)
controls = grep("Extraction.Control",pData(lungData)$SampleType)
lungTrim = lungData[,-controls]
rareFeatures = which(rowSums(MRcounts(lungTrim)>0)<10)
lungTrim = lungTrim[-rareFeatures,]
lungp = cumNormStat(lungTrim,pFlag=TRUE,main="Trimmed lung data")
```

    ## Default value being used.

![](Metagenome_Seq_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
lungTrim = cumNorm(lungTrim,p=lungp)
```

After defining an appropriate model matrix, there are optional entries
for *fitZig*, including configurations by zigControl. This example
includes body site as covariates, and we want to test differentially
abundant bacteria between smokers and non-smokers.

``` r
smokingStatus = pData(lungTrim)$SmokingStatus
bodySite = pData(lungTrim)$SampleType
normFactor = normFactors(lungTrim)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(~smokingStatus+bodySite + normFactor)
settings = zigControl(maxit=10,verbose=TRUE)
fit = fitZig(obj = lungTrim,mod=mod,useCSSoffset = FALSE, 
             control=settings)
```

    ## it= 0, nll=88.42, log10(eps+1)=Inf, stillActive=1029
    ## it= 1, nll=93.56, log10(eps+1)=0.06, stillActive=261
    ## it= 2, nll=93.46, log10(eps+1)=0.05, stillActive=120
    ## it= 3, nll=93.80, log10(eps+1)=0.05, stillActive=22
    ## it= 4, nll=93.94, log10(eps+1)=0.03, stillActive=3
    ## it= 5, nll=93.93, log10(eps+1)=0.00, stillActive=1
    ## it= 6, nll=93.90, log10(eps+1)=0.00, stillActive=1
    ## it= 7, nll=93.87, log10(eps+1)=0.00, stillActive=1
    ## it= 8, nll=93.86, log10(eps+1)=0.00, stillActive=1
    ## it= 9, nll=93.85, log10(eps+1)=0.00, stillActive=1

Running fitZig by default, the covariate adjustment must be added to the
design matrix, as it is crucial for contrast.

``` r
settings = zigControl(maxit=1,verbose=FALSE)
mod = model.matrix(~bodySite)
colnames(mod) = levels(bodySite)
res = fitZig(obj = lungTrim,mod=mod,control=settings)
zigFit = slot(res,"fit")
finalMod = slot(res,"fit")$design
contrast.matrix = makeContrasts(BAL.A-BAL.B,OW-PSB,levels=finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2)
```

    ##       BAL.A...BAL.B  OW...PSB   AveExpr         F      P.Value  adj.P.Val
    ## 18531    0.37318792  2.075648 0.7343081 12.715105 5.359780e-05 0.02813711
    ## 6291    -0.10695735  1.658829 0.4671470 12.956898 5.482439e-05 0.02813711
    ## 37977   -0.37995461  2.174071 0.4526060 12.528733 8.203239e-05 0.02813711
    ## 6901     0.17344138  1.466113 0.2435881 12.018652 1.335806e-04 0.03212047
    ## 40291    0.06892926  1.700238 0.2195735 11.803380 1.560761e-04 0.03212047
    ## 36117   -0.28665883  2.233996 0.4084024 10.571931 3.012092e-04 0.05013569
    ## 7343    -0.22859078  1.559465 0.3116465 10.090602 3.931844e-04 0.05013569
    ## 7342     0.59882970  1.902346 0.5334647  9.410984 4.901651e-04 0.05013569
    ## 1727     1.09837459 -2.160466 0.7780167  9.346013 5.027597e-04 0.05013569
    ## 40329   -0.07145998  1.481582 0.2475735  9.700136 5.259032e-04 0.05013569

To consider characteristics in abundance, the option `MRcoefs` can be
used in the MR tables in a specific group.

``` r
taxa = 
  sapply(strsplit(as.character(fData(lungTrim)$taxa),split=";"),
         function(i){i[length(i)]})
head(MRcoefs(fit,taxa=taxa,coef=2))
```

    ##                                   smokingStatusSmoker      pvalues   adjPvalues
    ## Neisseria polysaccharea                     -4.031612 3.927097e-11 2.959194e-08
    ## Neisseria meningitidis                      -3.958899 5.751592e-11 2.959194e-08
    ## Prevotella intermedia                       -2.927686 4.339587e-09 8.930871e-07
    ## Porphyromonas sp. UQD 414                   -2.675306 1.788697e-07 1.357269e-05
    ## Prevotella paludivivens                      2.575672 1.360718e-07 1.272890e-05
    ## Leptotrichia sp. oral clone FP036            2.574172 3.544957e-04 1.414122e-03

Looking for this output, we can observe two **Neisseria**, two
**Prevotella**, a **Leptotrichia**, and **Porphyromonas** are
differentially abundant.

The `coef` parameter refers to the coefficient of interest to be tested.
We can see the previous model using the `fitLogNormal` parameter, which
provides the p-value resolution for ten.

``` r
coeffOfInterest = 2
res = fitLogNormal(obj = lungTrim, mod = mod, useCSSoffset = FALSE,
B = 10, coef = coeffOfInterest)
adjustedPvalues = p.adjust(res$p, method = "fdr")
foldChange = abs(res$fit$coef[, coeffOfInterest])
sigList = which(adjustedPvalues <= 0.05)
sigList = sigList[order(foldChange[sigList])]
head(taxa[sigList])
```

    ## [1] "Megasphaera micronuciformis"         "Prevotella genomosp. C2"            
    ## [3] "Anaeroglobus geminatus"              "Veillonella montpellierensis"       
    ## [5] "Capnocytophaga sp. oral clone DZ074" "Macrococcus caseolyticus"

Making adjustments on Pvalue, we can observe other bacterial species of
the pulmonary flora.

## 3.1.Presence Absence

The presence-absence test’s idea is to assess whether a particular
characteristic is in a greater/lesser proportion between groups of
individuals. The `fitPA` parameter can be used to calculate the
presence-absence for each organism.

``` r
classes = pData(mouseData)$diet
res = fitPA(mouseData[1:5,],cl=classes)
classes = pData(mouseData)$diet
res = fitDO(mouseData[1:100,],cl=classes,norm=FALSE,log=FALSE)
head(res)
```

    ##                         oddsRatio      lower    upper   pvalues adjPvalues
    ## Prevotellaceae:1              Inf 0.01630496      Inf 1.0000000  1.0000000
    ## Lachnospiraceae:1             Inf 0.01630496      Inf 1.0000000  1.0000000
    ## Unclassified-Screened:1       Inf 0.01630496      Inf 1.0000000  1.0000000
    ## Clostridiales:1                 0 0.00000000 24.77661 0.3884892  0.7470946
    ## Clostridiales:2               Inf 0.01630496      Inf 1.0000000  1.0000000
    ## Firmicutes:1                    0 0.00000000 24.77661 0.3884892  0.7470946

# 3.2.Feature correlations

`CorrelationTest` can be used to correlate the classes of
microorganisms.

``` r
cors = correlationTest(mouseData[55:60, ], norm = FALSE, log = FALSE)
head(cors)
```

    ##                                     correlation            p
    ## Clostridiales:11-Lachnospiraceae:35 -0.02205882 7.965979e-01
    ## Clostridiales:11-Coprobacillus:3    -0.01701180 8.424431e-01
    ## Clostridiales:11-Lactobacillales:3  -0.01264304 8.825644e-01
    ## Clostridiales:11-Enterococcaceae:3   0.57315130 1.663001e-13
    ## Clostridiales:11-Enterococcaceae:4  -0.01264304 8.825644e-01
    ## Lachnospiraceae:35-Coprobacillus:3   0.24572606 3.548360e-03

## 3.3.Unique OTUs or features

The `uniqueFeatures` function is used To track missing resources in any
number of classes; The `uniqueFeatures` function provides a table of
resource ids, the number of positive resources, and readings for each
group.

``` r
cl = pData(mouseData)[["diet"]]
uniqueFeatures(mouseData,cl,nsamples = 10,nreads = 100)
```

    ##                      featureIndices Samp. in BK Samp. in Western Reads in BK
    ## Enterococcaceae:28             2458           0               36           0
    ## Lachnospiraceae:1453           2826          16                0         192
    ## Firmicutes:367                 4165           0               32           0
    ## Prevotellaceae:143             7030          50                0         109
    ## Lachnospiraceae:4122           7844          15                0         505
    ## Enterococcus:182               8384           0               34           0
    ## Lachnospiraceae:4347           8668          50                0         130
    ## Prevotella:79                  8749          50                0         100
    ## Prevotella:81                  8994          71                0         370
    ## Prevotellaceae:433             9223          53                0         154
    ##                      Reads in Western
    ## Enterococcaceae:28                123
    ## Lachnospiraceae:1453                0
    ## Firmicutes:367                    112
    ## Prevotellaceae:143                  0
    ## Lachnospiraceae:4122                0
    ## Enterococcus:182                  163
    ## Lachnospiraceae:4347                0
    ## Prevotella:79                       0
    ## Prevotella:81                       0
    ## Prevotellaceae:433                  0

## 3.4 Aggregating counts

The `ggTax` function in an MRexperiment can be used to aggregate the
counting matrix (normalized or not) at the user’s desired level

``` r
obj = aggTax(mouseData,lvl='phylum',out='matrix')
head(obj[1:5,1:5])
```

    ##                PM1:20080107 PM1:20080108 PM1:20080114 PM1:20071211 PM1:20080121
    ## Actinobacteria            0            3            2           37            0
    ## Bacteroidetes           486          921         1103          607          818
    ## Cyanobacteria             0            0            0            0            0
    ## Firmicutes              455          922         1637          772         1254
    ## NA                        5           25            5            8            4

Sample aggregation can also be done using the `aggregateBySample` or
`aggsamp` function, selecting the column of interest.

``` r
obj = aggSamp(mouseData,fct='mouseID',out='matrix')
head(obj[1:5,1:5])
```

    ##                         PM1 PM10       PM11 PM12 PM2
    ## Prevotellaceae:1          0    0 0.00000000    0   0
    ## Lachnospiraceae:1         0    0 0.00000000    0   0
    ## Unclassified-Screened:1   0    0 0.08333333    0   0
    ## Clostridiales:1           0    0 0.00000000    0   0
    ## Clostridiales:2           0    0 0.00000000    0   0

## **4.Results**

To visualize the analyzed data, metagenomeSeq has several plotting
functions.

To access the abundance heatmap we use the `plotMRheatmap` function:

``` r
trials = pData(mouseData)$diet
heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
```

``` r
trials = pData(mouseData)$diet
heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj=mouseData,n=200,cexRow = 0.4,cexCol = 0.4,trace="none",
                col = heatmapCols,ColSideColors = heatmapColColors)
```

![**Figure1.** Taxonomic comparison of all samples: Red values indicate
counts close to zero. The Row color labels indicate OTU taxonomic class;
column color labels indicate diet (green = high fat, yellow = low
fat)](Metagenome_Seq_files/figure-gfm/unnamed-chunk-38-1.png)

For basic correlation structures we use `plotCorr` function

``` r
plotCorr(obj = mouseData, n = 200, cexRow = 0.25, cexCol = 0.25,
trace = "none", dendrogram = "none", col = heatmapCols)
```

![**Figure2.** Taxonomic
correlation](Metagenome_Seq_files/figure-gfm/unnamed-chunk-39-1.png)

For principal coordinate analyses (PcoA) we use `plotOrd` ,and for
rarefaction effects `plotRare`

``` r
cl = factor(pData(mouseData)$diet)
plotOrd(mouseData,tran=TRUE,usePCA=FALSE,useDist=TRUE,bg=cl,pch=21)
```

![**Figure3.** principal coordinate analyses (PCoA) comparison from
groups submitted to different
diets)](Metagenome_Seq_files/figure-gfm/unnamed-chunk-40-1.png)

``` r
res = plotRare(mouseData,cl=cl,pch=21,bg=cl)
tmp=lapply(levels(cl), function(lv) 
  lm(res[,"ident"]~res[,"libSize"]-1, subset=cl==lv))
for(i in 1:length(levels(cl))){
   abline(tmp[[i]], col=i)
}
legend("topleft", c("Diet 1","Diet 2"), text.col=c(1,2),box.col=NA)
```

![**Figure4.** Rarefaction
effect](Metagenome_Seq_files/figure-gfm/unnamed-chunk-41-1.png)

## Feature specific

In the MetagenomeSeq package, `plotOTU` was used to plot the *Neisseria
meningitidis* normalized log (cpt) present in the 779th line of the
lungTrim counting matrix. And `plotGenus` was used to plot the
normalized log (cpt) of all *Neisseria meningitidis* annotated.

``` r
head(MRtable(fit,coef=2,taxa=1:length(fData(lungTrim)$taxa)))
```

    ##     +samples in group 0 +samples in group 1 counts in group 0 counts in group 1
    ## 63                   24                   6              1538                11
    ## 779                  23                   7              1512                22
    ## 358                  24                   1               390                 1
    ## 499                  21                   2               326                 2
    ## 25                   15                  26               162              1893
    ## 928                   2                  11                 4                91
    ##     smokingStatusSmoker      pvalues   adjPvalues
    ## 63            -4.031612 3.927097e-11 2.959194e-08
    ## 779           -3.958899 5.751592e-11 2.959194e-08
    ## 358           -2.927686 4.339587e-09 8.930871e-07
    ## 499           -2.675306 1.788697e-07 1.357269e-05
    ## 25             2.575672 1.360718e-07 1.272890e-05
    ## 928            2.574172 3.544957e-04 1.414122e-03

``` r
patients=sapply(strsplit(rownames(pData(lungTrim)),split="_"),
          function(i){
            i[3]
          })
pData(lungTrim)$patients=patients
classIndex=list(smoker=which(pData(lungTrim)$SmokingStatus=="Smoker"))
classIndex$nonsmoker=which(pData(lungTrim)$SmokingStatus=="NonSmoker")
otu = 779
plotOTU(lungTrim,otu=otu,classIndex,main="Neisseria meningitidis")
```

![**Figure5.** Abundance between two groups of
comparison](Metagenome_Seq_files/figure-gfm/unnamed-chunk-42-1.png)

``` r
x = fData(lungTrim)$taxa[otu]
otulist = grep(x,fData(lungTrim)$taxa)
```

``` r
plotGenus(lungTrim, otulist, classIndex, labs = FALSE, main = "Neisseria meningitidis")
lablist <- c("S", "NS")
axis(1, at = seq(1, 6, by = 1), labels = rep(lablist, times = 3))
```

![**Figure6**. Multiple OTU abundances in groups of
comparison](Metagenome_Seq_files/figure-gfm/unnamed-chunk-43-1.png)

``` r
classIndex = list(Western = which(pData(mouseData)$diet == "Western"))
classIndex$BK = which(pData(mouseData)$diet == "BK")
otuIndex = 8770
dates = pData(mouseData)$date
plotFeature(mouseData, norm = FALSE, log = FALSE, otuIndex, classIndex,
col = dates, sortby = dates, ylab = "Raw reads")
```

![**Figure7.** Raw reads
abundances](Metagenome_Seq_files/figure-gfm/unnamed-chunk-44-1.png)![**Figure7.**
Raw reads
abundances](Metagenome_Seq_files/figure-gfm/unnamed-chunk-44-2.png)

# **4.Discussion**

## Data Mouse

Analyzing data mouse from MetagenomeSeq package, we can see that mice
supplemented with a low-fat diet have more abundance of Lachnospiraceae
family. In contrast, mice supplemented with a high-fat diet have
abundance of Lactococcus, followed by Lachnospiraceae (see Figure 1).
All families have a similar correlation on the diagonal line(Figure 2).

The Analises of PcoA compares the variance in 2 different groups of
diet. The diet1 is represented on the Top on the left side, and Diet2 is
on the top too, on the right side. Analyzing CP1 (on the x-axis), We can
clearly see a variance in diet2, creating a separation from Diet1.
Analyzing Cp2 (Y-axis), Diet1 tend to be close together and we don’t see
much variance with Diet2. (Figure 3).

In Figure 4, the population of microorganisms from Diet2 needs a small
number of readings to obtain high associated characteristics; In
contrast, the group of microorganisms from Diet1 needs many readings to
achieve high associated resources (Figure 4).

## Data from Human lung microbiome

In the human lung microbiome, **Neisseria meningitidis** was
differentially more abundant in nonsmokers (Figure5 and 6).

# **5.Conclusion**

MetagenomeSeq is a useful package and is important to understand
variances in a community of microorganisms due to phenotypic
differences. Using MetagemeSeq we can do differential abundance analyses
to evaluate microorganisms’ similarities and genetic diversity.
Additionally, we can identify species involved in various diseases, and
we can outline preventive measures and treatment strategies para for
diseases.

## **References**

1.Calle, M. L. (2019) ‘Statistical Analysis of Metagenomics Data’,
Genomics & informatics, 17(1), pp. e6-e6.

2.Schulz, F., Andreani, J., Francis, R., Boudjemaa, H., Bou Khalil, J.
Y., Lee, J., La Scola, B. and Woyke, T. (2020) ‘Advantages and Limits of
Metagenomic Assembly and Binning of a Giant Virus’, mSystems, 5(3),
pp. e00048-20.

3.Sleator, R. D., Shortall, C. and Hill, C. (2008) ‘Metagenomics’,
Letters in applied microbiology, 47(5), pp. 361-366.
