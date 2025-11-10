# All-stage OA classification
All-stage OA classification identified molecular subtypes of OA and further validated certain molecular subtype as therapeutic subtypes via clinical trial. The code available in this repository could reproduce the results in our paper: 1. Quality check, 2. All-stage OA classification, 3. Transcriptome characteristics, 4. OAClassifier, 5.Clinical phenotypes, 6.Drug repurposing, 7.Mediation analysis, 8.Proteomics. 

## 1. Quality check
Quality Check (QC) is essential for ensuring the integrity and reliability of the data. This step involves detecting and removal outliers using PCA, and assessing the distribution of gene expression data to determine the quality of the dataset.

## 2. All-stage OA classification
Spectra clustering was used to classify all-stage OA patients ,including wild,middle and severe, based on ssGSEA scores.

## 3. Transcriptome characteristics
Characterize OA subtype based on subtype-specific pathways related to OA pathology,drivers or symptoms.

## 4. OAClassifier
To extrapolate the all-stage OA classification system to other datasets and compare different OA classification systems, we developed aOAClassifier. This tool is designed to identify all-stage OA subtypes from validation and external datasets based on the random forest algorithm.

## 5. Clinical phenotypes
Firstly, we analyzed the differences in basic information (including age, sex, affected side, and OA stage based on the Kellgren-Lawrence (KL) grade) among the subtypes. Additionally, we identified modifiable and recognizable risk factors (e.g., obesity and trauma) specific to certain subtype, which allows for the implementation of targeted prevention—such as weight loss—to prevent the development of obesity-related OA subtypes.

## 6. Drug repurposing
Subtype-specific drug candidates were obtained by integrating subtype-specific differentially expressed genes (DEGs) and predicted drug-induced gene expression profiles [available at https://drive.google.com/uc?export=download&id=1clsvmAZPVeQMAwxuBRRxTJaxEGXLqTXa; Pham et al., 2021, Nature Machine Intelligence] via computational drug repurposing algorithms [Chen et al., 2018, Gastroenterology].

## 7.Mediation analysis
The Expecto [Zhou J, et al., 2018. Nature Genetics] was used for the identification of subtype-specific transcriptional effect variants. After prediction, the PRDM6 eQTL was predicted to be specific to the ECM Degradation Mixed Pain subtype. Mendelian randomization reveals PRDM6 eQTLs have a significantly detrimental effect on chronic knee pain. Subsequnent mediation analysis showed that vitamin D supplements have causally associated with PRDM6 eQTLs and chronic knee pain.

## 8.Proteomics
Change of protein quantification values in synovial fluid of the individuals at baseline and hyaluromic acid treatment.

## 7. Contact

**Xianan Mo** < 12318017@zju.edu.cn >

Department of School of Medicine, Zhejiang University University, China

