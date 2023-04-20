# *Fgf15* FACS

Analyzing genes changed with diet in *Fgf15*-positive cells from ileum.

## Files
```
+-- data
    +-- Run_2425
        +-- Sample_113433
            +-- ReadsPerGene.out.tab
        +-- Sample_113434
            +-- ReadsPerGene.out.tab
        +-- ...
        +-- Run_2425_seeley.csv
    +-- Run_2643
        +-- Sample_123057
            +-- ReadsPerGene.out.tab
        +-- Sample_123058
            +-- ReadsPerGene.out.tab
        +-- ...
        +--
        +-- Run_2643_seeley.csv
    +-- genes.csv
+-- scripts
    +-- analysis.R
    +-- mapping.sh
+-- environment.yml
+-- README.md
```


## Pipeline
### Setup
1. Create `anaconda` environment: `conda env create -f environment.yml`
2. Activate `anaconda` environment: `conda activate Fgf15`

### Analysis
1. Map FASTQ files: `bash scripts/map.sh` (note this does not work in the GitHub repo because FASTQ files are not included).
2. Analyze count files: `Rscript scripts/analyze.R`

## Results
All results found in the **results** folder
* **results.csv**: numeric results from analysis pipeline
* **kegg.csv**: KEGG enrichment results
* **plots/kegg.png**: figure from the paper
