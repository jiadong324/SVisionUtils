
<div align=left><img width=30% height=30% src="https://github.com/xjtu-omics/SVision/blob/master/supports/svision-logo.png"/></div>

This repository provides support for SVision downstream CSV filter and analysis.


### Usage

#### Dependencies

Please install pandas, numpy and [intervaltree](https://pypi.org/project/intervaltree/). 

The call set for the paper is ./supports/HG00733.svision.s5.graph.vcf.

#### Prepare config file

**NOTE:** Default values in the config file are used to produce results in the paper

The config file requires: 

1. Chromosomes of interest. Default value constains the autosomes.
2. The path to bedtools.
3. The path to RepeatMasker and Tandem Repeat Finder annotated human reference genome GRCh38. Please download [TRF](https://drive.google.com/file/d/17w_aKnAsU1dFxSmwQ3PfIArENjMOiZys/view?usp=sharing) 
   and [RMSK](https://drive.google.com/file/d/1QCpQLIEP-b0ApL2HwPdl1ehD6thxyZpq/view?usp=sharing) in BED format.
   
4. Regions to exclude in the filter. A [BED](https://github.com/jiadong324/SVisionUtils/blob/master/supports/grch38.exclude_regions_cen.bed) file is avaiable.
5. The path to reference genome used in SV detection. 

#### Run filter


```
python FilterMain.py supports/HG00733.svision.s5.graph.vcf -g ./supports/HG00733.graph_exactly_match.txt -w ./output_dir -i 0,3
```

This will generate three files:

*prefix*.filtered.vcf: SVision discoveries filtered by low mapping quality regions, gaps and centromeres.

*prefix*.Raw-CSVs.tsv: SVision CSVs filtered by graph structures.

*prefix*.HQ-CSVs.tsv: CSVs additionally filtered by tandem repeats.

