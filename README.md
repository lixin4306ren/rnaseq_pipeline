# 基于genome alignment的RNA-Seq分析流程
## 安装
此流程依赖于array suite，安装array suite参考此[网页](http://www.arrayserver.com/wiki/index.php?title=Oshell#Overview)。

## 使用
主程序为`Run_OSA_all_modules.pl`，主要参数用法：
```
perl Run_OSA_all_modules.pl sample.list [output_dir] [database_dir] [genome_version] [gene_model_version]
```

`sample.list`为所有样品目录路径，目录下为fastq.gz文件，格式如下：  
```
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S1
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S2
```
`genome_version`默认基因组版本`Human.B37.3`  
`gene_model_version`默认版本为`RefGene`  
基因组和基因注释版本对应信息参考此[网页](http://www.arrayserver.com/wiki/index.php?title=A_list_of_compiled_genome_and_gene_model_from_OmicSoft)


## 举例
如果是`human`  
```
perl ~/scripts/rnaseq_pipeline/Run_OSA_all_modules.pl sample.list
```
如果是其他物种，需指明所用genome版本和gene model版本，例如`mouse`  
```
perl ~/scripts/rnaseq_pipeline/Run_OSA_all_modules.pl sample.list ./ /home/others/xli/data/genome Mouse.mm10 Ensembl.R80
```
会生成`align.sh`文件和响应的`*.oscript`文件。用`qsub`提交`align.sh`运行alignment步骤。完成后所有结果在`output_dir/sample_name`目录。



# kmer-based alignment free的RNA-Seq分析流程
本流程基于软件[kallisto](https://pachterlab.github.io/kallisto/)。
## index genome
可以此[链接](http://bio.math.berkeley.edu/kallisto/transcriptomes/)下载常用的cDNA序列。
```
kallisto index -i transcripts.idx transcripts.fasta.gz
```
## quantification
```
kallisto quant -i ~/data/cDNA/Homo_sapiens.GRCh38.rel79.cdna.all.idx -o S1 -t 2 -b 100 ~/raw_data/2016_08_01_ZY/RNA-Seq/S1/S1.R1.fastq.gz ~/raw_data/2016_08_01_ZY/RNA-Seq/S1/S1.R2.fastq.gz
```
此步会生成`abundance.tsv`文件，包含各个转录本的表达量。

## 生成gene level的FPKM
使用R package [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
```
library(tximport)
dir<-"./sleuth"                   #sleuth目录下保存上步各个样品的kallisto运行结果
samplePath<-"./sample.list"       #所有样品list
read.table(samplePath,header=T)->sampleInfor
files <- file.path(dir,  sampleInfor[,1], "abundance.tsv")
names(files) <-  as.vector(sampleInfor[,1])
load(t2gPath)
txi <- tximport(files, type = "kallisto", tx2gene = t2g, reader = read_tsv)
save(txi,file="txi.rda")          #txi.rda保存结果用于后续基因表达差异分析
```

# 使用rmarkdown生成html报告
以`/home/others/xli/RNA-Seq2`目录下数据为例。
安装`knitr`包：
```
install.packages('knitr', dependencies = TRUE)
```
安装其他相关包：
```
install.packages("ggplot2")
install.packages("readr")
install.packages("stringr")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("DT")

source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
biocLite("DESeq2")
biocLite("genefilter")
biocLite("AnnotationDbi")
biocLite("org.Mm.eg.db") #小鼠基因注释文件
biocLite("EnrichmentBrowser")
biocLite("ReportingTools")
biocLite("KEGGgraph")
biocLite("org.Mm.eg.db")
```

准备样品信息文件`sample2.list`，格式如下：
```
sample 	condition
lane1-1	T
lane1-2	T
lane1-3	T
lane1-4	W
lane1-5	W
lane1-6	W
```

运行rmarkdown报告模版：
```
/opt/apps/R/bin/R -e "rmarkdown::render('rnaseq_report.Rmd')"
```

