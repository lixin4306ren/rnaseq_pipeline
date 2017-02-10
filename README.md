# 基于genome alignment的RNA-Seq分析流程
## 安装
此流程依赖于array suite，安装array suite参考此[网页](http://www.arrayserver.com/wiki/index.php?title=Oshell#Overview)。

## 主要参数
主程序为`RNA_seq_align_all_module.pl`，用`-h`显示所有参数如下：
```
        RNA-Seq pipeline

        Usage: RNA_seq_align_all_module.pl <options>

                -sample_list        sample list

                -out_dir            output result folder

                -ref_dir            reference OSA base path

                -ref                reference name

                -gene_model         gene_model name

                -gtf                gene gtf file

                -gtf_exon           exon gff file for DEXSeq

                -step               align or merge_rmdup or count or count_exon or tdf

                -h or -help         Show Help , have a choice
```


`sample.list`为所有样品目录路径，目录下为fastq.gz文件，格式如下：  
```
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S1
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S2
```
`-ref`默认基因组版本Human.B37.3  
`-gene_model`默认版本为RefGene  
基因组和基因注释版本对应信息参考此[网页](http://www.arrayserver.com/wiki/index.php?title=A_list_of_compiled_genome_and_gene_model_from_OmicSoft)


## alignment
如何是human  
```
```
如果是其他物种，需指明所用genome版本和gene model版本，例如mouse  
```
perl ~/scripts/rnaseq_pipeline/Run_OSA_all_modules.pl sample.list ./ /home/others/xli/data/genome Mouse.mm10 Ensembl.R80
```
生成`align.sh`文件和响应的`*.oscript`文件。用`qsub`提交`align.sh`运行alignment步骤。完成后所有结果在`/home/others/xli/RNA-Seq/sample_name`目录。



## kmer-baased alignment free的RNA-Seq分析流程
本流程基于软件[kallisto](https://pachterlab.github.io/kallisto/)。
### index genome
可以此[链接](http://bio.math.berkeley.edu/kallisto/transcriptomes/)下载常用的cDNA序列。
```
kallisto index -i transcripts.idx transcripts.fasta.gz
```
### quantification
```
kallisto quant -i ~/data/cDNA/Homo_sapiens.GRCh38.rel79.cdna.all.idx -o S1 -t 2 -b 100 ~/raw_data/2016_08_01_ZY/RNA-Seq/S1/S1.R1.fastq.gz ~/raw_data/2016_08_01_ZY/RNA-Seq/S1/S1.R2.fastq.gz
```



