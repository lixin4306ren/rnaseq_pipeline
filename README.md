# RNA-Seq分析流程
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
```
perl ~/scripts/rnaseq_pipeline/RNA_seq_align_all_module.pl -sample_list sample.list -out_dir /home/others/xli/RNA-Seq -step align
```
生成`align.sh`文件和响应的`*.oscript`文件。用`qsub`提交`align.sh`运行alignment步骤。完成后所有结果在`/home/others/xli/RNA-Seq/sample_name`目录。

## 合并bam文件去duplicates

## 生成可供IGV broswer浏览的tdf文件

