# RNA-Seq分析流程
## 安装
此流程用OSA进行比对，安装OSA参考此[网页](http://www.arrayserver.com/wiki/index.php?title=Oshell#Overview)。

sample.list格式为所有样品目录路径，目录下为fastq文件
```
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S1
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S2
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S3
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S4
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S5
/home/others/xli/raw_data/2016_08_01_ZY/RNA-Seq/S6
```

nohup /usr/local/bin/bcl2fastq --runfolder-dir /home/yliu/data/Xin/2016_08_01_ZY/ --output-dir 2016_08_01_ZY_100bp --sample-sheet 2016_08_01_ZY.csv --use-bases-mask Y100n51,I6,Y100n51 

mono /home/others/xli/soft/oshell/16-8-11/oshell.exe --runscript /home/others/xli/data/genome /home/others/xli/RNA-Seq/Alignment.oscript /tmp/ /opt/apps/mono/bin/mono


perl ~/scripts/RNA_seq_align_all_module.pl -sample_list sample.list -step align

perl ~/scripts/RNA_seq_align.pl -sample_list sample.list -ref_dir ~/amber3/no_back_up/data/genome/hg19_OSA/ -ref Human.B37 -gene_model RefGene -step align

perl ~/scripts/RNA_seq_align.pl -sample_list sample.list -ref_dir ~/amber3/no_back_up/data/genome/galgal4/for_OSA/ -ref gal4 -gene_model gal4gene -step merge
perl ~/scripts/RNA_seq_align.pl -sample_list sample.list -gtf ~/amber3/no_back_up/data/gene/chicken/Gallus_gallus.Galgal4.71.gtf -step count
perl ~/scripts/RNA_seq_align.pl -sample_list sample.list -gtf_exon ~/amber3/no_back_up/data/gene/chicken/Gallus_gallus.Galgal4.71.for.dexseq.gff -step count_exon
perl ~/scripts/RNA_seq_align.pl -sample_list sample.list -ref_dir ~/amber3/no_back_up/data/genome/galgal4/for_OSA/ -ref gal4 -gene_model gal4gene -step align_sum
perl ~/scripts/RNA_seq_align.pl -sample_list sample.list -ref_dir ~/amber3/no_back_up/data/genome/galgal4/for_OSA/ -step tdf

