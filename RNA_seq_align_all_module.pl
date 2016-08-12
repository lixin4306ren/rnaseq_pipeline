#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;

my ($sample_list,$ref,$help,$step,$ref_dir,$gene_model,$gtf,$gtf_exon,$out_dir);

my $Function='RNA-Seq pipeline';
#my $bash="#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -q batch\n";

GetOptions(
        "sample_list:s"=>\$sample_list,
	"out_dir:s"=>\$out_dir,
        "ref_dir:s"=>\$ref_dir,
		"ref:s"=>\$ref,
		"gene_model:s"=>\$gene_model,
		"gtf:s"=>\$gtf,
		"gtf_exon:s"=>\$gtf_exon,
        "help"=>\$help,
        "step:s"=>\$step,

);

if(!defined($sample_list) ||!defined($step)||defined($help) ){

        Usage();

}
if(!defined $out_dir){
	$out_dir=`pwd`;
}

if (!defined $ref) {
$ref_dir="/home/others/xli/data/genome";
$ref="Human.B37.3";
$gene_model="RefGene";
}

if ($step eq 'align') {
open O,">align.sh"||die;
chomp($out_dir);
my $cmd="perl ~xli/scripts/Run_OSA_all_modules.pl $sample_list $out_dir $ref_dir $ref $gene_model";
my $cmd2=`$cmd`;
print O $cmd2,"\n";

#my $tmp_cmd="perl /home/jhmi/xinli/scripts/qsub_sge.pl --resource cegs,mf=30G,h_vmem=40G align.sh";
#system($tmp_cmd);

}
elsif($step eq 'align_sum'){
my %pair_mapped_read;
my %single_mapped_read;
my %total_read;
my %anno_read;
my %total_read2;
my %total_gene_exp;
open IN,$sample_list;
while (<IN>) {
	chomp;
	opendir DIR,$_;
	my $sample_name=(split /\//,$_)[-1];
	my $count_file=$sample_name.".count";
	open TMP,$count_file||die;
	while (my $line=<TMP>) {
		chomp $line;
		#print "$line\n";exit;
		my @infor=split /\s+/,$line;
		if ($line=~/^EN/) {
			$anno_read{$sample_name}+=$infor[1];
			if ($infor[1]>=10) {
				$total_gene_exp{$sample_name}++;
			}
		}
		$total_read2{$sample_name}+=$infor[1];
	}
	close TMP;

	#print "$sample_name\n";exit;
	while (my $file=readdir(DIR)) {
		if ($file=~/.*_R1_.*fastq.gz$/ or $file=~/\.fastq$/) {
			my $res_dir=(split /\//,$file)[-1];
			my $res_file=$res_dir."/".$res_dir.".AlignmentSummary.txt";
			open TMP,$res_file||die;
			while (my $line=<TMP>) {
				my @infor=split /\s+/, $line;
				if ($line=~/^Total read\#/) {
					$total_read{$sample_name}+=$infor[2];
				}elsif($line=~/^Uniquely paired read\#/){
					$pair_mapped_read{$sample_name}+=$infor[3];
				}
				elsif($line=~/Uniquely mapped read1\#/ or $line=~/Uniquely mapped read2\#/){
					$single_mapped_read{$sample_name}+=$infor[3];
				}
			}
			close TMP;
		}

	}
	closedir DIR;
}
close IN;

foreach my $key (keys %total_read) {
	#print $anno_read{$key};exit;
#	if ($total_read2{$key} ==($pair_mapped_read{$key}/2+$single_mapped_read{$key})) {
#		print "$key\n";
#	}
#	else{
#		print $total_read2{$key}, "\t",($pair_mapped_read{$key}/2+$single_mapped_read{$key}),"\n";
#	}
	print "$key\t$total_read{$key}\t",$pair_mapped_read{$key}+$single_mapped_read{$key},"\t",($pair_mapped_read{$key}+$single_mapped_read{$key})/$total_read{$key},"\t",$anno_read{$key}/($pair_mapped_read{$key}/2+$single_mapped_read{$key}),"\t$total_gene_exp{$key}\n";

}
}
elsif($step eq 'merge'){

open IN,$sample_list;
open O,">merge.sh"||die;
while (<IN>) {
	chomp;
	opendir DIR,$_;
	my $sample_name=(split /\//,$_)[-1];
	my $cmd= "samtools merge $sample_name.sort.bam ";
	while (my $file=readdir(DIR)) {
		if ($file=~/.*_R1_.*fastq.gz$/) {
			my $bam=$file;
			$bam=~s/\.fastq\.gz//;
			$bam=~s/R1_//;
			$bam.=".bam";
			$cmd.=" $file/$bam";
		}
	}
	closedir DIR;
#	print O $bash,"\n";
	print O $cmd,"\n";
	#$cmd.="|samtools rmdup - $sample_name.sort.rmdup.bam\n";
}
close IN;
close O;
#my $tmp_cmd="perl /home/jhmi/xinli/scripts/qsub_sge.pl merge.sh";
#system($tmp_cmd);

}


elsif($step eq 'count'){

open IN,$sample_list;
open O,">count.sh"||die;

#samtools sort -n 41_1.sort.bam 41_1.sort.byname
#samtools view 41_1.sort.byname.bam|htseq-count -q -s reverse  - ~/amber3/no_back_up/data/gene/Homo_sapiens.GRCh37.66.main.gtf > 41_1.sort.byname.bam.count
while (<IN>) {
	chomp;
	my $sample_name=(split /\//,$_)[-1];
	my $cmd= "samtools sort -n $sample_name.sort.bam $sample_name.sort.bam.byname && samtools view $sample_name.sort.bam.byname.bam|htseq-count -q -s reverse  - $gtf > $sample_name.count";
#	print O $bash,"\n";
	print O $cmd,"\n";
}
close IN;
close O;
#my $tmp_cmd="perl /home/jhmi/xinli/scripts/qsub_sge.pl --resource gwas,mf=50G,h_vmem=60G count.sh";
#system($tmp_cmd);

}
elsif($step eq 'count_exon'){

open IN,$sample_list;
open O,">count_exon.sh"||die;

#samtools sort -n 41_1.sort.bam 41_1.sort.byname
#samtools view 41_1.sort.byname.bam|htseq-count -q -s reverse  - ~/amber3/no_back_up/data/gene/Homo_sapiens.GRCh37.66.main.gtf > 41_1.sort.byname.bam.count
while (<IN>) {
	chomp;
	my $sample_name=(split /\//,$_)[-1];
	my $cmd= "samtools view $sample_name.sort.bam.byname.bam|python /home/jhmi/xinli/bin/dexseq_count.py -p yes -s reverse $gtf_exon - $sample_name.exon.count";
#	print O $bash,"\n";
	print O $cmd,"\n";
}
close IN;
close O;
#my $tmp_cmd="perl /home/jhmi/xinli/scripts/qsub_sge.pl --resource gwas,mf=50G,h_vmem=60G count.sh";
#system($tmp_cmd);

}
elsif($step eq 'tdf'){
open IN,$sample_list;
open O,">tdf.sh"||die;

#~/soft/IGVTools/igvtools count -z 5 -w 25 -e 100 5SFM_2.sort.bam 5SFM_2.sort.bam.tdf hg19
while (<IN>) {
	chomp;
	my $sample_name=(split /\//,$_)[-1];
	my $cmd= "/opt/apps/IGVTools/igvtools count -z 5 -w 25 -e 100 $sample_name.sort.bam $sample_name.sort.bam.tdf hg19";
#	print O $bash,"\n";
	print O $cmd,"\n";
}
close IN;
close O;
#my $tmp_cmd="perl /home/jhmi/xinli/scripts/qsub_sge.pl count.sh";
#system($tmp_cmd);

}

sub Usage {
    print << "    Usage";

        $Function

        Usage: $0 <options>

                -sample_list        sample list

                -ref_dir            reference OSA base path

                -ref                reference name
				
                -gene_model         gene_model name

                -gtf                gene gtf file

                -gtf_exon           exon gff file for DEXSeq

                -step               align or merge_rmdup or count or count_exon or tdf

                -h or -help         Show Help , have a choice

    Usage
        exit;

}
