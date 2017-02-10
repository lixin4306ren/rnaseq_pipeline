#!/usr/bin/perl
use strict;
use warnings;

#my $out_put_dir;
my $sample_list=$ARGV[0];
my $out_put_dir=$ARGV[1]||"./";
my $ref_dir=$ARGV[2]||"/home/others/xli/data/genome";
my $ref=$ARGV[3]||"Human.B37.3";
my $gene_model=$ARGV[4]||"RefGene";
my $pwd=`pwd`;chomp $pwd;

open IN,$sample_list;
while (<IN>) {
        chomp;
        opendir DIR,$_;
        my $fastq_files="";
        my $sample_name=(split /\//,$_)[-1];
        while (my $file=readdir(DIR)) {
	#print "$file\n";
                if ($file=~/R1.*fastq.gz$/) {
			my $pair=$file;
                        $pair=~s/R1/R2/;
                        $fastq_files.="$_/$file\n$_/$pair\n";
                }

        }
        closedir DIR;


print "mono /home/others/xli/soft/oshell/16-8-11/oshell.exe --runscript $ref_dir $pwd/$sample_name.oscript ./ /opt/apps/mono/bin/mono\n";

open O,">$pwd/$sample_name.oscript"||die;

chomp($fastq_files);

my $slash="\\\\";
my $back="\b";

print O <<"EOT";
Begin Macro;
 
\@ThreadNum\@ 6;
\@ProjectName\@ $sample_name;
\@ProjectFolder\@ \"$out_put_dir\";
 
\@FileNames\@
"
$fastq_files
";
 
\@CompressionMethod\@ Gzip;
\@Gzip\@ True;
\@PairedEnd\@ True;
\@ReferenceName\@ $ref;
\@GeneModelName\@ $gene_model;
 
End;
EOT



print O <<'EOT2';
#region Oshell run
//Create the OmicSoft project enviroment
Begin NewProject;
File "@ProjectFolder@/@ProjectName@.osprj";
Options /Distributed=true;
End;
 
#region Raw data QC
//QC: Raw sequence, using NgsQCWizard
Begin NgsQCWizard /Namespace=NgsLib;
Files
"@FileNames@";
Options /ThreadNumberPerJob=@ThreadNum@ /FileFormat=AUTO /CompressionMethod=@CompressionMethod@ /QualityEncoding=Automatic /PreviewMode=false 
/BasicStatistics=true /BaseDistribution=true /QualityBoxPlot=true /KMerAnalysis=true /SequenceDuplication=true;
Output rawseq_qc;
End;

 
// save the OmicSoft project enviroment after each major step
Begin SaveProject;
Project @ProjectName@;
File "@ProjectFolder@/@ProjectName@.osprj";
End;
#endregion
 
#region RNA-Seq alignment to human reference using OSA
Begin MapRnaSeqReadsToGenome /Namespace=NgsLib;
Files 
"@FileNames@";
Reference @ReferenceName@;
GeneModel @GeneModelName@;
Trimming  /Mode=TrimByQuality /ReadTrimQuality=2;
Options  /BamSubFolder=primary_alignment /ParallelJobNumber=1 /PairedEnd=@PairedEnd@ /FileFormat=AUTO /AutoPenalty=True 
/FixedPenalty=2 /Greedy=false /IndelPenalty=2 /DetectIndels=False /MaxMiddleInsertionSize=10 /MaxMiddleDeletionSize=10 
/MaxEndInsertionSize=10 /MaxEndDeletionSize=10 /MinDistalEndSize=3 /ExcludeNonUniqueMapping=True /ReportCutoff=10 
/WriteReadsInSeparateFiles=True /OutputFolder="" /GenerateSamFiles=True /ThreadNumberPerJob=@ThreadNum@ 
/InsertSizeStandardDeviation=40 /ExpectedInsertSize=300 /InsertOnSameStrand=False /InsertOnDifferentStrand=True 
/QualityEncoding=Automatic /CompressionMethod=@CompressionMethod@ /Gzip=@Gzip@ /SearchNovelExonJunction=True /ExcludeUnmappedInBam=True;
Output primary_alignment;
End;
 
Begin SaveProject;
Project @ProjectName@;
File "@ProjectFolder@/@ProjectName@.osprj";
End;
 
#endregion
 
#region Alignment QC
//Alignment QC Metrics using RnaSeqQCMetrics
Begin RnaSeqQCMetrics /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
GeneModel @GeneModelName@;
Metrics Alignment,Flag,Profile,Source,InsertSize,Duplication,Coverage,Strand;
Options /ThreadNumberPerJob=@ThreadNum@ /ExcludeFailedAlignments=true /ExcludeSecondaryAlignments=true 
/ExcludeMultiReads=false /ExcludeSingletons=false;
Output alignment_qc;
End;
 
//RNA-Seq 5'->3' trend using SummarizeRnaSeqTrend53
Begin SummarizeRnaSeqTrend53 /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
GeneModel @GeneModelName@;
Options /OutputFolder="@ProjectFolder@/@ProjectName@/primary_alignment" /ThreadNumberPerJob=@ThreadNum@ 
/BinNumber=100 /TranscriptLengthBins=500, 1000, 2000, 3000, 4000, 5000 /ExcludeGenesWithMultipleIsoforms=true 
/ScaleCoverage=true /ReportTranscriptData=false;
Output qc_trend53;
End;
 
Begin SaveProject;
Project @ProjectName@;
File "@ProjectFolder@/@ProjectName@.osprj";
End;
#endregion
 
#region Quantification at gene, isoform, exon and exon-junction level
//gene level rsem estimation
Begin ReportGeneTranscriptCounts /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
GeneModel @GeneModelName@;
Options /ExpressionMeasurement=RPKM /ThreadNumberPerJob=@ThreadNum@ /Add1=False /CountFragments=True /ExcludeMultiReads=False 
/UseEffectiveTranscriptLength=True /CountStrandedReads=True /CountReverseStrandedReads=True;
Output quantification_gene;
End;
 
//transcript level rsem estimation
Begin ReportGeneTranscriptCounts /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
GeneModel @GeneModelName@;
Options /ExpressionMeasurement=RPKM_Transcript /ThreadNumberPerJob=@ThreadNum@ /Add1=False /CountFragments=True /ExcludeMultiReads=False 
/UseEffectiveTranscriptLength=True /CountStrandedReads=True /CountReverseStrandedReads=True;
Output quantification_tx;
End;
 
//quantify the known exon and exon junctions based on a gene model
Begin ReportExonCounts /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
GeneModel @GeneModelName@;
Options /ThreadNumberPerJob=@ThreadNum@ /ReportExonCounts=True /ReportExonJunctionCounts=True /ExcludeSingletons=False /Rpkm=True 
/CountStrandedReads=True /CountReverseStrandedReads=True;
Output quantification_known;
End;
 
//quantify the each exon junction based on alignment and annotate with gene model info (known or novel exon junctions)
Begin ReportExonJunctions /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
Options /ExcludeMultiReads=False /ExcludeSingletons=False /NMCutoff=4 /JunctionOverhangCutoff=4 /GenerateReport=True 
/GenerateBedFile=False /BedFileOutputFolder="" /MinimalHit=1 /ThreadNumberPerJob=@ThreadNum@ 
/OutputFolder="@ProjectFolder@/@ProjectName@/primary_alignment";
Output quantification_jxn;
End;
 
Begin SaveProject;
Project @ProjectName@;
File "@ProjectFolder@/@ProjectName@.osprj";
End;
#endregion
 
#region Fusion gene detection
//Detect fusion based on inter-transcript read pairs (PE)
Begin ReportPairedEndFusionGenes /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
GeneModel @GeneModelName@;
Options /GenerateData=True /OutputFusionReads=True /MinimalHit=2 /FusionReportCutoff=1 /ThreadNumberPerJob=@ThreadNum@ 
/FilterBy=DefaultList /DefaultFilterListVersion=v1 /FilterGeneListFileName="" /FilterGeneFamilyFileName="" 
/GenerateTableland=True /OutputFolder="@ProjectFolder@/@ProjectName@/fusion_PEReads";
Output FusionPE;
End;
 
//Detect fusion based alignment of junction spanning reads (SE)
Begin MapFusionReads /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
Reference @ReferenceName@;
GeneModel @GeneModelName@;
Trimming /Mode=TrimByQuality /ReadTrimQuality=2;
Options /FusionVersion=2 /ParallelJobNumber=1 /PairedEnd=False /RnaMode=True /FileFormat=BAM /AutoPenalty=True /FixedPenalty=2 
/OutputFolder="@ProjectFolder@/@ProjectName@/fusion_SE_alignment" /MaxMiddleInsertionSize=10 
/ThreadNumberPerJob=@ThreadNum@ /QualityEncoding=Automatic /CompressionMethod=None /Gzip=False /MinimalFusionAlignmentLength=25 
/FilterUnlikelyFusionReads=True /FullLengthPenaltyProportion=8 /OutputFusionReads=True /MinimalHit=4 /MinimalFusionSpan=5000 
/FusionReportCutoff=1 /NonCanonicalSpliceJunctionPenalty=2 /FilterBy=DefaultList /DefaultFilterListVersion=v1 
/FilterGeneListFileName="" /FilterGeneFamilyFileName="" /GenerateTableland=True;
Output FusionSE;
End;
#endregion
 
#region Mutation detection
//detect mutation and annotate with dbSNP; will download dbSNP from omicsoft at first time 
Begin SummarizeMutation /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\primary_alignment;
Options /BaseQualityCutoff=20 /MapQualityCutoff=20 /MinimalIndelSize=1 /ExcludeSingletons=False /ExcludeMultiReads=False 
/LeftExclusion=5 /RightExclusion=5 /ThreadNumberPerJob=@ThreadNum@ /MinimalTotalHit=20 /MinimalMutationHit=5 /MinimalMutationFrequency=0.15 
/ExcludeNonMutantSites=True /GenerateSummarizedReport=True /GenerateIndividualReport=False /GenerateTableland=True 
/MaxFrequencyCutoff=0.15 /DbsnpVersion=v135 /OutputFolder="@ProjectFolder@/@ProjectName@/mutation";
Output MutationSummary;
End;
 
//Annotate Mutation 
Begin AnnotateMutation /Namespace=NgsLib;
Project @ProjectName@;
Data @ProjectName@\\MutationSummary.MutationReport;
ID ID;
Chromosome Chromosome;
Position Position;
Mutation Mutation;
Other (default);
Options /ReferenceLibraryID=@ReferenceName@ /GeneModelID=@GeneModelName@ /DbsnpVersion=v135 /GenerateClusteringFlag=True 
/ClusteringFlagWindowSize=10 /GenerateTableland=True /AnnotateLongestTranscriptOnly=False /AnnotateFunctionalMutation=False 
/AnnotateSomaticMutation=False;
Output MutationAnnotated.MutationReport;
End;
 
Begin SaveProject;
Project @ProjectName@;
File "@ProjectFolder@/@ProjectName@.osprj";
End;
 
#endregion
 
Begin ExportView;
Project @ProjectName@;
OutputFolder "@ProjectFolder@/@ProjectName@/ExportedResults";
End;
 
Begin CloseProject;
Project @ProjectName@;
End;

EOT2



close O;

}
close IN;


