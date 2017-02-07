#!/usr/bin/perl -w 

use strict; 
use warnings;
use Devel::Size qw(size total_size);
use List::Util 'max';
use List::Util qw(sum);

do "/home/hansenlo/Code/usefulFunctions.pl";
die $@ if $@;


my $expFastq=shift; #the file containing the fastq reads for the experiment
my $expKmerCounts=shift; #the file containing the counts of how often each kmer occurs
my $controlFile=shift; #the file containing the control reads if present


#given a contig file and a cutoff indicating how much to start aligning from each end
#iteratively increase the size of the end alignments until you get a unique alignment or you reach
#a preset size limit
sub iterativeAlignmentLocal($ $ $ $ $ $);


#given a set of contig edges that align uniquely 
#extend the contig edges until they reach a cutoff number of mismatches with the reference
#arguments are contig edges cutoff reference genome and contig seq
sub extendAlignmentsLocal($ $ $ $);

#given a contig alingnment file and a contig fasta file select those contigs that align with substantial soft or hard clipping or do not align at all
#output is a fasta file also given is the cutoff number of bases that need to be soft or hard clipped.  
sub filterAlignmentsLocal($ $ $);


#Program to take as input a file name containing a read library and running all commands on it

my($cmd, $line, $i, $header, $outputUnique, $outputOverlapEdena, $contigEdgesOutput, $alignmentOutput, $tempIn, $tempOut);
my($allReads, @splitLine, $path, $outputKmerLib, $dataset, $blastDB, $tempOutunMappPoorMap);

@splitLine=split(/\//, $expFastq);

print "$splitLine[$#splitLine]\n";

$splitLine[$#splitLine]=~m/(.*).fastq$/;

$header="unique_".$1;

$dataset=$1;

#print "string is $1\n";
#exit();


$expFastq=~m/(.*)\/(.*).fastq$/;

$path=$1."/";



#$blastDB="/data6/Bowtie2Index/elegansCe10";
#$blastDB="/data/Bowtie2Index/elegansCe10";

#$blastDB="/data/Bowtie2Index/hg19";


#$blastDB="/data/Bowtie2Index/hg19GATK";

$blastDB="/data/BwaIndex/allChrhg19InOrder.fa";
#$blastDB="/data/BwaIndex/allChr_cElegans10.fa";
#$blastDB="/data6/sukrit/081216_MiSeq_MMB1newdel_genomeSeq/MappingToReference/MMB1genomeCIRC84.fasta";


if(1==1)
{
#compiling the code
#$cmd="g++ findDiff_BitString_HashTableVer.cpp spooky.cpp -g -O3 -o findDiff_BitString_HashTableVer -std=c++11 -Wno-deprecated -fopenmp";

    $cmd="make";

#$cmd="g++ findDiff_BitString_HashTableVer_FilterOnControl.cpp spooky.cpp -g -O3 -o temp -std=c++0x -Wno-deprecated";


#$cmd="g++ findDiff_BitString.cpp spooky.cpp -g -O3 -o findDiff_BitString -std=c++0x -Wno-deprecated";


print "#####cmd is $cmd\n";
system($cmd)==0
    or die "system $cmd failed\n";

}
$outputUnique=$header.".fastq";



if(1==1)
{

#getting those reads that passed the filter running the code
#$tempIn="/data/HEM0013-4/".$file;

#$tempIn=$file;

#$outputKmerCount="/data4/HEM0013_158/uniqueKmerCount_unique_HEM0013-158-MC_allReads.dat";

#$cmd="time ./findDiff_HashTable_filterReads $tempIn 20 0 $outputKmerCount > "."/home/hansenlo/SeqDiff/Results/".$outputUnique;

#$cmd="time ./findDiff_BitString_HashTableVer $controlFile $file  31 0 $cutoff temp.dat > "."/home/hansenlo/SeqDiff/Results/".$outputUnique;

$cmd="time ./variantFinder $expKmerCounts $expFastq 45 > temp.dat";

#$cmd="time ./temp $controlFile $file  31 0 $cutoff temp.dat > "."/data7/SeqDiffResults/Results/".$outputUnique;;


#$cmd="time ./findDiff_BitString $controlFile $file 20  0 $cutoff temp.dat > "."/data5/SeqDiffResults/Results/".$outputUnique;;


#$cmd="time ./findDiff_BitString $controlFile $file  18 0 $cutoff temp.dat > "."/data4/".$outputUnique;

#deleting old debugging file if present
#system("rm debuggingMatrix.dat");
    #or die "system rm debuggingMatrix.dat failed\n";



print "#####cmd is $cmd\n";
system($cmd)==0
    or die "system $cmd failed\n";



}


#copying sequence to appropriate places
if(1==1)
{
#move contigs to appropriate directory
$cmd="cp /home/hansenlo/SeqDiff/gitHubProject/SeqDiff/contigs.fa /data5/SeqDiffResults/Results/".$header."_contigs.fasta";
print "####cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";



#move unique reads to the appropriat directory 
$cmd="cp /data/uniqueReads.fastq /data5/SeqDiffResults/Results/".$header."_unique.fastq";
print "####cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";
}



#finding the overlap between unique reads
$outputOverlapEdena="/data5/SeqDiffResults/Results/".$header.".ovl";

if(1==2)
{

$cmd="/home/hansenlo/bin/edena -nThreads 14 -minOverlap 25 -singleEnd /data5/SeqDiffResults/Results/".$outputUnique." -p $header";

print "#####cmd is $cmd\n";
system($cmd)==0
    or die "system $cmd failed\n";


#moving output file to correct dirctory
$cmd="mv ".$header.".ovl"." /data5/SeqDiffResults/Results/";
system($cmd)==0
    or die "system $cmd failed\n";
}

#using velvet
if(1==2)
{
   

    $cmd="velveth /home/hansenlo/SeqDiff/Temp/ 31 -fastq -short /data5/SeqDiffResults/Results/".$outputUnique;
    
#running Andys reads through the assembler
 #$cmd="velveth /home/hansenlo/SeqDiff/Temp/ 31 -fasta -short /home/hansenlo/AndysStuff/mergedAndysFiles.fasta";

    print "####cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";

    #open(OUT, ">temp.dat");
    #for($i=1; $i<50; $i++)
    #{
	#$cutoff=$i;

#	$cmd="velvetg /home/hansenlo/SeqDiff/Temp/ -read_trkg yes -cov_cutoff $cutoff -min_contig_lgth 101";
 

	print "####cmd is $cmd\n\n";
	system($cmd)==0
	    or die "system $cmd failed\n";


	#listing how many contigs there are
	$cmd="grep \">\" Temp/contigs.fa | wc -l";
	#print "#######cmd is $cmd\n\n";

	print "######Number of contigs is ###########\n\n\n\n\n\n";

	my $numContigs=`$cmd`;

	chomp($numContigs);
	print "numContigs is $numContigs\n";

    print "\n\n\n\n\n\n\n\n";



    $cmd="cp /home/hansenlo/SeqDiff/Temp/contigs.fa /data5/SeqDiffResults/Results/".$header."_contigs.fasta";
    print "####cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";

	#print OUT "$cutoff\t$numContigs\n";

    #}
    #close(OUT);

}


if(1==2)
{

#####assembling reads into contigs
#using ednea
$cmd="/home/hansenlo/bin/edena -c 90 -minCoverage 2 -e ".$outputOverlapEdena." -p $header";

#using velvet
#$cmd="velvetg /home/hansenlo/SeqDiff/Temp -cov_cutoff 5";
#$cmd="velvetg /home/hansenlo/SeqDiff/Temp";



$cmd="tr ‘[:lower:]‘ ‘[:upper:]‘ <  /data5/SeqDiffResults/Results/".$outputUnique." > /data5/SeqDiffResults/Results/temp.fastq";


print "####cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

$cmd="fastq_to_fasta -n -i  /data5/SeqDiffResults/Results/temp.fastq > /data5/SeqDiffResults/Results/foo.dat";

print "####cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

$cmd="time /home/hansenlo/bin/CAP3/cap3 /data5/SeqDiffResults/Results/foo.dat -p 99  > temp.dat";
print "####cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


if(1==2)
{
$cmd="mv *contig* /data5/SeqDiffResults/Results/";
system($cmd)==0
    or die "system $cmd failed\n";

$cmd="mv *log* /data5/SeqDiffResults/Results/";
system($cmd)==0
    or die "system $cmd failed\n";
}
}



$contigEdgesOutput="/data5/SeqDiffResults/Results/".$header."_contigEdges.fasta";

#$contigEdgesOutput="/data5/SeqDiffResults/Results/".$header."RefGenomeOnly_contigEdges.fasta";
if(1==2)
{
#####getting the edges of contigs
#$contigEdgesOutput="/data5/SeqDiffResults/Results/".$header."_contigEdges.fasta";

#$cmd="time perl VariantPerlScripts.pl "."/data5/SeqDiffResults/Results/".$header."_contigs.fasta"." 30 > $contigEdgesOutput";

#CAP3 command
#$cmd="time perl VariantPerlScripts.pl "."/data5/SeqDiffResults/Results/foo.dat.cap.contigs"." 30 > $contigEdgesOutput";

#velvet command temp.dat is a dummy file extract flag means use the extract fasta function
$cmd="time perl VariantPerlScripts.pl "."/home/hansenlo/SeqDiff/gitHubProject/SeqDiff/contigs.fa"." 40 temp.dat extract  > $contigEdgesOutput";

print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#copy contigs to Results folder
#$cmd="cp /home/hansenlo/SeqDiff/Temp/contigs.fa "."/data5/SeqDiffResults/Results/".$header."_contigs.fasta";
#print "#######cmd is $cmd\n\n";
#system($cmd)==0
#    or die "system $cmd failed\n";

}


if(1==1)
{


#aligning all unique reads to artifical chromosome
#$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header.".sam";



##############################aligning all unique reads##################
#aligning all unique reads 
$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header.".sam";

#$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."_refGenomeOnly.sam";

#$blastDB="/home/hansenlo/Genomes/Bowtie2Index/hg19";


#$blastDB="/data/BowtieIndex/elegansCe10";


#$blastDB="/home/gabdank/genomes_fasta/elegans_W235_index";

#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB "."/data5/SeqDiffResults/Results/"."$outputUnique --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_uniqueReads.fastq "."-S $alignmentOutput";

#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB "."/data/uniqueReads.fastq --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_uniqueReads.fastq "."-S $alignmentOutput";

$cmd="time /data/bin/bwa-master/bwa mem -t 20 $blastDB"." /data/uniqueReads.fastq > $alignmentOutput";


#fasta command DONT FORGET
#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB -f "."/data5/SeqDiffResults/Results/"."$outputUnique --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_uniqueReads.fastq "."-S $alignmentOutput";


#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 4 -x /home/hansenlo/Genomes/Bowtie2Index/seqA "."/data5/SeqDiffResults/Results/"."$outputUnique --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_uniqueReads.fastq "."-S $alignmentOutput";


print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".bam";

#converting alignments to bam format
$cmd="samtools view -bS $alignmentOutput > ".$tempOut;
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#sorting the aligments
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".sorted";

#sorting bam files
$cmd="samtools sort -m 40000000000 $tempIn $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating an index
$tempIn=$tempOut.".bam";

$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating a unique alignment
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".sorted.unique.bam";
$cmd="samtools view -bq 15 $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#index the unique alignment 
$tempIn=$tempOut;
$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#convert to bed file
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".sorted.unique.bed";
$cmd="bamToBed -i $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";
}



######################Aligning the contigs#####################
#$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."refGenomeOnly_contigEdges.sam";

$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."_contigs.sam";

my $contigs="/data5/SeqDiffResults/Results/".$header."_contigs.fasta";

if(1==1)
{

#$cmd="rm /data5/SeqDiffResults/Results/Alignment/".$header."*contigs*";  ##############Be careful with THIS!!!!! make sure you are not deleting the wrong files!!!!

print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#$blastDB="/home/hansenlo/Genomes/Bowtie2Index/hg19";

#$blastDB="/data/BowtieIndex/elegansCe10";

#$blastDB="/home/gabdank/genomes_fasta/elegans_W235_index";

#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2   -p 20 -x $blastDB -f $contigs --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta "."-S $alignmentOutput";

$cmd="time /data/bin/bwa-master/bwa mem -t 20 $blastDB"." $contigs > $alignmentOutput";


#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 --mp 10 --rdg 3,2 --rfg 3,2  -p 20 -x $blastDB -f $contigs --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta "."-S $alignmentOutput";



#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 4 -x /home/hansenlo/Genomes/Bowtie2Index/seqA -f $contigEdgesOutput --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_ContigEdges.fasta "."-S $alignmentOutput";


print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#getting bam files for contig edges alignment
my $contigheader=$header."_contigs";
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader.".bam";

#converting alignments to bam format
$cmd="samtools view -bS $alignmentOutput > ".$tempOut;
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#sorting the aligments
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader.".sorted";

#sorting bam files
$cmd="samtools sort $tempIn $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating and index
$tempIn=$tempOut.".bam";

$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


##Getting the unmapped contigs
my $tempOutUnMapped="/data5/SeqDiffResults/Results/Alignment/".$contigheader."_unmapped.sam";

$cmd="time samtools view -f 4 $tempIn > $tempOutUnMapped";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


my $tempOutPoorMap="/data5/SeqDiffResults/Results/Alignment/".$contigheader."_poorlyMapped.sam";
#$cmd="time samtools view -h $tempIn".'| awk \'{if(\$1 ~ /^@/) {print \$0} else if(\$5<8) {print \$0}}\' |'."samtools view -hSo $tempOutPoorMap -";

$cmd="time samtools view -h $tempIn".'| awk \'{if($5<8) {print $0}}\' |'."samtools view -hSo $tempOutPoorMap -";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#Converting to fasta format and concentanating files
#$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader."_temp.sam"
$cmd="time samtools bam2fq $tempOutUnMapped > temp.fa";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

$cmd="time samtools bam2fq $tempOutPoorMap > foo.fa";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#concatenating files
$tempOutunMappPoorMap="/data5/SeqDiffResults/Results/Alignment/".$contigheader."_unMapped_poorlyMapped.fasta";
$cmd="time cat temp.fa foo.fa > $tempOutunMappPoorMap";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#removing intermediate files
$cmd="rm foo.fa temp.fa";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";





#creating a unique alignment
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader.".sorted.unique.bam";
$cmd="samtools view -bq 20 $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating and index
$tempIn=$tempOut;

$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#convert to bed file
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader.".sorted.unique.bed";
$cmd="bamToBed -i $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

}


#######################now working with contig edges###################

    my $alignmentFile=shift;
    my $cutoff=shift;
    my $contigFile=shift;

print "starting to filter the alignments \n";

#&filterAlignments($alignmentOutput, 0.2, $contigs);


print "Starting iterative alignment\n";
#iterativeAlignment($tempOutunMappPoorMap, 20, 80,  $blastDB, $header);

#&iterativeAlignment("/data5/SeqDiffResults/Results/Alignment/unique_platinumChr21_plusUnmapped_contigs_unMapped_poorlyMapped.fasta", 20, 80,  $blastDB, $header);

#&iterativeAlignment("/data5/SeqDiffResults/Results/unique_platinumChr21_plusUnmapped_contigs.fasta", 20, 80,  $blastDB, $header);

#&iterativeAlignment("/data5/SeqDiffResults/Results/unique_allPlatinum_contigs.fasta", 20, 80,  $blastDB, $header);

my $mappingQualityCutoff=35;
#&iterativeAlignment("clippedContigs.fa", 20, 60,  $blastDB, $header, $mappingQualityCutoff);

#&iterativeAlignment($contigs, 20, 180,  $blastDB, $header, $mappingQualityCutoff);


print "Starting to extend alignments\n";

#&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/human19/allChrhg19InOrder.fa", "test.dat");


#&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/human19/allChrhg19InOrder.fa", "/data5/SeqDiffResults/Results/Alignment/unique_platinumChr21_plusUnmapped_contigs_unMapped_poorlyMapped.fasta");

#&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/human19/allChrhg19InOrder.fa", "/data5/SeqDiffResults/Results/unique_platinumChr21_plusUnmapped_contigs.fasta");

#&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/human19/allChrhg19InOrder.fa", "/data5/SeqDiffResults/Results/unique_allPlatinum_contigs.fasta");

#&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/human19/allChrhg19InOrder.fa", $contigs);

#&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/cElegans10/allChr.fa", $contigs);

#&extendAlignments("mappedReads_IDs.dat", 2, "/data6/sukrit/081216_MiSeq_MMB1newdel_genomeSeq/MappingToReference/MMB1genomeCIRC84.fasta", $contigs);




if(1==2)
{


    #&iterativeAlignment("/data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta", 20, 80,  $blastDB, $header);

    #&extendAlignments("mappedReads_IDs.dat", 2, "/data6/Genomes/cElegans10/allChr.fa", "/data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta");

    #&extendAlignments("mappedReads_IDs.dat", 2, "/data/Genomes/human19/allChrhg19InOrder.fa", "/data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta");


    #my $edges=shift;
    #my $cutoff=shift;
    #my $ref=shift;
    #my $contig=shift;


    #my $contigs=shift;
    #my $startSize=shift;
    #my $endSize=shift;
    #my $blastDB=shift;
    #my $header=shift;



#aligning the contig edges
#$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."refGenomeOnly_contigEdges.sam";

$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."_contigEdges.sam";

#$blastDB="/home/hansenlo/Genomes/Bowtie2Index/hg19";

#$blastDB="/data/BowtieIndex/elegansCe10";

#$blastDB="/home/gabdank/genomes_fasta/elegans_W235_index";

$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 --local -x $blastDB -f $contigEdgesOutput --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_ContigEdges.fasta "."-S $alignmentOutput";

$cmd="time /data/bin/bwa-master/bwa mem -t 20 $blastDB"." $contigEdgesOutput > $alignmentOutput";


#$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes/hg19 --outSAMunmapped Within --readFilesIn $contigEdgesOutput --runThreadN 20 --outFileNamePrefix Test/test --outFilterMismatchNmax 6"; 

#my $StarOutput="/home/hansenlo/SeqDiff/Test/testAligned.out.sam";


#system("cp $StarOutput $alignmentOutput"); 



#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 4 -x /home/hansenlo/Genomes/Bowtie2Index/seqA -f $contigEdgesOutput --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_ContigEdges.fasta "."-S $alignmentOutput";


print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#getting bam files for contig edges alignment
my $edgeheader=$header."_ContigEdges";
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$edgeheader.".bam";

#converting alignments to bam format
$cmd="samtools view -bS $alignmentOutput > ".$tempOut;
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#sorting the aligments
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$edgeheader.".sorted";

#sorting bam files
$cmd="samtools sort $tempIn $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating and index
$tempIn=$tempOut.".bam";


#convert to bed file all Alignments
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$edgeheader.".sorted.bed";
$cmd=" /data/bin/Bedtools/bedtools2/bin/bamToBed -i $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating a unique alignment
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$edgeheader.".sorted.unique.bam";
$cmd="samtools view -bq 40 $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#convert to bed file
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$edgeheader.".sorted.unique.bed";
$cmd=" /data/bin/Bedtools/bedtools2/bin/bamToBed -i $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";
}








if(1==2)
{

    $controlFile=$expFastq;

@splitLine=split(/\//,$controlFile);

print "$splitLine[$#splitLine]\n";

$splitLine[$#splitLine]=~m/(.*).fastq$/;

##############################aligning all Control reads##################
#$controlFile=~m/(.*).fastq/;
$header=$1;

#$header="allmut_debug";

#aligning all control reads 
$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header.".sam";
#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 4 -x /home/hansenlo/Genomes/Bowtie2Index/hg19 "."/data/HEM0013-4/"."$controlFile --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable.fastq "."-S $alignmentOutput";

#$blastDB="/home/hansenlo/Genomes/Bowtie2Index/hg19";

#$blastDB="/data/BowtieIndex/elegansCe10";


$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB "."$controlFile --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable.fastq "."-S $alignmentOutput";

$cmd="time /data/bin/bwa-master/bwa mem -t 20 $blastDB"." $controlFile > $alignmentOutput";



#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 4 -x /home/hansenlo/Genomes/Bowtie2Index/seqA "."-f /data5/SeqDiffResults/Results/"."foo.dat.cap.contigs --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable.fastq "."-S $alignmentOutput";

print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".bam";

#converting alignments to bam format
$cmd="samtools view -bS $alignmentOutput > ".$tempOut;
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#sorting the aligments
$tempIn=$tempOut;
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".sorted";

#sorting bam files
$cmd="samtools sort $tempIn $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating an index
$tempIn=$tempOut.".bam";

$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#creating a unique alignment
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$header.".sorted.unique.bam";
$cmd="samtools view -bq 15 $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#index the unique alignment 
$tempIn=$tempOut;
$cmd="samtools index $tempIn";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";
}



sub iterativeAlignmentLocal($ $ $ $ $ $)
{
    my $contigs=shift;
    my $startSize=shift;
    my $endSize=shift;
    my $blastDB=shift;
    my $header=shift;
    my $mappingCutoff=shift;

    my(@samOutput, $flag, $size, $i, @readsToReAlign, %contigs, $contigID, $locationContigs, $j);
    my($seq, $ctr, $start, $end, $strand, %multi, $key, $k, $l, $rowSize);


    $ctr=0;
    
    system("rm mappedReads_IDs.dat");
    open(wasMapped, ">mappedReads_IDs.dat");


    #system("cp $contigs needToBeMapped.fasta");
    #open(needMapped, ">needToBeMapped.fasta");
    

    open(willNotMap, ">finalMultiMapping.bed");

    open(mappedSam, ">mappedEnds.sam");


    %contigs=&readFasta($contigs);


    #$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/temp_contigEdges.sam";


    $alignmentOutput="/home/hansenlo/SeqDiff/Test/testAligned.out.sam";


    #$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."_contigEdges.sam";


    #$locationContigs="/data5/SeqDiffResults/Results/".$header."_contigs.fasta";


    	#getting contig edges
    $cmd="time perl VariantPerlScripts.pl $contigs $startSize temp.dat extract  > needToBeMapped.fasta";

    #$cmd="time perl VariantPerlScripts.pl $contigs 100 temp.dat extract  > needToBeMapped.fasta";

    


    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";
		

    #exit;

    my ($needsMappedCtr);

    $needsMappedCtr=10;
    for($j=$startSize+1; $j<=$endSize; $j++)
    {


	#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB -f needToBeMapped.fasta "." -S $alignmentOutput";

	#only unique alignments
	#$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes --readFilesIn needToBeMapped.fasta --runThreadN 20 --outFileNamePrefix Test/test --outFilterMultimapNmax 1 --outFilterMismatchNmax 6"; 

	#$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes --outSAMunmapped Within --readFilesIn needToBeMapped.fasta --runThreadN 20 --outFileNamePrefix Test/test --outFilterMismatchNmax 6"; 

	#$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes/hg19 --outSAMunmapped Within --readFilesIn needToBeMapped.fasta --runThreadN 20 --outFileNamePrefix Test/test --outFilterMismatchNmax 6"; 


#	$cmd="time /data/bin/bwa-master/bwa mem -t 20 $blastDB"." needToBeMapped.fasta > $alignmentOutput"; #this command will not work with reads shorter than about 70 bps

	#if entire set of contigs map then exit out of loop 
	if($needsMappedCtr==0)
	{
	    last;
	}


#bwa commands to run alignment on short reads and covert to sam format
	$cmd="time /data/bin/bwa-master/bwa aln -t 20 $blastDB needToBeMapped.fasta > reads.sai";
	print "#######cmd is $cmd\n\n";
	system($cmd)==0
	    or die "system $cmd failed\n";

	$cmd="time /data/bin/bwa-master/bwa samse $blastDB reads.sai needToBeMapped.fasta > $alignmentOutput";
	print "#######cmd is $cmd\n\n";
	system($cmd)==0
	    or die "system $cmd failed\n";

	
	#exit;
	
	open(my $needMapped, ">needToBeMapped.fasta");


	@samOutput=&parseSamFile($alignmentOutput);

	$size=@samOutput;


	#debugging code
	if(1==2)
	{
	    for($i=0; $i<$size; $i++)
	    {
		my $foo=@{$samOutput[$i]};
	  
		print "$foo\t";
		for($k=0; $k<11; $k++)
		{
		    print "$samOutput[$i][$k]\t";

		}

		print "\n";
	    }
    
	    exit;
	}


	$needsMappedCtr=0;

	for($i=0; $i<$size; $i++)
	{

	#skip header lines
	    if($samOutput[$i][0]=~m/^\@/)
	    {
		#print out the header for the sam file
		if($j==($startSize+1))
		{
	
		    $rowSize=@{$samOutput[$i]};
		    for($l=0; $l<$rowSize; $l++)
		    {
			if($l<($rowSize-1))
			{
			    print mappedSam "$samOutput[$i][$l]\t";
			}else
			{
			    print mappedSam "$samOutput[$i][$l]";
			}
		    }

		    print mappedSam "\n";
		}
		
		next;
	    }
    

    #check flag to see if read is mapped
    #3rd bit indicates if read is mapped 
    #if bit is set read is not mapped
    #integer 4 represents a number with 3rd bit set
	    if(($samOutput[$i][1] & 4)==4)
	    {
		if($j==$startSize+1)
		{
		    $ctr++;
		}

		if($j>$startSize+1)
		{
		    print willNotMap "chrI\t1000\t2000\t$samOutput[$i][0]\t0\t+\n";
		}

	    #read is not mapped so don't put an entry in the mapped reads bed file
		next;	
	    }
    

	    



    #checking mapping quality to see if unique
	    #also checking to make sure no soft clipping by looking for S in Cigar string
#bowtie mapping quality	    
#if($samOutput[$i][4]>15)
	    if($samOutput[$i][4]>35 && $samOutput[$i][5]!~/S/)
	    {
	    


		$start=$samOutput[$i][3]-1; #start of contig edge
		$end=$samOutput[$i][3]+$j-1; #new end of contig edge

		if($samOutput[$i][1]==0)
		{
		    $strand="+";

		}elsif($samOutput[$i][1]==16)
		{
		    $strand="-";

		}

		


		if($samOutput[$i][0]=~m/^142231\_/)
		{
		    my $temp=1;

		}


		$rowSize=@{$samOutput[$i]};
		#print out the alignment record for the contig edge that was mapped
		for($l=0; $l<$rowSize; $l++)
		{
		    if($l<($rowSize-1))
		    {
			print mappedSam "$samOutput[$i][$l]\t";
		    }else
		    {
			print mappedSam "$samOutput[$i][$l]";
		    }
		}

		print mappedSam "\n";

		
		
		$end=$end-1; #bed format is zero based but same format is 1 based
		print wasMapped "$samOutput[$i][2]\t$start\t$end\t$samOutput[$i][0]\t$samOutput[$i][4]\t$strand\n";
	    }elsif(($samOutput[$i][1] & 256)!=256)
	    {
		
		$samOutput[$i][0]=~m/(.*)\_\d+\_(firstPart|endPart)/;

		
		$contigID=$1;

		print $needMapped ">$samOutput[$i][0]\n";
		
		$needsMappedCtr=$needsMappedCtr+1;

		#extract first part of contig copy to needs mapped file
		if($samOutput[$i][0]=~m/firstPart/)
		{
		    $seq=substr($contigs{$contigID}, 0, $j); 
		    print $needMapped "$seq\n";
       
	  
		}else
		{
		    $seq=substr($contigs{$contigID}, -1*$j); 
		    print $needMapped "$seq\n";
		}

	    }


	}

	close $needMapped;
    
	#last;
    }

    close willNotMap;
    close wasMapped;

    #adding final set of non unique mappings
    open(noMap, ">>finalMultiMapping.bed");

    
    %multi=readFasta("needToBeMapped.fasta");

    foreach $key (keys %multi)
    {
	
	print noMap "chrI\t1000\t2000\t$key\t0\t+\n";
    }


    print "will not map is $ctr\n";


    close mappedSam;


#converting alignments to bam format
    $cmd="samtools view -bS mappedEnds.sam > mappedEnds.bam";
    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";


#sorting the aligments
#sorting bam files
    $cmd="samtools sort -m 40000000000 mappedEnds.bam mappedEnds.sorted";
    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";

#creating an index
    $cmd="samtools index mappedEnds.sorted.bam";
    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";

#moving bam file to correct directory
    $cmd="mv mappedEnds.sorted.bam\*  /data5/SeqDiffResults/Results/Alignment";
    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";

    #removing intermediate files
    $cmd="rm mappedEnds.bam";
    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";

}


sub extendAlignmentsLocal($ $ $ $)
{
    my $edges=shift;
    my $cutoff=shift;
    my $ref=shift;
    my $contig=shift;


    my (%refGenome, @contigEdges, $size, $i, $length, $contigID, $refString, @refSeq, $contigString, %contigSeq, @splitContig);
    my ($sizeEdge, $mismatches, $posCtr, $index, $flag, $sizeContig, $ctr, $newEnd, $j, $newStart);
    
    %refGenome=&readFasta($ref);
    
    @contigEdges=&parseBedfile($edges);

    %contigSeq=&readFasta($contig);


    $size=@contigEdges;


    open(my $extended, ">extendedContigEdges.bed");


    for($i=0; $i<$size; $i++)
    {

	#	$contigID=~m/^(\d+)\_/;


	if($contigEdges[$i][0] eq "chrM")
	{
	    next;
	}


	$contigEdges[$i][3]=~m/(.*)\_(\d+)\_(firstPart|endPart)/;

	$contigID=$1;
	#$length=$2;
	$mismatches=0;
	$posCtr=0;

    
	if($contigEdges[$i][5] eq "+")
	{
		    #$refString=substr($refGenome{$contigEdges[$i][0]}, $contigEdges[$i][1], $length);	
	    $contigString=$contigSeq{$contigID};
	}

	if($contigEdges[$i][5] eq "-")
	{
	    $contigString=&reverse_complement_IUPAC($contigSeq{$contigID});
	}


			
	$length=length($contigString);

       

	if(($contigEdges[$i][3]=~m/firstPart/ && $contigEdges[$i][5] eq "+") || ($contigEdges[$i][3]=~m/endPart/ && $contigEdges[$i][5] eq "-") )
	{	


	    $refString=uc(substr($refGenome{$contigEdges[$i][0]}, $contigEdges[$i][1], $length));	


#	    if(!defined($refGenome{$contigEdges[$i][0]}))
#	    {
#		my $temp=0;
#	    }


	    @refSeq=split(//, $refString);


	    
	    #my $temp=substr($refGenome{$contigEdges[$i][0]}, ($contigEdges[$i][2]-1), $length);	

	    $sizeEdge=$contigEdges[$i][2]-$contigEdges[$i][1];
	    
#	    print "\n\n\n";
#	    print "$refString\n";
#	    print "$contigString\n";
	    

	    #my $debug= (" " x ($sizeEdge-2)).$temp;

#	    print "$debug\n";

#	    print "\n\n\n";
	    
	    @splitContig=split(//, $contigString);

	   #  116156_41
	   

	    #if($contigEdges[$i][3]=~m/133774\_/)
	    #{
	#	my $foo=1;

	 #   }

	    $index=$sizeEdge-1;
	    $ctr=0;
	    $mismatches=0;
	    $posCtr=1;

	    $sizeContig=length($contigString);

	    my $ctrAdvanced=0;

	    for($j=$index; $j<$sizeContig; $j++)
	    {
	       		
		if($splitContig[$j] ne $refSeq[$j])
		{
		    $mismatches+=1;   		   
		}


		if($contigEdges[$i][3]=~m/116167\_/)
		{
		    my $foo=1;

		}

		
		

		$ctr++;

		if($mismatches==2 && ($posCtr==2))
		{
		    #$newEnd=$contigEdges[$i][1]+$j-2;
 		    $newEnd=$contigEdges[$i][2]+$ctrAdvanced-2;
		    print $extended "$contigEdges[$i][0]\t$contigEdges[$i][1]\t$newEnd\t$contigEdges[$i][3]\t$contigEdges[$i][4]\t$contigEdges[$i][5]\n";
		
		    $posCtr=0;
		    $mismatches=0;
		    last;
		}


		if($mismatches>=2 && ($posCtr==3))
		{
		    $newEnd=$contigEdges[$i][2]+$ctrAdvanced-3;
		    #$newEnd=$contigEdges[$i][1]+$j-3;
		    print $extended "$contigEdges[$i][0]\t$contigEdges[$i][1]\t$newEnd\t$contigEdges[$i][3]\t$contigEdges[$i][4]\t$contigEdges[$i][5]\n";

		    $posCtr=0;
		    $mismatches=0;
		    last;

		}


		if($mismatches>0)
		{
		    $posCtr+=1;
		    
		    if($posCtr>3)
		    {
			$posCtr=0;
			$mismatches=0;
		    }

		}

		
		$ctrAdvanced++;
	    }	
	
	}



	if(($contigEdges[$i][3]=~m/endPart/ && $contigEdges[$i][5] eq "+") || ($contigEdges[$i][3]=~m/firstPart/ && $contigEdges[$i][5] eq "-") )
	{	

	    $refString=uc(substr($refGenome{$contigEdges[$i][0]}, $contigEdges[$i][2]-$length, $length));	

	    @refSeq=split(//, $refString);
       	    

		if($contigEdges[$i][5] eq "+")
		{
		    #$refString=substr($refGenome{$contigEdges[$i][0]}, $contigEdges[$i][1], $length);	
		    $contigString=$contigSeq{$contigID};
		}

		if($contigEdges[$i][5] eq "-")
		{
		    $contigString=&reverse_complement_IUPAC($contigSeq{$contigID});
		}


	    $sizeEdge=$contigEdges[$i][2]-$contigEdges[$i][1];

	    
#	    print "\n\n\n";
#	    print "$refString\n";
#	    print "$contigString\n";
	    

	    #my $debug=(" " x ($length-($sizeEdge))).substr($refString, ($length-($sizeEdge)));

#	    print "$debug\n";

#	    print "\n\n\n";
	    
	    @splitContig=split(//, $contigString);

#	    if($contigEdges[$i][3]=~m/116167\_/)
#	    {
#		my $foo=1;

#	    }

	    $index=$length-$sizeEdge;
	    $ctr=0;
	    $mismatches=0;
	    $posCtr=1;

	    $sizeContig=length($contigString);

	    my $ctrAdvanced=0;
	    for($j=$index; $j>=0; $j--)
	    {

		if(!defined($splitContig[$j]) || !defined($refSeq[$j]))
		{
		    my $temp=1;
		}

	       		
		if($splitContig[$j] ne $refSeq[$j])
		{
		    $mismatches+=1;   		   
		}

		if(!defined($splitContig[$j]) || !defined($refSeq[$j]))
		{
		    my $temp=1;
		    
		}




		$ctr++;

		#also print out length of contig edges

		if($mismatches==2 && ($posCtr==2))
		{
		    #$newStart=$contigEdges[$i][2]-($length-$j-1);
		    $newStart=$contigEdges[$i][1]-$ctrAdvanced+2;
		    
		    print $extended "$contigEdges[$i][0]\t$newStart\t$contigEdges[$i][2]\t$contigEdges[$i][3]\t$contigEdges[$i][4]\t$contigEdges[$i][5]\n";
		
		    $posCtr=0;
		    $mismatches=0;
		    last;
		}


		if($mismatches>=2 && ($posCtr==3))
		{
		    #$newStart=$contigEdges[$i][2]-($length-$j-2);
		    $newStart=$contigEdges[$i][1]-$ctrAdvanced+3;
		    
		    print $extended "$contigEdges[$i][0]\t$newStart\t$contigEdges[$i][2]\t$contigEdges[$i][3]\t$contigEdges[$i][4]\t$contigEdges[$i][5]\n";
		
		
		    $posCtr=0;
		    $mismatches=0;
		    last;

		}


		if($mismatches>0)
		{
		    $posCtr+=1;
		    
		    if($posCtr>3)
		    {
			$posCtr=0;
			$mismatches=0;
		    }

		}

		
		$ctrAdvanced++;

	    }	
	
	}





    }

}

sub readFastaContig($)
{

    my $file=shift;
 
    
    open(SEQ, "$file")
	|| die "can't open $file file";




    my($line, $id, $seq, %seqResults, $ctr);

    $ctr=0;
    while($line=<SEQ>)
    {
	chomp($line);

	if($line =~m /^\>.*/)
	{

	    if($ctr>0)
	    {
		$seqResults{$id}=$seq;

		 #print "$id\n";
		
	    }
	 
	   
	    $line =~m /^\>(\d+)_.*/; 
	    $id=$1; 
	    $seq="";
	    $ctr++;
	   
	   
	    
	}
	else
	{
	    $seq=$seq.$line;
	}
    }

    
    $seqResults{$id}=$seq;

    close(SEQ);

    return(%seqResults);

}

sub filterAlignmentsLocal($ $ $)
{

    my $alignmentFile=shift;
    my $cutoff=shift;
    my $contigFile=shift;

    


    my(@samOutput, %contigs, $size, %hardClipped, $i, $cigarString, $sizeClippingFirst, $sizeClippingLast, $sizeContig, $fractionClipped, @allMatches, $sizeClipping);
    my($maxSclipping, $maxHclipping, $sizeSContig, $sizeHContig, $maxClipped, $contigSize);

    open(my $clipped, ">clippedContigs.fa");


    @samOutput=&parseSamFile($alignmentFile);

    %contigs=&readFasta($contigFile);

    $size=@samOutput;
    
    for($i=0; $i<$size; $i++)
    {
	if($samOutput[$i][0]=~m/^\@/)
	{
	    next;
	}
	
#	$sizeClipping=-1;
	$maxSclipping=-1;
	$maxHclipping=-1;
	$sizeSContig=-1;
	$sizeHContig=-1;
	$sizeContig=-1;
	$maxClipped=-1;

	$cigarString=$samOutput[$i][5];

	if($cigarString=~m/(\d+)S/)
	{
	    @allMatches=$cigarString=~m/(\d+)S/g; #getting all soft clipped bases from both ends 

	    $maxSclipping=max(@allMatches); #getting the max amount of soft clipped bases

	    #$sizeClipping=$1;
	    $sizeSContig=length($samOutput[$i][9]); #getting the length of the soft clipped contig

	   
	    
	}
	
	if($cigarString=~m/(\d+)H/)
	{
	    #$sizeClipping=$1;
	    
	    @allMatches=$cigarString=~m/(\d+)H/g; #getting all hard clipped bases from both ends 

	    $maxHclipping=max(@allMatches); #getting the max amount of hard clipped bases

	    #$sizeClipping=$1;
	    $sizeHContig=length($samOutput[$i][9]) + sum(@allMatches); #getting the length of the hard clipped contig


	}

	if($maxHclipping<$maxSclipping)
	{
	    $maxClipped=$maxSclipping;
	}else
	{
	    $maxClipped=$maxHclipping;
	}
	
	if($sizeHContig<$sizeSContig)
	{
	    $sizeContig=$sizeSContig;

	}else
	{
	    $sizeContig=$sizeHContig;

	}




	    if(!defined($sizeClipping))
	    {
		my $temp;
	    }

	    if(($maxClipped/$sizeContig)>=$cutoff && $maxClipped > 0)
	    {
		print $clipped ">".$samOutput[$i][0]."\n";
		print $clipped "$contigs{$samOutput[$i][0]}\n";

		
	    }


    }

    

}



