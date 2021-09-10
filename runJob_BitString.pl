#!/usr/bin/perl -w 

use strict; 
use warnings;
use Devel::Size qw(size total_size);

do "/home/hansenlo/Code/usefulFunctions.pl";
die $@ if $@;


my $expFastq=shift; #the file containing the fastq reads for the experiment
my $expKmerCounts=shift; #the file containing the counts of how often each kmer occurs
my $controlFile=shift; #the file containing the control reads if present


#given a contig file and a cutoff indicating how much to start aligning from each end
#iteratively increase the size of the end alignments until you get a unique alignment or you reach
#a preset size limit
sub iterativeAlignment($ $ $ $ $);


#given a set of contig edges that align uniquely 
#extend the contig edges until they reach a cutoff number of mismatches with the reference
#arguments are contig edges cutoff reference genome and contig seq
sub extendAlignments($ $ $ $);

#Program to take as input a file name containing a read library and running all commands on it

my($cmd, $line, $i, $header, $outputUnique, $outputOverlapEdena, $contigEdgesOutput, $alignmentOutput, $tempIn, $tempOut);
my($allReads, @splitLine, $path, $outputKmerLib, $dataset, $blastDB);

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


$blastDB="/data/Bowtie2Index/hg19GATK";


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



if(1==2)
{

#getting those reads that passed the filter running the code
#$tempIn="/data/HEM0013-4/".$file;

#$tempIn=$file;

#$outputKmerCount="/data4/HEM0013_158/uniqueKmerCount_unique_HEM0013-158-MC_allReads.dat";

#$cmd="time ./findDiff_HashTable_filterReads $tempIn 20 0 $outputKmerCount > "."/home/hansenlo/SeqDiff/Results/".$outputUnique;

#$cmd="time ./findDiff_BitString_HashTableVer $controlFile $file  31 0 $cutoff temp.dat > "."/home/hansenlo/SeqDiff/Results/".$outputUnique;

$cmd="time ./variantFinder $expKmerCounts $expFastq 45";

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


#copying files to appropriate places
if(1==2)
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


if(1==2)
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

$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB "."/data/uniqueReads.fastq --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_uniqueReads.fastq "."-S $alignmentOutput";




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



if(1==1)
{
######################Aligning the contigs#####################
#$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."refGenomeOnly_contigEdges.sam";

$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."_contigs_bowtie2.sam";

my $contigs="/data5/SeqDiffResults/Results/".$header."_contigs.fasta";


#$cmd="rm /data5/SeqDiffResults/Results/Alignment/".$header."*contigs*";  ##############Be careful with THIS!!!!! make sure you are not deleting the wrong files!!!!

print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";


#$blastDB="/home/hansenlo/Genomes/Bowtie2Index/hg19";

#$blastDB="/data/BowtieIndex/elegansCe10";

#$blastDB="/home/gabdank/genomes_fasta/elegans_W235_index";

$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2   -p 20 -x $blastDB -f $contigs --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta "."-S $alignmentOutput";

#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 --mp 10 --rdg 3,2 --rfg 3,2  -p 20 -x $blastDB -f $contigs --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_Contigs.fasta "."-S $alignmentOutput";



#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 4 -x /home/hansenlo/Genomes/Bowtie2Index/seqA -f $contigEdgesOutput --un /data5/SeqDiffResults/Results/Alignment/".$header."_unAlignable_ContigEdges.fasta "."-S $alignmentOutput";


print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";

#getting bam files for contig edges alignment
my $contigheader=$header."_contigs_bowtie2";
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

#convert to bed file all alignments
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader.".sorted.bed";
$cmd="bamToBed -i $tempIn > $tempOut";
print "#######cmd is $cmd\n\n";
system($cmd)==0
    or die "system $cmd failed\n";



#creating a unique alignment
$tempOut="/data5/SeqDiffResults/Results/Alignment/".$contigheader.".sorted.unique.bam";
$cmd="samtools view -bq 7 $tempIn > $tempOut";
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

sub iterativeAlignment($ $ $ $ $)
{
    my $contigs=shift;
    my $startSize=shift;
    my $endSize=shift;
    my $blastDB=shift;
    my $header=shift;


    my(@samOutput, $flag, $size, $i, @readsToReAlign, %contigs, $contigID, $locationContigs, $j);
    my($seq, $ctr, $start, $end, $strand, %multi, $key, $k);


    $ctr=0;
    
    system("rm mappedReads_IDs.dat");
    open(wasMapped, ">mappedReads_IDs.dat");


    system("cp $contigs needToBeMapped.fasta");
    #open(needMapped, ">needToBeMapped.fasta");
    

    open(willNotMap, ">finalMultiMapping.bed");


    %contigs=&readFasta($contigs);


    #$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/temp_contigEdges.sam";


    $alignmentOutput="/home/hansenlo/SeqDiff/Test/testAligned.out.sam";


    #$alignmentOutput="/data5/SeqDiffResults/Results/Alignment/".$header."_contigEdges.sam";


    #$locationContigs="/data5/SeqDiffResults/Results/".$header."_contigs.fasta";


    	#getting contig edges
    $cmd="time perl VariantPerlScripts.pl $contigs $startSize temp.dat extract  > needToBeMapped.fasta";

    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";
		

    #exit;



    for($j=$startSize+1; $j<=$endSize; $j++)
    {


	#$cmd="time /home/hansenlo/bin/bowtie2-2.0.0-beta7/bowtie2 -p 20 -x $blastDB -f needToBeMapped.fasta "." -S $alignmentOutput";

	#only unique alignments
	#$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes --readFilesIn needToBeMapped.fasta --runThreadN 20 --outFileNamePrefix Test/test --outFilterMultimapNmax 1 --outFilterMismatchNmax 6"; 

	#$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes --outSAMunmapped Within --readFilesIn needToBeMapped.fasta --runThreadN 20 --outFileNamePrefix Test/test --outFilterMismatchNmax 6"; 

	$cmd="/data/bin/STAR_2.3.0e.Linux_x86_64/STAR --genomeDir /data/StarGenomes/hg19 --outSAMunmapped Within --readFilesIn needToBeMapped.fasta --runThreadN 20 --outFileNamePrefix Test/test --outFilterMismatchNmax 6"; 


	
	print "#######cmd is $cmd\n\n";
	system($cmd)==0
	    or die "system $cmd failed\n";

	

	
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


	for($i=0; $i<$size; $i++)
	{

	#skip header lines
	    if($samOutput[$i][0]=~m/^\@/)
	    {

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
	    if($samOutput[$i][4]==255 && $samOutput[$i][5]!~/S/)
	    {
	    


		$start=$samOutput[$i][3]-1;
		$end=$samOutput[$i][3]+$j-1;

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

		
		

		print wasMapped "$samOutput[$i][2]\t$start\t$end\t$samOutput[$i][0]\t$samOutput[$i][4]\t$strand\n";
	    }elsif(($samOutput[$i][1] & 256)!=256)
	    {
		
		$samOutput[$i][0]=~m/(.*)\_(\d+)\_(firstPart|endPart)/;

		$contigID=$1;

		print $needMapped ">$samOutput[$i][0]\n";
		
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


    #adding final set of non unique mappings
    open(noMap, ">>finalMultiMapping.bed");

    
    %multi=readFasta("needToBeMapped.fasta");

    foreach $key (keys %multi)
    {
	
	print noMap "chrI\t1000\t2000\t$key\t0\t+\n";
    }


    print "will not map is $ctr\n";

}


sub extendAlignments($ $ $ $)
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

	$contigEdges[$i][3]=~m/(.*)\_(\d+)\_(firstPart|endPart)/;

	$contigID=$1;
	$length=$2;
	$mismatches=0;
	$posCtr=0;

    


	if(($contigEdges[$i][3]=~m/firstPart/ && $contigEdges[$i][5] eq "+") || ($contigEdges[$i][3]=~m/endPart/ && $contigEdges[$i][5] eq "-") )
	{	


	    $refString=uc(substr($refGenome{$contigEdges[$i][0]}, $contigEdges[$i][2]-2, $length));	

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





	    my $temp=substr($refGenome{$contigEdges[$i][0]}, ($contigEdges[$i][2]-2), $length);	

	    $sizeEdge=$contigEdges[$i][2]-$contigEdges[$i][1];

	    
#	    print "\n\n\n";
#	    print "$refString\n";
#	    print "$contigString\n";
	    

	    my $debug= (" " x ($sizeEdge-2)).$temp;

#	    print "$debug\n";

#	    print "\n\n\n";
	    
	    @splitContig=split(//, $contigString);

	    #if($contigEdges[$i][3]=~m/133774\_/)
	    #{
	#	my $foo=1;

	 #   }

	    $index=$sizeEdge-2;
	    $ctr=0;
	    $mismatches=0;
	    $posCtr=1;

	    $sizeContig=length($contigString);

	    for($j=$index; $j<$sizeContig; $j++)
	    {
	       		
		if($splitContig[$j] ne $refSeq[$ctr])
		{
		    $mismatches+=1;   		   
		}

		$ctr++;

		if($mismatches==2 && ($posCtr==2))
		{
		    $newEnd=$contigEdges[$i][1]+$j-2;
		    print $extended "$contigEdges[$i][0]\t$contigEdges[$i][1]\t$newEnd\t$contigEdges[$i][3]\t$contigEdges[$i][4]\t$contigEdges[$i][5]\n";
		
		    $posCtr=0;
		    $mismatches=0;
		    last;
		}


		if($mismatches>=2 && ($posCtr==3))
		{
		    $newEnd=$contigEdges[$i][1]+$j-3;
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

		

	    }	
	
	}



	if(($contigEdges[$i][3]=~m/endPart/ && $contigEdges[$i][5] eq "+") || ($contigEdges[$i][3]=~m/firstPart/ && $contigEdges[$i][5] eq "-") )
	{	

	    $refString=uc(substr($refGenome{$contigEdges[$i][0]}, $contigEdges[$i][2]-$length-1, $length));	

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
	    

	    my $debug=(" " x ($length-($sizeEdge))).substr($refString, ($length-($sizeEdge)));

#	    print "$debug\n";

#	    print "\n\n\n";
	    
	    @splitContig=split(//, $contigString);

	    #if($contigEdges[$i][3]=~m/133774\_/)
	    #{
	#	my $foo=1;

	 #   }

	    $index=$length-$sizeEdge;
	    $ctr=0;
	    $mismatches=0;
	    $posCtr=1;

	    $sizeContig=length($contigString);

	    for($j=$index; $j>=0; $j--)
	    {
	       		
		if($splitContig[$j] ne $refSeq[$j])
		{
		    $mismatches+=1;   		   
		}

		$ctr++;

		if($mismatches==2 && ($posCtr==2))
		{
		    $newStart=$contigEdges[$i][2]-($length-$j-1);
		    print $extended "$contigEdges[$i][0]\t$newStart\t$contigEdges[$i][2]\t$contigEdges[$i][3]\t$contigEdges[$i][4]\t$contigEdges[$i][5]\n";
		
		    $posCtr=0;
		    $mismatches=0;
		    last;
		}


		if($mismatches>=2 && ($posCtr==3))
		{
		    $newStart=$contigEdges[$i][2]-($length-$j-2);
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

		

	    }	
	
	}





    }

}
