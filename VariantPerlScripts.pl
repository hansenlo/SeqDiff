#!/usr/bin/perl -w 

use strict; 
use warnings;
use Devel::Size qw(size total_size);
use Bio::Cigar;


do "/home/hansenlo/Code/usefulFunctions.pl";
die $@ if $@;

#subroutine to extract the first and last lines of a sequence from a fasta file
#inputs to function are name of fasta file and integer representing the number
#of bps to extract from beginning and end of sequence
sub extractFasta($ $);

#given a file contianing the alignment for the ends of contigs  find contigs in which the aligned start and end positions are very different 
#than what they should be input is file name containing aligned ends of contigs in bed format and an integer representing 
#how different aligned start and end positions should differ to be included input also includes file of contig sequence
sub filterContigs($ $ $);


#given a file contianing the alignment for the ends of contigs  find contigs in which the aligned start and end positions are very different 
#than what they should be input is file name containing aligned ends of contigs in bed format and an integer representing 
#how different aligned start and end positions should differ 
#last argument is the path to a fasta containg the reference genome variants are called againist
sub filterContigsVer2($ $ $ $);


#given a file containing the alignment of contigs find the SNP in the contig sequence and report location and SNP
sub getSNPS($ $);

#given two strings will calculate the number of mismatches between them
sub calculateHammingDistance($ $);

#function to take as input a vcf file and merge the same variants if they meet a certain cutoff number of reads
#last argument is the name of the file to print merged calls to 
sub mergeVariants($ $ $);

#annoying have to do this because contig names are non standard
#will remove everything but the first number in the contig id and 
#this is what will be use as the key
sub readFastaContig($);


#given a Sam file identify SNPS insertions and deletions
##Important mapping quality is assumed to be in bowtie format!!! i.e. 0-42  
sub getVar_MDstring($);

my $inputFile=shift;
#my $contigSeq=shift;

my $cutoff=shift;

my $contigSeq=shift;

my $flag=shift;

my $genome=shift;
#my  $avg=`samtools mpileup -ABr chr21:16442562-16442689 -d 1000000 /data3/GenomeInABottle/chr21_notAligned.sorted.bam |  awk \'{sum+=\$4} END { print sum/NR}\'`;

#my $test=$avg/2;

#print "average is $test\n";


#exit;

#my $genome="/data/Genomes/human19/allChrhg19InOrder.fa";

# $genome="/data/Genomes/entireHuman19Broad/hs37d5_chrAdded.fa";


#my $genome="/data6/sukrit/081216_MiSeq_MMB1newdel_genomeSeq/MappingToReference/MMB1genomeCIRC84.fasta";

#my $genome="/data/Genomes/cElegans10/allChr.fa";



if($flag eq 'filter')
{

    &filterContigsVer2($inputFile, $contigSeq, $cutoff, $genome);


    #&filterContigs($inputFile, $contigSeq, $cutoff);

    #&getSNPS_MDstring($inputFile);
}

if($flag eq 'extract')
{
    &extractFasta($inputFile, $cutoff);
}


if($flag eq 'SNP')
{

&getSNPS($inputFile, $contigSeq);

}


if($flag eq 'cigar')
{

    &getVar_MDstring($inputFile);

    #merging indels
    &mergeVariants("temp.vcf", 6, "mergedIndels.vcf");

    #merging Snps
    &mergeVariants("temp2.vcf", 6,"mergedSNPs.vcf");


}



#subroutine given two strings will return number of mismatches between them
sub calculateHammingDistance($ $)
{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub extractFasta($ $)
{
    my $fasta=shift;
    my $numNuc=shift; #number of nucleotides to take for start and end of string

    
    
    my (%Seq, $i, $key, $fastaSeq, $chunk, $length, $middle, $mid);

    %Seq=&readFasta($fasta);
    
#loop through the fasta file
    
    foreach $key(keys %Seq){
	
	$fastaSeq=$Seq{$key};

	$key=~s/\s/_/g;
	
	#$chunk=substr($fastaSeq, 0, $numNuc);
	

	$length=length($fastaSeq);

	if($numNuc>$length)
	{
	    $numNuc=$length-1;
	}


	$mid=int($length/2);

	$chunk=substr($fastaSeq, int($mid-($numNuc/2)), $numNuc);

	#print ">".$key."_".$length."_middle\n";
	#print "$chunk\n";

	$chunk=substr($fastaSeq, 0, $numNuc);
	
	print ">".$key."_".$length."_firstPart\n";
	
	print "$chunk\n";
	
	$chunk=substr($fastaSeq, ($numNuc)*-1);
	
	#get end of fasta sequence
	print ">".$key."_".$length."_endPart\n";
	print "$chunk\n";
    }
}



   
sub filterContigs($ $ $)
{
    my $contigsBed=shift;
    my $seq=shift;
    my $cutoff=shift;


    my (@regions, $numContigs, $i, $j, $diff, $header, $line, %firstPart, %endPart, $key, $deletion);
    my ($index, $start, $end, $secondKey, $secondIndex, @temp, $size, $dev, %refSeq, %contigSeqs, $referenceSeq);
    my ($contigID, $refSubSeq, $contigSeq, $matches, $min, $temp, $numCol, @merged, $printString, $varSeq);

    my(@refString, @contigString, $chr, $null, $compareSeq, $minimum, $minIndex, $distance, $refSeq);

    open(OUT, ">inversions.bed");


    #header for vcf file
    print "##fileformat=VCFv4.1\n";
    print "##INFO=<ID=Contig, Number=1, Type=String, Description=\"Contig id the variant was derived from\">\n";
    print "##INFO=<ID=Ends, Number=2, Type=Integer, Description=\"The mapping quality of the two ends of the contig\">\n";
    print "##INFO=<ID=Length, Number=1, Type=Integer, Description=\"The length of the indel\">\n";
    print "##INFO=<ID=Type, Number=1, Type=String, Description=\"The type of indel\">\n";
    print "##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"The type of indel for indels larger than 15 bps\">\n";
    print "##INFO=<ID=SVLEN, Number=1, Type=Integer, Description=\"The size of the indel\">\n";
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";




    @regions=&parseBedfile($contigsBed);

    #human
    %refSeq=&readFasta("/data/Genomes/human19/allChrhg19InOrder.fa");

    #worm
    #%refSeq=&readFasta("/home/gabdank/AF_SOL_490/elegans.genome.fa");
    #%refSeq=&readFasta("/data6/Genomes/cElegans10/allChr.fa");



    %contigSeqs=&readFastaContig($seq);


    $numContigs=@regions;

    #store the start and end of each contig in seperate hash tables
    foreach($i=0; $i<$numContigs; $i++)
    {


	if($regions[$i][0]=~m/chrUn/)
	{
	    next;
	}


	if($regions[$i][0]=~m/random/)
	{
	    next;
	} 

	if($regions[$i][3]=~m/firstPart/)
	{
	    $firstPart{$regions[$i][3]}=$i;
	}elsif($regions[$i][3]=~m/endPart/)
	{
	    $endPart{$regions[$i][3]}=$i;
	}
	
    }

    #pair the second part of each contig with the first part
    #calcuate distance between the two
    foreach $key(keys %firstPart){


	#print "$key\n\n";

	$secondKey=$key;
	$secondKey=~s/firstPart/endPart/;

	$index=$firstPart{$key};


	#get the size the contig is supposed to be
	#$key=~m/size=(\d+)/;
	$key=~m/(.*)\_(\d+)\_(firstPart|endPart)/;

	$contigID=$1;
	$size=$2;

	$contigID=~m/^(\d+)\_/;

	$contigID=$1;
 
#	my $numUnderscore = ($contigID =~ tr/\_//); #counting the number of matches for and underscore
	
#	if($numUnderscore>1)
#	{
#	    $contigID=~s/\_\d+//;
#	}


	#if the end of the contig is not uniquely alignable throw the contig away
	if(exists $endPart{$secondKey})
	{
	    $secondIndex=$endPart{$secondKey};
	}else
	{
	    next;
	}


	#if start and end of contig aligne to different chromosomes throw away

	if($regions[$index][0] ne $regions[$secondIndex][0])
	{
	    next;

	}


	####Uncomment this after processing NIST data
	#if start and end of contig align to different strands keep
#	if($regions[$index][5] ne $regions[$secondIndex][5])
#	{
#	     @temp=@regions[$index, $secondIndex];

	    
#	    &printBed(\@temp);
#	     next;
#	}

	
	#next;

	#calculate size difference between aligned start and end of contigs
	if($regions[$index][5] eq "+")
	{
	    $diff=$regions[$secondIndex][2]-$regions[$index][1];
	     
	    @temp=@regions[$index, $secondIndex];

	    
	    $dev=abs($diff-$size);
	    

	    #compare contig sequence with refernce sequence
	    #if there are many differences print to a file
#	    if(exists($refSeq{$regions[$secondIndex][0]}))
#	    {
		#print "$key\t\t$regions[$secondIndex][0]\t\t$regions[$index][1]\t\t$size\n";


#		if($key eq "NODE_8233_length_78_cov_5.320513_108_firstPart")
#		{
#		    $temp=1;
#		}

#		$refSubSeq=uc(substr($refSeq{$regions[$secondIndex][0]}, $regions[$index][1], $size));
#		$contigSeq=uc($contigSeqs{$contigID});
#	    }else
#	    {
#		next;
#	    }

#	    $matches = ($refSubSeq ^ $contigSeq) =~ tr/\0//;
	    
	   
	    if($dev>=$cutoff)
	    {
	    #contig likely contains an inversion mutation
		if(abs($size-$matches)==$cutoff)
		{
		    
		    if(1==2)
		    {
		    #calculate the position of the mismatch
		    #########Only use when $cutoff is set to one#########
			@refString=split(//, $refSubSeq);
			@contigString=split(//,$contigSeq);

			for($i=0; $i<=$#refString; $i++)
			{
		
			    if($refString[$i] ne $contigString[$i])
			    {
				$chr=$regions[$secondIndex][0];
				$start=$regions[$index][1]+$i;
				$end=$regions[$index][1]+$i+1;
				
				print OUT "$chr\t$start\t$end\n";
				
			    }
			    
			    
			}

		    
		    }

                  ##############end of block#############

		    push @{$temp[0]}, abs($size-$matches);
		    push @{$temp[1]}, abs($size-$matches);

		    for($i=0; $i<=$#temp; $i++)
		    {
	
			$numCol=@{$temp[$i]};

			for($j=0; $j<$numCol; $j++)
			{
			      
			    if($j<($numCol-1))
			    {
				#print OUT "$temp[$i][$j]\t"; 
			    }
			    else
			    {
				#print OUT "$temp[$i][$j]\n";
			    }
			}
			
		    }



		    #print "$contigID\n";
		    #&printBed(\@temp);
		    next;

		   
		}
	    }
	    #next;


	    ###indels start code#####################
	    ##printing out indels for pos strand
	    if($dev>=$cutoff)
	    {
		push @{$temp[0]}, ($diff-$size)*-1;
		push @{$temp[1]}, ($diff-$size)*-1;
		#&printBed(\@temp);

		
#		if(abs(($diff-$size)*-1)<16)
#		{
		     #@merged=([$temp[0][0], $temp[0][1], $temp[1][2], $temp[0][3], $temp[0][4], $temp[1][4], $temp[0][5], $temp[0][6]]);
		

		#if(abs(($diff-$size)*-1)<16)
		#{
		 #   &printBed(\@temp);
		#}
		#&printBed(\@merged);
#		}		

				
		if(1==1)
		{

		    if($dev>16)
		    {
			next;
		    }

	   		
		#calculating exact indel position moving from one side of the contig###
		    if(exists($refSeq{$regions[$secondIndex][0]}))
		    {
		    #start at one end of contig go until hit first missmatch 
		    #$refSubSeq=uc(substr($refSeq{$regions[$secondIndex][0]}, $regions[$secondIndex][1], $size));
		    
			#deletion occured
			$deletion=0;
			if(($diff-$size)>0)
			{
			    $deletion=1;
			}

		    #deletion
			if($deletion)
			{
			    $refSubSeq=uc(substr($refSeq{$regions[$index][0]}, $regions[$index][1], $size+$dev));		    
			}else #insertion
			{
			    $refSubSeq=uc(substr($refSeq{$regions[$index][0]}, $regions[$index][1], $size-$dev));		    
			}


			if($regions[$index][3] eq "209397_GTGGGAGTTTGAAGATGCATCTCAGAAGGAG_70")
			{
			    $temp=1;
			}




			$contigSeq=uc($contigSeqs{$contigID});
		    }else
		    {
			next;
		    }
	    

		    #$refSubSeq=&reverse_complement_IUPAC($refSubSeq);

		    $null="Z"x$dev;
		

		#means deletion occured
		    if($deletion)
		    {
			$compareSeq=$contigSeq;
			$refSeq=$refSubSeq;
		    }
		    else #means insertion occured
		    {
			$compareSeq=$refSubSeq;
			$refSeq=$contigSeq;
		    }




		    $minimum=10000000000000000;
		    $minIndex=-1;
		    for($i=0; $i<=length($compareSeq); $i++)
		    {
			if($deletion || $dev <16)
			{
			    $temp=$compareSeq;
			    substr($temp, $i, 0)=$null;    
		    #print "$temp\n";



			    $distance=calculateHammingDistance($temp, $refSeq);
			    
			    if($distance<$minimum)
			    {
				$minimum=$distance;
				$minIndex=$i;
			    }

			}

		    }
		
		    substr($compareSeq, $minIndex, 0)=$null;
		$refSubSeq=$refSeq;
		$contigSeq=$compareSeq;

		    #$start=$regions[$index][1]+$minIndex-($dev+10);
		    #$end=$start+($dev+10);

		    $start=$regions[$index][1]+$minIndex;
		    $end=$start+$dev;


		    if($dev<16)
		    {
			#print regions as vcf format
		    
			#print "$regions[$secondIndex][0]}\t$start\t$end\n";

			#$temp=($diff-$size)*-1;
			$temp=abs(($diff-$size));
			#print "$regions[$secondIndex][0]\t$start\t$end\t$regions[$index][4]\t$regions[$secondIndex][4]\t$temp\n";

			if($deletion)
			{
			    $varSeq=substr($refSeq, $minIndex-1, $dev+1);
			    $referenceSeq=substr($refSeq, $minIndex-1, 1);
			    $printString="$regions[$secondIndex][0]\t$start\t"."."."\t$varSeq\t$referenceSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];Length=$temp;Type=deletion";
			}else
			{
			    $varSeq=substr($refSeq, $minIndex-1, $dev+1);
			    $referenceSeq=substr($refSeq, $minIndex-1, 1);
			    $printString="$regions[$secondIndex][0]\t$start\t"."."."\t$referenceSeq\t$varSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];Length=$temp;Type=deletion";
			}
			
		    

		}else
		{

		    $temp=abs(($diff-$size));
		    
		    if($deletion)
		    {
			$varSeq=substr($refSeq, $minIndex-1, $dev+1);
			$referenceSeq=substr($refSeq, $minIndex-1, 1);
			$printString="$regions[$secondIndex][0]\t$start\t"."."."\t.\t.\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];SVLEN=$temp;SVTYPE=deletion";
		    }else
		    {
			#$varSeq=substr($refSeq, $minIndex-1, $dev+1);
			#$referenceSeq=substr($refSeq, $minIndex-1, 1);
			#$printString="$regions[$secondIndex][0]\t$start\t"."."."\t.\t.\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];SVLEN=$temp;SVTYPE=insertion";
		    
			$printString="$regions[$secondIndex][0]\t$regions[$index][1]\t"."."."\t.\t.\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];SVLEN=$temp;SVTYPE=insertion";
		    
		    }


		    


		}

		    print "$printString\n";
		###################################################

		}








		next;
	    }

	}

	if($regions[$index][5] eq "-")
	{
	    $diff=$regions[$index][2]-$regions[$secondIndex][1];

	    $dev=abs($diff-$size);

	    @temp=@regions[$index, $secondIndex];

	    #compare contig sequence with refernce sequence
	    #if there are many differences print to a file
	    if(exists($refSeq{$regions[$secondIndex][0]}))
	    {
		$refSubSeq=uc(substr($refSeq{$regions[$secondIndex][0]}, $regions[$secondIndex][1], $size));
		$contigSeq=uc($contigSeqs{$contigID});

		if(!exists($contigSeqs{$contigID}))
		{
		    my $temp=1
		}


	    }else
	    {
		next;
	    }


	    $refSubSeq=&reverse_complement_IUPAC($refSubSeq);

	    $matches = ($refSubSeq ^ $contigSeq) =~ tr/\0//;
	    

	    if($dev>=$cutoff)
	    {
	    #contig likely contains an inversion mutation
		if(abs($size-$matches)==$cutoff)
		{

		    if(1==2)
		    {
		     #calculate the position of the mismatch
		    #########Only use when $cutoff is set to one#########
			@refString=split(//, $refSubSeq);
			@contigString=split(//,$contigSeq);
			
			for($i=0; $i<=$#refString; $i++)
			{
			    
			    if($refString[$i] ne $contigString[$i])
			    {
				$chr=$regions[$secondIndex][0];
				$start=$regions[$secondIndex][1]+(length($refSubSeq)-$i-1);
				$end=$regions[$secondIndex][1]+(length($refSubSeq)-$i);
				
				print OUT "$chr\t$start\t$end\n";
			       
				
			    }
			    
			    
			}

			
		    }

                  ##############end of block#############

		    

		   
		    push @{$temp[0]}, abs($size-$matches);
		    push @{$temp[1]}, abs($size-$matches);
		    
		    for($i=0; $i<=$#temp; $i++)
		      {
	
			  $numCol=@{$temp[$i]};

			  for($j=0; $j<$numCol; $j++)
			  {
			      
			      if($j<($numCol-1))
			      {
				  #print OUT "$temp[$i][$j]\t"; 
			      }
			      else
			      {
				  #print OUT "$temp[$i][$j]\n";
			      }
			  }
		    
		      }

		    #print "$contigID\n";
		    
		    #&printBed(\@temp);
		    next;

		}
	    }

	    #next;

#printing out indels for negative strand
	    if($dev>=$cutoff)
	    {
		#@temp=@regions[$secondIndex, $index];

		#push @{$temp[0]}, ($diff-$size)*-1;
		#push @{$temp[1]}, ($diff-$size)*-1;

		
	if(1==1)
	{
	    if($dev>16)
	    {
		next;
	    }


	   		
		#calculating exact indel position moving from one side of the contig###
		if(exists($refSeq{$regions[$secondIndex][0]}))
		{
		    #start at one end of contig go until hit first missmatch 
		    #$refSubSeq=uc(substr($refSeq{$regions[$secondIndex][0]}, $regions[$secondIndex][1], $size));
		    
		    #set flag for deltion
		    $deletion=0;
		    if(($diff-$size)>0)
		    {
			$deletion=1;
		    }



		    #deletion
		    if($deletion)
		    {
			$refSubSeq=uc(substr($refSeq{$regions[$secondIndex][0]}, $regions[$secondIndex][1], $size+$dev));		    
		    }else #insertion
		    {
			$refSubSeq=uc(substr($refSeq{$regions[$secondIndex][0]}, $regions[$secondIndex][1], $size-$dev));		    
		    }




			$contigSeq=uc($contigSeqs{$contigID});
		}else
		{
		    next;
		}
	    
		if($regions[$index][3]=~m/^510\_TACAATGGTGTGGAAGCTCTATGTGTCGATA\_35.*/)
		{
		    $temp=1;
		}



		$refSubSeq=&reverse_complement_IUPAC($refSubSeq);

		$null="Z"x$dev;
		
		#means deletion occured
		if($deletion)
		{
		    $compareSeq=$contigSeq;
		    $refSeq=$refSubSeq;
		}
		else #means insertion occured
		{
		    $compareSeq=$refSubSeq;
		    $refSeq=$contigSeq;
		}

		$minimum=10000000000000000;
		$minIndex=-1;
		#for($i=0; $i<=length($compareSeq); $i++)
	    for($i=length($compareSeq); $i>=0; $i--)
	    {
		if($deletion || $dev <16)
		{
		    $temp=$compareSeq;
		    substr($temp, $i, 0)=$null;    
		    #print "$temp\n";
		    

		    
		    $distance=calculateHammingDistance($temp, $refSeq);
		    
		    if($distance<$minimum)
		    {
			$minimum=$distance;
			$minIndex=$i;
		    }

		}
	    }
		
	    substr($compareSeq, $minIndex, 0)=$null;
		#$refSubSeq=$refSeq;
		#$contigSeq=$contigSeq;

	    $refSubSeq=$refSeq;
	    $contigSeq=$compareSeq;

	    if($deletion)
	    {
		#means deletion occured
		$start=$regions[$secondIndex][1]+$size-$minIndex;
		$end=$start+($dev);
	    }else
	    {
		$start=$regions[$secondIndex][1]+$size-$minIndex-$dev;
		$end=$start+($dev);
	    }


	     if($dev<16)
		    {
			#print regions as vcf format
		    
			#print "$regions[$secondIndex][0]}\t$start\t$end\n";

			#$temp=($diff-$size)*-1;
			$temp=abs(($diff-$size));
			#print "$regions[$secondIndex][0]\t$start\t$end\t$regions[$index][4]\t$regions[$secondIndex][4]\t$temp\n";

			if($deletion)
			{
			    $varSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, $dev+1));
			    $referenceSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, 1));
			    $printString="$regions[$secondIndex][0]\t$start\t"."."."\t$varSeq\t$referenceSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];Length=$temp;Type=deletion";
			}else
			{
			    #$varSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, $dev+1));
			    #$referenceSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, 1));
			    #$printString="$regions[$secondIndex][0]\t$start\t"."."."\t$referenceSeq\t$varSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];Length=$temp;Type=insertion";
			
			    $start=$regions[$secondIndex][1]+$size-$dev;
			    $printString="$regions[$secondIndex][0]\t$start\t"."."."\t$referenceSeq\t$varSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];Length=$temp;Type=insertion";
			

			    
			}
			
		    

		}else
		{

		    $temp=abs(($diff-$size));
		    
		    if($deletion)
		    {
			$varSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, $dev+1));
			$referenceSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, 1));
			$printString="$regions[$secondIndex][0]\t$start\t"."."."\t.\t.\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];SVLEN=$temp;SVTYPE=deletion";
		    }else
		    {
			$varSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, $dev+1));
			$referenceSeq=&reverse_complement_IUPAC(substr($refSeq, $minIndex, 1));
			$printString="$regions[$secondIndex][0]\t$start\t"."."."\t.\t.\t"."."."\tPASS\tcontig=$contigID;Ends=$regions[$index][4],$regions[$secondIndex][4];SVLEN=$temp;SVTYPE=insertion";
		    }


		    


		}

		    print "$printString\n";





	    #if($dev<16)
	    #{
	#	$temp=($diff-$size)*-1;
	#	print "$regions[$secondIndex][0]\t$start\t$end\t$regions[$index][4]\t$regions[$secondIndex][4]\t$temp\n";

	 #   }
		###################################################

	}
		


#		if(abs(($diff-$size)*-1)<16)
#		{
		     ##@merged=([$temp[0][0], $temp[0][1], $temp[1][2], $temp[0][3], $temp[0][4], $temp[1][4], $temp[0][5], $temp[0][6]]);
		
		#if(abs(($diff-$size)*-1)<16)
		#{
		#    push @{$temp[0]}, ($diff-$size)*-1;
		 #   push @{$temp[1]}, ($diff-$size)*-1;

		  #  &printBed(\@temp);
		#}
		     #&printBed(\@merged);
#		}		

		next;
	    }
	   
	}

	


    }
    

	
    &printBed(\@temp);

    #$temp=1;
}



sub getSNPS($ $)
{
    my $contigsBed=shift;
    my $seq=shift;
   

    my (@regions, $numContigs, $i, $j, $diff, $header, $line, %firstPart, %endPart, $key, @SNPs);
    my ($index, $start, $end, $secondKey, $secondIndex, @temp, $size, $dev, %refSeq, %contigSeqs);
    my ($contigID, $refSubSeq, $contigSeq, $matches, $min, $temp, $numCol, $refBase, $contigBase);

    my(@refString, @contigString, $chr);

    open(OUT, ">SNPs.bed");



    @regions=&parseBedfile($contigsBed);

    #human
    #%refSeq=&readFasta("/data/HumanGenomeHg19/allChrHg19.fa");

    #worm
    #%refSeq=&readFasta("/home/gabdank/AF_SOL_490/elegans.genome.fa");

    %refSeq=&readFasta("/data6/Genomes/cElegans10/allChr.fa");




    %contigSeqs=&readFasta($seq);


    $numContigs=@regions;

    my $ctr=0;
    for($i=0; $i<$#regions; $i++)
    {

	$size=$regions[$i][2]-$regions[$i][1];

	if($regions[$i][0]=~m/chrUn/)
	{
	    next;
	}


	if($regions[$i][0]=~m/random/)
	{
	    next;
	 } 



	$refSubSeq=uc(substr($refSeq{$regions[$i][0]}, $regions[$i][1], $size));
	
	#print($regions[$i][3]);
	$contigSeq=uc($contigSeqs{$regions[$i][3]});


	#if($regions[$i][3]=~m/^2081_/)
	#{
	 #   $temp=1;
	#}



	if($size!=length($contigSeq))
	{
	    next;
	}



	if($regions[$i][5] eq "-")
	{
	    $refSubSeq=&reverse_complement_IUPAC($refSubSeq);
	}

	$matches = ($refSubSeq ^ $contigSeq);
	
	#print($refSubSeq."\n");
	#print("contig seq is \n");
	#print($contigSeq."\n");

	while($matches =~ /[^\0]/g){
	    #print substr($refSubSeq,$-[0],1), ' ', substr($contigSeq,$-[0],1), ' ', $-[0], "\n";
	

	    $refBase=substr($refSubSeq,$-[0],1);
	    $contigBase=substr($contigSeq,$-[0],1);


	    if($regions[$i][5] eq "-")
	    {
		$start=$regions[$i][1]+length($refSubSeq)-$-[0];
		$end=$regions[$i][1]+length($refSubSeq)-$-[0];
		$refBase=&reverse_complement_IUPAC($refBase);
		$contigBase=&reverse_complement_IUPAC($contigBase);

	    }else
	    {

		$start=$regions[$i][1]+$-[0]+1;
		$end=$regions[$i][1]+$-[0]+1;
	    }

	    if($refBase eq "N" || $contigBase eq "N")
	    {
		next;
	    }

	    ##if($start==29832090)
	    ##{
	##	$temp=1;
		##print $regions[$i];
	  ##  }


	    push @{$SNPs[$ctr]}, $regions[$i][0], $start, $end, $refBase, $contigBase;

	    $ctr++;

	    #print($regions[$i][0]."\t".$start."\t".$end."\t".$refBase."\t".$contigBase."\n");
	}

	#exit();
    }
    
    my (@removeIndex, $ctrSame);

    $size=@SNPs;
    @temp=();

    $ctrSame=0;
    #run through SNPs remove those that are clustered in the same place likely an artifact
    for($i=1; $i<$size; $i++)
    {
	if(($SNPs[$i][0] eq $SNPs[$i-1][0]) && abs($SNPs[$i][1]-$SNPs[$i-1][1]+1)<=3)
	{
	    
		$ctrSame++;	
		push @temp, $i;
		push @temp, $i-1;
				

	}else
	{

	    if($ctrSame>3)
	    {
	
		push(@removeIndex, @temp);

		@temp=();

	    }

	    @temp=();
	    $ctrSame=0;
	}
	


    }


    $size=@removeIndex;
    for($i=0; $i<$size; $i++)
    {
	$SNPs[$removeIndex[$i]][0]=-1;	
    }


    #printing the remaining SNPs
    $size=@SNPs;
    for($i=0; $i<$size; $i++)
    {
	$numCol=@{$SNPs[$i]};

	for($j=0; $j<$numCol; $j++)
	{
	    if($SNPs[$i][0] ne -1)
	    {
	 
		if($j<($numCol-1))
		{
		    print "$SNPs[$i][$j]\t"; 
		}
		else
		{
		    print "$SNPs[$i][$j]\n";
		}
	    }

	}

    }
}


sub getVar_MDstring($)
{

    my $samFile=shift;


    my($cmd, $i, $line, $ref, $read, $operation, @data, $delFlag, $insFlag, $refSeq, $altSeq, $refSeqSNP, $altSeqSNP);
    my($delStart, $delEnd, $insStart, $insEnd, $size, $prevPos, $prevContig);

    $refSeq="";
    $altSeq="";

    $cmd="/data/bin/SamTSV/jvarkit/dist-1.128/sam2tsv -r /data/Genomes/human19/allChrhg19InOrder.fa  $samFile > temp.dat";

    #$cmd="/data/bin/SamTSV/jvarkit/dist-1.128/sam2tsv -r /data6/sukrit/081216_MiSeq_MMB1newdel_genomeSeq/MappingToReference/MMB1genomeCIRC84.fasta  $samFile > temp.dat";


    #$cmd="/data/bin/SamTSV/jvarkit/dist-1.128/sam2tsv -r /data/Genomes/cElegans10/allChr.fa  $samFile > temp.dat";



    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";



    open(FIN, "temp.dat")
	|| die "can't open temp.dat file";
    

    

    #vcf file containing snps
    open(SNP, ">SnpCalls.bed");	



    open(Indels, ">indelCalls.bed");	
  


if(1==1)
{
    open(Indels, ">indelCalls.vcf");	
    
    print Indels "##fileformat=VCFv4.2\n";
    print Indels "##reference=hg19\n";
    print Indels "##INFO=<ID=Contig,Number=1,Type=String,Description=\"Contig id the variant was derived from\">\n";
    print Indels "##INFO=<ID=Length,Number=1,Type=Integer,Description=\"The length of the indel\">\n";
    print Indels "##INFO=<ID=Type,Number=1,Type=String,Description=\"The type of the Variant\">\n";
    print Indels "##INFO=<ID=AlleleFraction,Number=1,Type=Float,Description=\"The number of reads in contig cluster variant was derived from divided by the coverage in a 200 bp window centered on variant\">\n";
    print Indels "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n";
    print Indels "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n";
 } 

if(1==1)
{
    open(SNP, ">SnpCalls.vcf");	



    print SNP "##fileformat=VCFv4.2\n";
    print SNP "##reference=hg19\n";
    print SNP "##INFO=<ID=Contig,Number=1,Type=String,Description=\"Contig id the variant was derived from\">\n";
    print SNP "##INFO=<ID=Type,Number=1,Type=String,Description=\"The type of the Variant\">\n";
    print SNP "##INFO=<ID=AlleleFraction,Number=1,Type=Float,Description=\"The number of reads in contig cluster variant was derived from divided by the coverage in a 200 bp window centered on variant\">\n";
    print SNP "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n";
    print SNP "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n";
}



  

{ my $ofh = select Indels;
	  $| = 1;
	  select $ofh;
}


    
    $delFlag=0;
    $insFlag=0;

    $prevPos=0;

    my(@contig, $j, $lineCtr);
    
    $lineCtr=0;

    $line=<FIN>; #reading in the header line
    push @contig, [split(/\t/, $line)];
    $prevContig="N";


    while($line=<FIN>)
    {
	$lineCtr++;

	chomp($line);


	###ONLY USE THIS FOR GENOME IN A BOTTLE###

#	if($line!~/chr21/) #only print out variants on chromosome 21
#	{
#	    next;
#	}

	$prevPos=0;






	#if($#contig>0)
	#{
	$prevContig=$contig[$#contig][0]; #storing the contig before whatever contig is on the next line
	#}

	push @contig, [split(/\t/, $line)]; # get next line of file

	if( ($contig[$#contig][0] ne $prevContig) ) #if on a different contig process the previous contig
	{

	    $delFlag=0;
	    $insFlag=0;

	    $prevPos=0;

	    #print "on contig $contig[$#contig][0] line counter is $lineCtr prevContig is $prevContig\n";

	    #for($j=0; $j<=$#contig; $j++)
	    for($j=($#contig * 0.1); $j<=($#contig-($#contig * 0.1)); $j++)
	    {
	
		#if($j<($#contig * 0.1) or ($j > $#contig-($#contig * 0.1))  ) #do not call variants based on the edges of contigs much more susceptible to alignment artifacts. 
		#{
		 #   next;

		#}
		


	
		#$prevPos=$data[6];

		#if enter contig in the middle of the variant go until the end of the variant before calling variants
		if(($contig[$j][6]!~m/^-?\d+\z/ || $contig[$j][3]!~m/^-?\d+\z/) && $prevPos==0)
		{
		    next;
		}

	      

		if($j>0 && $contig[$j-1][6]=~m/\d+/)
		{
		    $prevPos=$contig[$j-1][6];
		
		}elsif($contig[$j-1][6]=~m/\d+/)
		{
		    $prevPos=$contig[$j][6];
		}


		

		@data=@{$contig[$j]};

		
		#@data=split(/\t/, $line);

	#variant is a snp
		if((uc($data[4]) ne uc($data[7])) && $data[8] eq "M")
		{
		    if(uc($data[4]) eq "N")
		    {
			next;
		    }

		    if($data[0] eq "161419_12")
		    {

			my $temp=1;

		    }

	    
	#	    print SNP "$data[2]\t$data[6]\t$data[6]\t$data[4]\t$data[7]\t$data[0]\n"; 
#		    print SNP "$data[2]\t$data[6]\t$data[6]\t$data[4]\t$data[7]\n"; 

		    #print Indels "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n";

		    $altSeqSNP=uc($data[4]);
		    $refSeqSNP=uc($data[7]);

		    print SNP "$data[2]\t$data[6]\t.\t$refSeqSNP\t$altSeqSNP\t.\tPASS\tContig=$data[0];Type=SNP;AlleleFraction=.\tGT\t0/1\n";


	    
#		    print Indels "$data[2]\t$delStart\t.\t$refSeq\t$altSeq\t.\tPASS\tContig=$data[0];Length=$size;Type=DEL;AlleleFraction=$allelFraction\tGT\t0/1\n"; 



		}
	    
		if($data[8] eq "D")
		{
		    
		    #check this when done
		    if($data[0] eq "152695_49")
		    {
			my $temp=1;

		    }




		    if($delFlag==1)
		    {
			$refSeq=$refSeq.$data[7];
		    }

	    #first base in deleted region
		    if($delFlag==0)
		    {

		

			$delFlag=1;
			$delStart=$data[6];
			$delEnd=$data[6];
			$refSeq=$contig[$j-1][7].$data[7];
			$altSeq=$contig[$j-1][7];
			
		    }else
		    {
		#adding one to size of deleted region 
			$delEnd+=1;
		    }
	

		   

		}

		if($data[8] eq "I")
		{
		    if($insFlag==1)
		    {
			$altSeq=$altSeq.$data[4];
		    }

		    

		    


	    #first base in deleted region
		    if($insFlag==0)
		    {
			$insFlag=1;
			$insStart=$prevPos;
			$insEnd=$insStart+1;
			$refSeq=$contig[$j-1][7];
			$altSeq=$contig[$j-1][7].$data[4];

		    }else
		    {
		#adding one to size of deleted region 
			$insEnd+=1;
		    }
	
		}

	


		if($data[8] ne "I" && $data[8] ne "D")
		{
	   
		    if($delFlag==1)
		    {
			 #$size=($delEnd-$delStart)*-1; # length to be printed for bed format
			$size=($delEnd-$delStart);
			$delStart--;
			$delEnd--;
			#print Indels "$data[2]\t$delStart\t$delEnd\t$size\t$data[0]\n";

			#$data[2]=~s/chr//; #removing chr string from chromosome name in order to match genome in a bottle

			#print Indels "$data[2]\t$delStart\t$delEnd\t$size\n";
			
			#print Indels "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

			#$refSeq=$refSeq.$data[7];

			$refSeq=uc($refSeq);
			$altSeq=uc($altSeq);

			
			my $regionStart=$delStart-100;
			my $regionsEnd=$delStart+100;
			my $avg;
#  $avg=`samtools mpileup -ABr $data[2]:$regionStart-$regionsEnd -d 1000000 /data7/PlatinumAlignments/allPlatium.sorted.bam |  awk \'{sum+=\$4} END { print sum/NR}\'`;
			#my $avg=0.001; #temporary fix remove later

			$data[0]=~m/\_(\d+$)/;

			my $numReadsInCluster=$1;

			$avg=0.0000001;
			#my $allelFraction=$numReadsInCluster/$avg;
			#print "$allelFraction\n";
			
			my $allelFraction=".";



			#if(($numReadsInCluster/$avg)>0.6) #if alternative allel in cluster is greater than 50% of overall coverage call it a homozygote
			#{	
			    print Indels "$data[2]\t$delStart\t.\t$refSeq\t$altSeq\t.\tPASS\tContig=$data[0];Length=$size;Type=DEL;AlleleFraction=$allelFraction\tGT\t0/1\n"; 

			#}else
			#{
			 #   print Indels "$data[2]\t$delStart\t.\t$refSeq\t$altSeq\t.\tPASS\tContig=$data[0];Length=$size;Type=DEL;AlleleFraction=$allelFraction\tGT\t0/1\n"; 
			#}



			#if(($numReadsInCluster/$avg)>0.10)
			#{
			 #   print Indels "$data[2]\t$delStart\t.\t$refSeq\t$altSeq\t.\tPASS\tContig=$data[0];Length=$size;Type=DEL;AlleleFraction=$allelFraction\tGT\t0/1\n"; #vcf format
			#}
			    

			#print Indels "##INFO=<ID=Contig, Number=1, Type=String, Description=\"Contig id the variant was derived from\">\n";
			#print Indels "##INFO=<ID=Length, Number=1, Type=Integer, Description=\"The length of the indel\">\n";
			#print Indels "##INFO=<ID=Type, Number=1, Type=String, Description=\"The type of indel\">\n";

			$delFlag=0;
		    }
 

		    if($insFlag==1)
		    {
			$size=($insEnd-$insStart)*1;
			#print Indels "$data[2]\t$insStart\t$insEnd\t$size\t$data[0]\n";
	    
			#$data[2]=~s/chr//; #removing chr string from chromosome name in order to match genome in a bottle

			#$insStart--;  #do this for bed format
			#$insEnd--; #do this for bed format

			#print Indels "$data[2]\t$insStart\t$insEnd\t$size\n";
	    
			#$altSeq=$altSeq.$data[4];

			$altSeq=uc($altSeq);
			$refSeq=uc($refSeq);


			my $regionStart=$insStart-100;
			my $regionsEnd=$insStart+100;
			
			my $avg;
#my  $avg=`samtools mpileup -ABr $data[2]:$regionStart-$regionsEnd -d 1000000 /data7/PlatinumAlignments/allPlatium.sorted.bam |  awk \'{sum+=\$4} END { print sum/NR}\'`;
			#my $avg=0.001; #temporary fix remove later


			$data[0]=~m/\_(\d+$)/;

			my $numReadsInCluster=$1;
			
			$avg=0.0000001;

			#my $allelFraction=$numReadsInCluster/$avg;
			#print "$allelFraction\n";
			my $allelFraction=".";

			#if(($numReadsInCluster/$avg)>0.10)
			#{

			 #   if(($numReadsInCluster/$avg)>0.6) #if alternative allel in cluster is greater than 50% of overall coverage call it a homozygote
			  #  {	
			#	print Indels "$data[2]\t$insStart\t.\t$refSeq\t$altSeq\t.\tPASS\tContig=$data[0];Length=$size;Type=INS;AlleleFraction=$allelFraction\tGT\t0/1\n";
			 #   }else
			  #  {
				print Indels "$data[2]\t$insStart\t.\t$refSeq\t$altSeq\t.\tPASS\tContig=$data[0];Length=$size;Type=INS;AlleleFraction=$allelFraction\tGT\t0/1\n";
			   # }


			#}
			    

			$insFlag=0;

		    }


	     
		}
	 
	    }

	    @contig=(); #clear the contig array 
	    push @contig, [split(/\t/, $line)]; #add the current line which is the first base of the contig back in
	}

    }

    

    close(Indels);
    close(SNP);



$cmd="time /data/bin/vcflib/bin/vcfstreamsort -a indelCalls.vcf > foo.dat";
system($cmd)==0
    or die "system $cmd failed\n";

$cmd="time /data/bin/vcflib/bin/vcfstreamsort -a SnpCalls.vcf > foo2.dat";
system($cmd)==0
    or die "system $cmd failed\n";



#removing duplicate variant calls
#$cmd="time /data/bin/vcflib/bin/vcfuniq foo.dat > indelCalls.vcf";
#system($cmd)==0
#    or die "system $cmd failed\n";


if(1==2)
{

#if 2 variants have the same position keep only one of them 
$cmd="(head -n 7 indelCalls.vcf && tail -n +3 indelCalls.vcf | sort -u -k1,1 -k2,2) > foo.dat";
system($cmd)==0
    or die "system $cmd failed\n";


$cmd="cp foo.dat indelCalls.vcf";
system($cmd)==0
    or die "system $cmd failed\n";


}



$cmd="cat indelCalls.vcf | vcf-sort > temp.vcf";
system($cmd)==0
    or die "system $cmd failed\n";

$cmd="cat SnpCalls.vcf | vcf-sort > temp2.vcf";
system($cmd)==0
    or die "system $cmd failed\n";




$cmd="bgzip -c temp.vcf > indelCalls.vcf.gz";
system($cmd)==0
    or die "system $cmd failed\n";


$cmd="tabix -p vcf indelCalls.vcf.gz";
system($cmd)==0
    or die "system $cmd failed\n";







    if(1==2)
    {

     #My attempts at parsing md and cigar string will come back to

    my @sam=&parseSamFile($samFile);
    
    my ($i, $size, $ctr, $numFields, $md, $cigar, @characters, $parsedCigar);

    $size=@sam;

    for($i=0; $i<$size; $i++)
    {

	#read was not mapped skip to next
	if(($sam[$i][1] & 4)==4)
	{
	    next;
	}
    
	#if mapping quality is less than 20 than skip to next alignment 
	if($sam[$i][4]<20)
	{
	    next;
	}

	#getting the cigar string
	$cigar=$sam[$i][5];


	#look for the md string need to do this because md field is optional and don't know for sure which column it will be
	$ctr=11;
	
	$numFields=$#{$sam[$i]};
	while($ctr<=$numFields)
	{
	    if($sam[$i][$ctr]=~m/^MD/)
	    {
		$md=$sam[$i][$ctr];
		last;
	    }

	    
	    $ctr++;
	}
	
	$parsedCigar = Bio::Cigar->new($cigar);

	#$parsedCigar = Bio::Cigar->new("2M1D1M1I4M");

#	print "$cigar\n";

	

	
	$ctr=1;

#	@characters=($cigar=~/(\w)/g);

#	if(

	    
	



    }


    }    




}

#function to take as input a vcf file and merge the same variants if they meet a certain cutoff number of reads
sub mergeVariants($ $ $)
{
    my $vcfFile=shift;
    my $cutoff=shift;
    my $outputFile=shift;

    my @vcfMatrix=&parseVCFMatrix($vcfFile);

    my $rowNum=@vcfMatrix;

    my ($i, $sameVarFlag, $sameVarCtr, $covCounter, $j, $sameIndex, $z, $numColumns, $size, $regionsStart, $regionsEnd);
    my($contigID, $allelFraction);
    
    open(Indels, ">$outputFile");	

    


#or $vcfMatrix[$i][7]=~m/.*Contig\=\d+\_\d+\_\d+.*/

    $sameVarFlag=0;
    $sameVarCtr=0;
    for($i=0; $i<($rowNum-1); $i++)
    {


	$covCounter=0;
	if($vcfMatrix[$i][0]=~m/^\#/ )
	{

#	    or $vcfMatrix[$i][7]=~m/.*Contig\=\d+\_\d+\_.*/

	    #print Indels "$vcfMatrix[$i][0]\n"; #header lines have no tabs this should print them out in their entiriry. 
	    $sameVarFlag=0;
	    $sameVarCtr=0;
  

	    $size=@{$vcfMatrix[$i]};
	    #print out one copy of the many variant duplicates
	    for($z=0; $z<$size; $z++)
	    {
		if($z!=($size-1))
		{
		    print Indels "$vcfMatrix[$i][$z]\t";
		}else
		{
		    print Indels "$vcfMatrix[$i][$z]\n"; 
		}	
	    }


	    next;
	}

	if($vcfMatrix[$i][1]==9422116)
	{
	    my $temp=1;
	}
	



	#checking to see if Variants are the same
	if( ($vcfMatrix[$i][0] eq $vcfMatrix[$i+1][0] ) && ($vcfMatrix[$i][1]==$vcfMatrix[$i+1][1]) && ($vcfMatrix[$i][3] eq $vcfMatrix[$i+1][3]) && ($vcfMatrix[$i][4] eq $vcfMatrix[$i+1][4]) )
	{
	    $sameVarFlag=1;
	    $sameVarCtr++;

	}else
	{
	    if($sameVarFlag==1) #have reached the end of a stretch of similar variants
	    {
		$sameIndex=$i-$sameVarCtr;
		$covCounter=0;
		for($j=$sameIndex; $j<=($i); $j++) #sum up the count of reads for every cluster
		{
		    if($vcfMatrix[$j][7]=~/Contig\=\d+_(\d+)\;.*/)
		    {
			if($1>$covCounter)
			{
			    $covCounter=$1;

			    #$covCounter+=$1;
			}
		    }

		    if(!(defined($1)))
		    {
			my $temp=1;
		    }

		    

		}

		$sameVarFlag=0;
		$sameVarCtr=0;
		    

		$regionsStart=$vcfMatrix[$sameIndex][1]-100;
		$regionsEnd=$vcfMatrix[$sameIndex][1]+100;

	
		if($covCounter>0)
		{
		    #my  $avg=`samtools mpileup -ABr $vcfMatrix[$sameIndex][0]:$regionsStart-$regionsEnd -d 1000000 /data7/PlatinumAlignments/allPlatium.sorted.bam |  awk \'{sum+=\$4} END { print sum/NR}\'`;

		    my $avg=0.001;
		     $allelFraction=$covCounter/$avg;
		}else #if all the variants that make the same call come from an impure cluster than do not attempt to make any zygousity calls or allelfraction calculations
		{
		    	for($z=0; $z<$size; $z++)
			{
			    if($z!=($size-1))
			    {
				print Indels "$vcfMatrix[$sameIndex][$z]\t";
			    }else
			    {
				print Indels "$vcfMatrix[$sameIndex][$z]\n"; 
			    }	
			}

			next;


		    $allelFraction=-1;
		}


		#print "$allelFraction\n";

		#if($allelFraction>1)
		#{
		 #   my $temp=1;
		#}


		if($covCounter < $cutoff && $covCounter > 0 ) #not enough reads in cluster to meet the minimum coverage
		{
		    $sameVarFlag=0;
		    $sameVarCtr=0;
		    next;
		}

	
		
		if($allelFraction > 0.10) #remove any low allel fraction variants
		{
		    $vcfMatrix[$sameIndex][7]=~/Contig\=(\d+)_\d+/;
		    $contigID=$1;

		    if(!(defined($contigID)))
		    {
			my $temp=1;
		    }

			#replacing old allel faction and count of reads in the contig with the new counts
		    $vcfMatrix[$sameIndex][7]=~s/Contig\=\d+_\d+/Contig\=$contigID\_$covCounter/;	
		    $vcfMatrix[$sameIndex][7]=~s/AlleleFraction=\./AlleleFraction=$allelFraction/;



		    if($allelFraction>0.5) #if alternative allel in cluster is greater than 50% of overall coverage call it a homozygote
		    {	
			$size=@{$vcfMatrix[$sameIndex]};
			
			#print out one copy of the many variant duplicates
			for($z=0; $z<$size; $z++)
			{
			    if($z!=9)
			    {
				print Indels "$vcfMatrix[$sameIndex][$z]\t";
			    }

			    if($z==9) #reached genotype field
			    {
				print Indels "1/1\n"; 
			    }
			
			}
		    }else
		    {
		
			#heterozygout
			$size=@{$vcfMatrix[$sameIndex]};
			#print out one copy of the many variant duplicates
			for($z=0; $z<$size; $z++)
			{
			    if($z!=9)
			    {
				print Indels "$vcfMatrix[$sameIndex][$z]\t";
			    }

			    if($z==9) #reached genotype field
			    {
				print Indels "0/1\n"; 
			    }
			
			}
			
		       
		    }


		}
			    


	    }else
	    {
		
		##if only have one cluster that calls the variant make a decision on the basis of the the read count for that cluster
				

		    if($vcfMatrix[$i][7]=~/Contig\=\d+_(\d+)\;.*/)
		    {
			$covCounter=$1;

			$regionsStart=$vcfMatrix[$i][1]-100;
			$regionsEnd=$vcfMatrix[$i][1]+100;


			#my  $avg=`samtools mpileup -ABr $vcfMatrix[$i][0]:$regionsStart-$regionsEnd -d 1000000 /data7/PlatinumAlignments/allPlatium.sorted.bam |  awk \'{sum+=\$4} END { print sum/NR}\'`;
			my $avg=0.001;

			my $allelFraction=$covCounter/$avg;
			#print "$allelFraction\n";

		    


			if($covCounter < $cutoff) #not enough reads in cluster to meet the minimum coverage
			{
			    $sameVarFlag=0;
			    $sameVarCtr=0;
			    next;
			}




			if(($covCounter/$avg)>0.10) #remove any low allel fraction variants
			{
			    $vcfMatrix[$i][7]=~m/Contig\=(\d+)_\d+/;
			    $contigID=$1;

			#replacing old allel faction and count of reads in the contig with the new counts
			    $vcfMatrix[$i][7]=~s/Contig\=\d+_\d+/Contig\=$contigID\_$covCounter/;	
			    $vcfMatrix[$i][7]=~s/AlleleFraction=\./AlleleFraction=$allelFraction/;




	
			    if(($covCounter/$avg)>0.5) #if alternative allel in cluster is greater than 50% of overall coverage call it a homozygote
			    {	
				$size=@{$vcfMatrix[$i]};
				#print out one copy of the many variant duplicates
				for($z=0; $z<$size; $z++)
				{
				    if($z!=9)
				    {
					print Indels "$vcfMatrix[$i][$z]\t";
				    }

				    if($z==9) #reached genotype field
				    {
				    	print Indels "1/1\n"; 
				    }
			
				}
			    }else
			    {
		
			#heterozygout
				$size=@{$vcfMatrix[$i]};
			#print out one copy of the many variant duplicates
				for($z=0; $z<$size; $z++)
				{
				    if($z!=9)
				    {
					print Indels "$vcfMatrix[$i][$z]\t";
				    }

				    if($z==9) #reached genotype field
				    {
					print Indels "0/1\n"; 
				    }
			
				}
			
		       
			    }

			

		    

			}

		    }else #if variant is derived from a impure cluster do not attempt to make a zygosity call 
		    {
			for($z=0; $z<$size; $z++)
			{
			    if($z!=($size-1))
			    {
				print Indels "$vcfMatrix[$i][$z]\t";
			    }else
			    {
				print Indels "$vcfMatrix[$i][$z]\n"; 
			    }	
			}



		    }


	    }

	}



    }
   
    
    
    close(Indels);


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

sub filterContigsVer2($ $ $ $)
{
    my $contigsBed=shift;
    my $seq=shift;
    my $cutoff=shift;
    my $genomeFile=shift;

    my (@regions, $numContigs, $i, $j, $diff, $header, $line, %firstPart, %endPart, $key, $deletion);
    my ($index, $start, $end, $secondKey, $secondIndex, @temp, $size, $dev, %refSeq, %contigSeqs, $referenceSeq);
    my ($contigID, $refSubSeq, $contigSeq, $matches, $min, $temp, $numCol, @merged, $printString, $varSeq);

    my(@refString, @contigString, $chr, $null, $compareSeq, $minimum, $minIndex, $distance, $refSeq, $length, $endQuality, %genome);

#    open(OUT, ">inversions.bed");

    %contigSeqs=&readFasta($seq);
    %genome=&readFasta($genomeFile);


if(1==2)
{
    #header for vcf file
    print "##fileformat=VCFv4.1\n";
    print "##INFO=<ID=Contig, Number=1, Type=String, Description=\"Contig id the variant was derived from\">\n";
    print "##INFO=<ID=Ends, Number=2, Type=Integer, Description=\"The mapping quality of the two ends of the contig\">\n";
    print "##INFO=<ID=Length, Number=1, Type=Integer, Description=\"The length of the indel\">\n";
    print "##INFO=<ID=Type, Number=1, Type=String, Description=\"The type of indel\">\n";
    print "##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"The type of indel for indels larger than 15 bps\">\n";
    print "##INFO=<ID=SVLEN, Number=1, Type=Integer, Description=\"The size of the indel\">\n";
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}



#Justins GIB header
if(1==1)
{
    #header for vcf file
    print "##fileformat=VCFv4.1\n";
    print "##INFO=<ID=Contig, Number=1, Type=String, Description=\"Contig id the variant was derived from\">\n";
    print "##INFO=<ID=Ends, Number=2, Type=Integer, Description=\"The mapping quality of the two ends of the contig\">\n";
    print "##INFO=<ID=SVLEN, Number=1, Type=Integer, Description=\"The length of the structural Variant\">\n";
    print "##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"The structural variant type\">\n";
    print "##INFO=<ID=END, Number=1, Type=String, Description=\"The end location of the variant\">\n";
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}




    @regions=&parseBedfile($contigsBed);

    #human
    #%refSeq=&readFasta("/data/Genomes/human19/allChrhg19InOrder.fa");

    #worm
    #%refSeq=&readFasta("/home/gabdank/AF_SOL_490/elegans.genome.fa");
    #%refSeq=&readFasta("/data6/Genomes/cElegans10/allChr.fa");



    #%contigSeqs=&readFastaContig($seq);


    $numContigs=@regions;

    #store the start and end of each contig in seperate hash tables
    foreach($i=0; $i<$numContigs; $i++)
    {

	if($regions[$i][4]<30)
	{
	    next;
	}


	if($regions[$i][0]=~m/chrUn/)
	{
	    next;
	}


	if($regions[$i][0]=~m/random/)
	{
	    next;
	} 

	if($regions[$i][3]=~m/firstPart/)
	{
	    $firstPart{$regions[$i][3]}=$i;

	}elsif($regions[$i][3]=~m/endPart/)
	{
	    $endPart{$regions[$i][3]}=$i;
	}
	
    }

    #pair the second part of each contig with the first part
    #calcuate distance between the two
    foreach $key(keys %firstPart){


	#print "$key\n\n";

	$secondKey=$key;
	$secondKey=~s/firstPart/endPart/;

	$index=$firstPart{$key};





	#get the size the contig is supposed to be
	#$key=~m/size=(\d+)/;
	$key=~m/(.*)\_(\d+)\_(firstPart|endPart)/;

	$contigID=$1;
	$size=$2;
	
	if($contigID=~m/5851255/)
	{
	    my $temp=1;
	}

#	$contigID=~m/^(\d+)\_/;

#	$contigID=$1;
 
	
	$length=length($contigSeqs{$contigID});

	#if the end of the contig is not uniquely alignable throw the contig away
	if(exists $endPart{$secondKey})
	{
	    $secondIndex=$endPart{$secondKey};
	}else
	{
	    next;
	}


	#if start and end of contig aligne to different chromosomes throw away

	if($regions[$index][0] ne $regions[$secondIndex][0])
	{
	    next;

	}


	####Uncomment this after processing NIST data
	#if start and end of contig align to different strands keep
	if($regions[$index][5] ne $regions[$secondIndex][5])
	{
	    next;
	}	
	
#{
	#$diff=abs($regions[$index][1]-$regions[$secondIndex][2]);

	#push @{$regions[$index]}, $diff;   
	#push @{$regions[$secondIndex]}, $diff;   
	
	
	 #    @temp=@regions[$index, $secondIndex];

	    
	  #  &printBed(\@temp);
	   #  next;
#	}

	
	#next;

	


	#calculate size difference between aligned start and end of contigs
	if($regions[$index][5] eq "+")
	{
	    #$diff=$regions[$secondIndex][2]-$regions[$index][1];
	     
	    #$diff=abs($regions[$index][1]-$regions[$secondIndex][2]);
	    
	    $diff=$regions[$secondIndex][2]-$regions[$index][1];

	    if($contigID=~m/.*3967490.*/)
	    {
		my $temp=1;
	    }



	     
	    #@temp=@regions[$index, $secondIndex];
	    





	    $dev=abs($diff-$length);
	    $start=$regions[$index][2];
	    $end="$regions[$index][4],$regions[$secondIndex][4];";
	    $endQuality="$regions[$index][4],$regions[$secondIndex][4];";

	    #for bed format only
	    $start=$regions[$index][2];
	    #$end=$regions[$secondIndex][1];
	    $end=$start;


	}    

#$key=~m/(.*)\_(\d+)\_(firstPart|endPart)/;
	

	if($regions[$index][5] eq "-")
	{
	    if($contigID=~m/.*1175227.*/)
	    {
		my $temp=1;
	    }

	    #$diff=$regions[$secondIndex][2]-$regions[$index][1];
	     
	    #$diff=abs($regions[$index][2]-$regions[$secondIndex][1]);
	     
	    $diff=$regions[$index][2]-$regions[$secondIndex][1];
	    

	    #@temp=@regions[$index, $secondIndex];

	    
	    $dev=abs($diff-$length);
	    $start=$regions[$secondIndex][2];
	    $end="$regions[$index][4],$regions[$secondIndex][4];";
	    $endQuality="$regions[$index][4],$regions[$secondIndex][4];";

	    #for bed format only
	    $start=$regions[$secondIndex][2];
	    #$end=$regions[$index][1];
	    $end=$start;

	    
	}

	
	#if($dev>16)
	#{
	 #   next;
	#}

   
	    if($dev>=$cutoff)
	    {
		#$start=$start-5;
		#$end=$end+5;


		


	
		if($start>$end)
		{
		    my $temp=$start;
		    $start=$end;
		    $end=$start;

		    #if($dev<=3)
		    #{
		    #}

		}


		#uncomment if looking to compare variants across samples
		if(1==2)
		{
		    $start=$start-10;
		    $end=$end+10;
		}


	#deletion occured
		if(($diff-$length)>0)
		{
		    $varSeq=".";
		    $referenceSeq=".";







		   
			    
		    #$printString="$regions[$secondIndex][0]\t$start\t"."."."\t$varSeq\t$referenceSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$end;Length=$dev;Type=deletion";

		    #caculate ratio of heterzyoute versus homozygoute snps for the deleted regions

		    if(1==2)
		    {
			open(OUT, ">foo.bed");

			my $temp=$start+$dev;
		    #my $cmd="echo $regions[$index][0]\t$start\t$temp > foo.dat ";
		    
		    #system($cmd)==0
	#		or die "system $cmd failed\n";

			print OUT "$regions[$index][0]\t$start\t$temp\n";

			close(OUT);

		    #get the heterozygout and homozygout snp count for variants that fall in the deleted regions
			my $cmd="/data/bin/bedtools2/bin/intersectBed -u -a /data3/GenomeInABottle/ZgosityGenomeInBottleV1.19.bed -b foo.bed > foo2.bed";
			system($cmd)==0
			    or die "system $cmd failed\n";
    
			my @output=`wc -l < foo2.bed`;

			my ($hetCount, $homoCount, $total, $zygoRatio);
			
			if($output[0]>0)
			{
			    @output=`grep \"0/1\" foo2.bed | wc -l`;
			    $output[0]=~m/(\d+).*/;
			    $hetCount=$1;

			    @output=`grep \"1/1\" foo2.bed | wc -l`;
			    $output[0]=~m/(\d+).*/;
			    $homoCount=$1;


			    $total=0;
			    $total=$hetCount+$homoCount;

			    if($homoCount>0)
			    {
				$zygoRatio=$hetCount/$homoCount;
				$zygoRatio=$zygoRatio."/".$hetCount."/".$homoCount;

			    }else
			    {
				$zygoRatio=-1;
				$zygoRatio=$zygoRatio."/".$hetCount."/".$homoCount;


			    }
		       			
			}else
			{
			    $zygoRatio=-1;
			    $zygoRatio=$zygoRatio."/"."0"."/"."0";

			}

			print "$regions[$index][0]\t$start\t$end\t"."$contigID"."ZygosityRatio".$zygoRatio."deletionsize_"."$dev\n";
		    

		    }


#		    print "##fileformat=VCFv4.1\n";
#		    print "##INFO=<ID=Contig, Number=1, Type=String, Description=\"Contig id the variant was derived from\">\n";
#		    print "##INFO=<ID=Ends, Number=2, Type=Integer, Description=\"The mapping quality of the two ends of the contig\">\n";
#		    print "##INFO=<ID=Length, Number=1, Type=Integer, Description=\"The length of the indel\">\n";
#		    print "##INFO=<ID=Type, Number=1, Type=String, Description=\"The type of indel\">\n";
#		    print "##INFO=<ID=SVTYPE, Number=1, Type=String, Description=\"The type of indel for indels larger than 15 bps\">\n";
#		    print "##INFO=<ID=SVLEN, Number=1, Type=Integer, Description=\"The size of the indel\">\n";
#		    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
		
		    #$varSeq=".";
		    #$referenceSeq=".";

		    $start=$start-1;

		    $varSeq=uc(substr($genome{$regions[$secondIndex][0]}, $start-1, 1) );
		    $referenceSeq=uc(substr($genome{$regions[$secondIndex][0]}, $start-1, $dev));

		    if($dev>1000000)
		    {
			next;
		    }

		    #making sure the reference sequence has only allowed bases
#		    if($referenceSeq !~ m/A|C|G|T|a|c|g|t|N/)
#		    {
#			next;
#		    }


		    ##to filter SV variants that only occur on fasta records that hve chr in the title super dangerous be careful!!!
		    if($regions[$secondIndex][0] !~ m/chr/)
		    {
		#	next;
		    }


		    my $delEnd=$start+$dev;

		    #$printString="$regions[$secondIndex][0]\t$start\t"."."."\t$referenceSeq\t$varSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$endQuality;Length=$dev;Type=deletion\tGT\t1|1";

		    $printString="$regions[$secondIndex][0]\t$start\t"."."."\t$referenceSeq\t$varSeq\t"."."."\tPASS\tcontig=$contigID;SVLEN=$dev;SVTYPE=DEL;END=$delEnd\tGT\t".".";


		    #print "$dev\n";

		    print "$printString\n";
		    
		    #print "$regions[$index][0]\t$start\t$end\n";



		    #if($dev<300)
		    #{
			#print "$printString\n";
		
		    #}
		}

		#insertion occured
		if(($diff-$length)<0)
		{

		    $varSeq=".";
		    $referenceSeq=".";

		    #print "$regions[$index][0]\t$start\t$end\t"."$contigID"."insertionsize_"."$dev\n";

		    #$varSeq=uc(substr($genome{$regions[$secondIndex][0]}, $start, 1) );
		    #$referenceSeq=uc(substr($genome{$regions[$secondIndex][0]}, $start, $dev));

		    

		    #print "$regions[$index][0]\t$start\t$end\n";

		   # $printString="$regions[$secondIndex][0]\t$start\t"."."."\t$varSeq\t$referenceSeq\t"."."."\tPASS\tcontig=$contigID;Ends=$end;Length=$dev;Type=insertion\tGT\t1|1";

		    #print "$printString\n";

		}
		   
	    }
	
	    #next;

	
    }




}
