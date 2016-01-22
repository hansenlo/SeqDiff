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


#given a file containing the alignment of contigs find the SNP in the contig sequence and report location and SNP
sub getSNPS($ $);

#given two strings will calculate the number of mismatches between them
sub calculateHammingDistance($ $);




#given a Sam file identify SNPS insertions and deletions
##Important mapping quality is assumed to be in bowtie format!!! i.e. 0-42  
sub getVar_MDstring($);

my $inputFile=shift;
#my $contigSeq=shift;

my $cutoff=shift;

my $contigSeq=shift;

my $flag=shift;

if($flag eq 'filter')
{
    &filterContigs($inputFile, $contigSeq, $cutoff);

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
    #%refSeq=&readFasta("/data/HumanGenomeHg19/allChrHg19.fa");

    #worm
    #%refSeq=&readFasta("/home/gabdank/AF_SOL_490/elegans.genome.fa");
    %refSeq=&readFasta("/data6/Genomes/cElegans10/allChr.fa");



    %contigSeqs=&readFasta($seq);


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

				
		if(1==2)
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

		
	if(1==2)
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


    my($cmd, $i, $line, $ref, $read, $operation, @data, $delFlag, $insFlag);
    my($delStart, $delEnd, $insStart, $insEnd, $size, $prevPos);

    $cmd="/data/bin/SamTSV/jvarkit/dist-1.128/sam2tsv -r /data/Genomes/human19/allChrhg19InOrder.fa  $samFile > temp.dat";

    print "#######cmd is $cmd\n\n";
    system($cmd)==0
	or die "system $cmd failed\n";



    open(FIN, "temp.dat")
	|| die "can't open temp.dat file";
    

    

    #vcf file containing snps
    open(SNP, ">SnpCalls.bed");	

    open(Indels, ">indelCalls.bed");	
    
    $delFlag=0;
    $insFlag=0;

    $prevPos=0;

    while($line=<FIN>)
    {

	chomp($line);
	
	$prevPos=$data[6];
	@data=split(/\t/, $line);

	#variant is a snp
	if((uc($data[4]) ne uc($data[7])) && $data[8] eq "M")
	{
	    if(uc($data[4]) eq "N")
	    {
		next;
	    }

	    print SNP "$data[2]\t$data[6]\t$data[6]\t$data[4]\t$data[7]\n"; 

	}
	
	if($data[8] eq "D")
	{

	    #first base in deleted region
	    if($delFlag==0)
	    {
		$delFlag=1;
		$delStart=$data[6]-1;
		$delEnd=$data[6];

	    }else
	    {
		#adding one to size of deleted region 
		$delEnd+=1;
	    }
	


	}

	if($data[8] eq "I")
	{

	    #first base in deleted region
	    if($insFlag==0)
	    {
		$insFlag=1;
		$insStart=$prevPos-1;
		$insEnd=$insStart+1;

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
		$size=($delEnd-$delStart)*-1;
		$delStart--;
		$delEnd--;
		print Indels "$data[2]\t$delStart\t$delEnd\t$size\n";
		
		$delFlag=0;
	    }
 

	    if($insFlag==1)
	    {
		$size=($insEnd-$insStart)*1;
		print Indels "$data[2]\t$insStart\t$insEnd\t$size\n";
	    
		$insFlag=0;

	    }


	     
	}


	 
    }



    

    close(Indels);
    close(SNP);


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
