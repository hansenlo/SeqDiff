#!/bin/bash
#Shell script to run the genome in a bottle analysis intersecting snp calls with genome in a  bottle snp calls


CMD="sed -i 's/chr//g' SnpCalls.bed" #removing chr from the snp calls

eval $CMD


CMD="sed -i 's/chr//g' indelCalls.bed" #removing chr from the snp calls

eval $CMD


##############################snp calls

sort SnpCalls.bed | uniq > temp.bed # removing duplicate elements

mv temp.bed SnpCalls.bed

grep "^21" SnpCalls.bed > temp.bed #get only chr21

mv temp.bed SnpCalls.bed

# CMD="intersectBed -a SnpCalls.bed -b /data3/GenomeInABottle/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed > temp.bed"

CMD="intersectBed -a SnpCalls.bed -b /data3/GenomeInABottle/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed > temp.bed"


eval $CMD



mv temp.bed SnpCalls_NistRegions.bed

intersectBed -a SnpCalls_NistRegions.bed -b /data3/GenomeInABottle/allSNPs_only_NIST_v2.19.recode.vcf > SnpCalls_TruePositives.bed #find the snps that intersect

intersectBed -v -a SnpCalls_NistRegions.bed -b /data3/GenomeInABottle/allSNPs_only_NIST_v2.19.recode.vcf > SnpCalls_FalsePositives.bed #find the snps that do not intersect


intersectBed -v -b SnpCalls_NistRegions.bed -a /data3/GenomeInABottle/allSNPs_only_NIST_v2.19.recode.vcf > Snp_NistCalls_FalseNegatives.vcf #find the snps that intersect

sed -i -e 's/^/chr/' SnpCalls_FalsePositives.bed #add chr back to each line because db snp is in that format

intersectBed -a SnpCalls_FalsePositives.bed -b /home/hansenlo/SeqDiff/Annotations/dbSNP_138_UCSC.bed > SnpCalls_FalsePositives_intersectDBSnp.bed #find the snps that intersect


intersectBed -v -a SnpCalls_FalsePositives.bed -b /home/hansenlo/SeqDiff/Annotations/dbSNP_138_UCSC.bed > SnpCalls_FalsePositives_DoNot_intersectDBSnp.bed #find the snps that are not found in DBsnp



######################indel calls

sort indelCalls.bed | uniq > temp.bed # removing duplicate elements

mv temp.bed indelCalls.bed

grep "^21" indelCalls.bed > temp.bed #get only chr21

mv temp.bed indelCalls.bed

CMD="intersectBed -a indelCalls.bed -b /data3/GenomeInABottle/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed > temp.bed"

eval $CMD


mv temp.bed indelCalls_NistRegions.bed


intersectBed -a indelCalls_NistRegions.bed -b /data3/GenomeInABottle/allIndels_only_NIST_v2.19.recode.vcf > indelCalls_TruePositives.bed #find the snps that intersect

intersectBed -v -a indelCalls_NistRegions.bed -b /data3/GenomeInABottle/allIndels_only_NIST_v2.19.recode.vcf > indelCalls_FalsePositives.bed #find the snps that do not intersect


intersectBed -v -b indelCalls_NistRegions.bed -a /data3/GenomeInABottle/allIndels_only_NIST_v2.19.recode.vcf > indelCalls_NistCalls_FalseNegatives.vcf #find the nist indels that do not intersect

