'''
RNA seq variant calling pipeline. 
May 30 2018 Jongsoo Yoon 
https://software.broadinstitute.org/gatk/documentation/article.php?id=3891

Skipped step 4 (indel realignment) and step 5 (base recalibration) as 
we're only interested in RNA variant calling for DNA mutation validation purposes

Assumes one has STAR aligned RNA seq BAM file ready.
Will process STAR bam 
'''

import sys
import subprocess
import shlex
import argparse
import re
import os
import pysam

JAVA = '/home/users/cjyoon/tools/jdk1.8.0_171/bin/java' # need java version 1.8 to work, otherwise complains of malformed walker argument
GATK = '/home/users/tools/gatk/gatk-3.5/GenomeAnalysisTK.jar'
PICARD = '/home/users/tools/picard/dist/picard.jar'
REFERENCE = '/home/users/cjyoon/reference/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/ref_genome.fa' # STAR alignment reference (contains 'chr' in front of chromosome number
def argument_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_star_bam', required=True)
	args = vars(parser.parse_args())
	rna_bam = args['input_star_bam']
	return rna_bam

def is_bam_sorted(bamfile):
	'''checks if a bam file is sorted.'''
	bam = pysam.AlignmentFile(bamfile)
	try:
		if bam.header['HD']['SO']== 'coordinate':
			return True
		else:
			return False
	except KeyError:
		return False

def cleanup(list_of_files):
	'''delete list of files, usually intermediate files that aren't necessary for keeping'''
	for afile in list_of_files:
		if os.path.isfile(afile):
			cmd = 'rm -rf ' + afile
			subprocess.call(shlex.split(cmd))

	return 0


def main():
	rna_bam = argument_parser()

	# sort bam
	if is_bam_sorted(rna_bam):
		# if bam is already sorted, then skip
		print(f'{rna_bam} is already sorted')
		sorted_rna_bam = rna_bam
	else:
		sorted_rna_bam = re.sub(r'.bam$', '.sorted.bam', rna_bam)
		# sorting = f'{JAVA} -jar {PICARD} AddOrReplaceReadGroups I={rna_bam} O={sorted_rna_bam} SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample'
		sorting = f'samtools sort -o {sorted_rna_bam} {rna_bam}'
		print(sorting)
		subprocess.call(shlex.split(sorting))

	# replace header so that GATK can properly process the bam file
	reheadered_rna_bam = re.sub(r'.bam$', '.rg.bam', sorted_rna_bam)
	replace = f"samtools view -H {sorted_rna_bam} | sed 's/^@RG.*/@RG\tID:GRPundef\tSM:sample\tLB:library\tPL:platform/g'| samtools reheader - {sorted_rna_bam} > {reheadered_rna_bam}"
	print(replace)
	os.system(replace)
	

	# mark duplicate
	marked_bam = re.sub(r'.bam$', '.md.bam', reheadered_rna_bam)
	marking = f'java -jar {PICARD} MarkDuplicates I={reheadered_rna_bam} O={marked_bam} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics'
	print(marking)
	subprocess.call(shlex.split(marking))

	# SplitNCigarReads
	split_rna_bam = re.sub(r'.bam$', '.split.bam', marked_bam)
	splitNcigaring = f'{JAVA} -jar {GATK} -T SplitNCigarReads -R {REFERENCE} -I {marked_bam} -o {split_rna_bam} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
	print(splitNcigaring)
	subprocess.call(shlex.split(splitNcigaring))

	# RNA variant calling
	variantVCF = re.sub(r'.bam$', '.bam.vcf', split_rna_bam)
	variantCalling = f'{JAVA} -jar {GATK} -T HaplotypeCaller -R {REFERENCE} -I {split_rna_bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -o {variantVCF}'
	print(variantCalling)
	subprocess.call(shlex.split(variantCalling))

	print('RNA variant calling output ' + variantVCF)
	cleanup([sorted_rna_bam, marked_bam])
	print('Done cleaning up intermediate files')	
	return 0


if __name__=='__main__':
	main()
