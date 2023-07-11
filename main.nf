

process FASTQC {
	tag "FASTQC on $name using $task.cpus CPUs and $task.memory memory"

	input:
	tuple val(name), val(reads)

	output:
	path '*'

	script:
	"""
	fastqc ${reads}R1.fastq.gz ${reads}R2.fastq.gz -o ./
	"""
}

process CUTADAPT{
	tag "CUTADAPT on $name using $task.cpus CPUs and $task.memory memory" 
	publishDir "${params.outdir}/trimmed/${name}/", mode:'copy'
		
	input:
	tuple val(name), val(reads)
	
	output:
	tuple val("trimmed_${name}"), path("trimmed_*.fastq.gz")

	script:
	"""
	cutadapt -m 10 -q 20 -j 8 -o trimmed_${name}_R1.fastq.gz -p trimmed${name}_R2.fastq.gz ${reads}R1.fastq.gz ${reads}R2.fastq.gz
	"""
}


process BWA_ALIGN {
	tag "BWA_ALIGN on $name using $task.cpus CPUs and $task.memory memory"

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t ${task.cpus} ${params.refindex} ${reads}* > ${name}.sam
	samtools view -Sb ${name}.sam -o ${name}.bam
	"""
}


process SORT_INDEX {
  tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
  publishDir "${params.outdir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(bam)

	output:
	tuple val(name), path("${name}.sorted.bam"), path("${name}.sorted.bai")


	script:
	"""
	samtools sort ${name}.bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}

process DEDUPLICATE {
	tag "DEDUPLICATE on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(sortedBam), path(sortedBai)
	
	output:
		tuple val(name), path("${name}.deduplicated.bam"), path(sortedBai)
  path "*.markdup.txt"

	script:
	"""
	picard MarkDuplicates I=${sortedBam} O=${name}.markdup.bam M=${name}.markdup.txt
	picard BuildBamIndex INPUT=${name}.markdup.bam
 samtools view -b -F 0x400 -o ${name}.deduplicated.bam ${name}.markdup.bam
	"""
}

process QC_STATS {
	tag "QC_STATS on $name using $task.cpus CPUs and $task.memory memory"

	input:
	tuple val(name), path(bam), path(bai)

	output:
	path ("QC_results/*")
	
	script:
	"""
	samtools view -H $bam
	qualimap bamqc -bam $bam -gff ${params.covbed} -outdir ${name} -outfile ${name}.qualimap -outformat HTML
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard BedToIntervalList -I ${params.covbed} -O ${name}.interval_list -SD ${params.ref}/GRCh38-p10.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}/GRCh38-p10.fa -O ${name}.aln_metrics
	rm -r ?
	"""
}


process MULTIQC {
	tag "MultiQC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/multiqc/", mode:'copy'

	input:
	path('*')
		// tuple val(name), path("*")
		// path(markdupTxt)

	output:
	path "*"

	script:
	"""

	multiqc . -n ${name}MultiQC.html
	"""

}


workflow {

samplesList = channel.fromList(params.samples)
samplesList.view()

 fastqced =	FASTQC(samplesList)
	cutAdapted = CUTADAPT(samplesList)
	bam = BWA_ALIGN(cutAdapted)
 sortedBam =	SORT_INDEX(bam)
	deduplicatedBam = DEDUPLICATE(sortedBam)
 stats =	QC_STATS(sortedBam)
 MULTIQC(stats[0].mix(fastqced).mix(deduplicatedBam[1]).collect())
}