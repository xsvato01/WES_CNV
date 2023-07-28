

process FASTQC {
	tag "FASTQC on $name using $task.cpus CPUs and $task.memory memory"

	input:
	tuple val(name), val(reads), val(type)

	output:
	path '*'

	script:
	"""
	fastqc ${reads}_R1.fastq.gz ${reads}_R2.fastq.gz -o ./
	"""
}

process CUTADAPT{
	tag "CUTADAPT on $name using $task.cpus CPUs and $task.memory memory" 
	publishDir "${params.outdir}/trimmedFQ/", mode:'copy'
		
	input:
	tuple val(name), val(reads), val(type)
	
	output:
	tuple val("${name}"), path("*.fastq.gz"), val(type)

	script:
	"""
	cutadapt -m 10 -q 20 -j 8 -o ${name}_trim_R1.fastq.gz -p ${name}_trim_R2.fastq.gz ${reads}_R1.fastq.gz ${reads}_R2.fastq.gz
	"""
}


process BWA_ALIGN {
	tag "BWA_ALIGN on $name using $task.cpus CPUs and $task.memory memory"

	input:
	tuple val(name), path(reads), val(type)

	output:
	tuple val(name), path("${name}.bam"), val(type)

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
	tuple val(name), path(bam), val(type)

	output:
	tuple val(name), path("${name}.sorted.bam"), path("${name}.sorted.bai"), val(type)


	script:
	"""
	samtools sort ${name}.bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}

process DEDUPLICATE {
	tag "DEDUPLICATE on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/mapped/${name}", mode:'copy'

	input:
	tuple val(name), path(sortedBam), path(sortedBai), val(type)
	
	output:
		tuple val(name), path("${name}.deduplicated.bam"), path("${name}.deduplicated.bai"), val(type)
  path "*.txt"

	script:
	"""
	#picard UmiAwareMarkDuplicatesWithMateCigar I=${sortedBam} O=${name}.markdup.bam M=${name}.markdup.txt UMI_METRICS=${name}_umi_metrics.txt
	picard MarkDuplicates I=${sortedBam} O=${name}.markdup.bam M=${name}.markdup.txt
	umi_tools dedup --umi-separator=_ --output-stats=${name}_dedup_stats.txt --stdin=${sortedBam} --stdout=${name}.deduplicated.bam
	samtools index ${name}.deduplicated.bam ${name}.deduplicated.bai
 #samtools view -b -F 0x400 -o ${name}.deduplicated.bam ${name}.markdup.bam
	"""
}


process CNVKIT_COVERAGE {
	tag "CNVKIT_COVERAGE on $name using $task.cpus CPUs and $task.memory memory"
		publishDir "${params.outdir}/cnvkit/${name}", mode:'copy'


	input:
	tuple val(name), path(bam), path(bai), val(type)

	output:
	tuple val(name), path("*targetcoverage.cnn"), path("*antitargetcoverage.cnn"), val(type)
	
	script:
	"""
	cnvkit.py coverage ${bam} ${params.targetBed} -o ${name}.targetcoverage.cnn
	cnvkit.py coverage ${bam} ${params.antitargetBed} -o ${name}.antitargetcoverage.cnn
	"""
}

process CNVKIT_REFERENCE {
	tag "CNVKIT_REFERENCE using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/cnvkit/reference", mode:'copy'


	input:
	path "*"

	output:
 path "Reference.cnn"	

	script:
	"""
	cnvkit.py reference *coverage.cnn -f ${params.ref}/GRCh38-p10.fa -o Reference.cnn
	"""
}

process CNVKIT_TUMOR {
	tag "CNVKIT_TUMOR on $name using $task.cpus CPUs and $task.memory memory"
		publishDir "${params.outdir}/cnvkit/${name}", mode:'copy'


	input:
	tuple val(name), path(targetCov), path(antitargetCov)
	path "reference"	

	output:
 path "*"	

	script:
	"""
 cnvkit.py fix ${targetCov} ${antitargetCov} ${reference} -o ${name}.cnr
 cnvkit.py segment ${name}.cnr -o ${name}.cns

 cnvkit.py scatter ${name}.cnr -s ${name}.cns --y-max=3 --y-min=-3 -o ${name}-scatter.png
 cnvkit.py diagram ${name}.cnr -s ${name}.cns -o ${name}-diagram.pdf
	"""
}

process QC_STATS {
	tag "QC_STATS on $name using $task.cpus CPUs and $task.memory memory"

	input:
	tuple val(name), path(bam), path(bai)

	output:
	 path "*"
	
	script:
	"""
	samtools view -H $bam
	qualimap bamqc -bam $bam -gff ${params.covbed} -outdir ${name} -outfile ${name}.qualimap -outformat HTML
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard BedToIntervalList -I ${params.covbed} -O ${name}.interval_list -SD ${params.ref}/GRCh38-p10.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}/GRCh38-p10.fa -O ${name}.aln_metrics
	"""
}


process MULTIQC {
	tag "MultiQC on all samples using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/multiqc/", mode:'copy'

	input:
	path('*')
		// tuple val(name), path("*")
		// path(markdupTxt)

	output:
	path "*"

	script:
	"""
	multiqc . -n MultiQC.html
	"""

}


workflow {

samplesList = channel.fromList(params.samples)
 fastqced =	FASTQC(samplesList)
	cutAdapted = CUTADAPT(samplesList)
	bam = BWA_ALIGN(cutAdapted)
 sortedBam =	SORT_INDEX(bam)
	deduplicatedBam = DEDUPLICATE(sortedBam)
	coverage = CNVKIT_COVERAGE(deduplicatedBam[0])

 coverage
 .branch {
					Normal: it.type == "normal"
					Tumor: it.type == "tumor"
 	}.set{sorted}

cnvkitReference = CNVKIT_REFERENCE(sorted.Normal.collect())
tumorWithReference = sorted.Tumor.combine(cnvkitReference)
tumorWithReference.view()

 // stats =	QC_STATS(deduplicatedBam[0])
	// mqcChannel = stats.mix(cutAdapted).mix(deduplicatedBam[1]).collect()
	// mqcChannel.view()
 // MULTIQC(mqcChannel)
}