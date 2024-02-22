process FASTQC {
	tag "FASTQC on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "s_mem"
	label "s_cpu"


	input:
	val(sample)

	output:
	path '*'

	script:
	"""
	fastqc ${sample.name}_trim_R2.fastq.gz ${params.rawDataDir}/${sample.run}/raw_fastq/${sample.name}_R1.fastq.gz ${params.rawDataDir}/${sample.run}/raw_fastq/${sample.name}_R2.fastq.gz-o ./
	"""
}

process CUTADAPT{
	tag "CUTADAPT on $sample.name using $task.cpus CPUs and $task.memory memory" 
	// publishDir "${params.outdir}/trimmedFQ/", mode:'copy'
	label "m_mem"
	label "s_cpu"
		
	input:
	val(sample)
	
	output:
	tuple val(sample.name), path("*.fastq.gz"), val(sample.type)

	script:
	"""
	echo $sample
	cutadapt -m 10 -q 20 -j 8 -o ${sample.name}_trim_R1.fastq.gz -p ${sample.name}_trim_R2.fastq.gz ${params.rawDataDir}/${sample.run}/raw_fastq/${sample.name}_R1.fastq.gz ${params.rawDataDir}/${sample.run}/raw_fastq/${sample.name}_R2.fastq.gz
	"""
}


process BWA_ALIGN {
	tag "BWA_ALIGN on $name using $task.cpus CPUs and $task.memory memory"
	label "m_mem"
	label "l_cpu"

	input:
	tuple val(name), path(reads), val(type)

	output:
	tuple val(name), path("${name}.bam"), val(type)

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t ${task.cpus} ${params.refindex} ${reads}* | samtools view -Sb -o ${name}.bam
	"""
}


process SORT_INDEX {
 tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
  // publishDir "${params.outdir}/mapped/${name}", mode:'copy'
	label "s_mem"
	label "s_cpu"

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
	// publishDir "${params.outdir}/mapped/${name}", mode:'copy'
	label "m_mem"
	label "s_cpu"

	input:
	tuple val(name), path(sortedBam), path(sortedBai), val(type)
	
	output:
		tuple val(name), path("${name}.deduplicated.bam"), path("${name}.deduplicated.bai"), val(type)

	script:
	"""
	#picard UmiAwareMarkDuplicatesWithMateCigar I=${sortedBam} O=${name}.markdup.bam M=${name}.markdup.txt UMI_METRICS=${name}_umi_metrics.txt
	#picard MarkDuplicates I=${sortedBam} O=${name}.markdup.bam M=${name}.markdup.txt
	umi_tools dedup --umi-separator=_ --output-stats=${name}_dedup_stats.txt --stdin=${sortedBam} --stdout=${name}.deduplicated.bam
	samtools index ${name}.deduplicated.bam ${name}.deduplicated.bai
	"""
}


process CNVKIT_COVERAGE {
	tag "CNVKIT_COVERAGE on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/cnvkit/${type}/${name}", mode:'copy'
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), path(bam), path(bai), val(type)

	output:
	tuple val(name), path("${name}.targetcoverage.cnn"), path("${name}.antitargetcoverage.cnn"), val(type)
	
	script:
	"""
	cnvkit.py coverage ${bam} ${params.targetBedGeneNames} -o ${name}.targetcoverage.cnn
	cnvkit.py coverage ${bam} ${params.antitargetBed} -o ${name}.antitargetcoverage.cnn
	"""
}

process CNVKIT_REFERENCE {
	tag "CNVKIT_REFERENCE using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/cnvkit/reference", mode:'copy'
	label "s_mem"
	label "s_cpu"

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
	publishDir "${params.outdir}/cnvkit/patients/${name}", mode:'copy'
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), path(targetCov), path(antitargetCov), path(reference), path(vcf)

	output:
 path "*"	
	tuple val(name), path("${name}.cns")

	script:
	"""
 cnvkit.py fix ${targetCov} ${antitargetCov} ${reference} -o ${name}_WeirdChr.cnr
 cat ${name}_WeirdChr.cnr | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") {\$1=\$1;print \$0}' > ${name}.cnr
	
	#hmm-germline
	cnvkit.py segment ${name}.cnr -m cbs -v $vcf -o ${name}.cns 
	#cnvkit.py segmetrics -s ${name}.cn{s,r} --ci -o ${name}.segmetrics.cns
#	cnvkit.py call ${name}.segmetrics.cns --filter ci -m clonal -v $vcf -o ${name}.calls.segmetrics.cns

		cnvkit.py call ${name}.cns -v $vcf -m none --purity 1 -o ${name}.calls.segmetrics.cns


 #cat ${name}_WeirdChr.cnr | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") {\$1=\$1;print \$0}' > ${name}.cnr
 #cat ${name}.calls.segmetrics.cns | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") {\$1=\$1;print \$0}' > ${name}.cns

 cnvkit.py scatter ${name}.cnr -s ${name}.calls.segmetrics.cns -v $vcf --y-max=3 --y-min=-3 -o ${name}-scatter.pdf
 #cnvkit.py scatter ${name}.cnr -s ${name}.cns -v $vcf --y-max=3 --y-min=-3 -o ${name}-scatter-cbs-filtered.pdf
 cnvkit.py diagram ${name}.cnr -s ${name}.cns -o ${name}-diagram.pdf
	"""
}

process CNVKIT_PROCESS_TABLE{
	tag "CNVKIT_PROCESS_TABLE on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/cnvkit/patients/${name}", mode:'copy'
	label "s_mem"
	label "s_cpu"
	
	input:
	tuple val(name), path(cns)

	output:
 path "${name}.souhrn.txt"	

	script:
	"""
	head -n 1 ${cns} | tr -d '\\n' > ${name}.souhrn.txt
	bedtools map -c 4 -a <(bedtools sort -i <( tail -n +2 ${cns})) -b $params.GrCh38cytomap -o concat > ${name}Cytomap.txt
 bedtools map -c 4,4,5,6,7,8 -a ${name}Cytomap.txt -b $params.GrCh38CNV -o collapse,count,collapse,sum,sum,sum -g ${params.genomeLens}> souhrn.txt
 echo -e "\\tCytoCords\\tPubMedID\\tCitations_count\\tCnv_type\\tsamplesizeCount\\tobservedgainsCount\\tobservedlossesCount">> ${name}.souhrn.txt
	cat souhrn.txt >> ${name}.souhrn.txt
	"""
}

process QC_STATS {
	tag "QC_STATS on $name using $task.cpus CPUs and $task.memory memory"
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), path(bam), path(bai)

	output:
	path "*"
	
	script:
	"""
	samtools view -H $bam
	qualimap bamqc -bam $bam -gff ${params.covbed} --java-mem-size=6G -outdir ${name} -outfile ${name}.qualimap -outformat HTML
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard BedToIntervalList -I ${params.covbed} -O ${name}.interval_list -SD ${params.ref}/GRCh38-p10.dict
	picard CollectHsMetrics -I $bam -BAIT_INTERVALS ${name}.interval_list -TARGET_INTERVALS ${name}.interval_list -R ${params.ref}/GRCh38-p10.fa -O ${name}.aln_metrics
	"""
}


process HAPLOCALLER {
	tag "HaplotypeCaller on $name  using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/HaplotypeCaller/", mode:'copy'
	container "broadinstitute/gatk:4.1.3.0"
	label "s_mem"
	label "s_cpu"

	input:
	tuple val(name), path(bam), path(bai)

	output:
	tuple val(name), path("${name}.g.vcf.gz")

	script:
	"""
	gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ${params.ref}/GRCh38-p10.fa \
   -I $bam \
   -O ${name}.g.vcf.gz \
   -ERC GVCF
	"""

}

process MULTIQC {
	tag "MultiQC on all samples using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outdir}/multiqc/", mode:'copy'
	label "s_mem"
	label "s_cpu"

	input:
	path('*')

	output:
	path "*"

	script:
	"""
	multiqc . -n MultiQC.html
	"""

}


workflow {
 samplesList = channel.fromList(params.samples)
	cutAdapted = CUTADAPT(samplesList)
 fastqced =	FASTQC(cutAdapted)
	bam = BWA_ALIGN(cutAdapted)
 sortedBam =	SORT_INDEX(bam)
	deduplicatedBam = DEDUPLICATE(sortedBam)

	 deduplicatedBam.branch {  //this is for HAPLOCALLER process, to only call vars on tumor samples for later SNP allele frequencies
					Normal: it[3] == "normal"
					 return [it[0],it[1],it[2]] //without type
					Tumor: it[3] == "proband"
					 return [it[0],it[1],it[2]] //without type
 	}.set{dedupBam}

	coverage = CNVKIT_COVERAGE(deduplicatedBam[0])
 coverage.branch {
					Normal: it[3] == "normal"
					 return [it[0],it[1],it[2]] //without type
					Tumor: it[3] == "proband"
					 return [it[0],it[1],it[2]] //without type
 	}.set{coverage}

	 NormalOnlyCoveragePaths = coverage.Normal.map({return [it[1],it[2]]
	 })

 cnvkitReference = CNVKIT_REFERENCE(NormalOnlyCoveragePaths.collect())
 tumorWithReference = coverage.Tumor.combine(cnvkitReference)
	TumorVCFs = HAPLOCALLER(dedupBam.Tumor)
	CoverageVcfTumor = tumorWithReference.join(TumorVCFs)//.view()
 cnvOutput =	CNVKIT_TUMOR(CoverageVcfTumor)
 CNVKIT_PROCESS_TABLE(cnvOutput[1])
 stats =	QC_STATS(deduplicatedBam[0])
	mqcChannel = stats.mix(fastqced).collect()
  MULTIQC(mqcChannel)
}
