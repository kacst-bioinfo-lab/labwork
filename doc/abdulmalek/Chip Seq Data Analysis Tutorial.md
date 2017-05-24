 CRITICAL NOTE *We provide all code as Unix shell instructions that generally run on Unix/Linux distributions, allow line-by-line processing of large files and are typically very robust*

**You can download bash file from https://github.com/kacst-bioinfo-lab/labwork/blob/master/code/abdulmalek/ChipBash.sh
*We always suggest keeping large files in a compressed format (e.g., using gzip).*
## Pre-Processing:
**Sub-setting Large fastq.gz files with command line**
For sub-setting large files you can use the command line
*Step 0 : Subsetting Large Zipped Fastq files*

`for sample in chip_dmel input_dmel; do`
 `zcat ${sample}.fastq.gz | sed -n 1,4000000p | gzip > ${sample}_sub.fastq.gz`
`done`
The -n 1,4000000 will keep the fastq format intact; hence a single read in fastq file has 4 attributes

### Quality control with Fastx
*Assess the read quality by the average quality score. For raw sequencing files in fastq format, Fastx is used to do quality control.*
we can do this in terminal:

**Step 1: Sequence quality check**
`for sample in chip_dmel input_dmel; do
`    
### FASTX Statistics
    `fastx_quality_stats -i <(gunzip -c ${sample}_sub.fastq.gz) -o quality/${sample}_sub_stats.txt`
   
 ### FASTX quality score
       `fastq_quality_boxplot_graph.sh -i quality/${sample}_sub_stats.txt -o quality/${sample}_quality.png -t ${sample} `  
### FASTX nucleotide distribution
`    fastx_nucleotide_distribution_graph.sh -i quality/${sample}_sub_stats.txt -o quality/${sample}_nuc.png -t ${sample}
`    
#Remove intermediate file
`    rm quality/${sample}_sub_stats.txt
`
`done
`
### Step 2: Raw reads count
*Count the number of total and unique reads. In addition, it is worthwhile to check the number and
identity of the most abundant sequences, which might identify a high amount of linkers or other contaminants.*
`for sample in chip_dmel input_dmel; do
`    `echo -en $sample"\t"`
    # Number of unique reads and most repeated read
`    gunzip -c ${sample}_sub.fastq.gz | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max){max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'
`
`done`


###  Step 3: Read length check
`for sample in chip_dmel input_dmel; do
`    `echo -en $sample"\t"`
    # Read length
`    gunzip -c ${sample}_sub.fastq.gz | awk '((NR-2)%4==0){count[length($1)]++}END{for(len in count){print len}}'
`    # Truncate longer reads to 36bp (if necessary)
`    LEN=36
`    `gunzip -c ${sample}_sub.fastq.gz | awk -v LEN=$LEN '{if((NR-2)%2==0){print substr($1,1,LEN)}else{print $0}}' | gzip > ${sample}_36bp.fastq.gz`
`done`

### Read mapping
*Read mapping. Map reads uniquely to the reference genome allowing for mismatches. We recommend using the SAM output format and then convert the files to sorted BAM files (compressed binary version) and associated index (BAI) files using SAMTools.  Convert BAM files to BED files using BEDTools.*
 CRITICAL STEP *Reads from all compared samples should be mapped with the same settings in order to avoid bias in downstream analyses such as peak calling.*

`for sample in chip_dmel input_dmel; do
`        
`gunzip -c ${sample}_36bp.fastq.gz | bowtie2 -x bowtie_index/dm3 - > ${sample}.sam`

`done`

`for sample in chip_dmel input_dmel; do
`    
    # Convert file from SAM to BAM format
    `samtools view -b -S ${sample}.sam > ${sample}_nonSorted.bam`
     # Sort BAM file
     `samtools sort -m 50000000 ${sample}_nonSorted.bam ${sample}`
      # Create index file (BAI)
      `samtools index ${sample}.bam`
        # Remove intermediate files
         `rm ${sample}.sam ${sample}_nonSorted.bam`
`done`

### Step 5: Mapped reads counts 
*Count the number of mapped reads, unique read coordinates and the maximum of reads mapped to the same genomic position.*
`for sample in chip_dmel input_dmel; do
`    
    `echo -en $sample"\t"`
    # Number of raw reads
    `raw=$(samtools view ${sample}.bam | wc -l)`
    # Number of raw, unique and most repeated reads
`    bamToBed -i ${sample}.bam | awk -v RAW=$raw '{coordinates=$1":"$2"-"$3;total++;count[coordinates]++}END{for(coordinates in count){if(!max||count[coordinates]>max){max=count[coordinates];maxCoor= coordinates};if(count[coordinates]==1){unique++}};print RAW,total,total*100/RAW,unique,unique*100/total,maxCoor,count[maxCoor],count[maxCoor]*100/total}'
`    
     # Total and top 10 of non-mapped reads
`    samtools view -f 0x0004 ${sample}.bam | awk '{read=$10;total++;count[read]++}END{print "Total_non-mapped_reads",total;for(read in count){print read,count[read]+0}}' | sort -k2,2nr | head -11
`
`done`

### Step 6: Peak calling. 
*For each immunoprecipitation sample and its corresponding input control sample, call peaks using MACS, with a stringent FDR threshold (e.g., FDR ≤ 1%) to identify confident peaks and with the default P value (10  − 5 ) to identify regions with nonrandom enrichments.*
`for pair in chip_dmel-input_dmel; do`
    `echo -en $pair"\t"`
    `chip=$(echo $pair | sed 's/-.*//')`
   ` input=$(echo $pair | sed 's/.*-//')`
    `macs2 callpeak -t ${chip}.bam  -c ${input}.bam -n short -g dm -p 1e-2 --nomodel --extsize 36 --keep-dup all`
`done`

### Step 7: sort the narrowpeak file and convert to bed
*This step will facilitate the count table step,hence the counting would use the sorted bam file it is recommended to sort the peaks to speed up the counting.* 
`for sample in chip_dmel input_dmel; do` 
    `cat short_peaks.narrowPeak |  sort -k1,1 -k2,2n > narrowPeak.bed`
`done`

### Step 8 : Count Table
*will count how many reads are in those peaks using bedtools multicov*
**Count by bedtools**
`for sample in chip_dmel input_dmel;do` 
    `bedtools multicov -bams ${sample}.bam -bed narrowPeak.bed > counts_table.bed`
`done`
 CRITICAL NOTE * When dealing with merged bed file make sure to add the peak id to the fourth column . Here we do not need change the bed file as the peak id is the fourth column.*
To do so for merged file 
`cat merge.bed | awk '{$3=$3"\t""peak_"NR}1' OFS="\t" > bed_for_multicov.bed`
But we will skip this step

# Enter R
`R`

**You can download dmel_RChipseeker.R from https://github.com/kacst-bioinfo-lab/labwork/blob/master/code/abdulmalek/dmel_RChipseeker.R
**Then type** 
`source ("dmel_RChipseeker.R")`
