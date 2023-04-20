#!/usr/bin/bash

# Analyzing data for Intestine FGF15 FACS sorted

# Run FASTQC
mkdir FASTQC
SAMPLES=`ls -R -d seeley/* | grep -E Sample | sed -e 's/seeley\///g'`
for f in $SAMPLES; do
  cd seeley/${f}
  mkdir ../../FASTQC/${f}
  echo 'FASTQC on' $f
  zcat *.fastq.gz | fastqc -t 8 -o ../../FASTQC/$f stdin
  cd ../../
done

# Remove low quality reads (Phred 20)
for f in $SAMPLES; do
  date
  cd seeley/${f}
  echo 'Filtering' ${f}
  zcat *.fastq.gz | fastq_quality_filter -q 20 -z -o ${f}_filtered.fastq.gz
  cd ../../
done

# Map filtered reads to STAR
mkdir STAR_outs
for f in $SAMPLES; do
  date
  cd seeley/${f}
  echo 'Aligning' ${f}
  ~/STAR/bin/Linux_x86_64/STAR \
    --runThreadN 8 \
    --genomeDir ~/STAR/genome_mm_CreERT2_GfpL10a \
    --sjdbGTFfile ~/STAR/mm_CreERT2_GfpL10a.gtf \
    --readFilesIn ${f}_filtered.fastq.gz \
    --readFilesCommand zcat \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic
  cd ../../
done
