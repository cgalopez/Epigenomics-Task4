# Epigenomics tasks

## Task 4: EN‐TEx ATAC‐seq data: downstream analyses

Run the docker:
```bash
sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```

We download metadata file for ATAC-seq experiments:

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&biosample_ontology.term_name=sigmoid+colon&biosample_ontology.term_name=stomach&assembly=GRCh38&assay_title=ATAC-seq"
```

Create folder and download bigBed files:

```bash
mkdir data/bigBed.files data/bigWig.files

grep -F ATAC-seq metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
    wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done
```

Check integrity of bigBed files:

```bash
for file_type in bigBed bigWig; do

# retrieve original MD5 hash from the metadata
    ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt
    # compute MD5 hash on the downloaded files 
    cat data/"$file_type".files/md5sum.txt |\
    while read filename original_md5sum; do 
        md5sum data/"$file_type".files/"$filename"."$file_type" |\
        awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
    done > tmp 
    mv tmp data/"$file_type".files/md5sum.txt

    # make sure there are no files for which original and computed MD5 hashes differ
    awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

Create needed folders and convert bigBed to BED files:

```bash
mkdir analyses/peaks.analysis

mkdir data/bed.files

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
    bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
```

Create files containing peaks falling into promoters for each tissue:

```bash
mkdir annotation
cd annotation
cp ../ChIPseq/annotation/gencode.v24.protein.coding.non.redundant.TSS.bed/
cd ..

cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
    bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
    cut -f7 |\
    sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".ATAC-seq.txt
done

wc -l analyses/peaks.analysis/genes.with.peaks.stomach.ATAC-seq.txt
# 15029

wc -l analyses/peaks.analysis/genes.with.peaks.sigmoid_colon.ATAC-seq.txt
#14830
```
The number of peaks falling within promoters in sigmoid colon are 15,029 and in stomach are 14,830

We proceed to find the number of peaks falling outside gene body.

After downloading annotation file and uncompressing it

```bash
cd annotation
cp ../ChIP-seq/annotation/gencode.v24.primary_assembly.annotation.gtf
cd ..
```

Select entries corresponding to protein coding genes from annotation file, remove mitochondrial genes and move from 1-based to 0 -based coordinate system:

```bash
awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.gene.body.bed
```

Create files containing peaks falling outside gene body for each tissue by using bedtools intersect with option -v, which returns entries of first file which do not overlap in the second:

```bash
cut -f-2 analyses/bigBed.peaks.ids.txt |
while read filename tissue; do 
    bedtools subtract -a data/bed.files/"$filename".bed -b annotation/gencode.v24.gene.body.bed -A |
    sort -u > analyses/peaks.analysis/peaks.outside.gene.body."$tissue".ATAC-seq.txt
done

wc -l analyses/peaks.analysis/peaks.outside.gene.body.stomach.ATAC-seq.txt
#34537

wc -l analyses/peaks.analysis/peaks.outside.gene.body.sigmoid_colon.ATAC-seq.txt
#37035
```
The number of peaks falling outside gene body in sigmoid colon are 37,035 and in stomach are 34,537
