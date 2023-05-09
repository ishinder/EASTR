wd=$(pwd -P)
hisat_index=$wd/ref/hisat2/hisat_index
NCPU=$1
TEMPLOC="$wd"/tmp
ALIGNLOC="$wd"/BAM/original

mkdir -p "$TEMPLOC" "$ALIGNLOC"

for i in $(ls "$wd"/fastq/*_1.fastq); do
    name=$(basename "$i" | cut -d'_' -f1)
    echo "$name"
    if [ ! -f "$ALIGNLOC"/"${name}".bam ]; then
        hisat2 --rf --no-unal -p $NCPU \
                -x "$hisat_index" \
                -1 "$wd"/fastq/"${name}"_1.fastq \
                -2 "$wd"/fastq/"${name}"_2.fastq \
                -S "$TEMPLOC"/"${name}".sam
        
        samtools view -bS "$TEMPLOC"/"${name}".sam | \
        samtools sort -@ $NCPU - -o "$ALIGNLOC"/"${name}"_hisat.bam
        rm "$TEMPLOC"/"${name}".sam
    fi
done

rmdir tmp