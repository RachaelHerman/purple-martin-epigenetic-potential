#!/bin/bash
# Usage: ./submit_merge_bams.sh samples.txt

SAMPLE_LIST=$1

while read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=merge_${SAMPLE}
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
#SBATCH --mem=512G
#SBATCH --time=48:00:00

module load samtools

export JAVA_HOME=/path/to/openjdk17
export PATH=\$JAVA_HOME/bin:\$PATH

./merge_bams_script.sh $SAMPLE
EOT

done < "$SAMPLE_LIST"
