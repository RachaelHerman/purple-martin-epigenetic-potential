#!/bin/bash
BASE_DIR="/path/to/puma_wgs"
PIPELINE_SCRIPT="merge_map_clean_allsets.py"
BWA_INDEX="/path/to/GCA_022316685.1_UM_PSub_1_genomic"
LOG_DIR="$BASE_DIR/slurm_logs"
mkdir -p "$LOG_DIR"
NBTHREADS=40

export PATH=/path/to/adapterremoval-2.3.4/build:$PATH

SAMPLES=()
for SAMPLE in "$BASE_DIR"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE")
    if [ -d "$SAMPLE/$SAMPLE_NAME" ]; then
        SAMPLES+=("$SAMPLE")
    fi
done

NUM_SAMPLES=${#SAMPLES[@]}
SAMPLES_FILE="$LOG_DIR/samples.txt"
printf "%s\n" "${SAMPLES[@]}" > "$SAMPLES_FILE"

sbatch --array=1-${NUM_SAMPLES}%4 <<EOT
#!/bin/bash
#SBATCH --job-name=puma_map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$NBTHREADS
#SBATCH --mem=120G
#SBATCH --time=48:00:00

module load bwa
module load python/3.13
module load samtools/1.7

SAMPLE_DIR=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_FILE")
SAMPLE_NAME=\$(basename "\$SAMPLE_DIR")

python3 $PIPELINE_SCRIPT \$SAMPLE_NAME \$SAMPLE_DIR/\$SAMPLE_NAME $NBTHREADS $BWA_INDEX
EOT
