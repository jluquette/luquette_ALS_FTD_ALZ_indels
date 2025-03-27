## 01_run_cellbender_on_our_samples.sh
#!/bin/bash
#SBATCH --partition=bch-gpu               
#SBATCH --time=24:00:00                   
#SBATCH --job-name=run_cellbender_oursample 
#SBATCH --output=log_12/output_%A_%a.txt  
#SBATCH --nodes=1                         
#SBATCH --gres=gpu:1                      
#SBATCH --mem=30G                         
#SBATCH --array=0-11%1                    

# Load specific CUDA version
module load cuda/11.7
echo "Loaded CUDA version 11.7"

# Set environment variables for CUDA
export PATH=/usr/local/cuda-11.7/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.7/lib64:$LD_LIBRARY_PATH

# Activate conda environment and set library path
conda_base=$(conda info --base)
source $conda_base/etc/profile.d/conda.sh
conda activate cellBender_new
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Set TMPDIR
export TMPDIR=/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/02_cellbender_our_group/tmpfiles

# Print CUDA version to log
cuda_version=$(module list 2>&1 | grep cuda)
echo "Using CUDA version: $cuda_version"

# Check GPU status
nvidia-smi

# Activate Python environment and run Python code to check its tmp_dir and CUDA availability
python << EOF
import torch
print("Torch CUDA available:", torch.cuda.is_available())
import tempfile
# Create a temporary directory
with tempfile.TemporaryDirectory() as tmp_dir:
    print("This is new Temporary directory:", tmp_dir)
EOF

# Setup directories
cellranger_output_dir="/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/01_cell_ranger_output_our_group"
cellbender_output_dir="/lab-share/Gene-Lee-ANR-e2/KKS/01_Ambient_RNA/out/02_cellbender_our_group"

# Dynamically assign input and output directories based on the SLURM array task ID
srr_dirs=($(ls $cellranger_output_dir | grep "^\(.*-normal\|.*-TDP43\)$"))

# Check the number of srr_dirs
if [ ${#srr_dirs[@]} -ne 12 ]; then
    echo "Error: srr_dirs was set incorrectly"
    exit 1
fi

srr_dir="${srr_dirs[$SLURM_ARRAY_TASK_ID]}"
input_directory="${cellranger_output_dir}/${srr_dir}/outs/raw_feature_bc_matrix.h5"
output_directory="${cellbender_output_dir}/${srr_dir}"

# Check if input file exists
if [ ! -f "$input_directory" ]; then
    echo "Error: Input file for sample ${srr_dir} does not exist."
    echo "Skipping sample ${srr_dir}."
else
    # If the file exists, create output directory
    mkdir -p $output_directory
    # Start processing
    echo "Start processing sample ${srr_dir}"
    echo $(date "+%Y-%m-%d %H:%M:%S")
    # Code to check its tmp_dir
    echo "This is new Temporary directory: ${TMPDIR} for ${srr_dir}"
    ## Run cellbender
    cellbender remove-background \
        --cuda \
        --input $input_directory \
        --output ${output_directory}/output.h5
    echo "End processing sample ${srr_dir}"
    echo $(date "+%Y-%m-%d %H:%M:%S")
fi


