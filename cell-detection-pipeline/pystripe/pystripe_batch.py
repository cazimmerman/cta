import os,sys,subprocess

cwd = os.getcwd()
utils_fld = os.path.join(cwd,"utils")
sys.path.append(utils_fld)
from pipeline_utils import fast_scandir

fpath = os.path.join('/jukebox/YOUR_DATA',sys.argv[1],'rawdata/resolution_3.6x')
slurmpath = '/jukebox/YOUR_DIR/lightsheet-pipeline/pystripe/slurm_scripts/pystripe.sh'
slurmpath_flat = '/jukebox/YOUR_DIR/lightsheet-pipeline/pystripe/slurm_scripts/pystripe_flat.sh'

try:
    input_dir = fast_scandir(os.path.join(fpath,'Ex_488_Em_525_stitched'))[-1]
except:
    input_dir = os.path.join(fpath,'Ex_488_Em_525_stitched')
flat = os.path.join(fpath,'Ex_488_Em_525_stitched','flat.tiff')
output_dir = os.path.join(fpath,'Ex_488_Em_525_corrected')
if os.path.exists(flat):
    run_job = 'sbatch --export=ALL,input_dir=\'' + input_dir + '\',flat=\'' + flat + '\',corrected_dir=\'' + output_dir + '\' ' + slurmpath_flat
else:
    run_job = 'sbatch --export=ALL,input_dir=\'' + input_dir + '\',corrected_dir=\'' + output_dir + '\' ' + slurmpath
jobid = subprocess.check_output(run_job, shell=True)

try:
#    input_dir = fast_scandir(os.path.join(fpath,'Ex_647_Em_690_stitched'))[-1]
    input_dir = fast_scandir(os.path.join(fpath,'Ex_639_Em_690_stitched'))[-1]
except:
#    input_dir = os.path.join(fpath,'Ex_647_Em_690_stitched')
    input_dir = os.path.join(fpath,'Ex_639_Em_690_stitched')
flat = os.path.join(fpath,'Ex_639_Em_690_stitched','flat.tiff')
output_dir = os.path.join(fpath,'Ex_647_Em_690_corrected')
if os.path.exists(flat):
    run_job = 'sbatch --export=ALL,input_dir=\'' + input_dir + '\',flat=\'' + flat + '\',corrected_dir=\'' + output_dir + '\' ' + slurmpath_flat
else:
    run_job = 'sbatch --export=ALL,input_dir=\'' + input_dir + '\',corrected_dir=\'' + output_dir + '\' ' + slurmpath
jobid = subprocess.check_output(run_job, shell=True)