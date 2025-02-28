# General imports
import os,sys,glob,shutil,json
cwd = os.getcwd()
utils_fld = os.path.join(cwd,"utils")
sys.path.append(utils_fld)
from pipeline_utils import fast_scandir
# ClearMap2 imports
sys.path.append('/jukebox/YOUR_DIR/ClearMap2-master')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ParallelProcessing.BlockProcessing as bp

# Select datasets to analyze
if __name__ == "__main__":
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,imaging_request,'rawdata/resolution_3.6x')

	# Initialize ClearMap2 workspace object
	# directory = os.path.join(output_path,fpath)
	expression_raw = '/Ex_642_Em_2_corrected/Z<Z,4>.tif'
	ws = wsp.Workspace('CellMap',directory=dst_dir)
	ws.update(raw=expression_raw)
	ws.debug = False
	print()
	ws.info()

	# Finally, figure out how many blocks there so we can know how many array jobs to make
	print("Determining number of blocks")
	sys.stdout.flush()
	blocks = bp.split_into_blocks(ws.source('raw'),
                    processes='serial',
                    axes=[2],
                    size_min=10,
                    size_max=30,
                    overlap=5,
                    verbose=True)

	n_blocks = len(blocks)
	block_json_file = os.path.join(dst_dir,'block_processing_info.json')
	block_data = {
		'n_blocks':n_blocks
		}
	with open(block_json_file,'w') as outfile:
		json.dump(block_data,outfile)
	print(f"Wrote block info to file: {block_json_file}")
