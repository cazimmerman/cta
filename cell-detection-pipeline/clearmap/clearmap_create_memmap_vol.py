# General imports
import os,sys,glob,shutil,json, pickle
import numpy as np
from PIL import Image
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

# ClearMap2 imports
sys.path.append('/jukebox/YOUR_DIR/ClearMap2-master')
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ParallelProcessing.BlockProcessing as bp
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Utils.HierarchicalDict as hdict

def fast_scandir(dirname):
	""" gets all folders recursively """
	subfolders= [f.path for f in os.scandir(dirname) if f.is_dir()]
	for dirname in list(subfolders):
		subfolders.extend(fast_scandir(dirname))
	return subfolders

# Select datasets to analyze
if __name__ == "__main__":
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,imaging_request,'rawdata/resolution_3.6x')
	os.makedirs(dst_dir,exist_ok=True)
	print(f"Creating dst dir: {dst_dir}")

	# Initialize ClearMap2 workspace object
	# directory = os.path.join(output_path,fpath)
	expression_raw = '/Ex_642_Em_2_corrected/Z<Z,4>.tif'
	ws = wsp.Workspace('CellMap',directory=dst_dir)
	ws.update(raw=expression_raw)
	ws.debug = False
	print()
	ws.info()

	# Create combined npy volume -- this is CPU heavy
	source = ws.source('raw')
	sink = ws.filename('stitched')
	sys.stdout.flush()
	if not os.path.exists(sink):
		print("Converting raw files into single stitched.npy file")
		sys.stdout.flush()
		io.convert(source,sink,verbose=True)
		print("Successfully created stitched.npy file")
	else:
		print("Stitched.npy file already exists")

	print()
	print("Determining cell detection parameters...")
	print()
	xyz = np.around(np.shape(ws.source('stitched')))/2
	xyz = xyz.astype(int)
	ws.debug = 'background'
	slicing = (slice(xyz[0]-1000,xyz[0]+1000),slice(xyz[1]-1000,xyz[1]+1000),slice(xyz[2]-100,xyz[2]+100))
	ws.create_debug('stitched',slicing=slicing)
	im = Image.fromarray(np.amax(ws.source('stitched'),2))
	im.save(dst_dir+'/background_correction_volume_raw_MIP.tif')
	clearmap_params_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_parameter.p'
	with open(clearmap_params_file,'rb') as f:
	    cell_detection_parameter = pickle.load(f)
	print("median filtering...")
	data = ndimage.median_filter(ws.source('stitched'),3)
	print("background subtracting...")
	data = cells.remove_background(data,cell_detection_parameter['background_correction']['shape'],cell_detection_parameter['background_correction']['form'])
	im = Image.fromarray(np.amax(data,2))
	im.save(dst_dir+'/background_correction_volume_subtracted_MIP.tif')
	print("sigma clipping...")
	data = np.reshape(data,-1)
	mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=10, cenfunc='median', stdfunc='mad_std')
	print("finished!!!")
	print(" ")
	print("(mean, median, std)")
	print((mean, median, std))
	print(" ")
	os.remove(dst_dir+'/background_stitched.npy')
	cell_detection_parameter['shape_detection']['threshold'] = int(10*mean)
	clearmap_params_file = dst_dir+'/cell_detection_parameter.p'
	with open(clearmap_params_file, 'wb') as f:
	    pickle.dump(cell_detection_parameter,f)
	cell_detection_parameter = None
	with open(clearmap_params_file,'rb') as f:
	    cell_detection_parameter = pickle.load(f)
	print("path                    : " + clearmap_params_file)
	hdict.pprint(cell_detection_parameter)

	sys.stdout.flush()
