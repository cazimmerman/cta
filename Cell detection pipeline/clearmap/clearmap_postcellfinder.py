# General imports
import os,sys,pickle
import numpy as np
import pandas as pd
import json
import tifffile
from concurrent.futures import ProcessPoolExecutor
from brain_atlas_toolkit import graph_tools

# ClearMap2 imports
sys.path.append('/jukebox/YOUR_DIR/ClearMap2-master')
os.chmod("/jukebox/YOUR_DIR/ClearMap2-master/ClearMap/External/elastix/build/bin/transformix",0o777)
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Settings as settings
import ClearMap.Alignment.Resampling as res
import ClearMap.Alignment.Elastix as elx

from functools import partial

def get_count_and_volume(region_idx,intensity_idx,segment_name_dict,eroded_atlas_vol,segment_list,):
	"""
	---PURPOSE---
	Given a list of segments, get the number of counts in each brain region.
	---INPUT---
	region_idx        - 1D array containing the brain region ID in which each cell was detected
	segment_name_dict - dictionary mapping ids:segment names
	eroded_atlas_vol  - the eroded atlas volume
	segment_list      - list of integer ids
	---OUTPUT---
	count_dict   - dictionary where keys are segment ids and values are counts
	volume_dict  - dictionary where keys are segment ids and values are volume of that region in the eroded atlas
	"""
	count_dict = {}
	intensity_dict = {}
	volume_dict = {}
	print(f"parallel processing {segment_list}")
	for atlas_segment_id in segment_list:
		segment_name = segment_name_dict[int(atlas_segment_id)]
		# if int(atlas_segment_id) not in atlas_segments:
		# 	print('SKIPPING ' + segment_name + ' (' + atlas_segment_id + ')')
		# 	count_dict[segment_name] = 0
		# 	intensity_dict[segment_name] = 0
		# 	volume_dict[segment_name] = 0
		# else:
		#print('RUNNING ' + segment_name + ' (' + atlas_segment_id + ')' + str(np.shape(region_idx)))
		count_dict[segment_name] = np.count_nonzero(region_idx==int(atlas_segment_id))
		intensity_dict[segment_name] = np.sum(intensity_idx[np.where(region_idx==int(atlas_segment_id))])

			######
			# volume_dict[segment_name] = np.sum(eroded_atlas_vol == int(atlas_segment_id)) * np.power(0.02,3) # PMA
		volume_dict[segment_name] = np.sum(eroded_atlas_vol == int(atlas_segment_id)) * np.power(0.025,3) # Allen CCFv3
			######

	return count_dict,intensity_dict,volume_dict

if __name__ == "__main__":
	n_cores = os.cpu_count()
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")
	atlas = sys.argv[4].strip().rstrip("/")

	if atlas == 'Princeton':
		# Princeton Mouse Atlas
		eroded_atlas_file = '/jukebox/YOUR_DIR/clearmap-output/utilities/princeton-atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels.tif'
		segment_props_file = '/jukebox/YOUR_DIR/clearmap-output/utilities/princeton-atlas/PMA_16bit_hierarch_labels_segment_properties_info'
		ontology_json_file = '/jukebox/YOUR_DIR/clearmap-output/utilities/princeton-atlas/PMA_ontology.json'
	elif atlas == 'Allen':
		# Allen Mouse Brain Atlas
		eroded_atlas_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/modified-atlas/annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz.tif'
		segment_props_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/modified-atlas/allenatlas_2017_16bit_hierarch_labels_segment_properties_info_cz'
		ontology_json_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/modified-atlas/allen_hierarch_labels_fillmissing_cz.json'
	else:
		sys.exit(f"Atlas provided: {atlas} is not accepted. Must be one of ['Princeton','Allen']")

	request_name,sample_name = sample_dir.split('/')[-2:]
	workspace_dir =  os.path.join(output_rootpath,
		request_name,sample_name,imaging_request,
		"rawdata/resolution_3.6x")

	#  Set paths to Elastix transformation files and downsized files
	elastix_inverse_dir = os.path.join(workspace_dir,'elastix_inverse_transform')
	ch488_downsized_dir = os.path.join(workspace_dir,'Ex_488_Em_0_downsized')
	ch642_downsized_dir = os.path.join(workspace_dir,'Ex_642_Em_2_downsized')
	ch488_downsized_file = os.path.join(ch488_downsized_dir,'downsized_for_atlas_ch488.tif')
	ch642_downsized_file = os.path.join(ch642_downsized_dir,'downsized_for_atlas_ch642.tif')

	# Initialize ClearMap2 workspace object
	ws = wsp.Workspace('CellMap',directory=workspace_dir)
	ws.debug = False
	print()
	ws.info()

	# Load atlas files
	eroded_atlas_vol = np.array(tifffile.imread(eroded_atlas_file)).astype('uint16')
	atlas_segments = np.unique(eroded_atlas_vol)
	atlas_segments = np.array([x for x in atlas_segments if x!=0])

	with open(segment_props_file,'r') as infile:
		segment_props_dict = json.load(infile)

	with open(ontology_json_file,'r') as infile:
		ontology_dict = json.load(infile)

	import xml.etree.ElementTree as ET

## ORIGINAL CODE __________________________________________
	# fpath = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name,"points/cell_classification.xml")
	# if os.path.exists(fpath):
	# 	tree = ET.parse(fpath)
	# 	root = tree.getroot()
	# 	n_cells = len(root[1][2])
	# 	coordinates_keep = np.empty([n_cells-1, 3], dtype=int)
	# 	for idx1 in range(1,n_cells):
	# 		onecell = root[1][2][idx1]
	# 		for idx2,child in enumerate(onecell):
	# 			coordinates_keep[idx1-1,idx2] = int(child.text)
	# else:
	# 	fpath1 = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name+'-pt1',"points/cell_classification.xml")
	# 	tree = ET.parse(fpath1)
	# 	root = tree.getroot()
	# 	n_cells = len(root[1][2])
	# 	coordinates_keep1 = np.empty([n_cells-1, 3], dtype=int)
	# 	for idx1 in range(1,n_cells):
	# 		onecell = root[1][2][idx1]
	# 		for idx2,child in enumerate(onecell):
	# 			coordinates_keep1[idx1-1,idx2] = int(child.text)
	# 	fpath2 = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name+'-pt2',"points/cell_classification.xml")
	# 	tree = ET.parse(fpath2)
	# 	root = tree.getroot()
	# 	n_cells = len(root[1][2])
	# 	coordinates_keep2 = np.empty([n_cells-1, 3], dtype=int)
	# 	for idx1 in range(1,n_cells):
	# 		onecell = root[1][2][idx1]
	# 		for idx2,child in enumerate(onecell):
	# 			coordinates_keep2[idx1-1,idx2] = int(child.text)
	# 	print(fpath1)
	# 	print(coordinates_keep1.shape)
	# 	print(fpath2)
	# 	print(coordinates_keep2.shape)
	# 	print('concatenated')
	# 	coordinates_keep = np.array([*coordinates_keep1,*coordinates_keep2])
	# 	print(coordinates_keep.shape)
	# 	print()

## ORIGINAL CODE __________________________________________
		
## NEW CODE ROB 01/31/24 ------------------------------------------
	fpath = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name,"points/cell_classification.xml")
	if os.path.exists(fpath):
		tree = ET.parse(fpath)
		root = tree.getroot()
		n_cells = len(root[1][2])
		coordinates_keep = np.empty([n_cells-1, 3], dtype=int)
		for idx1 in range(1,n_cells):
			onecell = root[1][2][idx1]
			for idx2,child in enumerate(onecell):
				coordinates_keep[idx1-1,idx2] = int(child.text)
	else:
		fpath1 = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name+'-pt1',"points/cell_classification.xml")
		tree = ET.parse(fpath1)
		root = tree.getroot()
		n_cells = len(root[1][2])
		coordinates_keep1 = np.empty([n_cells-1, 3], dtype=int)
		for idx1 in range(1,n_cells):
			onecell = root[1][2][idx1]
			for idx2,child in enumerate(onecell):
				coordinates_keep1[idx1-1,idx2] = int(child.text)


		fpath2 = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name+'-pt2',"points/cell_classification.xml")
		tree = ET.parse(fpath2)
		root = tree.getroot()
		n_cells = len(root[1][2])
		coordinates_keep2 = np.empty([n_cells-1, 3], dtype=int)
		for idx1 in range(1,n_cells):
			onecell = root[1][2][idx1]
			for idx2,child in enumerate(onecell):
				coordinates_keep2[idx1-1,idx2] = int(child.text)

		fpath3 = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name+'-pt3',"points/cell_classification.xml")
		tree = ET.parse(fpath3)
		root = tree.getroot()
		n_cells = len(root[1][2])
		coordinates_keep3 = np.empty([n_cells-1, 3], dtype=int)
		for idx1 in range(1,n_cells):
			onecell = root[1][2][idx1]
			for idx2,child in enumerate(onecell):
				coordinates_keep3[idx1-1,idx2] = int(child.text)

		fpath4 = os.path.join("/jukebox/YOUR_DIR/clearmap-output/cellfinder",request_name,sample_name+'-pt4',"points/cell_classification.xml")
		tree = ET.parse(fpath4)
		root = tree.getroot()
		n_cells = len(root[1][2])
		coordinates_keep4 = np.empty([n_cells-1, 3], dtype=int)
		for idx1 in range(1,n_cells):
			onecell = root[1][2][idx1]
			for idx2,child in enumerate(onecell):
				coordinates_keep4[idx1-1,idx2] = int(child.text)


		print(fpath1)
		print(coordinates_keep1.shape)
		print(fpath2)
		print(coordinates_keep2.shape)
		print(fpath3)
		print(coordinates_keep3.shape)
		print(fpath4)
		print(coordinates_keep4.shape)

		print('concatenated')
		coordinates_keep = np.array([*coordinates_keep1,*coordinates_keep2,*coordinates_keep3,*coordinates_keep4])
		print(coordinates_keep.shape)
		print()



## NEW CODE ROB 01/31/24 ------------------------------------------


	raw_filtered_source = ws.source('cells', postfix='raw_filtered')
	coordinates_raw = np.hstack([raw_filtered_source[c][:,None] for c in 'xyz'])
	size = np.hstack([raw_filtered_source[c][:,None] for c in ['size']])
	region = np.hstack([raw_filtered_source[c][:,None] for c in ['region']])
	intensity_raw = np.hstack([raw_filtered_source[c][:,None] for c in ['intensity']])
	transformed_filtered_source = ws.source('cells', postfix='transformed_to_atlas')
	coordinates_transformed = np.hstack([transformed_filtered_source[c][:,None] for c in 'xyz'])

	xyz = np.asarray([(int(X[0]), int(X[1]), int(X[2])) for X in coordinates_transformed])
	for idx, val in enumerate(xyz):
		ID = eroded_atlas_vol[val[2],val[1],val[0]]
		region[idx] = ID

	X = coordinates_raw
	searched_values = coordinates_keep
	dims = X.max(0)+1
	X1D = np.ravel_multi_index(X.T,dims)
	searched_valuesID = np.ravel_multi_index(searched_values.T,dims)
	sidx = X1D.argsort()
	idx_keep = sidx[np.searchsorted(X1D,searched_valuesID,sorter=sidx)]

	coordinates_transformed_keep = coordinates_transformed[idx_keep,:]
	size_keep =  size[idx_keep]
	region_keep =  region[idx_keep]
	intensity_raw_keep =  intensity_raw[idx_keep]
	thresholds_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_filter.p'
	with open(thresholds_file,'rb') as f:
		thresholds_dict = pickle.load(f)

	sfx = str(thresholds_dict['size'][0]) + 'size_' + str(thresholds_dict['intensity'][0]) + "intensity"
	print(coordinates_raw.shape)
	print(coordinates_transformed.shape)
	print(coordinates_keep.shape)
	print(coordinates_transformed_keep.shape)
	print(thresholds_dict)
	print()

	# cells_filtered1 = cells.filter_cells(
	# 				 source = ws.filename('cells', postfix='transformed_to_atlas'),
	# 				 sink = ws.filename('cells', postfix='transformed_to_atlas_TEMP'),
	# 				 thresholds=thresholds_dict)
	# cells_filtered2 = cells.filter_cells(
	# 				 source = ws.filename('cells', postfix='transformed_to_atlas_classified'),
	# 				 sink = ws.filename('cells', postfix='transformed_to_atlas_classified_TEMP'),
	# 				 thresholds=thresholds_dict)
	# df = pd.DataFrame([np.size(ws.source('cells', postfix='transformed_to_atlas_TEMP')),np.size(ws.source('cells', postfix='transformed_to_atlas_classified_TEMP'))])
	# basename_csv = f"cellfinder_stats.csv"
	# savename_csv = os.path.join(workspace_dir,basename_csv)
	# df.to_csv(savename_csv,index=False)

	if 0 < thresholds_dict['intensity'][0] < 500:
		clearmap_params_file = os.path.join(output_rootpath,request_name,
			sample_name,imaging_request,'rawdata/resolution_3.6x')+'/cell_detection_parameter.p'
		with open(clearmap_params_file,'rb') as f:
			cell_detection_parameter = pickle.load(f)
		f = thresholds_dict['intensity'][0]/10
		thresholds_dict['intensity'] = (int(f*cell_detection_parameter['shape_detection']['threshold']),None)
		print('Resetting intensity threshold to '+str(f)+'x background: '+str(thresholds_dict['intensity'][0]))
		print(thresholds_dict)
		print()

	cells_to_save = np.hstack((coordinates_keep,size_keep,intensity_raw_keep,region_keep))
	header = ['x','y','z','size','intensity','region']
	dtypes = [int, int, int, int, float, int]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	savename = ws.filename('cells',postfix='raw_classified')
	io.write(savename,output_array)
	print(f'Saving raw cell detection results to: {savename}')
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print()

	cells_to_save = np.hstack((coordinates_transformed_keep,size_keep,intensity_raw_keep,region_keep))
	header = ['x','y','z','size','intensity','region']
	dtypes = [int, int, int, int, float, int]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	savename = ws.filename('cells',postfix='transformed_to_atlas_classified')
	io.write(savename,output_array)
	cells_filtered = cells.filter_cells(
					 source = ws.filename('cells', postfix='transformed_to_atlas_classified'),
					 sink = ws.filename('cells', postfix='transformed_to_atlas_classified_'+sfx),
					 thresholds=thresholds_dict)
	savename = ws.filename('cells',postfix='transformed_to_atlas_classified_'+sfx)
	io.write(savename,cells_filtered)
	print(f'Saving registered cell detection results to: {savename}')
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print()

	transformed_to_atlas_source = ws.source('cells', postfix='transformed_to_atlas_classified_'+sfx)
	coordinates = np.hstack([transformed_to_atlas_source[c][:,None] for c in 'xyz'])
	size = np.hstack([transformed_to_atlas_source[c][:,None] for c in ['size']])
	region = np.hstack([transformed_to_atlas_source[c][:,None] for c in ['region']])
	intensity_raw = np.hstack([transformed_to_atlas_source[c][:,None] for c in ['intensity']])

	coordinates = np.delete(coordinates,np.where(np.isnan(intensity_raw)),0)
	size = np.delete(size,np.where(np.isnan(intensity_raw)))
	region = np.delete(region,np.where(np.isnan(intensity_raw)))
	intensity_raw = np.delete(intensity_raw,np.where(np.isnan(intensity_raw)))

	print('Generating intensity and count weights...')

	######
	# hemisphere = coordinates[:,2]>=264 # PMA
	hemisphere = coordinates[:,2]>=228 # Allen CCFv3
	######

	intensity_sum_A = np.nansum(intensity_raw[np.where(hemisphere)])
	counts_sum_A = np.count_nonzero(hemisphere)
	intensity_sum_B = np.nansum(intensity_raw[np.where(~hemisphere)])
	counts_sum_B = np.count_nonzero(~hemisphere)
	intensity_weights = np.empty([len(intensity_raw),1], dtype=float)
	count_weights = np.empty([len(intensity_raw),1], dtype=float)
	for idx, val in enumerate(intensity_raw):
		if hemisphere[idx]:
			intensity_weights[idx] = (val/intensity_sum_A)*0.5
			count_weights[idx] = (1/counts_sum_A)*0.5
		elif ~hemisphere[idx]:
			intensity_weights[idx] = (val/intensity_sum_B)*0.5
			count_weights[idx] = (1/counts_sum_B)*0.5
		else:
			print('WARNING')
	print(np.sum(intensity_weights))
	print()
	print(np.sum(count_weights))
	print(np.sum(~hemisphere))
	print(np.sum(hemisphere))
	print('Done!!!')
	print()

	cells_to_save = np.hstack((coordinates,size.reshape((len(size),1)),intensity_raw.reshape((len(intensity_raw),1)),region.reshape((len(region),1)),hemisphere.reshape((len(hemisphere),1)),count_weights,intensity_weights))
	header = ['x','y','z','size','intensity','region','hemisphere','count-weight','intensity-weight']
	dtypes = [int, int, int, int, float, int, int, float, float]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	savename = ws.filename('cells',postfix='transformed_to_atlas_classified_'+sfx)
	io.write(savename,output_array)
	print(f'Saving final transformed/filtered cell detection results to: {savename}')



	transformed_to_atlas_classified = ws.source('cells', postfix='transformed_to_atlas_classified_'+sfx)
	var1 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['x']])
	var1 = [*var1,*var1]
	var2 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['y']])
	var2 = [*var2,*var2]
	var3 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['z']])
	var3b = 455-var3
	var3 = [*var3,*var3b]
	var4 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['size']])
	var4 = [*var4,*var4]
	var5 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['intensity']])
	var5 = [*var5,*var5]
	var6 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['region']])
	var6 = [*var6,*var6]
	var7 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['hemisphere']])
	var7 = [*var7,*var7]
	var8 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['count-weight']])
	var8 = [*var8,*var8]
	var9 = np.hstack([transformed_to_atlas_classified[c][:,None] for c in ['intensity-weight']])
	var9 = [*var9,*var9]
	cells_to_save = np.hstack((var1,var2,var3,var4,var5,var6,var7,var8,var9))
	header = ['x','y','z','size','intensity','region','hemisphere','count-weight','intensity-weight']
	dtypes = [int, int, int, int, float, int, int, float, float]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	savename = ws.filename('cells',postfix='transformed_to_atlas_classified_kde_'+sfx)
	io.write(savename,output_array)
	print(f'Saving final transformed/filtered cell detection results for KDE to: {savename}')
	print()
	print(ws.source('cells', postfix='transformed_to_atlas_classified'))
	print(ws.source('cells', postfix='transformed_to_atlas_classified_'+sfx))
	print(ws.source('cells', postfix='transformed_to_atlas_classified_kde_'+sfx))

	# remove unused voxels from eroded atlas volume
	# if thresholds_dict['x'][0] > 0:
	# 	eroded_atlas_vol[0:thresholds_dict['x'][0]-1,:,:] = 0
	# if thresholds_dict['x'][1] <= 351:
	# 	eroded_atlas_vol[thresholds_dict['x'][1]:-1,:,:] = 0
	# if thresholds_dict['y'][0] > 0:
	# 	eroded_atlas_vol[:,0:thresholds_dict['y'][0]-1,:] = 0
	# if thresholds_dict['y'][1] <= 639:
	# 	eroded_atlas_vol[:,thresholds_dict['y'][1]:-1,:] = 0
	# if thresholds_dict['z'][0] > 0:
	# 	eroded_atlas_vol[:,:,0:thresholds_dict['z'][0]-1] = 0
	# if thresholds_dict['z'][1] <= 539:
	# 	eroded_atlas_vol[:,:,thresholds_dict['z'][1]:-1] = 0
	# if thresholds_dict['region'][0] > 1:
	# 	eroded_atlas_vol[eroded_atlas_vol<thresholds_dict['region'][0]] = 0
	# if thresholds_dict['region'][1] <= 1329:
	# 	eroded_atlas_vol[eroded_atlas_vol>=thresholds_dict['region'][1]] = 0
	# rmv = np.argwhere(np.logical_and(eroded_atlas_vol[:,slice(0,109),:]>=914,eroded_atlas_vol[:,slice(0,109),:]<989))
	# for idx,val in enumerate(rmv):
	# 	eroded_atlas_vol[val[0],val[1],val[2]] = 0

	# if thresholds_dict['x'][0] > 0:
	# 	eroded_atlas_vol[0:thresholds_dict['x'][0]-1,:,:] = 0
	# if thresholds_dict['x'][1] < 320:
	# 	eroded_atlas_vol[thresholds_dict['x'][1]:-1,:,:] = 0
	# if thresholds_dict['y'][0] > 0:
	# 	eroded_atlas_vol[:,0:thresholds_dict['y'][0]-1,:] = 0
	# if thresholds_dict['y'][1] < 528:
	# 	eroded_atlas_vol[:,thresholds_dict['y'][1]:-1,:] = 0
	# if thresholds_dict['z'][0] > 0:
	# 	eroded_atlas_vol[:,:,0:thresholds_dict['z'][0]-1] = 0
	# if thresholds_dict['z'][1] < 456:
	# 	eroded_atlas_vol[:,:,thresholds_dict['z'][1]:-1] = 0
	# if thresholds_dict['region'][0] > 0:
	# 	eroded_atlas_vol[eroded_atlas_vol<thresholds_dict['region'][0]] = 0
	# if thresholds_dict['region'][1] < 1347:
	# 	eroded_atlas_vol[eroded_atlas_vol>=thresholds_dict['region'][1]] = 0
	# rmv = np.argwhere(np.logical_and(eroded_atlas_vol[:,slice(0,111),:]>=388,eroded_atlas_vol[:,slice(0,111),:]<=462))
	# for idx,val in enumerate(rmv):
	# 	eroded_atlas_vol[val[0],val[1],val[2]] = 0

	# Get cell counts for each brain region
	ids = segment_props_dict['inline']['ids']
	segment_names = segment_props_dict['inline']['properties'][0]['values']
	segment_name_dict = {int(ids[ii]):segment_names[ii].split(':')[1].strip() for ii in range(len(ids))}
	print()
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print()

	# Loop over all ids in the atlas JSON, not just the ones in the volume
	print("Determining counts and volumes for each region [whole brain]...")
	print()
	count_dict_names = {}
	intensity_dict_names = {}
	volume_dict_names = {}
	chunk_size = 25 # Each core runs get_count_and_volume() on this many different regions
	src = ws.source('cells', postfix='transformed_to_atlas_classified_'+sfx)
	region = np.hstack([src[c][:,None] for c in ['region']])
	intensity_raw = np.hstack([src[c][:,None] for c in ['intensity']])
	region_idx = region
	intensity_idx = intensity_raw
	chunked_segment_lists = [ids[i:i+chunk_size] for i in range(0,len(ids),chunk_size)]
	get_count_and_volume_par = partial(get_count_and_volume,
		region_idx,intensity_idx,segment_name_dict,eroded_atlas_vol)
	with ProcessPoolExecutor(max_workers=n_cores) as executor:
		for count_dict_i,intensity_dict_i,volume_dict_i in executor.map(get_count_and_volume_par,chunked_segment_lists):
			try:
				for key in count_dict_i:
					count_dict_names[key] = count_dict_i[key]
				for key in intensity_dict_i:
					intensity_dict_names[key] = intensity_dict_i[key]
				for key in volume_dict_i:
					volume_dict_names[key] = volume_dict_i[key]
			except Exception as exc:
				print(f'generated an exception: {exc}')
	sys.stdout.flush()
	print()
	print("Correcting counts [whole brain]...")

	# Correct counts by adding sum of counts in progeny regions
	ontology_graph = graph_tools.Graph(ontology_dict)
	corrected_count_dict = {}
	corrected_intensity_dict = {}
	corrected_volume_dict = {}
	for region in count_dict_names.keys():
		counts_region = count_dict_names[region]
		intensity_region = intensity_dict_names[region]
		volume_region = volume_dict_names[region]
		progeny = ontology_graph.get_progeny(region)
		if progeny != []:
			for prog in progeny:
				try:
					counts_region += count_dict_names[prog]
					intensity_region += intensity_dict_names[prog]
					volume_region += volume_dict_names[prog]
				except KeyError:
					print('error')
					continue
		corrected_count_dict[region] = counts_region
		corrected_intensity_dict[region] = intensity_region
		corrected_volume_dict[region] = volume_region
	print()
	print(corrected_count_dict)
	print()
	sys.stdout.flush()

	# Save region cell counts to region_cell_counts.csv
	df = pd.DataFrame([corrected_count_dict,corrected_intensity_dict,corrected_volume_dict])
	basename_csv = f"region_summary_statistics_classified_{sfx}.csv"
	savename_csv = os.path.join(workspace_dir,basename_csv)
	df.to_csv(savename_csv,index=False)
	print(f'Saving cell detection results broken down by brain region to: {savename_csv}')

	print()
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print()

	# Loop over all ids in the atlas JSON, not just the ones in the volume
	print("Determining counts and volumes for each region [hemisphere 0]...")
	print()
	count_dict_names = {}
	intensity_dict_names = {}
	volume_dict_names = {}
	chunk_size = 25 # Each core runs get_count_and_volume() on this many different regions
	src = ws.source('cells', postfix='transformed_to_atlas_classified_'+sfx)
	region = np.hstack([src[c][:,None] for c in ['region']])
	intensity_raw = np.hstack([src[c][:,None] for c in ['intensity']])
	region_idx_left = np.array(region)[np.where(~hemisphere)]
	intensity_idx_left = np.array(intensity_raw)[np.where(~hemisphere)]
	chunked_segment_lists = [ids[i:i+chunk_size] for i in range(0,len(ids),chunk_size)]
	get_count_and_volume_par = partial(get_count_and_volume,
		region_idx_left,intensity_idx_left,segment_name_dict,eroded_atlas_vol)
	with ProcessPoolExecutor(max_workers=n_cores) as executor:
		for count_dict_i,intensity_dict_i,volume_dict_i in executor.map(get_count_and_volume_par,chunked_segment_lists):
			try:
				for key in count_dict_i:
					count_dict_names[key] = count_dict_i[key]
				for key in intensity_dict_i:
					intensity_dict_names[key] = intensity_dict_i[key]
				for key in volume_dict_i:
					volume_dict_names[key] = volume_dict_i[key]
			except Exception as exc:
				print(f'generated an exception: {exc}')
	sys.stdout.flush()
	print()
	print("Correcting counts [hemisphere 0]...")

	# Correct counts by adding sum of counts in progeny regions
	ontology_graph = graph_tools.Graph(ontology_dict)
	corrected_count_dict = {}
	corrected_intensity_dict = {}
	corrected_volume_dict = {}
	for region in count_dict_names.keys():
		counts_region = count_dict_names[region]
		intensity_region = intensity_dict_names[region]
		volume_region = volume_dict_names[region]
		progeny = ontology_graph.get_progeny(region)
		if progeny != []:
			for prog in progeny:
				try:
					counts_region += count_dict_names[prog]
					intensity_region += intensity_dict_names[prog]
					volume_region += volume_dict_names[prog]
				except KeyError:
					continue
		corrected_count_dict[region] = counts_region
		corrected_intensity_dict[region] = intensity_region
		corrected_volume_dict[region] = volume_region
	print()
	print(corrected_count_dict)
	print()
	sys.stdout.flush()


	# Save region cell counts to region_cell_counts.csv
	df = pd.DataFrame([corrected_count_dict,corrected_intensity_dict,corrected_volume_dict])
	basename_csv = f"region_summary_statistics_classified_hemisphere0_{sfx}.csv"
	savename_csv = os.path.join(workspace_dir,basename_csv)
	df.to_csv(savename_csv,index=False)
	print(f'Saving cell detection results broken down by brain region to: {savename_csv}')

	print()
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print()

	# Loop over all ids in the atlas JSON, not just the ones in the volume
	print("Determining counts and volumes for each region [hemisphere 1]...")
	print()
	count_dict_names = {}
	intensity_dict_names = {}
	volume_dict_names = {}
	chunk_size = 25 # Each core runs get_count_and_volume() on this many different regions
	src = ws.source('cells', postfix='transformed_to_atlas_classified_'+sfx)
	region = np.hstack([src[c][:,None] for c in ['region']])
	intensity_raw = np.hstack([src[c][:,None] for c in ['intensity']])
	region_idx_right = np.array(region)[np.where(hemisphere)]
	intensity_idx_right = np.array(intensity_raw)[np.where(hemisphere)]
	chunked_segment_lists = [ids[i:i+chunk_size] for i in range(0,len(ids),chunk_size)]
	get_count_and_volume_par = partial(get_count_and_volume,
		region_idx_right,intensity_idx_right,segment_name_dict,eroded_atlas_vol)
	with ProcessPoolExecutor(max_workers=n_cores) as executor:
		for count_dict_i,intensity_dict_i,volume_dict_i in executor.map(get_count_and_volume_par,chunked_segment_lists):
			try:
				for key in count_dict_i:
					count_dict_names[key] = count_dict_i[key]
				for key in intensity_dict_i:
					intensity_dict_names[key] = intensity_dict_i[key]
				for key in volume_dict_i:
					volume_dict_names[key] = volume_dict_i[key]
			except Exception as exc:
				print(f'generated an exception: {exc}')
	sys.stdout.flush()
	print()
	print("Correcting counts [hemisphere 1]...")

	# Correct counts by adding sum of counts in progeny regions
	ontology_graph = graph_tools.Graph(ontology_dict)
	corrected_count_dict = {}
	corrected_intensity_dict = {}
	corrected_volume_dict = {}
	for region in count_dict_names.keys():
		counts_region = count_dict_names[region]
		intensity_region = intensity_dict_names[region]
		volume_region = volume_dict_names[region]
		progeny = ontology_graph.get_progeny(region)
		if progeny != []:
			for prog in progeny:
				try:
					counts_region += count_dict_names[prog]
					intensity_region += intensity_dict_names[prog]
					volume_region += volume_dict_names[prog]
				except KeyError:
					continue
		corrected_count_dict[region] = counts_region
		corrected_intensity_dict[region] = intensity_region
		corrected_volume_dict[region] = volume_region
	print()
	print(corrected_count_dict)
	print()
	sys.stdout.flush()


	# Save region cell counts to region_cell_counts.csv
	df = pd.DataFrame([corrected_count_dict,corrected_intensity_dict,corrected_volume_dict])
	basename_csv = f"region_summary_statistics_classified_hemisphere1_{sfx}.csv"
	savename_csv = os.path.join(workspace_dir,basename_csv)
	df.to_csv(savename_csv,index=False)
	print(f'Saving cell detection results broken down by brain region to: {savename_csv}')
	print()
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
