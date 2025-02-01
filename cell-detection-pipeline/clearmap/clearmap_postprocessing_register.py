# General imports
import os,sys,pickle
import numpy as np
import pandas as pd
import json
import tifffile
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor
from brain_atlas_toolkit import graph_tools

#Added by Rob
import itertools

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
		eroded_atlas_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/allen-atlas-cz/annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2.tif'
		segment_props_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/allen-atlas-cz/allenatlas_2017_16bit_hierarch_labels_segment_properties_info_cz'
		ontology_json_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/allen-atlas-cz/allen_hierarch_labels_fillmissing_cz.json'
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

	# Set thresholds for initial cell filtering
	init_thresholds_dict = {}
	init_thresholds_dict['size'] = (200,None) # previously 200
	cells_filtered = cells.filter_cells(
					 source = ws.filename('cells', postfix='raw'),
					 sink = ws.filename('cells', postfix='raw_filtered'),
					 thresholds=init_thresholds_dict)
	savename = ws.filename('cells',postfix='raw_filtered')
	io.write(savename,cells_filtered)

	# Load detected cells from cells_raw.npy
	raw_source = ws.source('cells', postfix='raw_filtered')
	size = np.hstack([raw_source[c][:,None] for c in ['size']])
	coordinates_raw = np.hstack([raw_source[c][:,None] for c in 'xyz'])

	# Swap axes to go from horizontal to sagittal orientation
	coordinates_raw_swapped_axes = np.zeros_like(coordinates_raw)
	coordinates_raw_swapped_axes[:,0] = coordinates_raw[:,2]
	coordinates_raw_swapped_axes[:,1] = coordinates_raw[:,1]
	coordinates_raw_swapped_axes[:,2] = coordinates_raw[:,0]

	# Transform cell coordinates from raw 642-space to downsampled 642-space
	coordinates_resampled = res.resample_points(
					coordinates_raw_swapped_axes, sink=None, orientation=None,
					source_shape=io.shape(ws.filename('stitched'))[::-1],
					sink_shape=io.shape(ch642_downsized_file))

	# Transform cell coordinates from 642-space to 488-space
	coordinates_aligned_to_488 = elx.transform_points(
					coordinates_resampled, sink=None,
					transform_directory=os.path.join(elastix_inverse_dir,
						'488_to_642'),
					temp_file='/tmp/elastix_input_pipeline_'+sample_name+'_'+imaging_request+'_step0.bin',
					result_directory='/tmp/elastix_output_pipeline_'+sample_name+'_'+imaging_request+'_step0')

	# Change permissions back since transformix alters permissions when it is run
	os.chmod(os.path.join(elastix_inverse_dir,
		'488_to_642','TransformParameters.0.txt'),0o777)
	os.chmod(os.path.join(elastix_inverse_dir,
		'488_to_642','TransformParameters.1.txt'),0o777)

	# Transform cell coordinates from 488-space to atlas-space
	coordinates_aligned_to_atlas = elx.transform_points(
					coordinates_aligned_to_488, sink=None,
					transform_directory=elastix_inverse_dir,
					binary=True, indices=False,
					temp_file='/tmp/elastix_input_pipeline_'+sample_name+'_'+imaging_request+'_step1.bin',
					result_directory='/tmp/elastix_output_pipeline_'+sample_name+'_'+imaging_request+'_step1')

	# Change permissions back since transformix alters permissions when it is run
	os.chmod(os.path.join(elastix_inverse_dir,'TransformParameters.0.txt'),0o777)
	os.chmod(os.path.join(elastix_inverse_dir,'TransformParameters.1.txt'),0o777)

	# Load atlas files
	eroded_atlas_vol = np.array(tifffile.imread(eroded_atlas_file)).astype('uint16')
	atlas_segments = np.unique(eroded_atlas_vol)
	atlas_segments = np.array([x for x in atlas_segments if x!=0])

	with open(segment_props_file,'r') as infile:
		segment_props_dict = json.load(infile)

	with open(ontology_json_file,'r') as infile:
		ontology_dict = json.load(infile)

	# remove unused voxels from eroded atlas volume
	# thresholds_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_filter.p'
	# with open(thresholds_file,'rb') as f:
	# 	thresholds_dict = pickle.load(f)

	#####
	# rmv = np.argwhere(np.logical_and(eroded_atlas_vol[:,slice(0,109),:]>=914,eroded_atlas_vol[:,slice(0,109),:]<989))
	# for idx,val in enumerate(rmv):
	# 	eroded_atlas_vol[val[0],val[1],val[2]] = 0
	# init_thresholds_dict['region'] = (2,1105)
	# rmv = np.argwhere(np.logical_and(eroded_atlas_vol[:,slice(0,111),:]>=388,eroded_atlas_vol[:,slice(0,111),:]<=462))
	# for idx,val in enumerate(rmv):
	# 	eroded_atlas_vol[val[0],val[1],val[2]] = 0
	init_thresholds_dict['region'] = (2,1115)
	#####

	# Record the brain region ID where a cell is detected, 0 if not in a region
	cell_regions = np.empty([len(coordinates_aligned_to_atlas), 1], dtype=int)
	xyz = np.asarray([(int(X[0]), int(X[1]), int(X[2])) for X in coordinates_aligned_to_atlas])
	source = ws.source('stitched')
	intensity = np.empty([len(coordinates_raw), 1], dtype=float)
	width = 1
	dist = 3
	clearmap_params_file = '/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_parameter.p'
	with open(clearmap_params_file,'rb') as f:
		cell_detection_parameter = pickle.load(f)
	from scipy import ndimage
	import ClearMap.Utils.HierarchicalDict as hdict
	for idx, val in enumerate(xyz):
		try:
			if np.amin(val)>=0:
				ID = eroded_atlas_vol[val[2],val[1],val[0]]
				cell_regions[idx] = ID
				if (cell_regions[idx]>=init_thresholds_dict['region'][0]) and (cell_regions[idx]<init_thresholds_dict['region'][1]):
					if (size[idx]>=init_thresholds_dict['size'][0]):
						neighbors = np.setdiff1d(np.argwhere(np.amax(np.abs(coordinates_raw[slice(idx-50,idx+51),:]-coordinates_raw[idx,:]),1)<=dist)+idx-50,idx)
						neighbors = neighbors[np.argwhere(size[neighbors]>=init_thresholds_dict['size'][0])]
						if neighbors.size>0:
							if (np.all(size[neighbors]<size[idx])) or ((np.any(size[neighbors]==size[idx])) and (np.all(neighbors>idx))):
								data = source[slice(coordinates_raw[idx,0]-22,coordinates_raw[idx,0]+23),slice(coordinates_raw[idx,1]-22,coordinates_raw[idx,1]+23),slice(coordinates_raw[idx,2]-2,coordinates_raw[idx,2]+3)]
								data = cells.remove_background(ndimage.median_filter(data,3),cell_detection_parameter['background_correction']['shape'],cell_detection_parameter['background_correction']['form'])
								intensity[idx] = np.mean(np.asarray(data[slice(21,24),slice(21,24),slice(1,4)]))
								#intensity[idx] = np.mean(np.asarray(source[slice(coordinates_raw[idx,0]-width,coordinates_raw[idx,0]+width),slice(coordinates_raw[idx,1]-width,coordinates_raw[idx,1]+width),slice(coordinates_raw[idx,2]-width,coordinates_raw[idx,2]+width)]))
							else:
								cell_regions[idx] = 0
								intensity[idx] = np.nan
						else:
							data = source[slice(coordinates_raw[idx,0]-22,coordinates_raw[idx,0]+23),slice(coordinates_raw[idx,1]-22,coordinates_raw[idx,1]+23),slice(coordinates_raw[idx,2]-2,coordinates_raw[idx,2]+3)]
							data = cells.remove_background(ndimage.median_filter(data,3),cell_detection_parameter['background_correction']['shape'],cell_detection_parameter['background_correction']['form'])
							intensity[idx] = np.mean(np.asarray(data[slice(21,24),slice(21,24),slice(1,4)]))
							#intensity[idx] = np.mean(np.asarray(source[slice(coordinates_raw[idx,0]-width,coordinates_raw[idx,0]+width),slice(coordinates_raw[idx,1]-width,coordinates_raw[idx,1]+width),slice(coordinates_raw[idx,2]-width,coordinates_raw[idx,2]+width)]))
					else:
						intensity[idx] = np.nan
				else:
					intensity[idx] = np.nan
			else:
				cell_regions[idx] = 0
				intensity[idx] = np.nan
		except Exception as e:
			cell_regions[idx] = 0
			intensity[idx] = np.nan
			pass

	# Add brain region ID to raw cell array
	cells_to_save = np.hstack((coordinates_raw,size,intensity,cell_regions))
	header = ['x','y','z','size','intensity','region']
	dtypes = [int, int, int, int, float, int]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	savename = ws.filename('cells',postfix='raw_filtered')
	io.write(savename,output_array)
	cells_filtered = cells.filter_cells(
					 source = ws.filename('cells', postfix='raw_filtered'),
					 sink = ws.filename('cells', postfix='raw_filtered'),
					 thresholds=init_thresholds_dict)
	savename = ws.filename('cells',postfix='raw_filtered')
	io.write(savename,cells_filtered)
	print()
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print(f'Saving filtered cell detection results to: {savename}')

## ORIGINAL CODE _________________

	# # save xml file to read into napari
	# data = ws.source('cells', postfix='raw_filtered')
	# cellArray = []
	# for i in data:
	#     cellArray.append(i)
	# root = ET.Element("CellCounter_Marker_File")
	# prop = ET.SubElement(root, "Image_Properties")
	# marker = ET.SubElement(root, "Marker_Data")
	# ET.SubElement(marker, "Current_Type").text = "1"
	# type = ET.SubElement(marker, "Marker_Type")
	# ET.SubElement(prop, "Image_Filename").text = "placeholder.tif"
	# ET.SubElement(type, "Type").text = "2"
	# for i in cellArray:
	#     point = ET.SubElement(type, "Marker")
	#     ET.SubElement(point, "MarkerX").text = str(i[0])
	#     ET.SubElement(point, "MarkerY").text = str(i[1])
	#     ET.SubElement(point, "MarkerZ").text = str(i[2])
	# tree = ET.ElementTree(root)
	# #ET.indent(tree, space="\t", level=0)
	# tree.write(os.path.join(workspace_dir,"cells_raw_filtered.xml"), encoding="utf-8", xml_declaration=True)


	## ORIGINAL CODE _________________


	## NEW CODE BELOW Rob 01/31/24    ---------------------------------------

# split cell list into 4 seperate xml files to read into napari/cellfinder - we split into 4 so that 1) We never timeout on cellfinder which can take days to run with a lot of cells and 2) We can possibly speed up the cellfinder step by parallelizing the classification into 4 batch jobs instead of 1.
	data = ws.source('cells', postfix='raw_filtered')
	cellArray = []
	for i in data:
	    cellArray.append(i)

	splitSize = int(np.ceil(len(cellArray)/4))
	splitCellArray = [list(filter(lambda x: x is not None, subList)) for subList in itertools.zip_longest(*([iter(cellArray)]*splitSize),fillvalue = None)]

# xml file generation
	for j,subList in enumerate(splitCellArray):
		root = ET.Element("CellCounter_Marker_File")
		prop = ET.SubElement(root, "Image_Properties")
		marker = ET.SubElement(root, "Marker_Data")
		ET.SubElement(marker, "Current_Type").text = "1"
		type = ET.SubElement(marker, "Marker_Type")
		ET.SubElement(prop, "Image_Filename").text = "placeholder.tif"
		ET.SubElement(type, "Type").text = "2"
		for i in subList:
			point = ET.SubElement(type, "Marker")
			ET.SubElement(point, "MarkerX").text = str(i[0])
			ET.SubElement(point, "MarkerY").text = str(i[1])
			ET.SubElement(point, "MarkerZ").text = str(i[2])
		tree = ET.ElementTree(root)
		#ET.indent(tree, space="\t", level=0)
		tree.write(os.path.join(workspace_dir,"cells_raw_filtered_pt" + str(j+1) + ".xml"), encoding="utf-8", xml_declaration=True)

	## NEW CODE ABOVE---------------------------------------



	# Add brain region ID to transformed cell array
	cells_to_save = np.hstack((coordinates_aligned_to_atlas,size,intensity,cell_regions))
	header = ['x','y','z','size','intensity','region']
	dtypes = [int, int, int, int, float, int]
	dt = {'names' : header, 'formats' : dtypes}
	output_array = np.zeros(len(cells_to_save), dtype=dt)
	for i,h in enumerate(header):
		output_array[h] = cells_to_save[:,i]
	savename = ws.filename('cells',postfix='transformed_to_atlas')
	io.write(savename,output_array)
	cells_filtered = cells.filter_cells(
					 source = ws.filename('cells', postfix='transformed_to_atlas'),
					 sink = ws.filename('cells', postfix='transformed_to_atlas'),
					 thresholds=init_thresholds_dict)
	savename = ws.filename('cells',postfix='transformed_to_atlas')
	io.write(savename,cells_filtered)
	print(f'Saving registered cell detection results to: {savename}')
	print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
	print()
