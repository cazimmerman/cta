# General imports
import os,sys,pickle
import numpy as np
import pandas as pd
import json
import tifffile
from concurrent.futures import ProcessPoolExecutor
from brain_atlas_toolkit import graph_tools

# ClearMap2 imports
sys.path.append('/jukebox/witten/Chris/python/ClearMap2-master')
os.chmod("/jukebox/witten/Chris/python/ClearMap2-master/ClearMap/External/elastix/build/bin/transformix",0o777)
import ClearMap.IO.Workspace as wsp
import ClearMap.IO.IO as io
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Settings as settings
import ClearMap.Alignment.Resampling as res
import ClearMap.Alignment.Elastix as elx

from functools import partial

def get_count_and_volume(segment_name_dict,eroded_atlas_vol,segment_list,):
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
		if int(atlas_segment_id) not in atlas_segments:
			count_dict[segment_name] = 0
			intensity_dict[segment_name] = 0
			volume_dict[segment_name] = 0
		else:
			count_dict[segment_name] = 1
			intensity_dict[segment_name] = 1
			volume_dict[segment_name] = 1
	return count_dict,intensity_dict,volume_dict

if __name__ == "__main__":
	n_cores = os.cpu_count()

	eroded_atlas_file = '/jukebox/witten/Chris/data/clearmap2/utilities/allen-atlas-cz/annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz.tif'
	segment_props_file = '/jukebox/witten/Chris/data/clearmap2/utilities/allen-atlas-cz/allenatlas_2017_16bit_hierarch_labels_segment_properties_info'
	ontology_json_file = '/jukebox/witten/Chris/data/clearmap2/utilities/allen-atlas-cz/allen_hierarch_labels_fillmissing.json'

	eroded_atlas_vol = np.array(tifffile.imread(eroded_atlas_file)).astype('uint16')
	atlas_segments = np.unique(eroded_atlas_vol)
	atlas_segments = np.array([x for x in atlas_segments if x!=0])

	with open(segment_props_file,'r') as infile:
		segment_props_dict = json.load(infile)

	with open(ontology_json_file,'r') as infile:
		ontology_dict = json.load(infile)

	ontology_graph = graph_tools.Graph(ontology_dict)
	ontology_graph.print_branch('root')
	out = ontology_graph.make_html_ontology()
	Func = open('/jukebox/witten/Chris/clearmap2/utilities/allen-atlas-cz/allen_CCFv3_ontology_cz/allen_CCFv3_annotation_cz.html',"w")
	Func.write(out)
	Func.close()
#	print(out)
	quit()

	ids = segment_props_dict['inline']['ids']
	segment_names = segment_props_dict['inline']['properties'][0]['values']
	segment_name_dict = {int(ids[ii]):segment_names[ii].split(':')[1].strip() for ii in range(len(ids))}
	# Loop over all ids in the atlas JSON, not just the ones in the volume
	count_dict_names = {}
	intensity_dict_names = {}
	volume_dict_names = {}
	chunk_size = 25 # Each core runs get_count_and_volume() on this many different regions


	chunked_segment_lists = [ids[i:i+chunk_size] for i in range(0,len(ids),chunk_size)]
	get_count_and_volume_par = partial(get_count_and_volume,
		segment_name_dict,eroded_atlas_vol)
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
	print("This is the thing...")

	print(type(count_dict_names.keys()))
	# Correct counts by adding sum of counts in progeny regions
	ontology_graph = graph_tools.Graph(ontology_dict)
	regionlist = list(count_dict_names.keys())
	for idxxxx,region in enumerate(regionlist):
		print('   ')
		print('   ')
		print(region)
		print('---------------------------')
		progeny = ontology_graph.get_progeny(region)
		x = np.empty([len(progeny)])
		if progeny != []:
			for idx,prog in enumerate(progeny):
				i = 0
				while i < len(regionlist):
					if prog == regionlist[i]:
						x[idx] = i
					i += 1
				print(prog)
		fname = '/jukebox/witten/Chris/clearmap2/utilities/allen-atlas-cz/allen_CCFv3_ontology_cz/region_' + str(idxxxx) + '.p'
		with open(fname, 'wb') as f:
		    pickle.dump(x,f)
		print(x)
	print("Corrected counts. Now saving out final CSV file")
	sys.stdout.flush()
