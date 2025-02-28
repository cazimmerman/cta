import tifffile, glob, os, sys
import pickle
import numpy as np
from scipy import stats,ndimage
from concurrent.futures import ProcessPoolExecutor
import math
import pandas as pd
from functools import partial, reduce
import edt
import matplotlib.pyplot as plt
import time
from conditions import conditions_dict

sys.path.append('/jukebox/witten/cz15/python/lightsheet-pipeline/kde-visualization/utils/') # to get precomputed_utils in path
from precomputed_annotations import make_precomputed_annotation_layer

layer_parent_dir = '/jukebox/LightSheetData/neuroglancer/public/lightserv/cz15/cfos_heatmaps/filtered_200px'

pma_annotation_file = '/jukebox/LightSheetData/atlas/annotation_sagittal_atlas_20um_16bit_hierarch_labels.tif'

input_dir = "/jukebox/witten/Chris/data/clearmap2"
save_dir = "/jukebox/wang/ahoag/for_cz/cfos_heatmaps/output_filtered_200px"

if __name__ == '__main__':
	start_time = time.time()
	print(sys.argv)
	step = sys.argv[1]
	condition_name = sys.argv[2]

	## Step 0: Save all points from this condition into a single dataframe,
	## include weight on points

	condition_dict = conditions_dict.get(condition_name)
	brains = condition_dict.get('brains')
	if step == 'step0':
		pts_condition_container = []
		for pth_dict in brains:
			points_pth = os.path.join(input_dir,
				pth_dict['request_name'],
				pth_dict['sample_name'],
				pth_dict['imaging_request'],
				'rawdata',
				'resolution_3.6x',
				'cells_transformed_to_atlas_filtered_200px.npy')

			df_pts = pd.DataFrame(np.load(points_pth))
			df_pts['weight'] = 1/len(df_pts)
			pts_condition_container.append(df_pts)

		df_all = pd.concat(pts_condition_container)[['x','y','z','weight']]
		savename_pts_condition = os.path.join(save_dir,
			f"{condition_name}_cells_transformed_to_atlas_filtered_200px.csv")
		df_all.to_csv(savename_pts_condition)
		print(f"Saved {savename_pts_condition}")

	## Step 1: Make precomputed layer of all cells in this condition
	if step == 'step1':

		# Define layer directory
		layer_name = f'{condition_name}_all_points'
		layer_dir = os.path.join(layer_parent_dir,layer_name)

		savename_pts_condition = os.path.join(save_dir,
			f"{condition_name}_cells_transformed_to_atlas_filtered_200px.csv")
		df_all = pd.read_csv(savename_pts_condition)

		# just get points
		coordinates = df_all[['x','y','z']].to_numpy().astype('float32')
		unique_coordinates = np.unique(coordinates,axis=0)
		make_precomputed_annotation_layer(
			unique_coordinates,
			layer_dir,
			grid_shape=[1,1,1],
         	chunk_size=[352,640,540],
         	dimensions_m=[2e-05,2e-05,2e-05],
         	limit=10000,
         	debug=False)

	print("--- Took %s seconds ---" % (time.time() - start_time))
