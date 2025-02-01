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
from taskqueue import LocalTaskQueue
import igneous.task_creation as tc
from cloudvolume.lib import touch

from conditions import conditions_dict

sys.path.append('/jukebox/witten/Chris/python/lightsheet-pipeline/kde-visualization/utils/') # to get precomputed_utils in path
from precomputed_utils import (make_info_file,
 calculate_chunks)


layer_parent_dir = '/jukebox/LightSheetData/neuroglancer/public/lightserv/cz15/cfos_heatmaps/cta-final/' # a location on bucket that Neuroglancer can see
pma_annotation_file = '/jukebox/witten/Chris/data/clearmap2/utilities/allen-atlas-cz/annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2_kde.tif'

# path to base folder where cell coordinates can be found
input_dir = "/jukebox/witten/Chris/data/clearmap2"

# where you want to save the results
save_dir = "/jukebox/witten/Chris/data/clearmap2/kde-visualization/cta-final/cgrp_pct_of_pbn/"

# atlas_dir contains the files:
# pma_positions_inbrain_20um.npy,
# pma_positions_20um.npy,
# pma_indices_inbrain_20um.npy
atlas_dir = "/jukebox/witten/Chris/python/lightsheet-pipeline/kde-visualization/atlas_resources/"


def upload_volume_slice_coronal(y,cloud_vol,progress_dir):
    """
    ---PURPOSE---
    Upload a single y plane from a numpy volume to a cloud volume object.
    Made to be run in parallel. Uses a global variable for the image_vol,
    the 3d volume from which we're writing. This is to enable parallel
    processing. If you pass this image_vol in as a parameter it can't
    process multiple planes in parallel.
    ---INPUT---
    z            0-indexed integer corresponding to the z plane to upload
    cloud_vol    CloudVolume object to which you want to upload the slice
    progress_dir Directory with list of files that get touched as slices are uploaded
    """

    print('Processing slice y=',y)
    z_dim,x_dim = image_vol[:,y,:].shape
    array = image_vol[:,y,].reshape((z_dim,1,x_dim)).T

    cloud_vol[:,y,:] = array
    touch(os.path.join(progress_dir, str(y)))
    return "Success"

if __name__ == '__main__':
	start_time = time.time()
	print(sys.argv)
	step = sys.argv[1]
	condition_name = sys.argv[2]
	chunk_size = int(sys.argv[3]) # how many kernel evaluations will take place in this array job, only used for steps
	try:
		chunk_index = int(os.environ.get('SLURM_ARRAY_TASK_ID'))
	except TypeError:
		chunk_index = None

	print("#step,condition_name,chunk_size,chunk_index")
	print(step,condition_name,chunk_size,chunk_index)

	## Step 0: Save all points in a single dataframe
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
				'cells_transformed_to_atlas_classified_kde_350size_0intensity.npy')

			df_pts = pd.DataFrame(np.load(points_pth))
			try:
				df_pts['weight'] = 1/pth_dict['pbn_fos']
				print('Weighting by PB Fos counts')
			except:
				df_pts['weight'] = 1/len(df_pts)
				print('Weighting by total Fos counts')
			pts_condition_container.append(df_pts)

		df_all = pd.concat(pts_condition_container)[['x','y','z','weight']]
		print(df_all)
		print(np.sum(df_all['weight']))
		savename_pts_condition = os.path.join(save_dir,
			f"{condition_name}_cells_transformed_to_atlas_classified.csv")
		df_all.to_csv(savename_pts_condition)
		print(f"Saved {savename_pts_condition}")

	## Step 1: Make the eroded array of positions we will pass to the kernel
	## Also, make the kernel
	if step == 'step1':


		# # Load PMA
		# print("Loading PMA")
		# pma_vol = tifffile.imread(pma_annotation_file)
		# pma_vol_bool = pma_vol.astype('bool')
        #
		# print("making EDT full res")
		# # make edt
		# anisotropy = (25,25,25) # resolution in microns in x,y,z according to edt docs
		# edt_pma_20um = edt.edt(
		#   pma_vol_bool,anisotropy=anisotropy,
		#   black_border=True, order='C',
		#   parallel=6 # number of threads, <= 0 sets to num cpu
		# )
        #
		# #zero out near edges in lower res edt
		# edge_mask_20um = edt_pma_20um < 1
        #
		# # swap axes of x and z in edge mask since meshgrid uses x,y,z order
		# edge_mask_20um_fixedorientation = np.swapaxes(edge_mask_20um,0,2)
        #
		# print("making meshgrid")
		# # Make meshgrid
		# xdim = pma_vol.shape[2] # horizontal
		# ydim = pma_vol.shape[1] # coronal
		# zdim = pma_vol.shape[0] # sagittal
        #
		# nx = xdim*1j
		# ny = ydim*1j
		# nz = zdim*1j
		# X, Y, Z = np.mgrid[0:xdim-1:nx, 0:ydim-1:ny, 0:zdim-1:nz]
        #
		# # Mask out X, Y and Z grids beyond the edges of the volume
		# X[edge_mask_20um_fixedorientation] = 0
		# Y[edge_mask_20um_fixedorientation] = 0
		# Z[edge_mask_20um_fixedorientation] = 0
        #
		# print("making positions array")
		# # Make the positions array that has the right format to pass to the kde kernel
		# positions = np.vstack([X.ravel(), Y.ravel(),Z.ravel()])
        #
		# # Save the positions array and good_indices
		# savename_positions = os.path.join(atlas_dir,'pma_positions_20um.npy')
		# np.save(savename_positions,positions)
		# print(f"Saved {savename_positions}")
        #
		# # Make a list of indices for which position is not [0,0,0]. These are
		# # the only indices of positions at which we need to compute the kernel
		# print("finding good indices inside brain")
		# savename_indices = os.path.join(atlas_dir,'pma_indices_inbrain_20um.npy')
		# if os.path.exists(savename_indices):
		# 	good_indices = np.load(savename_indices)
		# 	print(good_indices.shape)
		# else:
		# 	good_indices = np.array(
		# 		[ii for ii in range(positions.shape[1]) if any(positions[:,ii])])
        #
		# 	# Save the good indices
		# 	print(good_indices.shape)
		# 	np.save(savename_indices,np.array(good_indices))
		# 	print(f"Saved {savename_indices}")
        #
		# # Save the good positions array
		# print("Making good_positions array")
		# good_positions = np.take(positions,good_indices,axis=1)
        #
		# savename_good_positions = os.path.join(atlas_dir,'pma_positions_inbrain_20um.npy')
		# np.save(savename_good_positions,good_positions)
		# print(f"Saved {savename_good_positions}")



		# Make the kernel if it doesn't already exist
		savename_kernel = os.path.join(save_dir,f'{condition_name}_kernel_weighted.p')

		if not os.path.exists(savename_kernel):
			savename_pts_condition = os.path.join(save_dir,
				f"{condition_name}_cells_transformed_to_atlas_classified.csv")
			df_all = pd.read_csv(savename_pts_condition)

			values = np.vstack([df_all['x'],df_all['y'],df_all['z']])
			weights = df_all['weight']
			kernel = stats.gaussian_kde(values,weights=weights)

			# Save the kernel
			with open(savename_kernel,'wb') as outfile:
				pickle.dump(kernel,outfile)
			print(f"Saved {savename_kernel}")

	## Step 2: Make and save the individual array job kernel evaluations on the PMA grid
	if step == 'step2':
		# Load in saved cells
		savename_pts_condition = os.path.join(save_dir,
			f"{condition_name}_cells_transformed_to_atlas_classified.csv")
		df_all = pd.read_csv(savename_pts_condition)

		# load the kernel
		print("loading kernel")
		savename_kernel = os.path.join(save_dir,f'{condition_name}_kernel_weighted.p')
		with open(savename_kernel,'rb') as infile:
			kernel = pickle.load(infile)
		#kernel.set_bandwidth(kernel.factor*0.25)
		print('Estimated KDE bandwidth:')
		print(kernel.factor)
		#kernel.set_bandwidth(kernel.factor*(1/4))
		print('Resetting KDE bandwidth:')

		kernel.set_bandwidth(0.04) ### [standard] or 0.02 [halfbw]

		print(kernel.factor)
		# Load in good positions
		print("loading good positions")
		savename_good_positions = os.path.join(atlas_dir,'pma_positions_inbrain_20um.npy')
		good_positions = np.load(savename_good_positions)
		max_index = good_positions.shape[1] # the total number of positions over all array jobs
		print(good_positions.shape)
		# evaluate the kernel for this chunk
		print(f"Evaluating kde for chunk index: {chunk_index}")
		chunk_min = chunk_index*chunk_size
		chunk_max = min(chunk_min + chunk_size,max_index)
		print(f"Chunk: {chunk_min}-{chunk_max}")
		kde_estimate_array = kernel(good_positions[:,chunk_min:chunk_max])

		# Make chunk save dir if it does not exist
		#chunk_save_dir = os.path.join(save_dir,f'fullres_kde_chunks_{condition_name}_weighted_decr75pc_halfbw')
		chunk_save_dir = os.path.join(save_dir,f'fullres_kde_chunks_{condition_name}_weighted_decr75pc')
		os.makedirs(chunk_save_dir,exist_ok=True)

		# Save chunk file
		print("Saving kde chunk array")
		savename_kde_array = os.path.join(chunk_save_dir,f'chunk{chunk_index}_{chunk_min}-{chunk_max}.npy')
		np.save(savename_kde_array,kde_estimate_array)
		print(f"Saved {savename_kde_array}")

	## Step 3: Put together the individual kernel evals into a single final kde array and save it
	if step == "step3":
		# Load pma volume
		pma_vol = tifffile.imread(pma_annotation_file)

		# Load the positions array
		print("Loading positions")
		savename_positions = os.path.join(atlas_dir,'pma_positions_20um.npy')
		positions = np.load(savename_positions)

		# Load the good positions array (inside brain)
		print("Loading good_indices")
		savename_indices = os.path.join(atlas_dir,'pma_indices_inbrain_20um.npy')
		good_indices = np.load(savename_indices)

		# initialize array to store the full kde result, including zeros outside brain
		full_result = np.zeros(positions.shape[1])

		# initialize array to store kde result that we calculated inside brain
		good_result = np.zeros(len(good_indices))

		# Loop through our chunk files and create good_result array
		print("Combining chunk files into single large array")
		#chunk_dir = os.path.join(save_dir,f'fullres_kde_chunks_{condition_name}_weighted_decr75pc_halfbw')
		chunk_dir = os.path.join(save_dir,f'fullres_kde_chunks_{condition_name}_weighted_decr75pc')
		chunkfiles = sorted(glob.glob(chunk_dir + '/chunk*npy'),
			key=lambda x: int(x.split('_')[-2].split('chunk')[-1]))

		for ii,chunkfile in enumerate(chunkfiles):
			if ii % 10 == 0:
				print(f"{ii+1}/{len(chunkfiles)}")
			chunk_min = ii*chunk_size
			chunk_max = min(chunk_min + chunk_size,len(good_result))
			chunk_result = np.load(chunkfile)
			try:
				good_result[chunk_min:chunk_max] = chunk_result
			except:
				print("Issues with chunk:")
				print(ii,chunk_min,chunk_max)

		# Fill the final kde result array with the subset of kde calculations we just made
		full_result[good_indices]=good_result

		# # reshape to 3D
		kde3d_20um = np.reshape(full_result,pma_vol.shape[::-1]) # have to use reverse PMA since that was the meshgrid order

		print(np.sum(kde3d_20um))
		kde3d_20um = kde3d_20um / np.sum(kde3d_20um) * 100 # normalize to 100%
		print(np.sum(kde3d_20um))

		# condition_dict = conditions_dict.get(condition_name)
		# brains = condition_dict.get('brains')
		# x = 0
		# for pth_dict in brains:
		# 	x = x + pth_dict['pbn_fos']
		# print(x)
		# savename_pts_condition = os.path.join(save_dir,
		# 	f"{condition_name}_cells_transformed_to_atlas_classified.csv")
		# df_all = pd.read_csv(savename_pts_condition)
		# print(len(df_all['weight']))
		# kde3d_20um = kde3d_20um * len(df_all['weight'])/x
		# print(np.sum(kde3d_20um))

		# Swap orientation back to PMA since meshgrid swapped x and z
		kde3d_20um_fixedorientation = np.swapaxes(kde3d_20um,0,2) / np.power(0.025,3) # scale it to be per mm^3
		print(np.sum(kde3d_20um_fixedorientation))

		# # Save resulting kde array to .npy file
		#savename_kde_condition = os.path.join(save_dir,f'kde_20um_{condition_name}_weighted_decr75pc_halfbw.npy')
		savename_kde_condition = os.path.join(save_dir,f'kde_20um_{condition_name}_weighted_decr75pc.npy')
		np.save(savename_kde_condition,kde3d_20um_fixedorientation)
		print(f"Saved {savename_kde_condition}")

	## Step 4: Make precomputed layer of heatmap
	if step == "step4":
		# Load pma volume
		pma_vol = tifffile.imread(pma_annotation_file)

		# Make precomputed info file
		volume_size = pma_vol.shape[::-1] # z,y,x, -> x,y,z
		resolution = (25000,25000,25000) # nm
		layer_name = condition_name + '_celldensity_pma_20um_weighted_decr75pckernel'
		layer_dir = os.path.join(layer_parent_dir,layer_name)
		layer_type = 'image'
		data_type = 'float32'

		# Make coronal chunks
		cloud_vol = make_info_file(
			volume_size,
			resolution,
			layer_dir,
			layer_type,
			data_type,
			voxel_offset=[0,0,0],
			chunk_size=[1024,1,1024],
			description=f'Cfos density heatmap for condition: {condition_name}',
			commit=True)

		# Make precomputed layer
		# First, load completed kde file
		savename_kde_condition = os.path.join(save_dir,f'kde_20um_{condition_name}_weighted_decr75pc.npy')
		image_vol = np.load(savename_kde_condition).astype('float32')

		# Make progress dir
		progress_parent_dir = layer_parent_dir + '/progress_dirs'
		progress_dir = os.path.join(
			progress_parent_dir,f'progress_{layer_name}')
		os.makedirs(progress_dir,exist_ok=True) # unlike os.mkdir doesn't crash on prexisting

		done_files = set([ int(y) for y in os.listdir(progress_dir) ])
		all_files = set(range(cloud_vol.bounds.minpt.y, cloud_vol.bounds.maxpt.y))

		to_upload = [ int(y) for y in list(all_files.difference(done_files)) ]
		to_upload.sort()
		print("Remaining slices to upload are:",to_upload)

		upload_volume_slice_par = partial(
			upload_volume_slice_coronal,
			cloud_vol=cloud_vol,
			progress_dir=progress_dir)

		with ProcessPoolExecutor(max_workers=32) as executor:
			for res in executor.map(upload_volume_slice_par, to_upload):
				try:
					print(res)
				except Exception as exc:
					print(f'generated an exception: {exc}')

	## Step 5: Make log(kde_2/kde_1) array and precomputed layer for it
	if step == "step5":
		conditions = list(conditions_dict.keys())
		print(conditions)
		# Load kde 1
		savename_kde_condition1 = os.path.join(save_dir,f'kde_20um_{conditions[0]}_weighted_decr75pc.npy')
		kde_1 = np.load(savename_kde_condition1).astype('float32')

		# Load kde 2
		savename_kde_condition2 = os.path.join(save_dir,f'kde_20um_{conditions[1]}_weighted_decr75pc.npy')
		kde_2 = np.load(savename_kde_condition2).astype('float32')

		# # Make log(kde_2/kde_1)
		# image_vol = np.log(kde_2/kde_1)
		# np.nan_to_num(image_vol, copy=False, nan=0, posinf=0, neginf=-0)

		# Make log(kde_2/kde_1)
		L = 0.1
		image_vol = np.log((kde_2+L)/(kde_1+L))
		np.nan_to_num(image_vol, copy=False, nan=0, posinf=0, neginf=-0)

		# Save result
		savename_result = os.path.join(save_dir,f'log_kde_20um_{conditions[1]}_over_{conditions[0]}.npy')
		np.save(savename_result,image_vol)
		print(f"Saved {savename_result}")

		# Make info file
		volume_size = image_vol.shape[::-1] # z,y,x, -> x,y,z
		resolution = (25000,25000,25000) # nm

		layer_name = f'log_celldensity_20um_{conditions[1]}_over_{conditions[0]}'
		layer_dir = os.path.join(layer_parent_dir,layer_name)
		layer_type = 'image'
		data_type = 'float32'

		cloud_vol = make_info_file(
			volume_size,
			resolution,
			layer_dir,
			layer_type,
			data_type,
			voxel_offset=[0,0,0],
			chunk_size=[1024,1,1024],
			description=f'log(kde_{conditions[1]}_over_{conditions[0]})',
			commit=True)

		# Make progress dir
		progress_parent_dir = layer_parent_dir + '/progress_dirs'
		progress_dir = os.path.join(
			progress_parent_dir,f'progress_{layer_name}')
		os.makedirs(progress_dir,exist_ok=True) # unlike os.mkdir doesn't crash on prexisting

		done_files = set([ int(y) for y in os.listdir(progress_dir) ])
		all_files = set(range(cloud_vol.bounds.minpt.y, cloud_vol.bounds.maxpt.y))

		to_upload = [ int(y) for y in list(all_files.difference(done_files)) ]
		to_upload.sort()
		print("Remaining slices to upload are:",to_upload)

		upload_volume_slice_par = partial(
			upload_volume_slice_coronal,
			cloud_vol=cloud_vol,
			progress_dir=progress_dir)

		with ProcessPoolExecutor(max_workers=32) as executor:
			for res in executor.map(upload_volume_slice_par, to_upload):
				try:
					print(res)
				except Exception as exc:
					print(f'generated an exception: {exc}')

	## Step 6: Make difference map (kde_2-kde_1)
	## and precomputed layer for it
	if step == "step6":
		conditions = list(conditions_dict.keys())
		print(conditions)
		# Load kde 1
		savename_kde_condition1 = os.path.join(save_dir,f'kde_20um_{conditions[0]}_weighted_decr75pc.npy')
		kde_1 = np.load(savename_kde_condition1).astype('float32')

		# Load kde 2
		savename_kde_condition2 = os.path.join(save_dir,f'kde_20um_{conditions[1]}_weighted_decr75pc.npy')
		kde_2 = np.load(savename_kde_condition2).astype('float32')

		# Make difference map
		image_vol = (kde_2-kde_1)
		print(image_vol.shape)

		# Save result
		savename_result = os.path.join(save_dir,f'diff_kde_20um_{conditions[1]}_{conditions[0]}.npy')
		np.save(savename_result,image_vol)
		print(f"Saved {savename_result}")

		# Make info file
		volume_size = image_vol.shape[::-1] # z,y,x, -> x,y,z
		resolution = (25000,25000,25000) # nm

		layer_name = f'diff_celldensity_20um_{conditions[1]}_{conditions[0]}'
		layer_dir = os.path.join(layer_parent_dir,layer_name)
		layer_type = 'image'
		data_type = 'float32'

		cloud_vol = make_info_file(
			volume_size,
			resolution,
			layer_dir,
			layer_type,
			data_type,
			voxel_offset=[0,0,0],
			chunk_size=[1024,1,1024],
			description=layer_name,
			commit=True)

		# Make progress dir
		progress_parent_dir = layer_parent_dir + '/progress_dirs'
		progress_dir = os.path.join(
			progress_parent_dir,f'progress_{layer_name}')
		os.makedirs(progress_dir,exist_ok=True) # unlike os.mkdir doesn't crash on prexisting

		done_files = set([ int(y) for y in os.listdir(progress_dir) ])
		all_files = set(range(cloud_vol.bounds.minpt.y, cloud_vol.bounds.maxpt.y))

		to_upload = [ int(y) for y in list(all_files.difference(done_files)) ]
		to_upload.sort()
		print("Remaining slices to upload are:",to_upload)

		upload_volume_slice_par = partial(
			upload_volume_slice_coronal,
			cloud_vol=cloud_vol,
			progress_dir=progress_dir)

		with ProcessPoolExecutor(max_workers=32) as executor:
			for res in executor.map(upload_volume_slice_par, to_upload):
				try:
					print(res)
				except Exception as exc:
					print(f'generated an exception: {exc}')

	## Step 6: Mean map
	## and precomputed layer for it
	if step == "step7":
		conditions = list(conditions_dict.keys())
		print(conditions)
		# Load kde 1
		#savename_kde_condition1 = os.path.join(save_dir,f'kde_20um_{conditions[0]}_weighted_decr75pc_halfbw.npy')
		savename_kde_condition1 = os.path.join(save_dir,f'kde_20um_{conditions[0]}_weighted_decr75pc.npy')
		kde_1 = np.load(savename_kde_condition1).astype('float32')

		# Load kde 2
		#savename_kde_condition2 = os.path.join(save_dir,f'kde_20um_{conditions[1]}_weighted_decr75pc_halfbw.npy')
		savename_kde_condition2 = os.path.join(save_dir,f'kde_20um_{conditions[1]}_weighted_decr75pc.npy')
		kde_2 = np.load(savename_kde_condition2).astype('float32')

		# Make mean map
		image_vol = np.mean(np.array([kde_2,kde_1]),axis=0 )
		print(image_vol.shape)

		# Save result
		savename_result = os.path.join(save_dir,f'mean_kde_20um_{conditions[1]}_{conditions[0]}.npy')
		np.save(savename_result,image_vol)
		print(f"Saved {savename_result}")

		# Make info file
		volume_size = image_vol.shape[::-1] # z,y,x, -> x,y,z
		resolution = (25000,25000,25000) # nm

		layer_name = f'mean_celldensity_20um_{conditions[1]}_{conditions[0]}'
		layer_dir = os.path.join(layer_parent_dir,layer_name)
		layer_type = 'image'
		data_type = 'float32'

		cloud_vol = make_info_file(
			volume_size,
			resolution,
			layer_dir,
			layer_type,
			data_type,
			voxel_offset=[0,0,0],
			chunk_size=[1024,1,1024],
			description=layer_name,
			commit=True)

		# Make progress dir
		progress_parent_dir = layer_parent_dir + '/progress_dirs'
		progress_dir = os.path.join(
			progress_parent_dir,f'progress_{layer_name}')
		os.makedirs(progress_dir,exist_ok=True) # unlike os.mkdir doesn't crash on prexisting

		done_files = set([ int(y) for y in os.listdir(progress_dir) ])
		all_files = set(range(cloud_vol.bounds.minpt.y, cloud_vol.bounds.maxpt.y))

		to_upload = [ int(y) for y in list(all_files.difference(done_files)) ]
		to_upload.sort()
		print("Remaining slices to upload are:",to_upload)

		upload_volume_slice_par = partial(
			upload_volume_slice_coronal,
			cloud_vol=cloud_vol,
			progress_dir=progress_dir)

		with ProcessPoolExecutor(max_workers=32) as executor:
			for res in executor.map(upload_volume_slice_par, to_upload):
				try:
					print(res)
				except Exception as exc:
					print(f'generated an exception: {exc}')

	print("--- Took %s seconds ---" % (time.time() - start_time))
