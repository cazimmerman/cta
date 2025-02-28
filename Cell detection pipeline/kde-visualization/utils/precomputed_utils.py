from collections import defaultdict
from cloudvolume import CloudVolume
from cloudvolume.lib import touch
import os
import numpy as np


def make_info_file(volume_size,resolution,layer_dir,layer_type,data_type,
    voxel_offset=[0,0,0],chunk_size=[1024,1024,1],description='Dummy description',commit=True):
    """ 
    ---PURPOSE---
    Make the cloudvolume info file
    ---INPUT---
    volume_size     [Nx,Ny,Nz] in voxels, e.g. [2160,2560,1271]
    resolution      [size of x pix in nm,size of y pix in nm,size of z pix in nm], e.g. [5000,5000,10000]
    layer_dir       Where to save the precomputed layer
    layer_type      'image' or 'segmentation'
    data_type       'int16', 'float32', etc... (see neuroglancer source code: src/neuroglancer/util/data_type.ts for supported types)
    voxel_offset    x,y,z offset in voxels from the origin
    chunk_size      The size of each subvolume chunk in x,y,z.
    commit          if True, will write the info/provenance file to disk. 
                    if False, just creates it in memory
    """
    info = CloudVolume.create_new_info(
        num_channels = 1,
        layer_type = layer_type, # 'image' or 'segmentation'
        data_type = data_type, # 
        encoding = 'raw', # other options: 'jpeg', 'compressed_segmentation' (req. uint32 or uint64)
        resolution = resolution, # Size of X,Y,Z pixels in nanometers, 
        voxel_offset = voxel_offset, # values X,Y,Z values in voxels
        chunk_size = chunk_size, # rechunk of image X,Y,Z in voxels -- only used for downsampling task I think
        volume_size = volume_size, # X,Y,Z size in voxels
        )

    vol = CloudVolume(f'file://{layer_dir}', info=info)
    vol.provenance.description = description
    vol.provenance.owners = ['ahoag@princeton.edu'] # list of contact email addresses
    if commit:
        vol.commit_info() # generates info json file
        vol.commit_provenance() # generates provenance json file
        print("Created CloudVolume info file: ",vol.info_cloudpath)
    return vol

def upload_volume_slice(z,image_vol,cloud_vol,progress_dir):
    """ 
    ---PURPOSE---
    Upload a single z plane from a numpy volume to a cloud volume object.
    Made to be run in parallel.
    ---INPUT---
    z            0-indexed integer corresponding to the z plane to upload
    image_vol    3D numpy array from which you want to upload the slice
    cloud_vol    CloudVolume object to which you want to upload the slice 
    progress_dir Directory with list of files that get touched as slices are uploaded
    """ 

    print('Processing slice z=',z)
    y_dim,x_dim = image_vol[z].shape
    array = image_vol[z].reshape((1,y_dim,x_dim)).T 

    cloud_vol[:,:, z] = array
    touch(os.path.join(progress_dir, str(z)))
    return "Success"

def upload_volume_slice_coronal(y,image_vol,cloud_vol,progress_dir):
    """ 
    ---PURPOSE---
    Upload a single y plane from a numpy volume to a cloud volume object.
    Made to be run in parallel.
    ---INPUT---
    z            0-indexed integer corresponding to the z plane to upload
    image_vol    3D numpy array from which you want to upload the slice
    cloud_vol    CloudVolume object to which you want to upload the slice 
    progress_dir Directory with list of files that get touched as slices are uploaded
    """ 

    print('Processing slice y=',y)
    z_dim,x_dim = image_vol[:,y,:].shape
    array = image_vol[:,y,].reshape((z_dim,1,x_dim)).T  

    cloud_vol[:,y,:] = array
    touch(os.path.join(progress_dir, str(y)))
    return "Success"

def calculate_chunks(downsample, mip):
    """
    Chunks default to 64,64,64 so we want different chunks at 
    different resolutions
    """
    d = defaultdict(dict)
    d['full'][-1] = [1024,1024,1]
    d['full'][0] = [128,128,64]
    d['full'][1] = [128,128,64]
    d['full'][2] = [128,128,64]
    d['full'][3] = [128,128,64]
    d['full'][4] = [128,128,64]
    d['full'][5] = [64,64,64]
    d['full'][6] = [64,64,64]
    d['full'][7] = [64,64,64]
    d['full'][8] = [64,64,64]
    d['full'][9] = [64,64,64]

    try:
        result = d[downsample][mip]
    except:
        result = [64,64,64]
    return result

def calculate_factors(downsample, mip):
    """
    Scales get calculated by default by 2x2x1 downsampling
    """
    d = defaultdict(dict)
    d['full'][0] = [2,2,1]
    d['full'][1] = [2,2,2]
    d['full'][2] = [2,2,2]
    d['full'][3] = [2,2,2]
    d['full'][4] = [2,2,2]
    d['full'][5] = [2,2,2]
    d['full'][6] = [2,2,2]
    d['full'][7] = [2,2,2]
    d['full'][8] = [2,2,2]
    d['full'][9] = [2,2,2]

    try:
        result = d[downsample][mip]
    except:
        result = [2,2,1]
    return result

def determine_factors(chunk_size):
    # Determines how to calculate downsampling factors given an initial chunk size
    # The factors need to be defined in such a way that the chunks end up being 
    # maximally isotropic. 

    # Algorithm:
    # 1. Start by finding the largest number, N and its index in chunk_size, N_i. Downsample factor for N is D_N. Need to determine if 1 or 2.
    # 2. Then find the next largest number, M and its index in chunk_size, M_i. Downsample factor for N is D_M. Need to determine if 1 or 2.
    # 3. Now find the smallest number, P and its index in chunk_size, P_i. Downsample factor for P is D_P = 1
    # 4. Calculate X = round(N/M), Y = round(N/P).
    # 5. If X > 1, then let N = N/2. D_N = 2, D_M = 1. 
    #    Else, N and M are already within a factor of 2 and are isotropic. D_N=2, D_M=2. 
    #    If Y <= 1, N and P are already isotropic. Therefore chunk is isotropic -> D = [2,2,2]. Break out of while loop
    assert all([x>0 for x in chunk_size])
    factor_level_list = [[1,1,1]] # level 0 has [1,1,1] by default
    while True:
        D_N,D_M,D_P = 2,2,1 # initialize
        argsort = np.argsort(chunk_size)
        N_i = argsort[2]
        M_i = argsort[1]
        P_i = argsort[0]
        N = chunk_size[argsort[2]]
        M = chunk_size[argsort[1]]
        P = chunk_size[argsort[0]]
        X = round(N/M)
        Y = round(N/P)
        if X > 1:
            D_N, D_M = 2, 1
        elif Y <= 1: 
            D_P = 2
        factors = [0,0,0]
        factors[N_i] = D_N
        factors[M_i] = D_M
        factors[P_i] = D_P
        chunk_size = [0,0,0]
        chunk_size[N_i] = int(N/D_N)
        chunk_size[M_i] = int(M/D_M)
        chunk_size[P_i] = int(P/D_P)
        factor_level_list.append(factors)
        if factors == [2,2,2]:
            break
    # append a whole bunch more [2,2,2] levels
    for _ in range(20):
        factor_level_list.append([2,2,2])
    return factor_level_list
