# General imports
import pickle,sys

# ClearMap2 imports
sys.path.append('/jukebox/YOUR_DIR/ClearMap2-master')
import ClearMap.ImageProcessing.Experts.Cells as cells
import ClearMap.Utils.HierarchicalDict as hdict

# Original ClearMap paper (Renier et al., Cell, 2016)
# LaVision Ultramicroscope II
# 2x/0.5 NA objective at 0.8x zoom
# 4.06 x 4.06 x 4.5 µm per pixel

# PNI data
# Life Canvas SmartSPIM
# 3.6x/0.2 NA objective at 1.0x zoom
# 1.80 x 1.80 x 2.0 µm per pixel

# Axial scale factor: 2.26
# Voxel scale factor: 11.45

# Set ClearMap2 cell detection parameters
cell_detection_parameter = cells.default_cell_detection_parameter.copy()
cell_detection_parameter['iullumination_correction'] = None
cell_detection_parameter['background_correction']['shape'] = (21,21) # ClearMap paper: (7,7)
cell_detection_parameter['background_correction']['form'] = 'Disk'
cell_detection_parameter['background_correction']['save'] = False
cell_detection_parameter['equalization'] = None
cell_detection_parameter['dog_filter'] = None
cell_detection_parameter['maxima_detection']['h_max'] = None
cell_detection_parameter['maxima_detection']['shape'] = 11 # ClearMap paper: 5
cell_detection_parameter['maxima_detection']['threshold'] = None
cell_detection_parameter['maxima_detection']['valid'] = True
cell_detection_parameter['maxima_detection']['save'] = False
cell_detection_parameter['shape_detection']['threshold'] = None
cell_detection_parameter['shape_detection']['save'] = False
cell_detection_parameter['intensity_detection']['method'] = 'mean'
cell_detection_parameter['intensity_detection']['shape'] = None # ClearMap paper: 3
cell_detection_parameter['intensity_detection']['measure'] = ['background_correction']
cell_detection_parameter['verbose'] = True

# Save cell_detection_parameter.p
fname = '/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_parameter.p'
with open(fname, 'wb') as f:
    pickle.dump(cell_detection_parameter,f)
cell_detection_parameter = None

# Load and print cell_detection_parameter.p
with open(fname,'rb') as f:
    cell_detection_parameter = pickle.load(f)
print()
print('path                    : ' + fname)
hdict.pprint(cell_detection_parameter)
print()

# Set ClearMap2 cell detection filter thresholds
cell_detection_filter = {}

# PMA values
#cell_detection_filter['x'] = (0,352)
#cell_detection_filter['y'] = (59,510)
#cell_detection_filter['z'] = (0,540)
#cell_detection_filter['size'] = (400,None)
#cell_detection_filter['region'] = (2,1105)

# Allen CCFv3 values
cell_detection_filter['x'] = (0,320)
cell_detection_filter['y'] = (0,528)
cell_detection_filter['z'] = (0,456)
cell_detection_filter['size'] = (350,None)
cell_detection_filter['intensity'] = (0,None)
cell_detection_filter['region'] = (2,1028)

# Save cell_detection_filters.p
fname = '/jukebox/YOUR_DIR/lightsheet-pipeline/clearmap/cell_detection_filter.p'
with open(fname, 'wb') as f:
    pickle.dump(cell_detection_filter,f)
cell_detection_filter = None

# Load and print cell_detection_filters.p
with open(fname,'rb') as f:
    cell_detection_filter = pickle.load(f)
print()
print('path     : ' + fname)
hdict.pprint(cell_detection_filter)
print()
