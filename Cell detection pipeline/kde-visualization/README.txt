## README.txt

#### Installation/setup

## The conda environment you will need to run this pipeline can be installed by running:
## conda create --name <env> --file neuro_conda_env.txt
## where <env> is the name of the environment you want to use. I used "neuro" throughout this pipeline
## so if you change it to something else do a find and replace on neuro.
## Note that you might be able to use a different conda environment and just install the missing packages as needed manually using pip install

## You will need to copy the file:
/jukebox/wang/ahoag/pythonscripts/precomputed_utils.py
## to somewhere in your bucket folder and change the sys.path.append line at the top of:
scripts/make_heatmap_fullres_weighted_modkernel.py
## and the top of
scripts/make_cell_layers.py
## to point to the folder where you copied precomputed_utils.py
## This file has useful functions for making precomputed layers (neuroglancer-compatible data format)
## You may notice that the function: upload_volume_slice_coronal() from this file is also defined in
make_heatmap_fullres_weighted_modkernel.py
## That is because if this function is imported from an external file and used in a parallel process pool,
## the parallelization does not work! Not totally sure why but I think it has to do with sharing memory.

## Set the "atlas_dir" variable in:
make_heatmap_fullres_weighted_modkernel.py
## to the full path of the atlas_resources folder

#### Overview of what this pipeline does:

## Given an array of x,y,z coordinates from a treatment group of brains,
## the code creates a 3D kernel density estimate
## using scipy.stats.gaussian_kde(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
## The **creation** of the kernel is very fast even for millions of input coordinates.
## The **evaluation** of the kernel on the PMA grid is quite slow.
## The PMA grid has dimensions of (z,y,x): (540, 640, 352). Therefore the number of voxels is
## 540*640*352 = 121,651,200. A big fraction of these voxels are in the deadspace
## outside of the brain and can be ignored in the kernel evaluation.
## It turns out only around 41,000,000 voxels need to be evaluated.
## The first part of step 1 of the pipeline figures out the "good" voxels,
## i.e. the voxels inside the brain on which the kernel actually needs
## to be evaluated. Because the PMA never changes, this list of
## good voxels only needs to be computed once, so you
## will see some commented out lines in step 1 where I previously did this. You will
## never need to run those lines again unless you switch to a different atlas.
## The results from running this are saved in the folder atlas_resources/
## The second part of step 1 actually calculates the kernel so don't comment that part out.

## Even restricting to a smaller list of coordinates,
## evaluating the kernel is a pretty costly operation.
## Fortunately, it can be fully parallelized, and that is done via array jobs in step 2 of the pipeline (item 5 below).
## We define the chunk_size as the number of evaluations that happens in each array job
## If you use a chunk_size of 100000, you need 491 array jobs to evaluate
## the kernel on all of the good positions in the grid.
## Each chunk of 100000 evaluations takes between 45 minutes and 80 minutes (not sure why the large variance).
## You can submit all 491 array jobs at once via "sbatch --array=0-490 ... "
## (though they likely won't all run at once due to cluster resource restrictions)


## The main code that makes the kde is called:
scripts/make_heatmap_fullres_weighted_modkernel.py

## NOTE: The way I have been running this code is to use python on scotty
## or a computer with a bucket mount to do all steps except the kernel evaluations (step 2).
## For step 2, I use scripts/cfos_heatmap_pipeline.sh to submit sbatch array jobs to spock
## because this step requires heavy parallelization.
## If you really want to streamline the whole pipeline you could
## put all steps with dependencies into cfos_heatmap_pipeline, like what we do in the cell detection pipeline.
## I have some commented out code that could help you start to do that, though I don't have dependencies set up.

#### How to run the code:
## To make a cell heatmap for a group of brains in a single condition, do the following:

## 1)
## modify scripts/conditions.py. Make a conditions_dict whose keys are the names
## of the conditions you want to make heatmaps for. For each condition make a nested dictionary
## which has the key "brains" which maps to a list of dictionaries,
## one for each brain. Follow the format of the example dictionary in the file.

## 2)
## modify make_heatmap_fullres_weighted_modkernel.py to make sure input_dir, save_dir
## and atlas_dir point to the appropriate, existing folders.
## Also, this script is hardcoded to use the filtered_200px cell detection results.
## If you want to use a different filter size, you'll have to find and replace that in this file
## You will also have to create a new folder: output_filtered_{XYZ}px, where {XYZ} is the size of the new filter

## 3)
## Run step0 of the script with python for a single condition:
$ cd scripts
$ module load anacondapy/2020.11
$ conda activate neuro
$ python make_heatmap_fullres_weighted_modkernel.py step0 $condition_name 100000
## where $condition_name is one of the conditions that you listed in conditions.py
## This will create a CSV file in your save_dir containing all of the coordinates
## and weights from each of the brains in the condition cohort, based on what you provided
## in your conditions.py dictionary.
## The weights for each coordinate are 1/n_points_i where n_points_i is the
## number of coordinates in the brain from which the coordinate came.

## 4)
## Run step1 of the script with python for the same single conditon
$ module load anacondapy/2020.11 # don't need to do this if you already did it before
$ conda activate neuro # don't need to do this if you already did it before
$ python make_heatmap_fullres_weighted_modkernel.py step1 $condition_name 100000
## this generates the kernel for the list of coordinates and weights from all brains in this condition.

## 5)
## Run step 2 using cfos_heatmap_pipeline.sh
## modify cfos_heatmap_pipeline.sh to use the correct condition_name.
## I suggest keeping chunk_size at 100000. If you modify it, note that you will need
## more or less array jobs to complete the kernel evaluation.
## If you need to debug it, bring chunk_size down to 1000 or so and only run a few array jobs
## this can often help you diagnose the problem quickly without waiting hours.

## Note that in the python code of step2 is where I shrink the kernel size with this line:
>> kernel.set_bandwidth(kernel.factor*0.25)
## This does not normalize relative to other kernels, it just shrinks it to 1/4 of its original size.
## You could hardcode a bandwith here instead, which it sounds like what you might want to do to
## normalize the kernel size across conditions. Note that kernel.factor is just a number.
## The kernel covariance matrix sets the shape of the kernel
## (i.e. the elliptical gaussian's x,y,z standard deviations and its major axis angle)
## This can be inspected using:
>> kernel.covariance
## Unfortunately, the covariance cannot be modified in this implementation.
## I included a jupyter notebook:
notebooks/cell_heatmap_density_sandbox.ipynb
## which shows you how to load the kernel in python and see its properties like the
## bandwidth and the covariance. There is a lot of mess in there too, sorry.
## The scipy docs for gaussian_kde() are pretty good.


## 6)
## Run step 3 of the pipeline using python.
## Step 2 of the pipeline saved the kernel evalulations from each of the 491 chunks in separate files.
## Step 3 combines these into a single 3D array the same shape as the PMA,
## making sure to apply the evalulations to the "good" positions inside the brain.
## Outside of the brain the value of the kde is set to 0.
## The resulting final kde file is saved as {save_dir}/kde_20um_{condition_name}_weighted_decr75pc.npy

## 7)
## Run step 4 of the pipeline using python.
## This step makes a precomputed layer (Neuroglancer-compatible data format) of the heatmap
## and saves it to a place on bucket that Neuroglancer can see.

## 8)
## Re-run steps 3-7 in this README for a different condition you want to compare to your first condition

## 9)
## Depending on whether you want to make a log ratio or a simple difference map between the two conditions,
## run steps 5 or 6 in the main .py file, respectively, using python.
## This will create a precomputed layer for the log ratio or difference map, respectively.

## 10)
## Finally, if you want to see the cells from all brains in a given condition,
## run:
$ python make_cell_layers.py step0 $condition_name # not necessary to run if you have run step0 of the main pipeline
$ python make_cell_layers.py step1 $condition_name
## This will create a precomputed layer for the cells. Note that you will see all cells equally in the same color,
## regardless of what brain they came from. This script could be modified to make individual brain layers
## and then wyou could color them differently in neuroglancer using the GUI tools.

## 11)
## To make yourself a neuroglancer link which has the 3 panel layout with control, treatment,
## and treatment-control all with atlas boundaries overlaid in the Reds, Blues and RdBu colormaps like in the demo I sent you (the url for that demo is in the file scripts/neuroglancer_template_url.txt), run the code:
$ python neuroglancer_link_maker.py $control_name_short $treatment_name_short
## where $control_name_short and $treatment_name_short are the names you want displayed in the top layer panels.
## They do not have to be identical to the condition names in conditions.py (see the demo for example)
## This code will read your conditions.py file and use the first condition as the control condition
## and the second condition as the treatment condition and generate you a new link with the proper
## layers loaded in and with the names for the layers you supplied as command line arguments.
## Note that the link is set up to show the difference between condition 2 and condition 1
## in the third panel. If you want the log ratio, you will have to manually make that link.

## If the link maker script does not work or you want to try to load layers manually first go to:
https://neuroglancer-demo.appspot.com/
## Then find your layer of interest here:
https://lightsheetatlas.pni.princeton.edu/public/lightserv/cz15/cfos_heatmaps/filtered_200px/
## copy the url to the layer of interest and then load it into neuroglancer using the source box
## select or type in the precomputed:// source type and then enter the URL to the layer and click enter
## Once it loads, there should be a text box appearing at the bottom right of the screen saying
## "create as image layer" for image layers and saying something similar for other types of layers.
## Click that text box or simply hit enter again and it should load the layer. To get the colormap right,
## use the demo link in neuroglancer_template_url.txt and extract the code in the "shader" box
## in the "Rendering" tab of the image layers and paste that into the shader box of your new layer
## You may need to adjust the histogram in the Rendering tab to get the signal to show up properly.
