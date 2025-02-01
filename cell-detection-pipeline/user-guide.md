**These instructions are specific for Princeton Neuroscience Institute users and will need to be modified for external use**
**This code base is still under development**
**This document was written by Robert Fetcho with help from Christopher Zimmerman**

**Setting Up**

The workflow for the pipeline consists of utilizing three software packages (pystripe, ClearMap2 and BrainPipe) alongside the pipeline code itself (designed to run on spock utilizing slurm jobs) and requires the setup of 3 conda environments (lightsheet, ClearMap\_spock and cellfinder) to do all of this.    
Notes: 

1) This guide assumes you will be working from your user folder within jukebox(cup) directory \- if you are working from another location, you may need to find/replace additional sections of path names.   
2) This guide tries to make things easier by directly copying working code/packages from Witten lab that have already been edited/debugged to work with the pipeline. If you want to install packages directly from source (ie: git clone) to get the latest versions, be warned that this may involve some non-trivial troubleshooting to get things running smoothly (see end of this document for some notes regarding what code was changed in these packages)  
3) Paths are hardcoded into the scripts \- TODO \- fix this?.  
4) For each “package” of the pipeline, there is a folder in lightsheet-pipeline with scripts and also a folder with logs for outputs of each script, very important for any necessary debugging.  
5) **It’s recommended all of this be done on “spock” as opposed to “scotty” (at minimum, the environments be installed while on spock.** As it seems there are some issues related to using scotty and clearmaps compilation of C code.

1\)  Log into spock and navigate to your jukebox folder (ie: cd /jukebox/LAB/YOUR\_NAME)

2\) Make a new folder for all of this called “lightSheetPipe” (ie: mkdir lightSheetPipe) and cd to that dir.

3\) Copy the lightsheet-pipeline code from Rob’s folder to your dir (cp /jukebox/witten/Rob/lightsheet-pipeline-GENERAL.zip .)

4\) Unzip  and move contents into current directory   
unzip lightsheet-pipeline-GENERAL.zip  
mv lightsheet-pipeline-GENERAL/\* .  
rm \-r lightsheet-pipeline-GENERAL  
rm lightsheet-pipeline-GENERAL.zip

**Checkpoint:**

You should have a folder (/jukebox/LAB/YOUR\_NAME/lightSheetPipe) containing the following folders: BrainPipe, clearmap-output, ClearMap2-master, lightsheet-pipeline, pystripe.

6\) Change hardcoded paths in lightsheet-pipeline code. In addition to files within the lightsheet-pipeline folder there are 4 folders containing .py and .sh files that require hardcoded path name changes \- cellfinder, clearmap, elastix and pystripe.  There are 2 “types” of hardcoded paths \- paths to the DATA and paths to the PIPELINE.  Using your favorite editor (ie Visual Studio Code, you can find/replace all instances across a workplace/folder so it is very fast), for *EVERY* .py and .sh file in these folders, try to find and replace the following (replace for all instances found, but not every file has both or any, but just check them all):

1) PIPELINE path change:  
   

	Find: YOUR\_DIR   
	Replace: (Your actual directory)

	Example:   
	The path in the code is: /jukebox/YOUR\_DIR/lightsheet-pipeline  
Rob’s pipeline is in /jukebox/witten/Rob/lightSheetPipe/lightsheet-pipeline  
	So he would find: YOUR\_DIR and replace with: witten/Rob/lightSheetPipe

2) DATA path change:

	Find: YOUR\_DATA  
	Replace: (your data directory)

Example: 

	The path in the code is:   
/jukebox/YOUR\_DATA/$request\_name/$request\_name-$sample\_name/….  
Rob’s data is in /jukebox/SmartSPIMData/lightserv\_temp/rf6456/…..  
	So he would find: YOUR\_DATA and replace: SmartSPIMData/lightserv\_temp/rf6456

**SUGGESTED: Spot check if paths in code make sense (very easy to have extra ‘/’ for example)**

7\) Change email address in .sh files for notifications related to jobs. Look in all of the .sh files and find/replace YOU with your princeton email address.

Find: YOU  
	Replace: (your princeton email username)

Example: 

	The line in the .sh file is:   
\#SBATCH –mail-user=YOUR\_EMAIL@princeton.edu  
Rob’s email is rob.fetcho@princeton.edu  
	So he would find: YOUR\_EMAIL and replace: rob.fetcho

8\) Create and set up lightsheet conda environment

1) Navigate to the pipeline folder containing the yml files:  
   cd /jukebox/LAB/YOUR\_NAME/lightSheetPipe/lightsheet-pipeline/environment-yml-files/

2) Activate the cluster’s anacondapy module:  
   module load anacondapy/2020.11  
     
3) Create a new conda environment:  
   conda env create \--name lightsheet \--file=lightsheet.yml  
     
4) Activate environment and install pystripe:  
   conda activate lightsheet  
   pip install https://github.com/chunglabmit/pystripe/archive/master.zip  
     
5) Deactivate environment for now:  
   conda deactivate  
   

9\) Create ClearMap\_spock environment

1) Assuming you are still in the folder with the yml files and the anacondapy module is active, create a new environment:  
   conda env create \-f ClearMap\_spock.yml  
   

10\) Create and set up cellfinder environment

1) Assuming the anacondapy module is still active, create a new environment using the cellfinder.yml  
   conda env create \--name cellfinder \--file=cellfinder.yml

**Running the pipeline**

The pipeline consists of 5 (with an optional “2.5th” command) individual .sh commands that need to be run sequentially.  Each command takes identical input \- the specific data sample folder hierarchy:  
	command.sh “request\_name” “imaging\_request\_number “ “sample\_name”

For example: Chris has a cohort called zimmerman\_01 and samples named 111, 112, 113. His data folder hierarchy is: 

/jukebox/SmartSPIMData/lightserv\_temp/cz15/request\_name/request\_name\-sample\_name/imaging\_request\_number/rawdata/resolution\_3.6x/…..

Thus, the data for Chris’ sample 112 from the zimmerman\_01 cohort exists in the following directory:

/jukebox/SmartSPIMData/lightserv\_temp/cz15/zimmerman\_01/zimmerman\_01\-112/imaging\_request\_1/rawdata/resolution\_3.6x/…..

So, when Chris calls his pipeline commands, it will look like this:  
	  
	command.sh zimmerman\_01 1 112

Wherever your imaging data resides, you should format it with that triple hierarchy so that things run smoothly.

Additionally, you may need to change permissions of the .sh files to run them on spock. To do this use chmod:  
	chmod 777 command.sh

With that in mind

Running the pipeline:

Note: All of these commands  should be run from your lightsheet-pipeline folder (ie: /jukebox/witten/YOUR\_NAME/lightSheetPipe/lightsheet-pipeline). The environments used by each script are indicated below for knowledge sake, but are activated within the code, so you do not need to be in an active environment when running these commands.

1\) Striping artifact removal (ENVIRONMENT: lightsheet):  
	  
	./pystripe\_pipeline.sh request\_name imaging\_request\_number sample\_name

	Success? Check logs or look into raw data folder for ‘corrected’ files.

2\) Atlas registration (ENVIRONMENT: ClearMap\_spock):

	./elastix\_pipeline.sh request\_name imaging\_request\_number sample\_name  
	  
	Success? clearmap-output for sample should now have “downsized” imaging files and an   
elastix\_inverse\_transform folder with stuff in it.

2.5) (OPTIONAL) Cell detection parameter setup:

	The lightsheet-pipeline should already have a correctly set up “cell\_detection\_parameter.p” file; however, if you need to regenerate or generate with different settings (edited in clearmap\_parameters.py file), run the following:

	./clearmap\_parameters.sh 

3\) Candidate cell detection (ENVIRONMENT: ClearMap\_spock):

	./clearmap\_pipeline.sh request\_name imaging\_request\_number sample\_name

	Success? clearmap-output folder for sample should have a large (like 400gb) stitched.npy file and   
a cells\_blocks folder with a bunch of .p files in it. Should have a cells\_transformed\_to\_atlas.npy file in your clearmap\_output \- note it seems like it’s possible the job “Fails” in slurm but actually succeeded \- check the log/check for this .npy file, if it got through registering cell detection results than it is fine. Note: In clearmap\_postprocessing\_register.py there is an initial filter   
based on cell size (line 61+) that is currently at 200 but can be altered potentially.

4\) Candidate cell classification (ENVIRONMENT: cellfinder):

	./cellfinder\_pipeline.sh request\_name imaging\_request\_number sample\_name

	Success? clearmap-output folder should have a cellfinder directory, at end of folders for current  
	sampe should be a “cell\_classification.xml” file.

5\) Cleanup and post-processing (ENVIRONMENT: ClearMap\_spock):

	./clearmap\_postcellfinder.sh request\_name imaging\_request\_number sample\_name

	Success? Clearmap-output should have output files such as the ones mentioned below.

**Useful Output Files:**

‘cells\_transformed\_to\_atlas\_classified\_….npy’ : contains the atlas coordinates (x,y,z), brain region, and estimated size for every Fos\+ cell

‘Region\_summary\_statistics\_classified\_….npy’ : contains the total count of Fos\+ cells for every brain region in the atlas

**Additional Info/Points:**

\- Pipeline assumes the raw imaging files were stitched already using the “terastitcher” package in python.  
\- Kernel density estimate (KDE) heatmaps via scipy.stats python code is included in the pipeline folders (lightsheet-pipeline/kde-visualization) but is not user friendly right now (TODO)  
\- One option for modeling this fos data is a mixed-effects negative binomial regression across conditions using the ‘GLMMTMB’ package in R  
\- The ‘stitched.npy’ file created in the first step of the clearmap\_pipeline.sh is quite large (400gb), so may be wise to delete for each sample after finishing processing for the cohort.

**Potential Failure Points**  
\- *Registration of autofluorescence channel to atlas (spim\_inverse\_register.py during elastix\_pipeline)* : Fix \- register Fos channel directly to atlas. Contact Chris (or Rob)  
\- *Preparing ClearMap output for Cellfinder classification* *(clearmap\_postprocessing\_register.py during clearmap\_pipeline.sh)*:  May happen especially if run before atlas registration is complete. Fix \- Run python file again by running clearmap\_postprocessing\_registeronly.sh. \*\* Note: Rob had an error on one sample during clearmap\_pipeline.sh during postprocessing but it was not because atlas was not done registering \- unclear why but error was a failure to import the “Cells” module from ClearMap \- rerunning the clearmap\_postprocessing\_registeronly.sh it worked on the second attempt.  
\- *Cellfinder classification times out (cellfinder\_inference.sh*) : Happens when too many cell candidates (ie: \> 5million). Fix \- break the cell candidates into two separate files and run Cellfinder twice. Contact Chris (or Rob) if this happens \- NOTE: This should not really happen anymore, the default behavior now is to split the cell candidates into 4 files and run in parallel.

**Main changes made to package code for pipeline (to help you get things running if you download packages fresh from github)**

1\) Changes to ClearMap package

1) In file ClearMap2/ClearMap/ImageProcessing/Experts/Cells.py:

	Add to imports: from scipy import ndimage  
	Replace Line 368/369 with:  
	background \= wrap\_step('background\_correction', ndimage.median\_filter(corrected,3), remove\_background,remove\_previous\_result=True, \*\*default\_step\_params)

2) In file ClearMap2/ClearMap/ImageProcessing/IlluminationCorrection/py

	  
	Comment out line 35 (\#import matplotlib.pyplot as plt)

3) In file ClearMap2/ClearMap/IO/Workspace.py

	  
	Comment out line 33 (\#import ClearMap.Visualization.Plot3d as q\_plot\_3d)

2\) Changes to BrainPipe package

1) In file BrainPipe/tools/imageprocessing/preprocessing.py:  
     
   Change line 11 to just import tifffile

 

TODO:

Implement clear method for fixing failed registration

Implement automatic re-try of clearmap\_postprocessing if it fails

Mergeblocks can fail randomly too (kind of light clearmap\_postprocessing) \- need to check/rerun if this happens. 

Combine all 5 steps into running of just one script with all batch jobs based on dependency.