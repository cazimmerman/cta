# General imports
import os,sys,glob,shutil,json
cwd = os.getcwd()
utils_fld = os.path.join(cwd,"utils")
sys.path.append(utils_fld)
from pipeline_utils import fast_scandir

# Select datasets to analyze
if __name__ == "__main__":
	sample_dir = sys.argv[1].strip().rstrip("/")
	imaging_request = sys.argv[2].strip().rstrip("/")
	output_rootpath = sys.argv[3].strip().rstrip("/")

	request_name,sample_name = sample_dir.split('/')[-2:]
	src_dir = os.path.join(sample_dir,imaging_request,
		'rawdata/resolution_3.6x')
	dst_dir = os.path.join(output_rootpath,request_name,
		sample_name,imaging_request,'rawdata/resolution_3.6x')

	if not os.path.exists(dst_dir):
		os.makedirs(dst_dir)
		print(f"Creating dst dir: {dst_dir}")

	# Link and rename image files. Need to check if files are in old
	# old folder convention e.g. Ex_642_Em_2
	# or new folder convention e.g. Ex_647_Em_680
	old_src_642 = os.path.join(src_dir,'Ex_642_Em_2_corrected')
	new_src_642 = os.path.join(src_dir,'Ex_647_Em_690_corrected')
	if os.path.exists(old_src_642):
		try:
			src_642 = fast_scandir(old_src_642)[-1]
		except:
			src_642 = old_src_642
	elif os.path.exists(new_src_642):
		try:
			src_642 = fast_scandir(new_src_642)[-1]
		except:
			src_642 = new_src_642

	# print(src_642)
	dst_642 = os.path.join(dst_dir,'Ex_642_Em_2_corrected')

	os.makedirs(dst_642,exist_ok=True)

	src_files_642 = sorted(glob.glob(src_642 + '/*tif'))
	n_src_files_642 = len(src_files_642)
	print(f"have {n_src_files_642} ch642 corrected planes")
	print("Sym linking Ch 642 files if not done already")
	for ii,src in enumerate(src_files_642):
		dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
		dst = os.path.join(dst_642,dst_basename)
		if not os.path.exists(dst):
			os.symlink(src,dst)
	## Now 488
	old_src_488 = os.path.join(src_dir,'Ex_488_Em_0_corrected')
	new_src_488 = os.path.join(src_dir,'Ex_488_Em_525_corrected')
	if os.path.exists(old_src_488):
		try:
			src_488 = fast_scandir(old_src_488)[-1]
		except:
			src_488 = old_src_488
	elif os.path.exists(new_src_488):
		try:
			src_488 = fast_scandir(new_src_488)[-1]
		except:
			src_488 = new_src_488

	dst_488 = os.path.join(dst_dir,'Ex_488_Em_0_corrected')
	os.makedirs(dst_488,exist_ok=True)

	print("Sym linking Ch 488 files if not done already")
	src_files_488 = sorted(glob.glob(src_488 + '/*tif'))
	for ii,src in enumerate(src_files_488):
		dst_basename = 'Z'+f'{ii}'.zfill(4)+'.tif'
		dst = os.path.join(dst_488,dst_basename)
		if not os.path.exists(dst):
			os.symlink(src,dst)
