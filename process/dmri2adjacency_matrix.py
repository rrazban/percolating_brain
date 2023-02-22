"""
Process the diffusion MRI data from the UK Biobank, 
ABCD study, or dHCP into a weighted adjacency matrix
where elements represent tract length or density.

Closely follows the DIPY tutorial for generating 
the connectivity matrix
https://dipy.org/documentation/1.3.0./examples_built/streamlines_analysis/streamline_tools/

Note that to run this script, you must first acquire
dMRI scans from the respective dataset.

"""


import numpy as np
from scipy.ndimage import binary_dilation

from dipy.core.gradients import gradient_table
from dipy.data import get_fnames
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti_data
from dipy.direction import peaks
from dipy.reconst import shm
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion
from dipy.tracking.streamline import Streamlines
from dipy.tracking.streamline import length
from dipy.reconst.dti import TensorModel


from nilearn import image

from glob import glob 
from multiprocessing import Pool
import shutil, zipfile, os, pickle, tarfile 
from functools import partial
import datetime



def get_remaining_abcd(fnames, pre_done_already):
    done_already = []
    for filename in pre_done_already:
        pre_eid = filename.split('/')[-1]
        first = pre_eid[pre_eid.index('-')+1:pre_eid.index('_')]
        new_eid = pre_eid[pre_eid.index('_')+1:]
        second = new_eid[new_eid.index('-')+1:new_eid.index('_')]
        fname='{0}_{1}'.format(first, second)
        done_already.append(fname)
    
    not_present = []
    for x in fnames:
        almost = x.split('/')[-1]
        firstTest = almost.find('_')
        fname = almost[:almost.find('_', firstTest+1)]
        if fname not in done_already:
            not_present.append(x)
    return not_present




def get_dMRI(dataset):  #have following directories correctly set to location of output files and input files
    done_already = glob('/shared/datasets/public/{0}/derivatives/*_density.txt'.format(dataset))   #output data files
    done_already = [fname.split('/')[-1].replace('_density.txt', '') for fname in done_already]

    if dataset=='ukb':
        fnames = glob('/shared/datasets/public/ukb/dti_warped/*')   #dMRI data files 
        remaining = [fname for fname in fnames if fname.split('/')[-1].replace('.nii.gz', '') not in done_already]
    elif dataset=='abcd':
        fnames = glob('/shared/datasets/public/abcd/AWS_downloads/fmriresults01/*.tgz')
        remaining = get_remaining_abcd(fnames, done_already)
    elif dataset=='dhcp':
        fnames = glob('/shared/datasets/public/dhcp/rel3_dhcp_dmri_shard_pipeline/*/*/dwi/*dwi.nii.gz') 
        remaining = [fname for fname in fnames if fname.split('/')[-1].replace('.nii.gz', '') not in done_already]

    return remaining

def get_atlas_info():
    label_fname = 'atlas/Talairach.nii.gz'
    white_matter_label = 14	#14 for Talairach!! #1 for HarOx and DesKi

    print(label_fname)
    print('white matter label: {0}'.format(white_matter_label))
    labels = load_nifti_data(label_fname)
    labels_img = image.load_img(label_fname)
    return labels, labels_img, white_matter_label

def restructure(tract_structure, data):
    upps = np.triu(tract_structure)	#includes diagonal, set k =1 to not include
    indi_x, indi_y = np.nonzero(upps)

    fresh_structure = np.zeros((len(tract_structure), len(tract_structure)))    #maybe i can do this as a matrix operation to speed things up? nope
    for x,y in zip(indi_x, indi_y):
            tract_lengths = length(data[x,y])
            fresh_structure[x,y] = np.mean(tract_lengths)

    fresh_structure += fresh_structure.T 
    return fresh_structure

def readin_ukb(fname):
  subj = fname.split('/')[-1].replace('.nii.gz', '.zip')
  compressed = '/shared/datasets/public/ukb/dti_raw/%s' % subj
  uncompressed = '/scratch/rostam/%s' % subj

  if os.path.exists(uncompressed)==False:
    with zipfile.ZipFile(compressed, 'r') as zip_ref:
      zip_ref.extractall(uncompressed)	
									
  bval_fname, bvec_fname = \
    "%s/dMRI/dMRI/bvals" % uncompressed,\
    "%s/dMRI/dMRI/bvecs" % uncompressed	

  return fname, bval_fname, bvec_fname, uncompressed, subj

def readin_abcd(fnamegz):
  uncompressed = '/scratch/rostam/'
  tar = tarfile.open(fnamegz, "r:gz")
  tar.extractall(path=uncompressed)
  tar.close()

  almost = fnamegz.split('/')[-1] 

  firstTest = almost.find('_')
  sub = almost[:firstTest]
  ses = almost[firstTest+1:almost.find('_', firstTest+1)]

  DIR = "{0}/sub-{1}/ses-{2}/dwi".format(uncompressed, sub, ses)
  fname = glob('{0}/*.nii'.format(DIR))[0]

  subj = fname.split('/')[-1].replace('.nii', '.zip')
  
  bval_fname = fname.replace('.nii', '.bval')
  bvec_fname = fname.replace('.nii', '.bvec')
  scratch_file = "{0}/sub-{1}".format(uncompressed, sub)

  return fname, bval_fname, bvec_fname, scratch_file, subj

def readin_dhcp(fname):
  subj = fname.split('/')[-1].replace('.nii.gz', '.zip')
  
  bval_fname = fname.replace('.nii.gz', '.bval')
  bvec_fname = fname.replace('.nii.gz', '.bvec')
  scratch_file = 'none'

  return fname, bval_fname, bvec_fname, scratch_file, subj





def compute_features(fname, labels, labels_img, white_matter_label, dataset):
  if dataset=='ukb': 
    fname, bval_fname, bvec_fname, scratch_file, subj = readin_ukb(fname)
  elif dataset=='abcd':
    fname, bval_fname, bvec_fname, scratch_file, subj = readin_abcd(fname)
  elif dataset=='dhcp':
    fname, bval_fname, bvec_fname, scratch_file, subj = readin_dhcp(fname)


  img = image.load_img(fname) # warped
  img = image.resample_to_img(img, labels_img)
  data, affine = img.get_fdata(), img.affine

  bvals, bvecs = read_bvals_bvecs(bval_fname, bvec_fname)

  if dataset=='dhcp':   #nothing needs to be removed for dhcp since no tar files
      pass
  else:
      shutil.rmtree(scratch_file)  

  gtab = gradient_table(bvals, bvecs)
  white_matter = binary_dilation((labels == white_matter_label))

 # ten_model = TensorModel(gtab)  #can set csamodel = ten_model to check out dti method
#  csamodel=ten_model
  csamodel = shm.CsaOdfModel(gtab, 6)   #default
  csapeaks = peaks.peaks_from_model(model=csamodel,
                                    data=data,
                                    sphere=peaks.default_sphere,
                                    relative_peak_threshold=.8,
                                    min_separation_angle=45,
                                    mask=white_matter)


  affine = np.eye(4)
  seeds = utils.seeds_from_mask(white_matter, affine, density=1)

  stopping_criterion = BinaryStoppingCriterion(white_matter)

  streamline_generator = LocalTracking(csapeaks, stopping_criterion, seeds,
                                       affine=affine, step_size=0.5)

  streamlines = Streamlines(streamline_generator)
  streamlines = [sl for sl in streamlines if len(sl)>1]    #necessary when density>1. some streamlines do not go anywhere   #also rarely occurs for density=1

  M, grouping = utils.connectivity_matrix(streamlines, affine,
                                          labels.astype(np.int16),	
                                          return_mapping=True,
                                          mapping_as_streamlines=True)

  remove = [0,white_matter_label]
  for r in remove:
     M[r, :] = 0
     M[:, r] = 0

  out_dir = '/shared/datasets/public/{0}/derivatives'.format(dataset)   #have this set to where you want output to go 
  np.savetxt('{0}/{1}'.format(out_dir, subj.replace('.zip', '_density.txt')), M)
  M_lengths = restructure(M, grouping)
  np.savetxt('{0}/{1}'.format(out_dir, subj.replace('.zip', '_length.txt')), M_lengths)
  print('{0} is complete'.format(subj), flush=True)
#  out_fname = '{0}/{1}'.format(out_dir, subj.replace('.zip', '_grp.pkl'))  #gives all the tracts, just in case i want to do something else besides density or length #large memory
#  fhandler = open(out_fname, 'wb')      
#  pickle.dump(grouping, fhandler)   #hard to load this in parallel during analysis, hence i make the length file right now
#  fhandler.close()



if __name__ == '__main__':
    dataset = 'ukb' #ukb, abcd, dhcp

    dMRIs = get_dMRI(dataset)#[:5]
    print("Number of jobs: {0}".format(len(dMRIs)), flush=True) 
    start_time = datetime.datetime.now()
    print("Start time: {0}".format(start_time), flush=True)
    labels, labels_img, white_matter_label = get_atlas_info()

    for dMRI in dMRIs:
        compute_features(dMRI, labels, labels_img, white_matter_label, dataset) 

    #speed up processing by running in parallel
#    with Pool(5) as p:
 #       prod_x = partial(compute_features, labels=labels, labels_img=labels_img, white_matter_label=white_matter_label, dataset=dataset)
  #      p.map(prod_x, dMRIs)

    end_time = datetime.datetime.now()
    print("End time: {0}".format(end_time))
    print("Duration of processing: {0}".format(end_time-start_time))
