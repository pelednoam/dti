#!/usr/bin/env python
from __future__ import division
import os

def usage():
    print ('usage: mask_rois.py -s subject -m modality [--skip n] [--nomorph]\n'
        '\nstage 0: organize files for single shell reconstruction only'
        '\nstage 1: setup registration from anatomical to B0 space'
        '\nstage 2: morph white matter mask to B0 space (requires stage 1)'
        '\nstage 3: morph labels to subject and extract to B0 space'
        '\nstage 4: create mask in B0 space'
        '\nstage 5: run dsistudio pipeline and convert to alphabetical indices'
        '\nstage 7: do linearly corrected matrix creation')

#####################
def organize_singleshell(subj,preproc_dir,shell):
    print 'preparing file conversions for singleshell analysis at b=%s'%shell 
    import numpy as np
    import nibabel as nib
    import subprocess

    stem = 'diff_disco_eddy'
    new_stem = 'diff_singleshell_%s'%shell
    sample = 69 
    btable_file = 'btable.txt'
    new_btable_file = 'btable_singleshell_%s.txt'%shell

    if shell=='1k':
        ixes = np.arange(sample)
    elif shell=='3k':
        ixes = np.arange(sample, 2*sample)
    elif shell=='5k':
        ixes = np.arange(2*sample, 4*sample)
    elif shell=='10k':
        ixes = np.arange(4*sample, 8*sample)

    if os.path.exists(os.path.join(preproc_dir,'%s.nii'%stem)):
        gzip_cmd = ('mri_convert %s %s'%
            (os.path.join(preproc_dir,'%s.nii'%stem),
             os.path.join(preproc_dir,'%s.nii.gz'%stem)))
        print gzip_cmd
        m=subprocess.call(gzip_cmd,shell=True)

    a = nib.load(os.path.join(preproc_dir,'%s.nii.gz'%stem))
    img = nib.Nifti1Image(a.get_data()[:,:,:,ixes], a.get_affine(),
        a.get_header())
    nib.save(img, os.path.join(preproc_dir,'%s.nii.gz'%new_stem))

    b = np.loadtxt(os.path.join(preproc_dir,btable_file))
    np.savetxt(os.path.join(preproc_dir,new_btable_file), b[ixes])

    print 'done with file conversions for singleshell analysis'

######################
def setup_mprage(subj,nmr_dir,dwi_dir,preproc_dir,dset_dir,modality):
    import subprocess
    fsrecon_dir = os.path.join(dwi_dir,subj,'fsrecon')
    try:
        os.makedirs(dset_dir)
    except OSError:
        pass

    mprage_img = os.path.join(fsrecon_dir,'mri','orig.mgz')
    mprage_256 = os.path.join(dset_dir,'mprage_256.nii.gz')
    mprage_nii = os.path.join(dset_dir,'mprage.nii.gz')

    b0_img = os.path.join(dset_dir, 'b0.nii.gz')
    nodif_brain = os.path.join(dset_dir, 'nodif_brain.nii.gz')

    convert_cmd = 'mri_convert %s %s' % (mprage_img, mprage_256)
    print convert_cmd
    p=subprocess.call(convert_cmd,shell=True)

    if modality in ('dsi'):
        print preproc_dir
        b0conv_cmd = 'mri_convert %s %s' % (
            os.path.join(preproc_dir,'b0.nii'),
            b0_img)
        q=subprocess.call(b0conv_cmd,shell=True)

    elif modality in ('qball','qbi','singleshell','sgsh'):
        try:
            os.link( os.path.join(preproc_dir, 'b0_first.nii.gz'),
                     b0_img)
            print 'ln preproc/diff/b0_first.nii.gz b0.nii.gz'
        except OSError:
            pass
    elif modality == 'mask_only':
        pass

    #reorient_cmd = 'fslreorient2std %s %s' % (mprage_nii, mprage_nii)
    #print reorient_cmd
    #q=subprocess.call(reorient_cmd,shell=True)

    mprage_brain = os.path.join(dset_dir,'mpr_brain.nii.gz')

    bet_cmd = 'bet %s %s -f 0.40 -B' % (mprage_nii, mprage_brain)
    print bet_cmd
    r=subprocess.call(bet_cmd,shell=True)

    bet2_cmd = 'bet %s %s -f 0.15 -g 0 -m' % (b0_img, nodif_brain) 
    print bet2_cmd
    u=subprocess.call(bet2_cmd,shell=True)

    upsample_cmd = ('mri_convert -vs 1 1 1 %s %s' %
        ( b0_img,
          os.path.join(dset_dir, 'b0d.nii.gz')))

    matflirt_cmd = ('flirt -in %s -ref %s -omat %s' %
        ( mprage_256, b0_img,
          os.path.join(dset_dir, 'mpr2562b0.mat')))
    s=subprocess.call(matflirt_cmd,shell=True)

    flirt_cmd = ('flirt -in %s -ref %s -applyxfm -init %s -out %s' %
        ( mprage_256, b0_img,
          os.path.join(dset_dir, 'mpr2562b0.mat'),
          os.path.join(dset_dir, 'mprage_b0.nii.gz')))
    t=subprocess.call(flirt_cmd,shell=True)

    print 'Done with mprage reorientation and masking'

#####################
def morph_and_extract_labels(subj,nmr_dir,subjdir,parclist_dir,masks_dir,
    nomorph=False):
    import subprocess

    os.chdir(subjdir)

    parc_stem = 'wmreg'

    try:
        os.symlink( os.path.join(nmr_dir,'%s_FS'%subj),
                    os.path.join(subj,'fsrecon'))
    except OSError:
        pass

    morph_cmd = ('mri_label2label --srcsubject fsaverage5c '
        '--srclabel fsaverage5c/label/%s/%s.%s --trgsubject %s/fsrecon '
        '--trglabel %s/fsrecon/label/%s/%s.%s --regmethod surface --hemi %s')

    #extract_cmd = ('mri_label2vol --label %s/fsrecon/label/%s/%s.%s.label '
    #	'--temp %s/dset/mprage.nii.gz --surf white '
    #	'--o %s/masks/%s/%s-%s.nii.gz --regheader %s/fsrecon/mri/orig.mgz')

    #extract_cmd = ('mri_label2vol --label %s/fsrecon/label/%s/%s.%s.label '
    #	'--temp %s/dset/mprage_b0.nii.gz --surf white '
    #	'--o %s/masks/%s_b0/%s-%s.nii.gz --regheader %s/fsrecon/mri/orig.mgz')

    extract_cmd = ('mri_label2vol --label %s/fsrecon/label/%s/%s.%s.label '
        '--temp %s/dset/mprage_256.nii.gz --surf white '
        '--o %s/masks/%s_%s/%s-%s.nii.gz --regheader %s/fsrecon/mri/wm.mgz')

    flirt_cmd = ('flirt -in %s/masks/%s_%s/%s-%s.nii.gz -ref %s/dset/b0.nii.gz '
        '-applyxfm -init %s/dset/mpr2562b0.mat '
        '-out %s/masks/%s_%s/%s-%s.nii.gz')

    for parc in ('laus250',):
    #for parc in ('laus125','laus250','laus500'):	
        try:
            os.makedirs(os.path.join(masks_dir,'%s_%s'%(parc,parc_stem)))
        except OSError:
            pass

        try:
            os.mkdir(os.path.join(subj, 'fsrecon', 'label', parc))
        except OSError:
            pass

        entries=[]
        entry_file = os.path.join(parclist_dir,'%s_entries'%parc)
        with open(entry_file,'r') as fd:
            for ln in fd:
                entries.append(ln.strip())

        for hemi in ('lh','rh'):
            for entry in entries:
                'Now morphing and extracting %s'%subj
                m_cmd = morph_cmd % (parc,hemi,entry,subj,subj,parc,hemi,
                    entry,hemi)
                if not nomorph:
                    print m_cmd
                    p = subprocess.call(m_cmd,shell=True)

                print subj,parc,hemi,entry

                e_cmd=extract_cmd % (subj,parc,hemi,entry,subj,subj,parc,
                    parc_stem,entry,hemi,subj)
                print e_cmd
                q = subprocess.call(e_cmd,shell=True)

                f_cmd=flirt_cmd % (subj,parc,parc_stem,entry,hemi,subj,subj,
                    subj,parc,parc_stem,entry,hemi)
                print f_cmd
                r = subprocess.call(f_cmd,shell=True)

    print 'Done with morphing and extracting'

#####################
def bigmask_creation(subj,subjects_dir,masks_dir):
    import nibabel as nib
    import numpy as np

    #mangle_hemi = lambda x:'%s_%s'%(x[-9:-7],x[:-10])
    #demangle_hemi = lambda x:'%s-%s.nii.gz'%(x[3:],x[:2])

    #parc_orig = 'laus125'
    #parc = '%s_nonoverlapping' % parc_orig
    for parc in ( 'laus250_wmreg', ):

        #masks_dir=os.path.join(masks_top_lev_dir,parc)
        centroids_dir=os.path.join(masks_dir,parc)

        #copy parc_orig to parc if parc doesnt exist yet
        if not os.path.exists(masks_dir):
            raise ValueError('No such parcellation')
            #os.makedirs(masks_dir)
            #import shutil
            #shutil.copytree(centroids_dir,masks_dir,symlinks=True)

        mask_template=nib.load(os.path.join(centroids_dir,
            'bankssts_1-lh.nii.gz'))

        mask_shape = mask_template.shape
        big_mask=np.zeros((mask_shape[0],mask_shape[1],mask_shape[2]),dtype=int)
        binary_mask=big_mask.copy()

        print 'building mask in diffusion space with shape {0}'.format(mask_shape)

        #create centroids
        print "collecting centroids"
        centroids={}

        #make sure dirents are in alphabetical order
        dirents = []
        for hemi_path in ('lh.nii.gz','rh.nii.gz'):
            for dirent in sorted(os.listdir(centroids_dir)):
                if dirent.endswith(hemi_path):
                    dirents.append(dirent)

        for i,dirent in enumerate(dirents):
            print 'DBG: collecting centroid for subsequent mask', i, dirent
            try:
                r=nib.load(os.path.join(centroids_dir,dirent)).get_data()
                r[np.where(r)]=1
                centroids.update({i:np.mean(np.where(r),axis=1)})
                big_mask+=r*(i+1)
                binary_mask+=r
            except nib.spatialimages.ImageFileError:
                continue

        #binary_mask[binary_mask>0]=1
        print "all centroids collected"

        xi,yi,zi=np.where(binary_mask>1)

        for x,y,z in zip(xi,yi,zi):
            #print x,y,z
            vox=np.array((x,y,z))
            closest_centroid=0
            dist=np.inf
            nr_collisions=2	#divisor starts at 2
                            #ensuring that a random voxel will be chosen
            for i,centroid in centroids.items():
                #print vox,centroid
                cent_dist=np.linalg.norm(vox-centroid)
                if cent_dist<dist:
                    dist=cent_dist
                    closest_centroid=i+1

                if cent_dist==dist:
                    if np.random.random() < 1./nr_collisions:
                        dist=cent_dist
                        closest_centroid=i+1
                    nr_collisions+=1
                    
            big_mask[x,y,z]=closest_centroid

        img=nib.Nifti1Image(big_mask,mask_template.get_affine(),
            mask_template.get_header())
        nib.save(img,os.path.join(masks_dir,'bigmask_%s.nii.gz'%parc))

        print '%s bigmask saved'%parc

def white_matter_mask(subj,subjdir,dset_dir):
    import subprocess

    convert_cmd = ('mri_convert %s %s' % (
        os.path.join( subjdir, subj, 'fsrecon', 'mri', 'wm.mgz'),
        os.path.join( dset_dir, 'wm_256.nii.gz')))
    print convert_cmd
    p=subprocess.call(convert_cmd,shell=True)

    flirt_cmd = ('flirt -in %s -ref %s -applyxfm -init %s -out %s' % (
        os.path.join(dset_dir, 'wm_256.nii.gz'),
        os.path.join(dset_dir, 'b0.nii.gz'),
        os.path.join(dset_dir, 'mpr2562b0.mat'),
        os.path.join(dset_dir, 'wm_b0.nii.gz')))
    print flirt_cmd
    q=subprocess.call(flirt_cmd,shell=True)

def run_dsistudio(subj, subjdir, masks_dir, preproc_dir, dset_dir, modality,
        shell='10k', nr_fibers=20000, trk_output='tracts20k', 
        suppress_output=False, preproc_stem=None):
    import subprocess

    if modality in ('qball','qbi'):
        stem = 'diff_disco_eddy'
        methodnr = '4'
        methodname = 'gqi'
        param0 = '1.25'
        param1 = None
        btable_string = 'btable.txt'
    elif modality in ('singleshell','sgsh'):
        stem = 'diff_singleshell_%s'%shell
        methodnr = '3'
        methodname = 'qbi'
        param0 = '0.006'
        param1 = '8'
        btable_string = 'btable_singleshell_%s.txt'%shell
        trk_output = '%s_sgsh_%s'%(trk_output, shell)
    elif modality in ('dsi',):
        stem = 'data_raw_disco_eddy-vol_moco_clean_ordered'
        methodnr = '0'
        methodname = 'dsi'
        param0 = '17'
        param1 = None
        btable_string = 'dsi515_b_table.txt'
    elif modality in ('dti',):
        #stem = 'dti_unprocessed'
        stem = 'dti'
        methodnr = 1 
        methodname = 'dti'
        param0 = None
        param1 = None
        btable_string = 'btable.txt'
    else:
        sys.stderr.write('Unrecognized modality %s\n'%modality)

    if preproc_stem is not None:
        stem = preproc_stem

    try:
        os.mkdir(preproc_dir)
    except OSError:
        pass

    if os.path.exists(os.path.join(preproc_dir,'%s.nii'%stem)):
        gzip_cmd = ('mri_convert %s %s'%
            (os.path.join(preproc_dir,'%s.nii'%stem),
             os.path.join(preproc_dir,'%s.nii.gz'%stem)))
        print gzip_cmd
        m=subprocess.call(gzip_cmd,shell=True)

    if modality in ('qball','qbi','singleshell','sgsh'):
        src_cmd=('dsi_studio --action=src --source=%s --b_table=%s --output=%s'
            % (
            os.path.join(preproc_dir,'%s.nii.gz' % stem),
            os.path.join(preproc_dir,btable_string),
            os.path.join(preproc_dir,'%s.src.gz' % stem)))
    elif modality in ('dsi',):
        src_cmd=('dsi_studio --action=src --source=%s --b_table=%s --output=%s'
            % (
            os.path.join(preproc_dir,'%s.nii.gz' % stem),
            os.path.join(subjdir, btable_string),
            os.path.join(preproc_dir,'%s.src.gz' % stem)))
    elif modality in ('dti',):
        src_cmd=('dsi_studio --action=src --source=%s --output=%s --b_table=%s' % (
            os.path.join(preproc_dir, 'dti.nii.gz'),
            os.path.join(preproc_dir, '%s.src.gz' % stem),
            os.path.join(preproc_dir,btable_string)
))
    print src_cmd
    p=subprocess.call(src_cmd,shell=True)

    rec_cmd=('dsi_studio --action=rec --source=%s --thread=10 --method=%s '
        '--num_fiber=8 --odf_order=8 --record_odf=1 --mask=%s' % (
        os.path.join(preproc_dir,'%s.src.gz' % stem),
        methodnr,
        os.path.join(dset_dir,'wm_b0.nii.gz')))
    if param0 is not None:
        rec_cmd += ' --param0=%s'%param0
    if param1 is not None:
        rec_cmd += ' --param1=%s'%param1
    print rec_cmd
    q=subprocess.call(rec_cmd,shell=True)

    trk_cmd=('dsi_studio --action=trk --source=%s --fiber_count=%i '
        '--output=%s' % (
        os.path.join(preproc_dir, 
            '%s.src.gz.%s%s.%sfib.gz' % (
            stem, 
            '' if methodname=='dti' else 'odf8.f8rec.',
            methodname, 
            reduce(lambda x,y:x+y,['sh%s.'%p if i==0 else'%s.'%p for i,p in
                enumerate((param0,param1)[::-1])if p is not None],'') )),
            nr_fibers,
        #os.path.join(subjdir,subj,'%s.trk'%trk_output),
        os.path.join(dset_dir, 'tracks','%s.%s'%(trk_output, '%s'))))
    print trk_cmd
    r1=subprocess.call(trk_cmd%'txt',shell=True)
    r2=subprocess.call(trk_cmd%'trk',shell=True)

    for parc in ('laus250',):

        ana_cmd=('dsi_studio --action=ana --source=%s --tract=%s '
            '--export=connectivity --roi=%s' % (
            os.path.join(preproc_dir,
                '%s.src.gz.%s%s.%sfib.gz' % (
                stem, 
                '' if methodname=='dti' else 'odf8.f8rec.',
                methodname,
                reduce(lambda x,y:x+y,['sh%s.'%p if i==0 else'%s.'%p for i,p in
                    enumerate((param0,param1)[::-1])if p is not None],'') )),
            #os.path.join(subjdir,subj,'%s.trk'%trk_output),
            os.path.join(dset_dir, 'tracks','%s.txt'%trk_output),
            os.path.join(masks_dir,'bigmask_%s_wmreg.nii.gz'%parc)))
        print ana_cmd
        s=subprocess.call(ana_cmd,shell=True)

        #==================================================================
        #put the resulting matrix back in alphabetical order

        from scipy import io as sio
        import numpy as np

        #fname=os.path.join(subjdir,subj, '%s.trk.connectivity.mat'%trk_output)
        fname=os.path.join(dset_dir, 'tracks', '%s.txt.connectivity.mat'%trk_output)
        
        adj = sio.loadmat(fname)
        labnam = adj['name']
        med_length = adj['tract_median_length']

        adj = adj['connectivity']

        adj = adj/med_length
        adj[np.isnan(adj)]=0

        labnam = ''.join(map(chr,labnam.flat)).split('\n')[:-1]

        stem = len(os.path.commonprefix(labnam))
        ord = np.array( map( lambda x:int(x[stem:])-1, labnam))
        #dsistudio indexing to python indexing

        keys={}
        for i,k in enumerate(ord):
            keys.update({k:i})
        ord = map(keys.get,xrange(len(ord)))

        try:
            adj = adj[np.ix_(ord,ord)]
        except IndexError:
            print ord
            raise IndexError('The mask has too few ROIs')
        np.save(os.path.join(dset_dir, 'tracks', '%s_alph_normed_%s.npy'%
            (trk_output, parc)),adj)

def linearly_corrected_matrix_creation(subj, subjdir, nmr_dir, parclist_dir, 
        masks_dir, shell=None, trk_stem='tracts20k'):
    import subprocess
    import numpy as np
    import nibabel as nib
    import os.path as op
    print 'moving on to harmonic mean algorithm'

    parc='laus250'
    parc_stem='wmreg'

    tracts_stem = trk_stem if shell is None else '%s_sgsh_%s'%(trk_stem, shell)

    #extract surf array 
    if not os.path.exists(os.path.join( subjdir, subj, 'label',
            'surf_%s.table.npy'%parc)):
        os.environ['SUBJECTS_DIR']=nmr_dir

        for hemi in ('lh','rh'):
            labels = []

            if not op.exists(op.join( subjdir, subj, 'label', '%s.%s.annot'%(hemi,parc))):
                label2annot_cmd = ('mris_label2annot --s %s --h %s --ldir %s '
                        '--a %s --ctab %s --no-unknown'%(
                    '%s'%subj,
                    hemi,
                    os.path.join(subjdir, subj, 'label', parc),
                    parc,
                    os.path.join(parclist_dir, '%s_ctab_%s.txt'%(parc, hemi))))
                print label2annot_cmd
                w = subprocess.call(label2annot_cmd, shell=True)

                #mris_anatomical_stats is needed to get the surface area of each
                #label which is not the same as the number of vertices
            if not op.exists(op.join(subjdir, subj, 'label', '%s.%s.table' % (hemi, parc))):
                extract_sa_cmd = 'mris_anatomical_stats -a %s -f %s %s %s'%(
                    os.path.join(subjdir, subj, 'label',
                        '%s.%s.annot'%(hemi, parc)),
                    os.path.join(subjdir, subj, 'label',
                        '%s.%s.table'%(hemi, parc)),
                    '%s'%subj,
                    hemi)
                print extract_sa_cmd
                p = subprocess.call(extract_sa_cmd, shell=True)

        #build surf array
        lh_table = open(os.path.join( subjdir, subj, 'fsrecon', 'label',
            'lh.%s.table'%parc))
        rh_table = open(os.path.join( subjdir, subj, 'fsrecon', 'label',
            'rh.%s.table'%parc))

        labels_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp_all.txt' % parc)
        labels = open(labels_cmp_fname)

        a=[] 
        for ln in labels:
            a.append(ln.strip())

        surfArray = [0]*len(a)

        for i,r in enumerate(a):
            if r[:2]=='rh':
                for line in rh_table:
                    if filter(None, line.split(" "))[0].startswith(r[3:]):
                        surfArray[i] = filter(None, line.split(" "))[2]
            if r[:2]=='lh':
                for line in lh_table:
                    if filter(None, line.split(" "))[0].startswith(r[3:]):
                        print line
                        surfArray[i] = filter(None, line.split(" "))[2]
            rh_table.seek(0,0)
            lh_table.seek(0,0)
        np.save( 
            os.path.join( subjdir, subj, 'fsrecon', 'label', 
                'surf_%s.table.npy'%parc), surfArray)

        for file in (lh_table, rh_table, labels):
            file.close()

    #calculate fiber lengths and then add values to matrix in
    #single iteration over track file
    roi_mask = nib.load(os.path.join( masks_dir, 'bigmask_%s_%s.nii.gz'
        % (parc, parc_stem))).get_data()

    nr_rois = int(np.max(roi_mask))
    connection_matrix = np.zeros((nr_rois, nr_rois))

    def euclidean_track_length(xyz_array):
        '''see docstring in connectomemapper util.length'''

        xyz = np.asarray(xyz_array)
        if xyz.shape[0] < 2:
            return 0
        dists = np.sqrt((np.diff(xyz, axis=0)**2).sum(axis=1))
        return np.sum(dists)

    def fiber_endpoints(fiber, voxel_size, roi_mask, fov):
        '''
        returns a 2x1 matrix with indices of the masks containing the
        start and end points of this fiber.
    
        adapted from cmp.connectionmatrix.creatematrix
        '''
        f = fiber[0]

        start_xyz = np.array(f[0,:], dtype=int)
        end_xyz = np.array(f[-1,:], dtype=int)

        #account for the voxel size in the trk file which is
        #some random number based on the trk file format.
        #also account for the difference in FOV between the trk
        #file and the ROImask image which is independent.
        #print fov, roi_mask.shape, voxel_size

       # print voxel_size

        voxel_size = (2,2,2)
       # voxel_size = (1,1,2)

        voxel_size = voxel_size * fov / roi_mask.shape

        start_voxel = start_xyz // voxel_size
        end_voxel = end_xyz // voxel_size

        #DSI studio prints out voxel size as 1,1,1 no matter what (I think?)
        #start_voxel = start_xyz // (voxel_size*2)
        #end_voxel = end_xyz // (voxel_size*2)

       # print fov, roi_mask.shape
       # print start_xyz, end_xyz, voxel_size, start_voxel, end_voxel

       # dims_to_flip = (0,1,)
       # dims_to_flip = (1,2,)
        dims_to_flip = (1,)

        for i in dims_to_flip:
            start_voxel[i] = roi_mask.shape[i] - start_voxel[i] - 1
            end_voxel[i] = roi_mask.shape[i] - end_voxel[i] -1

#        print roi_mask[start_voxel[0], start_voxel[1], start_voxel[2]]
#        print [start_voxel[i] for i in xrange(3)]
#        print [[start_voxel[i]] for i in xrange(3)]
#        print roi_mask[[[start_voxel[i]] for i in xrange(3)]]

        start_roi = [roi_mask[[[start_voxel[i]] for i in xrange(3)]]]
        end_roi = [roi_mask[[[end_voxel[i]] for i in xrange(3)]]]

        #print start_roi, end_roi

        #unbox from list comprehension
        return int(start_roi[0])-1, int(end_roi[0])-1
    
    #TRACKVIS WAY

    low_cutoff, high_cutoff = 20, 500

    #if not os.path.exists(os.path.join(subjdir, subj, 'fiber_lengths.npy')):
#    tracks_file = os.path.join(subjdir, subj, '%s.trk'%tracts_stem)
#    streams, hdr = nib.trackvis.read(tracks_file, as_generator=True)
    #if int(hdr['n_count']) != nr_fibers:
    #    raise Exception('The number of fibers is is not 20000')
#    voxel_size = np.array(hdr['voxel_size'], dtype=int)
#    nr_fibers = hdr['dim']

    #low_cutoff, high_cutoff = -np.inf, np.inf

#    for i,fib in enumerate(streams):
#        fiber_length = euclidean_track_length(fib[0])
#
#        #apply fiber cutoffs -- length filtering as in CMP
#        #spline filtering should be done earlier if needed
#        if not low_cutoff < fiber_length < high_cutoff:
#            continue
#
#        start_roi, end_roi = fiber_endpoints(fib, voxel_size, roi_mask, 
#            hdr['dim'])
#
#        if start_roi==-1 or end_roi==-1:
#            continue
#        connection_matrix[start_roi, end_roi] += 1./fiber_length
#    #np.save(os.path.join(subjdir, subj, 'fiber_lengths.npy'), fiber_lengths)

    #TEXT FILE WAY
    text_file = os.path.join(subjdir, subj, '%s.txt'%tracts_stem)

    with open(text_file) as fd:
        for ln in fd:
            tract = np.reshape( np.array( ln.split(), dtype=float), (-1, 3) )
            fiber_length = euclidean_track_length(tract)

            if not low_cutoff < fiber_length < high_cutoff:
                continue

            #get endpoints manually
            #no dealing with stupid trackvis space
            start_voxel = np.floor(tract[0])
            end_voxel = np.floor(tract[-1])

            #flip dimensions
            dims_to_flip = (1,)

            for i in dims_to_flip:
                start_voxel[i] = roi_mask.shape[i] - start_voxel[i] - 1
                end_voxel[i] = roi_mask.shape[i] - end_voxel[i] -1

            start_roi = [roi_mask[[[start_voxel[i]] for i in xrange(3)]]]
            start_roi = int(start_roi[0])-1
            end_roi = [roi_mask[[[end_voxel[i]] for i in xrange(3)]]]
            end_roi = int(end_roi[0])-1

            if start_roi==-1 or end_roi==-1:
                continue

            connection_matrix[start_roi, end_roi] += 1./fiber_length
    
    connection_matrix = connection_matrix+connection_matrix.T
    print 'finished evaluating fibers'

    #normalize by surface area 


    roisl, roisr = [],[]
    print "didn't find these ROIs"
    with open( os.path.join( parclist_dir, '%s_entries'%parc)) as fd:
        for ln in fd:
            for hemi, rois in zip(('lh','rh'),(roisl, roisr)):
                roi = '%s-%s'%(ln.strip(), hemi)
                if not os.path.exists(os.path.join(masks_dir, 
                        '%s_%s'%(parc,parc_stem), '%s.nii.gz'%roi)):
                    print masks_dir, '%s_%s'%(parc, parc_stem), roi
                    continue
                rois.append(roi)
    rois = roisl+roisr

    roi_ord = []
    mangle_hemi = lambda s:'%s-%s'%(s[3:],s[:2])
    with open( os.path.join( parclist_dir, '%s_cmp.txt'%parc )) as fd:
        for ln in fd:
            roi_ord.append(mangle_hemi(ln.strip()))
    surfs_table = np.load(os.path.join(subjdir, subj, 'fsrecon', 'label',
        'surf_%s.table.npy'%parc))

    sa = [int(surfs_table[roi_ord.index(roi)]) for roi in rois]
    
    normalization = ( np.tile(sa, (np.size(sa), 1)) 
        + np.tile(sa, (np.size(sa), 1)).T )

    print 'finished normalizing'

    connection_matrix = connection_matrix * 2 / normalization 

    #np.save( os.path.join(subjdir, subj, 'linear_harmonic_connections.npy'), connection_matrix)
    #np.save( os.path.join(subjdir, subj, 'linear_normalization.npy'), normalization)

    np.fill_diagonal(connection_matrix, 0)

    np.save( os.path.join(subjdir, subj, 'linear_harmonic_%s_alph_%s.npy'%
        (parc, trk_stem)), connection_matrix)

def main():
    import sys
    import getopt

    if os.environ['SUBJECTS_DIR']:
        subjdir=os.environ['SUBJECTS_DIR']+'/'
    else:
        subjdir='local_mount/space/truffles/1/users/rlaplant/data/qball'

    subjdir = '/cluster/neuromind/dwakeman/tsc_pilot/subjects'
    os.environ['SUBJECTS_DIR'] = subjdir

    if len(sys.argv)==1:
        usage()
        return

    if not os.environ.has_key('DSI_STUDIO_FA_TEMPLATE'):
        sys.stderr.write("Source DSIstudio first!\n")
        sys.exit(85)
    if not os.environ.has_key('FREESURFER_HOME'):
        sys.stderr.write("Source Freesurfer first!\n")
        sys.exit(86)

    try:
        opts,args=getopt.getopt(sys.argv[1:],'s:m:',['skip=','nomorph',
            'shell=', 'trkstem=', 'nfibers=', 'no-fsdir', 'sbx=', 'dset=',
            'preproc=','stem=','nmr=','nmrdir='])
    except getopt.GetoptError as e:
        usage()
        sys.exit(85)

    skip_stages=[]
    nomorph=False
    shell=None
    nr_fibers=20000
    dsistudiostem='tracts20k'
    no_separate_freesurfer_dir=False
    preproc_dir = None
    preproc_stem = None
    nmr_dir = None
    dset_name = 'dset'

    for opt,arg in opts:
        if opt in ['--skip']:
            skip_stages.append(int(arg))
        if opt in ['-m']:
            modality=arg
        if opt in ['-s']:
            subj=arg
        if opt in ['--shell']:
            shell=arg
        if opt in ['--nomorph']:
            nomorph=True
        if opt in ['--nfibers']:
            nr_fibers=int(arg)
        #if opt in ['--suppress-dsistudio-output']:
        if opt in ['--trkstem']:
            dsistudiostem=arg
        if opt in ['--no-fsdir']:
            no_separate_freesurfer_dir=True
        if opt in ['--sbx']:
            preproc_dir = os.path.join(subjdir, subj, arg)
        if opt in ['--preproc','--stem']:
            preproc_stem = arg
        if opt in ['--nmr','--nmrdir']:
            nmr_dir = arg
        if opt in ['--dset']:
            dset_name = arg

    if nr_fibers != 20000 and dsistudiostem == 'tracts20k':
        raise ValueError("Must change DSIstudio save stem from default")

    if nmr_dir is None:
        #nmr_dir = '/local_mount/space/truffles/1/users/rlaplant/data/hcpnmr/'
        nmr_dir = '/local_mount/space/truffles/2/users/rlaplant/data/nmrdb/'
    if no_separate_freesurfer_dir: 
        nmr_dir = subjdir 
    #nmr_dir = '/local_mount/space/truffles/1/users/rlaplant/data/epilepsy/'
    dwi_dir = subjdir

    masks_dir = os.path.join(dwi_dir, subj, 'masks')

    dset_dir = os.path.join(dwi_dir, subj, 'diff', dset_name)
    if preproc_dir is not None:
        pass
    elif modality in ('qball','qbi','singleshell','sgsh','dti'):
        preproc_dir = os.path.join(dset_dir, 'preproc')
    elif modality in ('dsi'):
        preproc_dir = os.path.join(dwi_dir, subj, 'sbx')
    elif modality == 'mask_only':
        preproc_dir = None

    dsi_dir = '/local_mount/space/truffles/1/users/rlaplant/data/dsi'

    if modality in ('singleshell','sgsh') and shell is None:
        sys.stderr.write('Specify shell as 1k, 3k, 5k, 10k only')
        sys.exit(85)

    # parclist_dir = '/autofs/cluster/neuromind/rlaplant/mridat/parclists/'
    parclist_dir = '/cluster/neuromind/npeled/parclists'

    if not 0 in skip_stages:
        if modality in ('singleshell','sgsh'):
            organize_singleshell(subj, preproc_dir, shell)
    if not 1 in skip_stages:
        setup_mprage(subj, nmr_dir, dwi_dir, preproc_dir, dset_dir, modality)
    if not 2 in skip_stages:
        white_matter_mask(subj, subjdir, dset_dir)
    if not 3 in skip_stages:
        morph_and_extract_labels(subj, nmr_dir, subjdir, parclist_dir,
            masks_dir, nomorph=nomorph)
    if not 4 in skip_stages:
        bigmask_creation(subj, subjdir, masks_dir)
    if not 5 in skip_stages:
        run_dsistudio(subj, subjdir, masks_dir, preproc_dir, dset_dir,
            modality, shell=shell, nr_fibers=nr_fibers, 
            trk_output=dsistudiostem, preproc_stem=preproc_stem)
    if not 7 in skip_stages:
        linearly_corrected_matrix_creation(subj, subjdir, nmr_dir, parclist_dir,
         masks_dir, trk_stem=dsistudiostem)
    
if __name__=='__main__':
    main()
