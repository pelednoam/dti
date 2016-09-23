#!/usr/bin/env python
from __future__ import division
import sys
import os
import os.path as op
import glob
import traceback
import errno
import numpy as np
from scipy import io as sio
import nibabel as nib
import subprocess
import shutil
import logging
import mne

logging.basicConfig(filename='log.log',level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


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

    if op.exists(op.join(preproc_dir,'%s.nii'%stem)):
        gzip_cmd = ('mri_convert %s %s'%
            (op.join(preproc_dir,'%s.nii'%stem),
             op.join(preproc_dir,'%s.nii.gz'%stem)))
        print gzip_cmd
        m=subprocess.call(gzip_cmd,shell=True)

    a = nib.load(op.join(preproc_dir,'%s.nii.gz'%stem))
    img = nib.Nifti1Image(a.get_data()[:,:,:,ixes], a.get_affine(),
        a.get_header())
    nib.save(img, op.join(preproc_dir,'%s.nii.gz'%new_stem))

    b = np.loadtxt(op.join(preproc_dir,btable_file))
    np.savetxt(op.join(preproc_dir,new_btable_file), b[ixes])

    print 'done with file conversions for singleshell analysis'

######################
def setup_mprage(subj,nmr_dir,dwi_dir,preproc_dir,dset_dir,modality):
    #todo: send a parameter where the fsrecon folder
    fsrecon_dir = op.join(dwi_dir,subj) #,'fsrecon')
    make_dir(dset_dir)
    make_dir(preproc_dir)
    mprage_img = op.join(fsrecon_dir,'mri','orig.mgz')
    if not op.isfile(mprage_img):
        raise Exception('No mprage_img! ({})'.format(mprage_img))
    mprage_256 = op.join(dset_dir,'mprage_256.nii.gz')
    mprage_nii = op.join(dset_dir,'mprage.nii.gz')
    if not op.isfile(mprage_nii):
        raise Exception('No mprage_nii! ({})'.format(mprage_nii))

    b0_img = op.join(dset_dir, 'b0.nii.gz')
    nodif_brain = op.join(dset_dir, 'nodif_brain.nii.gz')

    mri_convert(mprage_img, mprage_256)
    # convert_cmd = 'mri_convert %s %s' % (mprage_img, mprage_256)
    # if not op.isfile(mprage_256):
    #     print convert_cmd
    #     p=subprocess.call(convert_cmd,shell=True)

    if modality in ('dsi'):
        mri_convert(op.join(preproc_dir,'b0.nii'), b0_img)
        # print preproc_dir
        # b0conv_cmd = 'mri_convert %s %s' % (
        #     op.join(preproc_dir,'b0.nii'),
        #     b0_img)
        # q=subprocess.call(b0conv_cmd,shell=True)

    elif modality in ('qball','qbi','singleshell','sgsh'):
        create_link(op.join(dset_dir, b0_img), op.join(preproc_dir, 'b0_first.nii.gz'))
        # try:
        #     os.link(op.join(dset_dir, b0_img), op.join(preproc_dir, 'b0_first.nii.gz'))
        #     print 'ln preproc/diff/b0_first.nii.gz b0.nii.gz'
        # except OSError as err:
        #     if err.errno != errno.EEXIST:
        #         print(traceback.format_exc())

    elif modality == 'mask_only':
        pass

    #reorient_cmd = 'fslreorient2std %s %s' % (mprage_nii, mprage_nii)
    #print reorient_cmd
    #q=subprocess.call(reorient_cmd,shell=True)

    mprage_brain = op.join(dset_dir, 'mpr_brain.nii.gz')
    if not op.isfile(mprage_brain):
        # BET (Brain Extraction Tool) deletes non-brain tissue from an image of the whole head
        bet_cmd = 'bet %s %s -f 0.40 -B' % (mprage_nii, mprage_brain)
        run_script_and_check_output(bet_cmd, 'bet')

    if not op.isfile(nodif_brain):
        bet2_cmd = 'bet %s %s -f 0.15 -g 0 -m' % (b0_img, nodif_brain)
        run_script_and_check_output(bet2_cmd, 'bet')

    upsample_cmd = ('mri_convert -vs 1 1 1 %s %s' %
        ( b0_img,
          op.join(dset_dir, 'b0d.nii.gz')))

    mpr2562b0_fname = op.join(dset_dir, 'mpr2562b0.mat')
    if not op.isfile(mpr2562b0_fname):
        # FLIRT (FMRIB's Linear Image Registration Tool) is a fully automated robust and accurate tool for
        # linear (affine) intra- and inter-modal brain image registration.
        matflirt_cmd = ('flirt -in %s -ref %s -omat %s' %
            (mprage_256, b0_img, mpr2562b0_fname))
        run_script_and_check_output(matflirt_cmd, 'flirt', output_fname=mpr2562b0_fname)

    mprage_b0_fname = op.join(dset_dir, 'mprage_b0.nii.gz')
    if not op.isfile(mprage_b0_fname):
        flirt_cmd = ('flirt -in %s -ref %s -applyxfm -init %s -out %s' %
            (mprage_256, b0_img, mpr2562b0_fname, mprage_b0_fname))
        run_script_and_check_output(flirt_cmd, 'flirt', output_fname=mprage_b0_fname)

    print 'Done with mprage reorientation and masking'


def get_hemi_indifferent_roi(roi):
    return roi.replace('-rh', '').replace('-lh', '').replace('rh-', '').replace('lh-', '').\
        replace('.rh', '').replace('.lh', '').replace('rh.', '').replace('lh.', '').\
        replace('Right-', '').replace('Left-', '').replace('-Right', '').replace('-Left', '').\
        replace('Right.', '').replace('Left.', '').replace('.Right', '').replace('.Left', '').\
        replace('right-', '').replace('left-', '').replace('-right', '').replace('-left', '').\
        replace('right.', '').replace('left.', '').replace('.right', '').replace('.left', '')


def get_roi_hemi(roi):
    l_roi = roi.lower()
    for hemi, big_hemi in zip(['rh', 'lh'], ['right', 'left']):
        for h in [hemi, big_hemi]:
            for d in ['-', '_', '.']:
                if '%s%s' % (d, h) in l_roi:
                    return roi[:-3], hemi
                if '%s%s' % (h, d) in l_roi:
                    return roi[3:], hemi
    if roi != 'delete':
        print("Can't find hemi in %s!" % roi)
    return None, None

#####################
def morph_and_extract_labels(subj,nmr_dir,subjdir,parclist_dir,masks_dir,dset_dir,
    nomorph=False, n_jobs=1):

    os.chdir(subjdir)

    parc_stem = 'wmreg'

    # try:
    #     os.symlink( op.join(nmr_dir, subj), #'%s_FS'%subj),
    #                 op.join(subj,'fsrecon'))
    # except OSError:
    #     print(traceback.format_exc())

    morph_cmd = ('mri_label2label --srcsubject fsaverage5c '
        '--srclabel fsaverage5c/label/%s/%s.%s --trgsubject %s '
        '--trglabel %s/label/%s/%s.%s --regmethod surface --hemi %s')

    #extract_cmd = ('mri_label2vol --label %s/fsrecon/label/%s/%s.%s.label '
    #	'--temp %s/dset/mprage.nii.gz --surf white '
    #	'--o %s/masks/%s/%s-%s.nii.gz --regheader %s/fsrecon/mri/orig.mgz')

    #extract_cmd = ('mri_label2vol --label %s/fsrecon/label/%s/%s.%s.label '
    #	'--temp %s/dset/mprage_b0.nii.gz --surf white '
    #	'--o %s/masks/%s_b0/%s-%s.nii.gz --regheader %s/fsrecon/mri/orig.mgz')

    extract_cmd = ('mri_label2vol --label %s/label/%s/%s.%s.label '
        '--temp %s/mprage_256.nii.gz --surf white '
        '--o %s/masks/%s_%s/%s-%s.nii.gz --regheader %s/mri/wm.mgz')

    flirt_cmd = ('flirt -in %s/masks/%s_%s/%s-%s.nii.gz -ref %s/b0.nii.gz '
        '-applyxfm -init %s/mpr2562b0.mat '
        '-out %s/masks/%s_%s_bspace/%s-%s.nii.gz')

    for parc in ('laus250',):
    #for parc in ('laus125','laus250','laus500'):
        make_dir(op.join(masks_dir,'%s_%s'%(parc,parc_stem)))
        make_dir(op.join(masks_dir, '%s_%s_bspace' % (parc,parc_stem)))
        make_dir(op.join(subj, 'label', parc))

        entries=[]
        # entry_file = op.join(parclist_dir,'%s_entries'%parc)
        entry_file = op.join(parclist_dir, 'orders', '%s_cmp.txt'%parc)
        with open(entry_file,'r') as fd:
            for ln in fd:
                if 'delete' not in ln:
                    entries.append(ln.strip())

        # for label_fname in glob.glob('fsaverage5c/label/%s/*.label' % parc):
        #     entry = op.basename(op.split(label_fname)[1])
        #     entry = get_hemi_indifferent_roi(entry).replace('.label', '')
        labels_chunks = chunks(entries, len(entries) / n_jobs)
        params = [(subj, parc, parc_stem, labels_chunk, subjdir, dset_dir, masks_dir, morph_cmd, extract_cmd, flirt_cmd, nomorph)
                  for labels_chunk in labels_chunks]
        results = run_parallel(_morph_and_extract_label, params, n_jobs)
        labels_num = 0
        for result in results:
            for err in result:
                if err != '':
                    print(err)
                else:
                    labels_num += 1
    print('Done with morphing and extracting. {} labels were created'.format(labels_num))


def chunks(l, n):
    n = int(max(1, n))
    return [l[i:i + n] for i in range(0, len(l), n)]


def run_parallel(func, params, njobs=1):
    import multiprocessing
    if njobs == 1:
        results = [func(p) for p in params]
    else:
        pool = multiprocessing.Pool(processes=njobs)
        results = pool.map(func, params)
        pool.close()
    return results


def _morph_and_extract_label(p):
    subj, parc, parc_stem, labels_chunk, subjdir, dset_dir, masks_dir, morph_cmd, extract_cmd, flirt_cmd, nomorph = p
    errors = []
    for entry in labels_chunk:
        # print('Morphing and extracting label {}'.format(entry))
        entry, hemi = get_roi_hemi(entry)
        # hemi = entry[:2]
        # entry = entry[3:]

        if op.isdir(op.join(masks_dir, '{}_{}_all'.format(parc, parc_stem))):
            shutil.copy(op.join(masks_dir, '{}_{}_all'.format(parc, parc_stem), '{}-{}.nii.gz'.format(entry, hemi)),
                        op.join(masks_dir, '{}_{}'.format(parc, parc_stem), '{}-{}.nii.gz'.format(entry, hemi)))
            errors.append('')
            continue

        source_label_fname = op.join(subjdir, 'fsaverage5c', 'label', parc, '{}.{}.label'.format(hemi, entry))
        if not op.isfile(source_label_fname):
            logging.warning('{}.{} not in fsaverage5c/label/{}'.format(hemi, entry, parc))
            errors.append('')
            continue

        # 'Now morphing and extracting %s'%subj
        m_cmd = morph_cmd % (parc, hemi, entry, subj, subj, parc, hemi, entry, hemi)
        target_label_fname = op.join(subjdir, subj, 'label', parc, '{}.{}.label'.format(hemi, entry))
        if not nomorph:
            err = run_script_and_check_output(m_cmd, 'mri_label2label', raise_error=False, output_fname=target_label_fname)
            if err != '':
                errors.append(err)
                continue

        e_cmd = extract_cmd % (subj, parc, hemi, entry, dset_dir, subj, parc,
                               parc_stem, entry, hemi, subj)
        target_label_fname = op.join(masks_dir, '{}_{}'.format(parc, parc_stem), '{}-{}.nii.gz'.format(entry, hemi))
        err = run_script_and_check_output(e_cmd, 'mri_label2vol', raise_error=False, output_fname=target_label_fname)
        if err != '':
            errors.append(err)
            continue

        f_cmd = flirt_cmd % (subj, parc, parc_stem, entry, hemi, dset_dir, dset_dir,
                             subj, parc, parc_stem, entry, hemi)
        output_fname = op.join(masks_dir, '{}_{}_bspace'.format(parc, parc_stem), '{}-{}.nii.gz'.format(entry, hemi))
        err = run_script_and_check_output(f_cmd, 'flirt', raise_error=False, output_fname=output_fname)
        if err != '':
            errors.append(err)
            continue
        errors.append('')
    return errors


def run_script_and_check_output(cmd, script_name, done_str='done', do_print_cmd=True, raise_error=True, output_fname=''):
    if output_fname == '':
        if script_name == 'bet':
            output_fname = cmd.split(' ')[2]
        elif script_name == 'mri_convert':
            output_fname = cmd.split(' ')[2]
        elif '>' in cmd and len(cmd.split('>')) == 2:
            output_fname = cmd.split('>')[1].strip()
        else:
            output_name_arg = get_opt_val(cmd, '--output')
            if output_name_arg != '':
                output_fname = output_name_arg
    if output_fname != '' and op.isfile(output_fname):
        return ''
    q = subprocess.check_output('{} | tee /dev/stderr'.format(cmd), shell=True)
    if do_print_cmd:
        print cmd
    err = ''
    if output_fname != '':
        if not op.isfile(output_fname):
            err = 'output file {} was not created!'.format(output_fname)
    elif done_str != '' and not done_str in q and not done_str.title() in q:
        err = 'Error running the {} script! {}'.format(script_name, cmd)
        print(q)
    if err != '':
        logging.error(err)
        if raise_error:
            raise Exception(err)
        else:
            print(err)
    return err


def get_dirents(centroids_dir, only_namebases=False):
    dirents = []
    for hemi_path in ('lh.nii.gz','rh.nii.gz'):
        for dirent in sorted(os.listdir(centroids_dir)):
            if dirent.endswith(hemi_path):
                if only_namebases:
                    dirents.append(dirent[:-(len('nii.gz') + 1)])
                else:
                    dirents.append(dirent)
    return dirents

#####################
def bigmask_creation(subj,subjects_dir,masks_dir):
    #mangle_hemi = lambda x:'%s_%s'%(x[-9:-7],x[:-10])
    #demangle_hemi = lambda x:'%s-%s.nii.gz'%(x[3:],x[:2])

    #parc_orig = 'laus125'
    #parc = '%s_nonoverlapping' % parc_orig
    parcs = ['laus250_wmreg']
    for parc in parcs:
        output_fname = op.join(masks_dir,'bigmask_%s.nii.gz'%parc)
        if op.isfile(output_fname):
            continue
        #masks_dir=op.join(masks_top_lev_dir,parc)
        centroids_dir = op.join(masks_dir, '{}_bspace'.format(parc))

        #copy parc_orig to parc if parc doesnt exist yet
        if not op.exists(masks_dir):
            raise ValueError('No such parcellation')
            #os.makedirs(masks_dir)
            #shutil.copytree(centroids_dir,masks_dir,symlinks=True)

        mask_template = nib.load(op.join(centroids_dir, 'bankssts_1-lh.nii.gz'))

        mask_shape = mask_template.shape
        big_mask = np.zeros((mask_shape[0],mask_shape[1],mask_shape[2]),dtype=np.float32)
        binary_mask = big_mask.copy()

        print('building mask in diffusion space with shape {0}'.format(mask_shape))

        #create centroids
        print("collecting centroids")
        centroids = {}

        #make sure dirents are in alphabetical order
        dirents = get_dirents(centroids_dir)

        for i, dirent in enumerate(dirents):
            logging.debug('collecting centroid for subsequent mask {} {}'.format(i, dirent))
            try:
                r=nib.load(op.join(centroids_dir,dirent)).get_data()
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
            closest_centroid = 0
            dist = np.inf
            nr_collisions = 2 	#divisor starts at 2
                            #ensuring that a random voxel will be chosen
            for i,centroid in centroids.items():
                #print vox,centroid
                cent_dist=np.linalg.norm(vox-centroid)
                if cent_dist < dist:
                    dist = cent_dist
                    closest_centroid=i+1

                if cent_dist == dist:
                    if np.random.random() < 1./nr_collisions:
                        dist=cent_dist
                        closest_centroid=i+1
                    nr_collisions+=1
                    
            big_mask[x,y,z]=closest_centroid

        img = nib.Nifti1Image(big_mask,mask_template.get_affine(),
            mask_template.get_header())
        nib.save(img, output_fname)

        print '%s bigmask saved'%parc

    success = np.all(np.array([op.isfile(op.join(masks_dir, 'bigmask_%s.nii.gz' % parc)) for parc in parcs]))
    if not success:
        logging.error('Not all the big masks were created!')
        print('Not all the big masks were created!')
    else:
        print('All the big masks were created!')
        logging.info('All the big masks were created!')


def white_matter_mask(subj,subjdir,dset_dir):
    convert_cmd = ('mri_convert %s %s' % (
        op.join( subjdir, subj, 'mri', 'wm.mgz'),
        op.join( dset_dir, 'wm_256.nii.gz')))
    if not op.isfile(op.join(dset_dir, 'wm_256.nii.gz')):
        run_script_and_check_output(convert_cmd, 'mri_convert')

    flirt_cmd = ('flirt -in %s -ref %s -applyxfm -init %s -out %s' % (
        op.join(dset_dir, 'wm_256.nii.gz'),
        op.join(dset_dir, 'b0.nii.gz'),
        op.join(dset_dir, 'mpr2562b0.mat'),
        op.join(dset_dir, 'wm_b0.nii.gz')))
    output_fname = op.join(dset_dir, 'wm_b0.nii.gz')
    if not op.isfile(output_fname):
        run_script_and_check_output(flirt_cmd, 'flirt', output_fname=output_fname)


def run_dsistudio(subj, subjdir, masks_dir, preproc_dir, dset_dir, modality,
        shell='10k', nr_fibers=20000, trk_output='tracts20k', 
        suppress_output=False, preproc_stem=None):

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
        logging.error('Unrecognized modality %s!'%modality)
        raise Exception('Unrecognized modality %s!'%modality)

    if preproc_stem is not None:
        stem = preproc_stem

    make_dir(preproc_dir)

    if op.isfile(op.join(dset_dir,'%s.nii.gz'%stem)) and not op.islink(op.join(preproc_dir, '%s.nii.gz'%stem)):
        create_link(op.join(dset_dir, '%s.nii.gz'%stem), op.join(preproc_dir,'%s.nii.gz'%stem))
    elif op.exists(op.join(preproc_dir,'%s.nii'%stem)) and not op.isfile(op.join(preproc_dir,'%s.nii.gz'%stem)):
        gzip_cmd = ('mri_convert %s %s'%
            (op.join(preproc_dir,'%s.nii'%stem),
             op.join(preproc_dir,'%s.nii.gz'%stem)))
        run_script_and_check_output(gzip_cmd, 'mri_convert')

    if not op.isfile(op.join(preproc_dir,btable_string)):
        bvals = np.genfromtxt(op.join(dset_dir, 'bvals'))
        bvals = bvals.reshape((len(bvals), 1))
        bvecs = np.genfromtxt(op.join(dset_dir, 'bvecs'))
        btable = np.hstack((bvals, bvecs))
        np.savetxt(op.join(preproc_dir,btable_string), btable)

    if modality in ('qball','qbi','singleshell','sgsh'):
        src_cmd=('dsi_studio --action=src --source=%s --b_table=%s --output=%s'
            % (
            op.join(preproc_dir,'%s.nii.gz' % stem),
            op.join(preproc_dir,btable_string),
            op.join(preproc_dir,'%s.src.gz' % stem)))
    elif modality in ('dsi',):
        src_cmd=('dsi_studio --action=src --source=%s --b_table=%s --output=%s'
            % (
            op.join(preproc_dir,'%s.nii.gz' % stem),
            op.join(subjdir, btable_string),
            op.join(preproc_dir,'%s.src.gz' % stem)))
    elif modality in ('dti',):
        src_cmd=('dsi_studio --action=src --source=%s --output=%s --b_table=%s' % (
            op.join(preproc_dir, 'dti.nii.gz'),
            op.join(preproc_dir, '%s.src.gz' % stem),
            op.join(preproc_dir,btable_string)))
    if not op.isfile(op.join(preproc_dir, '%s.src.gz' % stem)):
        run_script_and_check_output(src_cmd, 'dsi_studio')

    # if not op.isfile(op.join(preproc_dir, '%s.src.gz' % stem)):
    #     raise Exception("Can't find the dsi_studio output!")

    rec_cmd=('dsi_studio --action=rec --source=%s --thread=10 --method=%s '
        '--num_fiber=8 --odf_order=8 --record_odf=1 --mask=%s' % (
        op.join(preproc_dir,'%s.src.gz' % stem),
        methodnr,
        op.join(dset_dir,'wm_b0.nii.gz')))
    if param0 is not None:
        rec_cmd += ' --param0=%s'%param0
    if param1 is not None:
        rec_cmd += ' --param1=%s'%param1
    source_output_fname = '%s.src.gz.%s%s.%sfib.gz' % (
        stem, '' if methodname=='dti' else 'odf8.f8rec.', methodname,
        reduce(lambda x,y:x+y,['sh%s.'%p if i==0 else'%s.'%p for i,p in enumerate(
            (param0,param1)[::-1])if p is not None],''))
    run_script_and_check_output(rec_cmd, 'dsi_studio', output_fname=op.join(preproc_dir, source_output_fname))

    trk_cmd=('dsi_studio --action=trk --source=%s --fiber_count=%i --output=%s' % (
        op.join(preproc_dir, source_output_fname), nr_fibers,
        op.join(dset_dir, 'tracks','%s.%s'%(trk_output, '%s'))))
    # print trk_cmd
    if not op.isfile(op.join(dset_dir, 'tracks', 'tracts{}k.txt'.format(int(nr_fibers/1000)))):
        run_script_and_check_output(trk_cmd%'txt', 'dsi_studio')
    if not op.isfile(op.join(dset_dir, 'tracks', 'tracts{}k.trk'.format(int(nr_fibers / 1000)))):
        run_script_and_check_output(trk_cmd%'trk', 'dsi_studio')

    parcs = ['laus250']
    for parc in parcs:
        ana_cmd=('dsi_studio --action=ana --source=%s --tract=%s '
            '--export=connectivity --roi=%s' % (
            op.join(preproc_dir, source_output_fname),
            #op.join(subjdir,subj,'%s.trk'%trk_output),
            op.join(dset_dir, 'tracks','%s.txt'%trk_output),
            op.join(masks_dir,'bigmask_%s_wmreg.nii.gz'%parc)))
        connectivity_output_fname = op.join(dset_dir, 'tracks', '%s.txt.connectivity.mat'%trk_output)
        run_script_and_check_output(ana_cmd, 'dsi_studio', output_fname=connectivity_output_fname)


        #==================================================================
        #put the resulting matrix back in alphabetical order

        alph_normed_fname = op.join(dset_dir, 'tracks', '%s_alph_normed_%s.npy'%
            (trk_output, parc))

        if not op.isfile(alph_normed_fname):
            fname=op.join(dset_dir, 'tracks', '%s.txt.connectivity.mat'%trk_output)
            adj = sio.loadmat(fname)
            labnam = adj['name']
            med_length = adj['tract_median_length']

            adj = adj['connectivity']

            adj = adj/med_length
            print('This RuntimeWarning is ok')
            adj[np.isnan(adj)]=0

            labnam = ''.join(map(chr,labnam.flat)).split('\n')[:-1]

            stem = len(op.commonprefix(labnam))
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
            np.save(alph_normed_fname, adj)

    success = np.all(np.array([op.isfile(op.join(dset_dir, 'tracks', '{}_alph_normed_{}.npy'.format(
        trk_output, parc))) for parc in parcs]))
    if not success:
        logging.error('Not all the connectivity matrices were created!')
        print('Not all the connectivity matrices were created!')
    else:
        print('All the connectivity matrices were created!')
        logging.info('All the connectivity matrices were created!')
        print('Finish with run_dsistudio!')


def create_ctab_files(subj, subjdir, parc, parclist_dir):
    from mne.label import _read_annot
    import csv
    from collections import defaultdict
    annot_fname = op.join(subjdir, subj, 'label', '{}.{}.annot'.format('{hemi}', parc))
    if not op.isfile(annot_fname.format(hemi='rh')) or not op.isfile(annot_fname.format(hemi='lh')):
        labels = defaultdict(list)
        labels_files = glob.glob(os.path.join(subjdir, subj, 'label', parc, '*.label'))
        if len(labels_files) == 0:
            raise Exception('create_ctab_files: No labels files!')
        for label_file in labels_files:
            label = mne.read_label(label_file)
            if 'unknown' in label.name:
                continue
            # print(label.name)
            label.name = label.name if isinstance(label.name, str) else label.name.astype(str)
            labels[label.hemi].append(label)
        for hemi in ['rh', 'lh']:
            labels[hemi].sort(key=lambda l: l.name)
            mne.write_labels_to_annot(subject=subj, labels=labels[hemi], parc=parc, subjects_dir=subjdir,
                                      hemi=hemi, overwrite=True)
    if not op.isfile(annot_fname.format(hemi='rh')) or not op.isfile(annot_fname.format(hemi='lh')):
        logging.error('Not all the annotation files were created!')
        raise Exception('Not all the annotation files were created!')

    for hemi in ['rh', 'lh']:
        new_lut_fname = op.join(subjdir, subj, 'label', '{}_ctab_{}.txt'.format(parc, hemi))
        if op.isfile(new_lut_fname):
            continue
        lut_new = []
        _, ctab, names = _read_annot(annot_fname.format(hemi=hemi))
        names = [name.astype(str) for name in names]
        for index, (label, cval) in enumerate(zip(names, ctab)):
            r, g, b, a, _ = cval
            lut_new.append([index, label, r, g, b, a])
        with open(new_lut_fname, 'w') as fp:
            csv_writer = csv.writer(fp, delimiter='\t')
            csv_writer.writerows(lut_new)
        if not op.isfile(new_lut_fname):
            err_msg = 'The ctab file for hemi {} was not created! ({})'.format(hemi, new_lut_fname)
            logging.error(err_msg)
            raise Exception(err_msg)


def _calc_label_stat(p):
    subjdir, subj, parc, labels_fnames = p
    res = []
    make_dir(op.join(subjdir, subj, 'label', '%s_stat' % parc))
    for label_fname in labels_fnames:
        label_name = op.basename(label_fname)
        _, hemi = get_roi_hemi(label_name)
        output_fname = op.join(subjdir, subj, 'label', '%s_stat' % parc,
                               '%s.%s.%s.out.table' % (hemi, parc, label_name))
        table_fname = op.join(subjdir, subj, 'label', '%s_stat' % parc, '%s.%s.%s.table' % (hemi, parc, label_name))
        annot_fname = op.join(subjdir, subj, 'label', '%s.%s.annot' % (hemi, parc))
        extract_sa_cmd = 'mris_anatomical_stats -a %s -f %s -l %s %s %s > %s' % (
            annot_fname, table_fname, label_fname, subj, hemi, output_fname)
        run_script_and_check_output(extract_sa_cmd, 'mris_anatomical_stats', output_fname=table_fname)
        with open(output_fname, 'r') as stat_file:
            line = next(stat_file)
            while 'total surface area' not in line:
                line = next(stat_file)
            total_surface_area = int(line.split(' ')[-2])
        res.append((label_name, total_surface_area))
    return res


def linearly_corrected_matrix_creation(subj, subjdir, dset_dir, parclist_dir,
        masks_dir, shell=None, trk_stem='tracts20k', parc='laus250', parc_stem='wmreg', n_jobs=1):

    def create_orders_files():
        labels_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp.txt' % parc)
        labels_alph_fname = op.join(parclist_dir, 'orders', '%s_alph.txt' % parc)
        if not op.isfile(labels_cmp_fname) or not op.isfile(labels_alph_fname):
            print('create order files!')
            import getlausords
            getlausords.create_cmp_alph_files(subjdir, parclist_dir, parc)

    def get_order_labels_names():
        labels_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp_all.txt' % parc)
        create_orders_files()
        with open(labels_cmp_fname) as labels:
            labels_names = []
            for ln in labels:
                labels_names.append(ln.strip())
        return labels_names

    def call_calc_label_stat_parallel():
        surf_table_fname = op.join(subjdir, subj, 'label', 'surf_%s.table.npy' % parc)
        calc_labels_total_surface_area(get_order_labels_names())
        if not op.isfile(surf_table_fname):
            ordered_labels_names = get_order_labels_names()
            ordered_labels_names = [get_roi_hemi(l) for l in ordered_labels_names]
            labels_names = glob.glob(op.join(subjdir, subj, 'label', parc, '*.label'))
            labels_chunks = chunks(labels_names, len(labels_names) / n_jobs)
            params = [(subjdir, subj, parc, labels_fnames_chunk) for labels_fnames_chunk in labels_chunks]
            results = run_parallel(_calc_label_stat, params, n_jobs)
            surfArray = [0] * len(ordered_labels_names)
            for res in results:
                for label_name, total_surface_area in res:
                    name, hemi = get_roi_hemi(op.splitext(label_name)[0])
                    index = ordered_labels_names.index((name, hemi))
                    if index >= 0:
                        surfArray[index] = total_surface_area
                    else:
                        print("Can't find {} {}".format(name, hemi))
            np.save(surf_table_fname, surfArray)

    def calc_labels_total_surface_area():
        output_fname = op.join(subjdir, subj, 'label', 'surf_%s.table.npy' % parc)
        lh_table = open(op.join(subjdir, subj, 'label', 'lh.%s.table' % parc))
        rh_table = open(op.join(subjdir, subj, 'label', 'rh.%s.table' % parc))
        labels_names = get_order_labels_names()
        print('len of surfArray = {}'.format(len(labels_names)))
        surfArray = [0] * len(labels_names)

        for i, r in enumerate(labels_names):
            for table, hemi in zip([rh_table, lh_table], ['rh', 'lh']):
                label_no_hemi, label_hemi = get_roi_hemi(r)
                if label_no_hemi is None:
                    continue
                if label_hemi == hemi:
                    for linue_num, line in enumerate(table):
                        if label_no_hemi in line:
                            # if filter(None, line.split(" "))[0][3:].startswith(r[3:]):
                            next_line = next(table)
                            while 'number of vertices' not in next_line:
                                if 'total surface area' in next_line:
                                    raise Exception("Can't find number of vertices for {} in {}".format(
                                        label_no_hemi, table.name))
                                next_line = next(table)
                            # surfArray[i] = int(filter(None, line.split(" "))[2])
                            surfArray[i] = int(next_line.split('=')[-1].strip())
                            # print '%s: %s' % (r, surfArray[i])
                            break
            rh_table.seek(0, 0)
            lh_table.seek(0, 0)

        surfArray = np.array(surfArray)
        check_surface_statistics(surfArray, labels_names)
        np.save(output_fname, surfArray)
        for file in (lh_table, rh_table):
            file.close()
        return surfArray

    def check_surface_statistics(surfArray, labels_names):
        labels = {hemi: [get_hemi_indifferent_roi(l.name) for l in mne.read_labels_from_annot(
            subjects_dir=subjdir, subject=subj, parc=parc, hemi=hemi)] for hemi in ['rh', 'lh']}
        no_labels_data = np.array(labels_names)[np.where(surfArray == 0)[0]]
        for no_label in no_labels_data:
            if no_label == 'delete':
                continue
            no_label_name, hemi = get_roi_hemi(no_label)
            if no_label_name in labels[hemi]:
                logging.error('No statistics found for {}!'.format(no_label))
                raise Exception('No statistics found for {}!'.format(no_label))



    print 'moving on to harmonic mean algorithm'
    linear_harmonic_output_fname = op.join(dset_dir, 'tracks', 'linear_harmonic_%s_alph_%s.npy' % (parc, trk_stem))
    if op.isfile(linear_harmonic_output_fname):
        return

    tracts_stem = trk_stem if shell is None else '%s_sgsh_%s'%(trk_stem, shell)

    #extract surf array 
    if not op.exists(op.join(subjdir, subj, 'label', 'surf_%s.table.npy'%parc)):
        # os.environ['SUBJECTS_DIR']=nmr_dir
        create_ctab_files(subj, subjdir, parc, parclist_dir)
        for hemi in ('lh','rh'):

            if not op.exists(op.join( subjdir, subj, 'label', '%s.%s.annot'%(hemi,parc))):
                label2annot_cmd = ('mris_label2annot --s %s --h %s --ldir %s '
                        '--a %s --ctab %s --no-unknown'%(subj,hemi,
                    op.join(subjdir, subj, 'label', parc), parc,
                    op.join(parclist_dir, '%s_ctab_%s.txt'%(parc, hemi))))
                run_script_and_check_output(label2annot_cmd, 'mris_label2annot')

                #mris_anatomical_stats is needed to get the surface area of each
                #label which is not the same as the number of vertices
            if not op.exists(op.join( subjdir, subj, 'label', '%s.%s.table'%(hemi,parc))):
                extract_sa_cmd = 'mris_anatomical_stats -a %s -f %s %s %s > %s'%(
                    op.join(subjdir, subj, 'label', '%s.%s.annot'%(hemi, parc)),
                    op.join(subjdir, subj, 'label', '%s.%s.table'%(hemi, parc)),
                    '%s'%subj, hemi,
                    op.join(subjdir, subj, 'label', '%s.%s.table' % (hemi, parc)))
                run_script_and_check_output(extract_sa_cmd, 'mris_anatomical_stats')

        # if not op.isfile(surf_table_fname):
    surfArray = calc_labels_total_surface_area()

    #calculate fiber lengths and then add values to matrix in
    #single iteration over track file
    roi_mask = nib.load(op.join( masks_dir, 'bigmask_%s_%s.nii.gz'
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

    #if not op.exists(op.join(subjdir, subj, 'fiber_lengths.npy')):
#    tracks_file = op.join(subjdir, subj, '%s.trk'%tracts_stem)
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
#    #np.save(op.join(subjdir, subj, 'fiber_lengths.npy'), fiber_lengths)

    #TEXT FILE WAY
    text_file = op.join(dset_dir, 'tracks', '%s.txt'%tracts_stem)

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

            if start_roi < len(connection_matrix) and end_roi < len(connection_matrix):
                connection_matrix[start_roi, end_roi] += 1./fiber_length
            # else:
            #     print('problem with {}'.format(ln))
    
    connection_matrix = connection_matrix+connection_matrix.T
    print 'finished evaluating fibers'

    #normalize by surface area

    roisl, roisr = [],[]
    print "didn't find these ROIs"
    labels = {hemi: [get_hemi_indifferent_roi(l.name) for l in mne.read_labels_from_annot(
        subjects_dir=subjdir, subject=subj, parc=parc, hemi=hemi)] for hemi in ['rh', 'lh']}
    with open( op.join( parclist_dir, '%s_entries'%parc)) as fd:
        for ln in fd:
            for hemi, rois in zip(('lh','rh'),(roisl, roisr)):
                roi = '%s-%s'%(ln.strip(), hemi)
                if not op.exists(op.join(masks_dir, '%s_%s'%(parc,parc_stem), '%s.nii.gz'%roi)):
                    if roi in labels[hemi]:
                        print masks_dir, '%s_%s'%(parc, parc_stem), roi
                    continue
                rois.append(roi)
    rois = roisl+roisr

    # create_orders_files()
    roi_ord = []
    mangle_hemi = lambda s:'%s-%s'%(s[3:],s[:2])
    # centroids_dir=op.join(masks_dir, '{}_{}'.format(parc, parc_stem))
    # dirents = get_dirents(centroids_dir, only_namebases=True)
    # import glob
    # masks = glob.glob(masks_dir, '{}_{}'.format(parc, parc_stem), '*.nii.gz')
    # for mask_fname in masks:
    #     roi = namebase(mask_fname)
    #     roi_ord.append(roi)
    with open( op.join( parclist_dir, 'orders', '%s_cmp_all.txt'%parc )) as fd:
        for ln in fd:
            roi, hemi = get_roi_hemi(ln.strip())
            roi_ord.append('%s-%s' % (roi, hemi))

    surfs_table = np.load(op.join(subjdir, subj, 'label', 'surf_%s.table.npy'%parc))
    print('surfs_table len: {}'.format(len(surfs_table)))
    sa = []
    zeros_norm_rois = []
    zeros_norm_rois_indices = []
    for roi in rois:
        if roi not in roi_ord:
            print('{} not in list!'.format(roi))
        else:
            # if dirents.index(roi) >= len(surfs_table):
            #     print('index of {} is bigger then the surfs table len ({})!'.format(roi, len(surfs_table)))
            # else:
            norm_val = int(surfs_table[roi_ord.index(roi)])
            if norm_val == 0:
                print('norm val zero! {}'.format(roi))
                zeros_norm_rois.append(roi)
                zeros_norm_rois_indices.append(roi_ord.index(roi))
            sa.append(norm_val)
    # sa = [int(surfs_table[roi_ord.index(roi)]) for roi in rois]

    normalization = ( np.tile(sa, (np.size(sa), 1))
        + np.tile(sa, (np.size(sa), 1)).T )

    if normalization.shape != connection_matrix.shape:
        print('normalization.shape', normalization.shape)
        print('connection_matrix.shape', connection_matrix.shape)
        err_msg = 'Normaliztion matrix and connections matrix are not in the same shape!'
        logging.error(err_msg)
        raise Exception(err_msg)

    print 'finished normalizing'

    connection_matrix = connection_matrix * 2 / normalization 

    #np.save( op.join(subjdir, subj, 'linear_harmonic_connections.npy'), connection_matrix)
    #np.save( op.join(subjdir, subj, 'linear_normalization.npy'), normalization)

    np.fill_diagonal(connection_matrix, 0)
    np.save(linear_harmonic_output_fname, connection_matrix)

    labels_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp.txt' % parc)
    labels_all_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp_all.txt' % parc)
    if not op.isfile(labels_all_cmp_fname):
        shutil.copy(labels_cmp_fname, labels_all_cmp_fname)

    # with open(labels_all_cmp_fname, 'r') as lb_all_file:
    #     all_labels = lb_all_file.readlines()
    # with open(labels_cmp_fname, 'w') as lb_file:
    #     for index, label in enumerate(all_labels):
    #         if index in zeros_norm_rois_indices:
    #             print('write delete in the cmp file! {}, {}'.format(label, index))
    #             lb_file.write('delete\n')
    #         else:
    #             lb_file.write('{}'.format(label))


def namebase(file_name):
    return os.path.splitext(os.path.basename(file_name))[0]


def make_dir(fol):
    if not os.path.isdir(fol):
        os.makedirs(fol)
    return fol


def comp_adj_and_cmp_all(subjdir, subj, parc, parc_stem, dset_dir, parclist_dir, masks_dir, trk_stem='tracts20k'):
    labels_order_fname = op.join(parclist_dir, 'order_%s_alph.txt'%parc)
    new_labels_order_fname = op.join(subjdir, subj, 'label', 'order_%s_alph.txt'%parc)
    if not op.isfile(new_labels_order_fname):
        shutil.copy(labels_order_fname, new_labels_order_fname)
    labels_cmp_all = np.genfromtxt(new_labels_order_fname, dtype='str')
    labels_cmp_all = [s for s in labels_cmp_all if s != 'delete']
    # labels_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp.txt' % parc)
    # labels_cmp = np.genfromtxt(labels_cmp_fname, dtype='str')
    # labels_cmp = [s for s in labels_cmp if s != 'delete']
    connection_matrix_fname = op.join(dset_dir, 'tracks', 'linear_harmonic_%s_alph_%s.npy'%(parc, trk_stem))
    connection_matrix = np.load(connection_matrix_fname)
    centroids_dir = op.join(masks_dir, '{}_{}'.format(parc, parc_stem))
    dirents = get_dirents(centroids_dir, True)
    dirents = ['{}_{}'.format(d[-2:], d[:-3]) for d in dirents]
    dirents_minus_labels_cmp_all = set(dirents) - set(labels_cmp_all)
    if len(dirents_minus_labels_cmp_all) > 0:
        err_msg = 'labels in dirent and not in parclist!'
        logging.error(err_msg)
        raise Exception(err_msg)

    labels_in_connections_matrix = ['{1}_{0}'.format(*get_roi_hemi(op.basename(l)[:-len('.nii.gz')]))  for l in
                                    glob.glob(op.join(masks_dir, '{}_{}_bspace'.format(parc, parc_stem), '*.nii.gz'))]
    labels_cmp_all = ['{1}_{0}'.format(*get_roi_hemi(l)) for l in labels_cmp_all]
    labels_to_remove = set(labels_cmp_all) - set(labels_in_connections_matrix)
    if len(labels_to_remove) > 0:
        with open(new_labels_order_fname, 'r') as f:
            labels_cmp_data = f.read()
        for label_to_remove in labels_to_remove:
            labels_cmp_data = labels_cmp_data.replace(label_to_remove, 'delete')
        with open(new_labels_order_fname, 'w') as f:
            f.write(labels_cmp_data)
        labels_cmp_all = np.genfromtxt(new_labels_order_fname, dtype='str')
        labels_cmp_all = [s for s in labels_cmp_all if s != 'delete']
    if connection_matrix.shape[0] != len(labels_cmp_all):
        err_msg = 'connection_matrix and labels_cmp not in the same lengths! {} {}'.format(
            connection_matrix.shape[0], len(labels_cmp_all))
        logging.error(err_msg)
        raise Exception(err_msg)
    else:
        print('Connection_matrix and labels order have the same length! Hooray!')
        print('When loading the matrix in cvu, for adj matrix choose:')
        print(connection_matrix_fname)
        print('And for labels order:')
        print(new_labels_order_fname)


def comp_masks_and_cmp(parc, parc_stem, masks_dir, parclist_dir):
    labels_cmp_fname = op.join(parclist_dir, 'orders', '%s_cmp.txt' % parc)
    with open(labels_cmp_fname) as lb_file:
        labels = []
        for ln in lb_file:
            if 'delete' not in ln:
                labels.append(ln.strip())

    centroids_dir=op.join(masks_dir, '{}_{}'.format(parc, parc_stem))
    dirents = get_dirents(centroids_dir, True)
    dirents = ['{}_{}'.format(d[-2:], d[:-3]) for d in dirents]
    print(set(dirents) - set(labels))
    print(set(labels) - set(dirents))
    print(len(dirents), len(labels))


def get_bs_fols(subject, subjdir):
    diff_fol = op.join(subjdir, subject, 'diff')
    b_folders = glob.glob(op.join(diff_fol, '*7t_set*'))  + glob.glob(op.join(diff_fol, 'b*_set*'))
    b_folders = [fol for fol in b_folders if '_eddy.' not in fol]
    unique_bs = set([op.basename(fol).split('_')[0] for fol in b_folders])
    return b_folders, unique_bs


def create_dir_structure(subject, subjdir):
    b_folders, unique_bs = get_bs_fols(subject, subjdir)
    if not op.isfile(op.join(subjdir, subject, 'mri', 'orig.nii.gz')):
        mri_convert(op.join(subjdir, subject, 'mri', 'orig.mgz'),
                    op.join(subjdir, subject, 'mri', 'orig.nii.gz'))
    if not op.isfile(op.join(subjdir, subject, 'mri', 'wm.nii.gz')):
        mri_convert(op.join(subjdir, subject, 'mri', 'wm.mgz'),
                    op.join(subjdir, subject, 'mri', 'wm.nii.gz'))
    for b_val in unique_bs:
        logging.debug('creating files for {}'.format(b_val))
        dest_fol = op.join(subjdir, subject, 'diff', '{}_dest'.format(b_val))
        make_dir(dest_fol)
        create_link(op.join(subjdir, subject, 'diff', '{}_b0'.format(b_val), '{}_b0.nii.gz'.format(b_val)),
                    op.join(dest_fol, 'b0.nii.gz'.format(b_val)))
        create_link(op.join(subjdir, subject, 'mri', 'orig.nii.gz'),
                    op.join(dest_fol, 'mprage.nii.gz'))
        create_link(op.join(subjdir, subject, 'mri', 'orig.nii.gz'),
                    op.join(dest_fol, 'mprage_256.nii.gz'))
        create_link(op.join(subjdir, subject, 'mri', 'wm.nii.gz'),
                    op.join(dest_fol, 'wm_256.nii.gz'))
        b_fols = [fol for fol in b_folders if op.basename(fol).split('_')[0] == b_val]
        bvals, bvecs, data_images = [], [], []
        for b_fol in b_fols:
            b_fol_name = op.basename(b_fol)
            logging.debug('Loading data for {}'.format(b_fol_name))
            if not op.isfile(op.join(dest_fol, 'bvals')):
                bvals.append(np.genfromtxt(op.join(b_fol, '{}.bvals'.format(b_fol_name))))
            if not op.isfile(op.join(dest_fol, 'bvecs')):
                bvecs.append(np.genfromtxt(op.join(b_fol, '{}.bvecs'.format(b_fol_name))))
            if not op.isfile(op.join(dest_fol, 'data.nii.gz')):
                data_images.append(nib.load(op.join(b_fol, '{}.nii.gz'.format(b_fol_name))))
        # assert(data_image.shape[3] == bvecs.shape[0] == bvals.shape[0])
        if not op.isfile(op.join(dest_fol, 'bvals')):
            bvals = np.hstack(bvals)
            assert (bvals.ndim == 1)
            np.savetxt(op.join(dest_fol, 'bvals'), bvals, fmt='%.6f')
        if not op.isfile(op.join(dest_fol, 'bvecs')):
            bvecs = np.vstack(bvecs)
            assert (bvecs.shape[1] == 3)
            np.savetxt(op.join(dest_fol, 'bvecs'), bvecs, fmt='%.14f')
        if not op.isfile(op.join(dest_fol, 'data.nii.gz')):
            data_image = nib.concat_images(data_images, axis=3)
            nib.save(data_image, op.join(dest_fol, 'data.nii.gz'))


def mri_convert(org_fname, dest_fname):
    cmd = 'mri_convert {} {}'.format(org_fname, dest_fname)
    run_script_and_check_output(cmd, 'mri_convert')


def create_link(real_file, link_fname, overwrite=True):
    if overwrite and op.islink(link_fname):
        os.remove(link_fname)
    if not op.islink(link_fname):
        os.symlink(real_file, link_fname)
    if not op.islink(link_fname):
        err_msg = "Link wasn't created! Real file: {}, link name: {}".format(real_file, link_fname)
        logging.error(err_msg)
        raise Exception(err_msg)


def get_opt_val(cmd, opt='--output'):
    opt_res = [s for s in cmd.split(' ') if s.startswith(opt)]
    if len(opt_res) == 1:
        return opt_res[0][len(opt) + 1:]
    else:
        return ''


def get_n_jobs(n_jobs):
    import multiprocessing
    cpu_num = multiprocessing.cpu_count()
    n_jobs = int(n_jobs)
    if n_jobs > cpu_num:
        n_jobs = cpu_num
    elif n_jobs < 0:
        n_jobs = cpu_num + n_jobs
    return n_jobs


def main():
    import sys
    import getopt

    # if os.environ['SUBJECTS_DIR']:
    #     subjdir=os.environ['SUBJECTS_DIR']+'/'
    # else:
    #     subjdir='local_mount/space/truffles/1/users/rlaplant/data/qball'

    if len(sys.argv)==1:
        usage()
        return

    try:
        opts,args=getopt.getopt(sys.argv[1:],'s:m:',['skip=','nomorph',
            'shell=', 'trkstem=', 'nfibers=', 'no-fsdir', 'sbx=', 'n_jobs=',
            'preproc=','stem=','nmr=','nmrdir=', 'dset=', 'subjdir=', 'parclists='])
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
    subjdir = None
    parclist_dir = None
    n_jobs = -1

    for opt,arg in opts:
        if opt in ['--skip']:
            skip_stages.append(int(arg))
        if opt in ['--subjdir']:
            subjdir = arg
            os.environ['SUBJECTS_DIR'] = subjdir
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
            preproc_dir = op.join(subjdir, subj, arg)
        if opt in ['--preproc','--stem']:
            preproc_stem = arg
        if opt in ['--nmr','--nmrdir']:
            nmr_dir = arg
        if opt in ['--dset']:
            dset_name = arg
        if opt in ['--parclists']:
            parclist_dir = arg
        if opt in ['--n_jobs']:
            n_jobs = arg

    n_jobs = get_n_jobs(n_jobs)
    if nr_fibers != 20000 and dsistudiostem == 'tracts20k':
        raise ValueError("Must change DSIstudio save stem from default")

    if subjdir is None:
        if os.environ.get('SUBJECTS_DIR', None):
            subjdir = os.environ['SUBJECTS_DIR'] + '/'
        else:
            raise Exception('No subjdir!')

    if nmr_dir is None and subjdir is None:
        #nmr_dir = '/local_mount/space/truffles/1/users/rlaplant/data/hcpnmr/'
        nmr_dir = '/local_mount/space/truffles/2/users/rlaplant/data/nmrdb/'
    if no_separate_freesurfer_dir: 
        nmr_dir = subjdir
    if nmr_dir is None:
        raise Exception('No nmr_dir!')
    #nmr_dir = '/local_mount/space/truffles/1/users/rlaplant/data/epilepsy/'
    dwi_dir = subjdir

    masks_dir = op.join(dwi_dir, subj, 'masks')

    dset_dir = op.join(dwi_dir, subj, 'diff', dset_name)
    create_dir_structure(subj, subjdir)
    if not op.isdir(dset_dir):
        raise Exception('No dset dir, or dset dir does not exist!')

    if preproc_dir is not None:
        pass
    elif modality in ('qball','qbi','singleshell','sgsh','dti'):
        # preproc_dir = op.join(dwi_dir, subj, 'diff', 'preproc')
        preproc_dir = op.join(dset_dir, 'preproc')
    elif modality in ('dsi'):
        preproc_dir = op.join(dwi_dir, subj, 'sbx')
    elif modality == 'mask_only':
        preproc_dir = None

    for fname in ['b0.nii.gz', 'bvals', 'bvecs', 'data.nii.gz', 'mprage.nii.gz', 'mprage_256.nii.gz', 'wm_256.nii.gz']:
        if not op.isfile(op.join(dset_dir, fname)):
            raise Exception('{} does not exist in dest dir!'.format(fname))
    # dsi_dir = '/local_mount/space/truffles/1/users/rlaplant/data/dsi'

    if modality in ('singleshell','sgsh') and shell is None:
        sys.stderr.write('Specify shell as 1k, 3k, 5k, 10k only')
        sys.exit(85)

    if parclist_dir is None:
        parclist_dir = '/autofs/cluster/neuromind/npeled/parclists'
        if not op.isdir(parclist_dir):
            raise Exception('No parclist_dir!')

    # if not os.environ.has_key('DSI_STUDIO_FA_TEMPLATE'):
    #     sys.stderr.write("Source DSIstudio first!\n")
    #     sys.exit(85)
    if not os.environ.has_key('FREESURFER_HOME'):
        sys.stderr.write("Source Freesurfer first!\n")
        sys.exit(86)
    try:
        import nibabel
    except:
        raise Exception('No nibabel!')


    if not 0 in skip_stages:
        if modality in ('singleshell','sgsh'):
            organize_singleshell(subj, preproc_dir, shell)
    if not 1 in skip_stages:
        setup_mprage(subj, nmr_dir, dwi_dir, preproc_dir, dset_dir, modality)
    if not 2 in skip_stages:
        white_matter_mask(subj, subjdir, dset_dir)
    if not 3 in skip_stages:
        morph_and_extract_labels(subj, nmr_dir, subjdir, parclist_dir,
            masks_dir, dset_dir, nomorph=nomorph, n_jobs=n_jobs)
    if not 4 in skip_stages:
        bigmask_creation(subj, subjdir, masks_dir)
    if not 5 in skip_stages:
        run_dsistudio(subj, subjdir, masks_dir, preproc_dir, dset_dir,
            modality, shell=shell, nr_fibers=nr_fibers, 
            trk_output=dsistudiostem, preproc_stem=preproc_stem)
    if not 7 in skip_stages:
        linearly_corrected_matrix_creation(subj, subjdir, dset_dir, parclist_dir, masks_dir, trk_stem=dsistudiostem, n_jobs=n_jobs)
    if not 8 in skip_stages:
        parc = 'laus250'
        parc_stem = 'wmreg'
        # comp_masks_and_cmp(parc, parc_stem, masks_dir, parclist_dir)
        comp_adj_and_cmp_all(subjdir, subj, parc, parc_stem, dset_dir, parclist_dir, masks_dir)



if __name__=='__main__':
    # '-s nmr01002 -m qball --skip 0 --skip 1 --skip 2 --skip 3 --skip 4 --skip 5 --no-fsdir --dset b15_dest --subjdir /cluster/neuromind/dwakeman/tsc_pilot/subjects --parclists /cluster/neuromind/npeled/parclists --stem data'
    main()
