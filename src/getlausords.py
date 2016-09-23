import mne
import os.path as op

def mangle_hemi(s):
	return '%s_%s'%(s[-2:],s[:-3])


def create_cmp_alph_files(subjects_dir, parclists_dir, prc='laus250'):
    #prc='myaparc_%s' % nr
    import re
    num = int(''.join(re.findall('\d', prc)))
    txtf = 'labels_%s.txt' % num
    cvunm = prc

    with open(op.join(parclists_dir, txtf), 'r') as fd:
        cmpv=fd.readline()

    cmpv=cmpv.replace('.','_').replace('\'','').strip().split(', ')
    print cmpv


    prl=mne.read_labels_from_annot(subjects_dir=subjects_dir,subject='fsaverage5c',parc=prc,hemi='lh')
    prr=mne.read_labels_from_annot(subjects_dir=subjects_dir,subject='fsaverage5c',parc=prc,hemi='rh')
    pr=prl + prr
    pr=map(lambda lb:mangle_hemi(lb.name),pr)
    with open(op.join(parclists_dir, 'orders', '%s_alph.txt' % cvunm),'w') as fd:
        i=0
        for u in pr:
            if u in cmpv:
                fd.write('%s\n'%u)
                i+=1
            else:
                fd.write('delete\n')

    print i

    with open(op.join(parclists_dir, 'orders', '%s_cmp.txt' % cvunm),'w') as fd:
        i=0
        for u in cmpv:
            if u in pr:
                fd.write('%s\n'%u)
                i+=1
            else:
                fd.write('delete\n')
    print i

if __name__ == '__main__':
    subjects_dir = '/cluster/neuromind/npeled/subjects'
    create_cmp_alph_files(subjects_dir)