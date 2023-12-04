import os

inp_dir = '/media/asangani1/CDK7_project/2023-11-28/pre-analysis/'

def create_bam_dir(inp_dir):
    try:
        os.system('mkdir ' + os.path.join(inp_dir,'bam_files'))
    except:
        pass
    samples = os.listdir(inp_dir)
    samples = [x for x in samples if '_files' not in x and 'spiker' not in x]
    for sample in samples:
        os.system('ln -s ' + os.path.join(inp_dir,sample,'bowtie2_dros','aligned.primary.rmdup.bam') + ' ' +os.path.join(inp_dir,'bam_files',sample+'.bam'))

    return os.path.join(inp_dir,'bam_files')
    

def split_files(inp_dir,bam_dir);
    try:
        os.system('mkdir ' + os.path.join(inp_dir,'spiker'))
    files = os.listdir(bam_dir)
    for file in files:
        os.system('split_bam.py --threads 8 -i ')


if __name__=='__main__':
    bam_dir = create_bam_dir(inp_dir)
