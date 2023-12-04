#!/home/bol7asan/anaconda3/envs/bokeh/bin/python

import glob
import os
from bokeh.palettes import Category20_20 as colors
from datetime import date



today = date.today()

def hex_to_rgb(hexcode, alpha=0.1):
    hexcode = hexcode.lstrip('#')
    rgbcode = ','.join([str(int(hexcode[i:i+2], 16)) for i in (0, 2, 4)])
    return str(rgbcode)


bwfiles=[os.path.basename(x) for x in glob.glob(snakemake.params.bwDir+'/*.bw')]

#Create trackDb folder
fp=open(snakemake.output[0],'w')
for _nn, bwfile in enumerate(bwfiles):
    fp.write('track '+bwfile.split('.bw')[0]+'\n')
    fp.write('bigDataUrl '+bwfile+'\n')
    fp.write('shortLabel '+bwfile+'\n')
    fp.write('longLabel '+bwfile+'\n')
    fp.write('autoScale on \n')
    fp.write('type bigWig \n')
    fp.write('visibility full \n')
    fp.write('shortLabel '+bwfile+'\n')
    fp.write('color '+hex_to_rgb(colors[_nn])+'\n\n')
fp.close()


#Create hub file 
fp=open(snakemake.output[1],'w')
fp.write('hub ChIP-Seq_' + str(today) +'\n')
fp.write('shortLabel ChIP-Seq hub\n')
fp.write('longLabel ChIP-Seq hub\n')
fp.write('genomesFile genomes.txt\n')
fp.write('email mohammed.alhusayan@pennmedicine.upenn.edu\n')
fp.write('descriptionUrl none')
fp.close()

#Create genomes file
fp=open(snakemake.output[2],'w')
fp.write('genome ' + snakemake.params.build + '\n')
fp.write('trackDb ' + snakemake.params.build + '/trackDb.txt')
fp.close()
