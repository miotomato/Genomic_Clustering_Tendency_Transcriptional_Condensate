'''
Created on 06/17/2022
@authors: Zhenjia Wang <zhenjia@virginia.edu>, Chongzhi Zang <zang@virginia.edu>, 

'''

import os,re,argparse,sys
import bisect
# from IOparser_BedBam import get_tag_regions
#from BART.OptValidator import
import numpy as np
import pandas as pd


plus = re.compile('\+')
minus = re.compile('\-')
      


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY', 'chrM'];

mm10_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chrX','chrY', 'chrM']


	        	    
		
def get_count_on_mapID(start,end,positions,scores,mid,expand,half_life):
    '''
    Count the tags/positions on DHS
    '''
    #if is_list_sorted(positions)==0:
    #    positions.sort()
    if start < end and mid:
        middle = int((start+end)*0.5)
        s = bisect.bisect_left(positions,middle-expand)
        e = bisect.bisect_right(positions,middle+expand)        
        
        rp_score = 0.0
        for ii in range(s,e):
            x = abs(positions[ii]-middle)
            y = 2**(-1*x/half_life)*scores[ii]
            rp_score +=y
            
        return rp_score
        
    # elif start < end:
    #     s = bisect.bisect_left(positions,start-expand)
    #     e = bisect.bisect_right(positions,end+expand)
    #     return e-s, end-start+expand*2
            
    else:
        return 0




def get_bdg_regions(bdgfile,chroms,resolution):
    ''' 
    Get score from bedGraph files 
    '''  
    infile = open(bdgfile,'r')
    try:
        line = infile.readline()
    except:
        sys.stderr.write('Not a valid bedGraph format of file: {} ! \n\n'.format(bdgfile))
        sys.exit(1)
    position_list, score_list = {},{}
    while line:
        line = line.strip().split() 
        if line[0] in chroms and len(line)>=3:
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            score = float(line[3])
            # ==== ignore zero-scored regions
            if score!=0:
                # ==== empty list for new chrom
                if chrom not in position_list:
                    position_list[chrom]=[]
                    score_list[chrom]=[]
                # ==== evenly distributed score within start-end region
                for ii in np.arange(start,end,resolution):
                    position_list[chrom].append(ii)
                    score_list[chrom].append(score)
        else:
            pass
        line = infile.readline()
    infile.close()
    
    # ==== sort score by positions
    sorted_position_list, sorted_score_list = {},{}
    for chrom in position_list:
        sorted_position_list[chrom]=[]
        sorted_score_list[chrom]=[]
        for (x,y) in sorted(zip(position_list[chrom],score_list[chrom])):
            sorted_position_list[chrom].append(x)
            sorted_score_list[chrom].append(y)
    
    return sorted_position_list, sorted_score_list



def read_count_on_mapfile(afile,bdgfile,species,resolution,mid,expand,half_life):

    # ==== read score from bedGraph file
    chroms = hg38_chroms if species=='hg38' else mm10_chroms
    sorted_position_list, sorted_score_list = get_bdg_regions(bdgfile,chroms,resolution)
    
    counting = {}   
    mapfile = open(afile,'r')
    line = mapfile.readline()
    while line:
        line = line.strip().split()
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        try:
            map_id = line[3]
        except:
            map_id = '{}_{}_{}'.format(chrom,start,end)
        if chrom in sorted_position_list:
            rp_score = get_count_on_mapID(start,end,sorted_position_list[chrom],sorted_score_list[chrom],mid,expand,half_life)
            counting[map_id]= np.round(rp_score,2)
        else:
            counting[map_id]= 0
        line = mapfile.readline()
    mapfile.close()
    return counting


def main(afile,bdgfile,outfile,species,resolution,mid,expand,half_life):

    counting = read_count_on_mapfile(afile,bdgfile,species,resolution,mid,expand,half_life)
    
    with open(outfile,'w') as outf:
        outf.write('ID\tRP\n')
        for i in counting:
            outf.write('{}\t{}\n'.format(i,counting[i]))
    
        # outf.write('ID\tChr\tLeft\tRight\tRP\n')
        # for i in counting:
        #     chrom,start,end = i.split('_')[0],int(i.split('_')[1]),int(i.split('_')[2])
        #     middle = int((start+end)*0.5)
        #     outf.write('{}\t{}\t{}\t{}\t{}\n'.format(i,chrom,middle-expand,middle+expand,counting[i]))
    
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--afile', action = 'store', type = str,dest = 'afile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    parser.add_argument('-b', '--bdgfile', action = 'store', type = str,dest = 'bdgfile', help = 'bedGraph format', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile to write the the RPKM for each region', metavar = '<file>',required=True)
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
   
    parser.add_argument('-r', '--resolution', action = 'store', type = int,dest = 'resolution', help = 'resolution for the bedGraph file', default=50)
    # parser.add_argument('-g', '--fragmentsize', action = 'store', type = int,dest = 'fragmentsize', help = 'fragmentsize for the shift of reads. Default: 147', metavar = '<int>',default=147)
    parser.add_argument('-m', '--mid', action = 'store_true', dest = 'mid', help = 'whether to use middle side for expansion. Default: False',default=True)
    parser.add_argument('-e', '--expand', action = 'store', type = int,dest = 'expand', help = 'expand of regions for RP calculation. Default: 100k', metavar = '<int>',default=100000)
    parser.add_argument('-u', '--half_life', action = 'store', type = int,dest = 'half_life', help = 'expand of regions for RP calculation. Default: 10k', metavar = '<int>',default=10000)
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.afile,args.bdgfile,args.outfile,args.species,args.resolution,args.mid,args.expand,args.half_life)
