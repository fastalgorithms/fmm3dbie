from numpy import *
from pylab import *

com_str = '--------'

# list of subroutines
subs_list = []

# file name in which subroutine can be found
files_list = []

# name of raw files where subroutine documentation should be output
out_list = []


subs_list.append('subroutine lpcomp_helm_comb_dir')
files_list.append('../src/helm_wrappers/helm_comb_dir.f')
out_list.append('raws/lpcomp-helm-comb-dir.raw')

for i in range(size(subs_list)):
    f = open(files_list[i],'r')
    g = open(out_list[i],'w')
    istart = 0
    iend = 0
    for line in f:
        if(line.find(subs_list[i]) !=-1):
            istart = 1
        if(istart == 1 and line.find(com_str) !=-1):
            istart = 2
            iend = 1
            continue
        if(istart == 2 and iend >0):
            if(line.find(com_str) !=-1):
                iend = 0
                g.close()
                f.close()
                break
            if(iend > 0):
                if(len(line)<3):
                    g.writelines('\n')
                g.writelines(line[3::])

            
                

        

            
            
            

