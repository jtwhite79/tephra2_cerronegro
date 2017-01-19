import os
import sys
import multiprocessing as mp
import numpy as np
import pylab

#exe = os.path.join("..","exe","tephra2")
exe = "tephra2"

def stage_a(obs,out_a,verbose):
    #os.system("../exe/tephra2 conf/stage_A.conf "+obs+" wind/stage_A.wind 1>"+out_a+" 2>nul")       
    conf = os.path.join("conf","stage_A.conf")
    wind = os.path.join("wind","stage_A.wind")
    
    if verbose:
        cmd_line = exe+" "+conf+" "+obs+" "+wind+" 1>"+out_a    
    else:   
        cmd_line = exe+" "+conf+" "+obs+" "+wind+" 1>"+out_a+" 2>nula"
    print cmd_line
    os.system(cmd_line)
            
    return
def stage_b(obs,out_b,verbose):
    #os.system("../exe/tephra2 conf/stage_B.conf "+obs+" wind/stage_B.wind 1>"+out_b+" 2>nul")    
    conf = os.path.join("conf","stage_B.conf")
    wind = os.path.join("wind","stage_B.wind")    
    if verbose:
        cmd_line = exe+" "+conf+" "+obs+" "+wind+" 1>"+out_b    
    else:   
        cmd_line = exe+" "+conf+" "+obs+" "+wind+" 1>"+out_b+" 2>nulb"
    print cmd_line
    os.system(cmd_line)
    return


def reduce(bin_idxs,out):
    reduced = []
    for ibin,idxs in enumerate(bin_idxs):
        o = out[:,idxs].sum(axis=1)
        reduced.append(o)
    return np.array(reduced).T        

def main(obs=os.path.join("base","cerronegro.obs"),\
    out_a=os.path.join("out","stage_A.out"),\
    out_b=os.path.join("out","stage_B.out"),\
    out_t=os.path.join("out","total.out"),\
    verbose=False,parallel=True):
    if parallel:    
        p_a = mp.Process(target=stage_a,args=(obs,out_a,verbose))
        p_b = mp.Process(target=stage_b,args=(obs,out_b,verbose))
        p_a.start()
        p_b.start()
        p_a.join()
        p_b.join()
    else:
        stage_a(obs,out_a,verbose)
        stage_b(obs,out_b,verbose)
        
    
    obs_bins = np.arange(-4.0,3.5,0.5)
    obin_centers = (obs_bins[:-1] + obs_bins[1:]) / 2.0
    bin_labels = ["[-999:"+str(obs_bins.min())+')']
    for lbin,ubin in zip(obs_bins[:-1],obs_bins[1:]):
        bin_labels.append('['+str(lbin)+":"+str(ubin)+')')
    bin_labels.append(':'+str(obs_bins.max())+'+')

    
    out_files = [out_a,out_b]
    
    f = open(out_files[0])
    header = f.readline()
    f.close()
    grainbins = header.strip().split()[4:]
    #ubnd,lbnd = [],[]
    bin_idxs = []
    for ib in range(len(obs_bins)+1):
        bin_idxs.append([])
    for igb,gb in enumerate(grainbins):
        raw = gb.split(':')
        l,u = float(raw[0][1:]),float(raw[1][:-1])
        center = (u+l)/2.0
        #lbnd.append(l)
        #ubnd.append(u)
        if center < obs_bins[0]:
            bin_idxs[0].append(igb)
        elif center > obs_bins[-1]:
            bin_idxs[-1].append(igb) 
        else:
        #--find the near bin center
            dist = np.abs(obin_centers - center)
            idx = np.argmin(dist)
            bin_idxs[idx+1].append(igb)
              
    #nbins = 1 #count the last greater than bin
    #for ibin,obin in enumerate(obs_bins):
    #    idxs = bin_idxs[ibin]
    #    nbins += len(idxs) 

    # for ibin,idxs in enumerate(bin_idxs):
    #     print bin_labels[ibin],idxs
    # return    

    f = open(out_files[0])
    header = f.readline()
    tot = np.loadtxt(f) 

    f.close()                          
    #for ibin in range(4,tot.shape[1]):
    #    tot[:,ibin] /= tot[:,3]


    reduced = reduce(bin_idxs,tot[:,4:])     
    tot = np.hstack((tot[:,:4],reduced)) 
     

    #--resave outfile[0] with the reduced grainsize columns
    f = open(out_files[0],'w')
    f.write("x  y  elev mass "+" ".join(bin_labels)+"\n")
    np.savetxt(f,tot,fmt="%20.8E",delimiter='')

    #--process any remaining out files    
    for out_file in out_files[1:]:  
            
        f = open(out_file)
        header = f.readline()
        out = np.loadtxt(f)        
        f.close()
        #for ibin in range(4,out.shape[1]):
        #    out[:,ibin] /= out[:,3]        
        reduced = reduce(bin_idxs,out[:,4:])
        
        out = np.hstack((out[:,:4],reduced)) 
        #for row in out:
         #   print row[3],row[4:].sum()
        f = open(out_file,'w')
        f.write("x  y elev mass "+" ".join(bin_labels)+"\n")
        np.savetxt(f,out,fmt="%20.8E",delimiter='')

        
        #--add data columns (including the total mass)
        tot[:,3:] += out[:,3:]

    f = open(out_t,'w')
    f.write("x  y elev  mass "+" ".join(bin_labels)+"\n")
    np.savetxt(f,tot,fmt="%20.8e",delimiter='')
    f.close()
    #for row in tot:
     #  print row[3],row[4:].sum()
    return   

if __name__ == "__main__":
    main(verbose=True,parallel=False)    