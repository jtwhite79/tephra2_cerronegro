
import os
import shutil
import multiprocessing as mp
import numpy as np
import matplotlib as mpl
import pylab

#--globals
wdir = ''
case = 'cerronegro_2reg'


#map_xlim = [515000,537000]
#map_ylim = [1373000,1392000]
#grid_step = 200
map_xlim = [500000, 545000]
map_ylim = [1370000,1400000]
grid_step = 100

obs_bins = np.arange(-4.0,3.5,0.5)
obin_centers = (obs_bins[:-1] + obs_bins[1:]) / 2.0
vent_x,vent_y = 532596.0,1381860.0
obs = os.path.join('background', 'cerronegro.xy')
obs_xys = np.loadtxt(obs, usecols=[0,1,2])

a_loc,b_loc = (0.15,0.65,0.8,0.25),(0.15,0.4,0.8,0.25)
t_loc,cb_loc = (0.15,0.15,0.8,0.25),(0.25,0.06,0.5,0.025)

min_mass = 1.0



def  resize(x,y,data):
    unique_x = np.sort(np.unique(x))
    unique_y = np.sort(np.unique(y))
    nrow,ncol = unique_y.shape[0],unique_x.shape[0]
    i_idx,j_idx = {},{}
    for j,ux in enumerate(unique_x):
        j_idx[ux] = j
    for i,uy in enumerate(unique_y):
        i_idx[uy] = i
    arr = np.zeros((nrow,ncol))-999
    for x,y,val in zip(x,y,data):
        
        arr[i_idx[y],j_idx[x]] = val
    return arr


def reduce(bin_idxs,out):
    reduced = []
    for ibin,idxs in enumerate(bin_idxs):
        o = out[:,idxs].sum(axis=1)
        reduced.append(o)
    return np.array(reduced)        


def load_rei():
    #--get the number of observation from the pst - use to cut out regul obs
    f = open(os.path.join(case,case+".pst"))
    for iline in range(4):
        line = f.readline()
    f.close()    
    nobs = int(line.strip().split()[1])
    #--load the residuals from each iteration
    if 'reg' in case:
        skip = 7
    else:
        skip = 5
    rei_dtype = np.dtype([('name', 'a12'), ('measured', np.float64),\
                     ('modeled', np.float64)])
    files = os.listdir(wdir+case+'/')
    rei_arrays = {}
    for f in files:
        raw = f.split('.')
        if len(raw) == 3 and raw[0] == case and raw[1] == 'rei':
            iter = int(raw[-1])
            rei = np.genfromtxt(os.path.join(wdir,case,f), skiprows=skip,\
             usecols=[0, 2, 3], dtype=rei_dtype)
            rei_arrays[iter] = rei[:nobs]
    if len(rei_arrays.keys()) == 0:
        if case+".rei" in files:
            rei = np.genfromtxt(os.path.join(wdir, case,case+".rei"), skiprows=skip,\
            usecols=[0, 2, 3], dtype=rei_dtype)
            rei_arrays[0] = rei[:nobs]            
    return rei_arrays


def plot_residual_grainhist():  
    rei_arrays = load_rei()
    iters = rei_arrays.keys()
    iters.sort()
    #--parse site names and grain size bins - keep as strings
    grain_bins,site_names = [],[]
    site_dict,grain_dict = {},{}
    for i,obs_name in enumerate(rei_arrays[iters[0]]["name"]):        
        sname,gbin = obs_name.split('_')
        gbin = gbin[2:]
        if "imas" not in obs_name:
            if sname not in site_names:
                site_names.append(sname)
                site_dict[sname] = [i]
            else:
                site_dict[sname].append(i)      
            if gbin not in grain_bins:
                grain_bins.append(gbin)
                grain_dict[gbin] = [i]
            else:
                grain_dict[gbin].append(i)

    #--plot grainsize historgrams for each site
    left = np.linspace(0, len(grain_bins)*2, len(grain_bins))
    width = 0.8
    for iiter in iters:
        rei = rei_arrays[iiter]
        fig = pylab.figure(figsize=(24,20))
        ifig,iax = 1,1
        for isite,[sname,idxs] in enumerate(site_dict.iteritems()):
            if isite % 24 == 0 and isite > 0:
                plt_name = case+"_grainhist_{1:02d}_iter{0:02d}.png".format(iiter,ifig)
                plt_path = os.path.join(wdir,case,"png",plt_name)
                pylab.savefig(plt_path)
                pylab.close(fig)
                fig = pylab.figure(figsize=(24,20))

                ifig += 1
                iax = 1
               
            ax = pylab.subplot(6,4,iax)
            ax.bar(left,rei['measured'][idxs],color='g',ec='none',label='measured')
            ax.bar(left+width,rei['modeled'][idxs],color='m',ec='none',label='modeled')
            ax.set_xticks(left+width)
            ax.set_xticklabels(grain_bins)
            ax.legend()
            ax.text(0.1,0.85,"sample "+str(sname),ha="left",\
                  transform=ax.transAxes)
            pylab.xticks(rotation=45)
            ax.set_ylim(0,500)
            iax += 1




def run_grid():
    import run_tephra2_stages as run
    #--write the grid xy file
    xs = np.linspace(map_xlim[0],map_xlim[1],grid_step)
    ys = np.linspace(map_ylim[0],map_ylim[1],grid_step)
    f = open("cerronegro_grid.xy",'w')
    for x in xs:
        for y in ys:
            f.write("{0:15.6E}  {1:15.6E}  0.0  \n".format(x,y))
    f.close()
    run.main(obs="cerronegro_grid.xy",out_a=os.path.join("out","stage_A_grid.out"),\
        out_b=os.path.join("out","stage_B_grid.out"),\
        out_t=os.path.join("out","total_grid.out"))
  

def plot_grid(a_array,b_array,t_array,vmin=min_mass,vmax=None):
    if vmin is None:
        vmin = min(t_array.min(),b_array.min(),a_array.min())
    if vmax is None:
        vmax = max(t_array.max(),b_array.max(),a_array.max())

    fig = pylab.figure(figsize=(5,10))    
    ax_a = pylab.axes(a_loc)
    pylab.yticks(rotation=45)
    pylab.xticks(rotation=45)
    ax_b = pylab.axes(b_loc)
    pylab.yticks(rotation=45)
    pylab.xticks(rotation=45)
    ax_t = pylab.axes(t_loc)
    pylab.yticks(rotation=45)
    pylab.xticks(rotation=45)
    ax_cb = pylab.axes(cb_loc)
    pylab.xticks(rotation=45)

    ax_a.imshow(np.flipud(a_array),extent=[map_xlim[0],map_xlim[1],\
               map_ylim[0],map_ylim[1]],\
               alpha=0.5,vmin=vmin,vmax=vmax)
    ax_a.scatter(obs_xys[:,1],obs_xys[:,2],marker='.')
    ax_a.scatter([vent_x],[vent_y],marker='^',s=200,color='k')
    ax_a.set_xlim(map_xlim)
    ax_a.set_ylim(map_ylim)
    ax_a.set_xticklabels([])
    
    ax_b.imshow(np.flipud(b_array),extent=[map_xlim[0],map_xlim[1],\
               map_ylim[0],map_ylim[1]],\
               alpha=0.5,vmin=vmin,vmax=vmax)
    ax_b.scatter(obs_xys[:,1],obs_xys[:,2],marker='.')
    ax_b.scatter([vent_x],[vent_y],marker='^',s=200,color='k')
    ax_b.set_xlim(map_xlim)
    ax_b.set_ylim(map_ylim)
    ax_b.set_xticklabels([])

    im = ax_t.imshow(np.flipud(t_array),extent=[map_xlim[0],map_xlim[1],\
               map_ylim[0],map_ylim[1]],\
               alpha=0.5,vmin=vmin,vmax=vmax)
    ax_t.scatter(obs_xys[:,1],obs_xys[:,2],marker='.')
    ax_t.scatter([vent_x],[vent_y],marker='^',s=200,color='k')
    ax_t.set_xlim(map_xlim)
    ax_t.set_ylim(map_ylim)

    cb = pylab.colorbar(im,cax=ax_cb,orientation='horizontal')
    ax_a.text(0.05,0.9,"Stage A",ha="left",\
              transform=ax_a.transAxes)
    ax_b.text(0.05,0.9,"Stage B",ha="left",\
              transform=ax_b.transAxes)
    ax_t.text(0.05,0.9,"Total",ha="left",\
              transform=ax_t.transAxes)
    return fig,ax_a,ax_b,ax_t



def run_plot_grid_maps():
    run_grid()
    
    fa= open(os.path.join("out","stage_A_grid.out"))
    fb= open(os.path.join("out","stage_B_grid.out"))
    ft= open(os.path.join("out","total_grid.out"))
    header = fa.readline()
    grainbins = header.strip().split()[4:]
    a_data = np.loadtxt(fa,skiprows=0)
    b_data = np.loadtxt(fb,skiprows=1)
    tot_data = np.loadtxt(ft,skiprows=1)
    x = np.sort(np.unique(a_data[:,0]))
    y = np.sort(np.unique(a_data[:,1]))
    X,Y = np.meshgrid(x,y)
    nrow,ncol = y.shape[0],x.shape[0]
    
    a_array = resize(a_data[:,0],a_data[:,1],a_data[:,3])
    b_array = resize(b_data[:,0],b_data[:,1],b_data[:,3])
    t_array = resize(tot_data[:,0],tot_data[:,1],tot_data[:,3])
    a_array = np.ma.masked_where(a_array<min_mass,a_array)
    b_array = np.ma.masked_where(b_array<min_mass,b_array)
    t_array = np.ma.masked_where(t_array<min_mass,t_array)
    fig,ax_a,ax_b,ax_t = plot_grid(a_array,b_array,t_array)
    ax_a.set_title("isomass")
    plt_name = case+"_grid_isomass.png"
    plt_path = os.path.join(wdir,case,"png",plt_name)
    pylab.savefig(plt_path)
    pylab.close(fig)

    cmap = pylab.get_cmap("jet")
    for ibin,gbin in enumerate(grainbins):
        a_array = resize(a_data[:,0],a_data[:,1],a_data[:,4+ibin])
        b_array = resize(b_data[:,0],b_data[:,1],b_data[:,4+ibin])
        t_array = resize(tot_data[:,0],tot_data[:,1],tot_data[:,4+ibin])
        a_array = np.ma.masked_where(a_array<min_mass,a_array)
        b_array = np.ma.masked_where(b_array<min_mass,b_array)
        t_array = np.ma.masked_where(t_array<min_mass,t_array)
        fig,ax_a,ax_b,ax_t = plot_grid(a_array,b_array,t_array)
        ax_a.set_title(gbin)  
        plt_name = case+"_grid_grainmap_"+gbin+".png" 
        plt_path = os.path.join(wdir,case,"png",plt_name)
        pylab.savefig(plt_path)
        pylab.close(fig)
    return



def plot_residual(rei): 
    resid = rei["measured"] - rei["modeled"]
    abs_resid = np.abs(resid)

    fig = pylab.figure(figsize=(5,15))
    ax_meas = pylab.axes(a_loc)
    pylab.yticks(rotation=45)
    pylab.xticks(rotation=45)
    ax_mod = pylab.axes(b_loc)
    pylab.yticks(rotation=45)
    pylab.xticks(rotation=45)
    ax_resid = pylab.axes(t_loc)
    pylab.yticks(rotation=45)
    pylab.xticks(rotation=45)
    
    ab_color = 'm'
    a_color = 'g'
    b_color = 'c'
    m_color = 'b'
    colors = []

    for snum,x,y in obs_xys:
        name_str = ""
        tot_resid = 0 #since some locations have A and B
        for oname,res in zip(rei["name"],resid):
            sname = oname.split('_')[0]
            #print sname
            try:
                snumm = int(sname[:2])
                if snum == snumm:
                    name_str += sname.upper()+' '        
                    tot_resid += res     
            except:
                continue
        #print name_str
        if "A" in name_str and "B" in name_str:
            colors.append(ab_color)
            
        elif "A" in name_str and "B" not in name_str:
            colors.append(a_color)
        elif "B" in name_str and "A" not in name_str:
            colors.append(b_color) 
        else:
            colors.append(m_color)          
        ax_mod.text(x,y,str(int(snum)),ha="center",va="center")            
        ax_meas.text(x,y,str(int(snum)),ha="center",va="center") 
        if tot_resid > 0:
            label = '+'
        else:
            label = '-'               
        ax_resid.text(x,y,label,ha="center",va="center")   

    norm_meas = rei["measured"]/rei["measured"].max()
    norm_mod = rei["modeled"]/rei["measured"].max()
    norm_resid = abs_resid / rei["measured"].max()
    min_size = 0
    meas_size = 500*norm_meas + min_size
    mod_size = 500*norm_mod + min_size
    resid_size = 500*norm_resid + min_size
    ax_meas.scatter(obs_xys[:,1],obs_xys[:,2],marker="o",s=meas_size,\
        color=colors,alpha=0.5)
    ax_mod.scatter(obs_xys[:,1],obs_xys[:,2],marker="o",s=mod_size,\
        color=colors,alpha=0.5)
    ax_resid.scatter(obs_xys[:,1],obs_xys[:,2],marker="o",s=resid_size,\
        color=colors,alpha=0.5)

    vent_x,vent_y = [532596.0],[1381860.0]
    ax_meas.scatter(vent_x,vent_y,marker='^',s=200,color='k')
    ax_mod.scatter(vent_x,vent_y,marker='^',s=200,color='k')
    ax_resid.scatter(vent_x,vent_y,marker='^',s=200,color='k')
    
    #--some trickery to get the legend symbols
    xlim,ylim = ax_meas.get_xlim(),ax_meas.get_ylim()
    ax_meas.scatter([0],[0],marker='o',s=70,color=ab_color,label="A+B")
    ax_meas.scatter([0],[0],marker='o',s=70,color=a_color,label="A")
    ax_meas.scatter([0],[0],marker='o',s=70,color=b_color,label="B")
    ax_meas.scatter([0],[0],marker='o',s=70,color=m_color,label="Mix")
    ax_meas.set_xlim(xlim)
    ax_meas.set_ylim(ylim)
    ax_mod.set_xlim(xlim)
    ax_mod.set_ylim(ylim)
    ax_resid.set_xlim(xlim)
    ax_resid.set_ylim(ylim)
    ax_meas.ticklabel_format(useOffset=False)
    ax_mod.ticklabel_format(useOffset=False)
    ax_resid.ticklabel_format(useOffset=False)
    ax_meas.set_ylabel("Northing ($m$)")
    ax_mod.set_ylabel("Northing ($m$)")
    ax_resid.set_ylabel("Northing ($m$)")
    ax_resid.set_xlabel("Easting ($m$)")
    ax_meas.set_xticklabels([])
    ax_mod.set_xticklabels([])
    ax_meas.grid()
    ax_mod.grid()
    ax_resid.grid()
    ax_meas.legend(loc=2,scatterpoints=1,frameon=False)
    
    ax_meas.text(0.75,0.925,"Measured",ha="left",\
              transform=ax_meas.transAxes)
    mmax = max(rei["measured"].max(),rei["modeled"].max())
    ax_meas.text(0.45,0.075,"Max={0:10.3E} $kg/m^2$".format(mmax),ha="left",\
              transform=ax_meas.transAxes)
    ax_mod.text(0.75,0.925,"Modeled",ha="left",\
              transform=ax_mod.transAxes)
    ax_resid.text(0.75,0.925,"Residual",ha="left",\
              transform=ax_resid.transAxes)
    ax_resid.text(0.45,0.075,"RSS={0:10.3E}".format(np.sum(resid**2)),ha="left",\
              transform=ax_resid.transAxes)
    #print np.dot(resid,resid),np.sum(resid**2)
    return fig,ax_meas,ax_mod,ax_resid



def plot_residual_maps():
    rei_arrays = load_rei()
    iters = rei_arrays.keys()
    iters.sort()
    #--get indices of individual grain size bins and isomass
    otype_dict = {}
    for oidx,oname in enumerate(rei_arrays[iters[0]]['name']):
        #print oname
        otype = oname.split('_')[1]
        if otype not in otype_dict.keys():
            otype_dict[otype] = [oidx]
        else:
            otype_dict[otype].append(oidx)

    for iiter in iters:
        for otype,idxs in otype_dict.iteritems():
            rei = rei_arrays[iiter][idxs]
            #print len(idxs),rei.shape
            fig,ax_meas,ax_mod,ax_resid = plot_residual(rei)
            ax_meas.set_title("iteration: "+str(iiter))
            plt_name = case+"_residual_"+otype+"_iter{0:03d}.png".format(iiter)
            plt_path = os.path.join(wdir,case,"png",plt_name)
            pylab.savefig(plt_path)
            pylab.close(fig)



def plot_residual_loe():
    rei_arrays = load_rei()
    iters = rei_arrays.keys()
    iters.sort()
    #--get indices of individual grain size bins and isomass
    otype_dict = {}
    for oidx,oname in enumerate(rei_arrays[iters[0]]['name']):
        #print oname
        otype = oname.split('_')[1]
        if otype not in otype_dict.keys():
            otype_dict[otype] = [oidx]
        else:
            otype_dict[otype].append(oidx)

    for iiter in iters:
        for otype,idxs in otype_dict.iteritems():
            rei = rei_arrays[iiter][idxs]
            fig = pylab.figure(figsize=(10,10))
            ax = pylab.subplot(111)
            ax.scatter(np.sqrt(rei["measured"]),np.sqrt(rei["modeled"]),color='b',marker='.')
            xmin,xmax = np.sqrt(rei["measured"].min()),np.sqrt(rei["measured"].max())
            ax.plot([xmin,xmax],[xmin,xmax],'k-') 
            ax.set_ylim(xmin,xmax)
            ax.set_xlim(xmin,xmax)
            ax.set_xlabel("Measured ($kg/m^2$)")
            ax.set_ylabel("Modeled ($kg/m^2$)")
            ax.set_title(otype) 
            plt_name = case+"_loe_"+otype+"_iter{0:03d}.png".format(iiter)
            plt_path = os.path.join(wdir,case,"png",plt_name)  
            pylab.savefig(plt_path)
            pylab.close(fig)


def plot_pars():
    #--get a list of parameter files by iteration
    #-- and load the parameter values into a dict
    par_values = {}
    case_dir = os.path.join(wdir,case)
    files = os.listdir(case_dir)
    iters,par_files  = [],[]
    for f in files:
        raw = f.split('.')
        if len(raw) == 3 and raw[1] == 'par':
            par_files.append(f)
            iiter = int(raw[2])
            iters.append(iiter)
            pvals = np.loadtxt(os.path.join(case_dir,f),usecols=[1],skiprows=1)
            par_values[iiter] = pvals
    
    #--load parameter initial values and bounds
    pst_file = open(os.path.join(wdir,case,case+'.pst'),'r')
    init_pvals,lbnd,ubnd = {},{},{}
    while True:
        line = pst_file.readline()
        if len(line) == 0:
            raise Exception()
        if line.startswith("* parameter data"):
            while True:
                line = pst_file.readline()
                if line.startswith('*'):
                    break
                raw = line.strip().split()
                pname,pval = raw[0],float(raw[3])    
                lb,ub = float(raw[4]),float(raw[5])
                if 'fixed' not in raw[1].lower():
                    pname = raw[0]
                    init_pvals[pname] = pval
                    lbnd[pname] = lb
                    ubnd[pname] = ub
            break      
    pst_file.close() 
    iters.sort()
            
    #--load one wind file to get the par names and find nwind
    pfile = open(os.path.join(case_dir,par_files[0]))
    pfile.readline()
    pnames = []
    nwind = 0
    for line in pfile:
        raw = line.strip().split()
        pname,pval = raw[0],float(raw[1])
        if "ws_a" in pname:
            nwind += 1
        pnames.append(pname)
  
    #--scalar mappable for color dependence
    norm = pylab.Normalize(vmin=0,vmax=nwind)
    sm = pylab.cm.ScalarMappable(cmap=pylab.get_cmap("summer"),norm=norm)
    colors = []
    for iwind in range(nwind):
        colors.append(sm.to_rgba(iwind))



    for iiter in iters:
        pvals = par_values[iiter]

        fig = pylab.figure(figsize=(10,10))
        ax_winda = pylab.subplot2grid((2,2),(0,0),polar=True)
        ax_windb = pylab.subplot2grid((2,2),(0,1),polar=True)
        ax_tbl = pylab.subplot2grid((2,2),(1,0),colspan=2)
        ax_winda.set_title("Stage A Wind")
        ax_windb.set_title("Stage B Wind")
        ax_tbl.set_title("Adjustable Parameter Values")

        ax_windb.set_theta_direction(-1)
        ax_windb.set_theta_zero_location('N')
        ax_windb.set_xticklabels(['', '', 'W', '', '', '', 'E', ''])
        ax_winda.set_theta_direction(-1)
        ax_winda.set_theta_zero_location('N')
        ax_winda.set_xticklabels(['', '', 'W', '', '', '', 'E', ''])

        #--stage a wind
        mags = []
        wind_idxs = []
        for iwind in range(nwind):
            wd,ws = 'wd_a_{0:02d}'.format(iwind),'ws_a_{0:02d}'.format(iwind)
            iwd,iws = pnames.index(wd),pnames.index(ws)
            wind_idxs.extend([iwd,iws])
            twd = np.pi/180*(pvals[iwd])
            tws = pvals[iws]
            #print iiter,iwind,twd,tws,pvals[iwd]
            mags.append(tws)
            ax_winda.annotate('', xy=(twd,tws),  xycoords='data',
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='->',color=colors[iwind]))       

        #--stage b wind        
        wind_idxs = []
        for iwind in range(nwind):
            wd,ws = 'wd_b_{0:02d}'.format(iwind),'ws_b_{0:02d}'.format(iwind)
            iwd,iws = pnames.index(wd),pnames.index(ws)
            wind_idxs.extend([iwd,iws])
            twd = np.pi/180*(pvals[iwd])
            tws = pvals[iws]
            #print iiter,iwind,twd,tws,pvals[iwd]
            mags.append(tws)
            ax_windb.annotate('', xy=(twd,tws),  xycoords='data',
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='->',color=colors[iwind]))      
        max_ws = max(mags)
        #--add the upper and lower bound    
        ax_windb.annotate('', xy=(np.pi/180*(lbnd[wd]),max_ws),  xycoords='data',
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='-',linestyle='dashed',color='0.1'))
        ax_windb.annotate('', xy=(np.pi/180*(ubnd[wd]),max_ws),  xycoords='data',
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='-',linestyle='dashed',color='0.1'))
        ax_winda.annotate('', xy=(np.pi/180*(lbnd[wd]),max_ws),  xycoords='data',
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='-',linestyle='dashed',color='0.1'))
        ax_winda.annotate('', xy=(np.pi/180*(ubnd[wd]),max_ws),  xycoords='data',
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='-',linestyle='dashed',color='0.1'))
        ax_windb.set_rmax(max_ws)     
        ax_winda.set_rmax(max_ws)

        #--boring parameter table
        ax_tbl.set_axis_off()
        row_labels,row_idxs = [],[]
        for iname,pname in enumerate(pnames):
            if not pname.startswith("wd") and not pname.startswith("ws") \
                and pname in init_pvals.keys():
                row_labels.append(pname)
                row_idxs.append(iname)

        col_labels = ['Lower Bound','Upper Bound','Initial Value','Fit Value']
        row_labels.sort()
        cell_data = []
        cell_fmt = '{0:15.3G}'
        for par in row_labels:
            idx = pnames.index(par)
            par_data = [cell_fmt.format(lbnd[par]),\
                        cell_fmt.format(ubnd[par]),\
                cell_fmt.format(init_pvals[par]),\
                cell_fmt.format(pvals[idx])]
            cell_data.append(par_data)
        t = ax_tbl.table(loc='center',cellText=cell_data,\
            rowLabels=row_labels,colLabels=col_labels,\
            colWidths=[0.25,0.25,0.25,0.25])
        t.set_fontsize(13)
        table_props = t.properties()
        table_cells = table_props['child_artists']
        for cell in table_cells: cell.set_height(0.075)
        ax_tbl.text(0.1,-0.1,"Note: alpha = 3.0",transform = ax_tbl.transAxes)

        plt_name = os.path.join(case_dir,"png",case+\
                "_parameters_iter{0:02d}.png".format(iiter))    
        pylab.savefig(plt_name)    
        



def run_pest():
    case_dir = os.path.join(case)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.mkdir(case_dir)
    os.mkdir(os.path.join(case,"png"))
    shutil.copy2(case+".pst",os.path.join(case,case+'.pst'))
   
    exe = os.path.join("..","exe","pest")
    pst_path = os.path.join(case,case)
    os.system(exe+" "+pst_path)


def clean_figs():
    fig_dir = os.path.join(case,"png")
    files = os.listdir(fig_dir)
    for f in files:
        if f.startswith(case):
            os.remove(os.path.join(fig_dir,f))


def main():
    
    run_pest()
    clean_figs()
    
    #plot_pars()
    #plot_residual_loe()
    #plot_residual_grainhist()
    #plot_residual_maps()
    #run_plot_grid_maps()
    
    p_hist = mp.Process(target=plot_residual_grainhist)
    p_rmaps = mp.Process(target=plot_residual_maps)
    p_loe = mp.Process(target=plot_residual_loe)
    p_gmaps = mp.Process(target=run_plot_grid_maps)
    p_hist.start()
    p_rmaps.start()
    p_loe.start()
    p_gmaps.start()
    return


if __name__ == '__main__':
    main()

