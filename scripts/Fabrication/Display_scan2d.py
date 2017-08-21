from analysis.lib.m2 import m2
reload(m2)
import numpy as np
from matplotlib import pyplot as plt
import os
import msvcrt

from analysis.lib.fitting import fit
from analysis.lib.fitting import common

class DisplayScan(m2.M2Analysis):

    def get_data(self):
        self.zfocus = self.f['instrument_settings']['master_of_space'].attrs['z']
        self.keyword = self.f['instrument_settings']['setup_controller'].attrs['keyword']
        self.xvalues = self.f['x'].value
        self.yvalues = self.f['y'].value
        self.countrates = self.f['countrate'].value
        '''self.fitresultx = self.f['fit_result']['x'].value
        self.fitresulty = self.f['fit_result']['y'].value'''
        return self.xvalues,self.yvalues,self.countrates

    def plot_data(self,title,save=True, **kw):
        vmax = kw.pop('vmax', np.amax(self.countrates))
        vmin = kw.pop('vmin', np.amin(self.countrates))

        use_save_location = kw.pop('use_save_location',False)
        save_location = kw.pop('save_location',None)
        colormap=kw.pop('colormap','afmhot')
        grid = kw.pop('grid',False)
        fig = plt.figure(**kw)
        ax = fig.add_subplot(1,1,1)

        
        #colors=ax.pcolormesh(self.xvalues,self.yvalues,self.countrates, vmin=vmin,vmax=vmax,cmap=colormap)
        colors=ax.imshow(self.countrates, cmap=colormap, interpolation='none',
                        vmin=vmin,vmax=vmax, 
                        extent=[np.amin(self.xvalues),np.amax(self.xvalues),np.amin(self.yvalues),np.amax(self.yvalues)], origin='lower')

        #ax.set_xlim([np.amin(self.xvalues),np.amax(self.xvalues)])
        #ax.set_ylim([np.amin(self.yvalues),np.amax(self.yvalues)])
        #ax.set_title('z = {:.2f}'.format(self.zfocus)+'_'+self.keyword+'_\n' + title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        if grid:
            ax.grid(color='r',linestyle='-')
        #ax3 = fig.add_subplot(1,1,1)
        # ax3.axis('off')
        plt.colorbar(colors,ax = ax, fraction=0.046, pad=0.04)
        #fig.colorbar(colors,ax=ax)
        plt.show()


        if save:
            try:
                fig.savefig(
                    os.path.join(self.folder,colormap+'_scan2d.png'))
                if use_save_location and (save_location!=None):
                    min_x = min(self.xvalues)
                    min_y = min(self.yvalues)
                    depth = self.folder[-4:]
                    print '%d,%d,%.2f,%s'%(min_x,min_y,self.zfocus,depth)
                    fig.savefig(
                        os.path.join(save_location, '%d,%d,%.2f%s.png'%(min_x,min_y,self.zfocus,depth)))
            except:
                print 'Figure has not been saved.'
        return fig,ax

    def find_NV_locations(self,g_sigma = 0.2,**kw):
        plot_NV_zoom=kw.pop('plot_NV_zoom',False)
        save_NVs = kw.pop('save_NVs',True)
        print_update = kw.pop('print_update',False)

        # guess parameters for the fit
        g_offset=np.average(self.countrates)
        g_A = self.countrates.max()
        g_sigmax = g_sigma
        g_sigmay = g_sigma
        g_theta=0

        #arrays of guesses for x0 and y0 - a fit attempt will be done for each
        xstepsize = g_sigmax*2 #the stepsize is 2*sigma, to ensure each NV is found
        ystepsize = g_sigmay*2
        g_x0s = np.arange(self.xvalues[0]+g_sigmax,self.xvalues[-1],xstepsize)
        g_y0s = np.arange(self.yvalues[0]+g_sigmay,self.yvalues[-1],ystepsize)

        #the fit is performed on a subset of the data around the intial guesses for x0 and y0
        patchsizex = int(len(self.xvalues)/abs(max(self.xvalues)-min(self.xvalues))*g_sigmax*3) #number of indices around g_x0
        patchsizey = int(len(self.yvalues)/abs(max(self.yvalues)-min(self.yvalues))*g_sigmay*3) #number of indices around g_x0
        # print patchsizex
        # print patchsizey

        #initialise some arrays or dictionaries to store results
        results=np.zeros((len(g_x0s),len(g_y0s)))
        x0s=np.array([0])
        y0s=np.array([0])
        self.NVs = {}
        NV_nr=0
        f=0

        total_nr_fit_attempts = len(g_x0s)*len(g_y0s)
        
        for i,g_x0 in enumerate(g_x0s):       

            i_min,i_max,x_around_x0 = self.crop_array_around_value(self.xvalues,g_x0,patchsizex)
            
            for j,g_y0 in enumerate(g_y0s):

                j_min,j_max,y_around_y0 = self.crop_array_around_value(self.yvalues,g_y0,patchsizey)
                #define meshgrid for x,y cropped around x0,y0
                meshx,meshy = np.meshgrid(x_around_x0,y_around_y0)
                data_x0y0 = self.countrates[j_min:j_max,i_min:i_max]
                
                #some progress indication for the impatient
                if print_update:
                    if f%100==0:
                        print f,'out of ',total_nr_fit_attempts,'attempts to fit completed'
                        print 'found', NV_nr, 'NVs so far'
                f+=1    
                
                #the actual fitting
                results[i,j],fitres = self.get_fit(data_x0y0,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigmax,g_sigmay,g_theta)
                #results[i,j],fitres = self.get_circular_fit(data,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigma) #general fit works better

                if results[i,j] == 1: #an NV was found
                    #first check if it is a unique NV centre, or a previously found one.
                    x0 = fitres['params_dict']['x0']
                    min_dx0 = x0s[np.argmin(abs(x0s-x0))]-x0 #minimal distance to all previously found x0s

                    y0 = fitres['params_dict']['y0']
                    min_dy0 = y0s[np.argmin(abs(y0s-y0))]-y0 #minimal distance to all previously found y0s                  

                    if ((abs(min_dx0)>g_sigmax) or (abs(min_dy0)>g_sigmay)): 
                        #if in either direction the NV is more than sigma away from a previously found one, it's unique
                        NV_nr += 1
                        'found NV! ', 'number', NV_nr
                        self.NVs['NV_'+str(NV_nr)] = fitres

                        if plot_NV_zoom:                      
                            extent = [meshx[0,0],meshx[-1,-1],meshy[0,0],meshy[-1,-1]]
                            fig,ax = plt.subplots()
                            ax.imshow(data_x0y0,extent = extent,origin='lower',interpolation='None')
                            plotname = 'NV_%d_x_%.2f_y_%.2f_y_'%(NV_nr,x0,y0)
                            title = self.folder+'\n'+plotname
                            ax.set_title(title)
                            try:
                                plt.savefig(os.path.join(self.folder,plotname+'.png'))
                            except:
                                print 'could not save fig'

                        # print 'x range',x_around_x0[0],g_x0,x_around_x0[-1] #for debugging
                        # print 'y range',y_around_y0[0],g_y0,y_around_y0[-1] #for debugging

                        #update x0s and y0s with our new NV
                        x0s=np.append(x0s,x0)
                        y0s=np.append(y0s,y0)

        if save_NVs:
            self.save_NVs_dict()

        return self.NVs


    def save_NVs_dict(self):
        f = self.analysis_h5data()
        if not 'NVs' in f:
            f.create_group('NVs')
        g = f['/NVs']

        for NV in self.NVs:
            if not NV in g:
                g.create_group(NV)

            for param in self.NVs[NV]['params_dict']:
                if param in g[NV]:
                    del g[NV][param]
                g[NV][param] = self.NVs[NV]['params_dict'][param]

            for error in self.NVs[NV]['error_dict']:
                if 'u_'+error in g[NV]:
                    del g[NV]['u_'+error]
                g[NV]['u_'+error]  = self.NVs[NV]['error_dict'][error]

            if 'reduced_chisq' in g[NV]:
                del g[NV]['reduced_chisq']
            g[NV]['reduced_chisq'] = self.NVs[NV]['reduced_chisq']                
            # for dictentry in NV:
            #     print  NVs[NV]
            #     if not dictentry in h:
            #         h.create_group(dictentry)
            #     g[NV][dictentry] = NVs[NV][dictentry]

        f.close()


    def plot_NV_locations(self, **kw):
        print_NV_data = kw.pop('print_NV_data',False)
        extent = (min(self.xvalues),max(self.xvalues), min(self.yvalues),max(self.yvalues))

        fig,ax = plt.subplots()
        c = ax.imshow(self.countrates, extent=extent,origin='lower',interpolation='None')

        for NV in self.NVs:
            fitres=self.NVs[NV]
            x0=fitres['params_dict']['x0']
            y0=fitres['params_dict']['y0']
            sigmax=fitres['params_dict']['sigmax']
            sigmay=fitres['params_dict']['sigmay']

            if print_NV_data:
                print 10*'*'
                print NV
                print 'x0 = %.2f um, y0 = %.2f um'%(x0,y0)
                print 'sigma_x = %.2f um, sigma_y = %.2f um'%(sigmax,sigmay)
                print 10*'*'
            ax.plot(x0,y0,'+',mew=2,markersize=20, label=NV)
        ax.legend()

        plt.colorbar(c)            
        plt.show()  

        plt.close()
        return ax

    def zoom_around_NV(self):
        for i in self.NVs:
            nr= i[3:]
            x0 = self.NVs[i]['params_dict']['x0']
            y0 = self.NVs[i]['params_dict']['y0']

            i_min,i_max,x_around_x0 = self.crop_array_around_value(self.xvalues,x0,11)

            j_min,j_max,y_around_y0 = self.crop_array_around_value(self.yvalues,y0,11)

            meshx,meshy = np.meshgrid(x_around_x0,y_around_y0)
            data_x0y0 = self.countrates[j_min:j_max,i_min:i_max]


            extent = [meshx[0,0],meshx[-1,-1],meshy[0,0],meshy[-1,-1]]
            fig,ax = plt.subplots()
            plotname = '_zoom_NV_%s_x_%.1f_y_%.1f'%(nr,x0,y0)
            ax.set_title(self.folder+'\n'+plotname)
            c=ax.imshow(data_x0y0,extent = extent,origin='lower',interpolation='None')
            plt.colorbar(c)
            plt.savefig(os.path.join(self.folder,plotname+'.png'))

            plt.close()

    def get_good_NVs_only(self):
        f = self.analysis_h5data()
        g = f['/NVs']

        parameters = ['A','a','x0','y0','sigmay','sigmax']
        reject=0
        good_NVs=np.array([])

        for NV in g:
            for p in parameters:
                if abs(g[NV]['u_'+p].value)>abs(0.5*g[NV][p].value):
                    # print p,g[NV][p].value,g[NV]['u_'+p].value
                    # print 'NV %s rejected'%(NV)
                    reject=1
                    break

            if 'reject' in g[NV]:
                del g[NV]['reject']
            g[NV]['reject'] = reject                
             
            if reject == 0:
                good_NVs = np.append(good_NVs,int(NV[3:]))
            reject=0

        f.close()
        return good_NVs

    def get_NV_information(self,NV_nr):
        f = self.analysis_h5data()
        g = f['/NVs']
        NV_dict = {}
        parameters = ['x0','y0','sigmay','sigmax']
        for p in parameters:
            NV_dict[p] = ( g['NV_'+str(int(NV_nr))][p].value, g['NV_'+str(int(NV_nr))]['u_'+p].value )

        f.close() 
        return NV_dict


    ##########################################################
    ##### auxiliary functions. maybe not necessarily within the class.
    ##########################################################
    def plot_fit2d(self,fitres,meshx,meshy):
        """
        General function that plots result from a 2d fit
        """
        extent = [meshx[0,0],meshx[-1,-1],meshy[0,0],meshy[-1,-1]] #we need this for plotting
        fig,ax1 = plt.subplots()
        ax1.imshow(fitres['fitfunc'](meshx,meshy),extent = extent,origin='lower')
        plt.show()


    def get_fit(self,data,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigmax,g_sigmay,g_theta,**kw):
        """
        fit a single 2d scan with a general 2d gaussian
        """
        do_plot = kw.pop('do_plot',False)
        
        fitres = fit.fit2d((meshx,meshy),data,common.fit_2d_gaussian, g_offset,g_A, g_x0, g_y0, g_sigmax, g_sigmay,g_theta, ret = True )
        if type(fitres)!=dict:  #the fit failed; no NV
            return 0,0

        sigmax  = fitres['params_dict']['sigmax']
        sigmay  = fitres['params_dict']['sigmay']
        A = fitres['params_dict']['A']

        parameters = ['x0','y0','sigmay','sigmax']
        for p in parameters:
            if abs(fitres['error_dict'][p])>abs(1.5*fitres['params_dict'][p]):
                return 0,0

        #rejection conditions
        ellipticity_condition = ((sigmax < sigmay*0.7) or (sigmax > sigmay*1.3)) #NV cannot be too elliptic
        amplitude_condition = A<(0.2*g_offset) #amplitude (S/N) cannot be too small
        size_condition = ((sigmax>0.5) or (sigmay>0.5)) #NV cannot have a diameter larger than 0.5um.
        conditions = (ellipticity_condition or amplitude_condition) or size_condition
        
        if conditions:
            return 0,0      
        # print fitres['params_dict']
        
        if do_plot:
            self.plot_fit2d(fitres,meshx,meshy)
        
        return 1, fitres


    def get_circular_fit(self,data,meshx,meshy,g_offset,g_A,g_x0,g_y0,g_sigma,**kw):
        """
        fit a single 2d scan with a circular 2d gaussian
        """
        do_plot = kw.pop('do_plot',False)

        fitres = fit.fit2d((meshx,meshy),data,common.fit_2d_gaussian_circular, g_offset,g_A, g_x0, g_y0, g_sigma, ret = True )
        if type(fitres)!=dict: #the fit failed; no NV
            return 0,0

        sigma  = fitres['params_dict']['sigma']
        A = fitres['params_dict']['A']

        #rejection conditions
        amplitude_condition = A<(0.2*g_offset) #amplitude cannot be too small
        size_condition = sigma>0.5 #NV cannot have a diameter larger than 0.5um.
        conditions = amplitude_condition or size_condition
        
        if conditions: 
            return 0,0

        if do_plot:
            plot_fit2d(fitres,meshx,meshy)
        
        return 1, fitres



    def crop_array_around_value(self,array,value,patchsize):
        """
        Generally useful function that crops an array to 2*patchsize around value
        """
        i_closest = np.argmin(abs(array-value))
        if (i_closest-patchsize)>0:
            i_min = i_closest-patchsize
        else:
            i_min = 0
        if (i_closest-patchsize)<len(array):
            i_max = i_closest+patchsize
        else:
            i_max = -1
        
        cropped_array = array[i_min:i_max]
        return i_min,i_max,cropped_array

def display_scan(older_than = None, nr_plots =1):
    """function that calls DisplayScan
    PARAMETERS:
    older_than - timestamp that identifies the first data (default = None; gives latest)
    nr_plot - give the total number of plots you want to show. default = 1."""

    for j in np.arange(nr_plots):
        timestamp,folder = tb.latest_data(contains='scan2d',older_than = older_than, return_timestamp=True)
        older_than=timestamp

        title = folder
        a = DisplayScan(folder)
        a.get_data()
        a.plot_data(title)
        a.f.close()



def table_of_good_NVs(folder):
    a = DisplayScan(folder)
    a.get_data() 
    goodNV = a.get_good_NVs_only()


class DisplayScanFlim(DisplayScan):

    def get_flim_data(self):
        self.flim_data = self.f['flim_data'].value
        #self.flim_data_syncs = self.f['flim_data_syncs'].value

    def plot_flim_data(self,title,save=True, **kw):
        self.zfocus = self.f['instrument_settings']['master_of_space'].attrs['z']
        self.xvalues = self.f['x'].value
        self.yvalues = self.f['y'].value
        self.keyword = self.f['instrument_settings']['setup_controller'].attrs['keyword']
        if kw.pop('raw', False):
            flim_data =  self.flim_data
            title=title+'_raw'
        else:
            flim_data =  self.flim_data.astype(np.float)/self.flim_data_syncs.reshape(np.append(self.flim_data.shape[:2],1)) * 1e4
        
        rmin = kw.pop('rmin', 0)
        rmax = kw.pop('rmax', len(flim_data[0,0]))
        print rmax
        title=title+'range: {} - {}'.format(rmin,rmax)

        title = title + ''
        self.flim_plot_data = np.sum(flim_data[:,:,rmin:rmax],axis=2)
        vmax = kw.pop('vmax', np.amax(self.flim_plot_data))
        vmin = kw.pop('vmin', np.amin(self.flim_plot_data))

        use_save_location = kw.pop('use_save_location',False)
        save_location = kw.pop('save_location',None)
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
        colorname='gist_earth'
        colors=ax.imshow(self.flim_plot_data, cmap=colorname, interpolation='none', vmin=vmin,vmax=vmax, extent=[np.amin(self.xvalues),np.amax(self.xvalues),np.amin(self.yvalues),np.amax(self.yvalues)], origin='lower')

        #ax.set_title('z = {:.2f}'.format(self.zfocus)+'_'+self.keyword+'_\n' + title)
        ax.set_xlabel('x', fontsize = '18' )
        ax.set_ylabel('y', fontsize = '18')
        ax3 = fig.add_subplot(1,1,1)
        #ax3,axis('on')
        #ax3.axis('off')
        plt.colorbar(colors,ax = ax, fraction=0.046, pad=0.04)
        #fig.colorbar(colors,ax=ax)
        plt.show()


        if save:
            try:
                fig.savefig(
                    os.path.join(self.folder,colorname+'_scan2dflim.png'))
            except:
                print 'Figure has not been saved.'

    def plot_flim_hist(self, title,save=True, **kw):

        if kw.pop('raw', False):
            flim_data =  self.flim_data_raw
        else:
            flim_data =  self.flim_data
        xmin = kw.pop('xmin', 0)
        xmax = kw.pop('xmax', len(flim_data[0]))
        ymin = kw.pop('ymin', 0)
        ymax = kw.pop('ymax', len(flim_data[:,1]))

        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(1,1,1)
        ax.semilogy(np.sum(flim_data[ymin:ymax,xmin:xmax,:],axis=(0,1)))
        ax.set_xlabel('bin')
        ax.set_ylabel(self.f.attrs['flim_units'],fontsize = '14', fontweight = 'bold')

        if save:
            try:
                fig.savefig(
                    os.path.join(self.folder,'_histflim.png'))
            except:
                print 'Figure has not been saved.'
        return ax

