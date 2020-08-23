"""
Created on Tue Jan  8 20:35:23 2019

@author: roger
"""
############ ler os resultados ############
dire = r'E:\Documents\EM USO\INPE\2018\00_Dissertacao\03_Analitico\Qaa_rogerio\qaapy\Results_OLCI'
######### Juntar os dados #################

def merge_param_dict(chdir,param):
    import os
    import glob
    param_dict = {}
    lista = glob.glob(chdir + '\\*.pickle')
    for i in lista:
        barcut = i.split('\\')[-1]
        icut = barcut.split('.')
        df = pd.read_pickle(i)
        param_dict[icut[0]] = df[param]
    return param_dict

########## Reordenar elementos ###########################
    
#import collections
#
#def move_element(odict, thekey, newpos):
#    odict[thekey] = odict.pop(thekey)
#    i = 0
#    for key, value in odict.items():
#        if key != thekey and i >= newpos:
#            odict[key] = odict.pop(key)
#        i += 1
#    return odict



from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['arial']})
rc('text', usetex=True)

############### Plotar ##################################

Plt_atALL(at_Dict,a_ACS_water,709)
Plt_atALL(adg_dict,data_obj.aCDM,709)
Plt_atALL(aphy_dict,data_obj.aphy,709)

def Plt_atALL(at_dict,acs,maxi):
    '''
    Function to plot lines from X,Y data
     Especific for ROGERIO's DATA
    '''
    #acs.index.name = 'Wavelength'
    #with plt.style.context(('fivethirtyeight')):
    with plt.style.context(('seaborn')):
        fig, axes = plt.subplots(nrows=3, ncols=3)#, sharex=True, sharey=True)
        
        #axes.text(0.1,0.5,"hahaha")
        #fig.suptitle('{} - Lago Curuai - Lb0 {} - eq {} '.format(ano,lb0,eq))
        #fig.add_subplot(111)
        axes[0,0].plot(acs, lw=1)
        axes[0,0].set_xlim((380,720))
        axes[0,0].set_ylim((-3,15))
        axes[0,0].set_title('a) ' + 'Medido', fontsize=20)
        axes[0,0].xaxis.set_tick_params(labelsize=14)
        axes[0,0].yaxis.set_tick_params(labelsize=14)
        
        #axes[0,0].set_ylabel('at (m-1)')
        #axes[0,0].set_xlabel('Wavelenght (nm)')
        #acs.plot(ylim=(0,22),legend=False,title='Total absorption ACS',ax=axes[0,0],lw=1)
        
        #xes = [(1,1),(1,0),(0,1),(2,0),(2,1),(3,0),(3,1)]
        xes = [(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]
        abcd = ['b) ','c) ','d) ','e) ','f) ','g) ','h) ','i) ']
        #xes = [(0,1),(1,0),(1,1)]
        for x,i,abc in zip(at_dict,xes,abcd):
            
            nome = x.split('_')[2]            
            axes[i].plot(at_dict[x].loc[:maxi],lw=1)
            axes[i].set_xlim((380,720))
            axes[i].set_ylim((-3,15))
            axes[i].set_title(abc + nome, fontsize=20)
            axes[i].xaxis.set_tick_params(labelsize=12)
            axes[i].yaxis.set_tick_params(labelsize=12)
            if i == (1,0):
                #axes[i].set_ylabel(r'$ \alpha_{t} \ (m^{-1})$', fontsize=35)
                #axes[i].set_ylabel(r'$ \alpha_{CDM} \ (m^{-1})$', fontsize=35)
                axes[i].set_ylabel(r'$ \alpha_{\varphi} \ (m^{-1})$', fontsize=35)
            if i == (2,1):
                #axes[i].set_xlabel('Wavelenght (nm)', fontsize=30)
                axes[i].set_xlabel('Comprimento de Onda (nm)', fontsize=35)
            #at_dict[x].loc[:max].plot(ylim=(0,22),legend=False,title='Total absorption {}'.format(x),ax=axes[i],lw=1)
            #plt.ylabel("at (m-1)")
        plt.subplots_adjust(top=0.946,bottom=0.097,left=0.078,right=0.977,hspace=0.3,wspace=0.2) 
        #fig.delaxes(axes[2,2])    
        #fig.delaxes(axes[2,1])    
        #fig.delaxes(axes[2,0])   
        plt.show()


###############################################################################
############## Gráfico para plotar a dispersão de cada modelo X dado real
############################################################################        


Plt_scater_ats(a_ACS_water,at_Dict,list(r_ref.keys()))
Plt_scater_ats(data_obj.aCDM,adg_dict,list(r_ref.keys()))
Plt_scater_ats(data_obj.aphy,aphy_dict,list(r_ref.keys()))


def Plt_scater_ats(acs,at_dict,wvl):
    '''
    Function to plot a scatter from X,Y data
     Especific for ROGERIO's DATA
    '''
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    
    
    with plt.style.context(('fivethirtyeight')):
        
        fig, axes = plt.subplots(nrows=3, ncols=3)
        a_ACS_corr = acs
        
        xes = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1)]
        abcd = ['a) ','b) ','c) ','d) ','e) ','f) ','g) ', 'h) ']
    
        for ats,i,abc in zip(at_dict,xes,abcd):
            nome = ats.split('_')[2]
            axes[i].set_title(abc + nome, fontsize=20)
            axes[i].xaxis.set_tick_params(labelsize=14)
            axes[i].yaxis.set_tick_params(labelsize=14)
            axes[i].set_xlim((-1,13))
            axes[i].set_ylim((-2.5,13))
            axes[i].arrow(0,0,25,25,color='blue')
            if i == (1,0):
                #axes[i].set_ylabel(r'$ \alpha_t \ estimado \  (m^{-1})$', fontsize=35)
                #axes[i].set_ylabel(r'$ \alpha_{CDM} \ estimado \  (m^{-1})$', fontsize=35)
                axes[i].set_ylabel(r'$ \alpha_{\varphi} \ estimado \  (m^{-1})$', fontsize=35)
            if i == (2,1):
                #axes[i].set_xlabel('Wavelenght (nm)', fontsize=30)
                #axes[i].set_xlabel(r'$ \alpha_t \ medido \ (m^{-1})$', fontsize=35)
                #axes[i].set_xlabel(r'$ \alpha_{CDM}  \ medido \ (m^{-1})$', fontsize=35)
                axes[i].set_xlabel(r'$ \alpha_{\varphi}  \ medido \ (m^{-1})$', fontsize=35)
                
            for wl,cor in zip(wvl[1:11],colors):
                lx = axes[i].scatter(a_ACS_corr.loc[wl],at_dict[ats].loc[wl],c=cor, s=15)
        axes[i].legend(loc='lower right',bbox_to_anchor=(2, -0.35),title='Bandas (nm):', fancybox = True)
                #leg = fig.legend(lx,(wvl[1:11]), loc='lower right', fontsize=20, title=r'$Bandas \  \lambda_0$',) #title_fontsize=20)
                #leg = fig.legend(lx, (413, 443, 490, 510, 560, 620, 665, 674, 681, 709), loc='lower right', fontsize=20, title=r'$Bandas \  \lambda_0$',) #title_fontsize=20)
                #leg.get_title().set_fontsize(20)
        plt.subplots_adjust(top=0.946,bottom=0.097,left=0.078,right=0.977,hspace=0.3,wspace=0.2) 
        fig.delaxes(axes[2,2])
        plt.show()
        
        
        
        
###############################################################################
############# Gráfico para plotar as stats de cada modelo (MAPE,NRMSE)
###############################################################################        
        
def Plt_stats_ats(acs,at_dict,param):
    '''
    Function to plot a scatter from X,Y data
     Especific for ROGERIO's DATA
    '''
    
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    
    
    with plt.style.context(('seaborn')):
    #with plt.style.context(('fivethirtyeight')):
        from qaapy.utils.stats import Stats_QAA_OLCI
        fig, axes = plt.subplots(nrows=3, ncols=3)
        a_ACS_corr = acs
        
        xes = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1)]
        abcd = ['a) ','b) ','c) ','d) ','e) ','f) ','g) ', 'h) ']
    
        for ats,i,abc in zip(at_dict,xes,abcd):
            nome = ats.split('_')[2]
            stat =  Stats_QAA_OLCI(acs,at_dict[ats].loc[:750])
            
            axes[i].plot(stat[param[0]],lw=1, c= '#122918',linestyle= '-')
            axes[i].plot(stat[param[1]]*100,lw=1, c= '#ab0c0c',linestyle= '--')
            
            
            axes[i].set_title(abc + nome, fontsize=20)
            axes[i].xaxis.set_tick_params(labelsize=14)
            axes[i].yaxis.set_tick_params(labelsize=14)
            axes[i].set_xlim((380,720))
            axes[i].set_ylim((0,150))
            
            #axes[i].arrow(0,0,25,25,color='blue')
            if i == (1,0):
                axes[i].set_ylabel(r'$ Mape / NRMSE (\%)$', fontsize=35)
            if i == (2,1):
                #axes[i].set_xlabel('Wavelenght (nm)', fontsize=35)
                axes[i].set_xlabel('Comprimento de Onda (nm)', fontsize=35)
                #axes[i].set_xlabel(r'$ \alpha_t \ medido (m^{-1})$', fontsize=30)
                
                
        leg = axes[i].legend(loc='lower right',bbox_to_anchor=(2, 0.2), fancybox = True, prop={'size': 20})
                #leg = fig.legend(lx,(wvl[1:11]), loc='lower right', fontsize=20, title=r'$Bandas \  \lambda_0$',) #title_fontsize=20)
                #leg = fig.legend(lx, (413, 443, 490, 510, 560, 620, 665, 674, 681, 709), loc='lower right', fontsize=20, title=r'$Bandas \  \lambda_0$',) #title_fontsize=20)
        leg.get_title().set_fontsize(20)
        plt.subplots_adjust(top=0.946,bottom=0.097,left=0.078,right=0.977,hspace=0.3,wspace=0.2) 
        fig.delaxes(axes[2,2])
        plt.show()        
        
        
Plt_stats_ats(data_obj.aCDM,adg_dict,['Mape','NRMSE']) 
Plt_stats_ats(data_obj.aphy,aphy_dict,['Mape','NRMSE']) 
            
        
        
        
        
        

###############################################################################
        ########### Plota as estatisticas relacionadas ao at para cada lambda0
        ########### SCRIPT: for_qaapy_T1
        #######################################################################
def Plt_stats_at(Dict,stat1,stat2,stat3,eq):
    from matplotlib import pyplot as plt
    with plt.style.context(('fivethirtyeight')):
        fig, axes = plt.subplots(nrows=2, ncols=2)
        #fig.title('2017 - Lago Curuai - Lb0 {} - eq {} '.format(r_ref['lb0'],lago))
        fig.suptitle('Equação {}'.format(eq)) # or plt.suptitle('Main title')
        
        #for nome in list(Dict.index):
        #plt.subplot(2, 2, 1)
        wl = '$ Comprimento \ de \  Onda \ (nm) $'
        l1 = axes[0,0].plot(Dict[stat1],lw=1)
        axes[0,0].set_title('a) ' + Dict[stat1].name)
        axes[0,0].grid(True)
        #axes[0,0].set_xlim((380,800))
        axes[0,0].set_ylim((0,100))
        axes[0,0].set_ylabel('(%)', fontsize=20)
        axes[0,0].set_xlabel( wl, fontsize=20)
        
        
        axes[0,1].plot(Dict[stat2],lw=1)
        axes[0,1].set_title('b) ' + Dict[stat2].name)
        axes[0,1].grid(True)
        axes[0,1].set_ylim((0,1))
        #axes[0,1].set_ylabel('at ', fontsize=30)
        axes[0,1].set_xlabel(wl, fontsize=20)
        
        axes[1,0].plot(Dict[stat3],lw=1)
        axes[1,0].set_title('c) ' + Dict[stat3].name)
        axes[1,0].grid(True)
        axes[1,0].set_ylim((0,1))
        axes[1,0].set_xlabel(wl, fontsize=20)

        leg = fig.legend(l1,('560 nm', '674 nm', '709 nm ', '754 nm'), loc='lower right', fontsize=20, title=r'$Bandas \  \lambda_0$',) #title_fontsize=20)
        leg.get_title().set_fontsize(20)
        fig.delaxes(axes[1,1])    
        plt.subplots_adjust(top=0.9,bottom=0.07,left=0.07,right=0.95,hspace=0.4,wspace=0.2)        

    plt.show()    




















