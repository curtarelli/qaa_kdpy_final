import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def Plt_QAA(X,Y,wvl,ano,lb0,eq,R2 = 'N'):
    '''
    Function to plot a scatter from X,Y data
     Especific for ROGERIO's DATA
    '''
    at_corr = X.loc[:750,:]
    fig, ax = plt.subplots()
    a_ACS_corr = Y
    i = 0
    for wl in wvl[1:11]:
        ax.scatter(a_ACS_corr.loc[wl],at_corr.loc[wl])
        ax.legend()
        if R2 == 'Y':
            df_temp = pd.DataFrame({'predito': at_corr.loc[wl], 'medido': a_ACS_corr.loc[wl]})
            df_corr_r2 = (df_temp.corr().iloc[1].iloc[0]) ** 2
            plt.text(10 , 24 - i, '%d R2 = %0.2f' % (wl,df_corr_r2))
            i += 0.7
        
    plt.xlim(0, 25)
    plt.ylim(0, 25)
    plt.title('{} - Lago Curuai - Lb0 {} - eq {} '.format(ano,lb0,eq))
    plt.xlabel("at ACS + water (<0.7m)")
    plt.ylabel("at QAA (m-1)")
    plt.arrow(0,0,25,25,color='blue')
    plt.show()
    ### R2
    
    ########### UNDER CONSTRUCTION ########################
    
def Plt_stats(Dict,stat1,stat2,stat3,ano,lb0,eq):
    with plt.style.context(('seaborn')):
        fig = plt.figure()
        #fig.title('2017 - Lago Curuai - Lb0 {} - eq {} '.format(r_ref['lb0'],lago))
        fig.suptitle('{} - Lago Curuai - {}nm - eq {} '.format(ano,lb0,eq)) # or plt.suptitle('Main title')
        
        #for nome in list(Dict.index):
        plt.subplot(2, 2, 1)
        plt.plot(Dict[stat1])
        plt.legend()
        plt.title(Dict[stat1].name)
        plt.xlabel('Wavelength')
        plt.ylim(0, 100)
        plt.grid(True)
        
        plt.subplot(2, 2, 2)
        plt.plot(Dict[stat2])
        plt.legend()
        plt.title(Dict[stat2].name)
        plt.xlabel('Wavelength')
        plt.ylim(0, 1)
        plt.grid(True)
        
        plt.subplot(2, 2, 3)
        plt.plot(Dict[stat3])
        plt.legend()
        plt.title(Dict[stat3].name)
        plt.xlabel('Wavelength')
        plt.ylim(0, 1)
        plt.grid(True)
    plt.show()    
    
def Plt_merged_stats(Dict,stat1,stat2,stat3,ano,lb0,eq):
    
    with plt.style.context(('default')):
        fig = plt.figure()
        #fig.title('2017 - Lago Curuai - Lb0 {} - eq {} '.format(r_ref['lb0'],lago))
        fig.suptitle('{} - Lago Curuai - {}nm - eq {} '.format(ano,lb0,eq)) # or plt.suptitle('Main title')
        
        ##
        plt.subplot(2, 2, 1)
        plt.plot(Dict[stat1])
        plt.legend(Dict[stat1].columns)
        plt.title(stat1)
        plt.xlabel('Wavelength')
        plt.grid(True)
        ##
        plt.subplot(2, 2, 2)
        plt.plot(Dict[stat2])
        plt.legend(Dict[stat2].columns)
        plt.title(stat2)
        plt.xlabel('Wavelength')
        plt.grid(True)
        ##
        plt.subplot(2, 2, 3)
        plt.plot(Dict[stat3])
        plt.legend(Dict[stat3].columns)
        plt.title(stat3)
        plt.xlabel('Wavelength')
        plt.grid(True)
    plt.show()
    
def Plt_merged_lb0(Dict,stat1,stat2,stat3):
    
    with plt.style.context(('seaborn')):
        fig = plt.figure()
        #fig.title('2017 - Lago Curuai - Lb0 {} - eq {} '.format(r_ref['lb0'],lago))
        #fig.suptitle('{} - Lago Curuai - {}nm - eq {} '.format(ano,lb0,eq)) # or plt.suptitle('Main title')
        
        ##
        plt.subplot(2, 2, 1)
        plt.plot(Dict[stat1])
        plt.legend(Dict[stat1].columns)
        plt.title(stat1)
        plt.xlabel('Wavelength')
        plt.grid(True)
        ##
        plt.subplot(2, 2, 2)
        plt.plot(Dict[stat2])
        plt.legend(Dict[stat2].columns)
        plt.title(stat2)
        plt.xlabel('Wavelength')
        plt.grid(True)
        ##
        plt.subplot(2, 2, 3)
        plt.plot(Dict[stat3])
        plt.legend(Dict[stat3].columns)
        plt.title(stat3)
        plt.xlabel('Wavelength')
        plt.grid(True)
    plt.show() 

def Plt_atNacs(at,acs,ano,lb0,eq):
    '''
    Function to plot lines from X,Y data
     Especific for ROGERIO's DATA
    '''
    acs.index.name = 'Wavelength'
    with plt.style.context(('seaborn')):
        fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
        fig.suptitle('{} - Lago Curuai - Lb0 {} - eq {} '.format(ano,lb0,eq))
        at.plot(ylim=(0,22),legend=False,title='Total absorption QAA',ax=axes[0])
        plt.ylabel("at (m-1)")
        acs.plot(ylim=(0,22),legend=False,title='Total absorption ACS',ax=axes[1])
        plt.ylabel("at (m-1)")
        plt.show()


def Plt_atALL(at_dict,acs,max):
    '''
    Function to plot lines from X,Y data
     Especific for ROGERIO's DATA
    '''
    acs.index.name = 'Wavelength'
    with plt.style.context(('seaborn')):
        fig, axes = plt.subplots(nrows=4, ncols=2) 
        #fig.suptitle('{} - Lago Curuai - Lb0 {} - eq {} '.format(ano,lb0,eq))
        ###
        acs.plot(ylim=(0,22),legend=False,title='Total absorption ACS',ax=axes[0,0])
        plt.ylabel("at (m-1)")
        
        xes = [(1,1),(1,0),(0,1),(2,0),(2,1),(3,0),(3,1)]
        #xes = [(0,1),(1,0),(1,1)]
        for x,i in zip(at_dict,xes):

            at_dict[x].loc[:max].plot(ylim=(0,22),legend=False,title='Total absorption {}'.format(x),ax=axes[i])
            plt.ylabel("at (m-1)")
            
        plt.show()
#            ###
#            at.plot(ylim=(0,22),legend=False,title='Total absorption QAA',ax=axes[0,1])
#            plt.ylabel("at (m-1)")
#            ###
#            at.plot(ylim=(0,22),legend=False,title='Total absorption QAA',ax=axes[0,2])
#            plt.ylabel("at (m-1)")
#            ###
#            at.plot(ylim=(0,22),legend=False,title='Total absorption QAA',ax=axes[1,0])
#            plt.ylabel("at (m-1)")
#            ###
#            at.plot(ylim=(0,22),legend=False,title='Total absorption QAA',ax=axes[1,1])
#            plt.ylabel("at (m-1)")        
#            ###




    