###############################################################################
#############           QAA - Monte Carlo - lee 2013 - Main App          #################
###############################################################################
######################## IMPORTS ##############################################
###############################################################################

import pandas as pd
import numpy as np


def ramdom_bands(wvl_ref,r_ref,lb0,n):
    
    import numpy as np
    import random
    s = random.sample(wvl_ref['referencia'], n)
    X = random.randint(1, 5)
    #eq_stdd = np.log10( (r_ref[s[0]] +r_ref[s[1]] ) / ( r_ref['r_lb0'] + ( X * ((r_ref[s[2]] ** 2) / r_ref[s[3]] )) )) 
    eq_stdd=[]
    return eq_stdd, s, X



def Stats_MC_QAA(Real,Predicted,nm):
    import pandas as pd
    '''
    Function to retrieve statistics (MAPE,R2,RMSE, NRMSE)
    ....
    '''
    Df_stats = pd.DataFrame()
    for i in np.arange(nm):
            Real_C = Real
            Predicted_C = Predicted[:351,:,i]
            try:
                Df_stats[('Mape_' + str(i))] = (np.mean((np.abs(Real_C.sub(Predicted_C).div(Real_C))),axis=1)) * 100
            except:
                pass
#            try:
#                Df_stats['R2_' + str(i)] = (df_temp.corr().iloc[1].iloc[0]) ** 2
#                Df_stats['R2_' + str(i)] = np.corrcoef(Real_C, Predicted_C,  rowvar=True)[0,1] ** 2
#            except:
#                pass
    return Df_stats


def Stats_MC_QAA_old(Real,Predicted,nm):
    import pandas as pd
    '''
    Function to retrieve statistics (MAPE,R2,RMSE, NRMSE)
    ....
    '''
    def MAPE(y_true, y_pred):
        y_true, y_pred = np.array(y_true), np.array(y_pred)
        return np.mean(np.abs((y_true - y_pred) / y_true)) * 100
    Df_stats = pd.DataFrame()
    for i in np.arange(nm):
        for Wl,wx in zip(range(400,751),range(351)): 
            Real_C = Real
            Predicted_C = Predicted[:351,:,i]
            try:
                Df_stats.at[Wl,('Mape_' + str(i))] = MAPE(Real_C.loc[Wl], Predicted_C[wx,:])
                #Df_stats.at[Wl,'Mape'] = (sum(abs((Real_C.loc[Wl] - Predicted_C.loc[Wl]) / Real_C.loc[Wl])) / len(Real_C.loc[Wl]))  * 100
            except:
                pass
            try:
                df_temp = pd.DataFrame({'predito': Predicted_C[wx,:], 'medido': Real_C.loc[Wl]})
                Df_stats.at[Wl,'R2_' + str(i)] = (df_temp.corr().iloc[1].iloc[0]) ** 2
            except:
                pass
    return Df_stats