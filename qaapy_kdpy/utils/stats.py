'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é um pacote de ferramentas uteis para calculo de parâmetros estatísticos
com base nos dados de entrada e/ou resultados dos algorítmos QAA e de estimativa de Kd.
'''
##  Importanto pacotes básicos necessários para as funções.
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error as mse
from matplotlib import pyplot as plt

###############################################################################
########################## Base statistic #####################################
###############################################################################
def mean_abs_perc_error(y_true, y_pred):
    '''
    Função para o calculo do Erro Absoluto Médio Percentual (Mean Absolute 
    Percentual Error - MAPE) com base em series de conjuntos amostrais de
    referência e simulados fornecidos pelo usuário.
    
    ----------
    Parameters
    ----------
    y_true [Series]
        Conjunto amostral de referência, "verdade" amostral.
        
        Ex.: Valores medidos de campo, valores de análise laboratorial.
        
    y_pred [Series]
        Conjunto amostral predito, valor estimado.

        Ex.: Valores simulados por modelo com uso de dados e/ou imagem.
    
    -------
    Returns
    -------
    MAPE [Value]
        Valor calculado para o Erro Absoluto Médio Percentual com base nos
        valores de referência e preditos fornecidos pelo usuário.
    '''
    y_true, y_pred = np.array(y_true), np.array(y_pred)
    return np.mean(np.abs((y_true - y_pred) / y_true)) * 100

###############################################################################
def correlation_coef(y_true,y_pred):
    '''
    Função para o calculo do Coeficiente de Correlação (r) com base em series
    de conjuntos amostrais de referência e simulados fornecidos pelo usuário.
    
    ----------
    Parameters
    ----------
    y_true [Series]
        Conjunto amostral de referência, "verdade" amostral.
        
        Ex.: Valores medidos de campo, valores de análise laboratorial.
        
    y_pred [Series]
        Conjunto amostral predito, valor estimado.

        Ex.: Valores simulados por modelo com uso de dados e/ou imagem.
    
    -------
    Returns
    -------
    r [Value]
        Valor calculado para o Coeficiente de Correlação (r) com base nos
        valores de referência e preditos fornecidos pelo usuário.
    '''
    df_temp = pd.DataFrame({'Real': y_true, 'predicted': y_pred})
    return (df_temp.corr().iloc[1].iloc[0])

###############################################################################
def determination_coef(y_true,y_pred):
    '''
    Função para o calculo do Coeficiente de Determinação (R²) com base em series
    de conjuntos amostrais de referência e simulados fornecidos pelo usuário.
    
    ----------
    Parameters
    ----------
    y_true [Series]
        Conjunto amostral de referência, "verdade" amostral.
        
        Ex.: Valores medidos de campo, valores de análise laboratorial.
        
    y_pred [Series]
        Conjunto amostral predito, valor estimado.

        Ex.: Valores simulados por modelo com uso de dados e/ou imagem.
    
    -------
    Returns
    -------
    R² [Value]
        Valor calculado para o Coeficiente de Determinação (R²) com base nos
        valores de referência e preditos fornecidos pelo usuário.
    '''
    df_temp = pd.DataFrame({'Real': y_true, 'predicted': y_pred})
    return (df_temp.corr().iloc[1].iloc[0]) ** 2

###############################################################################
def bias_error(y_true,y_pred):
    '''
    Função para o calculo do Enviesamento (Bias) com base em series de conjuntos
    amostrais de referência e simulados fornecidos pelo usuário.
    
    ----------
    Parameters
    ----------
    y_true [Series]
        Conjunto amostral de referência, "verdade" amostral.
        
        Ex.: Valores medidos de campo, valores de análise laboratorial.
        
    y_pred [Series]
        Conjunto amostral predito, valor estimado.

        Ex.: Valores simulados por modelo com uso de dados e/ou imagem.
    
    -------
    Returns
    -------
    Bias [Value]
        Valor calculado para o Enviesamento (Bias) com base nos valores de
        referência e preditos fornecidos pelo usuário.
    '''
    return np.mean(y_pred - y_true)
    
def root_mean_sqrt_error(y_true,y_pred):
    '''
    Função para o calculo da Raiz do Erro Quadrático Médio (Root-Mean-Square
    Error - RMSE) com base em series de conjuntos amostrais de referência e
    simulados fornecidos pelo usuário.
    
    ----------
    Parameters
    ----------
    y_true [Series]
        Conjunto amostral de referência, "verdade" amostral.
        
        Ex.: Valores medidos de campo, valores de análise laboratorial.
        
    y_pred [Series]
        Conjunto amostral predito, valor estimado.

        Ex.: Valores simulados por modelo com uso de dados e/ou imagem.
    
    -------
    Returns
    -------
    RMSE [Value]
        Valor calculado para a Raiz do Erro Quadrático Médio (Root-Mean-Square
        Error - RMSE) com base nos valores de referência e preditos fornecidos
        pelo usuário.
    '''
    return np.sqrt(mse(y_true, y_pred))

###############################################################################
########################    TESTES INCOMPLETOS   ##############################
###############################################################################
#def Stats_QAA(Real,Predicted):
#    '''
#    Function to retrieve statistics (MAPE,R2,RMSE, NRMSE)
#    ....
#    '''
#    def MAPE(y_true, y_pred):
#        from sklearn.metrics import mean_squared_error as mse
#        y_true, y_pred = np.array(y_true), np.array(y_pred)
#        return np.mean(np.abs((y_true - y_pred) / y_true)) * 100
#    Df_stats = pd.DataFrame()
#    for Wl in range(400,751):
#        Real_C = Real
#        Predicted_C = Predicted.loc[:750,:]
#        try:
#            Df_stats.at[Wl,'Mape'] = MAPE(Real_C.loc[Wl], Predicted_C.loc[Wl])
#            #Df_stats.at[Wl,'Mape'] = (sum(abs((Real_C.loc[Wl] - Predicted_C.loc[Wl]) / Real_C.loc[Wl])) / len(Real_C.loc[Wl]))  * 100
#        except:
#            pass
#        try:
#            df_temp = pd.DataFrame({'predito': Predicted_C.loc[Wl], 'medido': Real_C.loc[Wl]})
#            df_corr_r2 = (df_temp.corr().iloc[1].iloc[0]) ** 2
#            Df_stats.at[Wl,'R2'] =  df_corr_r2
#        except:
#            pass
#        try:
#            Df_stats.at[Wl,'MSE'] = mse(Real_C.loc[Wl],Predicted_C.loc[Wl])
#        except:
#            pass
#        Df_stats.at[Wl,'RMSE'] = np.sqrt(Df_stats.loc[Wl,'MSE'])
#        Df_stats.at[Wl,'NRMSE'] = Df_stats.loc[Wl,'RMSE'] / (Real_C.loc[Wl].max() - Real_C.loc[Wl].min())
#    return Df_stats

###############################################################################
def Stats_QAA_OLCI(Real,Predicted):
    '''
    Function to retrieve statistics (MAPE,R2,RMSE, NRMSE,Bias)
    '''
    from sklearn.metrics import mean_squared_error
    Df_stats = pd.DataFrame()
    for Wl in list(Real.index):
        _real_input = Real.loc[Wl]
        _predicted_Input = Predicted.loc[Wl]
        try:
            Df_stats.at[Wl,'Mape'] = mean_abs_perc_error(_real_input, _predicted_Input)
            #Df_stats.at[Wl,'Mape'] = (sum(abs((Real_C.loc[Wl] - Predicted_C.loc[Wl]) / Real_C.loc[Wl])) / len(Real_C.loc[Wl]))  * 100
        except:
            print('erro Mape')
            #pass
        try:
            Df_stats.at[Wl,'R2'] =  correlation_coef(_real_input,_predicted_Input)
        except:
            print('erro R2')
            #pass
        try:
            Df_stats.at[Wl,'MSE'] = mean_squared_error(_real_input, _predicted_Input)
        except:
            print('erro MSE')
            #pass
        try:
            Df_stats.at[Wl,'RMSE'] = np.sqrt(mean_squared_error(_real_input, _predicted_Input))
            Df_stats.at[Wl,'NRMSE'] = np.sqrt(mean_squared_error(_real_input, _predicted_Input)) / (_real_input.max() - _real_input.min())
        except:
            print('erro RMSE')
            #pass
        try:
            Df_stats.at[Wl,'Bias'] = bias_error(_real_input, _predicted_Input)
        except:
            print('erro Bias')
            #pass
        
    return Df_stats

###############################################################################
def Stats_QAA_Chl_OLCI(Real,Predicted,name='data'):
    '''
    Function to retrieve statistics (MAPE,R2,RMSE, NRMSE,Bias)
    '''
    import pandas as pd
    from sklearn.metrics import mean_squared_error
    
    Df_stats = pd.DataFrame()
    _real_input = Real
    _predicted_Input = Predicted
    
    try:
        Df_stats.at[name,'Mape'] = mean_abs_perc_error(_real_input, _predicted_Input)
        #Df_stats.at[Wl,'Mape'] = (sum(abs((Real_C.loc[Wl] - Predicted_C.loc[Wl]) / Real_C.loc[Wl])) / len(Real_C.loc[Wl]))  * 100
    except:
        print('erro Mape')
        #pass
    try:
        Df_stats.at[name,'R2'] =  correlation_coef(_real_input,_predicted_Input)
    except:
        print('erro R2')
        #pass
    try:
        Df_stats.at[name,'MSE'] = mean_squared_error(_real_input, _predicted_Input)
    except:
        print('erro MSE')
        #pass
    try:
        Df_stats.at[name,'RMSE'] = np.sqrt(mean_squared_error(_real_input, _predicted_Input))
        Df_stats.at[name,'NRMSE'] = np.sqrt(mean_squared_error(_real_input, _predicted_Input)) / (_real_input.max() - _real_input.min())
    except:
        print('erro RMSE')
        #pass
    try:
        Df_stats.at[name,'Bias'] = bias_error(_real_input, _predicted_Input)
    except:
        print('erro Bias')
        #pass
    return Df_stats

###############################################################################
def merge_opt_stats(D1_N,D2_Y):
    Dict = {}
    for stats in D1_N.columns:
        statNY =pd.DataFrame()
        statNY[ stats + '_optN'] = D1_N[stats]
        statNY[ stats + '_optY'] = D2_Y[stats]
        Dict[stats] = statNY
    return Dict

###############################################################################
def merge_lb0_stats(D1,lb1,D2,lb2):
    Dict = {}
    for stats in D1.columns:
        statNY =pd.DataFrame()
        statNY[ stats + '_' + str(lb1)] = D1[stats]
        statNY[ stats + '_' + str(lb2)] = D2[stats]
        Dict[stats] = statNY
    return Dict    