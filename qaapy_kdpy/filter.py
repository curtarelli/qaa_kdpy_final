'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Funções para filtrar e padronizar os dados do QAA para certas operações.
'''
##  Importando pacotes básicos necessários para este script.
import pandas as pd
import numpy as np

###############################################################################
def filter_wvl(rrs, wvl_ref):
    '''
    Função que filtra os dados de reflectância de sensoriamento remoto subsuperficial
    com base em lista/serie de bandas passada pelo usuário.
    
    ----------
    Parameters
    ----------
    rrs [Data Frame]
        Dados de reflectância de sensoriamento remoto subsuperficial (rrs).
    wvl_ref [Series]
        Série contendo bandas a serem filtradas e padronizadas nos dados.

    -------
    Returns
    -------
    rrs_ref [Data Frame]
        Dados de reflectância de sensoriamento remoro subsuperficial fitlrados
        para os comprimentos de onda selecionados pelo usuário.
    '''
    rrs_ref = {}
    for wlr in wvl_ref['referencia']:
        rrs_ref[wlr] = rrs.loc[wlr]
    rrs_ref['r_lb0'] = rrs.loc[wvl_ref['lb0']]

    return rrs_ref

###############################################################################
def ACS_add_water(acs, aw):
    '''
    Função para adicionar o valor do coeficiente de absorção espectral da água
    ao dado de absorção medido pelo ACS.

    ----------
    Parameters
    ----------
    acs [Data Frame]
        Dados de coeficiente de absorção espectral medidos pelo ACS.
    aw [Series]
        Dados de coeficiente de absorção.

    -------
    Returns
    -------
    df [Data Frame]
        Coeficientes de absorção espectral total na coluna d'água.
    '''
    aw_cut = aw.loc[:750,:][0]
    df = pd.DataFrame()
    for i in acs.columns:
        df[i] = acs[i] + aw_cut
    return df    

###############################################################################
##  TODO: Transformar em um filtro geral!
def filter_aw(aw, wvl, filter_in):
    '''
    Filtro para padronizar dados de absorção espectral da água.
    
    ----------
    Parameters
    ----------
    aw [Data Frame]
        Coeficiente de absorção espectral da água.
    wvl [Series / List]
        Comprimentos de onda.
    filter_in [Series / List]
        Comprimentos de onda de referência para o filtro.

    -------
    Returns
    -------
    aw_filtered [Data Frame]
        Coeficientes de absorção espectral da água filtrados.
    '''
    return aw.loc[708][0]
    
    return aw.iloc[(wvl[wvl == filter_in].index[0])][0]

###############################################################################
##  TODO: Transformar em um filtro geral!
def filter_a(a_data, wvl, filter_in, aw):
    '''
    Filtro para padronizar dados de absorção do ACS e adicionar absorção da 
    água (aw Pope and Fry 1997).
    
    ----------
    Parameters
    ----------
    a_data [Data Frame]
        Coeficiente de absorção ACS a serem filtrados.
    wvl [Series / List]
        Comprimentos de onda.
    filter_in [Series / List]
        Comprimentos de onda de referência para o filtro.
    aw [Series]
        Coeficiente de abosorçãoespectral da água.

    -------
    Returns
    -------
    a_filtered [Data Frame]
        Coeficientes de absorção total (ACS + aw) filtrados.
    '''
    return a_data.iloc[(wvl[wvl == filter_in].index[0]), :] + aw

###############################################################################
##  TODO: Transformar em um filtro geral!
def filter_at(at, wvl, filter_in):
    '''
    Filtro para padronizar dados de absorção total.
    
    ----------
    Parameters
    ----------
    at [Data Frame]
        Coeficiente de absorção total.
    wvl [Series / List]
        Comprimentos de onda.
    filter_in [Series / List]
        Comprimentos de onda de referência para o filtro.

    -------
    Returns
    -------
    at_filtered [Data Frame]
        Coeficientes de absorção total filtrados.
    '''
    return at.iloc[(wvl[wvl == filter_in].index[0])]