'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

RESUMO: Este script contém as funções necessárias para a realização dos testes
de ajustes dos passos empíricos do modelo QAA (Passos 2 e 4). O teste conduzido
é baseado na metodologia apresentada no seguinte artigo:
    
X. Liu et al. (2019) - "Remote Sensing of Secchi Depth in Highly Turbid Lake
Waters and Its Application with MERIS Data"
'''
##  Importando pacotes básicos para o código.
import numpy as np

###############################################################################
def calc_aRrs_ratio(R_ref, R_ref1, R_ref2, R_ref3):
    '''
    Esta função aplica a razão de três bandas proposta para o Passo 2 do modelo
    QAA. Esta função usa como valores de entrada a Rrs medida de campo e três 
    comprimentos de onda de referência de escolha do usuário.
    
    [ Razão = Rrs(ref1) / ( Rrs(ref2) + Rrs(ref3) ) ]
    
    ----------
    Parameters
    ----------
    R_ref [Dictionary]
        Reflectância de sensoriamento remoto superficial (Rrs) filtrada para 
        as bandas de referência.
    R_ref1 [Value]
        Comprimento de onda de referência 1 para Rrs, por padrão se usa o lb0
        (Comprimento de onda incial usado no modelo QAA).
    R_ref2 [Value]
        Comprimento de onda de referência 2 pra Rrs.
    R_ref3 [Value]
        Comprimento de onda de referência 3 pra Rrs.

    -------
    Returns
    -------
    Rrs_ratio [Series]
        Cálculo da razão de três bandas do modelo de absorção em lb0 (QAA - Passo 2).
    '''
    return R_ref[R_ref1] / (R_ref[R_ref2] + R_ref[R_ref3])

###############################################################################
def calc_bb_ref(u, wvl_ref, a):
    '''
    Esta função realiza o cálculo do retro espalhamento total na coluna
    d'água (bbp) com base nos dados de absorção total na coluna d'água
    
    [Ex.: absorção ACS].    
    
    ----------
    Parameters
    ----------
    u [Data Frame]
        Parâmetro u (QAA - Passo 1) calculado para cada estação amostral e por
        comprimento de onda.
    wvl_ref [Series]
        Faixa de comprimentos de onda utilizados nos cálculos.
    a [Data Frame]
        Absorção total na coluna d'água para cada estação amostral por comprimento de onda.

    -------
    Returns
    -------
    bbp [Data Frame]
        Retroespalhamento do particulado presente na coluna d'água com base na
        absorção total na coluna d'água.
    '''    
    return (( u.loc[wvl_ref] * a.loc[wvl_ref] ) / (1 - u.loc[wvl_ref]))

###############################################################################
def calc_bbp_ref(u, bbw, wvl_ref, a):
    '''
    Esta função realiza o cálculo do retro espalhamento do particulado na coluna
    d'água (bbp) com base nos dados de absorção total na coluna d'água
    
    [Ex.: absorção ACS].    
    
    ----------
    Parameters
    ----------
    u [Data Frame]
        Parâmetro u (QAA - Passo 1) calculado para cada estação amostral e por
        comprimento de onda.
    bbw [Data Frame]
        Retroespalhamento da água por comprimento de onda.
    wvl_ref [Series]
        Faixa de comprimentos de onda utilizados nos cálculos.
    a [Data Frame]
        Absorção total na coluna d'água para cada estação amostral por comprimento de onda.

    -------
    Returns
    -------
    bbp [Data Frame]
        Retroespalhamento do particulado presente na coluna d'água com base na
        absorção total na coluna d'água.
    '''    
    return (( u.loc[wvl_ref] * a.loc[wvl_ref] ) / (1 - u.loc[wvl_ref])) - bbw.loc[wvl_ref].values

###############################################################################
def calc_bbp_lin(wvl_ref, bbp, lb0):
    '''
    Esta função realiza os cálculos necessários para a regressão linear usada
    para estimar o decaimento (slope) do retroespalhamento do particulado presente
    na coluna d'água.
    
    ----------
    Parameters
    ----------
    wvl_ref [Series]
        Faixa de comprimentos de onda utilizados nos cálculos.
    bbp [Data Frame]
        Retroespalhamento do particulado presente na coluna d'água.
    lb0 [Value]
        Comprimento de onda incial usado no modelo QAA (Lambda Zero).

    -------
    Returns
    -------
    bbp_ln [Dictionary]
        Dicionário contendo um Data Frame do logarítmo natural de bbp espectral
        por estação amostral, uma Série contendo o logarítmo natural de bbp em
        lb0 e uma Série do logarítmo natural da razão entre lb0 e comprimento
        de onda.
    '''
    ratio_wl = np.log(float(lb0) / wvl_ref)
    ratio_wl = ratio_wl.rename(wvl_ref)
    ratio_wl = ratio_wl.rename(lb0, axis = 'columns')
    
    return {'ln_bbp': np.log(bbp),
            'ln_bbp_lb0': np.log(bbp.loc[lb0]),
            'ln_ratio': ratio_wl}

###############################################################################
def calc_nrrs_ratio(r_ref, r_ref1, r_ref2):
    '''
    Esta função aplica a razão de duas bandas proposta para o Passo 4 do modelo
    QAA. Esta função usa como valores de entrada a rrs medida de campo e dois 
    comprimentos de onda de referência de escolha do usuário.
    
    [ Razão = exp( rrs(ref1) / rrs(ref2) ) ]
    
    ----------
    Parameters
    ----------
    r_ref [Dictionary]
        Reflectância de sensoriamento remoto subsuperficial (rrs) filtrada para 
        as bandas de referência.
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.

    -------
    Returns
    -------
    rrs_ratio [Series]
        Cálculo da razão de duas bandas do modelo do slope de decaimento do bbp.
    '''
    return np.exp(r_ref[r_ref1] / r_ref[r_ref2])    

###############################################################################
def calc_nrrs_ratio_neg(r_ref, r_ref1, r_ref2):
    '''
    Esta função aplica a razão negativa de duas bandas proposta para o Passo 4
    do modelo QAA. Esta função usa como valores de entrada a rrs medida de campo
    e dois comprimentos de onda de referência de escolha do usuário.
    
    [ Razão = exp( - rrs(ref1) / rrs(ref2) ) ]
    
    ----------
    Parameters
    ----------
    r_ref [Dictionary]
        Reflectância de sensoriamento remoto subsuperficial (rrs) filtrada para 
        as bandas de referência.
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.
        
    -------
    Returns
    -------
    rrs_ratio [Series]
        Cálculo da razão de duas bandas do modelo do slope de decaimento do bbp.
    '''
    return np.exp(- r_ref[r_ref1] / r_ref[r_ref2])
