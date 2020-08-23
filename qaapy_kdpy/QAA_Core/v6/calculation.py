'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é um pacote de funções "exclusivas" para uso do modelo QAA v6.
'''
###############################################################################
def calc_alb0_v6(wvl_ref_lb0, r_ref, r_ref1, r_ref2, aw):
    '''
    Função para calcular A no comprimento de onda inicial (Lb0).
    
    ----------
    Parameters
    ----------
    wvl_ref_lb0 [Value]
        comprimento de onda de referência utilizado como Lambda 0.
    r_ref [Data Frame]
        Reflectância de sensoriamento remoto subsuperficial (rrs) filtrada em
        bandas de referência.
    r_ref1 : TYPE
        Comprimento de onda de referência 1 para rrs.
    r_ref2 : TYPE
        Comprimento de onda de referência 2 pra rrs.
    aw : TYPE
        Coeficiente de absorção da água.

    -------
    Returns
    -------
    a_lb0 [Series]
        Coeficiente de absorção total para o comprimento de onda lb0.
    '''
    return aw.loc[wvl_ref_lb0].values + 0.39 * (r_ref[wvl_ref_lb0] / (r_ref[r_ref1] + r_ref[r_ref2])) ** 1.14

###############################################################################
def calc_alb0_v6_492(wvl_ref_lb0, r_ref, r_ref1, r_ref2, aw):
    '''
    Função para calcular o coeficiente de absorção no através do ajuste da função
    a(lb0) para o comprimento de onda inicial (Lb0) de 492 nm com base em bandas
    do sensor Snetinel 2A - MSI.
    
    ----------
    Parameters
    ----------
    wvl_ref_lb0 [Value]
        comprimento de onda de referência utilizado como Lambda 0.
    r_ref [Data Frame]
        Reflectância de sensoriamento remoto subsuperficial (rrs) filtrada em
        bandas de referência.
    r_ref1 : TYPE
        Comprimento de onda de referência 1 para rrs.
    r_ref2 : TYPE
        Comprimento de onda de referência 2 pra rrs.
    aw : TYPE
        Coeficiente de absorção da água.

    -------
    Returns
    -------
    a_lb0 [Series]
        Coeficiente de absorção total para o comprimento de onda lb0.
    '''
    return aw.loc[wvl_ref_lb0].values + 0.4639 * (r_ref[wvl_ref_lb0] / (r_ref[r_ref1] + r_ref[r_ref2])) ** -1.3466

###############################################################################
def calc_alb0_v6_560(wvl_ref_lb0, r_ref, r_ref1, r_ref2, aw):
    '''
    Função para calcular o coeficiente de absorção no através do ajuste da função
    a(lb0) para o comprimento de onda inicial (Lb0) de 560 nm com base em bandas
    do sensor Snetinel 2A - MSI.
    
    ----------
    Parameters
    ----------
    wvl_ref_lb0 [Value]
        comprimento de onda de referência utilizado como Lambda 0.
    r_ref [Data Frame]
        Reflectância de sensoriamento remoto subsuperficial (rrs) filtrada em
        bandas de referência.
    r_ref1 : TYPE
        Comprimento de onda de referência 1 para rrs.
    r_ref2 : TYPE
        Comprimento de onda de referência 2 pra rrs.
    aw : TYPE
        Coeficiente de absorção da água.

    -------
    Returns
    -------
    a_lb0 [Series]
        Coeficiente de absorção total para o comprimento de onda lb0.
    '''
    return aw.loc[wvl_ref_lb0].values + 0.4310 * (r_ref[wvl_ref_lb0] / (r_ref[r_ref1] + r_ref[r_ref2])) ** -1.4408

###############################################################################
def calc_alb0_v6_aw(wvl_ref_lb0, r_ref, r_ref1, r_ref2, aw):
    '''
    Função para calcular Aw no comprimento de onda inicial (Lb0).
    
    ----------
    Parameters
    ----------
    wvl_ref_lb0 [Value]
        comprimento de onda de referência utilizado como Lambda 0.
    r_ref [Data Frame]
        Reflectância de sensoriamento remoto subsuperficial (rrs) filtrada em
        bandas de referência.
    r_ref1 : TYPE
        Comprimento de onda de referência 1 para rrs.
    r_ref2 : TYPE
        Comprimento de onda de referência 2 pra rrs.
    aw : TYPE
        Coeficiente de absorção da água.

    -------
    Returns
    -------
    aw_lb0 [Series]
        Coeficiente de absorção total para o comprimento de onda lb0.
    '''
    return aw.loc[wvl_ref_lb0].values

###############################################################################
def calc_bbp_lb0_v6(u, bbw, lb0, aref):
    '''
    Função para calcular o retroespalhamento do particulado na coluna d'água no
    comprimento de onda inicial usado no modelo QAA (Lambda Zero)
    
    ----------
    Parameters
    ----------
    u [Data Frame]
        Parâmetro u (QAA - Passo 1) calculado para cada estação amostral e por
        comprimento de onda.
    bbw [Data Frame]
        Retroespalhamento da água por comprimento de onda.
    lb0 [Value]
        Comprimento de onda incial usado no modelo QAA (Lambda Zero).
    aref [Series]
        Absorção referência no comprimento de onda incial (Lambda Zero) para a
        coluna d'água em cada estação amostral.
    
    -------
    Returns
    -------
    bbp_lb0 [Series]
        Retroespalhamento do partículado no comprimentod e onda inicial (Lamda Zero)
        para a coluna d'água em cada estação amostral.
    '''
    return (( u.loc[lb0] * aref ) / (1 - u.loc[lb0])) - bbw.loc[lb0].values
