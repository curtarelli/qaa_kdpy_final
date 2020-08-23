'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é um pacote de funções base usada nos modelos de QAA v5 e outras
versões, para os passos que forem identicos.

Versões com passos distintos usam pacote próprio "exclusivo"
'''
##  Importanto pacotes básicos necessários para as funções.
import numpy as np
import pandas as pd
import math as mt

###############################################################################
############################## FASE 0 e 1 #####################################
###############################################################################
# Função para realizar o cálculo de U
calc_u = lambda g0, g1, rrs: (-g0 + np.sqrt(np.power(g0, 2) + 4 * g1 * rrs)) / (2 * g1)

###############################################################################
def calc_bbp_lb0(u, bbw, lb0, aref):
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
        Retroespalhamento do partículado no comprimento de onda inicial (Lamda Zero)
        para a coluna d'água em cada estação amostral.
    '''
    return ((u.loc[lb0] * aref) / (1 - u.loc[lb0])) - bbw.loc[lb0].values

###############################################################################
# Função para realizar o cálculo de N
def calc_n(r_ref1, r_ref2):
    '''
    Função que aplica a equação para estimativa de decaimento do retroespalhamento
    do material particulado na coluna d'água (slope bbp) com base em duas bandas
    de rrs escolhidas pelo usuário, conforme Passo 4 do modelo QAA.
    
    MODELO = Padrão proposto por Lee et al. 2013.
    
    ----------
    Parameters
    ----------
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.
        
    -------
    Returns
    -------
    n [Series]
        Taxa de decaimento do retroespalhamento do particulado na coluna d'água
        (slope bbp).
    '''
    return 2 * (1 - (1.2 * np.exp(-0.9 * (r_ref1 / r_ref2))))

###############################################################################
def calc_n_acs_492(r_ref1, r_ref2):
    '''
    Função que aplica a equação para estimativa de decaimento do retroespalhamento
    do material particulado na coluna d'água (slope bbp) com base na razão de 2
    bandas rrs usando as bandas 665 e 704 nm do sensor Sentienl 2A - MSI onde foi
    usado o comprimento de onda 492 nm para regressão e ajuste do slope do bbp,
    conforme Passo 4 do modelo QAA.
    
    Este modelo usa o valor de bb estimado através do dado de absorção do ACS.
    
    para o comprimento de onda inicial (Lb0) com base em bandas do sensor
    Sentinel 2A - MSI.
    
    MODELO = Ajustado com dados de ACS + aw (Três Marias 2013)
    
    ----------
    Parameters
    ----------
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.

    -------
    Returns
    -------
    n [Series]
        Taxa de decaimento do retroespalhamento do particulado na coluna d'água
        (slope bbp).
    '''    
    return (0.5248 * np.exp(r_ref1 / r_ref2)) - 1.1849

###############################################################################
def calc_n_acs_560(r_ref1, r_ref2):
    '''
    Função que aplica a equação para estimativa de decaimento do retroespalhamento
    do material particulado na coluna d'água (slope bbp) com base na razão de 2
    bandas rrs usando as bandas 665 e 704 nm do sensor Sentienl 2A - MSI onde foi
    usado o comprimento de onda 560 nm para regressão e ajuste do slope do bbp,
    conforme Passo 4 do modelo QAA.
    
    Este modelo usa o valor de bb estimado através do dado de absorção do ACS.
    
    para o comprimento de onda inicial (Lb0) com base em bandas do sensor
    Sentinel 2A - MSI.
    
    MODELO = Ajustado com dados de ACS + aw (Três Marias 2013)
    
    ----------
    Parameters
    ----------
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.

    -------
    Returns
    -------
    n [Series]
        Taxa de decaimento do retroespalhamento do particulado na coluna d'água
        (slope bbp).
    '''    
    return (0.3991 * np.exp(r_ref1 / r_ref2)) - 0.5983

###############################################################################
def calc_n_hydro_560(r_ref1, r_ref2):
    '''
    Função que aplica a equação para estimativa de decaimento do retroespalhamento
    do material particulado na coluna d'água (slope bbp) com base na razão de 2
    bandas rrs usando as bandas 443 e 492 nm do sensor Sentienl 2A - MSI onde foi
    usado o comprimento de onda 560 nm para regressão e ajuste do slope do bbp,
    conforme Passo 4 do modelo QAA.
    
    Este modelo usa o valor de bb medidos pelo Hydroscat.
    
    para o comprimento de onda inicial (Lb0) com base em bandas do sensor
    Sentinel 2A - MSI.
    
    MODELO = Ajustado com dados de bb do Hydroscat (Três Marias 2013)
    
    ----------
    Parameters
    ----------
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.

    -------
    Returns
    -------
    n [Series]
        Taxa de decaimento do retroespalhamento do particulado na coluna d'água
        (slope bbp).
    ''' 
    return (1.1363 * np.exp(r_ref1 / r_ref2)) - 1.0422

###############################################################################
def calc_n_hydro_704(r_ref1, r_ref2):
    '''
    Função que aplica a equação para estimativa de decaimento do retroespalhamento
    do material particulado na coluna d'água (slope bbp) com base na razão de 2
    bandas rrs usando as bandas 443 e 492 nm do sensor Sentienl 2A - MSI onde foi
    usado o comprimento de onda 704 nm para regressão e ajuste do slope do bbp,
    conforme Passo 4 do modelo QAA.
    
    Este modelo usa o valor de bb medidos pelo Hydroscat.
    
    para o comprimento de onda inicial (Lb0) com base em bandas do sensor
    Sentinel 2A - MSI.
    
    MODELO = Ajustado com dados de bb do Hydroscat (Três Marias 2013)
    
    ----------
    Parameters
    ----------
    r_ref1 [Value]
        Comprimento de onda de referência 1 pra rrs.
    r_ref2 [Value]
        Comprimento de onda de referência 2 pra rrs.

    -------
    Returns
    -------
    n [Series]
        Taxa de decaimento do retroespalhamento do particulado na coluna d'água
        (slope bbp).
    ''' 
    return (1.3065 * np.exp(r_ref1 / r_ref2)) - 1.2805

###############################################################################
def calc_bbp(wvl, n, bbp_lb0, lb0):
    '''
    Função para calcular o retroespalhamento do material particulado na coluna
    d'água (bbp) para todos os comprimentos de onda.

    ----------
    Parameters
    ----------
    wvl [Series / List]
        Comprimentos de onda a serem utilizados.
    n [Series / List]
        Taxa de decaimento do bbp (expoente da equação).
    bbp_lb0 [Series]
        Retroespalhamento do partículado no comprimento de onda inicial (Lamda Zero)
        para a coluna d'água em cada estação amostral.
    lb0 [Value]
        Comprimento de onda inicial (Lambda Zero).

    -------
    Returns
    -------
    bbp [Data Frame]
        Retroespalhamento do material particulado na coluna d'água simulado.
    '''
    bbp = pd.DataFrame()
    for i in bbp_lb0.index:
        bbp[i] = bbp_lb0.loc[i] * np.power(lb0 / wvl, n.loc[i])
    bbp.index = wvl
    return bbp

###############################################################################
def calc_bb(bbp, bbw):
    '''
    Função para calcular o retroespalhamento total na coluna d'água (bb) para
    todos comprimentos de onda.

    ----------
    Parameters
    ----------
    bbp [Data Frame]
        Retroespalhamento do material particulado na coluna d'água.
    bbw [Series / List]
        Retroespalhamento da água.

    -------
    Returns
    -------
    bb [Data Frame]
        Retroespalhamento total na coluna d'água.
    '''
    bb = pd.DataFrame()
    bb = bbw.values + bbp

    return bb
    
###############################################################################
def calc_a_total(bbp, bbw, u):
    '''
    Função para calcular o coeficiente de absorção total na coluna d'água para
    todos os comprimentos de onda.

    ----------
    Parameters
    ----------
    bbp [Data Frame]
        Retroespalhamento do material particulado na coluna d'água.
    bbw [Series / List]
        Retroespalhamento da água.
    u [Series / List]
        Parâmetro u calculado a partir de g0 e g1 - Passos 0 e 1 do modelo QAA.

    -------
    Returns
    -------
    at [Data Frame]
        Coeficientes de absorção total na coluna d'água para todas as estações.
    '''
    # Gerando o DataFrame com os zeros para otimizar o processo
    at = pd.DataFrame()
    for i in list(bbp.index):
            at[i] = ((1 - u.loc[i]) * (bbw.loc[i].values + bbp.loc[i]))/u.loc[i]
            #at_old[i] = ((1 - u.T[i])*(bbw.T[i].values + bbp.T[i]))/u.T[i]
    return at.transpose()

###############################################################################
def calc_a_total_nw(bbp, bbw, u):
    '''
    Função para calcular o coeficiente de absorção total na coluna d'água para
    todos os comprimentos de onda.
    
    Versão 2 = new = nw

    ----------
    Parameters
    ----------
    bbp [Data Frame]
        Retroespalhamento do material particulado na coluna d'água.
    bbw [Series / List]
        Retroespalhamento da água.
    u [Series / List]
        Parâmetro u calculado a partir de g0 e g1 - Passos 0 e 1 do modelo QAA.

    -------
    Returns
    -------
    at [Data Frame]
        Coeficientes de absorção total na coluna d'água para todas as estações.
    '''   
    # Gerando o DataFrame com os zeros para otimizar o processo
    at_nw = pd.DataFrame()
    at_nw = ((1 - u).multiply(bbw.values + bbp))/u
    return at_nw

###############################################################################
################################ FASE 2 #######################################
###############################################################################
'''
As docstrings dessas funções devem ser feitas de forma apadronizar com o resto
dos scripts.

Perguntar pro Rogério depois.
'''
calc_zeta = lambda  r_ref1, r_ref_lb0: 0.74 + (0.2 / (0.8 + (r_ref1 / r_ref_lb0)))

###############################################################################
calc_S = lambda r_ref1, r_ref55x: 0.015 + (0.002/(0.6+(r_ref1/r_ref55x)))

###############################################################################
calc_Xi = lambda S, wvl1, wvl2: np.exp(S * (wvl1 - wvl2))

###############################################################################
def calc_adg443(at, aw, Xi, zeta, wl_r2, wl_r1 = 443):

    at_r1 = at.loc[wl_r1]
    at_r2 = at.loc[wl_r2]
    aw_r1 = aw.loc[wl_r1].values
    aw_r2 = aw.loc[wl_r2].values
    
    return ((at_r2 - zeta*at_r1) - (aw_r2 - zeta * aw_r1)) / (Xi - zeta)

###############################################################################
def calc_adg_total(adg443, S, wl_rx, wl_ref):
    
    adg_total = pd.DataFrame()
    for i in wl_rx:
        adg_total[i] = adg443 * np.exp(-S *( i - wl_ref))
    
    return adg_total.transpose()

###############################################################################
def calc_aphy443(at, aw, Xi, zeta, wl_r2, wl_r1 = 443):

    at_r1 = at.loc[wl_r1]
    at_r2 = at.loc[wl_r2]
    aw_r1 = aw.loc[wl_r1].values
    aw_r2 = aw.loc[wl_r2].values
    
    return ((Xi * at_r1 - at_r2) - (Xi * aw_r2 - aw_r2)) / (Xi - zeta)

###############################################################################
def calc_aphy(at, aw, adg_T):
    
    aphy = pd.DataFrame()
    for i in list(at.index):   
        aphy[i] = at.loc[i] - adg_T.loc[i] -  aw.loc[i].values
    return aphy.transpose()

###############################################################################
################################ FASE Kd ######################################
###############################################################################
def calc_raz_bbw_bb(bbw, bb):
    '''
    Função que calcula a razão entre o retroespalhamento da água (bbw) e o
    retroesplhamento total para uso nas simulações de Kd.
    
    ----------
    Parameters
    ----------
    bbw [Series / List]
        Retroespalhamento da água.
    bb [Data Frame]
        Retroespalhamento total na coluna d'água.
    
    -------
    Returns
    -------
    raz_bbw_bb [Data Frame]
        Data Frame contendo resultado da razão entre o retroespalhamento da água
        e o retroespalhamento total.
    '''
    return bbw.values / bb

###############################################################################
def calc_theta_sun(date):
    '''
    Função que aplica função para calculo dos angulos solares zenitais para cada
    estação/coleta/medição. Tendo como base planilha com dados de coordenadas 
    geográficas, datae horário das coletas.
    
    Função pela data do dia juliano de declinação solar.

    ----------
    Parameters
    ----------
    date [Data Frame]
        Data Frame padronizado contendo nome da estação, coordenadas, data e
        horário das medições.

    -------
    Returns
    -------
    theta_s [Series / List]
        Angulo zenital solar para cada momento de aquisição (Estação).
    '''
    theta_s = pd.DataFrame()
    for i in date:
        lat = mt.radians(date.loc['LAT'][i])
        dd = date.loc['DATE'][i].hour + date.loc['DATE'][i].minute / 60
        h = mt.radians(abs((dd - 12) * 15))
        dy = date.loc['DATE'][i].timetuple().tm_yday
        delta = mt.radians(23.45 * mt.sin((360 / 365) * (dy - 80)))
        cos_theta_s = mt.degrees(mt.acos((mt.sin(lat) * mt.sin(delta)) + (mt.cos(lat) * mt.cos(delta) * mt.cos(h))))
        theta_s[i] = [cos_theta_s]
        
    theta_s.index = ['Theta_Sun']
    return theta_s

###############################################################################
def calc_theta_s(date):
    '''
    Função que aplica função para calculo dos angulos solares zenitais para cada
    estação/coleta/medição. Tendo como base planilha com dados de coordenadas 
    geográficas, datae horário das coletas.
    
    Função pela data do dia juliano de declinação solar.

    ----------
    Parameters
    ----------
    date [Data Frame]
        Data Frame padronizado contendo nome da estação, coordenadas, data e
        horário das medições.

    -------
    Returns
    -------
    theta_s [Series / List]
        Angulo zenital solar para cada momento de aquisição (Estação).
    '''
    theta_s = {
                }
    for i in date.iterrows():
        name = i[0]
        i = i[1]
        lat = mt.radians(i.LAT)
        dd = i.DATE.hour + i.DATE.minute / 60
        h = mt.radians(abs((dd - 12) * 15))
        dy = i.DATE.timetuple().tm_yday
        delta = mt.radians(23.45 * mt.sin((360 / 365) * (dy - 80)))
        cos_theta_s = mt.degrees(mt.acos((mt.sin(lat) * mt.sin(delta)) + (mt.cos(lat) * mt.cos(delta) * mt.cos(h))))
        theta_s[name] = [cos_theta_s]
    
    theta_s = pd.DataFrame(theta_s)
    theta_s.index = ['Theta_Sun']
    
    return theta_s

###############################################################################
def calc_kd_lee_2013(m, theta_s, at, bb, chi, raz_bbw_bb):
    '''
    Função para o calculo do coeficiente de atenualção difusa da irradiância 
    descendenten (Kd) segundo Lee et al. 2013.

    ----------
    Parameters
    ----------
    m [Series / List]
        Conjunto de constantes de ajuste m padrão do modelo de Kd Lee et al. 2013.
    theta_s [Series / List]
        Angulo zenital solar para cada momento de aquisição (Estação).
    at [Data Frame]
        Coeficientes de absorção total na coluna d'água para todas as estações.
    bb [Data Frame]
        Retroespalhamento total na coluna d'água.
    chi [Value]
        Parâmetro de ajuste "chi" padrão do modelo de Kd Lee et al. 2013.
    raz_bbw_bb [Data Frame]
        Data Frame contendo resultado da razão entre o retroespalhamento da água
        e o retroespalhamento total.

    -------
    Returns
    -------
    Kd [Data Frame]
        Coeficiente de atenualção difusa da irradiância descendente na coluna
        d'água.
    '''
    return ((1 + (m[0] * theta_s.values)) * at) + ((1 - (chi * raz_bbw_bb)) * m[1] * (1 - m[2] * np.exp(-m[3] * at)) * bb)
