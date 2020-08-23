'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Modelo QAA v6 integrado ao modelo de Kd com base no ajuste dos Passos 2 e 4
usando bandas do sensor Sentinel 2A/MSI.

Este modelo é separado em duas fases, na primeira é aplicado o modelo "Quasi
Analytical Algorithm" (QAA) de Lee et al. (2014) até seu Passo 6, onde é
calculada o coeficiente de absorção espectral total (at). Na segunda fase é
aplicado o modelo de Kd de Lee et al. (2013) com base nos dados de at simulados
pelo QAA v6.
'''
###############################################################################
#####################          QAA.V6 - model      ############################
###############################################################################
######################## IMPORTS ##############################################
###############################################################################

# Import basic modules from python
import pandas as pd
import os
import math as mt

from qaapy_kdpy.conversion import Rrs_to_rrs
from qaapy_kdpy.file import export_result
from qaapy_kdpy.filter import filter_wvl
from qaapy_kdpy.empirical_calc import *

# QAA - V5 modules
from qaapy_kdpy.QAA_Core.v5.calculation import *

# QAA - V6 modules
from qaapy_kdpy.QAA_Core.v6.calculation import *


def Qaa_v6_model_RUN(data_obj, ano, lb0):
    '''
    Função que aplica o modelo QAA v6 dado um objeto de Data Frames contendo todos
    os dados necessários, uma string com o nome/ano da campanha para salvar o arquivo
    dos resultados e um valor de comprimento de onda inicial (Lambda Zero).
    
    ----------
    Parameters
    ----------
    data_obj [Object]
        Objeto com Data Frames contendo dados usados nas funções do modelo QAA/Kd
        alocados pela função read_files através dos argumentos contendo as strings
        do caminho do diretório para cada planilha '.xlsx' de dados.
        
        O objeto de Data Frames deve conter no mínimo:
            [acs]       Absorção total na coluna d'água;
            
            [Rrs]       Reflectância de sensoriamento remoto superficial;
            
            [aw]        Coeficiente de absorção espectral da água;
            
            [bbw]       Coeficiente de retroespalhamento espectral da água;
            
            [coleta]    Planilha contendo Coordenadas, datas e horários das
                        coletas por estação;
                        
            [Kd_ZEU]    Kd medido com dados até a zona de atenação da luz até 1%;
            
            [Kd_3m]     Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]     Kd medido com dados até 6 metros de profundidade.
        
        * Os dados devem concordar em faixa espectral, bandas e nome das estações.
        
    ano [String]
        String contendo o nome ou ano da campanha, usado para salvar arquivo
        dicionário contendo todos os dados produtos no formato '.pickle'.
        
        Ex.: 'TRM2019', '2013', 'PACO_2011', etc...
        
    lb0 [Value]
        Comprimento de onda incial usado no modelo QAA (Lambda Zero).
        
        Ex.: 665, 560, 490, etc...
                
    -------
    Returns
    -------
    result_dict [Dictionary]
        O resultado desta função é um dicionário contendo todos os resultados e
        dados de entrada relevantes.
        
        Além disso o modelo salva um arquivo contendo todos os dados em formato
        '.pickle' na pasta '.\Results_Espectral'.
        
        Dados salvos:
            [dados/lb0]         Comprimento de onda incial (Lambda Zero);
            
            [dados/wl]          Comprimentos de onda utilizados;
            
            [dados/r_ref]       Dados de rrs;
            
            [dados/campanha]    Nome/ano da campanha;
            
            [dados/equação]     Equação utilizada no modelo em questão;
            
            [dados/g0]          Coeficiente g1 - Passo 1 QAA;
            
            [dados/g1]          Coeficiente g0 - Passo 1 QAA;
            
            [a_Lb]              Absorção em Lambda Zero;
            
            [bb]                Retroespalhamento total;
            
            [at]                Absorção total estimada pelo modelo QAA;
            
            [n]                 Decaimento do retroespalhamento do particulado;
            
            [bbp_lb0]           Retroespalhamento do particulado em Lambda Zero;
            
            [bbp]               Retroespalhamento do particulado;
            
            [Kd_QAA]            Kd estimado pelo modelo QAA / Kd;
            
            [Kd_ZEU]            Kd medido com dados até a zona de atenação da luz
                                até 1%;
                                
            [Kd_3m]             Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]             Kd medido com dados até 6 metros de profundidade.         
    '''
    ###############################################################################
    ########################### QAA SETTINGS - GERAL ##############################
    ###############################################################################
    # data_obj.Rrs = data_obj.Rrs.drop([833, 865, 945, 1374, 1614, 2202])
    # data_obj.aw = data_obj.aw.drop([833, 865, 945, 1374, 1614, 2202])
    # data_obj.bbw = data_obj.bbw.drop([833, 865, 945, 1374, 1614, 2202])
    
    # Transformação do Rrs
    rrs = Rrs_to_rrs(data_obj.Rrs)
    Rrs = data_obj.Rrs
    
    # Cálculando o U
    g0 = 0.089
    g1 = 0.1245
    
    u = calc_u(g0, g1, rrs)

    ##############################################################################################################
    # Filtrando os dados
    wvl_ref = {
            # 'referencia': [443, 492, 560, 665, 704, 741, 783], # Bandas do sensor Sentinel-2A
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779], # Bandas do sensor OLCI
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779,865,885,900], # Bandas do sensor OLCI
            # 'referencia': [400,411,443,490,560,620,667,674,681,709,754,761,779,865,555], # Bandas v6
            # 'referencia': list(range(400, 951, 1)), # Bandas v6
            'referencia': list(range(400, 751, 1)), # Bandas PAR v6
            'lb0': lb0
            }
    
    wl = pd.Series(wvl_ref['referencia'])
    
    r_ref = filter_wvl(rrs, wvl_ref)
    R_ref = filter_wvl(Rrs, wvl_ref)
    
    ###############################################################################
    ############################ FASE 1 - QAA.V6 ##################################
    ###############################################################################
    # QAA v6 considera at_lb0 = a(670)

    aw = data_obj.aw
    
    # a_lb0_v6_492 = calc_alb0_v6_492(492, R_ref, 665, 704, data_obj.aw); lago = 'Eq. QAA.v6 AJUSTE 492 nm'
    a_lb0_v6_560 = calc_alb0_v6_560(560, R_ref, 665, 704, data_obj.aw); lago = 'Eq. QAA.v6 AJUSTE 560 nm'
    # a_lb0_v6 = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 490, data_obj.aw); lago = 'Eq. QAA.v6'
    # a_lb0_v6 = aw.loc[lb0].values; lago = 'Eq. QAA.v6 aw'
    # a_lb0_v6_olci = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 490, data_obj.aw); lago = 'QAA.v6_OLCI'
    # a_lb0_v6_s2a = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 492, data_obj.aw); lago = 'QAA.v6_S2a'

    ##  bbp (lambda_0)
    # bbp_lb0_492 = calc_bbp_lb0_v6(u, data_obj.bbw, 492, a_lb0_v6_492)
    bbp_lb0_560 = calc_bbp_lb0_v6(u, data_obj.bbw, 560, a_lb0_v6_560)
    
    ## Cálculando N
    n_acs_492 = calc_n_acs_492(r_ref[665], r_ref[704])
    # n_acs_560 = calc_n_acs_560(r_ref[665], r_ref[704])
    # n_hydro_560 = calc_n_hydro_560(r_ref[443], r_ref[492])
    # n_hydro_704 = calc_n_hydro_704(r_ref[443], r_ref[492])

    # Realizando cálculo do BBP para todos os comprimentos de onda
    # bbp_492_acs_492 = calc_bbp(wl, n_acs_492, bbp_lb0_492, 492)
    # bbp_492_acs_560 = calc_bbp(wl, n_acs_560, bbp_lb0_492, 492)
    # bbp_492_hydro_560 = calc_bbp(wl, n_hydro_560, bbp_lb0_492, 492)
    # bbp_492_hydro_704 = calc_bbp(wl, n_hydro_704, bbp_lb0_492, 492)
    
    bbp_560_acs_492 = calc_bbp(wl, n_acs_492, bbp_lb0_560, 560)
    # bbp_560_acs_560 = calc_bbp(wl, n_acs_560, bbp_lb0_560, 560)
    # bbp_560_hydro_560 = calc_bbp(wl, n_hydro_560, bbp_lb0_560, 560)
    # bbp_560_hydro_704 = calc_bbp(wl, n_hydro_704, bbp_lb0_560, 560)

    
    # Realizando cálculo do BB para todos comprimentos de onda
    
    # bb_492_acs_492 = calc_bb(bbp_492_acs_492, data_obj.bbw)
    # bb_492_acs_560 = calc_bb(bbp_492_acs_560, data_obj.bbw)
    # bb_492_hydro_560 = calc_bb(bbp_492_hydro_560, data_obj.bbw)
    # bb_492_hydro_704 = calc_bb(bbp_492_hydro_704, data_obj.bbw)
    
    bb_560_acs_492 = calc_bb(bbp_560_acs_492, data_obj.bbw)
    # bb_560_acs_560 = calc_bb(bbp_560_acs_560, data_obj.bbw)
    # bb_560_hydro_560 = calc_bb(bbp_560_hydro_560, data_obj.bbw)
    # bb_560_hydro_704 = calc_bb(bbp_560_hydro_704, data_obj.bbw)
    
    # Realizando cálculo do atotal para todos os comprimentos de onda
    # at_492_acs_492 = calc_a_total(bbp_492_acs_492, data_obj.bbw, u)
    # at_492_acs_560 = calc_a_total(bbp_492_acs_560, data_obj.bbw, u)
    # at_492_hydro_560 = calc_a_total(bbp_492_hydro_560, data_obj.bbw, u)
    # at_492_hydro_704 = calc_a_total(bbp_492_hydro_704, data_obj.bbw, u)
    
    at_560_acs_492 = calc_a_total(bbp_560_acs_492, data_obj.bbw, u)
    # at_560_acs_560 = calc_a_total(bbp_560_acs_560, data_obj.bbw, u)
    # at_560_hydro_560 = calc_a_total(bbp_560_hydro_560, data_obj.bbw, u)
    # at_560_hydro_704 = calc_a_total(bbp_560_hydro_704, data_obj.bbw, u)
    
    ###############################################################################
    ############################ FASE 2 - Kd ######################################
    ###############################################################################
    
    m = [0.005, 4.26, 0.52, 10.54]
    
    chi = 0.265
    
    theta_s = calc_theta_s(data_obj.coleta)
    
    '''
    Renomear tudo de acordo!!!!!!
    Reescrever o modelo e automatizar etapa dos ajustes
    '''
    
    # raz_bbw_bb_492_acs_492 = calc_raz_bbw_bb(data_obj.bbw, bb_492_acs_492)
    # raz_bbw_bb_492_acs_560 = calc_raz_bbw_bb(data_obj.bbw, bb_492_acs_560)
    # raz_bbw_bb_492_hydro_560 = calc_raz_bbw_bb(data_obj.bbw, bb_492_hydro_560)
    # raz_bbw_bb_492_hydro_704 = calc_raz_bbw_bb(data_obj.bbw, bb_492_hydro_704)
    
    raz_bbw_bb_560_acs_492 = calc_raz_bbw_bb(data_obj.bbw, bb_560_acs_492)
    # raz_bbw_bb_560_acs_560 = calc_raz_bbw_bb(data_obj.bbw, bb_560_acs_560)
    # raz_bbw_bb_560_hydro_560 = calc_raz_bbw_bb(data_obj.bbw, bb_560_hydro_560)
    # raz_bbw_bb_560_hydro_704 = calc_raz_bbw_bb(data_obj.bbw, bb_560_hydro_704)
    
    # kd_qaa_492_acs_492 = calc_kd_lee_2013(m, theta_s, at_492_acs_492, bb_492_acs_492, chi, raz_bbw_bb_492_acs_492)
    # kd_qaa_492_acs_560 = calc_kd_lee_2013(m, theta_s, at_492_acs_560, bb_492_acs_560, chi, raz_bbw_bb_492_acs_560)
    # kd_qaa_492_hydro_560 = calc_kd_lee_2013(m, theta_s, at_492_hydro_560, bb_492_hydro_560, chi, raz_bbw_bb_492_hydro_560)
    # kd_qaa_492_hydro_704 = calc_kd_lee_2013(m, theta_s, at_492_hydro_704, bb_492_hydro_704, chi, raz_bbw_bb_492_hydro_704)
    
    kd_qaa_560_acs_492 = calc_kd_lee_2013(m, theta_s, at_560_acs_492, bb_560_acs_492, chi, raz_bbw_bb_560_acs_492)
    # kd_qaa_560_acs_560 = calc_kd_lee_2013(m, theta_s, at_560_acs_560, bb_560_acs_560, chi, raz_bbw_bb_560_acs_560)
    # kd_qaa_560_hydro_560 = calc_kd_lee_2013(m, theta_s, at_560_hydro_560, bb_560_hydro_560, chi, raz_bbw_bb_560_hydro_560)
    # kd_qaa_560_hydro_704 = calc_kd_lee_2013(m, theta_s, at_560_hydro_704, bb_560_hydro_704, chi, raz_bbw_bb_560_hydro_704)
    
    ###############################################################################
    ########################### Teste Kd - absorção insitu ########################
    ###############################################################################
    
    at_acs = data_obj.acs + aw.values
    bb_acs = calc_bb_ref(u, wl, at_acs)
    bb_bb = data_obj.bb
    
    dif = set(theta_s).difference(at_acs)
    theta_s_2 = theta_s.drop(dif, axis = 1)
    
    raz_bbw_bb_acs = calc_raz_bbw_bb(data_obj.bbw, bb_acs)
    raz_bbw_bb_bb = calc_raz_bbw_bb(data_obj.bbw, bb_bb)
    
    kd_acs = calc_kd_lee_2013(m, theta_s_2, at_acs, bb_acs, chi, raz_bbw_bb_acs)
    kd_bb = calc_kd_lee_2013(m, theta_s_2, at_acs, bb_bb, chi, raz_bbw_bb_bb)
    
    ###############################################################################
    ########################### Data Extraction ###################################
    ###############################################################################
    result_dict = {}
    result_dict['dados'] = {
        'lb0': wvl_ref['lb0'],
        'wl': list(r_ref.keys()),
        'r_ref': r_ref,
        'campanhas': ano,
        'equacao': lago,
        'chute_inicial':'não otimizado',
        'h0h1h2': 'não otimizado',
        'g0': g0,
        'g1': g1
    }
    # result_dict['a_lb0_492'] = a_lb0_v6_492
    result_dict['a_lb0_560'] = a_lb0_v6_560
    
    # result_dict['bb_492_acs_492'] = bb_492_acs_492
    # result_dict['bb_492_acs_560'] = bb_492_acs_560
    # result_dict['bb_492_hydro_560'] = bb_492_hydro_560
    # result_dict['bb_492_hydro_704'] = bb_492_hydro_704
    result_dict['bb_560_acs_492'] = bb_560_acs_492
    # result_dict['bb_560_acs_560'] = bb_560_acs_560
    # result_dict['bb_560_hydro_560'] = bb_560_hydro_560
    # result_dict['bb_560_hydro_704'] = bb_560_hydro_704
    
    # result_dict['at_492_acs_492'] = at_492_acs_492
    # result_dict['at_492_acs_560'] = at_492_acs_560
    # result_dict['at_492_hydro_560'] = at_492_hydro_560
    # result_dict['at_492_hydro_704'] = at_492_hydro_704
    result_dict['at_560_acs_492'] = at_560_acs_492
    # result_dict['at_560_acs_560'] = at_560_acs_560
    # result_dict['at_560_hydro_560'] = at_560_hydro_560
    # result_dict['at_560_hydro_704'] = at_560_hydro_704
    
    result_dict['n_acs_492'] = n_acs_492
    # result_dict['n_acs_560'] = n_acs_560
    # result_dict['n_hydro_560'] = n_hydro_560
    # result_dict['n_hydro_704'] = n_hydro_704
    
    # result_dict['bbp_lb0_492'] = bbp_lb0_492
    result_dict['bbp_lb0_560'] = bbp_lb0_560
    
    # result_dict['bbp_492_acs_492'] = bbp_492_acs_492
    # result_dict['bbp_492_acs_560'] = bbp_492_acs_560
    # result_dict['bbp_492_hydro_560'] = bbp_492_hydro_560
    # result_dict['bbp_492_hydro_704'] = bbp_492_hydro_704
    result_dict['bbp_560_acs_492'] = bbp_560_acs_492
    # result_dict['bbp_560_acs_560'] = bbp_560_acs_560
    # result_dict['bbp_560_hydro_560'] = bbp_560_hydro_560
    # result_dict['bbp_560_hydro_704'] = bbp_560_hydro_704
    
    # result_dict['kd_qaa_492_acs_492'] = kd_qaa_492_acs_492
    # result_dict['kd_qaa_492_acs_560'] = kd_qaa_492_acs_560
    # result_dict['kd_qaa_492_hydro_560'] = kd_qaa_492_hydro_560
    # result_dict['kd_qaa_492_hydro_704'] = kd_qaa_492_hydro_704
    result_dict['kd_qaa_560_acs_492'] = kd_qaa_560_acs_492
    # result_dict['kd_qaa_560_acs_560'] = kd_qaa_560_acs_560
    # result_dict['kd_qaa_560_hydro_560'] = kd_qaa_560_hydro_560
    # result_dict['kd_qaa_560_hydro_704'] = kd_qaa_560_hydro_704
    
    result_dict['kd_acs'] = kd_acs
    result_dict['kd_bb'] = kd_bb
    
    result_dict['Kd_ZEU'] = data_obj.kd_zeu
    result_dict['Kd_3m'] = data_obj.kd_3m
    result_dict['Kd_6m'] = data_obj.kd_6m

    os.makedirs('./Results_Espectral', exist_ok = True)
    export_result(result_dict, './Results_Espectral/QAAv6_' + ano + '_' + str(wvl_ref['lb0']) + '_Ajustado.pickle')

    return result_dict

def Qaa_v6_model_msi(Rrs, wvl, b1, b2, b3, b4, b5, lb0, aw, bbw, theta_s):
    '''
    Função que aplica o modelo QAA v6 dado um objeto de Data Frames contendo todos
    os dados necessários, uma string com o nome/ano da campanha para salvar o arquivo
    dos resultados e um valor de comprimento de onda inicial (Lambda Zero).
    
    ----------
    Parameters
    ----------
    data_obj [Object]
        Objeto com Data Frames contendo dados usados nas funções do modelo QAA/Kd
        alocados pela função read_files através dos argumentos contendo as strings
        do caminho do diretório para cada planilha '.xlsx' de dados.
        
        O objeto de Data Frames deve conter no mínimo:
            [acs]       Absorção total na coluna d'água;
            
            [Rrs]       Reflectância de sensoriamento remoto superficial;
            
            [aw]        Coeficiente de absorção espectral da água;
            
            [bbw]       Coeficiente de retroespalhamento espectral da água;
            
            [coleta]    Planilha contendo Coordenadas, datas e horários das
                        coletas por estação;
                        
            [Kd_ZEU]    Kd medido com dados até a zona de atenação da luz até 1%;
            
            [Kd_3m]     Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]     Kd medido com dados até 6 metros de profundidade.
        
        * Os dados devem concordar em faixa espectral, bandas e nome das estações.
        
    ano [String]
        String contendo o nome ou ano da campanha, usado para salvar arquivo
        dicionário contendo todos os dados produtos no formato '.pickle'.
        
        Ex.: 'TRM2019', '2013', 'PACO_2011', etc...
        
    lb0 [Value]
        Comprimento de onda incial usado no modelo QAA (Lambda Zero).
        
        Ex.: 665, 560, 490, etc...
                
    -------
    Returns
    -------
    result_dict [Dictionary]
        O resultado desta função é um dicionário contendo todos os resultados e
        dados de entrada relevantes.
        
        Além disso o modelo salva um arquivo contendo todos os dados em formato
        '.pickle' na pasta '.\Results_Espectral'.
        
        Dados salvos:
            [dados/lb0]         Comprimento de onda incial (Lambda Zero);
            
            [dados/wl]          Comprimentos de onda utilizados;
            
            [dados/r_ref]       Dados de rrs;
            
            [dados/campanha]    Nome/ano da campanha;
            
            [dados/equação]     Equação utilizada no modelo em questão;
            
            [dados/g0]          Coeficiente g1 - Passo 1 QAA;
            
            [dados/g1]          Coeficiente g0 - Passo 1 QAA;
            
            [a_Lb]              Absorção em Lambda Zero;
            
            [bb]                Retroespalhamento total;
            
            [at]                Absorção total estimada pelo modelo QAA;
            
            [n]                 Decaimento do retroespalhamento do particulado;
            
            [bbp_lb0]           Retroespalhamento do particulado em Lambda Zero;
            
            [bbp]               Retroespalhamento do particulado;
            
            [Kd_QAA]            Kd estimado pelo modelo QAA / Kd;
            
            [Kd_ZEU]            Kd medido com dados até a zona de atenação da luz
                                até 1%;
                                
            [Kd_3m]             Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]             Kd medido com dados até 6 metros de profundidade.         
    '''
    ###############################################################################
    ########################### QAA SETTINGS - GERAL ##############################
    ###############################################################################
    aw = aw.drop([741, 783])
    bbw = bbw.drop([741, 783])
    
    # Transformação do Rrs
    #   rrs = Rrs_to_rrs(Rrs)
    rrs = Rrs / (0.52 + 1.7 * Rrs)
    #   rrs_lb0 = Rrs_to_rrs(b3)
    rrs_lb0 = b3 / (0.52 + 1.7 * b3)
    
    rrs_b4 = b4 / (0.52 + 1.7 * b4)
    rrs_b5 = b5 / (0.52 + 1.7 * b5)
    
    # Cálculando o U
    g0 = 0.089
    g1 = 0.1245
    
    u = calc_u(g0, g1, rrs)
    u_lb0 = calc_u(g0, g1, rrs_lb0)
    
    ###############################################################################
    ############################ FASE 1 - QAA.V6 ##################################
    ###############################################################################    
    a_lb0_v6_560 = aw.loc[lb0].values + 0.4310 * (b3 / (b4 + b5)) ** -1.4408; lago = 'Eq. QAA.v6 AJUSTE 560 nm'

    ##  bbp (lambda_0)
    bbp_lb0_560 = (( u_lb0 * a_lb0_v6_560 ) / (1 - u_lb0)) - bbw.loc[lb0].values
    
    ## Cálculando N
    n_acs_492 = (0.5248 * np.exp(rrs_b4 / rrs_b5)) - 1.1849    
    
    bbp_560_acs_492 = bbp_lb0_560 * np.power(lb0 / wvl, n_acs_492)

    # Realizando cálculo do BB para todos comprimentos de onda    
    bb_560_acs_492 = bbw.loc[wvl].values + bbp_560_acs_492
    
    # Realizando cálculo do atotal para todos os comprimentos de onda
    at_560_acs_492 = ((1 - u) * (bbw.loc[wvl].values + bbp_560_acs_492)) / u
    
    ###############################################################################
    ############################ FASE 2 - Kd ######################################
    ###############################################################################
    
    m = [0.005, 4.26, 0.52, 10.54]
    
    chi = 0.265
    
    #theta_s = calc_theta_s(data_obj.coleta)
    
    raz_bbw_bb_560_acs_492 = bbw.loc[wvl].values / bb_560_acs_492
    
    kd_qaa_560_acs_492 = ((1 + (m[0] * theta_s)) * at_560_acs_492) + ((1 - (chi * raz_bbw_bb_560_acs_492)) * m[1] * (1 - m[2] * np.exp(-m[3] * at_560_acs_492)) * bb_560_acs_492)   
    
    ###############################################################################
    ########################### Data Extraction ###################################
    ###############################################################################
    
    return kd_qaa_560_acs_492

def Qaa_v6_model_RUN_original(data_obj, ano, lb0):
    '''
    Função que aplica o modelo QAA v6 dado um objeto de Data Frames contendo todos
    os dados necessários, uma string com o nome/ano da campanha para salvar o arquivo
    dos resultados e um valor de comprimento de onda inicial (Lambda Zero).
    
    ----------
    Parameters
    ----------
    data_obj [Object]
        Objeto com Data Frames contendo dados usados nas funções do modelo QAA/Kd
        alocados pela função read_files através dos argumentos contendo as strings
        do caminho do diretório para cada planilha '.xlsx' de dados.
        
        O objeto de Data Frames deve conter no mínimo:
            [acs]       Absorção total na coluna d'água;
            
            [Rrs]       Reflectância de sensoriamento remoto superficial;
            
            [aw]        Coeficiente de absorção espectral da água;
            
            [bbw]       Coeficiente de retroespalhamento espectral da água;
            
            [coleta]    Planilha contendo Coordenadas, datas e horários das
                        coletas por estação;
                        
            [Kd_ZEU]    Kd medido com dados até a zona de atenação da luz até 1%;
            
            [Kd_3m]     Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]     Kd medido com dados até 6 metros de profundidade.
        
        * Os dados devem concordar em faixa espectral, bandas e nome das estações.
        
    ano [String]
        String contendo o nome ou ano da campanha, usado para salvar arquivo
        dicionário contendo todos os dados produtos no formato '.pickle'.
        
        Ex.: 'TRM2019', '2013', 'PACO_2011', etc...
        
    lb0 [Value]
        Comprimento de onda incial usado no modelo QAA (Lambda Zero).
        
        Ex.: 665, 560, 490, etc...
                
    -------
    Returns
    -------
    result_dict [Dictionary]
        O resultado desta função é um dicionário contendo todos os resultados e
        dados de entrada relevantes.
        
        Além disso o modelo salva um arquivo contendo todos os dados em formato
        '.pickle' na pasta '.\Results_Espectral'.
        
        Dados salvos:
            [dados/lb0]         Comprimento de onda incial (Lambda Zero);
            
            [dados/wl]          Comprimentos de onda utilizados;
            
            [dados/r_ref]       Dados de rrs;
            
            [dados/campanha]    Nome/ano da campanha;
            
            [dados/equação]     Equação utilizada no modelo em questão;
            
            [dados/g0]          Coeficiente g1 - Passo 1 QAA;
            
            [dados/g1]          Coeficiente g0 - Passo 1 QAA;
            
            [a_Lb]              Absorção em Lambda Zero;
            
            [bb]                Retroespalhamento total;
            
            [at]                Absorção total estimada pelo modelo QAA;
            
            [n]                 Decaimento do retroespalhamento do particulado;
            
            [bbp_lb0]           Retroespalhamento do particulado em Lambda Zero;
            
            [bbp]               Retroespalhamento do particulado;
            
            [Kd_QAA]            Kd estimado pelo modelo QAA / Kd;
            
            [Kd_ZEU]            Kd medido com dados até a zona de atenação da luz
                                até 1%;
                                
            [Kd_3m]             Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]             Kd medido com dados até 6 metros de profundidade.         
    '''
    ###############################################################################
    ########################### QAA SETTINGS - GERAL ##############################
    ###############################################################################
    # data_obj.Rrs = data_obj.Rrs.drop([833, 865, 945, 1374, 1614, 2202])
    # data_obj.aw = data_obj.aw.drop([833, 865, 945, 1374, 1614, 2202])
    # data_obj.bbw = data_obj.bbw.drop([833, 865, 945, 1374, 1614, 2202])
    
    # Transformação do Rrs
    rrs = Rrs_to_rrs(data_obj.Rrs)
    Rrs = data_obj.Rrs
    
    # Cálculando o U
    g0 = 0.089
    g1 = 0.1245
    
    u = calc_u(g0, g1, rrs)

    ##############################################################################################################
    # Filtrando os dados
    wvl_ref = {
            # 'referencia': [443, 492, 560, 665, 704, 741, 783], # Bandas do sensor Sentinel-2A
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779], # Bandas do sensor OLCI
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779,865,885,900], # Bandas do sensor OLCI
            # 'referencia': [400,411,443,490,560,620,667,674,681,709,754,761,779,865,555], # Bandas v6
            # 'referencia': list(range(400, 951, 1)), # Bandas v6
            'referencia': list(range(400, 751, 1)), # Bandas PAR v6
            'lb0': lb0
            }
    
    wl = pd.Series(wvl_ref['referencia'])
    
    r_ref = filter_wvl(rrs, wvl_ref)
    R_ref = filter_wvl(Rrs, wvl_ref)
    
    ###############################################################################
    ############################ FASE 1 - QAA.V6 ##################################
    ###############################################################################
    # QAA v6 considera at_lb0 = a(670)

    aw = data_obj.aw
    
    # a_lb0_v6_492 = calc_alb0_v6_492(492, R_ref, 665, 704, data_obj.aw); lago = 'Eq. QAA.v6 AJUSTE 492 nm'
    # a_lb0_v6_560 = calc_alb0_v6_560(560, R_ref, 665, 704, data_obj.aw); lago = 'Eq. QAA.v6 AJUSTE 560 nm'
    # a_lb0_v6 = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 490, data_obj.aw); lago = 'Eq. QAA.v6'
    # a_lb0_v6 = aw.loc[lb0].values; lago = 'Eq. QAA.v6 aw'
    # a_lb0_v6_olci = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 490, data_obj.aw); lago = 'QAA.v6_OLCI'
    a_lb0_v6_s2a = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 492, data_obj.aw); lago = 'QAA.v6_S2a'

    ##  bbp (lambda_0)
    # bbp_lb0_492 = calc_bbp_lb0_v6(u, data_obj.bbw, 492, a_lb0_v6_492)
    # bbp_lb0_560 = calc_bbp_lb0_v6(u, data_obj.bbw, 560, a_lb0_v6_560)
    bbp_lb0_665 = calc_bbp_lb0_v6(u, data_obj.bbw, 665, a_lb0_v6_s2a)
    
    ## Cálculando N
    # n_acs_492 = calc_n_acs_492(r_ref[665], r_ref[704])
    # n_acs_560 = calc_n_acs_560(r_ref[665], r_ref[704])
    # n_hydro_560 = calc_n_hydro_560(r_ref[443], r_ref[492])
    # n_hydro_704 = calc_n_hydro_704(r_ref[443], r_ref[492])
    n_slope = calc_n(r_ref[443], r_ref[560])

    # Realizando cálculo do BBP para todos os comprimentos de onda
    # bbp_492_acs_492 = calc_bbp(wl, n_acs_492, bbp_lb0_492, 492)
    # bbp_492_acs_560 = calc_bbp(wl, n_acs_560, bbp_lb0_492, 492)
    # bbp_492_hydro_560 = calc_bbp(wl, n_hydro_560, bbp_lb0_492, 492)
    # bbp_492_hydro_704 = calc_bbp(wl, n_hydro_704, bbp_lb0_492, 492)
    # bbp_560_acs_492 = calc_bbp(wl, n_acs_492, bbp_lb0_560, 560)
    # bbp_560_acs_560 = calc_bbp(wl, n_acs_560, bbp_lb0_560, 560)
    # bbp_560_hydro_560 = calc_bbp(wl, n_hydro_560, bbp_lb0_560, 560)
    # bbp_560_hydro_704 = calc_bbp(wl, n_hydro_704, bbp_lb0_560, 560)
    bbp_665 = calc_bbp(wl, n_slope, bbp_lb0_665, 665)

    
    # Realizando cálculo do BB para todos comprimentos de onda
    # bb_492_acs_492 = calc_bb(bbp_492_acs_492, data_obj.bbw)
    # bb_492_acs_560 = calc_bb(bbp_492_acs_560, data_obj.bbw)
    # bb_492_hydro_560 = calc_bb(bbp_492_hydro_560, data_obj.bbw)
    # bb_492_hydro_704 = calc_bb(bbp_492_hydro_704, data_obj.bbw)
    # bb_560_acs_492 = calc_bb(bbp_560_acs_492, data_obj.bbw)
    # bb_560_acs_560 = calc_bb(bbp_560_acs_560, data_obj.bbw)
    # bb_560_hydro_560 = calc_bb(bbp_560_hydro_560, data_obj.bbw)
    # bb_560_hydro_704 = calc_bb(bbp_560_hydro_704, data_obj.bbw)
    bb_665 = calc_bb(bbp_665, data_obj.bbw)
    
    # Realizando cálculo do atotal para todos os comprimentos de onda
    # at_492_acs_492 = calc_a_total(bbp_492_acs_492, data_obj.bbw, u)
    # at_492_acs_560 = calc_a_total(bbp_492_acs_560, data_obj.bbw, u)
    # at_492_hydro_560 = calc_a_total(bbp_492_hydro_560, data_obj.bbw, u)
    # at_492_hydro_704 = calc_a_total(bbp_492_hydro_704, data_obj.bbw, u)
    # at_560_acs_492 = calc_a_total(bbp_560_acs_492, data_obj.bbw, u)
    # at_560_acs_560 = calc_a_total(bbp_560_acs_560, data_obj.bbw, u)
    # at_560_hydro_560 = calc_a_total(bbp_560_hydro_560, data_obj.bbw, u)
    # at_560_hydro_704 = calc_a_total(bbp_560_hydro_704, data_obj.bbw, u)
    at_665 = calc_a_total(bbp_665, data_obj.bbw, u)
    
    ###############################################################################
    ############################ FASE 2 - Kd ######################################
    ###############################################################################
    
    m = [0.005, 4.26, 0.52, 10.54]
    
    chi = 0.265
    
    theta_s = calc_theta_s(data_obj.coleta)
    
    '''
    Renomear tudo de acordo!!!!!!
    Reescrever o modelo e automatizar etapa dos ajustes
    '''
    
    # raz_bbw_bb_492_acs_492 = calc_raz_bbw_bb(data_obj.bbw, bb_492_acs_492)
    # raz_bbw_bb_492_acs_560 = calc_raz_bbw_bb(data_obj.bbw, bb_492_acs_560)
    # raz_bbw_bb_492_hydro_560 = calc_raz_bbw_bb(data_obj.bbw, bb_492_hydro_560)
    # raz_bbw_bb_492_hydro_704 = calc_raz_bbw_bb(data_obj.bbw, bb_492_hydro_704)
    # raz_bbw_bb_560_acs_492 = calc_raz_bbw_bb(data_obj.bbw, bb_560_acs_492)
    # raz_bbw_bb_560_acs_560 = calc_raz_bbw_bb(data_obj.bbw, bb_560_acs_560)
    # raz_bbw_bb_560_hydro_560 = calc_raz_bbw_bb(data_obj.bbw, bb_560_hydro_560)
    # raz_bbw_bb_560_hydro_704 = calc_raz_bbw_bb(data_obj.bbw, bb_560_hydro_704)
    raz_bbw_bb_665 = calc_raz_bbw_bb(data_obj.bbw, bb_665)
    
    # kd_qaa_492_acs_492 = calc_kd_lee_2013(m, theta_s, at_492_acs_492, bb_492_acs_492, chi, raz_bbw_bb_492_acs_492)
    # kd_qaa_492_acs_560 = calc_kd_lee_2013(m, theta_s, at_492_acs_560, bb_492_acs_560, chi, raz_bbw_bb_492_acs_560)
    # kd_qaa_492_hydro_560 = calc_kd_lee_2013(m, theta_s, at_492_hydro_560, bb_492_hydro_560, chi, raz_bbw_bb_492_hydro_560)
    # kd_qaa_492_hydro_704 = calc_kd_lee_2013(m, theta_s, at_492_hydro_704, bb_492_hydro_704, chi, raz_bbw_bb_492_hydro_704)
    # kd_qaa_560_acs_492 = calc_kd_lee_2013(m, theta_s, at_560_acs_492, bb_560_acs_492, chi, raz_bbw_bb_560_acs_492)
    # kd_qaa_560_acs_560 = calc_kd_lee_2013(m, theta_s, at_560_acs_560, bb_560_acs_560, chi, raz_bbw_bb_560_acs_560)
    # kd_qaa_560_hydro_560 = calc_kd_lee_2013(m, theta_s, at_560_hydro_560, bb_560_hydro_560, chi, raz_bbw_bb_560_hydro_560)
    # kd_qaa_560_hydro_704 = calc_kd_lee_2013(m, theta_s, at_560_hydro_704, bb_560_hydro_704, chi, raz_bbw_bb_560_hydro_704)
    kd_qaa_665 = calc_kd_lee_2013(m, theta_s, at_665, bb_665, chi, raz_bbw_bb_665)
    
    ###############################################################################
    ########################### Teste Kd - absorção insitu ########################
    ###############################################################################
    
    '''
    at_acs = data_obj.acs + aw.values
    bb_acs = calc_bb_ref(u, wl, at_acs)
    bb_bb = data_obj.bb
    
    dif = set(theta_s).difference(at_acs)
    theta_s_2 = theta_s.drop(dif, axis = 1)
    
    raz_bbw_bb_acs = calc_raz_bbw_bb(data_obj.bbw, bb_acs)
    raz_bbw_bb_bb = calc_raz_bbw_bb(data_obj.bbw, bb_bb)
    
    kd_acs = calc_kd_lee_2013(m, theta_s_2, at_acs, bb_acs, chi, raz_bbw_bb_acs)
    kd_bb = calc_kd_lee_2013(m, theta_s_2, at_acs, bb_bb, chi, raz_bbw_bb_bb)
    
    '''
    ###############################################################################
    ########################### Data Extraction ###################################
    ###############################################################################
    result_dict = {}
    result_dict['dados'] = {
        'lb0': wvl_ref['lb0'],
        'wl': list(r_ref.keys()),
        'r_ref': r_ref,
        'campanhas': ano,
        'equacao': lago,
        'chute_inicial':'não otimizado',
        'h0h1h2': 'não otimizado',
        'g0': g0,
        'g1': g1
    }
    # result_dict['a_lb0_492'] = a_lb0_v6_492
    # result_dict['a_lb0_560'] = a_lb0_v6_560
    result_dict['a_lb0_665'] = a_lb0_v6_s2a
    
    # result_dict['bb_492_acs_492'] = bb_492_acs_492
    # result_dict['bb_492_acs_560'] = bb_492_acs_560
    # result_dict['bb_492_hydro_560'] = bb_492_hydro_560
    # result_dict['bb_492_hydro_704'] = bb_492_hydro_704
    # result_dict['bb_560_acs_492'] = bb_560_acs_492
    # result_dict['bb_560_acs_560'] = bb_560_acs_560
    # result_dict['bb_560_hydro_560'] = bb_560_hydro_560
    # result_dict['bb_560_hydro_704'] = bb_560_hydro_704
    result_dict['bb_665'] = bb_665
    
    # result_dict['at_492_acs_492'] = at_492_acs_492
    # result_dict['at_492_acs_560'] = at_492_acs_560
    # result_dict['at_492_hydro_560'] = at_492_hydro_560
    # result_dict['at_492_hydro_704'] = at_492_hydro_704
    # result_dict['at_560_acs_492'] = at_560_acs_492
    # result_dict['at_560_acs_560'] = at_560_acs_560
    # result_dict['at_560_hydro_560'] = at_560_hydro_560
    # result_dict['at_560_hydro_704'] = at_560_hydro_704
    result_dict['at_665'] = at_665
    
    # result_dict['n_acs_492'] = n_acs_492
    # result_dict['n_acs_560'] = n_acs_560
    # result_dict['n_hydro_560'] = n_hydro_560
    # result_dict['n_hydro_704'] = n_hydro_704
    result_dict['n_665'] = n_slope
    
    # result_dict['bbp_lb0_492'] = bbp_lb0_492
    # result_dict['bbp_lb0_560'] = bbp_lb0_560
    result_dict['bbp_lb0_665'] = bbp_lb0_665
    
    # result_dict['bbp_492_acs_492'] = bbp_492_acs_492
    # result_dict['bbp_492_acs_560'] = bbp_492_acs_560
    # result_dict['bbp_492_hydro_560'] = bbp_492_hydro_560
    # result_dict['bbp_492_hydro_704'] = bbp_492_hydro_704
    # result_dict['bbp_560_acs_492'] = bbp_560_acs_492
    # result_dict['bbp_560_acs_560'] = bbp_560_acs_560
    # result_dict['bbp_560_hydro_560'] = bbp_560_hydro_560
    # result_dict['bbp_560_hydro_704'] = bbp_560_hydro_704
    result_dict['bbp_665'] = bbp_665
    
    # result_dict['kd_qaa_492_acs_492'] = kd_qaa_492_acs_492
    # result_dict['kd_qaa_492_acs_560'] = kd_qaa_492_acs_560
    # result_dict['kd_qaa_492_hydro_560'] = kd_qaa_492_hydro_560
    # result_dict['kd_qaa_492_hydro_704'] = kd_qaa_492_hydro_704
    # result_dict['kd_qaa_560_acs_492'] = kd_qaa_560_acs_492
    # result_dict['kd_qaa_560_acs_560'] = kd_qaa_560_acs_560
    # result_dict['kd_qaa_560_hydro_560'] = kd_qaa_560_hydro_560
    # result_dict['kd_qaa_560_hydro_704'] = kd_qaa_560_hydro_704
    result_dict['kd_qaa_665'] = kd_qaa_665
    
    # result_dict['kd_acs'] = kd_acs
    # result_dict['kd_bb'] = kd_bb
    
    result_dict['Kd_ZEU'] = data_obj.kd_zeu
    result_dict['Kd_3m'] = data_obj.kd_3m
    result_dict['Kd_6m'] = data_obj.kd_6m

    os.makedirs('./Results_Espectral', exist_ok = True)
    export_result(result_dict, './Results_Espectral/QAAv6_' + ano + '_' + str(wvl_ref['lb0']) + '_Ajustado.pickle')

    return result_dict

def Qaa_v6_model_msi_original(Rrs, wvl, b1, b2, b3, b4, b5, lb0, aw, bbw, theta_s):
    '''
    Função que aplica o modelo QAA v6 dado um objeto de Data Frames contendo todos
    os dados necessários, uma string com o nome/ano da campanha para salvar o arquivo
    dos resultados e um valor de comprimento de onda inicial (Lambda Zero).
    
    ----------
    Parameters
    ----------
    data_obj [Object]
        Objeto com Data Frames contendo dados usados nas funções do modelo QAA/Kd
        alocados pela função read_files através dos argumentos contendo as strings
        do caminho do diretório para cada planilha '.xlsx' de dados.
        
        O objeto de Data Frames deve conter no mínimo:
            [acs]       Absorção total na coluna d'água;
            
            [Rrs]       Reflectância de sensoriamento remoto superficial;
            
            [aw]        Coeficiente de absorção espectral da água;
            
            [bbw]       Coeficiente de retroespalhamento espectral da água;
            
            [coleta]    Planilha contendo Coordenadas, datas e horários das
                        coletas por estação;
                        
            [Kd_ZEU]    Kd medido com dados até a zona de atenação da luz até 1%;
            
            [Kd_3m]     Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]     Kd medido com dados até 6 metros de profundidade.
        
        * Os dados devem concordar em faixa espectral, bandas e nome das estações.
        
    ano [String]
        String contendo o nome ou ano da campanha, usado para salvar arquivo
        dicionário contendo todos os dados produtos no formato '.pickle'.
        
        Ex.: 'TRM2019', '2013', 'PACO_2011', etc...
        
    lb0 [Value]
        Comprimento de onda incial usado no modelo QAA (Lambda Zero).
        
        Ex.: 665, 560, 490, etc...
                
    -------
    Returns
    -------
    result_dict [Dictionary]
        O resultado desta função é um dicionário contendo todos os resultados e
        dados de entrada relevantes.
        
        Além disso o modelo salva um arquivo contendo todos os dados em formato
        '.pickle' na pasta '.\Results_Espectral'.
        
        Dados salvos:
            [dados/lb0]         Comprimento de onda incial (Lambda Zero);
            
            [dados/wl]          Comprimentos de onda utilizados;
            
            [dados/r_ref]       Dados de rrs;
            
            [dados/campanha]    Nome/ano da campanha;
            
            [dados/equação]     Equação utilizada no modelo em questão;
            
            [dados/g0]          Coeficiente g1 - Passo 1 QAA;
            
            [dados/g1]          Coeficiente g0 - Passo 1 QAA;
            
            [a_Lb]              Absorção em Lambda Zero;
            
            [bb]                Retroespalhamento total;
            
            [at]                Absorção total estimada pelo modelo QAA;
            
            [n]                 Decaimento do retroespalhamento do particulado;
            
            [bbp_lb0]           Retroespalhamento do particulado em Lambda Zero;
            
            [bbp]               Retroespalhamento do particulado;
            
            [Kd_QAA]            Kd estimado pelo modelo QAA / Kd;
            
            [Kd_ZEU]            Kd medido com dados até a zona de atenação da luz
                                até 1%;
                                
            [Kd_3m]             Kd medido com dados até 3 metros de profundidade;
            
            [Kd_6m]             Kd medido com dados até 6 metros de profundidade.         
    '''
    ###############################################################################
    ########################### QAA SETTINGS - GERAL ##############################
    ###############################################################################
    aw = aw.drop([741, 783])
    bbw = bbw.drop([741, 783])
    
    # Transformação do Rrs
    #   rrs = Rrs_to_rrs(Rrs)
    rrs = Rrs / (0.52 + 1.7 * Rrs)
    #   rrs_lb0 = Rrs_to_rrs(b3)
    rrs_lb0 = b4 / (0.52 + 1.7 * b4)
    
    rrs_b1 = b1 / (0.52 + 1.7 * b1)
    rrs_b3 = b3 / (0.52 + 1.7 * b3)
    
    # Cálculando o U
    g0 = 0.089
    g1 = 0.1245
    
    u = calc_u(g0, g1, rrs)
    u_lb0 = calc_u(g0, g1, rrs_lb0)
    
    ###############################################################################
    ############################ FASE 1 - QAA.V6 ##################################
    ###############################################################################    
    a_lb0_v6_665 = aw.loc[lb0].values + 0.39 * (b4 / (b1 + b2)) ** 1.14; lago = 'Eq. QAA.v6 Original 665 nm'

    ##  bbp (lambda_0)
    bbp_lb0_665 = (( u_lb0 * a_lb0_v6_665 ) / (1 - u_lb0)) - bbw.loc[lb0].values
    
    ## Cálculando N
    n_665 = 2 * (1 - 2* np.exp(-0.9 * (rrs_b1 / rrs_b3))) 
    
    bbp_665 = bbp_lb0_665 * np.power(lb0 / wvl, n_665)

    # Realizando cálculo do BB para todos comprimentos de onda    
    bb_665 = bbw.loc[wvl].values + bbp_665
    
    # Realizando cálculo do atotal para todos os comprimentos de onda
    at_665 = ((1 - u) * (bbw.loc[wvl].values + bbp_665)) / u
    
    ###############################################################################
    ############################ FASE 2 - Kd ######################################
    ###############################################################################
    
    m = [0.005, 4.26, 0.52, 10.54]
    
    chi = 0.265
    
    #theta_s = calc_theta_s(data_obj.coleta)
    
    raz_bbw_bb_665 = bbw.loc[wvl].values / bb_665
    
    kd_qaa_665 = ((1 + (m[0] * theta_s)) * at_665) + ((1 - (chi * raz_bbw_bb_665)) * m[1] * (1 - m[2] * np.exp(-m[3] * at_665)) * bb_665)
    
    ###############################################################################
    ########################### Data Extraction ###################################
    ###############################################################################
    
    return kd_qaa_665
