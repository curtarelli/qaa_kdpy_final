'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Modelo QAA v5 integrado ao modelo de Kd com base no ajuste dos Passos 2 e 4
usando bandas do sensor Sentinel 2A/MSI.

Este modelo é separado em duas fases, na primeira é aplicado o modelo "Quasi
Analytical Algorithm" (QAA) de Lee et al. (2013a) até seu Passo 6, onde é
calculada o coeficiente de absorção espectral total (at). Na segunda fase é
aplicado o modelo de Kd de Lee et al. (2013b) com base nos dados de at simulados
pelo QAA v5.
'''

###############################################################################
#############           QAA.V5 - lee 2013 - Model         #####################
###############################################################################
######################## IMPORTS ##############################################
###############################################################################

# Import basic modules from python
import numpy as np
import pandas as pd
import os

# Qaapy basic modules
from qaapy_kdpy.conversion import Rrs_to_rrs
from qaapy_kdpy.file import export_result
from qaapy_kdpy.optimizer import optimizer_aLb0
from qaapy_kdpy.filter import filter_wvl

# QAA - V5 modules
from qaapy_kdpy.QAA_Core.v5.calculation import *

def Qaa_v5_model_RUN(data_obj, ano, lb0):
    '''
    Função que aplica o modelo QAA v5 dado um objeto de Data Frames contendo todos
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
    ########################### QAA SETTINGS - GENERAL ############################
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
    
    ###############################################################################
    #ToDo: checar se isso é necessário aqui, talvez criar uma função que padronize isso fora dos OLCI_models.
    # Filtrando os dados
    wvl_ref = {
            # 'referencia': [443, 492, 560, 665, 704, 741, 783], # Bandas do sensor Sentinel-2A
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779], # Bandas do sensor OLCI
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779,865,885,900], # Bandas do sensor OLCI
            # 'referencia': [400,411,443,490,560,620,667,674,681,709,754,761,779,865,555], # Bandas v5
            # 'referencia': list(range(400, 951, 1)), # Bandas v5
            'referencia': list(range(400, 751, 1)), # Bandas PAR v5
            'lb0': lb0
            }

    wl = pd.Series(wvl_ref['referencia'])
    
    r_ref = filter_wvl(rrs, wvl_ref)
    # R_ref = filter_wvl(Rrs, wvl_ref)
    
    ###############################################################################
    ############################ FASE 1 - QAA.V5 - at #############################
    ###############################################################################
    
    p0 = [-1.146, -1.336,- 0.469] # Valores iniciais
    
    # eq_lee_v5_OLCI = np.log10( (r_ref[443] + r_ref[490]) / ( r_ref['r_lb0'] + 5 * (r_ref[665]/r_ref[510] * r_ref[665]))); lago = 'QAA.v5_OLCI'
    # eq_lee_v5_s2a = np.log10( (r_ref[443] + r_ref[492]) / ( r_ref['r_lb0'] + 5 * (r_ref[665]/r_ref[492] * r_ref[665]))); lago = 'QAA.v5_S2a'
    eq_lee_v5 = np.log10( (r_ref[443] + r_ref[490]) / (r_ref['r_lb0'] + 5 * (r_ref[670] / r_ref[490] * r_ref[670]))); lago = 'QAA.v5'

    aw = data_obj.aw

    #### ATENÇÃO : a agua é adicionada ao acs dentro da função "optimizer_aLb0" ##########
    # a_Lb0_optimized = optimizer_aLb0(eq_lee_v5, wvl_ref['lb0'], p0, data_obj.aw, data_obj.acs, NOToptimize = True); opt = '_OptN'
    # a_Lb0_optimized = optimizer_aLb0(eq_lee_v5, wvl_ref['lb0'], p0, data_obj.aw, data_obj.acs, NOToptimize = False); opt = '_OptY'
    a_Lb0_optimized = aw.loc[lb0].values; opt = '_OptAw'
    
    #### bbp (lambda_0)
    # bbp_lb0 = calc_bbp_lb0(u, data_obj.bbw, 555, wvl_ref['lb0'], a_Lb0_optimized['aref'])
    bbp_lb0 = calc_bbp_lb0(u, data_obj.bbw, wvl_ref['lb0'], a_Lb0_optimized)
    
    ## Cálculando N
    n = calc_n(r_ref[443], r_ref[555])
    
    # Realizando cálculo do BBP para todos os comprimentos de onda
    bbp = calc_bbp(wl, n, bbp_lb0, wvl_ref['lb0'])
    
    # Realizando cálculo do BB para todos comprimentos de onda
    bb = calc_bb(bbp, data_obj.bbw)
    
    # Realizando cálculo do atotal para todos os comprimentos de onda
    at = calc_a_total_nw(bbp, data_obj.bbw, u)
    
    ###############################################################################
    ############################ FASE 2 - Kd ######################################
    ###############################################################################
    
    m = [0.005, 4.26, 0.52, 10.54]
    
    chi = 0.265
    
    theta_s = calc_theta_s(data_obj.coleta)
    
    raz_bbw_bb = calc_raz_bbw_bb(data_obj.bbw, bb)
    
    kd_qaa = calc_kd_lee_2013(m, theta_s, at, bb, chi, raz_bbw_bb)
    
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
    #   'chute_inicial': a_Lb0_optimized['chute_inicial'],
    #    'h0h1h2': a_Lb0_optimized['h_otimizado'],
        'g0': g0,
        'g1': g1
    }
    # result_dict['a_Lb'] = a_Lb0_optimized['aref']
    result_dict['a_Lb'] = a_Lb0_optimized
    result_dict['bb'] = bb
    result_dict['at'] = at
    result_dict['n'] = n
    result_dict['bbp_lb0'] = bbp_lb0
    result_dict['bbp'] = bbp
    result_dict['Kd_QAA'] = kd_qaa
    result_dict['Kd_ZEU'] = data_obj.kd_zeu
    result_dict['Kd_3m'] = data_obj.kd_3m
    result_dict['Kd_6m'] = data_obj.kd_6m

    os.makedirs('./Results_Espectral', exist_ok = True)
    export_result(result_dict, './Results_Espectral/QAAv5_' + ano + '_' + str(wvl_ref['lb0']) + opt  + '.pickle')
    
    return result_dict