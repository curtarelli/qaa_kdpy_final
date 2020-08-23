###############################################################################
#####################          QAA.V6 - model      ############################
###############################################################################
######################## IMPORTS ##############################################
###############################################################################
#TODO: remover pacontes que não serão utilizados

# Import basic modules from python
import pandas as pd
import os

from qaapy_kdpy.conversion import Rrs_to_rrs
from qaapy_kdpy.file import export_result
from qaapy_kdpy.filter import filter_wvl

# QAA - V5 modules
from qaapy_kdpy.QAA_Core.v5.calculation import *

# QAA - V6 modules
from qaapy_kdpy.QAA_Core.v6.calculation import *


def Qaa_v6_model_RUN(data_obj, ano, lb0):
    """
    TODO: WRITE DESCRIPTION AND CREATE A PATERN TO THIS function
    
    """
    ###############################################################################
    ########################### QAA SETTINGS - GERAL ##############################
    ###############################################################################
    data_obj.Rrs = data_obj.Rrs.drop([833, 865, 945, 1374, 1614, 2202])
    data_obj.aw = data_obj.aw.drop([833, 865, 945, 1374, 1614, 2202])
    data_obj.bbw = data_obj.bbw.drop([833, 865, 945, 1374, 1614, 2202])
    
    # Transformação do Rrs
    rrs = Rrs_to_rrs(data_obj.Rrs)
    Rrs = data_obj.Rrs
    
    # Cálculando o U
    g0 = 0.089
    g1 = 0.1245
    
    u = calc_u(g0, g1, rrs)

    ##############################################################################################################
    #ToDo: checar se isso é necessário aqui, talvez criar uma função que padronize isso fora dos OLCI_models.
    # Filtrando os dados
    wvl_ref = {
            'referencia': [443, 492, 560, 665, 704, 741, 783], # Bandas do sensor Sentinel-2A
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779], # Bandas do sensor OLCI
            # 'referencia': [400,413,443,490,510,560,620,665,674,681,709,754,761,764,768,779,865,885,900], # Bandas do sensor OLCI
            # 'referencia': [400,411,443,490,560,620,667,674,681,709,754,761,779,865,555],# Bandas v5  
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
    
    # a_lb0_v6 = calc_alb0_v6(670, R_ref, 443, 490, data_obj.aw); lago = 'Eq. QAA.v6'
    # a_lb0_v6_olci = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 490, data_obj.aw); lago = 'QAA.v6_OLCI'
    a_lb0_v6_s2a = calc_alb0_v6(wvl_ref['lb0'], R_ref, 443, 492, data_obj.aw); lago = 'QAA.v6_S2a'
    
    #### bbp (lambda_0)
    bbp_lb0 = calc_bbp_lb0_v6(u, data_obj.bbw, 665, wvl_ref['lb0'], a_lb0_v6_s2a)
    
    ## Cálculando N
    n = calc_n(r_ref[443], r_ref[560])
    
    # Realizando cálculo do BBP para todos os comprimentos de onda
    bbp = calc_bbp(wl, n, bbp_lb0, wvl_ref['lb0'])
    
    # Realizando cálculo do BB para todos comprimentos de onda
    bb = calc_bb(bbp, data_obj.bbw)
    
    # Realizando cálculo do atotal para todos os comprimentos de onda
    at = calc_a_total(bbp, data_obj.bbw, u)
    
    ###############################################################################
    ############################ FASE 2 - Kd ######################################
    ###############################################################################
    
    m = [0.005, 4.26, 0.52, 10.54]
    
    chi = 0.265
    
    theta_10 = 10
    
    theta_15 = 15
    
    raz_bbw_bb = calc_raz_bbw_bb(data_obj.bbw, bb)
    
    Kd_10 = calc_kd_lee_2013(m, theta_10, at, bb, chi, raz_bbw_bb)
    
    Kd_15 = calc_kd_lee_2013(m, theta_15, at, bb, chi, raz_bbw_bb)

    ###############################################################################
    ########################### Data Extraction ###################################
    ###############################################################################
    result_dict = {}
    result_dict['dados'] = {
        'lb0': wvl_ref['lb0'],
        'r_ref': list(r_ref.keys()),
        'campanhas': ano,
        'equacao': lago,
        'chute_inicial':'não otimizado',
        'h0h1h2': 'não otimizado',
        'g0': g0,
        'g1': g1
    }
    result_dict['a_Lb'] = a_lb0_v6_s2a
    result_dict['bb'] = bb
    result_dict['at'] = at
    result_dict['n'] = n
    result_dict['bbp_lb0'] = bbp_lb0
    result_dict['bbp'] = bbp
    result_dict['Kd_10'] = Kd_10
    result_dict['Kd_15'] = Kd_15
    result_dict['Kd_ZEU'] = data_obj.Kd_ZEU
    result_dict['Kd_3m'] = data_obj.Kd_3m
    result_dict['Kd_6m'] =data_obj.Kd_6m

    os.makedirs('./Results_S2A', exist_ok = True)
    export_result(result_dict, './Results_S2A/QAAv6_' + ano + '_S2a_' + str(wvl_ref['lb0']) + '.pickle')

    return result_dict