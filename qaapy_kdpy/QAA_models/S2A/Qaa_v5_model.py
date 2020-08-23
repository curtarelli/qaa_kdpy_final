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
    """
    TODO: WRITE DESCRIPTION AND CREATE A PATERN TO THIS function
    
    """
    ###############################################################################
    ########################### QAA SETTINGS - GENERAL ############################
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
    
    ###############################################################################
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
    ############################ FASE 1 - QAA.V5 - at #############################
    ###############################################################################
    
    p0 = [-1.146, -1.336,- 0.469] # Valores iniciais
    
    # eq_lee_v5_OLCI = np.log10( (r_ref[443] + r_ref[490]) / ( r_ref['r_lb0'] + 5 * (r_ref[665]/r_ref[510] * r_ref[665]))); lago = 'QAA.v5_OLCI'
    eq_lee_v5_s2a = np.log10( (r_ref[443] + r_ref[492]) / ( r_ref['lb0'] + 5 * (r_ref[665]/r_ref[492] * r_ref[665]))); lago = 'QAA.v5_S2a'

    #### ATENÇÃO : a agua é adicionada ao acs dentro da função "optimizer_aLb0" ##########
    a_Lb0_optimized = optimizer_aLb0(eq_lee_v5_s2a, wvl_ref['lb0'], p0, data_obj.aw, data_obj.acs, NOToptimize = True); opt = '_OptN'
    # a_Lb0_optimized = optimizer_aLb0(eq_lee_v5_s2a, wvl_ref['lb0'], p0, data_obj.aw, data_obj.acs, NOToptimize = False); opt = '_OptY'
    
    #### bbp (lambda_0)
    bbp_lb0 = calc_bbp_lb0(u, data_obj.bbw, 560, wvl_ref['lb0'], a_Lb0_optimized['aref'])
    
    ## Cálculando N
    n = calc_n(r_ref[443], r_ref[560])
    
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
        'chute_inicial': a_Lb0_optimized['chute_inicial'],
        'h0h1h2': a_Lb0_optimized['h_otimizado'],
        'g0': g0,
        'g1': g1
    }
    result_dict['a_Lb'] = a_Lb0_optimized['aref']
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
    export_result(result_dict, './Results_S2A/QAAv5_' + ano + '_S2a_' + str(wvl_ref['lb0']) + opt  + '.pickle')
    
    return result_dict