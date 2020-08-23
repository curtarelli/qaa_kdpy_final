'''
Created on Thu Oct 31 15:04:10 2019

@author: Victor
'''
####################################################
##############  Argumentos - V5 S2A  ###############
####################################################

from qaapy_kdpy.definition import create_args
from qaapy_kdpy.file import read_files
from qaapy_kdpy.QAA_models.S2A.Qaa_v6_model import Qaa_v6_model_RUN

args = create_args()

args.acs = {
            'file_name': r'.\00_Dados\Dados_Victor\ACS\aTot_Sus_TRM_s2a.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.Rrs = {
            'file_name': r'.\00_Dados\Dados_Victor\Rrs\Rrs_TRM_s2a.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.aw = {
                'file_name': r'.\00_Dados\Dados_Victor\aw_pope97_s2a.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.bbw = {
                'file_name': r'.\00_Dados\Dados_Victor\bbw_zhang_s2a.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.Kd_ZEU = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_ZEU\Kd_ZEU_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.Kd_3m = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_3m\Kd_3m_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.Kd_6m = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_6m\Kd_6m_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj = read_files(args)

qaa_v6_s2a = Qaa_v6_model_RUN(data_obj, 'TRM_2019', 665)

qaa_v6_s2a = load_data('.\Results_S2A\QAAv6_TRM_2019_S2a_665.pickle')
