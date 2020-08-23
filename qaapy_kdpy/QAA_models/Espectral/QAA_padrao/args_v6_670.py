'''
Created on Thu Oct 31 15:04:10 2019

@author: Victor
'''
####################################################
##############  Argumentos - V5 S2A  ###############
####################################################

from qaapy_kdpy.definition import create_args
from qaapy_kdpy.file import read_files
from qaapy_kdpy.QAA_models.Espectral.Qaa_v6_model import Qaa_v6_model_RUN

args = create_args()

args.acs = {
            'file_name': r'.\00_Dados\Dados_Victor\ACS\aTot_Sus_TRM_PAR.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.Rrs = {
            'file_name': r'.\00_Dados\Dados_Victor\Rrs\Rrs_TRM_Optm_PAR.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.aw = {
                'file_name': r'.\00_Dados\Dados_Victor\aw_pope97_PAR.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.bbw = {
                'file_name': r'.\00_Dados\Dados_Victor\bbw_zhang_PAR.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.coleta = {
                'file_name': r'.\00_Dados\Dados_Victor\Coleta\TRM_2019_Datas_Optm.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.kd_zeu = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_ZEU\Kd_ZEU_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.kd_3m = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_3m\Kd_3m_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.kd_6m = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_6m\Kd_6m_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj = read_files(args)

qaa_v6_670 = Qaa_v6_model_RUN(data_obj, 'TRM_2019', 670)

qaa_v6_670 = load_data('.\Results_Espectral\QAAv6_TRM_2019_670.pickle')
