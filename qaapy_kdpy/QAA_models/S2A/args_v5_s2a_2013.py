'''
Created on Thu Oct 31 15:04:10 2019

@author: Victor
'''
####################################################
##############  Argumentos - V5 S2A  ###############
####################################################

from qaapy_kdpy.definition import create_args
from qaapy_kdpy.file import read_files, load_data
from qaapy_kdpy.QAA_models.S2A.Qaa_v5_model import Qaa_v5_model_RUN

args = create_args()

args.acs = {
            'file_name': r'.\00_Dados\Dados_Renato\ACS\aACS_TRM2013_s2a.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.Rrs = {
            'file_name': r'.\00_Dados\Dados_Renato\Rrs\Rrs_TRM2013_s2a.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.aw = {
                'file_name': r'.\00_Dados\Dados_Renato\aw_pope97_s2a.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.bbw = {
                'file_name': r'.\00_Dados\Dados_Renato\bbw_zhang_s2a.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.Kd_ZEU = {
                'file_name': r'.\00_Dados\Dados_Renato\Kd_ZEU\Kd_TRM_2013_ZEU.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.Kd_3m = {
                'file_name': r'.\00_Dados\Dados_Renato\Kd_3m\Kd_TRM_2013_3m.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.Kd_6m = {
                'file_name': r'.\00_Dados\Dados_Renato\Kd_6m\Kd_TRM_2013_6m.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj = read_files(args)

lb0 = 560
Opt = '_OptY'

qaa_v5_s2a_560 = Qaa_v5_model_RUN(data_obj, 'TRM_2013', lb0)

qaa_v5_s2a_560 = load_data('.\Results_S2A\QAAv5_TRM_2013_S2a_' + str(lb0) + Opt + '.pickle')
