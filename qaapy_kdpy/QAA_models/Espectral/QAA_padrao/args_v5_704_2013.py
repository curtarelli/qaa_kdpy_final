'''
Created on Thu Oct 31 15:04:10 2019

@author: Victor
'''
####################################################
##############  Argumentos - V5 S2A  ###############
####################################################

from qaapy_kdpy.definition import create_args
from qaapy_kdpy.file import read_files, load_data
from qaapy_kdpy.QAA_models.Espectral.Qaa_v5_model import Qaa_v5_model_RUN

args = create_args()

args.acs = {
            'file_name': r'.\00_Dados\Dados_Renato\ACS\aACS_TRM2013_PAR-2.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.Rrs = {
            'file_name': r'.\00_Dados\Dados_Renato\Rrs\Rrs_TRM2013_PAR.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.aw = {
                'file_name': r'.\00_Dados\Dados_Renato\aw_pope97_PAR.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.bbw = {
                'file_name': r'.\00_Dados\Dados_Renato\bbw_zhang_PAR.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.coleta = {
                'file_name': r'.\00_Dados\Dados_Renato\Coleta\TRM_2013_Datas-2.xlsx', 
                'sheet_name': 'Sheet1'
            }


args.kd_zeu = {
                'file_name': r'.\00_Dados\Dados_Renato\Kd_ZEU\Kd_TRM_2013_ZEU.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.kd_3m = {
                'file_name': r'.\00_Dados\Dados_Renato\Kd_3m\Kd_TRM_2013_3m.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.kd_6m = {
                'file_name': r'.\00_Dados\Dados_Renato\Kd_6m\Kd_TRM_2013_6m.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj = read_files(args)

lb0 = 704   
Opt = '_OptY'
OptAw = '_OptAw'

qaa_v5_704_2 = Qaa_v5_model_RUN(data_obj, 'TRM_2013', 704)

qaa_v5_704_2 = load_data('.\Results_Espectral\QAAv5_TRM_2013_' + str(lb0) + Opt + '.pickle')
qaa_v5_704_2aw = load_data('.\Results_Espectral\QAAv5_TRM_2013_' + str(lb0) + OptAw + '.pickle')
