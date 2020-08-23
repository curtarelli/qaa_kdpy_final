###############################################################################
##############################  PACOTES GERAIS  ###############################
###############################################################################
##  Importando os pacotes básicos necessários para o script
import numpy as np
import pandas as pd

###############################################################################
#############################  PREPARANDO DADOS  ##############################
###############################################################################
##  Importando pacotes do modelo QAA/Kd necessárias para o teste
from qaapy_kdpy.empirical_definition import create_args
from qaapy_kdpy.empirical_file import read_files
from qaapy_kdpy.conversion import Rrs_to_rrs
from qaapy_kdpy.filter import filter_wvl
from qaapy_kdpy.utils.stats import *

##  QAA - Módulos de cálculo da versão V5
from qaapy_kdpy.QAA_Core.v5.calculation import *

##  Criando argumentos e lendo arquivos utilizados no teste
args = create_args()

args.acs = {
            'file_name': r'.\00_Dados\Dados_Renato\ACS\aACS_TRM2013.xlsx', 
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

args.bb = {
                'file_name': r'.\00_Dados\Dados_Renato\bb\bb_TRM_2013.xlsx', 
                'sheet_name': 'Sheet1'
            }


args.bbp = {
                'file_name': r'.\00_Dados\Dados_Renato\bb\bbp_TRM_2013.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj = read_files(args)         ##  Leitura dos dados

##  Morel and Maritonera (2001) Kw:
Kw_490 = 0.01660
Kw_495 = 0.01885
Kw_492 = Kw_490 + ((492 - 490) * ((Kw_490 - Kw_495) / (490 - 495)))

Rrs = data_obj.Rrs                  ##  Reflectância SR superfície

##  Separação dos comprimentos de onda de referência para os cálculos em um
##  dicionário de parâmetros
wvl_ref = {
            'referencia': list(range(400, 751, 1)), ##  Bandas PAR
            'lb0': lb0                              ##  Banda inicial (lambda zero)
            }

##  Variáveis para chamada nos cálculos    
wl = pd.Series(wvl_ref['referencia'])

##  Variáveis de reflectâncias de SR para chamada nos cálculos
r_ref = filter_wvl(rrs, wvl_ref)
R_ref = filter_wvl(Rrs, wvl_ref)

###############################################################################
###############################  OPERACOES  ###################################
###############################################################################
##  Importando pacotes para realização do teste - Empirical Models
from qaapy_kdpy.empirical_calc import *
from qaapy_kdpy.utils.plot_Kd_campo_vs_qaa import plot_ajuste
from scipy.optimize import curve_fit

##  Razão de três bandas Rrs para cálculo de a(lb0)
Rrs_ratio = calc_aRrs_ratio(R_ref, lb0, 665, 704)
Rrs_ratio = Rrs_ratio.rename('Rrs_ratio', axis = 'columns')

##  Função de ajuste para o cálculo dos índices da função a(lb0) [multiplicador e expoente]
fun_a = lambda Rrs_ratio, a, b: aw.loc[lb0].values + (a * np.power(Rrs_ratio, b))
at_lb0 = at.loc[lb0]
a_fit = curve_fit(fun_a, Rrs_ratio, at_lb0, [0.39, 1.14], method= 'lm', maxfev = 100000)[0]

##  Parâmetros para cálculo de R²
a_res = at_lb0 - fun_a(Rrs_ratio, a_fit[0], a_fit[1])
a_ss_res = np.sum(a_res ** 2)
a_ss_tot = np.sum((at_lb0 - np.mean(at_lb0)) ** 2)
a_r2 = 1 - (a_ss_res / a_ss_tot)                    ##  Resultado R² para a(lb0)

##  Curva de ajuste continua para plotagem
a_x = pd.Series(np.arange(0, 11, 0.01))
a_ref = fun_a(a_x, a_fit[0], a_fit[1])
a_ref = a_ref.rename(a_x, axis = 'rows')

# a_eq = 'y = 0,0159 + 0,4639 * x^-1,3466'        ##  equação lb0 = 492
# a_r = 'R² = 0,9672 ; n = 22'                    ##  R² e n lb0 = 492
a_eq = 'y = 0,0619 + 0,4310 * x^-1,4408'      ##  equação lb0 = 560
a_r = 'R² = 0,9261 ; n = 22'                  ##  R² e n lb0 = 560

##  Plot da disperção de pontos e curva de ajuste
plot_ajuste(Rrs_ratio,
            at_lb0,
            a_ref,
            a_eq,
            a_r,
            [0, 3], [0, 1],
            'Rrs(' + str(lb0) + ') / Rrs(665) + Rrs(704) [-]',
            'a(' + str(lb0) + ') [m-1]',
            'Ajuste QAA v6 - Passo 2',
            'Ajuste',
            'Estações',
            1)
