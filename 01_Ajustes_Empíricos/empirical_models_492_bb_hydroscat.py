'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

RESUMO: Este script é um teste para os passos 2 e 4 do modelo QAA em sua versão
V6, com base em metodologia proposta por:
    
X. Liu et al. (2019) - "Remote Sensing of Secchi Depth in Highly Turbid Lake
Waters and Its Application with MERIS Data"

lb0 = 492 nm

O script aplica um ajuste exponencial para o PASSO 2 entre os dados de a(lb0) e
a razão de três bandas proposta para a equação. Seguido de um ajuste linear para
o PASSO 4 entre o slope do bbp medido como referência e a razão de duas bandas
proposta para a equação.

Esta versão usa o dado de bb medido com uso do HYDROSCAT (6 bandas PAR) para
as regressões lineares do PASSO 4.
'''
###############################################################################
##############################  PACOTES GERAIS  ###############################
###############################################################################
##  Importando os pacotes básicos necessários para o script
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

###############################################################################
#############################  PREPARANDO DADOS  ##############################
###############################################################################
##  Importando pacotes do modelo QAA/Kd necessárias para o teste
from qaapy_kdpy.empirical_definition import create_args
from qaapy_kdpy.empirical_file import read_files
from qaapy_kdpy.conversion import Rrs_to_rrs
from qaapy_kdpy.filter import filter_wvl

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
                'file_name': r'.\00_Dados\Dados_Renato\bb\bb_interp_TRM_2013.xlsx', 
                'sheet_name': 'Sheet1'
            }


args.bbp = {
                'file_name': r'.\00_Dados\Dados_Renato\bb\bbp_interp_TRM_2013.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj = read_files(args)         ##  Leitura dos dados

##  Dados alocados nas suas respectivas variáveis para fácil chamada nas funções
rrs = Rrs_to_rrs(data_obj.Rrs)      ##  Reflectância SR subsuperfície
Rrs = data_obj.Rrs                  ##  Reflectância SR superfície

acs = data_obj.acs                  ##  Absorção ACS
aw = data_obj.aw                    ##  Absorção da água
bbw = data_obj.bbw                  ##  Retroespalhamento da água

bb = data_obj.bb                    ##  Retroespalhamento total
bbp = data_obj.bbp                  ##  Retroespalhamento do material particulado

at = pd.DataFrame()                 ##  Dataframe de absorção total vazio
at = acs + aw.values                ##  Absorção total

##  Parâmetros QAA
g0 = 0.089
g1 = 0.1245

##  Calculo do parâmetro u (STEP 1 QAA)
u = calc_u(g0, g1, rrs)

##  Comprimento de onda inicial (referência - lambda zero)
lb0 = 492

##  Separação dos comprimentos de onda de referência para os cálculos em um
##  dicionário de parâmetros
wvl_ref = {
            'referencia': list(range(400, 751, 1)), ##  Bandas PAR
            'bbp': list(bbp.index),                 ##  Bandas PAR HYDROSCAT/bbp
            'lb0': lb0                              ##  Banda inicial (lambda zero)
            }

##  Variáveis para chamada nos cálculos
wl = pd.Series(wvl_ref['referencia'])
wl_bbp = pd.Series(wvl_ref['bbp'])

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

a_eq = 'y = 0,0159 + 0,4639 * x^-1,3466'        ##  equação lb0 = 492
a_r = 'R² = 0,9672 ; n = 22'                    ##  R² e n lb0 = 492
# a_eq = 'y = 0,0619 + 0,4310 * x^-1,4408'        ##  equação lb0 = 560
# a_r = 'R² = 0,9261 ; n = 22'                    ##  R² e n lb0 = 560

##  Plot da disperção de pontos e curva de ajuste
plot_ajuste(Rrs_ratio,
            at_lb0,
            a_ref,
            a_eq,
            a_r,
            [0, 2.5], [0, 2.5],
            'Rrs(' + str(lb0) + ') / Rrs(665) + Rrs(704) [-]',
            'a(' + str(lb0) + ') [m-1]',
            'Ajuste QAA v6 - Passo 2',
            'Ajuste',
            'Estações',
            1)

##  Cálculos de logarítmos naturais e razão de comprimentos de onda para
bbp_ln = calc_bbp_lin(wl_bbp, bbp, 704)

##  Criando DataFrames vazios para alocar os dados de slope do bbp e R² do ajuste
##  por estação amostral
n_data = pd.DataFrame()
n_r2_data = pd.DataFrame()

##  Função de ajuste para o cálculo do slope do bbp por estação amostral
for i in bbp_ln['ln_bbp']:
    bbp_lb0 = bbp_ln['ln_bbp_lb0'][i]
    bbp_ratio = bbp_ln['ln_ratio']
    fun_n = lambda bbp_ratio, a: (a * bbp_ratio) + bbp_lb0
    ln_bbp = bbp_ln['ln_bbp'][i]
    n_fit = curve_fit(fun_n, bbp_ratio, ln_bbp, [0.5], method= 'lm', maxfev = 100000)[0]
    n_data[i] = n_fit
    
    ##  Parâmetros para cálculo de R²
    n_res = ln_bbp - fun_n(bbp_ratio, n_fit)
    n_ss_res = np.sum(n_res ** 2)
    n_ss_tot = np.sum((ln_bbp - np.mean(ln_bbp)) ** 2)
    n_r2 = [1 - (n_ss_res / n_ss_tot)]
    n_r2_data[i] = n_r2

##  Transpondo os resultados para de colunas para linhas
n_data_t = n_data.transpose()

##  Razão de duas bandas rrs para cálculo do slope do bbp
rrs_ratio = calc_nrrs_ratio(r_ref, 443, 492)
rrs_ratio = rrs_ratio.rename('rrs_ratio', axis = 'columns')

##  Retirando outliers
n_data_t = n_data_t.drop(['P04', 'P12', 'P13'], axis = 0)
rrs_ratio = rrs_ratio.drop(['P04', 'P12', 'P13'], axis = 0)

##  Função de ajuste para o cálculo dos índices da função do slope do bbp
##  [ganho e offset]
fun_nn = lambda rrs_ratio, a, b: (a * rrs_ratio) + b
nn_fit = curve_fit(fun_nn, rrs_ratio, n_data_t[0], [2, -2.4], method = 'lm', maxfev = 100000)[0]

##  Parâmetros para cálculo de R²
nn_res = n_data_t[0] - fun_nn(rrs_ratio, nn_fit[0], nn_fit[1])
nn_ss_res = np.sum(nn_res ** 2)
nn_ss_tot = np.sum((n_data_t[0] - np.mean(n_data_t[0])) ** 2)
nn_r2 = 1 - (nn_ss_res / nn_ss_tot)

nn_x = pd.Series(np.arange(0, 11, 0.01))
nn_ref = fun_nn(nn_x, nn_fit[0], nn_fit[1])
nn_ref = nn_ref.rename(nn_x, axis = 'rows')

nn_eq = 'y = 1,3065x - 1,2805'        ##  equação [Slope do bbp em 704 nm]
nn_r = 'R² = 0,2422 ; n = 19'         ##  R² e n [Slope do bbp em 704 nm]
# nn_eq = 'y = 1,1363x - 1,0422'        ##  equação [Slope do bbp em 560 nm]
# nn_r = 'R² = 0,1955 ; n = 19'         ##  R² e n [Slope do bbp em 560 nm]

##  Plot da disperção de pontos e curva de ajuste
plot_ajuste(rrs_ratio,
            n_data_t,
            nn_ref,
            nn_eq,
            nn_r,
            [1.5, 2.5], [0.8, 2],
            'e^(rrs(443) / rrs(492)) [-]',
            'Slope bbp (n) [-]',
            'Ajuste QAA v6 - Passo 4',
            'Ajuste',
            'Estações',
            1)

###############################################################################
###############################  TESTE R²  ####################################
###############################################################################
from scipy.stats import linregress as lr

slope, intercept, r_value, p_value, str_err = lr(rrs_ratio, n_data_t[0])

r2_teste = r_value ** 2
