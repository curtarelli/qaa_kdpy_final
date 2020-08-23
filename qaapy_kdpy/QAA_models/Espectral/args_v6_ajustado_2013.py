'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é a aplicação do modelo QAA/Kd v6 aplicado para os dados de 2013 com 
base nos modelos de ajuste realizados para a banda lb0 em 492 nm.

Este modelo é parametrizado para aplicação em imagens Sentinel 2A MSI.

Ajuste do modelo de absorção (Passo 2) utiliza razão de três bandas em 492, 665
e 704 nm.

O ajutes do modelo do slope de bbp (Passo 4) utiliza a razão de duas bandas em 
665 e 704 nm.

O modelo usa a chamada de argumentos para alocar os dados de entrada e então
aplicar no modelo QAA/Kd. Ou entãoum arquivo já processado pode ser aberto.

Após isso alguns calculso estatísticos são realizados para os valores simulados
de Kd nas bandas Sentinel 2A MSI dentro da faixa de 400 a 750 nm.
'''
##  Importando pacotes básicos necessários para o código.
import pandas as pd
import numpy as np

##  Importando pacotes básicos para ajustes e regressões lineares.
from scipy.stats import linregress as lr
from scipy.optimize import curve_fit

##  Importando pacotes básicos para calculos estatísticos.
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse

##  Importando pacotes do modelo QAA/Kd utilizados no código.
from qaapy_kdpy.definition import create_args
from qaapy_kdpy.file import read_files, load_data
from qaapy_kdpy.QAA_models.Espectral.Qaa_v6_model import Qaa_v6_model_RUN
from qaapy_kdpy.empirical_calc import *
from qaapy_kdpy.utils.plot_Kd_campo_vs_qaa import plot_ajuste
from qaapy_kdpy.utils.stats import *

###############################################################################
############################  Argumentos - V6 S2A  ############################
###############################################################################
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
                'file_name': r'.\00_Dados\Dados_Renato\bb\bb_TRM_2013_extra.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.coleta = {
                'file_name': r'.\00_Dados\Dados_Renato\Coleta\TRM_2013_Datas.xlsx', 
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

qaa_v6_fit13_560 = Qaa_v6_model_RUN(data_obj, 'TRM_2013', 560)

qaa_v6_fit13_560 = load_data(r'.\Results_Espectral\QAAv6_TRM_2013_560_Ajustado.pickle')

###############################################################################
########################    Estatísticas (N = 22)   ###########################
###############################################################################
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = qaa_v6_fit13_560['Kd_ZEU'].loc[443]                    ##  Dado de referência (y_true)
kd_pred_443 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[443]       ##  Dados simulados (y_pred)

##  Ajuste da curva de regressão linear para 443 nm.
fit_443 = curve_fit(fun_reg, kd_ref_443, kd_pred_443, pred_ab, method = 'lm', maxfev = 100000)[0]

##  Calculo do coeficiente de determinação R²
res_443 = kd_pred_443 - fun_reg(kd_ref_443, fit_443[0], fit_443[1])
ss_res_443 = np.sum(res_443 ** 2)
ss_tot_443 = np.sum((kd_pred_443 - np.mean(kd_pred_443)) ** 2)
r2_443 = 1 - (ss_res_443 / ss_tot_443)

##  Calculo das estatísticas MAE e MAPE
mae_443 = mae(kd_ref_443, kd_pred_443)
mape_443 = mean_abs_perc_error(kd_ref_443, kd_pred_443)

## Calculo das estatisticas MSE e RMSE
mse_443 = mse(kd_ref_443, kd_pred_443)
rmse_443 = np.sqrt(mse_443)

## Regressão linear usando pacote do scipy para teste
sl_443, inter_443, r_443, p_443, stderr_443 = lr(kd_ref_443, kd_pred_443)

##  Curva de melhor ajuste para plotagem
ref_443 = fun_reg(x_fit, fit_443[0], fit_443[1])
ref_443 = ref_443.rename(x_fit, axis = 'rows')

##  Texto para plotagem com dados dos resultados estatísticos e da regressão
eq_443 = 'y = 0,7067x + 0,4500'
stat_443 = 'R² = 0,1796, MAPE = 0,0%; RMSE = 0,3560; n = 22'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
            [0, 1.6], [0, 1.6],
            'Kd ZEU (443 nm)',
            'Kd QAA (443 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações', 
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 492 nm.
kd_ref_492 = qaa_v6_fit13_560['Kd_ZEU'].loc[492]
kd_pred_492 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[492]

fit_492 = curve_fit(fun_reg, kd_ref_492, kd_pred_492, pred_ab, method = 'lm', maxfev = 100000)[0]

res_492 = kd_pred_492 - fun_reg(kd_ref_492, fit_492[0], fit_492[1])
ss_res_492 = np.sum(res_492 ** 2)
ss_tot_492 = np.sum((kd_pred_492 - np.mean(kd_pred_492)) ** 2)
r2_492 = 1 - (ss_res_492 / ss_tot_492)

mae_492 = mae(kd_ref_492, kd_pred_492)
mape_492 = mean_abs_perc_error(kd_ref_492, kd_pred_492)

mse_492 = mse(kd_ref_492, kd_pred_492)
rmse_492 = np.sqrt(mse_492)

sl_492, inter_492, r_492, p_492, stderr_492 = lr(kd_ref_492, kd_pred_492)

ref_492 = fun_reg(x_fit, fit_492[0], fit_492[1])
ref_492 = ref_492.rename(x_fit, axis = 'rows')

eq_492 = 'y = 1,6916x - 0,3091'
stat_492 = 'R² = 0,9680, MAPE = 28,93%; RMSE = 0,5855; n = 22'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
            [0, 1.5], [0, 1.5],
            'Kd ZEU (492 nm)',
            'Kd QAA (492 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 560 nm
kd_ref_560 = qaa_v6_fit13_560['Kd_ZEU'].loc[560]
kd_pred_560 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[560]

fit_560 = curve_fit(fun_reg, kd_ref_560, kd_pred_560, pred_ab, method = 'lm', maxfev = 100000)[0]

res_560 = kd_pred_560 - fun_reg(kd_ref_560, fit_560[0], fit_560[1])
ss_res_560 = np.sum(res_560** 2)
ss_tot_560 = np.sum((kd_pred_560 - np.mean(kd_pred_560)) ** 2)
r2_560 = 1 - (ss_res_560 / ss_tot_560)

mae_560 = mae(kd_ref_560, kd_pred_560)
mape_560 = mean_abs_perc_error(kd_ref_560, kd_pred_560)

mse_560 = mse(kd_ref_560, kd_pred_560)
rmse_560 = np.sqrt(mse_560)

sl_560, inter_560, r_560, p_560, stderr_560 = lr(kd_ref_560, kd_pred_560)

ref_560 = fun_reg(x_fit, fit_560[0], fit_560[1])
ref_560 = ref_560.rename(x_fit, axis = 'rows')

eq_560 = 'y = 2,0201x - 0,2956'
stat_560 = 'R² = 0,9697, MAPE = 33,06%; RMSE = 0,4763; n = 22'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
            [0, 1], [0, 1],
            'Kd ZEU (560 nm)',
            'Kd QAA (560 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 665 nm
kd_ref_665 = qaa_v6_fit13_560['Kd_ZEU'].loc[665]
kd_pred_665 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[665]

fit_665 = curve_fit(fun_reg, kd_ref_665, kd_pred_665, pred_ab, method = 'lm', maxfev = 100000)[0]

res_665 = kd_pred_665 - fun_reg(kd_ref_665, fit_665[0], fit_665[1])
ss_res_665 = np.sum(res_665 ** 2)
ss_tot_665 = np.sum((kd_pred_665 - np.mean(kd_pred_665)) ** 2)
r2_665 = 1 - (ss_res_665 / ss_tot_665)

mae_665 = mae(kd_ref_665, kd_pred_665)
mape_665 = mean_abs_perc_error(kd_ref_665, kd_pred_665)

mse_665 = mse(kd_ref_665, kd_pred_665)
rmse_665 = np.sqrt(mse_665)

sl_665, inter_665, r_665, p_665, stderr_665 = lr(kd_ref_665, kd_pred_665)

ref_665 = fun_reg(x_fit, fit_665[0], fit_665[1])
ref_665 = ref_665.rename(x_fit, axis = 'rows')

eq_665 = 'y = 2,1719x - 0,6640'
stat_665 = 'R² = 0,8882; MAPE = 32,65%; RMSE = 0,4818; n = 22'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0, 2], [0, 2],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = qaa_v6_fit13_560['Kd_ZEU'].loc[704]
kd_pred_704 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[704]

fit_704 = curve_fit(fun_reg, kd_ref_704, kd_pred_704, pred_ab, method = 'lm', maxfev = 100000)[0]

res_704 = kd_pred_704 - fun_reg(kd_ref_704, fit_704[0], fit_704[1])
ss_res_704 = np.sum(res_704 ** 2)
ss_tot_704 = np.sum((kd_pred_704 - np.mean(kd_pred_704)) ** 2)
r2_704 = 1 - (ss_res_704 / ss_tot_704)

mae_704 = mae(kd_ref_704, kd_pred_704)
mape_704 = mean_abs_perc_error(kd_ref_704, kd_pred_704)

mse_704 = mse(kd_ref_704, kd_pred_704)
rmse_704 = np.sqrt(mse_704)

sl_704, inter_704, r_704, p_704, stderr_704 = lr(kd_ref_704, kd_pred_704)

ref_704 = fun_reg(x_fit, fit_704[0], fit_704[1])
ref_704 = ref_704.rename(x_fit, axis = 'rows')

eq_704 = 'y = 1,8832x - 0,4573'
stat_704 = 'R² = 0,8407, MAPE = 40,01%; RMSE = 0,5424; n = 22'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0, 2], [0, 2],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
########################    Estatísticas (N = 20)   ###########################
###############################################################################
##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = qaa_v6_fit13_560['Kd_ZEU'].loc[443].drop(['P25', 'P26'])               ##  Dado de referência (y_true)
kd_pred_443 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[443].drop(['P25', 'P26'])  ##  Dados simulados (y_pred)

##  Ajuste da curva de regressão linear para 443 nm.
fit_443 = curve_fit(fun_reg, kd_ref_443, kd_pred_443, pred_ab, method = 'lm', maxfev = 100000)[0]

##  Calculo do coeficiente de determinação R²
res_443 = kd_pred_443 - fun_reg(kd_ref_443, fit_443[0], fit_443[1])
ss_res_443 = np.sum(res_443 ** 2)
ss_tot_443 = np.sum((kd_pred_443 - np.mean(kd_pred_443)) ** 2)
r2_443 = 1 - (ss_res_443 / ss_tot_443)

##  Calculo das estatísticas MAE e MAPE
mae_443 = mae(kd_ref_443, kd_pred_443)
mape_443 = mean_abs_perc_error(kd_ref_443, kd_pred_443)

## Calculo das estatisticas MSE e RMSE
mse_443 = mse(kd_ref_443, kd_pred_443)
rmse_443 = np.sqrt(mse_443)

## Regressão linear usando pacote do scipy para teste
sl_443, inter_443, r_443, p_443, stderr_443 = lr(kd_ref_443, kd_pred_443)

##  Curva de melhor ajuste para plotagem
ref_443 = fun_reg(x_fit, fit_443[0], fit_443[1])
ref_443 = ref_443.rename(x_fit, axis = 'rows')

##  Texto para plotagem com dados dos resultados estatísticos e da regressão
eq_443 = 'y = 0,9999x + 0,2710'
stat_443 = 'R² = 0,8111, MAPE = 30,29%; RMSE = 0,3110; n = 20'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
            [0, 1.8], [0, 1.8],
            'Kd ZEU (443 nm)',
            'Kd QAA (443 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 492 nm
kd_ref_492 = qaa_v6_fit13_560['Kd_ZEU'].loc[492].drop(['P25', 'P26'])
kd_pred_492 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[492].drop(['P25', 'P26'])

fit_492 = curve_fit(fun_reg, kd_ref_492, kd_pred_492, pred_ab, method = 'lm', maxfev = 100000)[0]

res_492 = kd_pred_492 - fun_reg(kd_ref_492, fit_492[0], fit_492[1])
ss_res_492 = np.sum(res_492 ** 2)
ss_tot_492 = np.sum((kd_pred_492 - np.mean(kd_pred_492)) ** 2)
r2_492 = 1 - (ss_res_492 / ss_tot_492)

mae_492 = mae(kd_ref_492, kd_pred_492)
mape_492 = mean_abs_perc_error(kd_ref_492, kd_pred_492)

mse_492 = mse(kd_ref_492, kd_pred_492)
rmse_492 = np.sqrt(mse_492)

sl_492, inter_492, r_492, p_492, stderr_492 = lr(kd_ref_492, kd_pred_492)

ref_492 = fun_reg(x_fit, fit_492[0], fit_492[1])
ref_492 = ref_492.rename(x_fit, axis = 'rows')

eq_492 = 'y = 1,0110x + 0,1531'
stat_492 = 'R² = 0,7809, MAPE = 25,83%; RMSE = 0,2024; n = 20'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
            [0, 1.5], [0, 1.5],
            'Kd ZEU (492 nm)',
            'Kd QAA (492 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 560 nm
kd_ref_560 = qaa_v6_fit13_560['Kd_ZEU'].loc[560].drop(['P25', 'P26'])
kd_pred_560 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[560].drop(['P25', 'P26'])

fit_560 = curve_fit(fun_reg, kd_ref_560, kd_pred_560, pred_ab, method = 'lm', maxfev = 100000)[0]

res_560 = kd_pred_560 - fun_reg(kd_ref_560, fit_560[0], fit_560[1])
ss_res_560 = np.sum(res_560** 2)
ss_tot_560 = np.sum((kd_pred_560 - np.mean(kd_pred_560)) ** 2)
r2_560 = 1 - (ss_res_560 / ss_tot_560)

mae_560 = mae(kd_ref_560, kd_pred_560)
mape_560 = mean_abs_perc_error(kd_ref_560, kd_pred_560)

mse_560 = mse(kd_ref_560, kd_pred_560)
rmse_560 = np.sqrt(mse_560)

sl_560, inter_560, r_560, p_560, stderr_560 = lr(kd_ref_560, kd_pred_560)

ref_560 = fun_reg(x_fit, fit_560[0], fit_560[1])
ref_560 = ref_560.rename(x_fit, axis = 'rows')

eq_560 = 'y = 1,2216x + 0,0155'
stat_560 = 'R² = 0,6415, MAPE = 26,96%; RMSE = 0,1451; n = 20'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
            [0, 1], [0, 1],
            'Kd ZEU (560 nm)',
            'Kd QAA (560 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 665 nm
kd_ref_665 = qaa_v6_fit13_560['Kd_ZEU'].loc[665].drop(['P25', 'P26'])
kd_pred_665 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[665].drop(['P25', 'P26'])

fit_665 = curve_fit(fun_reg, kd_ref_665, kd_pred_665, pred_ab, method = 'lm', maxfev = 100000)[0]

res_665 = kd_pred_665 - fun_reg(kd_ref_665, fit_665[0], fit_665[1])
ss_res_665 = np.sum(res_665 ** 2)
ss_tot_665 = np.sum((kd_pred_665 - np.mean(kd_pred_665)) ** 2)
r2_665 = 1 - (ss_res_665 / ss_tot_665)

mae_665 = mae(kd_ref_665, kd_pred_665)
mape_665 = mean_abs_perc_error(kd_ref_665, kd_pred_665)

mse_665 = mse(kd_ref_665, kd_pred_665)
rmse_665 = np.sqrt(mse_665)

sl_665, inter_665, r_665, p_665, stderr_665 = lr(kd_ref_665, kd_pred_665)

ref_665 = fun_reg(x_fit, fit_665[0], fit_665[1])
ref_665 = ref_665.rename(x_fit, axis = 'rows')

eq_665 = 'y = 0,9118x + 0,2599'
stat_665 = 'R² = 0,6055; MAPE = 26,47%; RMSE = 0,2148; n = 20'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0, 1.5], [0, 1.5],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = qaa_v6_fit13_560['Kd_ZEU'].loc[704].drop(['P25', 'P26'])
kd_pred_704 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].loc[704].drop(['P25', 'P26'])

fit_704 = curve_fit(fun_reg, kd_ref_704, kd_pred_704, pred_ab, method = 'lm', maxfev = 100000)[0]

res_704 = kd_pred_704 - fun_reg(kd_ref_704, fit_704[0], fit_704[1])
ss_res_704 = np.sum(res_704 ** 2)
ss_tot_704 = np.sum((kd_pred_704 - np.mean(kd_pred_704)) ** 2)
r2_704 = 1 - (ss_res_704 / ss_tot_704)

mae_704 = mae(kd_ref_704, kd_pred_704)
mape_704 = mean_abs_perc_error(kd_ref_704, kd_pred_704)

mse_704 = mse(kd_ref_704, kd_pred_704)
rmse_704 = np.sqrt(mse_704)

sl_704, inter_704, r_704, p_704, stderr_704 = lr(kd_ref_704, kd_pred_704)

ref_704 = fun_reg(x_fit, fit_704[0], fit_704[1])
ref_704 = ref_704.rename(x_fit, axis = 'rows')

eq_704 = 'y = 0,6606x + 0,6029'
stat_704 = 'R² = 0,5527, MAPE = 35,45%; RMSE = 0,3173; n = 20'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0, 2], [0, 2],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ref vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)
