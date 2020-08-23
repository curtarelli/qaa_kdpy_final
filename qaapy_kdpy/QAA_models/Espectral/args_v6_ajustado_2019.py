'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é a aplicação do modelo QAA/Kd v6 aplicado para os dados de 2019 com 
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
import matplotlib.pyplot as plt
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
from qaapy_kdpy.utils.stats import mean_abs_perc_error

###############################################################################
############################  Argumentos - V6 S2A  ############################
###############################################################################
args = create_args()

args.acs = {
            'file_name': r'.\00_Dados\Dados_Victor\ACS\aTot_Sus_TRM_PAR.xlsx', 
            'sheet_name': 'Sheet1'
            }

args.Rrs = {
            'file_name': r'.\00_Dados\Dados_Victor\Rrs\Rrs_TRM_PAR.xlsx', 
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

args.bb = {
                'file_name': r'.\00_Dados\Dados_Victor\bb\bb_TRM_2019_extra.xlsx', 
                'sheet_name': 'Sheet1'
            }

args.coleta = {
                'file_name': r'.\00_Dados\Dados_Victor\Coleta\TRM_2019_Datas.xlsx', 
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

##  Rodando o modelo QAA/Kd com base nos arquivos chamados pelos argumentos.
qaa_v6_fit19_560 = Qaa_v6_model_RUN(data_obj, 'TRM_2019', 560)

##  Abrindo umarquivo contendo dicionário de dados resultados da aplicação do QAA/Kd.
qaa_v6_fit19_560 = load_data(r'.\Results_Espectral\QAAv6_TRM_2019_560_Ajustado.pickle')

plt.figure(figsize = (10, 5))

plt.subplot(121)
plt.plot(qaa_v6_fit13_560['Kd_ZEU']['P15'], 'r-')
plt.plot(qaa_v6_fit13_560['kd_qaa_560_acs_492']['P15'], 'b-')
plt.xlabel('Comprimento de Onda [nm]', fontsize = 14)
plt.ylabel('Kd [m ^ -1]', fontsize = 14)
plt.xlim([400, 700])
plt.ylim([0, 2])
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
plt.title('2013 - P15', fontsize = 14)
plt.show()

plt.subplot(122)
plt.plot(qaa_v6_fit19_560['Kd_ZEU']['TRM09'], 'r-')
plt.plot(qaa_v6_fit19_560['kd_qaa_560_acs_492']['TRM09'], 'b-')
plt.xlabel('Comprimento de Onda [nm]', fontsize = 14)
plt.xlim([400, 700])
plt.ylim([0, 2])
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
plt.title('2019 - TRM09', fontsize = 14)
plt.legend(['Kd in-situ', 'Kd QAA'], loc = 1, fontsize = 14) 
plt.show()

###############################################################################
########################    Estatísticas (N = 20)   ###########################
###############################################################################
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = qaa_v6_fit19_560['Kd_ZEU'].loc[443]                    ##  Dado de referência (y_true)
kd_pred_443 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[443]       ##  Dados simulados (y_pred)

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
eq_443 = 'y = 0,7175x + 0,6611'
stat_443 = 'R² = 0,1784, MAPE = 50,00%; RMSE = 0,4542; n = 20'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
            [0.5, 1.5], [0.5, 1.5],
            'Kd ZEU (443 nm)',
            'Kd QAA (443 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 492 nm.
kd_ref_492 = qaa_v6_fit19_560['Kd_ZEU'].loc[492]
kd_pred_492 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[492]

##  Ajuste da curva de regressão linear para 492 nm.
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

eq_492 = 'y = 1,2146x - 0,0380'
stat_492 = 'R² = 0,7738, MAPE = 15,44%; RMSE = 0,1052; n = 20'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
            [0, 1], [0, 1],
            'Kd ZEU (492 nm)',
            'Kd QAA (492 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 560 nm
kd_ref_560 = qaa_v6_fit19_560['Kd_ZEU'].loc[560]
kd_pred_560 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[560]

##  Ajuste da curva de regressão linear para 560 nm.
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

eq_560 = 'y = 1,4833x - 0,1744'
stat_560 = 'R² = 0,7957, MAPE = 7,96%; RMSE = 0,0421; n = 20'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
            [0, 1], [0, 1],
            'Kd ZEU (560 nm)',
            'Kd QAA (560 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 665 nm
kd_ref_665 = qaa_v6_fit19_560['Kd_ZEU'].loc[665]
kd_pred_665 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[665]

##  Ajuste da curva de regressão linear para 665 nm.
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

eq_665 = 'y = 0,9680x - 0,0744'
stat_665 = 'R² = 0,2691; MAPE = 9,12%; RMSE = 0,0697; n = 20'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0.5, 1], [0.5, 1],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = qaa_v6_fit19_560['Kd_ZEU'].loc[704]
kd_pred_704 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[704]

##  Ajuste da curva de regressão linear para 704 nm.
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

eq_704 = 'y = 0,5363x + 0,5681'
stat_704 = 'R² = 0,2072; MAPE = 22,38%; RMSE = 0,1936; n = 20'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0.5, 1.5], [0.5, 1.5],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(221)
plt.text(0.1, 1.3, 'a)', fontsize = 12)
plt.text(0.75, 0.25, eq_443)
plt.text(0.75, 0.15, 'R² = 0,1784; MAPE = 50,00%')
plt.text(0.75, 0.05, 'RMSE = 0,4542; n = 20')
plt.plot(ref_443, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (443 nm)', fontsize = 14)
plt.ylabel('Kd QAA (443 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
#plt.title('443 nm', fontsize = 14)
plt.show()

plt.subplot(222)
plt.text(0.5, 1.3, 'b)', fontsize = 12)
plt.text(0.75, 0.25, eq_492)
plt.text(0.75, 0.15, 'R² = 0,7738; MAPE = 15,44%')
plt.text(0.75, 0.05, 'RMSE = 0,1052; n = 20')
plt.plot(ref_492, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (492 nm)', fontsize = 14)
plt.ylabel('Kd QAA (492 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
#plt.title('492 nm', fontsize = 14)
plt.legend(['Fit', '1:1', '2019'], loc = 2, fontsize = 12) 
plt.show()

plt.subplot(223)
plt.text(0.1, 1.3, 'c)', fontsize = 12)
plt.text(0.75, 0.25, eq_560)
plt.text(0.75, 0.15, 'R² = 0,7957; MAPE = 7,96%')
plt.text(0.75, 0.05, 'RMSE = 0,0421; n = 20')
plt.plot(ref_560, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (560 nm)', fontsize = 14)
plt.ylabel('Kd QAA (560 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
#plt.title('560 nm', fontsize = 14)
plt.show()

plt.subplot(224)
plt.text(0.1, 1.3, 'd)', fontsize = 12)
plt.text(0.75, 0.25, eq_665)
plt.text(0.75, 0.15, 'R² = 0,2691; MAPE = 9,12%')
plt.text(0.75, 0.05, 'RMSE = 0,0697; n = 20')
plt.plot(ref_665, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (665 nm)', fontsize = 14)
plt.ylabel('Kd QAA (665 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
#plt.title('665 nm', fontsize = 14)
plt.show()

###############################################################################
########################    Estatísticas (N = 16)   ###########################
###############################################################################
##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = qaa_v6_fit19_560['Kd_ZEU'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])  ##  Dado de referência (y_true)
kd_pred_443 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])  ##  Dados simulados (y_pred)

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
eq_443 = 'y = 0,9007x + 0,4646'
stat_443 = 'R² = 0,4304, MAPE = 45,17%; RMSE = 0,4001; n = 16'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
            [0.5, 1.5], [0.5, 1.5],
            'Kd ZEU (443 nm)',
            'Kd QAA (443 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 492 nm.
kd_ref_492 = qaa_v6_fit19_560['Kd_ZEU'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_492 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

##  Ajuste da curva de regressão linear para 492 nm.
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

eq_492 = 'y = 1,2397x - 0,0511'
stat_492 = 'R² = 0,8292, MAPE = 15,73%; RMSE = 0,1064; n = 16'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
            [0, 1], [0, 1],
            'Kd ZEU (492 nm)',
            'Kd QAA (492 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 560 nm
kd_ref_560 = qaa_v6_fit19_560['Kd_ZEU'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_560 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

##  Ajuste da curva de regressão linear para 560 nm.
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

eq_560 = 'y = 1,5397x - 0,1886'
stat_560 = 'R² = 0,8972, MAPE = 7,50%; RMSE = 0,0401; n = 16'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
            [0, 1], [0, 1],
            'Kd ZEU (560 nm)',
            'Kd QAA (560 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 665 nm
kd_ref_665 = qaa_v6_fit19_560['Kd_ZEU'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_665 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

##  Ajuste da curva de regressão linear para 665 nm.
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

eq_665 = 'y = 1,3349x - 0,1620'
stat_665 = 'R² = 0,4772; MAPE = 9,23%; RMSE = 0,0718; n = 16'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0.5, 1], [0.5, 1],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = qaa_v6_fit19_560['Kd_ZEU'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_704 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

##  Ajuste da curva de regressão linear para 704 nm.
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

eq_704 = 'y = 0,6809x + 0,4478'
stat_704 = 'R² = 0,3003, MAPE = 22,67%; RMSE = 0,1941; n = 16'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0.5, 1.5], [0.5, 1.5],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
############################    Estatísticas 2    #############################
###############################################################################

kd_ref_acs_443 = qaa_v6_fit19_560['kd_acs'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_bb_443 = qaa_v6_fit19_560['kd_bb'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_qaa_492_acs_492_443 = qaa_v6_fit19_560['kd_qaa_492_acs_492'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_acs_492_443 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_hydro_704_443 = qaa_v6_fit19_560['kd_qaa_560_hydro_704'].loc[443].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_ref_acs_492 = qaa_v6_fit19_560['kd_acs'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_bb_492 = qaa_v6_fit19_560['kd_bb'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_qaa_492_acs_492_492 = qaa_v6_fit19_560['kd_qaa_492_acs_492'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_acs_492_492 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_hydro_704_492 = qaa_v6_fit19_560['kd_qaa_560_hydro_704'].loc[492].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_ref_acs_560 = qaa_v6_fit19_560['kd_acs'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_bb_560 = qaa_v6_fit19_560['kd_bb'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_qaa_492_acs_492_560 = qaa_v6_fit19_560['kd_qaa_492_acs_492'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_acs_492_560 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_hydro_704_560 = qaa_v6_fit19_560['kd_qaa_560_hydro_704'].loc[560].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_ref_acs_665 = qaa_v6_fit19_560['kd_acs'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_bb_665 = qaa_v6_fit19_560['kd_bb'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_qaa_492_acs_492_665 = qaa_v6_fit19_560['kd_qaa_492_acs_492'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_acs_492_665 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_hydro_704_665 = qaa_v6_fit19_560['kd_qaa_560_hydro_704'].loc[665].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_ref_acs_704 = qaa_v6_fit19_560['kd_acs'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_bb_704 = qaa_v6_fit19_560['kd_bb'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

kd_qaa_492_acs_492_704 = qaa_v6_fit19_560['kd_qaa_492_acs_492'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_acs_492_704 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_qaa_560_hydro_704_704 = qaa_v6_fit19_560['kd_qaa_560_hydro_704'].loc[704].drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

mape_492_acs_492_acs_443 = mean_abs_perc_error(kd_ref_acs_443, kd_qaa_492_acs_492_443)
mape_492_acs_492_bb_443 = mean_abs_perc_error(kd_ref_bb_443, kd_qaa_492_acs_492_443)
mape_560_acs_492_acs_443 = mean_abs_perc_error(kd_ref_acs_443, kd_qaa_560_acs_492_443)
mape_560_acs_492_bb_443 = mean_abs_perc_error(kd_ref_bb_443, kd_qaa_560_acs_492_443)
mape_560_hydro_704_acs_443 = mean_abs_perc_error(kd_ref_acs_443, kd_qaa_560_hydro_704_443)
mape_560_hydro_704_bb_443 = mean_abs_perc_error(kd_ref_bb_443, kd_qaa_560_hydro_704_443)

mape_492_acs_492_acs_492 = mean_abs_perc_error(kd_ref_acs_492, kd_qaa_492_acs_492_492)
mape_492_acs_492_bb_492 = mean_abs_perc_error(kd_ref_bb_492, kd_qaa_492_acs_492_492)
mape_560_acs_492_acs_492 = mean_abs_perc_error(kd_ref_acs_492, kd_qaa_560_acs_492_492)
mape_560_acs_492_bb_492 = mean_abs_perc_error(kd_ref_bb_492, kd_qaa_560_acs_492_492)
mape_560_hydro_704_acs_492 = mean_abs_perc_error(kd_ref_acs_492, kd_qaa_560_hydro_704_492)
mape_560_hydro_704_bb_492 = mean_abs_perc_error(kd_ref_bb_492, kd_qaa_560_hydro_704_492)

mape_492_acs_492_acs_560 = mean_abs_perc_error(kd_ref_acs_560, kd_qaa_492_acs_492_560)
mape_492_acs_492_bb_560 = mean_abs_perc_error(kd_ref_bb_560, kd_qaa_492_acs_492_560)
mape_560_acs_492_acs_560 = mean_abs_perc_error(kd_ref_acs_560, kd_qaa_560_acs_492_560)
mape_560_acs_492_bb_560 = mean_abs_perc_error(kd_ref_bb_560, kd_qaa_560_acs_492_560)
mape_560_hydro_704_acs_560 = mean_abs_perc_error(kd_ref_acs_560, kd_qaa_560_hydro_704_560)
mape_560_hydro_704_bb_560 = mean_abs_perc_error(kd_ref_bb_560, kd_qaa_560_hydro_704_560)

mape_492_acs_492_acs_665 = mean_abs_perc_error(kd_ref_acs_665, kd_qaa_492_acs_492_665)
mape_492_acs_492_bb_665 = mean_abs_perc_error(kd_ref_bb_665, kd_qaa_492_acs_492_665)
mape_560_acs_492_acs_665 = mean_abs_perc_error(kd_ref_acs_665, kd_qaa_560_acs_492_665)
mape_560_acs_492_bb_665 = mean_abs_perc_error(kd_ref_bb_665, kd_qaa_560_acs_492_443)
mape_560_hydro_704_acs_665 = mean_abs_perc_error(kd_ref_acs_665, kd_qaa_560_hydro_704_665)
mape_560_hydro_704_bb_665 = mean_abs_perc_error(kd_ref_bb_665, kd_qaa_560_hydro_704_665)

mape_492_acs_492_acs_704 = mean_abs_perc_error(kd_ref_acs_704, kd_qaa_492_acs_492_704)
mape_492_acs_492_bb_704 = mean_abs_perc_error(kd_ref_bb_704, kd_qaa_492_acs_492_704)
mape_560_acs_492_acs_704 = mean_abs_perc_error(kd_ref_acs_704, kd_qaa_560_acs_492_704)
mape_560_acs_492_bb_704 = mean_abs_perc_error(kd_ref_bb_704, kd_qaa_560_acs_492_704)
mape_560_hydro_704_acs_704 = mean_abs_perc_error(kd_ref_acs_704, kd_qaa_560_hydro_704_704)
mape_560_hydro_704_bb_704 = mean_abs_perc_error(kd_ref_bb_704, kd_qaa_560_hydro_704_704)

kd_ref_443 = kd_ref_443.drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_492 = kd_ref_492.drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_560 = kd_ref_560.drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_665 = kd_ref_665.drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_ref_704 = kd_ref_704.drop(['TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

mape_ref_acs_443 = mean_abs_perc_error(kd_ref_443, kd_ref_acs_443)
mape_ref_acs_492 = mean_abs_perc_error(kd_ref_492, kd_ref_acs_492)
mape_ref_acs_560 = mean_abs_perc_error(kd_ref_560, kd_ref_acs_560)
mape_ref_acs_665 = mean_abs_perc_error(kd_ref_665, kd_ref_acs_665)
mape_ref_acs_704 = mean_abs_perc_error(kd_ref_704, kd_ref_acs_704)

###############################################################################
###########################  APLICAÇÃO NAS IMAGENS  ###########################
###############################################################################

###############################################################################
###############    Estatísticas Algoritmo Global (N = 100)   ##################
###############################################################################
##  Regressão, plots e estatísticas para algoritmo global
kd_ref_global = pd.concat([kd_ref_443.copy(), kd_ref_492.copy(), kd_ref_560.copy(), kd_ref_665.copy(), kd_ref_704.copy()], axis = 0)
kd_pred_global = pd.concat([kd_pred_443.copy(), kd_pred_492.copy(), kd_pred_560.copy(), kd_pred_665.copy(), kd_pred_704.copy()], axis = 0)

##  Ajuste da curva de regressão linear para 443 nm.
fit_global = curve_fit(fun_reg, kd_ref_global, kd_pred_global, pred_ab, method = 'lm', maxfev = 100000)[0]

##  Calculo do coeficiente de determinação R²
res_global = kd_pred_global - fun_reg(kd_ref_global, fit_global[0], fit_global[1])
ss_res_global = np.sum(res_global ** 2)
ss_tot_global = np.sum((kd_pred_global - np.mean(kd_pred_global)) ** 2)
r2_global = 1 - (ss_res_global / ss_tot_global)

##  Calculo das estatísticas MAE e MAPE
mae_global = mae(kd_ref_global, kd_pred_global)
mape_global = mean_abs_perc_error(kd_ref_global, kd_pred_global)

## Calculo das estatisticas MSE e RMSE
mse_global = mse(kd_ref_global, kd_pred_global)
rmse_global = np.sqrt(mse_global)

## Regressão linear usando pacote do scipy para teste
sl_global, inter_global, r_global, p_global, stderr_global = lr(kd_ref_global, kd_pred_global)

##  Curva de melhor ajuste para plotagem
ref_global = fun_reg(x_fit, fit_global[0], fit_global[1])
ref_global = ref_global.rename(x_fit, axis = 'rows')

##  Texto para plotagem com dados dos resultados estatísticos e da regressão
eq_global = 'y = 1.5621x - 0.2226'
stat_global = 'R² = 0.8281; MAPE = 20,98%; RMSE = 0.2287 m^-1; n = 100'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_global,
            kd_pred_global,
            ref_global,
            eq_global,
            stat_global,
            [0, 1.5], [0, 1.5],
            'Kd - Field (GLOBAL MSI)',
            'Kd - QAA (GLOBAL MSI)',
            'Kd Field vs. Kd QAA - Global MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################

kd_ref19_443 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[443]
kd_ref19_492 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[492]
kd_ref19_560 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[560]
kd_ref19_665 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[665]
kd_ref19_704 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[704]

kd_pred19_443 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[443]
kd_pred19_492 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[492]
kd_pred19_560 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[560]
kd_pred19_665 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[665]
kd_pred19_704 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[704]

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(111)
plt.text(1, 0.6, eq_global, fontsize = 10)
plt.text(1, 0.4, 'R² = 0.8281; MAPE = 20.98%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0.2287 m^-1; n = 100', fontsize = 10)
plt.plot(ref_global, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref19_443, kd_pred19_443, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref19_492, kd_pred19_492, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref19_560, kd_pred19_560, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref19_665, kd_pred19_665, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.scatter(kd_ref19_704, kd_pred19_704, marker = '^', facecolors = 'none', edgecolors = '#FF1493')
plt.xlabel('Kd - Field (GLOBAL)', fontsize = 14)
plt.ylabel('Kd - QAA (GLOBAL)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.5, 1, 1.5])
plt.yticks([0, 0.5, 1, 1.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.legend(['Fit', '1:1', '443 nm (2019)', '492 nm', '560 nm', '665 nm', '704 nm'], loc = 0, fontsize = 12) 
plt.show()
