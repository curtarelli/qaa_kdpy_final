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
from qaapy_kdpy.QAA_models.Espectral.Qaa_v6_model import Qaa_v6_model_RUN_original
from qaapy_kdpy.empirical_calc import *
from qaapy_kdpy.utils.plot_Kd_campo_vs_qaa import plot_ajuste, plot_kd_nxm
from qaapy_kdpy.utils.stats import *

###############################################################################
############################  Argumentos - V6 S2A  ############################
###############################################################################
args19 = create_args()

args19.acs = {
            'file_name': r'.\00_Dados\Dados_Victor\ACS\aTot_Sus_TRM_PAR.xlsx', 
            'sheet_name': 'Sheet1'
            }

args19.Rrs = {
            'file_name': r'.\00_Dados\Dados_Victor\Rrs\Rrs_TRM_PAR.xlsx', 
            'sheet_name': 'Sheet1'
            }

args19.aw = {
                'file_name': r'.\00_Dados\Dados_Victor\aw_pope97_PAR.xlsx', 
                'sheet_name': 'Sheet1'
            }

args19.bbw = {
                'file_name': r'.\00_Dados\Dados_Victor\bbw_zhang_PAR.xlsx', 
                'sheet_name': 'Sheet1'
            }

args19.bb = {
                'file_name': r'.\00_Dados\Dados_Victor\bb\bb_TRM_2019_extra.xlsx', 
                'sheet_name': 'Sheet1'
            }

args19.coleta = {
                'file_name': r'.\00_Dados\Dados_Victor\Coleta\TRM_2019_Datas.xlsx', 
                'sheet_name': 'Sheet1'
            }

args19.kd_zeu = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_ZEU\Kd_ZEU_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

args19.kd_3m = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_3m\Kd_3m_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

args19.kd_6m = {
                'file_name': r'.\00_Dados\Dados_Victor\Kd_6m\Kd_6m_TRM.xlsx', 
                'sheet_name': 'Sheet1'
            }

data_obj19 = read_files(args19)

qaa_v6_fit19_665 = Qaa_v6_model_RUN_original(data_obj19, 'original_TRM_2019', 665)

###############################################################################
#########################  QAA V6 - Campanhas 2013/19  ########################
###############################################################################

##  Abrindo arquivos contendo dicionários de dados resultados da aplicação do QAA/Kd;
qaa_v6_fit19_665 = load_data(r'.\Results_Espectral\QAAv6_original_TRM_2019_665_Ajustado.pickle')

kd_ref = qaa_v6_fit19_665['Kd_ZEU'].copy()
kd_pred = qaa_v6_fit19_665['kd_qaa_665'].copy()

kd_ref19 = qaa_v6_fit19_665['Kd_ZEU'].copy()
kd_pred19 = qaa_v6_fit19_665['kd_qaa_665'].copy()

###############################################################################

dir_o19 = r'D:\prog\master\qaapy_kdpy\Results_Espectral\2019_lb0-665_Kd'

plot_kd_nxm(qaa_v6_fit19_665, 'kd_qaa_665', 'Kd_ZEU', 4, 5, dir_o19, '.png')

###############################################################################
########################    Estatísticas (N = 40)   ###########################
###############################################################################
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443_original = kd_ref.loc[443]         ##  Dado de referência (y_true)
kd_pred_443_original = kd_pred.loc[443]        ##  Dados simulados (y_pred)

##  Ajuste da curva de regressão linear para 443 nm.
fit_443_original = curve_fit(fun_reg, kd_ref_443_original, kd_pred_443_original, pred_ab, method = 'lm', maxfev = 100000)[0]

##  Calculo do coeficiente de determinação R²
res_443_original = kd_pred_443_original - fun_reg(kd_ref_443_original, fit_443_original[0], fit_443_original[1])
ss_res_443_original = np.sum(res_443_original ** 2)
ss_tot_443_original = np.sum((kd_pred_443_original - np.mean(kd_pred_443_original)) ** 2)
r2_443_original = 1 - (ss_res_443_original / ss_tot_443_original)

##  Calculo das estatísticas MAE e MAPE
mae_443_original = mae(kd_ref_443_original, kd_pred_443_original)
mape_443_original = mean_abs_perc_error(kd_ref_443_original, kd_pred_443_original)

## Calculo das estatisticas MSE e RMSE
mse_443_original = mse(kd_ref_443_original, kd_pred_443_original)
rmse_443_original = np.sqrt(mse_443_original)

## Regressão linear usando pacote do scipy para teste
sl_443_original, inter_443_original, r_443_original, p_443_original, stderr_443_original = lr(kd_ref_443_original, kd_pred_443_original)

##  Curva de melhor ajuste para plotagem
ref_443_original = fun_reg(x_fit, fit_443_original[0], fit_443_original[1])
ref_443_original = ref_443_original.rename(x_fit, axis = 'rows')

##  Texto para plotagem com dados dos resultados estatísticos e da regressão
eq_443_original = 'y = 0.5579x + 0.2686'
stat_443_original = 'R² = 0.2695; MAPE = 15.96%; RMSE = 0.1603 m^-1; n = 20'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443_original,
            kd_pred_443_original,
            ref_443_original,
            eq_443_original,
            stat_443_original,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (443 nm)',
            'Kd QAA (443 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 492 nm.
kd_ref_492_original = kd_ref.loc[492]
kd_pred_492_original = kd_pred.loc[492]

##  Ajuste da curva de regressão linear para 492 nm.
fit_492_original = curve_fit(fun_reg, kd_ref_492_original, kd_pred_492_original, pred_ab, method = 'lm', maxfev = 100000)[0]

res_492_original = kd_pred_492_original - fun_reg(kd_ref_492_original, fit_492_original[0], fit_492_original[1])
ss_res_492_original = np.sum(res_492_original ** 2)
ss_tot_492_original = np.sum((kd_pred_492_original - np.mean(kd_pred_492_original)) ** 2)
r2_492_original = 1 - (ss_res_492_original / ss_tot_492_original)

mae_492_original = mae(kd_ref_492_original, kd_pred_492_original)
mape_492_original = mean_abs_perc_error(kd_ref_492_original, kd_pred_492_original)

mse_492_original = mse(kd_ref_492_original, kd_pred_492_original)
rmse_492_original = np.sqrt(mse_492_original)

sl_492_original, inter_492_original, r_492_original, p_492_original, stderr_492_original = lr(kd_ref_492_original, kd_pred_492_original)

ref_492_original = fun_reg(x_fit, fit_492_original[0], fit_492_original[1])
ref_492_original = ref_492_original.rename(x_fit, axis = 'rows')

eq_492_original = 'y = 0.9587x - 0.1072'
stat_492_original = 'R² = 0.7852; MAPE = 22.64%; RMSE = 0.1380 m^-1; n = 20'

plot_ajuste(kd_ref_492_original,
            kd_pred_492_original,
            ref_492_original,
            eq_492_original,
            stat_492_original,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (492 nm)',
            'Kd QAA (492 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 560 nm
kd_ref_560_original = kd_ref.loc[560]
kd_pred_560_original = kd_pred.loc[560]

##  Ajuste da curva de regressão linear para 560 nm.
fit_560_original = curve_fit(fun_reg, kd_ref_560_original, kd_pred_560_original, pred_ab, method = 'lm', maxfev = 100000)[0]

res_560_original = kd_pred_560_original - fun_reg(kd_ref_560_original, fit_560_original[0], fit_560_original[1])
ss_res_560_original = np.sum(res_560_original ** 2)
ss_tot_560_original = np.sum((kd_pred_560_original - np.mean(kd_pred_560_original)) ** 2)
r2_560_original = 1 - (ss_res_560_original / ss_tot_560_original)

mae_560_original = mae(kd_ref_560_original, kd_pred_560_original)
mape_560_original = mean_abs_perc_error(kd_ref_560_original, kd_pred_560_original)

mse_560_original = mse(kd_ref_560_original, kd_pred_560_original)
rmse_560_original = np.sqrt(mse_560_original)

sl_560_original, inter_560_original, r_560_original, p_560_original, stderr_560_original = lr(kd_ref_560_original, kd_pred_560_original)

ref_560_original = fun_reg(x_fit, fit_560_original[0], fit_560_original[1])
ref_560_original = ref_560_original.rename(x_fit, axis = 'rows')

eq_560_original = 'y = 1.3756x - 0.2085'
stat_560_original = 'R² = 0.8155; MAPE = 19.02%; RMSE = 0.0766 m^-1; n = 20'

plot_ajuste(kd_ref_560_original,
            kd_pred_560_original,
            ref_560_original,
            eq_560_original,
            stat_560_original,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (560 nm)',
            'Kd QAA (560 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 665 nm
kd_ref_665_original = kd_ref.loc[665]
kd_pred_665_original = kd_pred.loc[665]

##  Ajuste da curva de regressão linear para 665 nm.
fit_665_original = curve_fit(fun_reg, kd_ref_665_original, kd_pred_665_original, pred_ab, method = 'lm', maxfev = 100000)[0]

res_665_original = kd_pred_665_original - fun_reg(kd_ref_665_original, fit_665_original[0], fit_665_original[1])
ss_res_665_original = np.sum(res_665_original ** 2)
ss_tot_665_original = np.sum((kd_pred_665_original - np.mean(kd_pred_665_original)) ** 2)
r2_665_original = 1 - (ss_res_665_original / ss_tot_665_original)

mae_665_original = mae(kd_ref_665_original, kd_pred_665_original)
mape_665_original = mean_abs_perc_error(kd_ref_665_original, kd_pred_665_original)

mse_665_original = mse(kd_ref_665_original, kd_pred_665_original)
rmse_665_original = np.sqrt(mse_665_original)

sl_665_original, inter_665_original, r_665_original, p_665_original, stderr_665_original = lr(kd_ref_665_original, kd_pred_665_original)

ref_665_original = fun_reg(x_fit, fit_665_original[0], fit_665_original[1])
ref_665_original = ref_665_original.rename(x_fit, axis = 'rows')

eq_665_original = 'y = 1.3516x - 0.1565'
stat_665_original = 'R² = 0.3463; MAPE = 11.52%; RMSE = 0.0932 m^-1; n = 20'

plot_ajuste(kd_ref_665_original,
            kd_pred_665_original,
            ref_665_original,
            eq_665_original,
            stat_665_original,
            [0, 5], [0, 5],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704_original = kd_ref.loc[704]
kd_pred_704_original = kd_pred.loc[704]

##  Ajuste da curva de regressão linear para 704 nm.
fit_704_original = curve_fit(fun_reg, kd_ref_704_original, kd_pred_704_original, pred_ab, method = 'lm', maxfev = 100000)[0]

res_704_original = kd_pred_704_original - fun_reg(kd_ref_704_original, fit_704_original[0], fit_704_original[1])
ss_res_704_original = np.sum(res_704_original ** 2)
ss_tot_704_original = np.sum((kd_pred_704_original - np.mean(kd_pred_704_original)) ** 2)
r2_704_original = 1 - (ss_res_704_original / ss_tot_704_original)

mae_704_original = mae(kd_ref_704_original, kd_pred_704_original)
mape_704_original = mean_abs_perc_error(kd_ref_704_original, kd_pred_704_original)

mse_704_original = mse(kd_ref_704_original, kd_pred_704_original)
rmse_704_original = np.sqrt(mse_704_original)

sl_704_original, inter_704_original, r_704_original, p_704_original, stderr_704_original = lr(kd_ref_704_original, kd_pred_704_original)

ref_704_original = fun_reg(x_fit, fit_704_original[0], fit_704_original[1])
ref_704_original = ref_704_original.rename(x_fit, axis = 'rows')

eq_704_original = 'y = 1.0751x + 0.2423'
stat_704_original = 'R² = 0.3129; MAPE = 36.84%; RMSE = 0.3169 m^-1; n = 20'

plot_ajuste(kd_ref_704_original,
            kd_pred_704_original,
            ref_704_original,
            eq_704_original,
            stat_704_original,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################

kd_ref19_443_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[443]
kd_ref19_492_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[492]
kd_ref19_560_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[560]
kd_ref19_665_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[665]

kd_pred19_443_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[443]
kd_pred19_492_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[492]
kd_pred19_560_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[560]
kd_pred19_665_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[665]

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(221)
plt.text(0.1, 2.2, 'a)', fontsize = 12)
plt.text(1, 0.6, eq_443_original, fontsize = 10)
plt.text(1, 0.4, 'R² = 0.2695; MAPE = 15.96%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0.1603 m^-1; n = 20', fontsize = 10)
plt.plot(ref_443_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref19_443_original, kd_pred19_443_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd - Field (443 nm)', fontsize = 14)
plt.ylabel('Kd - QAA (443 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('443 nm', fontsize = 14)
plt.show()

plt.subplot(222)
plt.text(1, 2.2, 'b)', fontsize = 12)
plt.text(1, 0.6, eq_492_original, fontsize = 10)
plt.text(1, 0.4, 'R² = 0.7852; MAPE = 22.64%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0.1380 m^-1; n = 20', fontsize = 10)
plt.plot(ref_492_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref19_492_original, kd_pred19_492_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd - Field (492 nm)', fontsize = 14)
plt.ylabel('Kd - QAA (492 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('492 nm', fontsize = 14)
plt.legend(['Fit', '1:1', '2019'], loc = 2, fontsize = 12) 
plt.show()

plt.subplot(223)
plt.text(0.1, 2.2, 'c)', fontsize = 12)
plt.text(1, 0.6, eq_560_original, fontsize = 10)
plt.text(1, 0.4, 'R² = 0.8155; MAPE = 19.02%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0.0766 m^-1; n = 20', fontsize = 10)
plt.plot(ref_560_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref19_560_original, kd_pred19_560_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd - Field (560 nm)', fontsize = 14)
plt.ylabel('Kd - QAA (560 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('560 nm', fontsize = 14)
plt.show()

plt.subplot(224)
plt.text(0.1, 2.2, 'd)', fontsize = 12)
plt.text(1, 0.6, eq_665_original, fontsize = 10)
plt.text(1, 0.4, 'R² = 0.3463; MAPE = 11.52%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0.0932 m^-1; n = 20', fontsize = 10)
plt.plot(ref_665_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref19_665_original, kd_pred19_665_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd - Field (665 nm)', fontsize = 14)
plt.ylabel('Kd - QAA (665 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('665 nm', fontsize = 14)
plt.show()

###############################################################################
###############    Estatísticas Algoritmo Global (N = 100)   ##################
###############################################################################
##  Regressão, plots e estatísticas para algoritmo global
kd_ref_global_original = pd.concat([kd_ref_443_original.copy(), kd_ref_492_original.copy(), kd_ref_560_original.copy(), kd_ref_665_original.copy(), kd_ref_704_original.copy()], axis = 0)
kd_pred_global_original = pd.concat([kd_pred_443_original.copy(), kd_pred_492_original.copy(), kd_pred_560_original.copy(), kd_pred_665_original.copy(), kd_pred_704_original.copy()], axis = 0)

##  Ajuste da curva de regressão linear para 443 nm.
fit_global_original = curve_fit(fun_reg, kd_ref_global_original, kd_pred_global_original, pred_ab, method = 'lm', maxfev = 100000)[0]

##  Calculo do coeficiente de determinação R²
res_global_original = kd_pred_global_original - fun_reg(kd_ref_global_original, fit_global_original[0], fit_global_original[1])
ss_res_global_original = np.sum(res_global_original ** 2)
ss_tot_global_original = np.sum((kd_pred_global_original - np.mean(kd_pred_global_original)) ** 2)
r2_global_original = 1 - (ss_res_global_original / ss_tot_global_original)

##  Calculo das estatísticas MAE e MAPE
mae_global_original = mae(kd_ref_global_original, kd_pred_global_original)
mape_global_original = mean_abs_perc_error(kd_ref_global_original, kd_pred_global_original)

## Calculo das estatisticas MSE e RMSE
mse_global_original = mse(kd_ref_global_original, kd_pred_global_original)
rmse_global_original = np.sqrt(mse_global_original)

## Regressão linear usando pacote do scipy para teste
sl_global_original, inter_global_original, r_global_original, p_global_original, stderr_global_original = lr(kd_ref_global_original, kd_pred_global_original)

##  Curva de melhor ajuste para plotagem
ref_global_original = fun_reg(x_fit, fit_global_original[0], fit_global_original[1])
ref_global_original = ref_global_original.rename(x_fit, axis = 'rows')

##  Texto para plotagem com dados dos resultados estatísticos e da regressão
eq_global_original = 'y = 1.2978x - 0.1833'
stat_global_original = 'R² = 0.6819; MAPE = 21.20%; RMSE = 0.1787 m^-1; n = 100'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_global_original,
            kd_pred_global_original,
            ref_global_original,
            eq_global_original,
            stat_global_original,
            [0, 2.5], [0, 2.5],
            'Kd - Field (GLOBAL MSI)',
            'Kd - QAA (GLOBAL MSI)',
            'Kd Field vs. Kd QAA - Global MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################

kd_ref19_443_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[443]
kd_ref19_492_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[492]
kd_ref19_560_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[560]
kd_ref19_665_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[665]
kd_ref19_704_original = qaa_v6_fit19_665['Kd_ZEU'].copy().loc[704]

kd_pred19_443_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[443]
kd_pred19_492_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[492]
kd_pred19_560_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[560]
kd_pred19_665_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[665]
kd_pred19_704_original = qaa_v6_fit19_665['kd_qaa_665'].copy().loc[704]

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(111)
plt.text(1, 0.6, eq_global_original, fontsize = 10)
plt.text(1, 0.4, 'R² = 0.6819; MAPE = 21.20%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0.1787 m^-1; n = 100', fontsize = 10)
plt.plot(ref_global_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref19_443_original, kd_pred19_443_original, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref19_492_original, kd_pred19_492_original, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref19_560_original, kd_pred19_560_original, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref19_665_original, kd_pred19_665_original, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.scatter(kd_ref19_704_original, kd_pred19_704_original, marker = '^', facecolors = 'none', edgecolors = '#FF1493')
plt.xlabel('Kd - Field (GLOBAL)', fontsize = 14)
plt.ylabel('Kd - QAA (GLOBAL)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.legend(['Fit', '1:1', '443 nm (2019)', '492 nm', '560 nm', '665 nm', '704 nm'], loc = 4, fontsize = 12) 
plt.show()
