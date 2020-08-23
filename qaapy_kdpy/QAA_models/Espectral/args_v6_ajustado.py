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
from qaapy_kdpy.utils.plot_Kd_campo_vs_qaa import plot_ajuste, plot_kd_nxm
from qaapy_kdpy.utils.stats import mean_abs_perc_error

###############################################################################
#########################  QAA V6 - Campanhas 2013/19  ########################
###############################################################################

##  Abrindo arquivos contendo dicionários de dados resultados da aplicação do QAA/Kd;
qaa_v6_fit13_560 = load_data(r'.\Results_Espectral\QAAv6_TRM_2013_560_Ajustado.pickle')
qaa_v6_fit19_560 = load_data(r'.\Results_Espectral\QAAv6_TRM_2019_560_Ajustado.pickle')

kd_ref = pd.concat([qaa_v6_fit13_560['Kd_ZEU'].copy(), qaa_v6_fit19_560['Kd_ZEU'].copy()], axis = 1)
kd_pred = pd.concat([qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy(), qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy()], axis = 1)

kd_ref13 = qaa_v6_fit13_560['Kd_ZEU'].copy().drop(['P25', 'P26'], axis = 1)
kd_pred13 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().drop(['P25', 'P26'], axis = 1)

kd_ref19 = qaa_v6_fit19_560['Kd_ZEU'].copy()
kd_pred19 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy()

###############################################################################

##  Plotando resultados exemplo de cada uma das campanhas;
plt.figure(figsize = (10, 5))

plt.subplot(121)
plt.plot(qaa_v6_fit13_560['Kd_ZEU']['P15'], 'r-')
plt.plot(qaa_v6_fit13_560['kd_qaa_560_acs_492']['P15'], 'b-')
plt.text(405, 1.75, 'Chl-a = 7,11 ug/L ; TSS = 2,5 mg/L', fontsize = 14)
plt.text(405, 1.6, 'aCDOM(443) = 0,38 1/m', fontsize = 14)
plt.xlabel('Comprimento de Onda [nm]', fontsize = 14)
plt.ylabel('Kd [m ^ -1]', fontsize = 14)
plt.xlim([400, 700])
plt.ylim([0, 2])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], ['0,00', '0,25', '0,50', '0,75', '1,00', '1,25', '1,50', '1,75', '2,00'], fontsize = 14)
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
plt.title('2013 - P15', fontsize = 14)
#plt.legend(['Kd in-situ', 'Kd QAA'], loc = 1, fontsize = 14) 
plt.show()

plt.subplot(122)
plt.plot(qaa_v6_fit19_560['Kd_ZEU']['TRM09'], 'r-')
plt.plot(qaa_v6_fit19_560['kd_qaa_560_acs_492']['TRM09'], 'b-')
plt.text(480, 1.5, 'Chl-a = 4,24 ug/L ; TSS = 1,1 mg/L', fontsize = 14)
plt.text(480, 1.35, 'aCDOM(443) = 0,22 1/m', fontsize = 14)
plt.xlabel('Comprimento de Onda [nm]', fontsize = 14)
plt.xlim([400, 700])
plt.ylim([0, 2])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2], ['0,00', '0,25', '0,50', '0,75', '1,00', '1,25', '1,50', '1,75', '2,00'], fontsize = 14)
plt.rc('xtick', labelsize = 14) 
plt.rc('ytick', labelsize = 14) 
plt.title('2019 - TRM09', fontsize = 14)
plt.legend(['Kd in-situ', 'Kd QAA'], loc = 1, fontsize = 14) 
plt.show()

###############################################################################

dir_o13 = r'D:\prog\Master_Victor_Curtarelli\qaapy_kdpy\Results_Espectral\2013_lb0-560_Kd'
dir_o19 = r'D:\prog\Master_Victor_Curtarelli\qaapy_kdpy\Results_Espectral\2019_lb0-560_Kd'

plot_kd_nxm(qaa_v6_fit13_560, 'kd_qaa_560_acs_492', 'Kd_ZEU', 4, 5, dir_o13, '.png')
plot_kd_nxm(qaa_v6_fit19_560, 'kd_qaa_560_acs_492', 'Kd_ZEU', 4, 5, dir_o19, '.png')

###############################################################################
########################    Estatísticas (N = 42)   ###########################
###############################################################################
'''
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = kd_ref.loc[443].drop(['P25', 'P26'])         ##  Dado de referência (y_true)
kd_pred_443 = kd_pred.loc[443].drop(['P25', 'P26'])        ##  Dados simulados (y_pred)

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
eq_443 = 'y = 0,8904x + 0,4460'
stat_443 = 'R² = 0,6012, MAPE = 40,14%; RMSE = 0,3892; n = 40'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
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
kd_ref_492 = kd_ref.loc[492]
kd_pred_492 = kd_pred.loc[492]

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

eq_492 = 'y = 1,6880x - 0,3107'
stat_492 = 'R² = 0,9681, MAPE = 22,50%; RMSE = 0,4300; n = 42'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
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
kd_ref_560 = kd_ref.loc[560]
kd_pred_560 = kd_pred.loc[560]
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

eq_560 = 'y = 2,0469x - 0,3455'
stat_560 = 'R² = 0,9662, MAPE = 21,11%; RMSE = 0,3459; n = 42'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
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
kd_ref_665 = kd_ref.loc[665]
kd_pred_665 = kd_pred.loc[665]

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

eq_665 = 'y = 2,2164x - 0,7256'
stat_665 = 'R² = 0,8991; MAPE = 21,44%; RMSE = 0,3520; n = 42'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = kd_ref.loc[704]
kd_pred_704 = kd_pred.loc[704]

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

eq_704 = 'y = 1,9098x - 0,5255'
stat_704 = 'R² = 0,8381; MAPE = 31,62%; RMSE = 0,4147; n = 42'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)
'''

###############################################################################
########################    Estatísticas (N = 40)   ###########################
###############################################################################
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = kd_ref.loc[443].drop(['P25', 'P26'])         ##  Dado de referência (y_true)
kd_pred_443 = kd_pred.loc[443].drop(['P25', 'P26'])        ##  Dados simulados (y_pred)

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
eq_443 = 'y = 0,8904x + 0,4460'
stat_443 = 'R² = 0,6012, MAPE = 40,14%; RMSE = 0,3892 m^-1; n = 40'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
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
kd_ref_492 = kd_ref.loc[492].drop(['P25', 'P26'])
kd_pred_492 = kd_pred.loc[492].drop(['P25', 'P26'])

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

eq_492 = 'y = 1,0918x + 0,0655'
stat_492 = 'R² = 0,7917, MAPE = 20,63%; RMSE = 0,1613 m^-1; n = 40'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
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
kd_ref_560 = kd_ref.loc[560].drop(['P25', 'P26'])
kd_pred_560 = kd_pred.loc[560].drop(['P25', 'P26'])

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

eq_560 = 'y = 1,3702x - 0,0885'
stat_560 = 'R² = 0,6406, MAPE = 17,46%; RMSE = 0,1068 m^-1; n = 40'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
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
kd_ref_665 = kd_ref.loc[665].drop(['P25', 'P26'])
kd_pred_665 = kd_pred.loc[665].drop(['P25', 'P26'])

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

eq_665 = 'y = 1,2424x - 0,0491'
stat_665 = 'R² = 0,6310; MAPE = 17,80%; RMSE = 0,1597 m^-1; n = 40'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
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
kd_ref_704 = kd_ref.loc[704].drop(['P25', 'P26'])
kd_pred_704 = kd_pred.loc[704].drop(['P25', 'P26'])

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

eq_704 = 'y = 0,8022x + 0,4114'
stat_704 = 'R² = 0,4870; MAPE = 28,91%; RMSE = 0,2628 m^-1; n = 40'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################

kd_ref13_443 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[443].drop(['P25', 'P26'])
kd_ref13_492 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[492].drop(['P25', 'P26'])
kd_ref13_560 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[560].drop(['P25', 'P26'])
kd_ref13_665 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[665].drop(['P25', 'P26'])

kd_ref19_443 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[443]
kd_ref19_492 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[492]
kd_ref19_560 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[560]
kd_ref19_665 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[665]

kd_pred13_443 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[443].drop(['P25', 'P26'])
kd_pred13_492 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[492].drop(['P25', 'P26'])
kd_pred13_560 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[560].drop(['P25', 'P26'])
kd_pred13_665 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[665].drop(['P25', 'P26'])

kd_pred19_443 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[443]
kd_pred19_492 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[492]
kd_pred19_560 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[560]
kd_pred19_665 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[665]

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(221)
plt.text(0.1, 2.2, 'a)', fontsize = 12)
plt.text(1, 0.6, eq_443, fontsize = 10)
plt.text(1, 0.4, 'R² = 0,6012; MAPE = 40,14%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0,3892 m^-1; n = 40', fontsize = 10)
plt.plot(ref_443, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_443, kd_pred13_443, c = 'b', marker = 'x')
plt.scatter(kd_ref19_443, kd_pred19_443, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (443 nm)', fontsize = 14)
plt.ylabel('Kd QAA (443 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('443 nm', fontsize = 14)
plt.show()

plt.subplot(222)
plt.text(1, 2.2, 'b)', fontsize = 12)
plt.text(1, 0.6, eq_492, fontsize = 10)
plt.text(1, 0.4, 'R² = 0,7917; MAPE = 20,63%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0,1613 m^-1; n = 40', fontsize = 10)
plt.plot(ref_492, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_492, kd_pred13_492, c = 'b', marker = 'x')
plt.scatter(kd_ref19_492, kd_pred19_492, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (492 nm)', fontsize = 14)
plt.ylabel('Kd QAA (492 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('492 nm', fontsize = 14)
plt.legend(['Ajuste', '1:1', '2013', '2019'], loc = 2, fontsize = 12) 
plt.show()

plt.subplot(223)
plt.text(0.1, 2.2, 'c)', fontsize = 12)
plt.text(1, 0.6, eq_560, fontsize = 10)
plt.text(1, 0.4, 'R² = 0,6406; MAPE = 17,46%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0,1068 m^-1; n = 40', fontsize = 10)
plt.plot(ref_560, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_560, kd_pred13_560, c = 'b', marker = 'x')
plt.scatter(kd_ref19_560, kd_pred19_560, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (560 nm)', fontsize = 14)
plt.ylabel('Kd QAA (560 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('560 nm', fontsize = 14)
plt.show()

plt.subplot(224)
plt.text(0.1, 2.2, 'd)', fontsize = 12)
plt.text(1, 0.6, eq_665, fontsize = 10)
plt.text(1, 0.4, 'R² = 0,6310; MAPE = 17,80%', fontsize = 10)
plt.text(1, 0.2, 'RMSE = 0,1597 m^-1; n = 40', fontsize = 10)
plt.plot(ref_665, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_665, kd_pred13_665, c = 'b', marker = 'x')
plt.scatter(kd_ref19_665, kd_pred19_665, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd ZEU (665 nm)', fontsize = 14)
plt.ylabel('Kd QAA (665 nm)', fontsize = 14)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5], ['0,0', '0,5', '1,0', '1,5', '2,0', '2,5'])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('665 nm', fontsize = 14)
plt.show()

###############################################################################

x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(241)
plt.text(0.1, 2.3, 'a)', fontsize = 12)
plt.text(0.75, 0.45, eq_443, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.6012; MAPE = 40.14%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.3892 m^-1; n = 40', fontsize = 9)
plt.plot(ref_443, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_443, kd_pred13_443, c = 'b', marker = 'x')
plt.scatter(kd_ref19_443, kd_pred19_443, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(242)
plt.text(0.1, 2.3, 'b)', fontsize = 12)
plt.text(0.75, 0.45, eq_492, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.7917; MAPE = 20.63%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.1613 m^-1; n = 40', fontsize = 9)
plt.plot(ref_492, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_492, kd_pred13_492, c = 'b', marker = 'x')
plt.scatter(kd_ref19_492, kd_pred19_492, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(243)
plt.text(0.1, 2.3, 'c)', fontsize = 12)
plt.text(0.75, 0.45, eq_560, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.6406; MAPE = 17.46%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.1068 m^-1; n = 40', fontsize = 9)
plt.plot(ref_560, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_560, kd_pred13_560, c = 'b', marker = 'x')
plt.scatter(kd_ref19_560, kd_pred19_560, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(244)
plt.text(0.1, 2.3, 'd)', fontsize = 12)
plt.text(0.75, 0.45, eq_665, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.6310; MAPE = 17.80%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.1597 m^-1; n = 40', fontsize = 9)
plt.plot(ref_665, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_665, kd_pred13_665, c = 'b', marker = 'x')
plt.scatter(kd_ref19_665, kd_pred19_665, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.legend(['Fit', '1:1', '2013', '2019'], loc = 1, fontsize = 8) 
plt.show()

plt.subplot(245)
plt.text(0.1, 2.3, 'e)', fontsize = 12)
plt.text(0.75, 0.45, eq_443_original, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.7104; MAPE = 17.65%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.2042 m^-1; n = 40', fontsize = 9)
plt.plot(ref_443_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_443_original, kd_pred13_443_original, c = 'b', marker = 'x')
plt.scatter(kd_ref19_443_original, kd_pred19_443_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(246)
plt.text(0.1, 2.3, 'f)', fontsize = 12)
plt.text(0.75, 0.45, eq_492_original, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.8072; MAPE = 20.25%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.1440 m^-1; n = 40', fontsize = 9)
plt.plot(ref_492_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_492_original, kd_pred13_492_original, c = 'b', marker = 'x')
plt.scatter(kd_ref19_492_original, kd_pred19_492_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(247)
plt.text(0.1, 2.3, 'g)', fontsize = 12)
plt.text(0.75, 0.45, eq_560_original, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.6959; MAPE = 18.01%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.0769 m^-1; n = 40', fontsize = 9)
plt.plot(ref_560_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_560_original, kd_pred13_560_original, c = 'b', marker = 'x')
plt.scatter(kd_ref19_560_original, kd_pred19_560_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(248)
plt.text(0.1, 2.3, 'h)', fontsize = 12)
plt.text(0.75, 0.45, eq_665_original, fontsize = 9)
plt.text(0.75, 0.25, 'R² = 0.6454; MAPE = 12.48%', fontsize = 9)
plt.text(0.75, 0.05, 'RMSE = 0.1092 m^-1; n = 40', fontsize = 9)
plt.plot(ref_665_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_665_original, kd_pred13_665_original, c = 'b', marker = 'x')
plt.scatter(kd_ref19_665_original, kd_pred19_665_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

###############################################################################
###############    Estatísticas Algoritmo Global (N = 200)   ##################
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
eq_global = 'y = 1.2842x - 0.0231'
stat_global = 'R² = 0.8325; MAPE = 24.99 %; RMSE = 0.2381 m^-1; n = 200'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_global,
            kd_pred_global,
            ref_global,
            eq_global,
            stat_global,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (GLOBAL MSI)',
            'Kd QAA (GLOBAL MSI)',
            'Kd ZEU vs. Kd QAA - Global MSI',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
kd_ref13_443 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[443].drop(['P25', 'P26'])
kd_ref13_492 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[492].drop(['P25', 'P26'])
kd_ref13_560 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[560].drop(['P25', 'P26'])
kd_ref13_665 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[665].drop(['P25', 'P26'])
kd_ref13_704 = qaa_v6_fit13_560['Kd_ZEU'].copy().loc[704].drop(['P25', 'P26'])

kd_ref19_443 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[443]
kd_ref19_492 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[492]
kd_ref19_560 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[560]
kd_ref19_665 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[665]
kd_ref19_704 = qaa_v6_fit19_560['Kd_ZEU'].copy().loc[704]

kd_pred13_443 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[443].drop(['P25', 'P26'])
kd_pred13_492 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[492].drop(['P25', 'P26'])
kd_pred13_560 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[560].drop(['P25', 'P26'])
kd_pred13_665 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[665].drop(['P25', 'P26'])
kd_pred13_704 = qaa_v6_fit13_560['kd_qaa_560_acs_492'].copy().loc[704].drop(['P25', 'P26'])

kd_pred19_443 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[443]
kd_pred19_492 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[492]
kd_pred19_560 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[560]
kd_pred19_665 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[665]
kd_pred19_704 = qaa_v6_fit19_560['kd_qaa_560_acs_492'].copy().loc[704]

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(121)
plt.text(0.1, 2.3, 'a)', fontsize = 12)
plt.text(1.4, 0.6, eq_global, fontsize = 10)
plt.text(1.4, 0.4, 'R² = 0.8325; MAPE = 24.99%', fontsize = 10)
plt.text(1.4, 0.2, 'RMSE = 0.2381 m^-1; n = 200', fontsize = 10)
plt.plot(ref_global, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_443, kd_pred13_443, c = '#0000FF', marker = 'x')
plt.scatter(kd_ref13_492, kd_pred13_492, c = '#4682B4', marker = 'x')
plt.scatter(kd_ref13_560, kd_pred13_560, c = '#008000', marker = 'x')
plt.scatter(kd_ref13_665, kd_pred13_665, c = '#800000', marker = 'x')
plt.scatter(kd_ref13_704, kd_pred13_704, c = '#DC143C', marker = 'x')
plt.scatter(kd_ref19_443, kd_pred19_443, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref19_492, kd_pred19_492, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref19_560, kd_pred19_560, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref19_665, kd_pred19_665, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.scatter(kd_ref19_704, kd_pred19_704, marker = '^', facecolors = 'none', edgecolors = '#FF1493')
plt.xlabel('Kd-measured', fontsize = 12)
plt.ylabel('Kd-SA', fontsize = 12)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('Estatísticas Globais MSI', fontsize = 14)
plt.show()

plt.subplot(122)
plt.text(0.1, 2.3, 'b)', fontsize = 12)
plt.text(0.2, 2.3, eq_global_original, fontsize = 10)
plt.text(0.2, 2.1, 'R² = 0.6490; MAPE = 20.02%', fontsize = 10)
plt.text(0.2, 1.9, 'RMSE = 0.1802 m^-1; n = 200', fontsize = 10)
plt.plot(ref_global_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref13_443_original, kd_pred13_443_original, c = '#0000FF', marker = 'x')
plt.scatter(kd_ref13_492_original, kd_pred13_492_original, c = '#4682B4', marker = 'x')
plt.scatter(kd_ref13_560_original, kd_pred13_560_original, c = '#008000', marker = 'x')
plt.scatter(kd_ref13_665_original, kd_pred13_665_original, c = '#800000', marker = 'x')
plt.scatter(kd_ref13_704_original, kd_pred13_704_original, c = '#DC143C', marker = 'x')
plt.scatter(kd_ref19_443_original, kd_pred19_443_original, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref19_492_original, kd_pred19_492_original, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref19_560_original, kd_pred19_560_original, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref19_665_original, kd_pred19_665_original, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.scatter(kd_ref19_704_original, kd_pred19_704_original, marker = '^', facecolors = 'none', edgecolors = '#FF1493')
plt.xlabel('Kd-measured', fontsize = 12)
plt.xlim([0, 2.5])
plt.ylim([0, 2.5])
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5])
plt.yticks([])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.legend(['Fit', '1:1', '443 nm (2013)', '492 nm', '560 nm', '665 nm', '704 nm', '443 nm (2019)', '492 nm', '560 nm', '665 nm', '704 nm'], loc = 4, fontsize = 12) 
plt.show()

###############################################################################
########################    Estatísticas (N = 36)   ###########################
###############################################################################
'''
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = kd_ref.loc[443].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])         ##  Dado de referência (y_true)
kd_pred_443 = kd_pred.loc[443].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])        ##  Dados simulados (y_pred)

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
eq_443 = 'y = 0,9376x + 0,3772'
stat_443 = 'R² = 0,7192, MAPE = 36,90%; RMSE = 0,3534; n = 36'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
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
kd_ref_492 = kd_ref.loc[492].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_492 = kd_pred.loc[492].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

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

eq_492 = 'y = 1,0842x + 0,0749'
stat_492 = 'R² = 0,7959, MAPE = 21,34%; RMSE = 0,1667; n = 36'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
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
kd_ref_560 = kd_ref.loc[560].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_560 = kd_pred.loc[560].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

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

eq_560 = 'y = 1,3527x - 0,0736'
stat_560 = 'R² = 0,6566, MAPE = 18,31%; RMSE = 0,1114; n = 36'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
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
kd_ref_665 = kd_ref.loc[665].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_665 = kd_pred.loc[665].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

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

eq_665 = 'y = 1,2225x - 0,0256'
stat_665 = 'R² = 0,6479; MAPE = 18,81%; RMSE = 0,1671; n = 36'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (665 nm)',
            'Kd QAA (665 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = kd_ref.loc[704].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])
kd_pred_704 = kd_pred.loc[704].drop(['P25', 'P26', 'TRMK1', 'TRMK2', 'TRMK3', 'TRMK4'])

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

eq_704 = 'y = 0,8121x + 0,4101'
stat_704 = 'R² = 0,5116; MAPE = 29,77%; RMSE = 0,2695; n = 40'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0, 2.5], [0, 2.5],
            'Kd ZEU (704 nm)',
            'Kd QAA (704 nm)',
            'Kd ZEU vs. Kd QAA - lb0 = 560 nm',
            'Ajuste',
            'Estações',
            leg_loc = 0,
            x1 = True)

'''
