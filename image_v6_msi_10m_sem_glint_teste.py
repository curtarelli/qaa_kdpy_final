'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------


'''
##  Importando pacotes básicos necessários para o código
import gc
from osgeo import gdal
from gdalconst import *
import pandas as pd
import numpy as np
import math as mt
import datetime as dt
import scipy.ndimage
from matplotlib import pyplot as plt

from scipy.stats import linregress as lr
from scipy.optimize import curve_fit

##  Importando pacotes básicos para calculos estatísticos.
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_squared_error as mse

##  Modelo QAA/Kd para se aplicar às bandas do Sentinel 2A msi
from qaapy_kdpy.QAA_models.Espectral.Qaa_v6_model import Qaa_v6_model_msi

##  Funções auxiliares para carregamento dos dados
from qaapy_kdpy.definition import create_args
from qaapy_kdpy.utils.image_app import *
from qaapy_kdpy.file import load_data

from qaapy_kdpy.empirical_calc import *
from qaapy_kdpy.utils.plot_Kd_campo_vs_qaa import plot_ajuste
from qaapy_kdpy.utils.stats import mean_abs_perc_error

##  Diretório onde se encontram as imagens
dir_tif = r'.\00_Dados\Imagens'
dir_b11 = r'.\00_Dados\Imagens\utm_mask_T23KMV_20190701T131251_B11.tif'

##  Carregando as bandas de interesse (B1 a B5)
b1, b2, b3, b4, b5 = open_msi_b1_to_b5(dir_tif)

##  Carregando banda B11 para remoção de glint;
b11 = open_tif(dir_b11)

b1 = np.delete(b1, 0, axis = 0)
b2 = np.delete(b2, [0, 1, 2, 3], axis = 0)
b3 = np.delete(b3, [0, 1, 2, 3], axis = 0)
b4 = np.delete(b4, [0, 1, 2, 3], axis = 0)
b5 = np.delete(b5, [0, 1], axis = 0)
b11 = np.delete(b11, [0, 1], axis = 0)

##  Reamostragem espacial das bandas B1 (60m / 6) e B5 (20 m / 2) [RE = 10 metros];
##  mode = ['reflect', 'constant', 'nearest', 'mirror', 'wrap'];
b1 = scipy.ndimage.zoom(b1, 6, order = 0, mode = 'nearest')
b5 = scipy.ndimage.zoom(b5, 2, order = 0, mode = 'nearest')
b11 = scipy.ndimage.zoom(b11, 2, order = 0, mode = 'nearest')

b1 = b1 - b11
b2 = b2 - b11
b3 = b3 - b11
b4 = b4 - b11
b5 = b5 - b11

b1 = b1 / mt.pi
b2 = b2 / mt.pi
b3 = b3 / mt.pi
b4 = b4 / mt.pi
b5 = b5 / mt.pi

##  Carregando dados do coeficiente de absorção da água pura
aw = pd.read_excel(r'.\00_Dados\Dados_Victor\aw_pope97_s2a_PAR.xlsx',
                   sheet_name = 'Sheet1',
                   header = 0,
                   index_col = 0)
    
##  Carregando dados do coeficiente de retroespalhamento da água pura
bbw = pd.read_excel(r'.\00_Dados\Dados_Victor\bbw_zhang_s2a_PAR.xlsx',
                    sheet_name = 'Sheet1',
                    header = 0,
                    index_col = 0)

##  Latitude central do reservatório (imagem)
lat = mt.radians(-18.6087777777778)

##  Data da imagem com a hora de aquisição
date = dt.datetime(2019, 7, 1, 13, 12, 51)
##  Hora decimal da aquisição
dd = date.hour + (date.minute / 60) + (date.second / (60 * 60))
##  Dia juliano (dia absoluto no ano)
dy = date.timetuple().tm_yday
##  Calculo da elevação solar
h = mt.radians(abs((dd - 12) * 15))
##  Delta para calculo do ângulo solar zenital
delta = mt.radians(23.45 * mt.sin((360 / 365) * (dy - 80)))
##  Calculo do ângulo solar zenital
theta_s = mt.degrees(mt.acos((mt.sin(lat) * mt.sin(delta)) + (mt.cos(lat) * mt.cos(delta) * mt.cos(h))))

##  
ref_band = r'.\00_Dados\Imagens\utm_mask_T23KMV_20190701T131251_B02.tif'
dataset_de_referencia = gdal.Open(ref_band, GA_ReadOnly)

##  Aplicando o modelo QAA v6 nas imagens carregadas
kd_msi_b1 = Qaa_v6_model_msi(b1, 443, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b1, r'.\03_Results_image\QAAv6_msi_2019_10m_sem_glint_443nm_teste.tif', dataset_de_referencia)
del kd_msi_b1
gc.collect()

kd_msi_b2 = Qaa_v6_model_msi(b2, 492, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b2, r'.\03_Results_image\QAAv6_msi_2019_10m_sem_glint_492nm_teste.tif', dataset_de_referencia)
del kd_msi_b2
gc.collect()

kd_msi_b3 = Qaa_v6_model_msi(b3, 560, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b3, r'.\03_Results_image\QAAv6_msi_2019_10m_sem_glint_560nm_teste.tif', dataset_de_referencia)
del kd_msi_b3
gc.collect()

kd_msi_b4 = Qaa_v6_model_msi(b4, 665, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b4, r'.\03_Results_image\QAAv6_msi_2019_10m_sem_glint_665nm_teste.tif', dataset_de_referencia)
del kd_msi_b4
gc.collect()

kd_msi_b5 = Qaa_v6_model_msi(b5, 704, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b5, r'.\03_Results_image\QAAv6_msi_2019_10m_sem_glint_704nm_teste.tif', dataset_de_referencia)
del kd_msi_b5
gc.collect()

###############################################################################

plt.figure(figsize = (10, 10))
plt.plot()
plt.title('Kd B2')
plt.imshow(kd_msi_b2, cmap='gray')
plt.show()

plt.figure(figsize = (10, 10))
plt.plot()
plt.title('Kd B3')
plt.imshow(kd_msi_b3, cmap='gray')
plt.show()

###############################################################################
########################    Estatísticas (N = 20)   ###########################
###############################################################################
qaa_v6_fit19_560 = load_data(r'.\Results_Espectral\QAAv6_TRM_2019_560_Ajustado.pickle')

kd_ref = qaa_v6_fit19_560['Kd_ZEU'].copy()
kd_pred = pd.read_excel(r'.\03_Results_Image\TRM_2019_QAAv6_msi_teste.xlsx',
                        sheet_name = 'Sheet1',
                        header = 0,
                        index_col = 0)

###############################################################################
##  Criação de linha 1:1 para aplicação nos plots
pred_ab = [1, 0]
x_fit = pd.Series(np.arange(0, 11, 0.01))

##  Função de regressão linear para aplicar nos dados
fun_reg = lambda x, a, b: (a * x) + b

##  Regressão, plots e estatísticas para 443 nm
kd_ref_443 = kd_ref.loc[443]                     ##  Dado de referência (y_true)
kd_pred_443 = kd_pred['Kd_QAAv6_10m_s-glint_443']        ##  Dados simulados (y_pred)

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
eq_443 = '443 nm = 0.97x - 0.04'
stat_443 = 'R² = 0.44; MAPE = 13.49%; RMSE = 0.14 [${m}^{-1}$]; n = 20'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443,
            kd_pred_443,
            ref_443,
            eq_443,
            stat_443,
            [0, 2.5], [0, 2.5],
            'Kd-measured (443 nm)',
            'Kd-SA MSI (443 nm)',
            'Kd-measured vs. Kd-SA MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 492 nm.
kd_ref_492 = kd_ref.loc[492]
kd_pred_492 = kd_pred['Kd_QAAv6_10m_s-glint_492']

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

eq_492 = '492 nm = 1.08x + 0.03'
stat_492 = 'R² = 0.70; MAPE = 14.09%; RMSE = 0.09 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_492,
            kd_pred_492,
            ref_492,
            eq_492,
            stat_492,
            [0, 2.5], [0, 2.5],
            'Kd-measured (492 nm)',
            'Kd-SA MSI (492 nm)',
            'Kd-measured vs. Kd-SA MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 560 nm
kd_ref_560 = kd_ref.loc[560]
kd_pred_560 = kd_pred['Kd_QAAv6_10m_s-glint_560']

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

eq_560 = '560 nm = 0.93x + 0.17'
stat_560 = 'R² = 0.56; MAPE = 38.34%; RMSE = 0.15 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_560,
            kd_pred_560,
            ref_560,
            eq_560,
            stat_560,
            [0, 2.5], [0, 2.5],
            'Kd-measured (560 nm)',
            'Kd-SA MSI (560 nm)',
            'Kd-measured vs. Kd-SA MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 665 nm
kd_ref_665 = kd_ref.loc[665]
kd_pred_665 = kd_pred['Kd_QAAv6_10m_s-glint_665']

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

eq_665 = '665 nm = -0.35x + 1.15'
stat_665 = 'R² = 0.01; MAPE = 38.82%; RMSE = 0.28 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_665,
            kd_pred_665,
            ref_665,
            eq_665,
            stat_665,
            [0, 5], [0, 5],
            'Kd-measured (665 nm)',
            'Kd-SA MSI (665 nm)',
            'Kd-measured vs. Kd-SA MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################
##  Regressão, plots e estatísticas para 704 nm
kd_ref_704 = kd_ref.loc[704]
kd_pred_704 = kd_pred['Kd_QAAv6_10m_s-glint_704']

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

eq_704 = '704 nm = 0.27x + 0.89'
stat_704 = 'R² = 0.14; MAPE = 34.86%; RMSE = 0.29 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_704,
            kd_pred_704,
            ref_704,
            eq_704,
            stat_704,
            [0, 2.5], [0, 2.5],
            'Kd-measured (704 nm)',
            'Kd-SA MSI (704 nm)',
            'Kd-measured vs. Kd-SA MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################

x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(221)
plt.text(0.05, 1.4, 'a)', fontsize = 16)
plt.text(0.4, 0.25, eq_443, fontsize = 16)
plt.text(0.4, 0.15, 'R² = 0.44; MAPE = 13 %', fontsize = 16)
plt.text(0.4, 0.05, 'RMSE = 0.14 [1/m]; n = 20', fontsize = 16)
plt.plot(ref_443, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured [1/m]', fontsize = 16)
plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(222)
plt.text(0.05, 1.4, 'b)', fontsize = 16)
plt.text(0.4, 0.25, eq_492, fontsize = 16)
plt.text(0.4, 0.15, 'R² = 0.70; MAPE = 14 %', fontsize = 16)
plt.text(0.4, 0.05, 'RMSE = 0.10 [1/m]; n = 20', fontsize = 16)
plt.plot(ref_492, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured [1/m]', fontsize = 16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([])
plt.legend(['Fit', '1:1', 'MSI/19'], loc = 1, fontsize = 16) 
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(223)
plt.text(0.05, 1.4, 'c)', fontsize = 16)
plt.text(0.4, 0.25, eq_560, fontsize = 16)
plt.text(0.4, 0.15, 'R² = 0.56; MAPE = 38 %', fontsize = 16)
plt.text(0.4, 0.05, 'RMSE = 0.15 [1/m]; n = 20', fontsize = 16)
plt.plot(ref_560, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured [1/m]', fontsize =  16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.ylabel('Kd-SA MSI', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(224)
plt.text(0.05, 1.4, 'd)', fontsize = 16)
plt.text(0.4, 0.25, eq_665, fontsize = 16)
plt.text(0.4, 0.15, 'R² = 0.01; MAPE = 39 %', fontsize = 16)
plt.text(0.4, 0.05, 'RMSE = 0.28 [1/m]; n = 20', fontsize = 16)
plt.plot(ref_665, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured [1/m]', fontsize = 16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

###############################################################################

plt.figure(figsize = (10, 10))

plt.subplot(241)
plt.text(0.05, 1.4, 'a)', fontsize = 16)
plt.text(0.4, 0.25, eq_443, fontsize = 14)
plt.text(0.4, 0.15, 'R² = 0.44; MAPE = 13 %', fontsize = 14)
plt.text(0.4, 0.05, 'RMSE = 0.14 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_443, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured [1/m]', fontsize = 16)
plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(242)
plt.text(0.05, 1.4, 'b)', fontsize = 16)
plt.text(0.4, 0.25, eq_492, fontsize = 14)
plt.text(0.4, 0.15, 'R² = 0.70; MAPE = 14 %', fontsize = 14)
plt.text(0.4, 0.05, 'RMSE = 0.10 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_492, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured [1/m]', fontsize = 16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(243)
plt.text(0.05, 1.4, 'c)', fontsize = 16)
plt.text(0.4, 0.25, eq_560, fontsize = 14)
plt.text(0.4, 0.15, 'R² = 0.56; MAPE = 38 %', fontsize = 14)
plt.text(0.4, 0.05, 'RMSE = 0.15 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_560, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured [1/m]', fontsize =  16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)

plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(244)
plt.text(0.05, 1.4, 'd)', fontsize = 16)
plt.text(0.4, 0.25, eq_665, fontsize = 14)
plt.text(0.4, 0.15, 'R² = 0.01; MAPE = 39 %', fontsize = 14)
plt.text(0.4, 0.05, 'RMSE = 0.28 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_665, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured [1/m]', fontsize = 16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([])
plt.legend(['Fit', '1:1', 'MSI/19'], loc = 1, fontsize = 16) 
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(245)
plt.text(0.05, 1.4, 'e)', fontsize = 16)
plt.text(0.05, 1.25, eq_443_original, fontsize = 14)
plt.text(0.05, 1.15, 'R² = 0.57; MAPE = 88 %', fontsize = 14)
plt.text(0.05, 1.05, 'RMSE = 0.76 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_443_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443_original, kd_pred_443_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured [1/m]', fontsize = 16)
plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(246)
plt.text(0.05, 1.4, 'f)', fontsize = 16)
plt.text(0.05, 1.25, eq_492_original, fontsize = 14)
plt.text(0.05, 1.15, 'R² = 0.69; MAPE = 84 %', fontsize = 14)
plt.text(0.05, 1.05, 'RMSE = 0.50 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_492_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492_original, kd_pred_492_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured [1/m]', fontsize = 16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(247)
plt.text(0.05, 1.4, 'g)', fontsize = 16)
plt.text(0.05, 1.25, eq_560_original, fontsize = 14)
plt.text(0.05, 1.15, 'R² = 0.69; MAPE = 79 %', fontsize = 14)
plt.text(0.05, 1.05, 'RMSE = 0.30 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_560_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560_original, kd_pred_560_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured [1/m]', fontsize =  16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

plt.subplot(248)
plt.text(0.05, 1.4, 'h)', fontsize = 16)
plt.text(0.05, 1.25, eq_665_original, fontsize = 14)
plt.text(0.05, 1.15, 'R² = 0.32; MAPE = 72 %', fontsize = 14)
plt.text(0.05, 1.05, 'RMSE = 0.48 [1/m]; n = 20', fontsize = 14)
plt.plot(ref_665_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665_original, kd_pred_665_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured [1/m]', fontsize = 16)
#plt.ylabel('Kd-SA MSI [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.show()

###############################################################################

ref_443 = fun_reg(pd.Series(np.arange(0.7, 1.25, 0.01)), fit_443[0], fit_443[1])
ref_443 = ref_443.rename(pd.Series(np.arange(0.7, 1.25, 0.01)), axis = 'rows')

ref_492 = fun_reg(pd.Series(np.arange(0.45, 0.9, 0.01)), fit_492[0], fit_492[1])
ref_492 = ref_492.rename(pd.Series(np.arange(0.45, 0.9, 0.01)), axis = 'rows')

ref_560 = fun_reg(pd.Series(np.arange(0.3, 0.6, 0.01)), fit_560[0], fit_560[1])
ref_560 = ref_560.rename(pd.Series(np.arange(0.3, 0.6, 0.01)), axis = 'rows')

ref_665 = fun_reg(pd.Series(np.arange(0.6, 0.8, 0.01)), fit_665[0], fit_665[1])
ref_665 = ref_665.rename(pd.Series(np.arange(0.6, 0.8, 0.01)), axis = 'rows')

ref_443_original = fun_reg(pd.Series(np.arange(0.7, 1.25, 0.01)), fit_443_original[0], fit_443_original[1])
ref_443_original = ref_443_original.rename(pd.Series(np.arange(0.7, 1.25, 0.01)), axis = 'rows')

ref_492_original = fun_reg(pd.Series(np.arange(0.45, 0.9, 0.01)), fit_492_original[0], fit_492_original[1])
ref_492_original = ref_492_original.rename(pd.Series(np.arange(0.45, 0.9, 0.01)), axis = 'rows')

ref_560_original = fun_reg(pd.Series(np.arange(0.3, 0.6, 0.01)), fit_560_original[0], fit_560_original[1])
ref_560_original = ref_560_original.rename(pd.Series(np.arange(0.3, 0.6, 0.01)), axis = 'rows')

ref_665_original = fun_reg(pd.Series(np.arange(0.6, 0.8, 0.01)), fit_665_original[0], fit_665_original[1])
ref_665_original = ref_665_original.rename(pd.Series(np.arange(0.6, 0.8, 0.01)), axis = 'rows')

#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(111)
plt.text(1, 0.90, eq_492, fontsize = 16, color = '#191970')
plt.text(1, 0.82, 'R² = 0.70; MAPE = 14 %', fontsize = 16, color = '#191970')
plt.text(1, 0.74, 'RMSE = 0.10 [1/m]; n = 20', fontsize = 16, color = '#191970')
plt.text(1, 0.66, eq_492, fontsize = 16, color = '#00FFFF')
plt.text(1, 0.58, 'R² = 0.70; MAPE = 14 %', fontsize = 16, color = '#00FFFF')
plt.text(1, 0.50, 'RMSE = 0.10 [1/m]; n = 20', fontsize = 16, color = '#00FFFF')
plt.text(1, 0.42, eq_560, fontsize = 16, color = '#00FF00')
plt.text(1, 0.34, 'R² = 0.56; MAPE = 38 %', fontsize = 16, color = '#00FF00')
plt.text(1, 0.26, 'RMSE = 0.15 [1/m]; n = 20', fontsize = 16, color = '#00FF00')
plt.text(1, 0.18, eq_665, fontsize = 16, color = '#FF0000')
plt.text(1, 0.10, 'R² = 0.01; MAPE = 39 %', fontsize = 16, color = '#FF0000')
plt.text(1, 0.02, 'RMSE = 0.28 [1/m]; n = 20', fontsize = 16, color = '#FF0000')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.plot(ref_443, color = '#191970', ls = '--')
plt.plot(ref_492, color = '#00FFFF', ls = '--')
plt.plot(ref_560, color = '#00FF00', ls = '--')
plt.plot(ref_665, color = '#FF0000', ls = '--')
plt.xlabel('Kd-measured [1/m]', fontsize = 16)
plt.ylabel('Kd-SA [1/m]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)
plt.legend(['1:1', '443 nm', '492 nm', '560 nm', '665 nm'], loc = 0, fontsize = 16) 
#plt.title('Estatísticas Globais MSI', fontsize = 14)
plt.show()

###############################################################################

plt.figure(figsize = (10, 10))

plt.subplot(121)
plt.text(0.1, 1.4, 'a)', fontsize = 18)
plt.text(0.2, 1.4, eq_443, fontsize = 16, color = '#191970')
plt.text(0.2, 1.32, 'R² = 0.44; MAPE = 13.49%', fontsize = 16, color = '#191970')
plt.text(0.2, 1.24, 'RMSE = 0.14 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#191970')
plt.text(0.95, 0.66, eq_492, fontsize = 16, color = '#00FFFF')
plt.text(0.95, 0.58, 'R² = 0.70; MAPE = 14.09%', fontsize = 16, color = '#00FFFF')
plt.text(0.95, 0.50, 'RMSE = 0.10 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#00FFFF')
plt.text(0.95, 0.42, eq_560, fontsize = 16, color = '#00FF00')
plt.text(0.95, 0.34, 'R² = 0.56; MAPE = 38.32%', fontsize = 16, color = '#00FF00')
plt.text(0.95, 0.26, 'RMSE = 0.15 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#00FF00')
plt.text(0.95, 0.18, eq_665, fontsize = 16, color = '#FF0000')
plt.text(0.95, 0.10, 'R² = 0.01; MAPE = 38.82%', fontsize = 16, color = '#FF0000')
plt.text(0.95, 0.02, 'RMSE = 0.28 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#FF0000')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.plot(ref_443, color = '#191970', ls = '--')
plt.plot(ref_492, color = '#00FFFF', ls = '--')
plt.plot(ref_560, color = '#00FF00', ls = '--')
plt.plot(ref_665, color = '#FF0000', ls = '--')
plt.xlabel('Kd-measured [${m}^{-1}$]', fontsize = 18)
plt.ylabel('Kd-SA [${m}^{-1}$]', fontsize = 18)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
#plt.legend(['1:1', '443 nm', '492 nm', '560 nm', '665 nm'], loc = 0, fontsize = 16) 
#plt.title('Estatísticas Globais MSI', fontsize = 14)
plt.show()

plt.subplot(122)
plt.text(0.1, 1.4, 'b)', fontsize = 18)
plt.text(0.2, 1.4, eq_443_original, fontsize = 16, color = '#191970')
plt.text(0.2, 1.32, 'R² = 0.56; MAPE = 88.18%', fontsize = 16, color = '#191970')
plt.text(0.2, 1.24, 'RMSE = 0.76 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#191970')
plt.text(0.2, 1.16, eq_492_original, fontsize = 16, color = '#00FFFF')
plt.text(0.2, 1.08, 'R² = 0.69; MAPE = 84.36%', fontsize = 16, color = '#00FFFF')
plt.text(0.2, 1.00, 'RMSE = 0.50 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#00FFFF')
plt.text(0.95, 0.80, eq_560_original, fontsize = 16, color = '#00FF00')
plt.text(0.95, 0.72, 'R² = 0.69; MAPE = 78.75%', fontsize = 16, color = '#00FF00')
plt.text(0.95, 0.64, 'RMSE = 0.30 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#00FF00')
plt.text(0.95, 0.56, eq_665_original, fontsize = 16, color = '#FF0000')
plt.text(0.95, 0.48, 'R² = 0.32; MAPE = 72.18%', fontsize = 16, color = '#FF0000')
plt.text(0.95, 0.40, 'RMSE = 0.49 [${m}^{-1}$]; n = 20', fontsize = 16, color = '#FF0000')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443_original, kd_pred_443_original, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref_492_original, kd_pred_492_original, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref_560_original, kd_pred_560_original, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref_665_original, kd_pred_665_original, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.plot(ref_443_original, color = '#191970', ls = '--')
plt.plot(ref_492_original, color = '#00FFFF', ls = '--')
plt.plot(ref_560_original, color = '#00FF00', ls = '--')
plt.plot(ref_665_original, color = '#FF0000', ls = '--')
plt.xlabel('Kd-measured [${m}^{-1}$]', fontsize = 18)
#plt.ylabel('Kd-SA [${m}^{-1}$]', fontsize = 16)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)
plt.legend(['1:1', '443 nm', '492 nm', '560 nm', '665 nm'], loc = 1, fontsize = 16) 
#plt.title('Estatísticas Globais MSI', fontsize = 14)
plt.show()

###############################################################################

x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(221)
plt.text(0.05, 1.4, 'a)', fontsize = 12)
plt.text(0.4, 0.25, eq_443, fontsize = 12)
plt.text(0.4, 0.15, 'R² = 0.4359; MAPE = 13.49%', fontsize = 12)
plt.text(0.4, 0.05, 'RMSE = 0.1399 m^-1; n = 20', fontsize = 12)
plt.plot(ref_443, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
plt.ylabel('Kd-SA MSI', fontsize = 12)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.show()

plt.subplot(222)
plt.text(0.05, 1.4, 'b)', fontsize = 12)
plt.text(0.4, 0.25, eq_492, fontsize = 12)
plt.text(0.4, 0.15, 'R² = 0.6956; MAPE = 14.09%', fontsize = 12)
plt.text(0.4, 0.05, 'RMSE = 0.0951 m^-1; n = 20', fontsize = 12)
plt.plot(ref_492, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([])
plt.yticks([])
plt.legend(['Fit', '1:1', 'MSI/19'], loc = 1, fontsize = 12) 
plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.show()

plt.subplot(223)
plt.text(0.05, 1.4, 'c)', fontsize = 12)
plt.text(0.4, 0.25, eq_560, fontsize = 12)
plt.text(0.4, 0.15, 'R² = 0.5601; MAPE = 38.32%', fontsize = 12)
plt.text(0.4, 0.05, 'RMSE = 0.1458 m^-1; n = 20', fontsize = 12)
plt.plot(ref_560, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.ylabel('Kd-SA MSI', fontsize = 12)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.show()

plt.subplot(224)
plt.text(0.05, 1.4, 'd)', fontsize = 12)
plt.text(0.4, 0.25, eq_665, fontsize = 12)
plt.text(0.4, 0.15, 'R² = 0.0095; MAPE = 38.82%', fontsize = 12)
plt.text(0.4, 0.05, 'RMSE = 0.2774 m^-1; n = 20', fontsize = 12)
plt.plot(ref_665, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = 'c')
#plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)
plt.show()

plt.subplot(245)
plt.text(0.05, 1.4, 'e)', fontsize = 12)
plt.text(0.2, 1.4, eq_443_original, fontsize = 9)
plt.text(0.2, 1.3, 'R² = 0.5651; MAPE = 88.18%', fontsize = 9)
plt.text(0.2, 1.2, 'RMSE = 0.7633 m^-1; n = 20', fontsize = 9)
plt.plot(ref_443_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443_original, kd_pred_443_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(246)
plt.text(0.05, 1.4, 'f)', fontsize = 12)
plt.text(0.2, 1.4, eq_492_original, fontsize = 9)
plt.text(0.2, 1.3, 'R² = 0.6864; MAPE = 84.36%', fontsize = 9)
plt.text(0.2, 1.2, 'RMSE = 0.4990 m^-1; n = 20', fontsize = 9)
plt.plot(ref_492_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492_original, kd_pred_492_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(247)
plt.text(0.05, 1.4, 'g)', fontsize = 12)
plt.text(0.2, 1.4, eq_560_original, fontsize = 9)
plt.text(0.2, 1.3, 'R² = 0.6888, MAPE = 78.75%', fontsize = 9)
plt.text(0.2, 1.2, 'RMSE = 0.2958 m^-1; n = 20', fontsize = 9)
plt.plot(ref_560_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560_original, kd_pred_560_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()

plt.subplot(248)
plt.text(0.05, 1.4, 'h)', fontsize = 12)
plt.text(0.2, 1.4, eq_665_original, fontsize = 9)
plt.text(0.2, 1.3, 'R² = 0.3164; MAPE = 72.18%', fontsize = 9)
plt.text(0.2, 1.2, 'RMSE = 0.4787 m^-1; n = 20', fontsize = 9)
plt.plot(ref_665_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665_original, kd_pred_665_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured', fontsize = 9)
#plt.ylabel('Kd-SA MSI', fontsize = 9)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([])
plt.rc('xtick', labelsize = 9)
plt.rc('ytick', labelsize = 9)
plt.show()


###############################################################################
