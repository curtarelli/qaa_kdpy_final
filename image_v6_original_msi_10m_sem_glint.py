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
from qaapy_kdpy.QAA_models.Espectral.Qaa_v6_model import Qaa_v6_model_msi_original

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
kd_msi_b1 = Qaa_v6_model_msi_original(b1, 443, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b1, r'.\03_Results_image\QAAv6_msi_original_2019_10m_sem_glint_443nm.tif', dataset_de_referencia)
del kd_msi_b1
gc.collect()

kd_msi_b2 = Qaa_v6_model_msi_original(b2, 492, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b2, r'.\03_Results_image\QAAv6_msi_original_2019_10m_sem_glint_492nm.tif', dataset_de_referencia)
del kd_msi_b2
gc.collect()

kd_msi_b3 = Qaa_v6_model_msi_original(b3, 560, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b3, r'.\03_Results_image\QAAv6_msi_original_2019_10m_sem_glint_560nm.tif', dataset_de_referencia)
del kd_msi_b3
gc.collect()

kd_msi_b4 = Qaa_v6_model_msi_original(b4, 665, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b4, r'.\03_Results_image\QAAv6_msi_original_2019_10m_sem_glint_665nm.tif', dataset_de_referencia)
del kd_msi_b4
gc.collect()

kd_msi_b5 = Qaa_v6_model_msi_original(b5, 704, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
salvar_banda(kd_msi_b5, r'.\03_Results_image\QAAv6_msi_original_2019_10m_sem_glint_704nm.tif', dataset_de_referencia)
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
qaa_v6_fit19_665 = load_data(r'.\Results_Espectral\QAAv6_original_TRM_2019_665_Ajustado.pickle')

kd_ref = qaa_v6_fit19_665['Kd_ZEU'].copy()
kd_pred = pd.read_excel(r'.\03_Results_Image\TRM_2019_QAAv6_msi.xlsx',
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
kd_ref_443_original = kd_ref.loc[443]                     ##  Dado de referência (y_true)
kd_pred_443_original = kd_pred['Kd_QAAv6_original_10m_s-glint_443']        ##  Dados simulados (y_pred)

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
eq_443_original = '443 nm = 0.16x - 0.04'
stat_443_original = 'R² = 0.57; MAPE = 88.18%; RMSE = 0.76 [${m}^{-1}$]; n = 20'

##  Plotando dispesão dos dados e os resultados da regressão
plot_ajuste(kd_ref_443_original,
            kd_pred_443_original,
            ref_443_original,
            eq_443_original,
            stat_443_original,
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
kd_ref_492_original = kd_ref.loc[492]
kd_pred_492_original = kd_pred['Kd_QAAv6_original_10m_s-glint_492']

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

eq_492_original = '492 nm = 0.23x - 0.04'
stat_492_original = 'R² = 0.69; MAPE = 84.36%; RMSE = 0.50 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_492_original,
            kd_pred_492_original,
            ref_492_original,
            eq_492_original,
            stat_492_original,
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
kd_ref_560_original = kd_ref.loc[560]
kd_pred_560_original = kd_pred['Kd_QAAv6_original_10m_s-glint_560']

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

eq_560_original = '560 nm = 0.32x - 0.04'
stat_560_original = 'R² = 0.69; MAPE = 78.75%; RMSE = 0.30 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_560_original,
            kd_pred_560_original,
            ref_560_original,
            eq_560_original,
            stat_560_original,
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
kd_ref_665_original = kd_ref.loc[665]
kd_pred_665_original = kd_pred['Kd_QAAv6_original_10m_s-glint_665']

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

eq_665_original = '665 nm = 0.55x - 0.18'
stat_665_original = 'R² = 0.32; MAPE = 72.18%; RMSE = 0.48 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_665_original,
            kd_pred_665_original,
            ref_665_original,
            eq_665_original,
            stat_665_original,
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
kd_ref_704_original = kd_ref.loc[704]
kd_pred_704_original = kd_pred['Kd_QAAv6_original_10m_s-glint_704']

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

eq_704_original = '704 nm = 0.66x - 0,30'
stat_704_original = 'R² = 0.34; MAPE = 70.34%; RMSE = 0.59 [${m}^{-1}$]; n = 20'

plot_ajuste(kd_ref_704_original,
            kd_pred_704_original,
            ref_704_original,
            eq_704_original,
            stat_704_original,
            [0, 2.5], [0, 2.5],
            'Kd-measured (704 nm)',
            'Kd-SA MSI (704 nm)',
            'Kd-measured vs. Kd-SA MSI',
            'Fit',
            'Stations',
            leg_loc = 0,
            x1 = True)

###############################################################################
#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(221)
plt.text(0.1, 1.2, 'a)', fontsize = 12)
plt.text(0.2, 1.4, eq_443_original, fontsize = 10)
plt.text(0.2, 1.25, 'R² = 0.5651; MAPE = 88.18%', fontsize = 10)
plt.text(0.2, 1.1, 'RMSE = 0.7633 m^-1; n = 20', fontsize = 10)
plt.plot(ref_443_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443_original, kd_pred_443_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured (443 nm)', fontsize = 14)
plt.ylabel('Kd-SA MSI (443 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('443 nm', fontsize = 14)
plt.show()

plt.subplot(222)
plt.text(0.1, 1.2, 'b)', fontsize = 12)
plt.text(0.2, 1.4, eq_492_original, fontsize = 10)
plt.text(0.2, 1.25, 'R² = 0.6864; MAPE = 84.36%', fontsize = 10)
plt.text(0.2, 1.1, 'RMSE = 0.4990 m^-1; n = 20', fontsize = 10)
plt.plot(ref_492_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_492_original, kd_pred_492_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured (492 nm)', fontsize = 14)
plt.ylabel('Kd-SA MSI (492 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('492 nm', fontsize = 14)
plt.legend(['Fit', '1:1', 'MSI/19'], loc = 1, fontsize = 12) 
plt.show()

plt.subplot(223)
plt.text(0.1, 1.2, 'c)', fontsize = 12)
plt.text(0.2, 1.4, eq_560_original, fontsize = 10)
plt.text(0.2, 1.25, 'R² = 0.6888, MAPE = 78.75%', fontsize = 10)
plt.text(0.2, 1.1, 'RMSE = 0.2958 m^-1; n = 20', fontsize = 10)
plt.plot(ref_560_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_560_original, kd_pred_560_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured (560 nm)', fontsize = 14)
plt.ylabel('Kd-SA MSI (560 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('560 nm', fontsize = 14)
plt.show()

plt.subplot(224)
plt.text(0.1, 1.2, 'd)', fontsize = 12)
plt.text(0.2, 1.4, eq_665_original, fontsize = 10)
plt.text(0.2, 1.25, 'R² = 0.3164; MAPE = 72.18%', fontsize = 10)
plt.text(0.2, 1.1, 'RMSE = 0.4787 m^-1; n = 20', fontsize = 10)
plt.plot(ref_665_original, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_665_original, kd_pred_665_original, marker = '^', facecolors = 'none', edgecolors = 'c')
plt.xlabel('Kd-measured (665 nm)', fontsize = 14)
plt.ylabel('Kd-SA MSI (665 nm)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
#plt.title('665 nm', fontsize = 14)
plt.show()

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
eq_global = 'y = 0,9241x + 0,1886'
stat_global = 'R² = 0,5791; MAPE = 27,86 %; RMSE = 0,2047 m^-1; n = 100'

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
#plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
x = np.linspace(0, 10, 1000)

plt.figure(figsize = (10, 10))

plt.subplot(111)
plt.text(0.6, 0.3, eq_global, fontsize = 10)
plt.text(0.6, 0.2, 'R² = 0,5791; MAPE = 27,86%', fontsize = 10)
plt.text(0.6, 0.1, 'RMSE = 0,2047 m^-1; n = 100', fontsize = 10)
plt.plot(ref_global, 'r--')
plt.plot(x, x, color = 'k', ls = '--')
plt.scatter(kd_ref_443, kd_pred_443, marker = '^', facecolors = 'none', edgecolors = '#191970')
plt.scatter(kd_ref_492, kd_pred_492, marker = '^', facecolors = 'none', edgecolors = '#00FFFF')
plt.scatter(kd_ref_560, kd_pred_560, marker = '^', facecolors = 'none', edgecolors = '#00FF00')
plt.scatter(kd_ref_665, kd_pred_665, marker = '^', facecolors = 'none', edgecolors = '#FF0000')
plt.scatter(kd_ref_704, kd_pred_704, marker = '^', facecolors = 'none', edgecolors = '#FF1493')
plt.xlabel('Kd de campo (GLOBAL)', fontsize = 14)
plt.ylabel('Kd QAA MSI (GLOBAL)', fontsize = 14)
plt.xlim([0, 1.5])
plt.ylim([0, 1.5])
plt.xticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5], ['0,00', '0,25', '0,5', '0,75', '1,0', '1,25', '1,5'])
plt.yticks([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5], ['0,00', '0,25', '0,5', '0,75', '1,0', '1,25', '1,5'])
plt.rc('xtick', labelsize = 14)
plt.rc('ytick', labelsize = 14)
plt.legend(['Ajuste', '1:1', '443 nm', '492 nm', '560 nm', '665 nm', '704 nm' ], loc = 4, fontsize = 12) 
#plt.title('Estatísticas Globais MSI', fontsize = 14)
plt.show()

