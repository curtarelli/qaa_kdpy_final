'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------


'''
##  Importando pacotes básicos necessários para o código
import os
import gc
import glob
from osgeo import gdal
from gdalconst import *
import pandas as pd
import geopandas as gpd
import numpy as np
import numpy.ma as ma
import math as mt
import datetime as dt
import scipy.ndimage

from shapely.geometry import mapping
import rasterio as rio
from rasterio.plot import plotting_extent
from rasterio.mask import mask
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep

import seaborn as sns
from matplotlib import pyplot as plt

import xml.etree.ElementTree as ET
import requests

##  Modelo QAA/Kd para se aplicar às bandas do Sentinel 2A msi
from qaapy_kdpy.QAA_models.Espectral.Qaa_v6_model import Qaa_v6_model_msi

##  Funções auxiliares para carregamento dos dados
from qaapy_kdpy.utils.image_app import *

# Prettier plotting with seaborn
sns.set(font_scale=1.5)

##  Diretório onde se encontram as imagens
dir_i = r'C:\Imagens'
dir_a = r'D:\prog\master\qaapy_kdpy'

os.chdir(dir_i)
list_i = glob.glob('20*')

# i = list_i[2]

for i in list_i:
    print('Entrou em ' + i)
    dir_ii = dir_i + '\\' + i
    os.chdir(dir_ii)
    list_ii = glob.glob('*')
    
    m_set = np.full([10980, 10980], 0)
    
    count_ii = 0
    
    # ii = list_ii[0]
    
    for ii in list_ii:
        count_ii =+ 1
        
        dir_iii = dir_ii + '\\' + ii + '\\crop'
        
        os.chdir(dir_iii)
        print('Entrou em ' + i + '\\' + ii)
        
        dir_b11 = dir_iii + '\\' + glob.glob('*B11.tif')[0]
        dir_mask = dir_iii + '\\' + glob.glob('crop_L1C*')[0]
        
        date_ref = glob.glob('crop_L1C*')[0].split('.tif')[0].split('_')[-1]
        
        ##  Carregando as bandas de interesse (B1 a B5)
        b1, b2, b3, b4, b5 = open_msi_b1_to_b5(dir_iii)
        
        ##  Carregando banda B11 para remoção de glint;
        b11 = open_tif(dir_b11)
        
        ##  Carregando Mascara;
        bmask = open_tif(dir_mask)
        bmask[bmask == 5] = 0
        bmask = ma.masked_where(bmask > 0, bmask)
        
        # bmask = np.clip(bmask, 4, 5)
        
        '''
        plt.figure(figsize = (10, 10))
        plt.plot()
        plt.title('')
        plt.imshow(bmask, cmap='gray')
        plt.show()
        '''
        
        '''
        plt.figure(figsize = (10, 10))
        plt.plot()
        plt.title('')
        plt.imshow(b1, cmap='gray')
        plt.show()
        '''
        
        '''
        plt.figure(figsize = (10, 10))
        plt.plot()
        plt.title('')
        plt.imshow(b2, cmap='gray')
        plt.show()
        '''
        
        '''
        plt.figure(figsize = (10, 10))
        plt.plot()
        plt.title('')
        plt.imshow(b5, cmap='gray')
        plt.show()
        '''
        
        ##  Reamostragem espacial das bandas B1 (60m / 6) e B5 (20 m / 2) [RE = 10 metros];
        ##  mode = ['reflect', 'constant', 'nearest', 'mirror', 'wrap'];
        b1 = scipy.ndimage.zoom(b1, 6, order = 0, mode = 'nearest')
        b5 = scipy.ndimage.zoom(b5, 2, order = 0, mode = 'nearest')
        b11 = scipy.ndimage.zoom(b11, 2, order = 0, mode = 'nearest')
        bmask = scipy.ndimage.zoom(bmask, 2, order = 0, mode = 'nearest')
        print('Reamostragem ' + i + '\\' + ii + ' [OK!]')
        
        b1 = np.delete(b1, [0, 1, 2, 3, 4, 5, 6], axis = 0)
        b1 = np.delete(b1, [0, 1, 2, 3, 4], axis = 1)
        b5 = np.delete(b5, [0], axis = 0)
        b5 = np.delete(b5, [0], axis = 1)
        b11 = np.delete(b11, [0], axis = 0)
        b11 = np.delete(b11, [0], axis = 1)
        bmask = np.delete(bmask, [0], axis = 0)
        bmask = np.delete(bmask, [0], axis = 1)
        print('Ajuste pixels ' + i + '\\' + ii + ' [OK!]')
        
        b1 = b1 - b11
        b2 = b2 - b11
        b3 = b3 - b11
        b4 = b4 - b11
        b5 = b5 - b11
        print('Remoção de glint ' + i + '\\' + ii + ' [OK!]')
        
        b1 = b1 / mt.pi
        b2 = b2 / mt.pi
        b3 = b3 / mt.pi
        b4 = b4 / mt.pi
        b5 = b5 / mt.pi
        print('Cálculo Rrs ' + i + '\\' + ii + ' [OK!]')
        
        b1 = ma.masked_array(b1, bmask)
        b2 = ma.masked_array(b2, bmask)
        b3 = ma.masked_array(b3, bmask)
        b4 = ma.masked_array(b4, bmask)
        b5 = ma.masked_array(b5, bmask)
        print('Mascara de nuvens e sombras ' + i + '\\' + ii + ' [OK!]')
        
        '''
        plt.figure(figsize = (10, 10))
        plt.plot()
        plt.title('')
        plt.imshow(b1, cmap='gray')
        plt.show()
        '''
        
        os.chdir(dir_a)
        print('Entrando no diretório do QAA-Kd-Py')
        
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
        
        print('Carregado aw e bbw [OK!]')
              
        ##  Latitude central do reservatório (imagem)
        lat = mt.radians(-18.6087777777778)
        
        year = int(date_ref[:4])
        month = int(date_ref[4:6])
        day = int(date_ref[6:8])
        hour = int(date_ref[9:11])
        minute = int(date_ref[11:13])
        second = int(date_ref[13:])
        print('Separação de: lat, ano, mês, dia, hora, minuto e segundo de passagem do satélite [OK!]')
        
        ##  Data da imagem com a hora de aquisição
        date = dt.datetime(year, month, day, hour, minute, second)
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
        print('Cálculo do ângulo zenital solar [OK!]')
        
        ##
        kd_msi_b2 = Qaa_v6_model_msi(b2, 492, b1, b2, b3, b4, b5, 560, aw, bbw, theta_s)
        print('Aplicação do algoritmo de Kd [OK!]')
        
        m_set =+ kd_msi_b2
        print('Adição do resultado em Matriz para média mensal [OK!]')
        
        del kd_msi_b2
        gc.collect()
        print('Arquivo deletado da memória [OK!]')
        
    mean_set = m_set / count_ii
    print('Cálculo da matríz da média mensal de Kd para ' + i + ' [OK!]')
    
    os.chdir(dir_iii)
    
    dir_save = dir_ii + '\\kd_' + i + 'B02.tif'
    
    ref_band = dir_iii + '\\' + glob.glob('*B02.tif')[0]
    dataset_de_referencia = gdal.Open(ref_band, GA_ReadOnly)
    salvar_banda(mean_set, dir_save, dataset_de_referencia)
    print('Salvo arquivo .TIF com média mensal de Kd para ' + i + ' [OK!]')
    
    del mean_set
    gc.collect()
    print('Arquivo de média mensal deletado da memória [OK!]')
    