'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é um pacote de ferramentas uteis para aplicação do modelo QAA e Kd
nas imagens (MSI, e outras), extrair pontos de interesse para análises estatísticas
e plotagem dos resultados de Kd da imagem.
'''
##  Importando pacotes necessários para as funções
from osgeo import gdal
from gdalconst import *
import os
import glob

gdal.UseExceptions()

##  Função para abrir um arquivo ".tif"
def open_tif(file_name):
    '''
    Função para abrir único arquivo '.tif' como uma matriz (Array), alocando ele
    na variável escolhida pelo usuário.
    
    ----------
    Parameters
    ----------
    file_name [String]
        String contendo o caminho para o diretório da imagem, incluindo o nome da
        imagem e a extensão '.tif' no final.
        
    -------
    Returns
    -------
    dataset [Matriz]
        Retorna o arquivo lido na forma de matriz, útil para a manipulação e
        aplicação de cálculos.
    '''
    
    try:
        dataset = gdal.Open(file_name, GA_ReadOnly).ReadAsArray()
        print('Arquivo aberto com sucesso.')
        
        return dataset
    
    except:
        print('Erro na abertura do arquivo.')

##  Função para abrir todas as bandas msi inclusas em uma pasta (em batch)
def open_msi_batch(dir_tif):
    '''
    Função para abrir os arquivos '.TIF' de cada umas das bandas Sentinel 2 MSI
    na forma de matriz (Array), alocando elas nas variáveis escolhidas pelo
    usuário.
    
    ----------
    Parameters
    ----------
    dir_tif [String]
        String contendo o caminho para o diretório das imagens, sem conter a extensão
        dos arquivos.
    
    -------
    Returns
    -------
    dataset [Matriz]
        Retorna o arquivo lido na forma de matriz, útil para a manipulação e
        aplicação de cálculos.
    '''
    dir_i = os.getcwd()
    
    os.chdir(dir_tif)
    
    b1 = glob.glob('*B01.TIF')[0]
    b2 = glob.glob('*B02.TIF')[0]
    b3 = glob.glob('*B03.TIF')[0]
    b4 = glob.glob('*B04.TIF')[0]
    b5 = glob.glob('*B05.TIF')[0]
    b6 = glob.glob('*B06.TIF')[0]
    b7 = glob.glob('*B07.TIF')[0]
    b8 = glob.glob('*B08.TIF')[0]
    b9 = glob.glob('*B09.TIF')[0]
    b10 = glob.glob('*B10.TIF')[0]
    b11 = glob.glob('*B11.TIF')[0]
    b12 = glob.glob('*B12.TIF')[0]
    
    dataset_b1 = open_tif(b1)
    dataset_b2 = open_tif(b2)
    dataset_b3 = open_tif(b3)
    dataset_b4 = open_tif(b4)
    dataset_b5 = open_tif(b5)
    dataset_b6 = open_tif(b6)
    dataset_b7 = open_tif(b7)
    dataset_b8 = open_tif(b8)
    dataset_b9 = open_tif(b9)
    dataset_b10 = open_tif(b10)
    dataset_b11 = open_tif(b11)
    dataset_b12 = open_tif(b12)
    
    os.chdir(dir_i)
    
    return dataset_b1, dataset_b2, dataset_b3, dataset_b4, dataset_b5, dataset_b6, dataset_b7, dataset_b8, dataset_b9, dataset_b10, dataset_b11, dataset_b12

##  Função para abrir bandas msi B1 a B5 inclusas em uma pasta (em batch)
def open_msi_b1_to_b5(dir_tif):
    '''
    Função para abrir os arquivos '.TIF' de cada umas das bandas Sentinel 2 MSI
    usandas nos calculos (B1 a B5) na forma de matriz (Array), alocando elas nas
    variáveis escolhidas pelo usuário.
    
    ----------
    Parameters
    ----------
    dir_tif [String]
        String contendo o caminho para o diretório das imagens, sem conter a extensão
        dos arquivos.
    
    -------
    Returns
    -------
    dataset [Matriz]
        Retorna o arquivo lido na forma de matriz, útil para a manipulação e
        aplicação de cálculos.
    '''
    dir_i = os.getcwd()
    
    os.chdir(dir_tif)
    
    b1 = glob.glob('*B01.TIF')[0]
    b2 = glob.glob('*B02.TIF')[0]
    b3 = glob.glob('*B03.TIF')[0]
    b4 = glob.glob('*B04.TIF')[0]
    b5 = glob.glob('*B05.TIF')[0]
    
    dataset_b1 = open_tif(b1)
    dataset_b2 = open_tif(b2)
    dataset_b3 = open_tif(b3)
    dataset_b4 = open_tif(b4)
    dataset_b5 = open_tif(b5)
    
    os.chdir(dir_i)    
    
    return dataset_b1, dataset_b2, dataset_b3, dataset_b4, dataset_b5

def salvar_banda(matriz_de_pixels, nome_do_arquivo, dataset_de_referencia):
    # obter metadados
    linhas = dataset_de_referencia.RasterYSize
    colunas = dataset_de_referencia.RasterXSize
    bandas = 1
    # definir driver
    driver = gdal.GetDriverByName('GTiff')
    # copiar tipo de dados da banda já existente
    data_type = dataset_de_referencia.GetRasterBand(1).DataType
    # criar novo dataset
    dataset_output = driver.Create(nome_do_arquivo, colunas, linhas, bandas, data_type)
    # copiar informações espaciais da banda já existente
    dataset_output.SetGeoTransform(dataset_de_referencia.GetGeoTransform())
    # copiar informações de projeção
    dataset_output.SetProjection(dataset_de_referencia.GetProjectionRef())
    # escrever dados da matriz NumPy na banda
    dataset_output.GetRasterBand(1).WriteArray(matriz_de_pixels)
    # salvar valores
    dataset_output.FlushCache()
    # fechar dataset
    dataset_output = None
