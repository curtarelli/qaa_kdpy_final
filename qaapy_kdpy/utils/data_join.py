'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é um pacote de ferramentas para manipular dados em para organizar
banco de dados de entrada dos modelos QAA/Kd.
'''
##  Importanto pacotes básicos necessários para as funções.
import os
import glob
import pandas as pd
from pandas import ExcelWriter as Ew

###############################################################################
def data_join(dir_i, nome_o):
    '''
    Função para juntar diversas planilhas com estações separadas em apenas uma
    planilha contendo todos os dados compilados.
    
    ----------
    Parameters
    ----------
    dir_i [String]
        String contendo o diretório onde se encontram os dados.
    nome_o [String]
        Nome do arquivo de saida contendo todas as medidas em colunas distintas,
        contendo extensão '.xslx'.

    -------
    Returns
    -------
    df_o [Data Frame]
        Joined Data Frame saved as '.xlsx'.
    '''
    ##      
    dir_o = os.getcwd()
    
    os.chdir(dir_i)

    list_i = glob.glob('*.xlsx')
    
    df_o = pd.read_excel(list_i[0],
                          na_values = '',
                          sheet_name = 'Sheet1',
                          header = 0,
                          index_col = 0,
                          usecols = 0)
    
    for i in list_i:
        nome_i = i.split('.')[0]
        
        df_i = pd.read_excel(i,
                             na_values = '',
                             sheet_name = 'Sheet1',
                             header = 0,
                             index_col = 0)
        
        df_o[nome_i] = df_i.iloc[:, 0]
        
    writer = Ew(nome_o)
    df_o.to_excel(writer, 'Sheet1', index = True)
    writer.save()
    
    os.chdir(dir_o)
    
    return df_o

###############################################################################
def csv_to_ex(dir_i):
    '''
    Função para salvar arquivo CSV no dormato XLSX (Excel)
	
    ----------
    Parameters
    ----------
    dir_i [String]
        String contendo o diretório onde se encontram os dados.

    -------
    Returns
    -------
    df [Data Frame]
        Converted Data Frame saved as '.xlsx'.
    '''
    os.chdir(dir_i)
    
    list_i = glob.glob('*.csv')
    
    for i in list_i:
        nome = i.split('.')[0]
        
        df = pd.read_csv(i)
        
        writer = Ew(nome)
        df.to_excel(writer, 'Sheet1', index = True)
        writer.save()    
    
    os.chdir(r'..\..')
    
    return df
