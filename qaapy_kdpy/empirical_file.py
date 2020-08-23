'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Funções para a leitura dos dados de entrada a serem utilizados no modelo QAA,
salvar e carregar resultados. 
'''
##  Importando pacotes básicos necessários para este pacote.
import pickle as pk
import pandas as pd
from collections import namedtuple
from qaapy_kdpy.empirical_definition import create_obj_data

###############################################################################
def read_files(args: namedtuple):
    '''
    Função para realizar a leitura dos dados definidos pelo usuário na chamadas de
    argumentos.

    ----------
    Parameters
    ----------
    args [Object]
        Objeto Args com as definições de dados a serem carregados.
    
    ------
    Raises
    ------
    Exception
        Caso seja passado um objeto de argumentos errado com nome diferente do
        criado pela função qaapy_kdpy.create_args um aviso pede ao usuário para
        que corrija esse problema.
        
    -------
    Returns
    -------
    data_obj [Object]
        Objeto com Data Frames contendo dados usados nas funções do modelo QAA/Kd
        alocados pela função read_files através dos argumentos contendo as strings
        do caminho do diretório para cada planilha '.xlsx' de dados.
    '''
    # Verificando o tipo do argumento
    if args.__qualname__ != 'Args':
        raise Exception('É necessário que o argumento seja um objeto do tipo Args.' + \
                                ' Para isso utilize a função qaapy.create_args')
    else:
        data_object = create_obj_data()
    
    # Fazendo a leitura dos dados com as informações presente no 'args'
    data_object.acs = pd.read_excel(
                                        args.acs['file_name'], 
                                        sheet_name = args.acs['sheet_name'],
                                        index_col = 0)
    
    data_object.Rrs = pd.read_excel(
                                        args.Rrs['file_name'],
                                        sheet_name = args.Rrs['sheet_name'],
                                        header = 0,
                                        index_col = 0)
       
    data_object.aw = pd.read_excel(
                                        args.aw['file_name'],
                                        sheet_name = args.aw['sheet_name'],
                                        header = 0,
                                        index_col = 0)
    
    data_object.bbw = pd.read_excel(
                                        args.bbw['file_name'],
                                        sheet_name = args.bbw['sheet_name'],
                                        header = 0,
                                        index_col = 0)
    
    data_object.bb = pd.read_excel(
                                        args.bb['file_name'],
                                        sheet_name = args.bb['sheet_name'],
                                        header = 0,
                                        index_col = 0)
    
    data_object.bbp = pd.read_excel(
                                        args.bbp['file_name'],
                                        sheet_name = args.bbp['sheet_name'],
                                        header = 0,
                                        index_col = 0)
    
    return data_object

###############################################################################
def export_result(data_to_save, filename):
    '''
    Função para exportar os dados do dicionário em formato '.pickle' para o
    diretório escolhido pelo usuário.
    
    ----------
    Parameters
    ----------
    data_to_save [Dictionary]
        Dicionário de dados contendo infromações resultados e valores de entrada.
    filename [String]
        String do caminho do diretório para salvar os resultados contendo o nome
        desejado para o arquivo seguido do sufixo '.pickle'.
    
    -------
    Returns
    -------
    Salva o arquivo '.pickle' com todos os dados no diretório desejado.
    '''
    with open(filename, 'wb') as file:
        pk.dump(data_to_save, file)

###############################################################################
def load_data(filename):
    '''
    Função para carregar os dados '.pickle' após salvos, para isso carregue a 
    função em uma variável.

    ----------
    Parameters
    ----------
    filename [String]
        String do caminho do diretório para carregar os resultados contendo o nome
        desejado para o arquivo seguido do sufixo '.pickle'.
    
    -------
    Returns
    -------
    loaded_data [Dictionary]
        Dicionário de dados contendo infromações resultados e valores carregados.
    '''
    with open(filename, 'rb') as file:
        return pk.load(file)