'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Funções para carregar os argumentos contendo os caminhos para os diretórios,
arquivos e abas contendo dados a serem utilizados no modelo QAA. Além disso há
uma função para criar o objeto de dados.

Argumentos: acs, Rrs, aw, bbw, coleta, Kd_ZEU, Kd_3m e Kd_6m
'''
##  Importando pacotes básicos necessários para o script
from collections import namedtuple

# Função para criar o objeto de argumentos a ser utilizado
create_args = lambda: namedtuple('Args', ['acs', 'Rrs', 'aw', 'bbw', 'bb', 'coleta', 'Kd_ZEU', 'Kd_3m', 'Kd_6m'])

# Função para criar objeto com todos os dados definidos pelo usuário
create_obj_data = lambda: namedtuple('Data', ['acs', 'Rrs', 'aw', 'bbw', 'bb', 'coleta', 'Kd_ZEU', 'Kd_3m', 'Kd_6m'])
