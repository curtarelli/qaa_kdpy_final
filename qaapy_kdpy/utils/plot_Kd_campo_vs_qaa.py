'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Este script é um pacote de ferramentas uteis para plotagem de dados de entrada
e resultados dos algorítmos QAA e de estimativa de Kd.
'''
##  Importanto pacotes básicos necessários para as funções.
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

###############################################################################
def plot_kd_separado(qaa_data, kd_qaa, kd_campo, dir_o, formato):
    '''
    Função para plotar gráfico contendo dado espectral de campo e dado por banda
    calculado pelo QAA, por estação, para todas as estações contidas em uma
    pasta definida pelo usuário.
    
    ----------
    Parameters
    ----------
    qaa_data [Dictionary]
        Dicionário contendo os resultados do QAA, que pode ser um arquivo '.pickle'
        carregado.
    kd_qaa [String]
        Nome do resultado de Kd estimado pelo QAA. Por padrão 'Kd_QAA'.
    kd_campo [String]
        DNome do Kd de campo medido. Por padrão 'Kd_ZEU', 'Kd_3m' ou 'Kd_6m'.
    dir_o [String]
        Caminho do diretório para salvar as imagens.
    formato [String]
        String contendo o formato para salvar a image.
        
        Ex.: '.png', '.jpeg', '.tif', etc...

    -------
    Returns
    -------
    Plot ou imagem salva no diretório e formato escolhidos.
    '''
    dir_i = os.getcwd()
    
    os.chdir(dir_o)
    
    label_campo = kd_campo.split('_')[0] + ' ' + kd_campo.split('_')[1]
    label_qaa = kd_qaa.split('_')[0] + ' ' + kd_qaa.split('_')[1]
    
    for i in qaa_data[kd_qaa]:
        fig, ax = plt.subplots()
        
        ax.plot (qaa_data[kd_campo][i], 'r-', label = label_campo)
        ax.plot(qaa_data[kd_qaa][i], 'b-', label = label_qaa)
        
        plt.xlabel('Comprimento de Onda [nm]')
        plt.ylabel('Coef. Atenuação Difusa de Ed (Kd) [m-1]')
        plt.xlim(400, 700)
        plt.ylim(0, 2)
        plt.title(label_campo + ' Campo x ' + label_qaa + ' QAA - ' + i + ' (Lb0 = ' + str(qaa_data['dados']['lb0']) + ')')
        plt.legend()
        
        '''
        Nesta etapa deve-se escolher se será salva imagem automaticamente ou manualmente, para salvar automaticamente
        basta retirar os # antes das linhas de comando. Lembre-se de configurar o IDE para plotar imagens em Qt5
        '''
        fig.savefig('QAAv6_ajuste_2013_' + i + '_' + kd_campo + '_' + str(qaa_data['dados']['lb0']) + formato, dpi = 300)
    
    os.chdir(dir_i)

###############################################################################        
def plot_kd_nxm(qaa_data, kd_qaa, kd_campo, n_rows, n_cols, dir_o, formato):
    '''
    Plot em grade com todos os pontos, conforme numero de linhas e colunas
    selecionadas pelo usuário.
    
    ----------
    Parameters
    ----------
    qaa_data [Dictionary]
        Dicionário contendo os resultados do QAA, que pode ser um arquivo '.pickle'
        carregado.
    kd_qaa [String]
        Nome do resultado de Kd estimado pelo QAA. Por padrão 'Kd_QAA'.
    kd_campo [String]
        DNome do Kd de campo medido. Por padrão 'Kd_ZEU', 'Kd_3m' ou 'Kd_6m'.
    n_rows [Value]
        Número de linha para organizar o plot.
    n_cols [Value]
        Número de colunas para organizar o plot.
    dir_o [String]
        Caminho do diretório para salvar as imagens.
    formato [String]
        String contendo o formato para salvar a image.
        
        Ex.: '.png', '.jpeg', '.tif', etc...

    -------
    Returns
    -------
    Imagem salva no diretório e formato escolhidos.

    '''
    dir_i = os.getcwd()
    
    os.chdir(dir_o)
    
    label_campo = kd_campo.split('_')[0] + ' ' + kd_campo.split('_')[1]
    label_qaa = kd_qaa.split('_')[0] + ' ' + kd_qaa.split('_')[1]
    
    n = 1
    
    count_y = ((n_rows - 2) * n_cols) + 1
    count_x = n_rows * n_cols
    count_xy = ((n_rows - 1) * n_cols) + 1
    
    list_y = []
    list_x = []
    
    while count_y >= n_cols + 1:
        list_y.append(count_y)
        count_y -= n_cols

    while count_x > count_xy:
        list_x.append(count_x)
        count_x -= 1
    
    fig, ax = plt.subplots(nrows = n_rows, ncols = n_cols, 
                           sharex = True, sharey = True)
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.05, right=0.98, hspace=0.1, wspace=0.1)
    fig.text(0.5, 0.01, 'Comprimento de Onda [nm]', ha = 'center', fontsize = 12)
    fig.text(0.01, 0.5, 'Coef. de Atenuação de Ed (Kd) [m-1]', va = 'center', rotation = 'vertical', fontsize = 12)
    fig.text(0.5, 0.97, label_campo + ' Campo x ' + label_qaa + ' QAA' + ' (Lb0 = ' + str(qaa_data['dados']['lb0']) + ')', ha = 'center', fontsize = 12)

    
    for i in qaa_data[kd_qaa]:
        plt.subplot(n_rows, n_cols, n)
        plt.plot(qaa_data[kd_campo][i], 'r-', lw = 0.5, label = label_campo)
        plt.plot(qaa_data[kd_qaa][i], 'bx', ms = 2, mew = 0.5, label = label_qaa)
        plt.xlim(400, 700)
        plt.ylim(0, 2.5)
        plt.title(i, fontsize = 12, va = 'top')
        
        if n == 1:
            plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], fontsize = 10)
            plt.xticks([])
            plt.legend(fontsize = 10)
        
        elif n in list_y:
            plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], fontsize = 10)
            plt.xticks([])
            
        elif n in list_x:
            plt.xticks([400, 500, 600, 700], fontsize = 10)
            plt.yticks([])
            
        elif n == count_xy:
            plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5], fontsize = 10)
            plt.xticks([400, 500, 600, 700], fontsize = 10)

        else:
            plt.xticks([])
            plt.yticks([])
        
        n += 1

    '''
    Nesta etapa deve-se escolher se será salva imagem automaticamente ou manualmente, para salvar automaticamente
    basta retirar os # antes das linhas de comando. Lembre-se de configurar o IDE para plotar imagens em Qt5.
    
    Ao salvar via plot/Qt5 o usuário pode ajustar as propriedades do plot manualmente.
    '''
    # fig.savefig(qaa_data['dados']['equacao'].split('.')[0] + \
    #             qaa_data['dados']['equacao'].split('.')[1] + \
    #             '_TRM_' + kd_qaa + '_' + kd_campo + '_' + \
    #             str(qaa_data['dados']['lb0']) + formato, dpi = 300)
    
    os.chdir(dir_i)
    
###############################################################################
def plot_ajuste(x, y, ajuste, eq, a_r, xlims, ylims, xlabel, ylabel, title, leg_fit, leg_data, leg_loc, x1 = False):
    '''
    Plot para os ajuste dos Passo 2 e 4 do modelo QAA v6 proposto.
    
    ----------
    Parameters
    ----------
    x [Series]
        Valores para ajuste a serem alocados no eixo x.
    y [Series]
        Valores para ajuste a serem alocados no eixo y.
    ajuste [Data Frame]
        Valores da reta de melhor ajuste da função.
    eq [String]
        String contendo a equação de melhor ajuste segundo teste.
    r [String]
        String contendo valores de R² e n para o ajuste.
    xlims [List]
        Lista de valores para ajuste do eixo x no plot.
    ylims [List]
        Lista de valores para ajuste do eixo y no plot.
    xlabel [String]
        String com rótulo do eixo x.
    ylabel [String]
        String com rótulo do eixo y.
    title [String]
        String com Título do plot.
    leg_fit [String]
        String com texto da legenda da curva de ajuste.
    leg_data [String]
        String com texto da legenda da sipersão dos dados.

    -------
    Returns
    -------
    Plot sem salvar a imagem para ajuste fino e posterior salvamento.
    '''
    fig, ax = plt.subplots()
    fig.text(0.3, 0.2, eq)
    fig.text(0.3, 0.15, a_r)    
    
    if x1 == True:      
        plt.arrow(0, 0, 10, 10, color = 'k', ls = '--')
        
        '''
        A ideia que foi deixada para trás aqui era da adição de uma linha 1:1
        que permite o comando .legend().
        
        A função plot.arro() não possui o comando .legend().
        '''
        
        # plt.legend(['1:1'], loc = leg_loc)
        # fun_x1 = lambda x: x
        # x1 = pd.Series(np.arange(0, 11, 0.01))
        # y1 = fun_x1(x1).rename(x1, axis = 'rows')
        
        # ax.plot(y1, 'k--')
    
    ax.plot(ajuste, 'r--')
    ax.scatter(x, y, c = 'k')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.title(title)
    plt.legend([leg_fit, leg_data], loc = leg_loc)      
    