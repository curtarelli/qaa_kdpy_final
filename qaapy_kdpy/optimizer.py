'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Otimizador utilizados no QAA v5 - Passo 2.
'''
##  Importando pacotes básicos necessários para o script
import numpy as np
from scipy.optimize import curve_fit

###############################################################################
def optimizer_aLb0(eq, lb0, p0, aw, a_acs, NOToptimize = False):
    '''
    Função de otimização da banda inicial de abosorção (Lambda Zero).

    ----------
    Parameters
    ----------
    eq [Data Frame]
        Data Frame contendo resultados da equação de estimativa de absorção total
        com base nos parâmetros padrão establecidos para o QAA v5.
    lb0 [Value]
        Comprimento de onda inicial (Lambda Zero).
    p0 [List]
        Lista contendo os parâmetros padrão de ajuste da equação de absorção.
    aw [Series]
        Absorção da água.
    a_acs [Data Frame]
        Absorção total do ACS .
    NOToptimize [Bool]
        Argumento de decisão pela optimização ou não da equação. Por padrão False.

    -------
    Returns
    -------
    Optimized [Dictionary]
        Dicionário contendo nova absorçãocalculada, valores dos parâmetros da 
        equação optimizados e antigos.
    '''
    aw_Lb0 = aw.loc[lb0]
    try:
        at_Lb0 = a_acs.loc[lb0]     ## + aw_Lb0.values modificado por rogerio no dia 21/01
    except:
        pass
    
    if not NOToptimize:
        # Definindo função de ajuste
        fun = lambda eq, a, b, c: aw_Lb0.values + np.power(10, (a + b  * eq + c * (np.power(eq, 2))))
        h = curve_fit(fun, eq, at_Lb0, p0, method= 'lm', maxfev = 100000)[0] 
    else:
        h = p0
    
    aref1 = aw_Lb0.values + np.power(10,( h[0] + h[1] * eq + h[2] * np.power(eq, 2)))
        
    return {'aref': aref1, 'h_otimizado': h, 'chute_inicial': p0}
