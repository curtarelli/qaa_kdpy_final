'''
@author: Rogério Flores Jr.
@coauthor: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------

Função para conversão de valores de reflectância de sensoriamento remoto superficial
(Rrs) para subsuperficial (rrs).
'''
# Função para transformação do Rrs(acima da água) para rrs(abaixo da água)
Rrs_to_rrs = lambda Rrs: Rrs / (0.52 + 1.7 * Rrs)