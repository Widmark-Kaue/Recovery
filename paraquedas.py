# -*- coding: utf-8 -*-
"""
Created on Tue May 26 19:18:21 2020

@author: Widmark Kauê
"""
import numpy as np
import matplotlib.pylab as plt

#M = 3.437 - massa do foguete no lançamento
#Massa do propelente = 533 g 
m = 2.904 #kg - massa do foguete após a queima do propelente
Vc = 7.5 #m/s - Velocidade terminal utilizada para os cálculos
w = 10 #m/s - Velocidade terminal máxima requirida na competição


def diameter(m, Vc):
    p = 1.2225 #kg/m^3 - densidade do ar
    Cd = 1.5 #coeficiente de arrasto do paraquedas
    g = 9.81 #m/s^2 - acelaração da gravidade
    Num = 8*m*g
    Den = p*np.pi*Cd*(Vc**2)
    d = np.sqrt(Num/Den)
    return d

def grafico(m, Vc, w = 10):
    """
    Cálcula o diâmetro do paraquedas necessário para a velocide terminal de entrada e plota
    um gráfico com a relação Velocidade terminal por massa de projeto incluindo uma
    faixa vermelha, que varia + ou - 1m/s da velocidade terminal exigida pela competição, como
    range de imprevisibilidade.

    Obs:Considera somente um evento de recuperação. 
        
    Parameters
    ----------
    m : Float
        Massa estimada do projeto em kg.
    Vc : Float
        Velocidade terminal utilizada para o cálculo do diâmetro.
    w : Float, optional
        Velocidade Terminal requirida na competição. The default is 10.

    Returns
    -------
    d : Float
        Diâmetro do paraquedas em metros.

    """
    
    ##Constantes
    d = diameter(m, Vc) #resolvendo o cálculo para o diâmetro do paraquedas
    p = 1.2225 #kg/m^3 - densidade do ar
    Cd = 1.5 #coeficiente de arrasto do paraquedas
    g = 9.81 #m/s^2 - acelaração da gravidade
    A = np.pi*((d/2)**2) #Área de referência do escoamento
    W = np.arange(w - 1,w, 0.1) #Máxima velocidade terminal aceitada pelo projeto
    
    ##Função
    Den = A*Cd*p
    x = np.arange(m,2*m,0.1)
    Num = lambda x: 2*x*g
    y = lambda x: np.sqrt(Num(x)/Den)

    #Pontos a serem mostradas
    i = Den*(W**2)/(2*g)
    
    ##Plot
    plt.plot(x,y(x))
    plt.plot(i,y(i),'ro')
    plt.xlabel ('Massa (kg)')
    plt.ylabel ('Vel. Terminal (m/s)')
    
    plt.show()
    
    return d
    
    
    
    
    