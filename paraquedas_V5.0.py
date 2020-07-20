# -*- coding: utf-8 -*-
"""
Created on Tue May 26 19:18:21 2020

@author: Widmark Kauê - Gerente de Recovery, Kosmos Rocketry
"""
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pylab as plt

#M = 3.437 - massa do foguete no lançamento
#Massa do propelente = 533 g 
m = 2.904 #kg - massa do foguete após a queima do propelente
Vc = 7.5 #m/s - Velocidade terminal utilizada para os cálculos
w = 10 #m/s - Velocidade terminal máxima requirida na competição


def diameter(m, Vc, Cd = 1.5):
    """
    Cálcula o diâmetro nominal do paraquedas considerando o sistema paraquedas-carga de
    recuperação em regime estacionário.
    
    Parameters
    ----------
    m : Float
        Massa do foguete após a queima do propelente, kg.
    Vc : Float
        Velocidade terminal desejada, m/s.
    Cd : Float, optional
        Coeficiente de arrasto do paraquedas, adimensional. The default is 1.5.

    Returns
    -------
    d : Float
        Diâmetro nominal do paraquedas, metros.
    p : Float
        Densidade do ar, kg/m^3.
    Cd : Float
        Coeficente de Arrasto padrão do paraquedas, adimensional.
    g : Float
        Aceleração da gravidade, m/s^2.

    """
    p = 1.2225 #kg/m^3 - densidade do ar
    g = 9.81 #m/s^2 - aceleração da gravidade
    Num = 8*m*g
    Den = p*np.pi*Cd*(Vc**2)
    d = np.sqrt(Num/Den)
    
    #retorna o diâmetro do paraquedas e os parâmetro principais para os cálculos.
    return d, p, Cd, g

def componentes(d, caps):
    """
    Executa o dimensionamento dos componentes listados abaixo:
        Cordas de suspensão;
        Caps.
    
    Obs: Para o diâmetro do Caps está sendo considerado o valor de 20% do diâmetro nominal.
    
    Parameters
    ----------
    d : Float
        Diãmetro nominal do paraquedas.
    caps : Boolean
        Condição de existência de caps no projeto.


    Returns
    -------
    Cs : Float
        Comprimento das cordas de suspensão, metros.
    c : Float
        Diâmetro do caps do paraquedas, metros.
    
    """
    Cs = 1.2*d #Cálculo do comprimento dos cabos de suspensão
    if (caps == True):
        c = 0.2*d  #Cálculo do diâmetro do caps do paraquedas
    else:
        c = 0
    return Cs, c

def area(d, c, caps):
    """
    Calcula a área nominal (quantidade de tecido) de um paraquedas semi-elipsóide com 
    proporção de b/a = 0.707 utilizando o resultado de uma integral de superfície.
    
    Cálcula, também, a área de nominal (quantidade de tecido) necessária para o caps do paraquedas.
    Obs: Essa fórmula só é valida para paraquedas sem buracos.

    Parameters
    ----------
    d : Float
        Diâmetro nominal do paraquedas, metros.
    c : Float
        Diâmetro do caps.
   caps : Boolean
        Condição de existência de caps no projeto.

    Returns
    -------
    S : Float
        Área nominal do paraquedas, m^2.
    Sc: Float
        Área nominal do caps, m^2

    """
    ##Cálculo da área nominal do paraquedas
    a = 0.5*d
    b = 0.707*a
    e = (np.sqrt((a**2) - (b**2)))/a
    l = np.log((1 + e)/(1 - e))
    S = 0.5*(2*np.pi*a**2 + (np.pi*(b**2)/e)*l)
    
    #Cálculo da área do caps, se exixtir no projeto.
    if (caps == True):
        Sc = np.pi*(c**2)/4
    else:
        Sc = 0
        
    return S, Sc

    
    
def Fmax(m, Vc, S, D, Ap, v1 = 5, n0 = 14):
    """
    Calcula as esforços exigidos no paraquedas em dois momentos: na abertura dele e em regime 
    estacionário para o tipo de paraquedas ringslot.

    Parameters
    ----------
    m : Float
        Massa do foguete após a queima do propelente, kg.
    Vc : Float
        Velocidade terminal, m/s.
    S : Float
        Área de refeência do escoamento, m^2.
    D : Array de Floats
       Arrays com valores de diâmetro nominal, densidade do ar, coeficiente de arrasto e acelaração da gravidade.
    Ap : Float
        Apogeu estimado do foguete, metros.
    v1 : Float, optional
        Velocidade de alongamento das linhas, m/s.
    n0 : Int, optional
        Constante de cada tipo de paraquedas. The default is 14.

    Returns
    -------
    Fo : Float
        Força de abertura no apogeu, Newtons.
    Fd : Float
        Força de arrasto em regime estacionário, Newtons.
    Fd_h : Function
        Função Força de arrasto por altura de abertura do paraquedas.

    """
############################Força de arrasto####################################
    
    Fd = 0.5*D[2]*D[1]*(Vc**2)*S #Cálculo da força de arrasto em regime estacionário
    
    #####Função força de arrasto por altura#######
    V_h = lambda h: np.sqrt(2*D[3]*(h)) #Função velocidade em queda livre por distância do apogeu.
    Fd_h = lambda h :0.5*D[2]*D[1]*((V_h(h))**2)*S #Função força de arrasto por altura.
    
    #####Função força de arrasto por tempo#######
    V_t = lambda t: -D[3]*t #Função velocidade em queda livre por tempo.
    Fd_t = lambda t: 0.5*D[2]*D[1]*((V_t(t))**2)*S #Função força de arrasto por altura.
    
    
###################Força de abertura do paraquedas##############################
    
    #Obs: Utilza o Método de Pflanz para o cálculo do fator de choque.
    

    #####Constantes#######
    m = m*2.2046 #massa do sistema em lb
    S = S*10.7639 #área de referência do escoamento em ft^2
    v1 = v1*3.2808 #velocidade de alongamento das linhas ft/s
    d = D[0]*3.2808 #diâmetro nominal do paraquedas em ft
    p = D[1]*0.00194032 #densidade do ar em slug/ft^3
    g = D[3]*3.2808 #gravidade em ft/s^2
    
    #######Cálculos dos parâmetros#######
    q = 0.5*p*(v1**2) #cálculo da pressão dinâmica
    tf = (n0*d)/v1 #cálculo do tempo de inflação do paraquedas
    Num = 2*m
    Den = S*p*g*v1*tf
    A = Num/Den #cálculo do parâmetro balístico 
    
    #O fator de redução da força é encontrado utilizando um gráfico disponibilizado no livro do Knacke.
    #Dado nosso tipo de paraquedas, ringslot, vemos que ele se encaixa dentro da curve 1 para n = 1 representado no livro.
    Mat = np.loadtxt("OFRFxBP.txt") #Matrix com valores de x e y para n = 1 do gráfico.
    X1 = interpolate.interp1d(Mat[:,0],  Mat[:,1], kind = 'cubic') #Interpolação cúbica dos valores da matrix.
    x1 = X1(A) #Identificação do fator de choque do paraquedas.
    
    Fo = S*q*x1 #Cálculo da força de abertura do paraquedas em lbf.
    Fo = Fo*4.44822 #Força de abertura do paraquedas em Newtons.
    
    # #####Função força de abertura por tempo#######
    # Num1_t = lambda t: np.sqrt(2*t*g) 
    # #Obs: vai ser necessário utilizar um método numérico para determinar tf.
    # Den1_t = lambda t: (S*p*g*tf*Num1_t(t))/(Num)
    # v1_t = lambda t: (Num1_t(t))/(1 + Den1_t(t))

    # ##Parâmetros em função do tempo###
    # q_t = lambda t: 0.5*p*(v1_t(t)**2) #cálculo da pressão dinâmica em função do tempo.
    # Den_t = lambda t: S*p*g*v1_t(t)*tf 
    # A_t = lambda t: Num/Den_t(t) #cálculo do parâmetro balístico em função do tempo. 
    
    # x1_t = lambda t: X1(A_t(t)) #Identificação do fator de choque do paraquedas em função do tempo.    
    
    # Fo_t = lambda t : S*q_t(t)*x1_t(t) #Função força de abertura do paraquedas em lbf em função do tempo.
    # Fo1_t = lambda t:  Fo_t(t)*4.44822 #Força de abertura do paraquedas em Newtons.
    
    return Fo, Fd, Fd_h, Fd_t
    
    
def main(m, Vc, w = 10, Ap = 1000, Caps = True, Ret = False, Prt = True):
    """
    Cálcula o diâmetro do paraquedas necessário para a velocide terminal de entrada e plota
    um gráfico com a relação Velocidade terminal por massa de projeto incluindo uma
    faixa vermelha, que varia + ou - 1m/s da velocidade terminal exigida pela competição, como range de imprevisibilidade.

    Obs:Considera somente um evento de recuperação. 
        
    Parameters
    ----------
    m : Float
        Massa estimada do projeto, kg.
    Vc : Float
         Velocidade terminal desejada, m/s.
    w : Float, optional
        Velocidade Terminal requirida na competição. The default is 10.
    Ap : Float, optional
        Apogeu esperado do foguete, metros. The default is 1000.
    Caps : Boolean, optional
        Condição de falso ou verdadeiro para a existência do caps no paraquedas. The default is True.
    Ret : Boolean, optional
        Condição de falso ou verdadeiro para o retorno dos parâmetros. The default is False.
    Prt : Boolean, optional
        Condição de falso ou verdadeiro para o print dos valores. The default is True.
    
    Returns
    -------
    D[0] : Float
          Diâmetro nominal do paraquedas, metros.
    C[0] : Float
        Comprimento dos cabos de suspensão, metros.
    C[1] : Float
        Diâmetro do caps do paraquedas, metros.
    A[0] : Float
        Área nominal do velame, m^2.
    A[1] : Float
        Área nominal do caps, m^2.
    Fv[0] : Float
        Força de abertura do paraquedas utilizando o método de Pflanz retirado do livros o Knacke, Newtons.
    Fv[0] : Float
        Força de arrasto do paraquedas em regime estacionário, Newtons.
    Fd_t : Function
        Força de arrasto em função do tempo considerando o foguete em queda livre.
    Fd_h: Function
        Força de arrasto em função da altura considerando o foguete em queda livre.
    V_m: Function
        Velocidade terminal do foguete em função da massa considerando o diâmetro dimensionado.
    V_dm: Function
        Velocidade terminal do foguete em função da massa e do diâmetro.    
    """   
    
######################################Constantes############################################
    D = diameter(m, Vc) #resolvendo o cálculo para o diâmetro do paraquedas 
    #D: vetor com parâmetros de diâmetro, pressão, coeficiente de arrasto e aceleração da         gravidade nessa ordem.
    
    A0 = np.pi*((D[0]/2)**2) #Área de referência do escoamento para o paraquedas dimensionado.
    W = np.arange(w - 1,w, 0.1) #Máxima velocidade terminal aceitada pelo projeto
        
    
###############################Esforços no paraquedas###############################    
    Fv = Fmax(m, Vc, A0, D, Ap) #Vetor com os esforços sofridos pelo paraquedas
    if(Prt == True):
        print("Esforços no paraquedas")
        print("Força de abertura do paraquedas:",Fv[0], "Newtons")
        print("Força de arrasto em regime estacionário:",Fv[1], "Newtons")
        
    ###
    ##Função força de arrasto por distância do apogeu.
    Fd_h = Fv[2] #Recebe a função do vetor de forças
    H = np.arange(0,11)
    #Plot
    plt.title("Força de arrasto por Distância do apogeu")
    plt.xlabel ('Distância (m)')
    plt.ylabel ('Força de arrasto (N)')
    plt.grid(True)
    plt.plot(H,Fd_h(H),'r')
    
    plt.show()
    ###
    
    ##Função força de arrasto por tempo.
    Fd_t = Fv[3] #Recebe a função do vetor de forças
    T = np.arange(0, 6, 1) #Cria um vetor com valores de tempo de 0 até 6s.
    #Plot
    plt.title("Força de arrasto por tempo")
    plt.xlabel ('Tempo (s)')
    plt.ylabel ('Força de arrasto (N)')
    plt.grid(True)
    plt.plot(T,Fd_t(T),'g')
    
    plt.show()
    
    # ##Função força de abertura do paraquedas por tempo.
    # Fo_t = Fv[4]
    # #Plot
    # plt.title("Opening-Force x Time")
    # plt.xlabel ('Tempo (s)')
    # plt.ylabel ('Opening-Force (N)')
    # plt.grid(True)
    # plt.plot(T,Fo_t(T))
    
    # plt.show()
###############################Dimensionamento de componentes###############################
    C = componentes(D[0], Caps)
    if(Prt == True):
        print("\nDimensionamento dos componentes")
        print("Diâmetro nominal do velame:",D[0],"metros")
        print("Comprimento dos cabos de suspensão:",C[0], "metros")
        if (Caps == True):
            print("Diâmetro do caps:",C[1], "metros")
        
##########################Dimensionamento da quantidade de tecido###############################        
    A = area(D[0], C[1], caps = Caps)
    if (Prt == True):
        print("\nQuantidade de tecido necessária")
        print ("Área nominal do velame:", A[0], "m^2")
        if (Caps == True):
            print("Área nominal do caps:", A[1],"m^2")
        print("Área total:", A[0] + A[1],"m^2")
    
################################Função Vel. Terminal por massa##############################
    Den = A0*D[2]*D[1]
    x = np.arange(m,2*m,0.1) #range de massas com parada em duas vezes a massa do foguete.
    Num =  lambda x: 2*x*D[3]
    V_m = lambda x: np.sqrt(Num(x)/Den) #Função da vel. terminal por massa
    i = Den*(W**2)/(2*D[3]) #Pontos a serem mostradas do range de imprevisibilidade
    
    #Plot
    plt.title("Range de Massa ideal")
    plt.xlabel ('Massa (kg)')
    plt.ylabel ('Vel. Terminal (m/s)')
    plt.grid(True)
    plt.plot(x,V_m(x))
    plt.plot(i,V_m(i),'ro', label = "Range de imprevisibilidade.")
    plt.legend()
    
    plt.show()

#########################Função Vel. Terminal por massa e diâmetro##############################
    d = np.linspace(D[0] - 0.3, D[0] + 0.3, 7) #range de diâmetros entre + ou - 0.1 do diâmetro nominal. 
    d = d*100
    A0_d = lambda d: np.pi*((d/2)**2) #área de referência do escoamento em função do diâmetro 
    Den_d = lambda d: A0_d(d)*D[2]*D[1] 
    V_dm = lambda x, d: np.sqrt(Num(x)/Den_d(d)) #Função da vel. terminal por massa
    
    #Plot
    plt.title("Curvas de diâmetros")
    plt.xlabel ('Massa (kg)')
    plt.ylabel ('Vel. Terminal (m/s)')
    plt.grid(True)
    for j in range(len(d)):
        d[j] = int(d[j])/100
        plt.plot(x,V_dm(x, d[j]), label = str(d[j])) #Plota as funções para cada um dos diâmetros
    plt.legend()
    
    plt.show()
    
    if (Ret == True):
        return D[0], C[0], C[1], A[0], A[1], Fv[0], Fv[1], Fd_t, Fd_h, V_m, V_dm
    
    
    
    
    
    
    