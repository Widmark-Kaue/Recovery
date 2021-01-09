# -*- coding: utf-8 -*-
"""
Created on Tue May 26 19:18:21 2020

@author: Widmark Kauê - Gerente de Recovery, Kosmos Rocketry
"""
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import matplotlib.pylab as plt

#M = 3.437 - massa do foguete no lançamento
#Massa do propelente = 533 g 
m = 2.904 #kg - massa do foguete após a queima do propelente
Vc = 7.5 #m/s - Velocidade terminal utilizada para os cálculos
w = 10 #m/s - Velocidade terminal máxima requirida na competição
Cd_f = (0.46 + 0.5)/2 #Coeficiente de arrasto do foguete. Levei em consideração um valor médio extraido do Open Rocket
d_f = 0.078 #metros -  Diâmetro do corpo do foguete.
S_f = np.pi*(d_f**2)/4 #m^2 - área de referência de escoamento do foguete.

def diameter(m, Vc, Cd = 1.5):
    """
    Cálcula o diâmetro nominal do paraquedas considerando o sistema paraquedas-carga de
    recuperação em regime permanente.
    
    Obs.: Considera-se que a variação do coeficiente de arrasto não é tão alta com a
    área e é adotado um valor constante de 1.5 como é mostrado nas literaturas.
    
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
    
    Fd = 0.5*D[2]*D[1]*(Vc**2)*S #Cálculo da força de arrasto em regime permanen
    
    #####Função pico força de arrasto por altura#######
    V_h = lambda h: np.sqrt(2*D[3]*(h)) #Função velocidade em queda livre por distância do apogeu.
    Fd_h = lambda h :0.5*D[2]*D[1]*((V_h(h))**2)*S #Função pico força de arrasto por altura.
    
    #####Função pico força de arrasto por tempo#######
    Vp_t = lambda t: -D[3]*t #Função velocidade em queda livre por tempo.
    Fpd_t = lambda t: 0.5*D[2]*D[1]*((Vp_t(t))**2)*S #Função pico de força de arrasto por altura.
    
    
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
    Den = D[2]*S*p*g*v1*tf
    A = Num/Den #cálculo do parâmetro balístico 
    
    #O fator de redução da força é encontrado utilizando um gráfico disponibilizado no livro do Knacke.
    #Dado nosso tipo de paraquedas, ringslot, vemos que ele se encaixa dentro da curve 1 para n = 1 representado no livro.
    Mat = np.loadtxt("OFRFxBP.txt") #Matrix com valores de x e y para n = 1 do gráfico.
    X1 = interp1d(Mat[:,0],  Mat[:,1], kind = 'cubic') #Interpolação cúbica dos valores da matrix.
    x1 = X1(A) #Identificação do fator de choque do paraquedas.
    
    Fo = D[2]*S*q*x1 #Cálculo da força de abertura do paraquedas em lbf.
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
    
    return Fo, Fd, Fd_h, Fpd_t, V_h
    
def Voo(V_h, mt, D, A0, Cd_f, S_f, R, Ap):
    """
    Resolve os sistema de equações diferenciais da dinâmica de voo do paraquedas de forma que
    é possível modelar funções da velocidade, força de arrasto e altura pelo tempo em
    regime transiente.
    
    Parameters
    ----------
    V_h : Function
        Função da velocidade em queda livre pela altura.
    mt : Float
        Masso total do sistema.
    D : Array
        D[0] : Float
            Diâmetro nominal do paraquedas, metros.
        D[1] : Float
            Densidade do ar, Kg/m^3.
        D[2]: Float
            Coeficiente de arrasto do paraquedas, adimensional.
        D[3]: Float
            Aceleração da gravidade, m/s^2.
    A0 : Float
        Área de referência de escoamento do paraquedas, m^2.
    Cd_f : TYPE
        Coeficiente de arrasto do foguete, adimensional.
    S_f : Float
        Área de referência de escoamento do foguete, m^2.
    R : Float
        Altura de implementação do paraquedas (distância do apogeu), metros.
    Ap : Float
        Apogeu estimado do foguete, metros.

    Returns
    -------
    H_t : Function
        Função da altura pelo tempo.
    V_t : Function
        Função velocidade terminal pelo tempo.
    Fd_t : Function
        Função força de arrasto pelo tempo.

    """
    ##Tempo para a altura de implementação do paraquedas##
    h_t = lambda t:  Ap - 0.5*D[3]*(t**2) #função da altura pelo tempo em queda livre.
    B = lambda t: h_t(t) - R #função para aplicar na bissecção
    
    a1 = b1 =  1 #Definindo  primeiro e segundo ponto do intervalo da bisseção.
   
        ##Encontrando onde a função muda de sinal 
    while (B(b1)>0):
        b1 = b1 + 10
    
    #encontrando tempo para atingir altura de implementação do foguete.
    t1 = bisect(B, a1, b1) 
    
    ##Velocidade de descida em função do tempo##
    k = (-D[1]*((Cd_f*S_f) + (D[2]*A0)))/(2*mt) #Cálcula uma constante padrão durante todo processo de integração
    V0 = V_h(R) #Determina a velocidade inicial do sistema de acordo com ruido do sensor.
    
    dvdt = lambda V, t: (k*(V)**2 + D[3]) #monta o modelo utilizado para a integração.
    t0 = np.linspace(0,300,900) #determina um range de tempo de até 5 min.
    V = integrate.odeint(dvdt, V0, t0) #Realiza a integração e retorna valores de velocidade para o range de tempo determinado. 
    
    #Nota: a função odeint retorna um array de tamanho (300,1), por conta de como a função de interpolação funciona 
    #é necessário passar os valores como está abaixo.
    V_t = interp1d(t0, V[:,0], kind = 'cubic') #Utilizar os valores de V e t para interpolar (cúbica) uma função.
    
    
    ##Força de arrasto no paraquedas em função do tempo##
    Fd_t = lambda t: 0.5*D[2]*D[1]*(V_t(t)**2)*A0 
    
    ##Altura em que se encontra o sistema em função do tempo.##
    Hi = Ap - R #detrmina a altura de implementação do paraquedas em relação ao solo.
    
    #Nota:Como a função quad retorna  uma tupla, o valor da integral e o erro, não é possível realizar operações aritiméticas 
    #com ela. Dessa forma, é cáculado primeiramente a tupla e depois acessado somente o valor da integral. 
    Quad = lambda t: integrate.quad(V_t,0,t) #Cálculo da integral e retorno da tupla.
    H_t = lambda t: Hi - Quad(t)[0] #Monta a função da altura pelo tempo.
    
    return  h_t, H_t, V_t, Fd_t, t1
    
def main(m, Vc, Cd_f, S_f, w = 10, Ap = 1000, R = 500, Caps = True, Ret = False, Prt = True, Plt = True):
    """
    Cálcula o diâmetro nominal do paraquedas necessário para a massa e a velocidade terminal de
    entrada.
    
    Cálcula os esforços atuantes no paraquedas, sendo eles a força de arrasto em regime permanente
    e a força de abertura utilizando o método de Pflanz. Também plota gráficos para o pico de força
    de arrasto que o paraquedas teria que amortecer caso a implementação não ocorrese no apogeu, estes
    valores variando com o tempo e a distênci do apogeu considerando o sistema em queda livre.
    
    Resolve o sistema de equações diferênciais da dinâmica de voo do paraquedas e retorna funções
    de altura, velocidade e força de arrasto pelo tempo incluindo os regimes transiente e permanente.
    Além disso, encontra o tempo necessário para o sistema entrar em regime permanente e o tempo
    necessário para o sistema atingir o solo.
    
    Dimensiona os componentes do sitema do paraquedas, cordas de suspensão e caps, com base no
    diâmetro cálculado anteriormente.
    
    Dimensiona a quantidade de tecido necessária para a manufatura do paraquedas, separando em
    tecido para o velame e para o caps.
    
    Plota um gráfico com a relação Velocidade terminal por massa de projeto, para o diâmetro nominal
    estimado, incluindo uma faixa vermelha, que varia + ou - 1m/s da velocidade terminal exigida
    pela competição, como range de imprevisibilidade.
    
    Plota um gráfico com a relação Velocidade terminal por massa para várias diâmetros de
    paraquedas, variando 30 cm para mais ou para menos do diâmatro nominal estimado no
    código.

    Obs:Considera somente um evento de recuperação. 
        
    Parameters
    ----------
    m : Float
        Massa estimada do projeto, kg.
    Vc : Float
        Velocidade terminal desejada, m/s.
    Cd_f: Float
        Coeficiente de arrasto do foguete, adimensional.
    S_f: Float
        Área de referência de escoamento do foguete, m^2.
    w : Float, optional
        Velocidade Terminal requirida na competição. The default is 10.
    Ap : Float, optional
        Apogeu esperado do foguete, metros. The default is 1000.
    R: Float, optional.
        Altura de implementação do paraquedas (distância do apogeu), metros.
    Caps : Boolean, optional
        Condição de falso ou verdadeiro para a existência do caps no paraquedas. The default is True.
    Ret : Boolean, optional
        Condição de falso ou verdadeiro para o retorno dos parâmetros. The default is False.
    Prt : Boolean, optional
        Condição de falso ou verdadeiro para o print dos valores. The default is True.
    Plt : Boolean, optional
        Condição de falso ou verdadeiro para o plot dos gráficos. The default is True.
        
    Returns
    -------
    D[0] : Float
          Diâmetro nominal do paraquedas, metros.
    C : Array
        C[0] : Float
            Comprimento dos cabos de suspensão, metros.
        C[1] : Float
            Diâmetro do caps do paraquedas, metros.
    A : Array 
        A[0] : Float
            Área nominal do velame, m^2.
        A[1] : Float
            Área nominal do caps, m^2.
    Fv[0:2] : Array
            Fv[0] : Float
                Força de abertura do paraquedas utilizando o método de Pflanz retirado do livros o Knacke, Newtons.
            Fv[1] : Float
                Força de arrasto do paraquedas em regime permanente, Newtons.
    Fd_t : Function
        Força de arrasto em função do tempo considerando o foguete em queda livre.
    Fd_h: Function
        Força de arrasto em função da altura considerando o foguete em queda livre.
    V_m: Function
        Velocidade terminal do foguete em função da massa considerando o diâmetro dimensionado.
    V_dm: Function
        Velocidade terminal do foguete em função da massa e do diâmetro.
    h_t: Function
        Altura do sistema em queda livre em função do tempo.
    H_t: Function
        Altura do sistema em função do tempo.
    V_t: Function
        Velocidade do sistema em função do tempo.
    Fd_t: Function
        Força de arrasto do paraquedas em função do tempo.
    """   
    
######################################Constantes############################################
    D = diameter(m, Vc) #resolvendo o cálculo para o diâmetro do paraquedas 
    #D: vetor com parâmetros de diâmetro, densidade, coeficiente de arrasto e aceleração da         gravidade nessa ordem.
    
    A0 = np.pi*((D[0]/2)**2) #Área de referência do escoamento para o paraquedas dimensionado.
    W = np.arange(w - 1,w, 0.1) #Máxima velocidade terminal aceitada pelo projeto
        
    
###############################Esforços no paraquedas###################################    
    Fv = Fmax(m, Vc, A0, D, Ap) #Vetor com os esforços sofridos pelo paraquedas
    if(Prt == True):
        print("Esforços no paraquedas")
        print("Força de abertura do paraquedas:",Fv[0], "Newtons")
        print("Força de arrasto em regime permanente:",Fv[1], "Newtons")
        
    ###
    ##Função pico de força de arrasto por distância do apogeu.
    Fd_h = Fv[2] #Recebe a função do vetor de forças
    H = np.arange(0,11,0.5) #Cria um vetor com valores de distânica do apogeu de 0 até 11 metros.
    #Plot
    if (Plt == True):
        plt.title("Pico de Força de arrasto por Distância do apogeu")
        plt.xlabel ('Distância (m)')
        plt.ylabel ('Força de arrasto (N)')
        plt.grid(True)
        plt.plot(H,Fd_h(H),'r')
        
        plt.show()
        
    ###
    
    ##Função pico de força de arrasto por tempo.
    Fpd_t = Fv[3] #Recebe a função do vetor de forças
    T = np.arange(0, 6, 0.2) #Cria um vetor com valores de tempo de 0 até 6s.
    #Plot
    if (Plt == True):
        plt.title("Pico de Força de arrasto por tempo")
        plt.xlabel ('Tempo (s)')
        plt.ylabel ('Força de arrasto (N)')
        plt.grid(True)
        plt.plot(T,Fpd_t(T),'g')
        
        plt.show()
    
###############################Dinâmica de Voo do paraquedas###############################
    V_h = Fv[4] #Velocidade do sistema em queda livre em função da altura
    wid = Voo(V_h, m, D, A0, Cd_f, S_f, R, Ap) #vetor com funções da dinâmica de voo
    h_t = wid[0] #função da altura pelo tempo antes da implementação do paraquedas.
    H_t = wid[1] #função da altura pelo tempo com o paraquedas aberto.
    V_t = wid[2] #função da velocidade pelo tempo
    Fd_t = wid[3] #função da força de arrasto pelo tempo   
    tr = wid[4] #tempo para atingir a altura de implementação do paraquedas.
    ###
    
    #####Tempo para o sistema atingir o regime permanente e atingir o solo#####
    #Para achar o tempo em que o sistema entra em equilibrio  e atingi o solo irei utilizar o método da bissecção
    
    #Transforma o problema de achar os valores de tempo supracitados em um problema de raiz.
    T1 = lambda t: V_t(t) - Vc #tempo para o regime permanente
    T2 = lambda t: H_t(t) - 0 #tempo para atingir o solo
    
    a1 = a2 = 1 #Defininindo primeiro ponto do intervalo da bissecção.
    b1 = b2 = 1 #Definindo segundo ponto do intervalo da bisseção.
   
        ##Encontrando onde as funções mudam de sinal 
    while ((T1(b1) > 0) or (T2(b2) > 0)):
        
        if (T1(b1) > 0):    
            b1 = b1 + 10
        
        if (T2(b2) > 0):
            b2 = b2 + 60
           
    t1 = bisect(T1, a1, b1) #Aplicação da bissecção - tempo para o regime permanente
    t = np.linspace(0, t1 + 0.5*t1, 100) #Array de valores de tempos para plotagem de gráficos.
    
    t2 = bisect(T2, a2, b2) #Aplicação da bissecção - tempo para atingir o solo.
    t0 = np.linspace(0, t2+tr,100) #Array de valores de tempos para plotagem de gráficos.
    
    #Prints
    if(Prt ==  True):
        print("\nDinâmica de Voo")
        print("Tempo que o sistema leva para atingir o regime permanente:",t1,"s")
        print("Tempo que o sistema leva para atingir o solo:",t2+tr,"s")
     
    ###
    t1 = int(t1*100)/100 #truncando o valor do tempo para duas casas decimais.
    
    ##Função força de arrasto pelo tempo.
    #Plot
    if (Plt == True):
        plt.title("Força de arrasto por tempo")
        plt.xlabel ('Tempo (s)')
        plt.ylabel ('Força de arrasto (N)')
        plt.grid(True)
        plt.plot(t,Fd_t(t),'r')
        plt.plot(t1,Fd_t(t1),'o', label = str(t1))
        plt.legend(loc = 7)
        
        plt.show()
    
    ###
    
    #Função velocidade de descida pelo tempo.
    #Plot
    if (Plt == True):
        plt.title("Velocidade de descida por tempo ")
        plt.xlabel ('Tempo (s)')
        plt.ylabel ('Vel. de descida (m/s)')
        plt.grid(True)
        plt.plot(t,V_t(t),'g')
        plt.plot(t1,V_t(t1),'o', label = str(t1))
        plt.legend(loc = 7)
        
        plt.show()
    
    ###
    
    #Função altura do sistema pelo tempo.
    #Nota: Por conta de como a função quad funciona não é possível passar um vetor de uma vez para função.
    #Então,será passado valor por valor do array de tempo para a função e gerado um vetor de imagens para a plotagem.
    s = np.zeros(len(t0)) #array que vai receber as imagens.
    
    for i in range(len(t0)):
        if (t0[i]<tr):
            s[i] = h_t(t0[i]) #cálculo das imagens com o paraquedas fechado.
        else:
            s[i] = H_t(t0[i]- tr) #cálculo das imagens com o paraquedas aberto.
    
    #Plot
    if (Plt == True):
        plt.title("Altura por tempo")
        plt.xlabel ('Tempo (s)')
        plt.ylabel ('Altura (m)')
        plt.grid(True)
        plt.plot(t0,s,'b')
        plt.plot(t2+tr,H_t(t2),'ro', label = str(int(t2+tr)))
        plt.legend(loc = 7)
        
        plt.show()

###############################Dimensionamento de componentes###############################
    C = componentes(D[0], Caps) #vetor com o dimensionamento dos componentes.
    if(Prt == True):
        print("\nDimensionamento dos componentes")
        print("Diâmetro nominal do velame:",D[0],"metros")
        print("Comprimento dos cabos de suspensão:",C[0], "metros")
        if (Caps == True):
            print("Diâmetro do caps:",C[1], "metros")
        
##########################Dimensionamento da quantidade de tecido###############################        
    A = area(D[0], C[1], caps = Caps) #vetor com os valores de área nominal do velame e caps.
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
    if (Plt == True):
        plt.title("Range de Massa ideal")
        plt.xlabel ('Massa (kg)')
        plt.ylabel ('Vel. Terminal (m/s)')
        plt.grid(True)
        plt.plot(x,V_m(x))
        plt.plot(i,V_m(i),'ro', label = "Range de imprevisibilidade.")
        plt.legend()
        
        plt.show()

#########################Função Vel. Terminal por massa e diâmetro##############################
    d = np.linspace(D[0] - 0.3, D[0] + 0.3, 7) #range de diâmetros entre + ou - 0.3 m do diâmetro nominal. 
    d = d*100
    A0_d = lambda d: np.pi*((d/2)**2) #área de referência do escoamento em função do diâmetro 
    Den_d = lambda d: A0_d(d)*D[2]*D[1] 
    V_dm = lambda x, d: np.sqrt(Num(x)/Den_d(d)) #Função da vel. terminal por massa e diâmetro
    
    #Plot
    if (Plt == True):
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
        return D[0], C, A, Fv[0:2], Fpd_t, Fd_h, V_m, V_dm, h_t, H_t, V_t, Fd_t
    
    
    
    
    
    