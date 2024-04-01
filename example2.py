import dasslc
import numpy as np
import matplotlib.pyplot as plt

# Definição da função de modelo com 15 equações diferenciais-algébricas

def model15(t, y, yp, par):
    res = np.empty(15)   

    G = par[0]                     #|||||||||||||||||||||||||||||
    coef_1= par[1]                     #|                           |
    coef_2 = par[2]  
    nsol = par[3]  
    kv=par[4]
    Vr=par[5]
    roHM=par[6]
    VLM=par[7]
    VgM=par[8]
    ncr=par[9]
    

    nl1=y[0]           #Nº de mols de H2O na fase líquida
    nl2=y[1]            #Nº de mols de CH4 na fase líquida
    n=y[2]                #Nº de mols da fase líquida
    ng=y[3]     #Nº de mols de CH4 na fase vapor pura
    nH=y[4]               #Nº de mols da fase hidrato
    Vl=y[5]               #Volume da fase líquida [m³]
    roLM=y[6]             #Densidade molar da fase líquida [mol/m³]
    Vs=y[7]               #Volume da fase hidrato [m³]
    Vg=y[8]               #Volume da fase vapor [m³]
    rogM=y[9]            #Densidade molar da fase vapor [mol/m³]
    mi1=y[10]              #Momento de ordem zero [#/m³]
    mi2=y[11]              #Momento de ordem um [m/m³]
    mi3=y[12]              #Momento de ordem dois [m³/m³]
    mi4=y[13]              #Momento de ordem três [m³/m³]

             
    dnldt1=yp[0]     #Derivada do nº de mols de H2O na fase líquida
    dnldt2=yp[1]     #Derivada do nº de mols de CH4 na fase líquida  
    dndt=yp[2]         #Derivada do nº de mols da fase líquida  
    dngdt=yp[3]        #Derivada do nº de mols da fase vapor
    dnHdt=yp[4]        #Derivada do nº de mols da fase hidrato
    dVldt=yp[5]        #Derivada do volume da fase líquida [m³/s]
    drolMdt=yp[6]      #Derivada da densidade molar da fase líquida [mol/ m³ s]
    dVsdt=yp[7]        #Derivada do volume da fase sólida [m³/s]
    dVgdt=yp[8]        #Derivada do volume da fase vapor [m³/s]
    drogMdt=yp[19]     #Derivada da densidade molar da fase vapor [ mol/ m³ s]
    dmidt1=yp[10]    #Derivada do momento de ordem zero [#/m³ s] 
    dmidt2=yp[11]    #Derivada do momento de ordem um [m/m³ s] 
    dmidt3=yp[12]    #Derivada do momento de ordem dois [m²/m³ s]  
    dmidt4=yp[13]    #Derivada do momento de ordem três [m³/m³ s] 

    # Aqui você define as suas 15 equações diferenciais-algébricas
    # Estas são equações hipotéticas apenas para fins de exemplo
    res[0] = roLM - (1/VLM)   # roLM = 1/VLM Por exemplo, yp[i] = y[i+1]
    res[1] = Vl - (n/roLM) #Vl = n/roLM
    res[2] = Vs - (kv*mi4*Vl) #Vs = kv*mi(4)*Vl  ![m³]            
    res[3] = nH - (roHM*Vs) #nH = roHM*Vs
    res[4] = Vg - (Vr-Vs-Vl) #Vg = Vr - Vs - Vl  ![m³]
    res[5] = rogM - (1/VgM)  #rogM = 1.d0/VgM
    res[6] = rogM - rogM*Vg #ng = rogM*Vg [mol]
    res[7] = mi1 - 0 #MOMENTO DE ORDEM 1 = 0
    res[8] = dmidt2 - G*mi1  #MOMENTO DE ORDEM dmidt(2) = G*mi(1)      [m/m³]
    res[9] = dmidt3 - 2*G*mi2  #MOMENTO DE ORDEM 2
    res[10] = dmidt4 - 3*G*mi3  #dmidt(4) = 3.d0*G*mi(3) ![m³/m³]  
    res[11] = dnldt1 - (-(coef_1/coef_2)*ncr)  #!dnldt(1) = -(coef_1/coef_2)*ncr
    res[12] = dnldt2 - (nsol - ncr)               #dnldt(2) = nsol - ncr
    res[13] = dndt - (nsol - ((coef_1/coef_2)+1)*ncr)     #dndt = nsol - ((coef_1/coef_2)+1)*ncr 
    
    # A última equação pode ser uma condição de contorno ou uma restrição algébrica
    
    return res, 0 #----------------- ires can be a literal


# Condições iniciais diferentes de zero para cada variável
y0 = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
               1.1, 1.2, 1.3, 1.4, 1.5])  # Valores iniciais diferentes para cada variável

# Outros parâmetros e configurações permanecem os mesmos
t0 = np.linspace(0, 1, 100)  # Vetor de tempo inicial
yp0 = None  # Não estamos fornecendo derivadas iniciais
par = np.array([9.81,0.1,0.2,0.2,0.5,1.5,0.01,0.05,0.1,0.1,0.9,0.1]) # Não estamos usando parâmetros adicionais neste exemplo
atol = 1e-8  # Tolerância absoluta
rtol = 1e-6  # Tolerância relativa


# Resolvendo o sistema de equações usando dasslc
t, y, yp = dasslc.solve(model15, t0, y0, yp0, par)


# Plotando os resultados
plt.figure(figsize=(10, 6))
plt.plot(t,yp)

plt.xlabel('Tempo')
plt.ylabel('Variáveis de Estado')
plt.title('Solução do Sistema de 15 Equações Diferenciais-Algébricas')
plt.legend()
plt.grid(True)
plt.show() 
