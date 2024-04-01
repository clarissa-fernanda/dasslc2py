import dasslc
import numpy as np
import matplotlib.pyplot as plt

# Definição da função de modelo com 15 equações diferenciais-algébricas

def model15(t, y, yp,par):
    res = np.empty(5)   

 
    #VLM=par[0]
    #n=par[1]
    #coef_1= par[0]                     #|                           |
    #coef_2 = par[1]  
    #nsol = par[2]  
    #ncr=par[3]
    #kv=par[6]
    #ncr=par[4]
    #mi4=par[8]
    G=par[0]

    Vl=y[0]               #Volume da fase líquida [m³]
    roLM=y[1]             #Densidade molar da fase líquida [mol/m³]
    mi1=y[2]              #Momento de ordem zero [#/m³]
    mi2=y[3]              #Momento de ordem um [m/m³]
    mi3=y[4]              #Momento de ordem dois [m³/m³]
   
    #nl1=y[2]                #Nº de mols de H2O na fase líquida
    #nl2=y[3]            #Nº de mols de CH4 na fase líquida
    #n=y[4]                #Nº de mols da fase líquida
    #Vs=y[5]               #Volume da fase hidrato [m³]

    
    dVldt=yp[0]        #Derivada do volume da fase líquida [m³/s]
    droLMdt=yp[1]      #Derivada da densidade molar da fase líquida [mol/ m³ s]
    #dnldt1=yp[2]     #Derivada do nº de mols de H2O na fase líquida
    #dnldt2=yp[3]     #Derivada do nº de mols de CH4 na fase líquida  
    #dndt=yp[4]         #Derivada do nº de mols da fase líquida  
    #dnldt1 = -(coef_1/coef_2)*ncr
    dmi1dt=yp[2]    #Derivada do momento de ordem zero [#/m³ s] 
    dmi2dt=yp[3]    #Derivada do momento de ordem um [m/m³ s] 
    dmi3dt=yp[4]    #Derivada do momento de ordem dois [m²/m³ s]  
    #dmi4dt=yp[5]    #Derivada do momento de ordem três [m³/m³ s]
    
 
    # Aqui você define as suas 15 equações diferenciais-algébricas
    # Estas são equações hipotéticas apenas para fins de exemplo
    res[0] = roLM - droLMdt # roLM = 1/VLM Por exemplo, yp[i] = y[i+1]
    res[1] = Vl - dVldt #Vl = n/roLM
    res[2] = mi1 - dmi1dt #MOMENTO DE ORDEM 1 = 
    res[3] = dmi2dt - G*mi1  #MOMENTO DE ORDEM dmidt(2) = G*mi(1)      [m/m³]
    res[4] = dmi3dt - 2*G*mi2  #MOM2ENO DE ORDEM 2
    #res[5] = dmi4dt - 3*G*mi3  #dmi3dt4) = 3.d0*G*mi(3) ![m³/m³]  
    #res[2] = dnldt1-  -(coef_1/coef_2)*ncr  #!dnldt(1) = -(coef_1/coef_2)*ncr
    #res[3] = dnldt2 - (nsol - ncr)                #dnldt(2) = nsol - ncr
    #res[4] = n- dndt     #dndt = nsol - ((coef_1/coef_2)+1)*ncr 
    
    # A última equação pode ser uma condição de contorno ou uma restrição algébrica

    print(y)

    return res, 0 #----------------- ires can be a literal

# Condições iniciais diferentes de zero para cada variável
y0 = np.array([1,2,1,1,2])  # Valores iniciais diferentes para cada variável

# Outros parâmetros e configurações permanecem os mesmos
t0 = np.array([500])
yp0 = None  # Não estamos fornecendo derivadas iniciais
par = np.array([1]) # Não estamos usando parâmetros adicionais neste exemplo
atol = 1e-10  # Tolerância absoluta
rtol = 1e-10  # Tolerância relativa


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
