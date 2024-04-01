from PR import PR
import dasslc
import numpy as np
import matplotlib.pyplot as plt
import math

from liquid_density_rackett import liquid_density_rackett

# Definição da função de modelo com 15 equações diferenciais-algébricas

def model15(t, y, yp,par):
    res = np.empty(14)


    #VLM=par[0]
    #n=par[1]
    #coef_1= par[0]                     #|                           |
    #coef_2 = par[1]
    #nsol = par[2]
    #ncr=par[3]
    #kv=par[6]
    #ncr=par[4]
    #mi4=par[8] 
    Vr = 0.001 
    G, coef_1, coef_2, ncr, roHM, kv, VgM, Vr, Aint, kd, neq,VLM = par

    # constantes
    #ks = 4*pi
      #volume do reator
    #Rin = 0.054 # raio do reator
    #Parâmetros admensionais do sistema
    #Re = 7000        #Reynolds
    #We = 1       #Weber
    #Eu = 10 #Euler


    roLM, Vl, mi1, mi2, mi3, mi4, nl1, nl2, n, Vs, rogM, ng, nH, Vg= y


    (droLMdt, dVldt, dmi1dt, dmi2dt, dmi3dt, dmi4dt,
     dnldt1, dnldt2, dndt, dVsdt, drogMdt, dngdt, dnHdt, dVgdt) = yp


    #taxas iniciais

    #Hl = (Vl+Vs)/(pi*Rin*Rin)
    #H = Vr/(pi*Rin*Rin)

    ncr = coef_2*roHM*kv*(Vl*dmi4dt + mi4*dVldt)
    #Aint = ((1/H)+((1/Hl)*((Re**1.75)*(We**3)/(Eu**3))))   # TAXA MOLAR DE CRESCIMENTO DA FASE HIDRATO
    nsol = (Aint*kd/Vl)*(neq-nl2)    # taxa de solubilização de metano
    ne = rogM + Vg*drogMdt + nsol   # [mol de CH4/s] #TAXA MOLAR DE ENTRADA DE METANO NO REATOR

    n += nl1
    n += nl2
    xl = [nl / n for nl in [nl1, nl2]]

    Tp = 276.0  # Substitua pelo valor desejado
    P = 7.09 # Substitua pelo valor desejado
    Acentric = [0.344, 0.011]  # Substitua pelos valores desejados
    Tc = [647.14, 190.56]  # Substitua pelos valores desejados
    Pc = [22064000 , 4599000]  # Substitua pelos valores desejados
    Vc = [55.95, 98.60]  # Substitua pelos valores desejados
    

    # Aqui você define as suas 15 equações diferenciais-algébricas
    # Estas são equações hipotéticas apenas para fins de exemplo
    res[0] = roLM - 1/VLM
    # res[0] = roLM - 1/1
    res[1] = Vl - n/roLM #Vl = n/roLM
    res[2] = mi1 - dmi1dt #MOMENTO DE ORDEM 1 =
    res[3] = dmi2dt - G*mi1  #MOMENTO DE ORDEM dmidt(2) = G*mi(1)      [m/m³]
    res[4] = dmi3dt - 2*G*mi2  #MOM2ENO DE ORDEM 2
    res[5] = dmi4dt - 3*G*mi3  #dmi3dt4) = 3.d0*G*mi(3) ![m³/m³]
    res[6] = dnldt1 -(coef_1/coef_2)*ncr  #!dnldt(1) = -(coef_1/coef_2)*ncr
    res[7] = dnldt2 - (nsol - ncr)                #dnldt(2) = nsol - ncr
    res[8] = dndt - (nsol - ((coef_1/coef_2)+1)*ncr)    #dndt = nsol - ((coef_1/coef_2)+1)*ncr
    res[9] = Vs - (kv*mi4*Vl) #Vs = kv*mi(4)*Vl  ![m³]
    res[10] = rogM - (1/VgM)  #rogM = 1.d0/VgM
    res[11] = ng - rogM*Vg #ng = rogM*Vg [mol]
    res[12] = nH - (roHM*Vs) #nH = roHM*Vs
    res[13] = Vg - (Vr-Vs-Vl) #Vg = Vr - Vs - Vl  ![m³]


    # A última equação pode ser uma condição de contorno ou uma restrição algébrica

    return res, 0 #----------------- ires can be a literal

nl1 = 16.6
nl2 = 1e-5    
n=0
n += nl1
n += nl2
xl = [nl / n for nl in [nl1, nl2]]
Tp = 276.0  # Substitua pelo valor desejado
P = 7.09 # Substitua pelo valor desejado
Acentric = [0.344, 0.011]  # Substitua pelos valores desejados
Tc = [647.14, 190.56]  # Substitua pelos valores desejados
Pc = [22064000 , 4599000]  # Substitua pelos valores desejados
Vc = [55.95, 98.60]  # Substitua pelos valores desejados
Vr = 0.001    #volume do reator

# Call the LiquidDensityRackett function to get Vm
Vm = liquid_density_rackett(Tp, P, xl, Acentric, Tc, Pc, Vc)

# Volume molar da fase líquida [m³/mol]
VLM = Vm

# Densidade molar da fase líquida [mol/m³]
roLM = 1.0 / VLM  

# Volume da fase líquida [m³]
Vl = n / roLM   

# Concentração da fase líquida [mol/m³]
Cl = n / Vl  

# Para impressão na saída
xLL = xl[1]  # Assuming xL is a list or array in Python (indexing starts from 0)
print(VLM)
#C. I. CALCULADAS
dd = 10e-6  # Diameter in meters
VH1 = (4 * math.pi / 3) * (dd / 2) ** 3  # Volume of a particle with average size in cubic meters
Ad = 4 * math.pi * (dd / 2) ** 2  # Surface area of a particle with average size
kv = 4*math.pi/3

fH = 0.000001       #Fração relativa de hidrato VH/VW
MM1 = 18.015
MM2 = 16.043
ro1 = 997e3
VW = nl1*MM1/ro1   #Volume de água inicial [m³]
VHH = fH*VW            #Volume de Hidrato inicial [m³] 0.005% do Volume de água
mi4 = VHH/kv            #Momento de ordem 3 [m³/m³]
mi1 = VHH/VH1           #Momento de ordem 0 [#/m³]
mi2 = dd*mi1          #Momento de ordem 1 [m/m³]
mi3 = Ad*mi1          #Momento de ordem 2 [m²/m³]

#Volume da fase sólida
Vs = (kv*mi4)*Vl  #[m³]

root=0
xg = [0.0, 0.0]  # Fração molar da fase vapor [mol] - H2O, CH4
xg[1] = 1.0 - xg[0]  # CH4

# Propriedades volumétricas da fase vapor
Vg = Vr - Vs - Vl  # Volume da fase vapor [m³]

# Cálculos com o modelo Peng-Robinson (PR)
root = 0  # Vapor
Vol, FugCoef= PR(Tp, P, xg, Acentric, Tc, Pc, root)  # Substitua pela função de cálculo em Python

# Volume molar da fase vapor [m³/mol]
VgM = Vol
rogM = 1.0 / VgM  # Densidade molar da fase vapor [mol/m³]
 # Coeficiente de fugacidade do CH4 na fase vapor
print(FugCoef)

# Número de mols de CH4 na fase vapor
ng = rogM * Vg  # [mol]

# Imprimindo ou usando os resultados conforme necessário
print(Vg, ng, rogM)
print(roLM)
print(Vs)
print([mi1, mi2, mi3, mi4])

# Condições iniciais diferentes de zero para cada variável
# y0 = np.array([1, 2, 0, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1])
# y0 = np.array([1, 2, mi1, mi2, mi3, mi4, nl1, nl2, 1, 1, 1, 1, 1])  # Valores iniciais diferentes para cada variável
y0 = np.array([Vl,roLM, mi1, mi2, mi3, mi4, nl1, nl2, n, Vs, rogM, ng, 1, Vg])  # Valores iniciais diferentes para cada variável

# Outros parâmetros e configurações permanecem os mesmos
t0 = np.array([500])
yp0 = None  # Não estamos fornecendo derivadas iniciais
par = np.array([1,1,1,2,1,1,VgM,Vr,1,1,1,VLM]) # Não estamos usando parâmetros adicionais neste exemplo

atol = 1e-20
rtol = 1e-12

# Resolvendo o sistema de equações usando dasslc
t, y, yp = dasslc.solve(model15, t0, y0, yp0,par)


# Plotando os resultados
plt.figure(figsize=(10, 6))
plt.plot(t,yp)

plt.xlabel('Tempo')
plt.ylabel('Variáveis de Estado')
plt.title('Solução do Sistema de 15 Equações Diferenciais-Algébricas')
plt.legend()
plt.grid(True)
plt.savefig('index.png', dpi=300)
plt.show()
