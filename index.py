import dasslc
import numpy as np
import matplotlib.pyplot as plt

from liquid_density_rackett import liquid_density_rackett

# Definição da função de modelo com 15 equações diferenciais-algébricas

def model15(t, y, yp,par):
    res = np.empty(13)   

 
    #VLM=par[0]
    #n=par[1]
    #coef_1= par[0]                     #|                           |
    #coef_2 = par[1]  
    #nsol = par[2]  
    #ncr=par[3]
    #kv=par[6]
    #ncr=par[4]
    #mi4=par[8]
    G, coef_1, coef_2, ncr, roHM, kv, VgM, Vg, Aint, kd, neq = par



    Vl, roLM, mi1, mi2, mi3, mi4, nl1, nl2, n, Vs, rogM, ng, nH = y

    
    (dVldt, droLMdt, dmi1dt, dmi2dt, dmi3dt, dmi4dt, 
     dnldt1, dnldt2, dndt, dVsdt, drogMdt, dngdt, dnHdt) = yp
    
    
    #taxas iniciais     
    #ncr = coef_2*roHM*kv*(Vl*dmi4dt + mi4)
    #nsol = (Aint*kd/Vl)*(neq-nl2) # taxa de solubilização de metano 
    ncr = coef_2*roHM*kv*(Vl*dmi4dt + mi4*dVldt) # TAXA MOLAR DE CRESCIMENTO DA FASE HIDRATO  
    nsol = (Aint*kd/Vl)*(neq-nl2)
    ne = rogM + Vg*drogMdt + nsol   # [mol de CH4/s] #TAXA MOLAR DE ENTRADA DE METANO NO REATOR

    #Aint = ((1.d0/H)+((1.d0/Hl)*((Re**1.75)*(We**3.d0)/(Eu**3.d0))*((rosusM/rogM)**-2.67d0)))*Vr  

    n += nl1
    n += nl2
    xl = [nl / n for nl in [nl1, nl2]]
    Tp = 300.0  # Substitua pelo valor desejado
    P = 10.0  # Substitua pelo valor desejado
    Acentric = [0.1, 0.2]  # Substitua pelos valores desejados
    Tc = [400.0, 500.0]  # Substitua pelos valores desejados
    Pc = [10.0, 15.0]  # Substitua pelos valores desejados
    Vc = [0.1, 0.2]  # Substitua pelos valores desejados
    VLM = liquid_density_rackett(Tp, P, xl, Acentric, Tc, Pc, Vc)
    
    # Aqui você define as suas 15 equações diferenciais-algébricas
    # Estas são equações hipotéticas apenas para fins de exemplo
    res[0] = roLM - droLMdt
    # res[0] = roLM - 1/1
    res[1] = Vl - dVldt #Vl = n/roLM
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
    #res[13] = Vg - (Vr-Vs-Vl) #Vg = Vr - Vs - Vl  ![m³]
    
    
    # A última equação pode ser uma condição de contorno ou uma restrição algébrica

    return res, 0 #----------------- ires can be a literal

# Condições iniciais diferentes de zero para cada variável
y0 = np.array([1,2,0,1,2,1,1,1,1,1,1,1,1])  # Valores iniciais diferentes para cada variável

# Outros parâmetros e configurações permanecem os mesmos
t0 = np.array([500])
yp0 = None  # Não estamos fornecendo derivadas iniciais
par = np.array([1,1,1,2,1,1,1,1,1,1,1]) # Não estamos usando parâmetros adicionais neste exemplo


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
