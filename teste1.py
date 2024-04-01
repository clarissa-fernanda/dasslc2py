import dasslc
import numpy as np
import matplotlib.pyplot as plt

# Definição da função de modelo com 15 equações diferenciais-algébricas

def model15(t, y, yp):
    res = np.empty(2)   

 
    #VLM=par[0]
    #n=par[1]
    
    Vl=y[0]               #Volume da fase líquida [m³]
    roLM=y[1]             #Densidade molar da fase líquida [mol/m³]
 
    dVldt=yp[0]        #Derivada do volume da fase líquida [m³/s]
    droLMdt=yp[1]      #Derivada da densidade molar da fase líquida [mol/ m³ s]
    
    # Aqui você define as suas 15 equações diferenciais-algébricas
    # Estas são equações hipotéticas apenas para fins de exemplo
    res[0] = roLM - droLMdt # roLM = 1/VLM Por exemplo, yp[i] = y[i+1]
    res[1] = Vl - dVldt #Vl = n/roLM
    
    # A última equação pode ser uma condição de contorno ou uma restrição algébrica

    return res, 0 #----------------- ires can be a literal

# Condições iniciais diferentes de zero para cada variável
y0 = np.array([1,2])  # Valores iniciais diferentes para cada variável

# Outros parâmetros e configurações permanecem os mesmos
t0 = np.array([500])
yp0 = None  # Não estamos fornecendo derivadas iniciais
par = np.array([9.81,0.1]) # Não estamos usando parâmetros adicionais neste exemplo
atol = 1e-10  # Tolerância absoluta
rtol = 1e-10  # Tolerância relativa


# Resolvendo o sistema de equações usando dasslc
t, y, yp = dasslc.solve(model15, t0, y0, yp0)


# Plotando os resultados
plt.figure(figsize=(10, 6))
plt.plot(t,yp)

plt.xlabel('Tempo')
plt.ylabel('Variáveis de Estado')
plt.title('Solução do Sistema de 15 Equações Diferenciais-Algébricas')
plt.legend()
plt.grid(True)
plt.show() 
