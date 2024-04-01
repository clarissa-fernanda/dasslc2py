import numpy as np

def PR(Tp, P, x, Acentric, Tc, Pc, root):
    NC = len(x)
    a = np.zeros(NC)
    b = np.zeros(NC)
    ac = np.zeros(NC)
    bmix = 0.0
    amix = 0.0

    # Cálculo dos parâmetros a (energético) e b (volumétrico) por componente
    for i in range(NC):
        m = 0.37464 + (1.54226 - 0.26992*Acentric[i])*Acentric[i]
        alfa = (1.0 + m*(1.0 - np.sqrt(Tp/Tc[i])))**2
        ac[i] = 0.45724*Tc[i]**2/Pc[i]
        a[i] = ac[i]*alfa

    for i in range(NC):
        b[i] = 0.07780*Tc[i]/Pc[i]

    # Normalizing mole fractions
    norm = np.sum(x)
    xAux = x / norm

    # Parâmetros de interação binária kij
    BinIntCoef = np.zeros((NC, NC))

    # Cálculo dos parâmetros a (energético) e b (volumétrico) para a mistura
    for i in range(NC):
        for j in range(NC):
            amix += xAux[i]*xAux[j]*np.sqrt(a[i]*a[j])*(1.0 - BinIntCoef[i,j])

    for i in range(NC):
        bmix += xAux[i]*b[i]

    Samix = np.zeros(NC)
    for index in range(NC):
        for i in range(NC):
            Samix[index] += 2.0*xAux[i]*np.sqrt(a[index]*a[i])*(1.0 - BinIntCoef[index,i])

    # Calculo do volume molar da fase root
    sigma_eos = 1.0 + np.sqrt(2.0)
    epsilon_eos = 1.0 - np.sqrt(2.0)

    aux = P / (8.3144621*Tp)

    coefCubic = np.zeros(4)
    coefCubic[0] = 1.0
    coefCubic[1] = (sigma_eos + epsilon_eos - 1)*bmix - 1/aux
    coefCubic[2] = sigma_eos*epsilon_eos*bmix**2 - (1/aux + bmix)*(sigma_eos + epsilon_eos)*bmix + amix/P
    coefCubic[3] = -(1/aux + bmix)*sigma_eos*epsilon_eos*bmix**2 - bmix*amix/P

    roots = np.roots(coefCubic)
    V = roots[roots > 0.0]

    if len(V) == 1:
        Vol = V[0]
    elif len(V) == 2:
        if root == 1:
            Vol = min(V)
        elif root == 0:
            Vol = max(V)
    elif len(V) == 3:
        if root == 1:
            Vol = min(V)
        elif root == 0:
            Vol = max(V)

    # Calculo do coeficiente de fugacidade do componente index
    Z = Vol*aux
    FugCoef = np.zeros(NC)
    for index in range(NC):
        FugCoef[index] = np.exp((b[index]/bmix)*(Z-1.0)-np.log((Vol-bmix)*aux)-amix/(bmix*8.3144621*Tp*(epsilon_eos-sigma_eos))*(Samix[index]/amix-b[index]/bmix)*np.log((Vol+epsilon_eos*bmix)/(Vol+sigma_eos*bmix)))

    return Vol, FugCoef
