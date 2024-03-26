def liquid_density_rackett(TT, P, x, Acentric, Tc, Pc, Vc):
    R = 8.3144621  # Constante universal dos gases [m³ Pa / mol K = J/ mol K]
    NC = len(x)
    
    # Temperatura Reduzida da Mistura
    KK = []
    for i in range(NC):
        row = []
        for j in range(NC):
            kij = 1.0 - (8.0 * (Vc[i]*Vc[j])**(1.0/2.0)) / ((Vc[i]**(1.0/3.0) + Vc[j]**(1.0/3.0))**(3.0))
            row.append(kij)
        KK.append(row)
    
    Tcc = []
    for i in range(NC):
        row = []
        for j in range(NC):
            tcc_ij = (1.0 - KK[i][j]) * ((Tc[i]*Tc[j])**(1.0/2.0))
            row.append(tcc_ij)
        Tcc.append(row)
    
    # Parâmetro Fi
    SumFi = sum(x[i]*Vc[i] for i in range(NC))
    Fi = [(x[i]*Vc[i]) / SumFi for i in range(NC)]
    
    Tcm = 0.0
    for i in range(NC):
        for j in range(NC):
            Tcm += Fi[i]*Fi[j]*Tcc[i][j]
    
    # Temperatura Reduzida
    Trm = TT / Tcm
    
    # Volume molar da mistura
    ZRA = [0.29056 - (0.08775 * Acentric[i]) for i in range(NC)]
    
    ZRAm = 0.0
    for i in range(NC):
        ZRAm += x[i]*ZRA[i]
    
    SumTc = 0.0
    for i in range(NC):
        SumTc += x[i]*Tc[i]/(Pc[i])
    
    Vm = R * SumTc * (ZRAm**(1.0 + (((1.0 - Trm))**(2.0/7.0))))
    
    return Vm
