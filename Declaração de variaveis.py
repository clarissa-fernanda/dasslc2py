Declaração de variaveis 


DECLARAÇÃO DAS VÁRIAVEIS DE CHAMADA DAS SUBROTINAS

    REAL*8 :: Vol, FugCoef(NC), xLHchute(NC), xLeqH(NC)
    INTEGER :: root
    REAL*8 :: GAMA(NC), IND(NC)
    REAL*8 :: Vm
    !Vol - Volume da fase root 
    !FugCoef - Coeficiente de fugacidade dos componentes NC na fase root
    !xLeqH - Composição de equilíbrio Hidrato-líquido-Vapor a T e P
    !xLchute - Chute inicial para calculo da xLeq
    !root - Fase que se deseja calcular o volume ou coeficiente de fugacidade ( root =  0 - Fase vapor / root = 1 - Fase líquida)
    !IND - Indicador de componente  
    !GAMA - Coeficiente de atividade obtido pelo modelo NRTL
!   !Vm - Volume molar da fase líquida 

 !DECLARAÇÃO DE VARIAVEIS AUXILIARES
!    
    REAL*8 :: xx, SumCoef, xeq0, Nsus, xeq1(NC), tt
    INTEGER :: i
    !xx - Soma das frações molares
    !SumCoef - Soma dos coeficentes estequiometricos
    !xeq0 - Chute inicial da fração molar do CH4 no equilíbrio da interface G-L [mol]
    !Nsus - número de mols da suspensão nH+nl
    !xeq1(NC) - Vetor da composição de equilíbrio na interface para auxiliar o calculo de equilíbrio G-L
    !tt - variavél auxiliar para looping de calculo de equilpibrio G-L
!
!_________________________________________________________________________________________________________________      
    !DECLARAÇÃO DOS DADOS DE ENTRADA
!    
    REAL*8 :: MM(NC), R, pi, Na, kv, ks, vis(NC), roH, NM(NF+1), visL, ro(NC)
    REAL*8 :: Vr, Rin, H, Hl, dag, s, Npo
    REAL*8 :: D, Ea, deltaf
    REAL*8 :: Tc(NC), Pc(Nc), Acentric(NC), Vc(NC)
    REAL*8 :: Aij(NC,NC), Bij(NC,NC), ALPHA
    REAL*8 :: v0 , B(2)
    REAL*8 :: v2, C(3)
    REAL*8 :: P, Tp
    REAL*8 :: VW, VHH, VH1, dd, Ad, fH
!    
    !MM(NC) - massa molar do componente NC [g/mol]
    !R - constante universal dos gases [m³ Pa / mol K = J/ mol K]
    !pí - [rad]
    !Na - nº de avogadro [ moleculas/mol]
    !kv - parâmetro volumétrico de forma para esfera
    !Ks - parâmetros supericial de forma para esfera 
    !vis - viscosidade da água (Componente predominânte na fase líquida) [mol/ m s]
    !roH - densidade mássica da fase hidrato [g/m³]
    !NM(NF+1) - nº de moleculas do componente NF na fase hidrato
    !Vr - Volume do reator [m³]
    !Rin - Raio interno do reator [m]
    !H - Altura do reator [m]
    !Hl - Altura do coluna de líquido [m] 
    !dag - Diametro do agitador [m]
    !s - taxa de agitação [rps]
    !Npo - Número de potencia do agitador 
    !D - Difusividade do metano em água [m²/s]
    !Ea - Energia de ativação do processo de difusão do metano em água [J/mol]
    !Deltaf - Expessura do filme em torno da partícula [m] - Mi(2)/100 (Aproximação)
    !Tc(NC) - Temperatura crítica do componente NC [K]
    !Pc(NC) - Pressão critica do componente NC [Pa]
    !Vc(NC) - Volume critico do componente NC [m³/mol]
    !Acentric(NC) - Fator acentrico do componente NC
    !Aij - Parametro independente de T NRTL
    !Bij - Parametro dependente de T NRTL
    !ALPHA - Parametro interativo 
    !Parâmetros da correlação da constante de Henry a diluição infinita para chute inicial da composição de equilíbrio do metano na interface G-L
    !v0 - [m³/mol]
    !B(1) - [Admensional]
    !B(2) - [T]
    !Parâmetros da correlação da fase hipotética de metano puro
    !v2 - [m³/mol]
    !C(1) - [Admensional]
    !C(2) - [Admensional]
    !C(3) - [Admensional]
    !P - Pressão do sistema [Pa]
    !Tp - Temperatura do sistema [K]
    !VW - Volume de água do sistema [m³]
    !VHH - Volume de hidrato inicial do sistema (10% do volume de água) [m³]
    !VH1 - Volume de um cristal de hidrato de tamanho médio [m3]
    !dd - Tamanho médio do cristal de hidrato (10 micrometros)  [m]
    !Ad - Area superficial de um cristal de hidrato de tamanho médio [m²]
    !fH - fração de hidrato na água VH/VW [Admensional]   
!
_______________________________________________________________________________________    
    !DECLARAÇÃO DAS VARIAVEIS NÃO DIFERENCIAVEIS DO PROBLEMA
!  
    REAL*8 :: ne, neq
    REAL*8 :: xeq(NC), xL(NC), xG(NC), FiV(NC), FiL0(NC), Cl
    REAL*8 :: aeq(NF+1), ab(NF+1), Asup, coef(NF+1), GamaB(NF+1), GamaEq(NF+1), Kdd(NF+1), bm(NF+1), bb(NF+1), beq(NF+1), Da(NF+1)
    REAL*8 :: MMH, Ld, Aint, kd, rosusM, visusM, MMsus, dsus, Re, We, Eu
    REAL*8 :: roHM, VLM, VGM, TETA(NF+1), rr, YY(2)
!    
!    
    !ncr - Taxa molar de crescimento da fase hidrato [mol/s] Declarado VarRes
    !nsol - Taxa molar de solubilização do metano [mol/s] Declarado VarRes
    !ne - Taxa molar de entrada de metano no reator [mol/s]
    !neq - Número de mols de metano em equilíbrio na interface G-L [mol]
    !xeq(NC) - Fração molar do componente NC no equilíbrio da interface G-L [mol]
    !G - Taxa volumétrica de crescimento da fase hidrato [m³/mol] Declarado VarRes
    !Keq - Parâmetro do produto das atividades na condição de equilíbrio [admensional]
    !Kb - Parâmetro do produto das atividades na condição do bulk [admensional]
    !aeq(NF) - Atividade na condição de equilíbrio do formador NF [admensional]
    !ab(NF) - Atividade na condição do bulk do formador NF [admensional]
    !Asup - Aréa superficial total das partículas de hidrato [m²]
    !coef(NF+1) - Coeficiente estequiométrico da água (1) e dos componentes formadores de hidrato NF
    !Kdd - Paramêtro de difusão dependente da composição [mol/m² s]
    !bm - Composição média ponderada
    !bb - Composição ponderada na condição do bulk
    !beq - Composição ponderada na condição do equilibrio HL
    !DR - Parametro de relação Difusão-reação
    !MMH - Massa molar do hidrato [g/mol]
    !Aint - Aréa interfacial entre as fases vapor e líquida [m²]
    !rosusM - Densidade molar da suspensão [mol/m³]
    !visusM - Viscosidade molar da suspensão [mol/ m s]
    !MMsus - Massa molar da suspensão [g/mol]
    !dsus - Energia especifica superficial da suspensão [J/ m²]
    !Re - Reynolds
    !We - Weber
    !Eu - Euler
    !kd - Coeficiente de difusão mássica de metano na interfase G-L [m/s]
    !xL(NC) - Fração molar do componente NC na fase líquida
    !xG(NC) - Fração molar do componente NC na fase vapor (Considerada pura em metano xG(2)=1)
    !FiV(NC) - Coeficiente de fugacidade do componente NC na fase vapor  
    !Ld - Coeficiente de difusão global do sistema [m²/s]
    !GamaB(NC) - Coeficiente de atividade do componente NC na condição do bulk
    !GamaEq(NC) - Coeficiente de atividade do componente NC na condição de equilíbrio
    !FiL0(NC) - Coeficiente de fugacidade do componente NC puro na fase líquida
    !x0(NC) - Fração molar para cálculo da fugacidade do componente NC puro na fase líquida
    !roHM - Densidade molar da fase hidrato [mol/m³]
    !VLM - Volume molar da fase líquida [m³/mol]
    !VGM - Volume molar da fase vapor [m³/mol]
    !TETA(NF+1) - Fator de Ocupação do formador no Hidrato
    !rr - Coeficiente de reação global do sistema [mol/m² s]
    !Cl - concentração da fase líquida [mol/m³]
    !YY - Fator de ocupação de cada cavidade

# PARÂMETROS CALCULADOS
# ÁREA INTERFACIAL GÁS - LÍQUIDO
# Propriedades do sistema

# Altura do reator [m]
H = Vr / (pi * Rin * Rin)

# Altura da coluna de líquido [m]
Hl = (Vl + Vs) / (pi * Rin * Rin)

# Viscosidade da fase líquida [mol/s m]
visL = sum(xL[i] * vis[i] for i in range(1, NC + 1))

# Propriedades da suspensão
dsus = 1.0e-3 * ((2.848 * (roLM - rogM)**2 - 8.056 * (roLM - rogM) + (8.53 * (roLM - rogM) / (Tp / Tc[2])**0.3125))**4)
rosusM = (n * roLM + nH * roHM) / (n + nH)
visusM = visL * math.exp(2.5 * (Vs / Vr) / (1.0 - 1.45 * (Vs / Vr)))

# Massa molar da suspensão [g/mol]
Nsus = nH + n
MMsus = nH * MMH / Nsus
MMsus += sum(nl[i] * MM[i] / Nsus for i in range(1, NC + 1))

# Parâmetros adimensionais do sistema
Re = Dag**2 * rosusM * s / visusM  # Reynolds
We = Dag**3 * rosusM * s**2 / dsus  # Weber
Eu = P / (Dag**2 * rosusM * MMsus * s**2)  # Euler

# Área interfacial G-L
Aint = ((1.0 / H) + (1.0 / Hl) * (Re**1.75 * We**3 / Eu**3 * (rosusM / rogM)**-2.67)) * Vr  # [m²]

# COEFICIENTE DE DIFUSÃO MÁSSICA DE METANO NA INTERFACE G-L
kd = 0.5 * D * (Npo * s**3 * Dag**5 / (visusM / rosusM))**(1.0 / 4.0)  # [m/s]

_______________________________________________________________________________________________________________      
    !DECLARAÇÃO DAS VARIAVEIS DIFERENCIAVEIS DO PROBLEMA
!    
    REAL*8 :: nl(NC), n, ng, nH, Vl, roLM, Vs, Vg, rogM, mi4
    !nl - número de mols do componente NC na fase líquida [mol]
    !n - número de mols da fase líquida [mol]
    !ng - número de mols de metano na fase vapor [mol] (Fase vapor considerada pura)
    !nH - número de mols da fase hidrato [mol]
    !Vm - Volume molar da fase líquida [m³/mol]
    !Vl - volume da fase líquida [m³]
    !Vs - volume da fase hidrato [m³]
    !Vg - volume da fase vapor [m³]
    !roLM - densidade molar da fase líquida [mol/m³]
    !rogM - densidade molar da fase vapor [mol/m³]
    !mi(j) - momento de ordem j [#/m³], [m/m³], [m²/m³], [m³/m³]
!    
!_________________________________________________________________________________________________________________        
    !DECLARAÇÃO DAS DERIVADAS DAS VARIAVEIS DIFERENCIAVEIS
!    
    REAL*8 :: dnldt(NC), dndt, dngdt, dnHdt, dVldt, droLMdt, dVsdt, dVgdt, drogMdt, dmi4dt
!    
    !dnldt(NC) - Derivada temporal do número de mols do componente NL na fase líquida [mol/s]
    !dndt - Derivada temporal do número de mols da fase líquida [mol/s]
    !dngdt - Derivada temporal do número de mols da fase vapor [mol/s]
    !dnHdt - Derivada temporal do número de mols da fase hidrato [mol/s]
    !dVldt - Derivada temporal do volume da fase líquida [m³/s]
    !droLdt - Derivada temporal da densidade molar da fase líquida [mol/m³s]
    !dVsdt - Derivada temporal do volume da fase sólida [m³/s]
    !dVgdt - Derivada temporal do volume da fase vapor [m³/s]
    !drogdt -  Derivada temporal da densidade molar da fase vapor [mol/m³s]
    !dmidt(j) - Derivada temporal do momento de ordem j [#/m³s], [m/m³s], [m²/m³s], [m³/m³s]

#dados de entrada

import math

# Atribuindo valores dos dados de entrada do problema
# H2O - (NC = 1, index = 1), CH4 - (NC = 2, index = 2)
# Fase vapor - root = 0
# Fase líquida - root = 1

P = 7.09e6  # Pressão do sistema [Pa]
Tp = 276.0  # Temperatura do sistema [K]

# Massa molar do componente NC [g/mol]
MM = {1: 18.015, 2: 16.043}
R = 8.3144621  # Constante universal dos gases [m³ Pa / mol K = J/ mol K]
pi = 4.0 * math.atan(1.0)  # π [rad]
roH = 917.8  # Densidade mássica da fase hidrato [g/m³]
ro1= 997.3           #!Densidade da H2O [g/m³] 
ro2 = 0.6563          #!Densidade do CH4 [g/m³]
roH = 917.83           # !Densidade mássica da fase hidrato [g/m³]

# Viscosidade dos componentes [mol/ m s]
vis = {
    1: (1.0e3 / MM1) * math.exp(-3.7188 + (578.918 / (Tp - 137.546))),
    2: (1.0e3 / MM2) * math.exp(-25.5947 + (25392.0 / (Tp + 969.306)))
}
vis1 = (1./MM1)*math.exp(-3.71880 + (578.9180 / (Tp-137.5460)))#!Viscosidade da H2O [mol/ m s]
vis2 = (1/MM2*math.exp(-25.59470 + (25392.0 / (Tp+969.3060))) #!Viscosidade da CH4 [mol/ m s]

Na = 6.0221409e23  # Nº de Avogadro [moleculas/mol]

kv = 4.0 * pi / 3.0  # Parâmetro volumétrico de forma para esfera
ks = 4.0 * pi  # Parâmetro superficial de forma para esfera

# Nº de moléculas de cada componente
NM1= 184.0
NM2= 8.0

Vr = 1.5e-3  # Volume do reator [m³]
Rin = 0.054  # Raio interno do reator [m]
D = 0.0347e-4  # Difusividade do metano em água [m²/s]
Ea = 18.36e3  # Energia de ativação do processo de difusão do metano em água [J/mol]
Dag = (2.0 / 3.0) * Rin  # Diâmetro do agitador [m]
s = 400.0  # Taxa de agitação [rps]
Npo = 200.0  # Número de potência do agitador
Deltaf = 1.0e-6  # Espessura do filme em torno da partícula [m]

# Temperatura crítica do componente NC [K]
Tc = {1: 647.14, 2: 190.56}

# Pressão crítica do componente NC [Pa]
Pc = {1: 22064000.0, 2: 4599000.0}

# Fator acentrico do componente NC
Acentric = {1: 0.344, 2: 0.011}

# Vc(NC) - Volume crítico do componente NC [cm³/mol]
Vc = {1: 55.95, 2: 98.60}  # Multiplicado por 1e-6 para converter para m³/mol

# Parâmetros da correlação da constante de Henry a diluição infinita
v0 = 32.0e-12  # [m³/mol]
B = {1: 15.826277, 2: -1559.063}

# Parâmetros da correlação da fugacidade da fase líquida hipotética do CH4
v2 = 5.2e-5  # [m³/mol]
C = {1: 44.4895, 2: -1.4457, 3: -2.6880}

# Parâmetros NRTL
Aij = [[0.0, -208.76], [879.302, 0.0]]  # Linear
Bij = [[0.0, -41.8913], [-430.631, 0.0]]  # Não Linear

ALPHA = 0.00160748  # Fator de não aleatoriedade
rr = 2.2e-4  # Coeficiente de reação global do sistema [mol/m² s]

# Parametros da fase sólida para cálculo das C.I.
dd = 10.0e-6  # Diametro médio dos cristais de hidrato [m]
VH1 = (4.0 * pi / 3.0) * (dd / 2.0) ** 3.0  # Volume de uma partícula de tamanho médio [m³]
Ad = 4.0 * pi * (dd / 2.0) ** 2.0  # Área superficial de uma partícula de tamanho médio


# Condições iniciais das variáveis (C.I.)
nl1 = 16.6 #H20 na fase liquida
nl2= 1.0e-5 # Número de mols DE CH4 na fase líquida [mol]
fH = 0.000001  # Fração relativa de hidrato VH/VW

# C.I. Calculadas
VW = nl1* MM1 / ro[1]  # Volume de água inicial [m³]
VHH = fH * VW  # Volume de Hidrato inicial [m³] (0.005% do Volume de água)
mi4= VHH/kv            #!Momento de ordem 3 [m³/m³]
mi1 = VHH/VH1           #!Momento de ordem 0 [#/m³]   
mi2 = dd*mi1          #Momento de ordem 1 [m/m³]
mi3 = Ad*mi1          #Momento de ordem 2 [m²/m³]

# Condições iniciais CALCULADAS DAS VARIÁVEIS
# FASE LÍQUIDA

# Número de mols da fase líquida [mol]
n = 0.0
for i in range(1, NC + 1):
    n += nl[i]

# Fração molar do componente i na fase líquida
xl = {i: nl[i] / n for i in range(1, NC + 1)}

# Verificar se a soma das frações molares é 1
xx = 0.0
for i in range(1, NC + 1):
    xx += xl[i]

if xx = 1.0:
    print('Soma das frações molares da fase líquida diferente de 1 - VERIFICAR')  # Verificação
    exit()
    
G = ((Asup*rr)/(1.d0 + Da(1) + Da(2)))*abs(dlog(Kb/Keq))
# Propriedades volumétricas da fase líquida
# Supondo que você tenha uma função LiquidDensityRackett definida em algum lugar
# Passando parâmetros para a função
Vm = LiquidDensityRackett(Tp, P, xl, Acentric, Tc, Pc, Vc, NC)
VLM = Vm  # Volume molar da fase líquida [m³/mol]
roLM = 1.0 / VLM  # Densidade molar da fase líquida [mol/m³]
Vl = n / roLM  # Volume da fase líquida [m³]
Cl = n / Vl  # Concentração da fase líquida [mol/m³]
xLL = xl[2]  # Para impressão na saída
print(Vm, roLM, Vl)
# input()

#EQUILÍBRIO DA INTERFACE SÓLIDO-LÍQUIDO (S-L)
# EQUILÍBRIO DA INTERFACE SÓLIDO-LÍQUIDO (S-L)
# Condições de Equilíbrio
# Cálculo da composição de equilíbrio dado Tp e P
xLHchute = {1: xL[1], 2: 1.0 - xL[1]}  # H2O e CH4

# Chamar a sub-rotina (ou função) xLEq para calcular a composição de equilíbrio
xLeqH = xLEq(Tp, P, xLHchute, Acentric, Tc, Pc, NC, Aij, Bij, ALPHA, C, v2)
xLeqHH = xLeqH[2]  # CH4 para impressão na saída

# FASE SÓLIDA
# Volume da fase sólida
Vs = kv * mi[4] * Vl  # [m³]

# Estequiometria da fase sólida
# Fator de Ocupação
TETA = {1: 1.0, 2: (YY[1] + YY[2]) / 2.0}  # H2O sempre 100%, CH4

coef = {i: TETA[i] * NM[i] for i in range(1, NC + 1)}

# Massa molar da fase hidrato [g/mol]
SumCoef = sum(coef.values())
MMH = sum((coef[i] * MM[i]) / SumCoef for i in range(1, NC + 1))  # [g/mol]

# Densidade molar da fase sólida
roHM = roH / MMH  # [mol/m³]

# Número de mols da fase sólida
nH = roHM * Vs

# FASE VAPOR
# Fração molar da fase vapor - CH4 puro
xg = {1: 0.0, 2: 1.0}  # H2O e CH4

# Propriedades volumétricas da fase vapor
# Volume da fase vapor
Vg = Vr - Vs - Vl  # [m³]

## chamar PR PARA CALCULAR O VG, NG E ROG
call PR(Tp, P, xg, Acentric, Tc, Pc, NC, root, Vol, FugCoef)  #PengRobinson
!     
       VgM = Vol         #Volume molar da fase vapor [m³/mol]
       rogM = 1.d0/VgM   #Densidade molar da fase vapor [mol/m³] 
       FIV(2) = FugCoef(2)     #Coeficiente de fugacidade do CH4 na fase vapor
!       
     !Numero de mols de CH4 na fase vapor
       ng = rogM*Vg  ![mol]
       
       !Write(*,*)Vg, ng, rogM
        !pause

## CHAMAR CALCULO DE EQUILIBRIO PROGRAMA DO IURI
## CHAMAR NRTL
##


# Constants and initial conditions
coef1 = 1.0 
coef2= 0.5 # coef(1) and coef(2)
ncr = 20
nsol = 100  # Assuming nsol is constant for simplicity
ne = 120
kv = 0.1  # Example value for kv
roHM = 1.0  # Example density value for roHM
mi4 = 0.5  # Example value for mi4
Vl = 1.0  # Assuming Vl is constant for simplicity
dmi4dt = 0.2  # Example value for dmi4dt
dVldt 4=   # Example value for dVldt

#Equações Dasslc 


#balanço de massa 
dn_dt = nsol - ((coef[0] / coef[1]) + 1) * ncr
ng = rogM*Vg 
nH = roHM*Vs
dng_dt = ne - nsol
dnH_dt = kv * roHM * (Vl * dmidt4 + mi4 * dVldt)
dnLdt1 = - (coef[0] / coef[1]) * ncr #H2O
dnLdt2 = nsol - ncr #METANO
nl1 = 16.6 3 #H20 na fase liquida
nl2= 1.0e-5 # Número de mols DE CH4 na fase líquida [mol]


#volume das fases 
Vl = n / roLM  # Volume da fase líquida [m³]
roLM = 1.0 / VLM #  densidade da fase liquida
# Volume da fase sólida
Vs = kv * mi4 * Vl  # [m³]

# Volume da fase vapor
Vg = Vr - Vs - Vl  # [m³]
VgM = Vol         #Volume molar da fase vapor [m³/mol]
rogM = 1/VgM   #Densidade molar da fase vapor [mol/m³] 
# momentos 
VHH = fH*VW            #Volume de Hidrato inicial [m³] 
Ad = 4*pi*(dd/2)**2     #Diametro médio dos cristais de hidrato
dd = 10*10-6     #Diametro médio dos cristais de hidrato
mi4 = VHH/kv            #Momento de ordem 3 [m³/m³]
mi1 = VHH/VH1           #Momento de ordem 0 [#/m³]   
mi2 = dd*mi1        #Momento de ordem 1 [m/m³]
mi3 = Ad*mi1          #Momento de ordem 2 [m²/m³]


#DEFININDO DERIVADAS TEMPORAIS E TAXAS INICIAIS PARA O MÉTODO NUMÉRICO
#NA AUSÊNCIA DE ESTIMATIVAS INICIAIS, INSIRA UM VALOR APROXIMADO QUE A DASSL TENTA ESTIMÁ-LA
# taxas iniciais 
#taxas molares 
nsol = (Aint*kd/Vl)*(neq-nl(2)) # taxa de solubilização de metano 
ncr = coef(2)*roHM*kv*(Vl*dmi4dt + mi4*dVldt) # TAXA MOLAR DE CRESCIMENTO DA FASE HIDRATO  
ne = rogM*dVgdt + Vg*drogMdt + nsol   # [mol de CH4/s] #TAXA MOLAR DE ENTRADA DE METANO NO REATOR

#balanço populacional 
dmidt1 = 0       
dmidt(2) = G*mi1      #[m/m³]
dmidt(3) = 2.G*mi2 #[m²/m³] 
dmi4dt = 3*G*mi3 #[m³/m³]

# DERIVADAS TEMPORAIS DAS DENSIDADES MOLARES E DOS VOLUMES DAS FASES

droLMdt = 0
dVldt = ((nsol/roLM)-((SumCoef*kv*roHM/roLM)*((Vl*dmi4dt)))+((n/(roLM*roLM))*droLMdt))/(1+(SumCoef*kv*roHM*mi4/roLM))
dVsdt = kv*(mi4*dVldt+Vl*dmi4dt)
dVgdt = - dVsdt - dVldt
drogMdt = 0

#DERIVADAS TEMPORAIS DO BALANÇO DE MASSA   
dnldt(1) = -(coef(1)/coef(2))*ncr
dnldt(2) = (nsol - ncr)
dndt = nsol - ((coef(1)/coef(2))+1)*ncr
dngdt = ne - nsol
dnHdt = kv*roHM*((Vl*dmi4dt)+(mi4*dVldt))

return [dn_dt, dng_dt, dnH_dt, dnLdt1, dnLdt2]

# Define the system of differential equations

def system_of_eqs(t, y):
    dn_dt = nsol - ((coef[0] / coef[1]) + 1) * ncr
    dng_dt = ne - nsol
    dnH_dt = kv * roHM * (Vl * dmidt4 + mi4 * dVldt)
    dnLdt1 = - (coef[0] / coef[1]) * ncr
    dnLdt2 = nsol - ncr
    return [dn_dt, dng_dt, dnH_dt, dnLdt1, dnLdt2]


