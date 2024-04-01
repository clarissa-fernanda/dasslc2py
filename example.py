import dasslc
import numpy as np
import matplotlib.pyplot as plt

# Definindo as constantes
G = 9.81  # Ajuste conforme necessário
Aint = 1.0  # Ajuste conforme necessário
kd = 1.0  # Ajuste conforme necessário
Vl = 1.0  # Ajuste conforme necessário
SumCoef = 1.0  # Ajuste conforme necessário
kv = 1.0  # Ajuste conforme necessário
rogM = 1.0  # Ajuste conforme necessário
roHM = 1.0  # Ajuste conforme necessário
coef = np.array([1.0, 1.0])  # Ajuste conforme necessário

# Condições iniciais
mi_initial = np.array([0.0, 0.0, 0.0, 0.0])
nl_initial = np.array([0.0, 0.0])
Vl_initial = 0.0
Vs_initial = 0.0
Vg_initial = 0.0
rogM_initial = 0.0

y0 = np.concatenate((mi_initial, nl_initial, [Vl_initial, Vs_initial, Vg_initial, rogM_initial]))

# Definindo as equações diferenciais e algébricas
def model(t, y, yp):
    mi = y[:4]
    nl = y[4:6]
    Vl, Vs, Vg, rogM = y[6:]

    dmidt = np.array([
        0.0,
        G * mi[0],
        2.0 * G * mi[1],
        3.0 * G * mi[2]
    ])

    droLMdt = 0.0
    dVldt = ((nsol / roLM) - ((SumCoef * kv * roHM / roLM) * (Vl * dmidt[3])) + (
            (n / (roLM * roLM)) * droLMdt)) / (1 + SumCoef * kv * roHM * mi[3] / roLM)
    dVsdt = kv * (mi[3] * dVldt + Vl * dmidt[3])
    dVgdt = -dVsdt - dVldt
    drogMdt = 0.0

    ne = rogM * dVgdt + Vg * drogMdt + nsol
    ncr = coef[1] * roHM * dVsdt

    dnldt = np.array([
        -(coef[0] / coef[1]) * ncr,
        nsol - ncr
    ])

    dndt = nsol - ((coef[0] / coef[1]) + 1) * ncr
    dngdt = ne - nsol
    dnHdt = kv * roHM * ((Vl * dmidt[3]) + (mi[3] * dVldt))

    return np.concatenate((dmidt, dnldt, [dVldt, dVsdt, dVgdt, drogMdt, dnldt[0], dnldt[1], dndt, dngdt, dnHdt]))

# Condições do problema
t0 = np.linspace(0, 1000, 100)
neq=5
nsol = (Aint * kd / Vl) * (neq - nl[1])  # Ajuste conforme necessário
roLM = 1.0  # Ajuste conforme necessário
n = 1.0  # Ajuste conforme necessário

# Resolvendo o sistema
t, y, yp = dasslc.solve(model, t0, y0)

# Plotando os resultados se necessário
# (Lembre-se de ajustar conforme a quantidade de variáveis no seu sistema)
plt.plot(t, y[:, 0], label='mi1')
plt.plot(t, y[:, 1], label='mi2')
plt.plot(t, y[:, 2], label='mi3')
plt.plot(t, y[:, 3], label='mi4')
plt.plot(t, y[:, 6], label='Vl')
plt.plot(t, y[:, 7], label='Vs')
plt.plot(t, y[:, 8], label='Vg')
plt.xlabel('Tempo')
plt.ylabel('Variáveis do Modelo')
plt.legend()
plt.show()
