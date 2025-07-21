import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os

# 1) Crear carpeta de salida para figuras y CSV
os.makedirs("figuras", exist_ok=True)

# 2) Definir el modelo SIR
def sir_model(y, t, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

# 3) Datos de evolución (Semanas 1, 2 y 3)
#    Tabla 1: Susceptibles, Infectados y Muertos
#    (Observación: S permanece ~constante en los primeros 3 semanas)
ciudades = {
    "Norilsk":  {
        "N": 200000,
        "I1": 300,  "R1": 0,
        "I2": 1500, "R2": 300,
        "I3": 5000, "R3": 650
    },
    "Kayerkan": {
        "N":  80000,
        "I1": 10,   "R1": 0,
        "I2": 700,  "R2": 120,
        "I3": 2300, "R3": 1000
    },
    "Dudinka":  {
        "N":  60000,
        "I1": 20,   "R1": 0,
        "I2": 1200, "R2": 60,
        "I3": 3000, "R3": 900
    }
}

# 4) Parámetros de simulación
duracion = 10                  # simulación de 10 semanas
gamma = 1 / duracion           # tasa de recuperación
t = np.linspace(0, duracion, 400)

# 5) Listas para resultados
infectados_por_ciudad = []     # curvas para figura resumen
resultados = []                # datos para CSV

# 6) Bucle de simulación por ciudad
for ciudad, d in ciudades.items():
    N      = d["N"]
    I1     = d["I1"]           # infectados en semana 1
    I0     = d["I3"]           # infectados en semana 3 → condición inicial
    R0_ini = d["R3"]           # muertos/recuperados en semana 3
    S0     = N - I0 - R0_ini   # susceptibles en t=0 (semana 3)

    # 6.1) Estimación de R₀: crecimiento exponencial entre sem1 y sem3
    R0_est = round(1 + (I0 - I1) / I1, 2)
    beta   = R0_est * gamma

    # 6.2) Condiciones iniciales normalizadas
    y0 = [S0 / N, I0 / N, R0_ini / N]

    # 6.3) Integrar SIR
    sol = odeint(sir_model, y0, t, args=(beta, gamma))
    S, I, R = sol.T

    # 6.4) Curva absoluta de infectados para resumen
    infectados_por_ciudad.append(I * N)

    # 6.5) Proyección a semana 10
    infectados_10 = int(I[-1] * N)
    muertos_10    = int(R[-1] * N)
    prop_inf_10   = round(infectados_10 / N, 4)

    resultados.append({
        "Ciudad":                   ciudad,
        "R0_estimado":              R0_est,
        "Infectados_Semana_10":     infectados_10,
        "Proporcion_Infectados_10": prop_inf_10,
        "Muertos_Semana_10":        muertos_10
    })

    # 6.6) Graficar S, I, R por ciudad
    plt.figure(figsize=(8,5))
    plt.plot(t+3, S*N, 'b', label='Susceptibles')
    plt.plot(t+3, I*N, 'r', label='Infectados')
    plt.plot(t+3, R*N, 'g', label='Recuperados/Muertos')
    plt.title(f'Evolución SIR – {ciudad}')
    plt.xlabel('Semana')
    plt.ylabel('Personas')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'figuras/SIR_{ciudad}.png')
    plt.close()

# 7) Figura resumen: promedio, min y max de infectados
arr = np.array(infectados_por_ciudad)
mean_inf = arr.mean(axis=0)
min_inf  = arr.min(axis=0)
max_inf  = arr.max(axis=0)

plt.figure(figsize=(10,6))
plt.plot(t+3, mean_inf, 'k', lw=2, label='Promedio Infectados')
plt.plot(t+3, max_inf, 'r--', lw=1, label='Máximo')
plt.plot(t+3, min_inf, 'b--', lw=1, label='Mínimo')
plt.fill_between(t+3, min_inf, max_inf, color='gray', alpha=0.2, label='Rango [min–max]')
plt.title('Proyección Promedio de Infectados (Semanas 3–13)')
plt.xlabel('Semana')
plt.ylabel('Infectados')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figuras/Resumen_Proyeccion_Total.png")
plt.close()

# 8) Guardar CSV con resultados
df = pd.DataFrame(resultados)
csv_path = "figuras/proyeccion_ciudades.csv"
df.to_csv(csv_path, index=False)
print("✅ Resultados guardados en", csv_path)
print(df)
