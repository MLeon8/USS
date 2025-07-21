import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os

# 1) Crear carpeta de salida para las figuras y el CSV
os.makedirs("figuras", exist_ok=True)

# 2) Definir el modelo SIR (derivadas)
def sir_model(y, t, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

# 3) Datos de las ciudades (semana 1 vs. semana 3)
#    N   = población total
#    I1  = infectados en semana 1
#    I3  = infectados en semana 3 (usados como I0)
#    R3  = recuperados/muertos en semana 3 (usados como R0_ini)
ciudades = {
    "Norilsk":  {"N": 200000, "I1": 300,  "I3": 5000, "R3": 650},
    "Kayerkan": {"N":  80000, "I1": 10,   "I3": 2300, "R3": 1000},
    "Dudinka":  {"N":  60000, "I1": 20,   "I3": 3000, "R3": 900}
}

# 4) Parámetros de simulación
duracion = 10                # semanas totales (de la 3 a la 13 en calendario)
gamma    = 1 / duracion      # tasa de recuperación
t        = np.linspace(0, duracion, 400)  # discretización temporal fina

# 5) Variables para almacenar resultados
infectados_por_ciudad = []  # para el gráfico resumen
resultados = []             # para el CSV final

# 6) Bucle por cada ciudad
for ciudad, datos in ciudades.items():
    N      = datos["N"]
    I1     = datos["I1"]
    I0     = datos["I3"]     # infectados al inicio de la simulación
    R0_ini = datos["R3"]     # recuperados/muertos al inicio
    S0     = N - I0 - R0_ini  # susceptibles al inicio

    # 6.1) Estimar R0 con crecimiento simple: R0 ≈ 1 + (I3 − I1)/I1
    R0_est = 1 + (I0 - I1) / I1
    R0_est = round(R0_est, 2)

    # 6.2) Calcular beta para SIR: beta = R0 × gamma
    beta = R0_est * gamma

    # 6.3) Condiciones iniciales normalizadas [S, I, R]
    y0 = [S0 / N, I0 / N, R0_ini / N]

    # 6.4) Integrar el sistema de ODEs
    sol = odeint(sir_model, y0, t, args=(beta, gamma))
    S, I, R = sol.T

    # 6.5) Guardar la curva de infectados (absolutos) para el resumen
    infectados_por_ciudad.append(I * N)

    # 6.6) Extraer proyecciones al final de la simulación (semana 10)
    infectados_final     = int(I[-1] * N)
    muertos_final        = int(R[-1] * N)
    proporcion_infectada = round(infectados_final / N, 4)

    # 6.7) Acumular resultados para el CSV
    resultados.append({
        "Ciudad":                   ciudad,
        "R0_estimado":              R0_est,
        "Infectados_Semana_10":     infectados_final,
        "Proporcion_Infectados_10": proporcion_infectada,
        "Muertos_Semana_10":        muertos_final
    })

    # 6.8) Graficar SIR para cada ciudad
    plt.figure(figsize=(8, 5))
    plt.plot(t, S * N, 'b', label='Susceptibles')
    plt.plot(t, I * N, 'r', label='Infectados')
    plt.plot(t, R * N, 'g', label='Recuperados/Muertos')
    plt.title(f'Evolución SIR – {ciudad}')
    plt.xlabel('Semanas (desde semana 3)')
    plt.ylabel('Número de personas')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'figuras/SIR_{ciudad}.png')
    plt.close()

# 7) Construir figura resumen (promedio, min y max de infectados)
infectados_arr = np.array(infectados_por_ciudad)
inf_mean = infectados_arr.mean(axis=0)
inf_min  = infectados_arr.min(axis=0)
inf_max  = infectados_arr.max(axis=0)

plt.figure(figsize=(10, 6))
plt.plot(t, inf_mean, 'k', linewidth=2, label='Promedio Infectados')
plt.plot(t, inf_max, 'r--', linewidth=1, label='Máximo')
plt.plot(t, inf_min, 'b--', linewidth=1, label='Mínimo')
plt.fill_between(t, inf_min, inf_max, color='gray', alpha=0.2, label='Rango [min–max]')
plt.title('Proyección Promedio de Infectados (Semanas 3–13)')
plt.xlabel('Semanas (desde semana 3)')
plt.ylabel('Infectados estimados')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("figuras/Resumen_Proyeccion_Total.png")
plt.close()

# 8) Exportar resultados a CSV
df = pd.DataFrame(resultados)
df.to_csv("figuras/proyeccion_ciudades.csv", index=False)
print("✅ Resultados guardados en figuras/proyeccion_ciudades.csv")
print(df)
