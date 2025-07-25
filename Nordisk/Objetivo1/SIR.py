import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Datos iniciales
data = {
    "Norilsk": {
        "S": [200000, 200000, 200000],
        "I": [300, 1500, 5000],
        "D": [0, 300, 650],
        "N": 200300  # Población inicial (S + I)
    },
    "Kayerkan": {
        "S": [80000, 80000, 80000],
        "I": [10, 700, 2300],
        "D": [0, 120, 1000],
        "N": 80010
    },
    "Dudinka": {
        "S": [60000, 60000, 60000],
        "I": [20, 1200, 3000],
        "D": [0, 60, 900],
        "N": 60020
    }
}

# Cálculo básico de R0
def calcular_R0(ciudad):
    I1 = data[ciudad]["I"][0]
    I3 = data[ciudad]["I"][2]
    return 1 + (I3 - I1)/I1

# Proyección lineal simple a semana 10
def proyectar_semana10(ciudad):
    # Tasa de crecimiento semanal de infectados
    tasa_crecimiento = (data[ciudad]["I"][2] - data[ciudad]["I"][0])/2
    
    # Proyección lineal (simplificación inicial)
    I10 = data[ciudad]["I"][0] + tasa_crecimiento*9
    
    # Estimación de muertos basada en CFR promedio
    cfr = data[ciudad]["D"][2]/data[ciudad]["I"][2]
    D10 = I10 * cfr
    
    return int(I10), int(D10)

# Resultados básicos
resultados = {}
for ciudad in data.keys():
    R0 = calcular_R0(ciudad)
    I10, D10 = proyectar_semana10(ciudad)
    resultados[ciudad] = {"R0": R0, "I10": I10, "D10": D10}
    
    print(f"{ciudad}:")
    print(f"R0 estimado: {R0:.2f}")
    print(f"Proyección semana 10 - Infectados: {I10}, Muertos: {D10}")
    print("----------------------")

# Guardar resultados básicos
df_resultados = pd.DataFrame(resultados).T
df_resultados.to_csv("resultados_basicos.csv")

# Gráficos básicos de evolución
for ciudad in data.keys():
    plt.figure(figsize=(10, 5))
    
    # Datos observados
    semanas = [1, 2, 3]
    plt.plot(semanas, data[ciudad]["I"], 'ro-', label='Infectados observados')
    plt.plot(semanas, data[ciudad]["D"], 'ks-', label='Muertos observados')
    
    # Proyección lineal
    semanas_proy = [1, 2, 3, 10]
    infectados_proy = [data[ciudad]["I"][0], data[ciudad]["I"][1], 
                      data[ciudad]["I"][2], resultados[ciudad]["I10"]]
    muertos_proy = [data[ciudad]["D"][0], data[ciudad]["D"][1], 
                   data[ciudad]["D"][2], resultados[ciudad]["D10"]]
    
    plt.plot(semanas_proy, infectados_proy, 'r--', label='Proyección infectados')
    plt.plot(semanas_proy, muertos_proy, 'k--', label='Proyección muertos')
    
    plt.title(f"Evolución en {ciudad}\nR0 = {resultados[ciudad]['R0']:.2f}")
    plt.xlabel("Semanas")
    plt.ylabel("Número de casos")
    plt.legend()
    plt.grid()
    plt.savefig(f"evolucion_basica_{ciudad}.png", dpi=300)
    plt.close()

    # Continuación del script anterior - Parte de modelado avanzado
from scipy.integrate import odeint

def modelo_SIRD(y, t, N, beta, gamma, mu):
    S, I, R, D = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I - mu * I
    dRdt = gamma * I
    dDdt = mu * I
    return dSdt, dIdt, dRdt, dDdt

def ajustar_y_proyectar(ciudad):
    # Datos de la ciudad
    N = data[ciudad]["N"]
    I0 = data[ciudad]["I"][0]
    D0 = data[ciudad]["D"][0]
    S0 = N - I0
    R0 = 0
    y0 = S0, I0, R0, D0
    
    # Semanas observadas (1-3)
    t_obs = np.array([0, 1, 2])  # semanas desde inicio
    
    # Función para ajustar parámetros
    def error(params):
        beta, gamma, mu = params
        sol = odeint(modelo_SIRD, y0, t_obs, args=(N, beta, gamma, mu))
        I_modelo = sol[:, 1]
        D_modelo = sol[:, 3]
        error_I = np.sum((I_modelo - data[ciudad]["I"])**2)
        error_D = np.sum((D_modelo - data[ciudad]["D"])**2)
        return error_I + error_D
    
    # Ajuste de parámetros
    from scipy.optimize import minimize
    params0 = [0.5, 0.1, 0.01]  # Valores iniciales para beta, gamma, mu
    bounds = [(0, 10), (0, 1), (0, 0.5)]  # Rangos para los parámetros
    res = minimize(error, params0, bounds=bounds)
    beta, gamma, mu = res.x
    
    # Calcular R0 ajustado
    R0_ajustado = beta / (gamma + mu)
    
    # Proyección a 10 semanas
    t_proy = np.linspace(0, 9, 100)  # Semanas 1-10
    sol = odeint(modelo_SIRD, y0, t_proy, args=(N, beta, gamma, mu))
    S, I, R, D = sol.T
    
    # Resultados
    resultados_ajustados = {
        "R0_ajustado": R0_ajustado,
        "beta": beta,
        "gamma": gamma,
        "mu": mu,
        "I10": int(I[-1]),
        "D10": int(D[-1]),
        "S10": int(S[-1])
    }
    
    # Gráfico comparativo
    plt.figure(figsize=(12, 6))
    
    # Datos observados
    semanas_obs = [1, 2, 3]
    plt.scatter(semanas_obs, data[ciudad]["I"], color='red', label='Infectados observados')
    plt.scatter(semanas_obs, data[ciudad]["D"], color='black', label='Muertos observados')
    
    # Modelo ajustado
    semanas_modelo = t_proy + 1
    plt.plot(semanas_modelo, I, 'r-', label='Infectados modelados')
    plt.plot(semanas_modelo, D, 'k-', label='Muertos modelados')
    plt.plot(semanas_modelo, S, 'b-', label='Susceptibles')
    
    plt.title(f"Modelo SIRD ajustado para {ciudad}\n$R_0$={R0_ajustado:.2f}, β={beta:.3f}, γ={gamma:.3f}, μ={mu:.3f}")
    plt.xlabel("Semanas")
    plt.ylabel("Número de personas")
    plt.legend()
    plt.grid()
    plt.savefig(f"modelo_ajustado_{ciudad}.png", dpi=300)
    plt.close()
    
    return resultados_ajustados

# Aplicar modelo ajustado a todas las ciudades
resultados_ajustados = {}
for ciudad in data.keys():
    resultados_ajustados[ciudad] = ajustar_y_proyectar(ciudad)
    print("\n" + ciudad + " - Modelo ajustado:")
    print(f"R0 ajustado: {resultados_ajustados[ciudad]['R0_ajustado']:.2f}")
    print(f"Semana 10 - S: {resultados_ajustados[ciudad]['S10']}, I: {resultados_ajustados[ciudad]['I10']}, D: {resultados_ajustados[ciudad]['D10']}")