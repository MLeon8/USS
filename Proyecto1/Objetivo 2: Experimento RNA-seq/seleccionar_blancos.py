
#!/usr/bin/env python3
import pandas as pd

# Leer archivo de anotaciones estructurales
estructura_df = pd.read_csv("estructuras/resumen_estructuras.csv")
print("Estructuras de blancos:")
print(estructura_df)

# Leer predicción de sitios de regulación
sitios_df = pd.read_csv("sitios_regulacion/prediccion_sitios_blancos.csv")
print("Sitios regulatorios detectados:")
print(sitios_df)

# Leer resultados de docking
docking_df = pd.read_csv("drogas/docking_ligandos_vs_blancos.csv")
print("Resultados de docking:")
print(docking_df)
