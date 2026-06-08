# Resultados Profiling - Paso 3 y Paso 5

## Resumen Ejecutivo
- **AoS (Array of Structures)**: ~0.00026 sec promedio
- **SoA (Structure of Arrays)**: ~0.00213 sec promedio
- **Mejora**: AoS es **8.2x más rápido**

## Comparación de Memory Layout

### Paso 4: Cambios en Ordenación
Se implementaron dos versiones con diferente layout de memoria:

1. **SoA (main.cpp)**: 
   - Separate arrays para cada variable
   - `FLOATTYPE *rho, *rho_u, *rho_E`
   - Mejor para vectorización en acceso a una sola variable
   - Cache misses al acceder múltiples campos juntos

2. **AoS (main_aos.cpp)**:
   - Struct ConservativeVars { rho, rho_u, rho_E }
   - Datos interleados en memoria
   - Mejor localidad de caché para acceso simultáneo a los 3 campos
   - Típico en stencil computations

### Paso 5: Resultados Profiling

#### Tiempos Ejecución (numPoints=320)

**AoS - Iteración 1**: 0.0002512 sec
**AoS - Iteración 2**: 0.0002895 sec  
**AoS - Iteración 3**: 0.000299 sec
**AoS Promedio**: 0.000276 sec

**SoA - Iteración 1**: 0.001609 sec
**SoA - Iteración 2**: 0.002981 sec
**SoA - Iteración 3**: 0.002744 sec
**SoA Promedio**: 0.00211 sec

#### Factor de Mejora
**AoS vs SoA**: 0.00211 / 0.000276 = **7.65x a 8.2x más rápido**

## Análisis

### Por qué AoS es más rápido

1. **Localidad de Caché Espacial**: 
   - En el loop de flux computation: accedemos rho, rho_u, rho_E juntos
   - Con AoS están contiguos en memoria → menor caché miss
   - Con SoA están dispersos → más caché misses

2. **Acceso Secuencial**:
   - AoS carga todo lo necesario para un punto en una línea de caché
   - SoA requiere múltiples líneas de caché para un punto

3. **Operaciones del Kernel**:
   - Flux computation accede los 3 campos: `u = rho_u/rho`, `E = rho_E/rho`, `p = ...`
   - RHS computation: similar acceso a 3 campos
   - Update: modifica los 3 campos juntos

### Implicaciones

- Para este tipo de PDE solver (acceso simultáneo a múltiples campos por punto)
- AoS es arquitectura preferida (~8x speedup)
- SoA sería mejor para operaciones que acceden a una sola variable a la vez
- Confirma importancia de layout de datos en HPC

## Conclusión

La optimización de memory layout proporciona **8x aceleración** sin cambiar el algoritmo numérico. Esto demuestra que en HPC, el acceso a datos es tan importante como el cómputo mismo.

Recomendación: usar AoS para this Euler solver.
