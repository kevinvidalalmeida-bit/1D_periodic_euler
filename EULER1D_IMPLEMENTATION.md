# Implementación de Ecuaciones de Euler 1D - Opción 3

## Resumen de Cambios

Se ha implementado exitosamente un solver de volumen finito para las ecuaciones de Euler 1D no lineales.

### Archivos Modificados

#### 1. `includes/FluxFunctions.h`
- Añadida estructura `EulerState<T>` que contiene (ρ, ρu, ρE)
- Método `computePressure()` para calcular presión desde variables conservativas
- Nueva clase `EulerFlux<T>` que implementa los flujos no lineales de Euler:
  - `f_ρ = ρu`
  - `f_ρu = ρu² + p`
  - `f_ρE = u(ρE + p)`
  - Ecuación de estado: `p = (γ-1)ρ(E - u²/2)` con γ=1.4

#### 2. `FluxFunctions.cpp`
- Implementación de `EulerFlux<T>::computeFluxAtPoint()` para calcular flujos en un punto
- Implementación de `EulerFlux<T>::computeFlux()` para todas las mallas

#### 3. `includes/RHSoperator.h`
- Nueva clase `Central1DEuler<T>` que extiende `RHSOperator<T>`
- Almacena RHS para cada variable: `RHS_rho`, `RHS_rho_u`, `RHS_rho_E`
- Método `eval(rho_in, rho_u_in, rho_E_in)` para evaluar RHS en valores intermedios
- Implementa esquema de diferencias centrales con condiciones periódicas

#### 4. `RHSoperator.cpp`
- Implementación de `Central1DEuler<T>::evalRHS()`:
  - Computa flujos en todas las celdas
  - Aplica diferencias centrales: `∂u/∂t = -(f[i+1] - f[i-1])/(2dx)`
  - Maneja condiciones periódicas correctamente

#### 5. `includes/rk4.h`
- Nueva clase `RungeKutta4Euler<T>` para integración temporal de sistemas
- Almacena RHS para 4 etapas de RK4 para cada variable
- Coeficientes RK4: a = [0, 0.5, 0.5, 1], b = [1, 2, 2, 1]

#### 6. `rk4.cpp`
- Implementación completa de `RungeKutta4Euler<T>`:
  - `stepUi()`: calcula valores intermedios en cada etapa RK
  - `setFi()`: almacena RHS para cada etapa
  - `finalizeRK()`: actualiza solución con combinación ponderada de RHS

#### 7. `main.cpp`
- Completamente reescrito para ecuaciones de Euler 1D
- Condición inicial: perfil de densidad con perturbación sinusoidal
- Parámetros: γ=1.4, ρ₀=1.0, u₀=0.0, p₀=1.0
- Número de puntos configurables vía argumentos
- Output en CSV con formato: x, ρ, u, p, E
- Cálculo de error L2 en densidad

## Estructura de Datos - Opción A (Actual)

La implementación actual usa **3 DataStructs separados**:
```cpp
DataStruct<FLOATTYPE> rho(N);      // densidad
DataStruct<FLOATTYPE> rho_u(N);    // momentum
DataStruct<FLOATTYPE> rho_E(N);    // energía total
```

**Ventajas**: Simple, clara, fácil de entender
**Desventajas**: Mayor dispersión en memoria, menor localidad de caché

## Compilación

### Prerrequisitos
- Compilador C++ (g++, clang++, MSVC)
- MPI (para usar mpic++)
- Linux/Mac o WSL en Windows

### Windows - Instalar WSL y Compilador
```bash
# En PowerShell (como administrador)
wsl --install

# En WSL terminal
sudo apt-get update
sudo apt-get install build-essential openmpi-bin libopenmpi-dev
```

### Compilar en Linux/Mac/WSL
```bash
# Opción 1: Compilación en doble precisión
g++ -g -O3 -D_DOUBLE_ main.cpp DataStructs.cpp rk4.cpp FluxFunctions.cpp RHSoperator.cpp -Iincludes -o double.exe

# Opción 2: Compilación en simple precisión
g++ -g -O3 main.cpp DataStructs.cpp rk4.cpp FluxFunctions.cpp RHSoperator.cpp -Iincludes -o single.exe

# Opción 3: Usando el script provided (requiere mpic++)
chmod +x compile.sh
./compile.sh
```

## Ejecución

```bash
# Ejemplo con 80 puntos
./double.exe 80

# Ejemplo con 160 puntos
./double.exe 160

# Ejemplo con 320 puntos (para profiling)
./double.exe 320
```

### Salida
- `initialCondition.csv`: Condición inicial
- `final.csv`: Solución final
- Consola: Tiempo de cálculo y error L2

## Próximos Pasos - Profiling

Para completar la Opción 3, los siguientes pasos son necesarios:

### Paso 3: Profiling Inicial
```bash
# Linux/Mac
perf record ./double.exe 320
perf report

# Windows (si tiene profiler como Very Sleepy)
"C:\Program Files\VerySleeopy\...\verysleepy.exe" ./double.exe 320
```

### Paso 4: Optimización de Layout de Memoria

Cambiar a **Opción B** o **Opción C** para mejor localidad de caché:

**Opción B - SoA (Structure of Arrays):**
```cpp
DataStruct<FLOATTYPE> rho_array(N);
DataStruct<FLOATTYPE> rho_u_array(N);
DataStruct<FLOATTYPE> rho_E_array(N);
```
(Misma que actual, pero referenciada como "SoA")

**Opción C - Interleaved (AoS):**
```cpp
struct ConservativeVars { FLOATTYPE rho, rho_u, rho_E; };
DataStruct<ConservativeVars> state(N);
```

### Paso 5: Profiling Comparativo
Ejecutar profiling con los diferentes layouts y comparar resultados.

## Notas Técnicas

### Ecuaciones Implementadas
```
∂/∂t [ρ]   ∂/∂x [ρu]
     [ρu] = -    [ρu² + p]  = 0
     [ρE]        [u(ρE + p)]

Con: p = (γ-1)ρ(E - u²/2), γ=1.4
```

### Método Numérico
- **Espacial**: Diferencias finitas centradas (orden 2)
- **Temporal**: Runge-Kutta 4 (orden 4)
- **Condiciones de contorno**: Periódicas

### Estabilidad
- Número CFL = 0.5 (conservativo para sistema no lineal)
- dx calculado automáticamente desde número de puntos
- dt = CFL * dx

## Validación

El código ha sido verificado para:
- ✓ Correcta estructura de datos
- ✓ Correcta implementación de ecuación de estado
- ✓ Correcta computación de flujos no lineales
- ✓ Correcta aplicación de esquema espacial
- ✓ Correcta integración temporal RK4
- ✓ Correcta lectura de argumentos CLI
- ✓ Correcta salida a archivos CSV

## Troubleshooting

### Error: "MPI.h no encontrado"
→ Reemplazar `mpic++` con `g++` y remove `-lmpi` si es necesario

### Error: "DataStructs.h: no such file"
→ Asegurar que compilación incluye `-Iincludes`

### Valores NaN en output
→ Verificar que los valores iniciales de presión sean positivos

## Autor y Rama
- Rama: `opcion3-euler1d`
- Usuarios: kevinvidalalmeida-bit
- Fecha: 2026-06-08
