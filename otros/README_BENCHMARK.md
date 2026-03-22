# Modelo de Rashevsky - Benchmark de Metodos Numericos

Implementacion de tres metodos numericos (Euler, Taylor orden 2, Trapecio) en Python, C++ y Fortran con benchmark de 10^7 pasos.

## Uso

### Compilar y ejecutar todo

```bash
make -f Makefile.benchmark run
```

### Ver comparacion

```bash
make -f Makefile.benchmark compare
```

### Graficar resultados

```bash
python3 comparar.py
```

### Limpiar

```bash
make -f Makefile.benchmark clean
```

## Archivos

### Codigo fuente
- `benchmark.py` - Python
- `benchmark.cpp` - C++
- `benchmark.f90` - Fortran

### Resultados
- `benchmark_python.dat`
- `benchmark_cpp.dat`
- `benchmark_fortran.dat`
- `comparacion.png`

## Parametros

- b = 0.02 (tasa nacimiento)
- d = 0.015 (tasa mortalidad)
- r = 0.1 (proporcion fija)
- p0 = 0.01 (condicion inicial)
- t_max = 50 (tiempo maximo)
- n_steps = 10^7 (pasos de integracion)

## Metodos

1. Euler: p(n+1) = p(n) + h*f(p(n))
2. Taylor 2: p(n+1) = p(n) + h*f(p(n)) + (h^2/2)*f'(p(n))
3. Trapecio: [p(n) + (h/2)*f(p(n)) + (h*r*b/2)] / [1 + (h*r*b/2)]

donde f(p) = r*b*(1-p)
