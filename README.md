# Producto de matrices con MPI: distribución unidimensional

Implementación paralela en C del producto de matrices **C = α · A · B** mediante paso de mensajes con **MPI**, repartiendo el trabajo entre procesos con una distribución por bloques de filas.

Desarrollado como práctica de la asignatura *Arquitectura de Computadores* (Grado en Ingeniería Informática, Universidade da Coruña, curso 2025/2026).

## Planteamiento

El proceso 0 lee las dimensiones (*m*, *n*, *k*) y el factor de escala α, y genera las matrices completas. A partir de ahí:

1. **Difusión de parámetros** a todos los procesos con `MPI_Bcast`.
2. **Replicación de B** en todos los procesos, ya que cada uno la necesita entera para calcular sus filas de C.
3. **Distribución de A por bloques de filas** con `MPI_Scatterv`.
4. **Cómputo local**: cada proceso multiplica su bloque de A por B.
5. **Recolección de C** en el proceso 0 con `MPI_Gatherv`.

Al terminar, el proceso 0 repite el cálculo de forma secuencial y compara ambos resultados elemento a elemento.

## Aspectos destacados

- **Funciona con cualquier tamaño de matriz y cualquier número de procesos.** Cuando el número de filas no es divisible entre los procesos, el resto se reparte de una en una entre los primeros, de forma que la carga queda equilibrada con una diferencia máxima de una fila. Funciona incluso cuando hay más procesos que filas.
- **Todas las comunicaciones son operaciones colectivas**, sin envíos punto a punto. Esto simplifica el código y deja la optimización de las comunicaciones en manos de la implementación de MPI.
- Cálculo explícito de `sendcounts` y `displs` para las versiones vectoriales de las colectivas.
- Las matrices completas solo se reservan en el proceso raíz.

## Validación

Verificado frente a la versión secuencial con 1, 3, 7 y 12 procesos para una matriz de 10 filas, lo que incluye repartos no exactos y procesos sin filas asignadas.

## Compilación y ejecución

```bash
make                                    # compila
make run NP=7 ARGS="10 7 5 2.0"         # ejecuta con 7 procesos: m n k alfa
```

O directamente:

```bash
mpicc -Wall -O2 producto_matrices_1d.c -o producto_matrices_1d
mpirun -np 4 ./producto_matrices_1d <m> <n> <k> <alfa>
```

## Archivos

```
producto_matrices_1d.c   Implementación
Makefile
docs/enunciado.pdf       Enunciado de la práctica
```

Ver también: [SUMMA](https://github.com/Maveriick123/summa-mpi), la versión con distribución bidimensional que elimina la necesidad de replicar B.
