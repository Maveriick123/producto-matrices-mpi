/*
    * Implementación SPDM
    * m, n, k y alfa pasados por teclado
    * Resultados correctos para cualquier m, n y k
    * Se hace el producto de dos matrices de tal forma que: alfa * A(m*k) * B(k*n) = C(m*n)
    * Inicialización de A y B en proceso 0
    * Distribución y replicación de datos
    * Recolección de matriz C al final en proceso 0
*/

/*
    * mpicc -Wall producto_matrices_1d.c -o producto_matrices_1d
    * mpirun -np 7 -oversubscribe ./producto_matrices_1d 6 5 4 2.2
*/

#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

void productoMatrices(float alfa, float *A, float *B, float *res, int m, int n, int k){
    int i, j, z;
    // El orden de los bucles es importante para el correcto funcionamiento 
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            for(z = 0; z < k; z++){
                res[i*n+j] += alfa * A[i*k+z] * B[z*n+j];
            }
        }
    }
    return;
}

void printMatriz(float *M, int filas, int columnas){
    int i,j;
    for(i = 0; i < filas; i++){
        for(j = 0; j < columnas; j++){
            printf("%c%12.6f%c", (j == 0)?'\t': ' ', M[i*columnas + j], (j+1==columnas)?'\n': ' ');
        }
    }
    return;
}


int main(int argc, char *argv[]){
    int nprocs, myid;
    int m, n, k;
    float alfa;
    int i,j;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Proceso 0 lee los parametros de entrada

    if(!myid){ // Vamos a asumir que los parámetros de entrada son correctos.
        if(argc > 4){ // 4 parametros + nombre del ejecutable.
            m = atoi(argv[1]);
            n = atoi(argv[2]);
            k = atoi(argv[3]);
            alfa = atof(argv[4]);
        } else{
            printf("Número de parámetros incorrecto");
            MPI_Finalize();
            return 0;
        }
    }

    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alfa, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //Forma dinámica
    float *B = malloc(n*k*sizeof(float));
    float *A = NULL;
    float *C = NULL;

    // El proceso 0 inicializa a A y B

    if(!myid){
        A = malloc(m*k*sizeof(float));
        C = malloc(m*n*sizeof(float));

        // C es row-major-order
        // Hacemos acceso lineal a memoria contigua
        for(i = 0;i < m;i++){
            for(j = 0;j < k;j++){
                // 1.0f evita conversiones implicitas, más eficiente --> float nativo
                A[i*k+j] = i + j;
            }
        }

        for(i = 0; i < k; i++){
            for(j = 0; j < n; j++){
                B[i*n+j] = 1 + 2*i + j;
            }
        }
    }

    // Número de filas por procesador (distrbución por filas)

    // Ejemplo: 10 filas y 3 procesos
    /*
    * 10/3 = 3 filas de base para todos
    * 10%3 = 1 fila restante
    * En este caso solo el 0 asume una fila de más (...)
    */
    
    int local_m;
    int base = m / nprocs;
    int resto = m % nprocs;
    if(myid < resto) local_m = base + 1;
    else local_m = base;

    // Inicialización A y C locales

    float *localA = malloc(local_m * k * sizeof(float));

    float *localC = malloc(local_m * n * sizeof(float));

    // Para localC inicializamos todo a 0, ya que con malloc tenemos valores indefinidos.
    
    for(i = 0; i < local_m; i++){
        for(j = 0; j < n; j ++){
            localC[i*n+j] = 0.0;
        }
    }

    // Calculo de sendcounts y displs para MPI_Scatterv

    int *sendcounts = NULL;
    int *displs = NULL;

    if(!myid){
        int offset = 0;
        sendcounts = malloc(nprocs * sizeof(int));
        displs = malloc(nprocs * sizeof(int));
        for (i = 0; i < nprocs; i++){
            sendcounts[i] = (base + (i < resto)) * k; // Floats que recibe cada proceso
            displs[i] = offset;
            offset += sendcounts[i]; // Suma acumulada
        }
    }



    MPI_Bcast(B,k*n,MPI_FLOAT,0, MPI_COMM_WORLD);
    // "The send buffer is ignored for all non-root processes"
    MPI_Scatterv(A, sendcounts, displs, MPI_FLOAT, localA, local_m*k, MPI_FLOAT, 0, MPI_COMM_WORLD);

    productoMatrices(alfa, localA, B, localC, local_m, n, k);

    // Calculo de recvcounts y displs para MPI_Gatherv

    int *recvcounts = NULL;
    int *displs_g = NULL;
    if(!myid){
        recvcounts = malloc(nprocs * sizeof(int));
        displs_g = malloc(nprocs * sizeof(int));
        int offset = 0;
        for(i = 0; i < nprocs; i++){
            recvcounts[i] = (base + (i<resto)) * n;
            displs_g[i] = offset;
            offset += recvcounts[i];
        }
    }

    MPI_Gatherv(localC, local_m*n, MPI_FLOAT, C, recvcounts,displs_g, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Mostrar los resultados y comparar con programa secuencial:

    if(myid == 0){
        printf("\nResultado en Paralelo\n");
        printf("\tMatriz A\n");
        printMatriz(A, m, k);
        printf("\tMatriz B\n");
        printMatriz(B, k, n);
        printf("\tMatriz C\n");
        printMatriz(C, m, n);
        printf("\nComprobación en Secuencial\n");
        float *m_sec = calloc(m*n, sizeof(float));
        productoMatrices(alfa, A,B,m_sec,m,n,k);
        int cnt = 0;
        for (i = 0; i < m; i++){
            for(j = 0; j < n; j++){
                if(C[i*n+j]!=m_sec[i*n+j]){
                    printf("\t (*) Error en la fila %d columna %d, ya que %f != %f", i, j, C[i*n+j],m_sec[i*n+j]);
                    cnt++;
                }
            }
        }
        (cnt == 0)? printf("\t (*) Matrices idénticas\n\n"): printf("\t (*) Matrices distintas\n\n");
        free(m_sec);
        free(A);
        free(C);
    }

    free(B);
    free(localA);
    free(localC);
    fflush(stdout);
    MPI_Finalize();
    return 0;
}





