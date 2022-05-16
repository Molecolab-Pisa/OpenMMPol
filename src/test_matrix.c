#include <stdio.h>
#include <stdlib.h>
#include "openmmpol.h"

int main(int argc, char **argv){
    if( argc != 3 ){
        printf("Synthax:\n   $ test_matrix.out M N\n");
        return 0;
    }
    
    double *myvec;
    int32_t m, n;
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    printf("==== C SIDE ====\n");
    printf("m (row): %d, n(columns): %d\n", m, n);
    
    myvec = (double *) malloc(m * n * sizeof(double));

    for(int i=0; i < m; i++)
        for(int j=0; j < n; j++){
            // printf("%d %d %d\n", i, j, i*n+j);
            myvec[i*n+j] = i*n+j;
        }

    for(int i=0; i < m; i++){
        for(int j=0; j < n; j++)
            printf("%5.2f  ", myvec[i*n+j]);
        printf("\n");
    }
    printf("================\n");

    mb22_print_matrix(myvec, (int32_t) m, (int32_t) n, true);
    return 0;
}            
