#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if( argc < 2 || argc > 3){
        printf("Syntax expected\n");
        printf("    $ test_init.exe input_file.mmp [ scratch_file.rwf ]\n");
        return 0;
    }
    
    char infile_mmp[120], infile_rwf[120];
    strcpy(infile_mmp, argv[1]);
        
    if( argc == 3 )
        strcpy(infile_rwf, argv[2]);
    else{
        strcpy(infile_rwf, infile_mmp);
        int len = strlen(infile_rwf) - 3;
        sprintf(&(infile_rwf[len]), "rwf");
    }
        
    printf("MMP file: '%s'\n", infile_mmp);
    printf("RWF file: '%s'\n", infile_rwf);

    w_mmpol_init(infile_mmp);
    
    int mm_atoms = get_mm_atoms();
    double *_cmm = (double *) get_cmm();
    double **cmm = (double **) malloc(sizeof(double *) * mm_atoms);
    for(int i=0; i < mm_atoms; i++){
        cmm[i] = &(_cmm[3*i]);
    }

    for(int i=0; i<mm_atoms; i++){
        for(int j=0; j<3; j++)
            printf("%5.2f ", cmm[i][j]);
        printf("\n");
    }

    return 0;
}
