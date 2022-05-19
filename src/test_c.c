#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if( argc != 2){
        printf("Syntax expected\n");
        printf("    $ test_init.exe input_file.mmp\n");
        return 0;
    }
    
    char infile_mmp[120];
    strcpy(infile_mmp, argv[1]);
        
    printf("MMP file: '%s'\n", infile_mmp);

    w_mmpol_init(infile_mmp);

    write_hdf5("test.hdf5");

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
