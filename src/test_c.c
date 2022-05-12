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

    return 0;
}
