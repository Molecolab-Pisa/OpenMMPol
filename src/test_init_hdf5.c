#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 2 && argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <HDF5 FILE> [<OUTPUT FILE>]\n");
        return 0;
    }
    
    set_verbose(OMMP_VERBOSE_DEBUG);
    mmpol_init_hdf5(argv[1]);
    
    ommp_terminate();
    
    return 0;
}

