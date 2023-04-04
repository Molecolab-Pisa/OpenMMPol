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
    
    OMMP_SYSTEM_PRT my_system = ommp_init_hdf5(argv[1], "system");
    
    if(argc == 3)
        ommp_print_summary_to_file(my_system, argv[2]);
    else
        ommp_print_summary(my_system);
    
    ommp_terminate(my_system);
    
    return 0;
}

