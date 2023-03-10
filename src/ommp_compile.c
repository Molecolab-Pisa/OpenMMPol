#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <XYZ FILE> <PRM FILE> <OMMP HDF5 INPUT FILE>\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_LOW);
    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[2]);

    ommp_save_as_hdf5(my_system, argv[3], "system");
    
    ommp_terminate(my_system);
    
    return 0;
}
