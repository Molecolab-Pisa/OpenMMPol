#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 3){
        printf("Syntax expected\n");
        printf("    $ ommp_pp <JSON SI FILE> <OMMP HDF5 OUTPUT FILE>\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_LOW);
    OMMP_SYSTEM_PRT my_system;
    OMMP_QM_HELPER_PRT my_qmh;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    if(my_qmh != NULL){
        printf("Only MM part of the system could be saved in HDF5 format, please remove all the section referring to QM.\n");
        return 1;
    }

    ommp_save_as_hdf5(my_system, argv[2], "system");
    if(my_system != NULL) ommp_terminate(my_system);
    
    return 0;
}
