#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 4){
        printf("Syntax expected\n");
        printf("    $ xyz2mmp <INPUt XYZ FILE> <INPUT PRM FILE> <OUTPUT MMP FILE>\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[2]);
    
    ommp_save_mmp(my_system, argv[3], 3);
    
    ommp_terminate(my_system);
    
    return 0;
}
