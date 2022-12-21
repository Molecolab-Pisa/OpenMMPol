#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 3 && argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <XYZ FILE> <PRM FILE> [<OUTPUT FILE>]\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[2]);
    
    if(argc == 4)
        ommp_print_summary_to_file(my_system, argv[3]);
    else
        ommp_print_summary(my_system);
    
    ommp_terminate(my_system);
    
    return 0;
}
