#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 3 && argc != 2){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <JSON FILE> [<OUTPUT FILE>]\n");
        return 1;
    }
    
    OMMP_SYSTEM_PRT my_system;
    OMMP_QM_HELPER_PRT my_qmh;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    
    if(argc == 3)
        ommp_print_summary_to_file(my_system, argv[2]);
    else
        ommp_print_summary(my_system);
    
    if(my_qmh != NULL) ommp_terminate_qm_helper(my_qmh);
    if(my_system != NULL) ommp_terminate(my_system);
    
    return 0;
}
