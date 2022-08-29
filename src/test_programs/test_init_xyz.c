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
    ommp_init_xyz(argv[1], argv[2]);
    
    if(argc == 4)
        ommp_print_summary_to_file(argv[3]);
    else
        ommp_print_summary();
    
    ommp_terminate();
    
    return 0;
}
