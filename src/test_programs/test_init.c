#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 3 && argc != 2){
        printf("Syntax expected\n");
        printf("    $ test_init.exe <INPUT FILE> [<OUTPUT FILE>]\n");
        return 0;
    }
    
    char infile[120], outfile[120];
    strcpy(infile, argv[1]);
    
    if(argc == 3)
        strcpy(outfile, argv[2]);
        
    OMMP_SYSTEM_PRT my_system = ommp_init_mmp(infile);
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    if(argc == 3)
        ommp_print_summary_to_file(my_system, outfile);
    else
        ommp_print_summary(my_system);
    
    ommp_terminate(my_system);
    
    return 0;
}
