#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <INPUt XYZ FILE> <INPUT PRM FILE> <OUTPUT MMP FILE>\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    ommp_init_xyz(argv[1], argv[2]);
    
    ommp_save_mmp(argv[3], 3);
    
    ommp_terminate();
    
    return 0;
}
