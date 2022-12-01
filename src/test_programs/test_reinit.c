#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init.exe <MMP FILE 1> <MMP FILE 2>\n");
        return 0;
    }
    
    char infile1[120], infile2[120];
    strcpy(infile1, argv[1]);
    strcpy(infile2, argv[2]);
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    ommp_init_mmp(infile1);
    ommp_terminate();
    ommp_init_mmp(infile2);
    
    ommp_print_summary();
    
    ommp_terminate();
    
    return 0;
}
