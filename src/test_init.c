#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if( argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init.exe <INPUT FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    char infile[120], outfile[120];
    strcpy(infile, argv[1]);
    strcpy(outfile, argv[2]);
        
    printf("Input file: '%s'\n", infile);

    w_mmpol_init(infile);
    set_verbose(OMMP_VERBOSE_DEBUG);
    print_summary();
    
    return 0;
}
