#include <stdio.h>
#include <stdlib.h>
#include <cjson/cJSON.h>

#include <openmmpol.h>
#define MAX_STR_LEN 1024

void c_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Read the whole file content
    FILE *fp = fopen(json_file, "r");
    // Get file size ...
    long int fsize;
    for(fsize = 0; getc(fp) != EOF; fsize++);
    rewind(fp);

    char *file_content = malloc(fsize * sizeof(char));
    fread(file_content, sizeof(char), fsize, fp);
    
    // Parse the input json
    cJSON *input_json = cJSON_Parse(file_content);
    free(file_content);
    if(input_json == NULL){
        const char *err = cJSON_GetErrorPtr();
        char error_string[MAX_STR_LEN];
        if(err != NULL)
            sprintf(error_string, "Error before %s\n", err);
        else
            sprintf(error_string, "Unexpected error during JSON parsing\n");
        //ommp_fatal(...);
        printf("%s", error_string);
    }
    
    char *string = NULL;
    string = cJSON_Print(input_json);
    if (string == NULL)
        printf("Failed to print.\n");
    else
        printf("%s", string);

    return;
}

void ommp_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Just an interface function to expose same names and functionalities in C and Fortran
    c_smartinput(json_file, ommp_sys, ommp_qmh);
}