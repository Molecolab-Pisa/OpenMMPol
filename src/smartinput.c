#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cjson/cJSON.h>
#include <openssl/md5.h>

#include <openmmpol.h>
#define MAX_STR_LEN 1024

bool str_ends_with(const char *s, const char *e){
    size_t slen = strlen(s);
    size_t elen = strlen(e);
    if(elen > slen) return false;

    if(strcmp(&(s[slen-elen]), e) != 0) return false;
    return true;
}

bool md5_file_check(const char* my_file, const char *md5_sum){
    // This function verifies if the file at the address of my_file has 
    // the md5 sum defined by md5_sum. Return false if the check fail.

    unsigned char c[MD5_DIGEST_LENGTH];
    FILE *fp = fopen(my_file, "rb");
    MD5_CTX mdContext;

    if(fp == NULL) return false;

#define BUF_SIZE 1024
    unsigned char buf[BUF_SIZE];

    MD5_Init(&mdContext);
    for(size_t nrd; (nrd = fread(buf, 1, BUF_SIZE, fp)) > 0;)
        MD5_Update (&mdContext, buf, nrd);
    MD5_Final (c, &mdContext);
#undef BUF_SIZE

    unsigned int l;
    
    for(int i = 0; i < MD5_DIGEST_LENGTH; i++){
        sscanf(&(md5_sum[i*2]), "%2x", &l);
        if(l != c[i]) return false;
    }
    fclose(fp);
    return true;
}

bool check_file(cJSON *file_json, char *path){
    // It should be a structure
    if(!cJSON_IsObject(file_json)){
        return false;
    }
    
    cJSON *file_data = file_json->child;
    char *md5sum = NULL;
    path = NULL;

    while(file_data != NULL){
        if(strcmp(file_data->string, "path") == 0){
            if(cJSON_IsString(file_data))
                path = file_data->valuestring;
            else{
                printf("path should be a string\n");
                return false;
            }
        }
        else if(strcmp(file_data->string, "md5sum") == 0){
            if(cJSON_IsString(file_data))
                md5sum = file_data->valuestring;
            else{
                printf("md5sum should be a string\n");
                return false;
            }
        }
        else{
            printf("Unrecognized\n");
        }
        file_data = file_data->next;
    }
    // It should contain a path
    if(path == NULL){
        printf("File entry does not contain a path\n");
        return false;
    }

    // The path should exist
    FILE *fp = fopen(path, "r");
    if(fp == NULL){
        printf("File doesn't exist!\n");
        return false;
    }

    // It could contain an md5sum
    if(md5sum != NULL){
        // The md5 sum should be correct
        bool md5chk = md5_file_check(path, md5sum);
        if(!md5chk){
            printf("MD5 check failed!\n");
            return false;
        }
    }
    return true;
}

void c_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Read the whole file content
    FILE *fp = fopen(json_file, "r");
    // Get file size ...
    size_t fsize;
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
    
    cJSON *cur = input_json->child;
    char *path;
    while(cur != NULL){
        printf("%s\n", cur->string);
        if(str_ends_with(cur->string, "_file")){
            check_file(cur, path);
        } 
        cur = cur->next;
    }
    
    return;
}

void ommp_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Just an interface function to expose same names and functionalities in C and Fortran
    c_smartinput(json_file, ommp_sys, ommp_qmh);
}