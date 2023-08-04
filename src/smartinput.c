#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cjson/cJSON.h>
#include <openssl/md5.h>

#include <openmmpol.h>

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

bool check_file(cJSON *file_json, char **path){
    // It should be a structure
    if(!cJSON_IsObject(file_json)){
        return false;
    }
    
    cJSON *file_data = file_json->child;
    char *md5sum = NULL;
    char errstring[OMMP_STR_CHAR_MAX];
    *path = NULL;


    while(file_data != NULL){
        if(strcmp(file_data->string, "path") == 0){
            if(cJSON_IsString(file_data))
                *path = file_data->valuestring;
            else{
                ommp_message("File's path should be a string.", OMMP_VERBOSE_LOW, "SI file");
                return false;
            }
        }
        else if(strcmp(file_data->string, "md5sum") == 0){
            if(cJSON_IsString(file_data))
                md5sum = file_data->valuestring;
            else{
                ommp_message("File's md5sum should be a string.", OMMP_VERBOSE_LOW, "SI file");
                return false;
            }
        }
        else{
            sprintf(errstring, "Unrecognized file attribute %s", file_data->string);
            ommp_message(errstring, OMMP_VERBOSE_LOW, "SIfile");
        }
        file_data = file_data->next;
    }
    // It should contain a path
    if(*path == NULL){
        ommp_message("File entry does not contain a path.", OMMP_VERBOSE_LOW, "SI file");
        return false;
    }

    // The path should exist
    FILE *fp = fopen(*path, "r");
    if(fp == NULL){
        sprintf(errstring, "File %s does not exist.", *path);
        ommp_message(errstring, OMMP_VERBOSE_LOW, "SI file");
        return false;
    }

    // It could contain an md5sum
    if(md5sum != NULL){
        // The md5 sum should be correct
        bool md5chk = md5_file_check(*path, md5sum);
        if(!md5chk){
            sprintf(errstring, "File %s does not correspond to md5sum %s.", *path, md5sum);
            ommp_message(errstring, OMMP_VERBOSE_LOW, "SI file");
            return false;
        }
    }
    return true;
}

void c_smartinput(const char *json_file, OMMP_SYSTEM_PRT *ommp_sys, OMMP_QM_HELPER_PRT *ommp_qmh){
    char msg[OMMP_STR_CHAR_MAX];
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    // Read the whole file content
    FILE *fp = fopen(json_file, "r");
    if(fp == NULL){
        sprintf(msg, "Unable to open %s", json_file);
        ommp_fatal(msg);
    }
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
        if(err != NULL)
            sprintf(msg, "Error before %s.", err);
        else
            sprintf(msg, "Unexpected error during JSON parsing.");
        ommp_fatal(msg);
    }
    
    cJSON *cur = input_json->child;
    char *path, *xyz_path = NULL, *prm_path = NULL, *hdf5_file = NULL, *mmpol_file = NULL;
    while(cur != NULL){
        sprintf(msg, "Parsing JSON element \"%s\".", cur->string);
        ommp_message(msg, OMMP_VERBOSE_DEBUG, "SI");
        if(str_ends_with(cur->string, "_file")){
            if(! check_file(cur, &path)){
                sprintf(msg, "File check on \"%s\" has failed.", cur->string);
                ommp_fatal(msg);
            };
        }

        if(strcmp(cur->string, "xyz_file") == 0){
            xyz_path = path;
        }
        else if(strcmp(cur->string, "prm_file") == 0){
            prm_path = path;
        }
        else{
            sprintf(msg, "Unrecognized JSON element \"%s\".", cur->string);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
            
        }

        cur = cur->next;
    }

    if(xyz_path != NULL){
        ommp_message("Trying initialization from Tinker .xyz file.", OMMP_VERBOSE_LOW, "SI");
        if(prm_path == NULL)
            ommp_fatal("xyz_file set but prm_file is missing.");
        if(hdf5_file != NULL)
            ommp_fatal("xyz_file set but also hdf5_file is set, this is ambiguous.");
        if(mmpol_file != NULL)
            ommp_fatal("xyz_file set but also mmpol_file is set, this is ambiguous.");

        *ommp_sys = ommp_init_xyz(xyz_path, prm_path);
    }
    
    return;
}

void ommp_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Just an interface function to expose same names and functionalities in C and Fortran
    c_smartinput(json_file, ommp_sys, ommp_qmh);
}
