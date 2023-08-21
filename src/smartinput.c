#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cjson/cJSON.h>
#include <openssl/md5.h>

#include <openmmpol.h>

typedef struct _semversion{
    int major;
    int minor;
    int patch;
    int ncommit;
    char commit[8];
    bool clean;
} semversion;

bool str_ends_with(const char *s, const char *e){
    size_t slen = strlen(s);
    size_t elen = strlen(e);
    if(elen > slen) return false;

    if(strcmp(&(s[slen-elen]), e) != 0) return false;
    return true;
}

semversion str_to_semversion(char *strin){
    semversion v;
    v.major = -1;
    v.minor = -1;
    v.patch = -1;
    v.ncommit = -1;
    strcpy(v.commit, "        ");
    v.clean = true;

    semversion v_err = v;

    char msg[OMMP_STR_CHAR_MAX];
    char str[256];
    strcpy(str, strin);

    char *major = strtok(str, ".");
    bool major_ok = (major != NULL);
    int i;

    sprintf(msg, "Major \"%s\".", major);
    ommp_message(msg, OMMP_VERBOSE_DEBUG, "SemVersDB");

    for(i = 0; major[i] != '\0' && major_ok; i++)
        major_ok &= isdigit(major[i]);
    major_ok &= (i>0);

    if(!major_ok){
        sprintf(msg, "Malformed major version in \"%s\".", strin);
        ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
        return v_err;
    }

    v.major = atoi(major);
    
    char *minor = strtok(NULL, ".");
    bool minor_ok = (minor != NULL);
    
    sprintf(msg, "Minor \"%s\".", minor);
    ommp_message(msg, OMMP_VERBOSE_DEBUG, "SemVersDB");

    for(i = 0; minor[i] != '\0' && minor_ok; i++)
        minor_ok &= isdigit(minor[i]);
    minor_ok &= (i>0);

    if(!minor_ok){
        sprintf(msg, "Malformed minor version in \"%s\".", strin);
        ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
        return v_err;
    }

    v.minor = atoi(minor);

    char *patchplus = strtok(NULL, "");
    char *patch = strtok(patchplus, "+");
    char *plus = strtok(NULL, "+");

    bool patch_ok = (patch != NULL);

    sprintf(msg, "Patch \"%s\".", patch);
    ommp_message(msg, OMMP_VERBOSE_DEBUG, "SemVersDB");

    for(i = 0; patch[i] != '\0' && patch_ok; i++)
        patch_ok &= isdigit(patch[i]);
    patch_ok &= (i>0);

    if(!patch_ok){
        sprintf(msg, "Malformed patch version in \"%s\".", strin);
        ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
        return v_err;
    }

    v.patch = atoi(patch);

    if(plus != NULL){
        if(strtok(NULL, "+") != NULL){
            sprintf(msg, "Malformed pre-release string in \"%s\".", strin);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
            return v_err;
        }

        sprintf(msg, "Pre-Release \"%s\".", plus);
        ommp_message(msg, OMMP_VERBOSE_DEBUG, "SemVersDB");

        char *ncommits = strtok(plus, ".");
        bool ncommits_ok = (ncommits != NULL);

        if(ncommits[0] != 'r') ncommits_ok = false;

        for(i = 1; ncommits[i] != '\0' && ncommits_ok; i++)
            ncommits_ok &= isdigit(ncommits[i]);
        ncommits_ok &= (i>0);
    
        if(!ncommits_ok){
            sprintf(msg, "Malformed pre-release string in \"%s\".", strin);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
            return v_err;
        }

        v.ncommit = atoi(&(ncommits[1]));
        
        char *commithash = strtok(NULL, ".");
        bool commithash_ok = (commithash != NULL);

        for(i = 0; commithash[i] != '\0' && commithash_ok; i++);
        commithash_ok &= (i == 8);
    
        if(!commithash_ok){
            sprintf(msg, "Malformed pre-release string in \"%s\".", strin);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
            return v_err;
        }

        strcpy(v.commit, commithash);

        char *clean = strtok(NULL, ".");
        if(clean != NULL){
            if(strcmp(clean, "dirty") == 0) 
                v.clean = false;
            else{
                sprintf(msg, "Malformed pre-release string in \"%s\".", strin);
                ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
                return v_err;
            }
        }
    }
    else{
        v.ncommit = 0;
    }
    return v;
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

bool check_version(char *verstr){
    semversion vreq = str_to_semversion(verstr);
    
    char *ommp_vstr[256];
    sprintf(ommp_vstr, "%s", OMMP_VERSION_STRING);
    semversion vommp = str_to_semversion(ommp_vstr);
    
    if(vreq.major < vommp.major) return true;
    if(vreq.major > vommp.major) return false;
    
    if(vreq.minor < vommp.minor) return true;
    if(vreq.minor > vommp.minor) return false;
    
    if(vreq.patch < vommp.patch) return true;
    if(vreq.patch > vommp.patch) return false;
    if(vreq.ncommit == 0 && vommp.ncommit > 0) return true;
    if(vreq.ncommit > 0 && vommp.ncommit > 0){
        if(vreq.ncommit == vommp.ncommit && 
           strcmp(vreq.commit, vommp.commit) == 0)
            return true;
        else
            return false;
    }

    return true;
}

void c_smartinput(const char *json_file, OMMP_SYSTEM_PRT *ommp_sys, OMMP_QM_HELPER_PRT *ommp_qmh){
    char msg[OMMP_STR_CHAR_MAX];
    
    // TODO Just for debugging
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);

    // Read the whole file content
    sprintf(msg, "Parsing JSON file \"%s\".", json_file);
    ommp_message(msg, OMMP_VERBOSE_LOW, "SI");

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
    char *json_name, *json_description;
    int32_t req_verbosity = OMMP_VERBOSE_DEFAULT;

    while(cur != NULL){
        sprintf(msg, "Parsing JSON element \"%s\".", cur->string);
        ommp_message(msg, OMMP_VERBOSE_DEBUG, "SIDB");
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
        else if(strcmp(cur->string, "version") == 0){
            if(!check_version(cur->valuestring)){
                sprintf(msg, "Required version (%s) is not compatible with OMMP version (%s).", cur->valuestring, OMMP_VERSION_STRING);
                ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
                ommp_fatal("Cannot complete Smart Input initialization.");
            }
                sprintf(msg, "Required version (%s) is compatible with OMMP version (%s).", cur->valuestring, OMMP_VERSION_STRING);
                ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
        }
        else if(strcmp(cur->string, "verbosity") == 0){
            if(strcmp(cur->valuestring, "none") == 0)
                req_verbosity = OMMP_VERBOSE_NONE;
            else if(strcmp(cur->valuestring, "low") == 0)
                req_verbosity = OMMP_VERBOSE_LOW;
            else if(strcmp(cur->valuestring, "high") == 0)
                req_verbosity = OMMP_VERBOSE_HIGH;
            else if(strcmp(cur->valuestring, "debug") == 0)
                req_verbosity = OMMP_VERBOSE_DEBUG;
        }
        else if(strcmp(cur->string, "name") == 0){
            json_name = cur->valuestring;
        }
        else if(strcmp(cur->string, "description") == 0){
            json_description = cur->valuestring;
        }
        else{
            sprintf(msg, "Unrecognized JSON element \"%s\".", cur->string);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
        }

        cur = cur->next;
    }

    ommp_set_verbose(req_verbosity);
    
    sprintf(msg, "Smart Input Name: %s", json_name);
    ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
    
    sprintf(msg, "Smart Input Description: %s", json_description);
    ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
    
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
    else if(mmpol_file != NULL){
        if(prm_path == NULL)
            ommp_fatal("prm_file set but it is not needed for mmpol input, this is ambiguous.");
        if(hdf5_file != NULL)
            ommp_fatal("mmpol_file set but also hdf5_file is set, this is ambiguous.");

        *ommp_sys = ommp_init_mmp(mmpol_file);
    }
    else if(hdf5_file != NULL){
        if(prm_path == NULL)
            ommp_fatal("prm_file set but it is not needed for mmpol input, this is ambiguous.");

        ommp_fatal("input from hdf5_file is still not implemented.");
    }
    else{
        ommp_fatal("No input for MM system found in Smart Input file, set one of xyz_file+prm_file, mmpol_file, hdf5_file");
    }
    
    return;
}

void ommp_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Just an interface function to expose same names and functionalities in C and Fortran
    c_smartinput(json_file, ommp_sys, ommp_qmh);
}
