#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cjson/cJSON.h>
#include <openssl/md5.h>

#include <openmmpol.h>

#define NELEM 11
char ELEMENTS[NELEM][2] = {"X", 
                           "H", "He", 
                           "Li", "Be", "B", "C", "N", "O", "F", "Ne"};

int element_to_Z(char *elem){
    for(int i=0; i<NELEM; i++)
        if(strcmp(elem, ELEMENTS[i]) == 0) return i;
    return -1;
}

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
        major_ok = (major_ok && isdigit(major[i]));
    major_ok = (major_ok && i>0);

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
        minor_ok = (minor_ok && isdigit(minor[i]));
    minor_ok = (minor_ok && i>0);

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
        patch_ok = (patch_ok && isdigit(patch[i]));
    patch_ok = (patch_ok && i>0);

    if(!patch_ok){
        sprintf(msg, "Malformed patch version in \"%s\".", strin);
        ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
        return v_err;
    }

    v.patch = atoi(patch);

    if(plus != NULL){
        if(strtok(NULL, "+") != NULL || plus[0] != 'r'){
            sprintf(msg, "Malformed pre-release string in \"%s\".", strin);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
            return v_err;
        }

        sprintf(msg, "Pre-Release \"%s\".", plus);
        ommp_message(msg, OMMP_VERBOSE_DEBUG, "SemVersDB");

        plus = &(plus[1]);
        char *ncommits = strtok(plus, ".");
        bool ncommits_ok = (ncommits != NULL);
        if(ncommits_ok && strcmp(ncommits, "dirty") == 0){
            v.clean = false;
            v.ncommit = 0;
        }
        else{
            for(i = 0; ncommits[i] != '\0' && ncommits_ok; i++)
                ncommits_ok = (ncommits_ok && isdigit(ncommits[i]));
            ncommits_ok = (ncommits_ok && i>0);
        
            if(!ncommits_ok){
                sprintf(msg, "Malformed pre-release string in \"%s\".", strin);
                ommp_message(msg, OMMP_VERBOSE_LOW, "SemVers");
                return v_err;
            }

            v.ncommit = atoi(&(ncommits[0]));
            
            char *commithash = strtok(NULL, ".");
            bool commithash_ok = (commithash != NULL);

            for(i = 0; commithash[i] != '\0' && commithash_ok; i++);
            commithash_ok = (commithash_ok && i == 8);
        
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

bool check_file(cJSON *file_json, char **path, char *outmode){
    // It should be a structure
    if(!cJSON_IsObject(file_json)){
        ommp_message("Unexpected errror: qm_json is not a cJSON object.", OMMP_VERBOSE_LOW, "SI file");
        return false;
    }
    
    cJSON *file_data = file_json->child;
    char *md5sum = NULL;
    char errstring[OMMP_STR_CHAR_MAX];
    char accessmode = 'r';
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
        else if(strcmp(file_data->string, "mode") == 0){
            // Default mode is read, for output either append or write or overwrite should be specified.
            if(strcmp(file_data->valuestring, "read") == 0){
                accessmode = 'r';
            }
            else if(strcmp(file_data->valuestring, "write") == 0){
                accessmode = 'w';
            }
            else{
                sprintf(errstring, "Unrecognized file access mode %s; assuming read.", file_data->valuestring);
                ommp_message(errstring, OMMP_VERBOSE_LOW, "SIfile");
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

    FILE *fp = fopen(*path, "r");
    if(fp == NULL && accessmode == 'r'){
        sprintf(errstring, "File %s in read mode does not exist.", *path);
        ommp_message(errstring, OMMP_VERBOSE_LOW, "SI file");
        return false;
    }
    
    if(fp != NULL && accessmode == 'w'){
        sprintf(errstring, "File %s in write mode already exists, cannot overwrite.", *path);
        ommp_message(errstring, OMMP_VERBOSE_LOW, "SI file");
        return false;
    }
    if(fp != NULL) fclose(fp);
    *outmode = accessmode;

    // It could contain an md5sum
    if(md5sum != NULL){
        if(accessmode == 'w'){
            sprintf(errstring, "md5sum check does not make any sense for output files; skipping.");
            ommp_message(errstring, OMMP_VERBOSE_LOW, "SI file");
        }
        else{
            // The md5 sum should be correct
            bool md5chk = md5_file_check(*path, md5sum);
            if(!md5chk){
                sprintf(errstring, "File %s does not correspond to md5sum %s.", *path, md5sum);
                ommp_message(errstring, OMMP_VERBOSE_LOW, "SI file");
                return false;
            }
        }
    }
    return true;
}

bool check_version(char *verstr){
    semversion vreq = str_to_semversion(verstr);
    if(vreq.major < 0) return false;
    
    char *ommp_vstr[256];
    sprintf(ommp_vstr, "%s", OMMP_VERSION_STRING);
    semversion vommp = str_to_semversion(ommp_vstr);
    if(vommp.major < 0) return false;
    
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

bool smartinput_qm(cJSON *qm_json, OMMP_QM_HELPER_PRT *qmh){
    // It should be a structure
    if(!cJSON_IsObject(qm_json)){
        ommp_message("Unexpected errror: qm_json is not a cJSON object.", OMMP_VERBOSE_LOW, "SI QMH");
        return false;
    }
    
    cJSON *qm_data = qm_json->child;
    char errstring[OMMP_STR_CHAR_MAX];
    int32_t natoms = 0;
    int32_t *qmz = NULL, *qmt = NULL;
    double *coords = NULL, *nucq = NULL;
    char *xyz_path = NULL, *prm_path = NULL;
    char mode, *path;

    while(qm_data != NULL){
        if(str_ends_with(qm_data->string, "_file")){
            if(! check_file(qm_data, &path, &mode)){
                sprintf(errstring, "File check on \"%s\" has failed.", qm_data->string);
                ommp_fatal(errstring);
            };
            
            if(mode != 'r' && mode != 'w'){
                sprintf(errstring, "Unrecognized file access mode (unexpected error) '%c'.", mode);
                ommp_fatal(errstring);
            }
        }

        if(strcmp(qm_data->string, "xyz_file") == 0){
            if(mode != 'r'){
                sprintf(errstring, "xyz_file should be in read mode.");
                ommp_fatal(errstring);
            }
            xyz_path = path;
        }
        else if(strcmp(qm_data->string, "prm_file") == 0){
            if(mode != 'r'){
                sprintf(errstring, "prm_file should be in read mode.");
                ommp_fatal(errstring);
            }
            prm_path = path;
        }
        else if(strcmp(qm_data->string, "qm_atoms") == 0){
            if(cJSON_IsArray(qm_data)){
                cJSON *qmat_array = qm_data->child;
                int nat;

                // Check number of elements in array
                for(nat = 0; qmat_array != NULL; qmat_array = qmat_array->next) nat++;
                if(natoms == 0){
                    if(nat == 0){
                        ommp_message("qm_atoms length is zero, that is wired!", OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    natoms = nat;
                }
                else if(natoms != nat){
                    ommp_message("qm_atoms length does not match the length of other arrays in qm",  
                                 OMMP_VERBOSE_LOW, "SI QMH");
                    return false;
                }
                qmz = (int32_t *) malloc(sizeof(int32_t) * natoms);
                nucq = (double *) malloc(sizeof(double) * natoms);

                qmat_array = qm_data->child;
                for(int i=0; qmat_array != NULL; i++){
                    if(!cJSON_IsString(qmat_array)){
                        ommp_message("qm_atoms should be an array of strings.", OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    qmz[i] = element_to_Z(qmat_array->valuestring);
                    if(qmz[i] <= 0){
                        sprintf(errstring, "%s is not a recognized element.", qmat_array->valuestring);
                        ommp_message(errstring, OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    nucq[i] = (double) qmz[i];
                    qmat_array = qmat_array->next;
                }

            }
            else{
                ommp_message("qm_atoms should be an array of strings.", OMMP_VERBOSE_LOW, "SI QMH");
                return false;
            }
        }
        else if(strcmp(qm_data->string, "qm_coords") == 0){
            if(cJSON_IsArray(qm_data)){
                cJSON *qc_array = qm_data->child;
                int nat;

                // Check number of elements in array
                for(nat = 0; qc_array != NULL; qc_array = qc_array->next) nat++;
                if(natoms == 0){
                    if(nat == 0){
                        ommp_message("qm_coords length is zero, that is wired!", OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    natoms = nat;
                }
                else if(natoms != nat){
                    ommp_message("qm_coords length does not match the length of other arrays in qm",  
                                 OMMP_VERBOSE_LOW, "SI QMH");
                    return false;
                }
                coords = (double *) malloc(sizeof(double) * natoms * 3);

                qc_array = qm_data->child;
                for(int i=0; qc_array != NULL; i++){
                    if(!cJSON_IsArray(qc_array)){
                        ommp_message("qm_coords should be a 3 x natoms matrix of double.", 
                                     OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    cJSON *rowel = qc_array->child;
                    for(int j=0; j < 3; j++){
                        if(!cJSON_IsNumber(rowel) || (j==2 && rowel->next != NULL)){
                            ommp_message("qm_coords should be a 3 x natoms matrix of double.", 
                                        OMMP_VERBOSE_LOW, "SI QMH");
                            return false;
                        }
                        coords[3*i+j] = rowel->valuedouble * OMMP_ANG2AU;
                        rowel = rowel->next;
                    }
                    qc_array = qc_array->next;
                }
            }
            else{
                ommp_message("qm_coords should be a 3 x natoms matrix of double.", 
                                OMMP_VERBOSE_LOW, "SI QMH");
                return false;
            }
        }
        else if(strcmp(qm_data->string, "qm_atom_types") == 0){
            if(cJSON_IsArray(qm_data)){
                cJSON *qmat_array = qm_data->child;
                int nat;

                // Check number of elements in array
                for(nat = 0; qmat_array != NULL; qmat_array = qmat_array->next) nat++;
                if(natoms == 0){
                    if(nat == 0){
                        ommp_message("qm_atom_types length is zero, that is wired!", OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    natoms = nat;
                }
                else if(natoms != nat){
                    ommp_message("qm_atom_types length does not match the length of other arrays in qm",  
                                 OMMP_VERBOSE_LOW, "SI QMH");
                    return false;
                }
                qmt = (int32_t *) malloc(sizeof(int32_t) * natoms);

                qmat_array = qm_data->child;
                for(int i=0; qmat_array != NULL; i++){
                    if(!cJSON_IsNumber(qmat_array)){
                        ommp_message("qm_atom_types should be an array of integers.", OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    qmt[i] = qmat_array->valueint;
                    if(qmt[i] <= 0){
                        sprintf(errstring, "Atom types should be positive (found %d).", qmat_array->valueint);
                        ommp_message(errstring, OMMP_VERBOSE_LOW, "SI QMH");
                        return false;
                    }
                    qmat_array = qmat_array->next;
                }

            }
            else{
                ommp_message("qm_atom_types should be an array of integers.", OMMP_VERBOSE_LOW, "SI QMH");
                return false;
            }
        }
        else{
            sprintf(errstring, "Unrecognized qm attribute %s", qm_data->string);
            ommp_message(errstring, OMMP_VERBOSE_LOW, "SI QMH");
        }
        qm_data = qm_data->next;
    }
    
    if(qmz != NULL && nucq != NULL){
        if(xyz_path != NULL){
            ommp_message("Both qm_atoms, qm_coords and xyz_file are set, this is ambiguous.", 
                         OMMP_VERBOSE_LOW, "SI QMH");
            return false;
        }
        if(coords == NULL){
            ommp_message("Both qm_atoms and qm_coords should be set in a correct input.",
                         OMMP_VERBOSE_LOW, "SI QMH");
            return false;
        }

        *qmh = ommp_init_qm_helper(natoms, coords, nucq, qmz);
    }
    else if(coords != NULL){
        ommp_message("Both qm_atoms and qm_coords should be set in a correct input.",
                     OMMP_VERBOSE_LOW, "SI QMH");
        return false;
    }
    else if(xyz_path != NULL){
        ommp_message("QM input from xyz is currently unsupported!",
                     OMMP_VERBOSE_LOW, "SI QMH");
        return false;

    }
    else{
        ommp_message("Either qm_atoms and qm_coords or xyz_file should be set in qm section!",
                     OMMP_VERBOSE_LOW, "SI QMH");
        return false;
    }

    if(qmt != NULL){
        ommp_qm_helper_set_attype(*qmh, qmt);
        if(prm_path != NULL){
            ommp_qm_helper_init_vdw_prm(*qmh, prm_path);
        }
        else{
            ommp_message("Since qm_atom_types is present but prm_file is not, atom types are set, while VdW are not.", OMMP_VERBOSE_LOW, "SI QMH");
        }
    }
    else if(prm_path != NULL){
        ommp_message("prm_file is only needed when qm_atom_types is set, your config does not make any sense.", OMMP_VERBOSE_LOW, "SI QMH");
        return false;
    }

    if(qmz != NULL) free(qmz);
    if(nucq != NULL) free(nucq);
    if(coords != NULL) free(coords); 

    return true;
}

void c_smartinput(const char *json_file, OMMP_SYSTEM_PRT *ommp_sys, OMMP_QM_HELPER_PRT *ommp_qmh){
    char msg[OMMP_STR_CHAR_MAX];
    
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
    char *path, *xyz_path = NULL, *prm_path = NULL,
         *hdf5_path = NULL, *mmpol_path = NULL,
         *output_path = NULL, mode;
    char *json_name, *json_description;
    int32_t req_verbosity = OMMP_VERBOSE_DEFAULT,
            req_solver = OMMP_SOLVER_DEFAULT,
            req_matv = OMMP_MATV_DEFAULT;
    
    int32_t nla = 0, *la_mm, *la_qm, *la_la, *la_ner;
    int32_t nfrozen=0, *frozenat;
    double *la_bl;
    *ommp_qmh = NULL;

    while(cur != NULL){
        sprintf(msg, "Parsing JSON element \"%s\".", cur->string);
        ommp_message(msg, OMMP_VERBOSE_DEBUG, "SIDB");
        if(str_ends_with(cur->string, "_file")){
            if(! check_file(cur, &path, &mode)){
                sprintf(msg, "File check on \"%s\" has failed.", cur->string);
                ommp_fatal(msg);
            };
            
            if(mode != 'r' && mode != 'w'){
                sprintf(msg, "Unrecognized file access mode (unexpected error) '%c'.", mode);
                ommp_fatal(msg);
            }
        }

        if(strcmp(cur->string, "xyz_file") == 0){
            if(mode != 'r'){
                sprintf(msg, "xyz_file should be in read mode.");
                ommp_fatal(msg);
            }
            xyz_path = path;
        }
        else if(strcmp(cur->string, "prm_file") == 0){
            if(mode != 'r'){
                sprintf(msg, "prm_file should be in read mode.");
                ommp_fatal(msg);
            }
            prm_path = path;
        }
        else if(strcmp(cur->string, "mmpol_file") == 0){
            if(mode != 'r'){
                sprintf(msg, "mmpol_file should be in read mode.");
                ommp_fatal(msg);
            }
            mmpol_path = path;
        }
        else if(strcmp(cur->string, "output_file") == 0){
            if(mode != 'w'){
                sprintf(msg, "output_file should be in write mode.");
                ommp_fatal(msg);
            }
            output_path = path;
        }
        else if(strcmp(cur->string, "qm") == 0){
            if(*ommp_qmh == NULL){
                if(!smartinput_qm(cur, ommp_qmh))
                    ommp_fatal("Error during creation of QM Helper object");
            }
            else
                ommp_fatal("Only a single qm section can be present in smart input");
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
        else if(strcmp(cur->string, "solver") == 0){
            if(strcmp(cur->valuestring, "default") == 0)
                req_solver = OMMP_SOLVER_DEFAULT;
            else if(strcmp(cur->valuestring, "conjugate gradient") == 0 || strcmp(cur->valuestring, "cg") == 0)
                req_solver = OMMP_SOLVER_CG;
            else if(strcmp(cur->valuestring, "inversion") == 0)
                req_solver = OMMP_SOLVER_INVERSION;
            else if(strcmp(cur->valuestring, "diis") == 0)
                req_solver = OMMP_SOLVER_DIIS;
            else{
                sprintf(msg, "Unrecognized option \"%s\" for solver; Available solvers are default, conjugate gradient, cg, inversion, diis.", cur->valuestring);
                ommp_fatal(msg);
            }
        }
        else if(strcmp(cur->string, "matrix_vector") == 0){
            if(strcmp(cur->valuestring, "default") == 0)
                req_matv = OMMP_MATV_DEFAULT;
            else if(strcmp(cur->valuestring, "direct") == 0)
                req_matv = OMMP_MATV_DIRECT;
            else if(strcmp(cur->valuestring, "incore") == 0)
                req_matv = OMMP_MATV_INCORE;
            else{
                sprintf(msg, "Unrecognized option \"%s\" for matrix_vector; Available solvers are default, direct, incore.", cur->valuestring);
                ommp_fatal(msg);
            }
        }
        else if(strcmp(cur->string, "frozen_atoms") == 0){
            if(!cJSON_IsArray(cur))
                ommp_fatal("frozen_atoms should be an array of integers!");
            if(nfrozen > 0) 
                ommp_fatal("Only a single frozen_atoms section should be present");
            cJSON *_arr = cur->child;
            for(nfrozen = 0; _arr != NULL; _arr = _arr->next){
                if(!cJSON_IsNumber(_arr))
                    ommp_fatal("frozen_atoms should be an array of integers!");
                nfrozen++;
            }
            
            if(nfrozen < 0)
                ommp_fatal("Wired error when parsing frozen_atoms");
            
            frozenat = (int32_t *) malloc(sizeof(int32_t) * nfrozen);
            
            _arr = cur->child;
            for(int i = 0; _arr != NULL; i++){
                frozenat[i] = _arr->valueint;
                _arr = _arr->next;
            }
        }
        else if(strcmp(cur->string, "link_atoms") == 0){
            if(!cJSON_IsArray(cur))
                ommp_fatal("link_atoms should be an array of structures!");
            if(nla > 0)
                ommp_fatal("Only a single link_atoms section should be present");
            cJSON *linkatom_arr = cur->child;
            for(nla = 0; linkatom_arr != NULL; linkatom_arr = linkatom_arr -> next){
                if(!cJSON_IsObject(linkatom_arr))
                    ommp_fatal("link_atoms should be an array of structures!");
                nla++;
            }

            if(nla < 0)
                ommp_fatal("Wired error when parsing link_atoms");
            la_mm = (int32_t *) malloc(sizeof(int32_t) * nla);
            la_qm = (int32_t *) malloc(sizeof(int32_t) * nla);
            la_la = (int32_t *) malloc(sizeof(int32_t) * nla);
            la_ner = (int32_t *) malloc(sizeof(int32_t) * nla);
            la_bl = (double *) malloc(sizeof(double) * nla);

            linkatom_arr = cur->child;
            for(int i=0; linkatom_arr != NULL; i++){
                la_mm[i] = la_qm[i] = la_la[i] = 0;
                la_bl[i] = -1.0;
                la_ner[i] = 0;
                cJSON *tmp;
                for(tmp = linkatom_arr->child; tmp != NULL; tmp = tmp->next){
                    if(strcmp(tmp->string, "MM_id") == 0 && cJSON_IsNumber(tmp))
                        la_mm[i] = tmp->valueint;
                    else if(strcmp(tmp->string, "QM_id") == 0 && cJSON_IsNumber(tmp))
                        la_qm[i] = tmp->valueint;
                    else if(strcmp(tmp->string, "LA_id") == 0 && cJSON_IsNumber(tmp))
                        la_la[i] = tmp->valueint;
                    else if(strcmp(tmp->string, "bond_length") == 0 && cJSON_IsNumber(tmp))
                        la_bl[i] = tmp->valuedouble * OMMP_ANG2AU;
                    else if(strcmp(tmp->string, "eel_remove") == 0 && cJSON_IsNumber(tmp))
                        la_ner[i] = tmp->valueint;
                    else{
                        sprintf(msg, "Unrecognized field %s in link atom %d.", tmp->string, i);
                        ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
                    }
                }
                if(la_mm[i] == 0){
                    sprintf(msg, "MM_id missing in link atom %d.", i);
                    ommp_fatal(msg);
                }
                if(la_qm[i] == 0){
                    sprintf(msg, "QM_id missing in link atom %d.", i);
                    ommp_fatal(msg);
                }
                if(la_la[i] == 0){
                    sprintf(msg, "LA_id missing in link atom %d.", i);
                    ommp_fatal(msg);
                }

                linkatom_arr = linkatom_arr -> next;
            }
        }
        else{
            sprintf(msg, "Unrecognized JSON element \"%s\".", cur->string);
            ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
        }

        cur = cur->next;
    }

    // Set verbosity
    ommp_set_verbose(req_verbosity);
    // Set output file
    if(output_path != NULL)
        ommp_set_outputfile(output_path);

    // Print information from JSON
    sprintf(msg, "Smart Input Name: %s", json_name);
    ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
    
    sprintf(msg, "Smart Input Description: %s", json_description);
    ommp_message(msg, OMMP_VERBOSE_LOW, "SI");
    
    // Input for MM
    if(xyz_path != NULL){
        ommp_message("Trying initialization from Tinker .xyz file.", OMMP_VERBOSE_LOW, "SI");
        if(prm_path == NULL)
            ommp_fatal("xyz_file set but prm_file is missing.");
        if(hdf5_path != NULL)
            ommp_fatal("xyz_file set but also hdf5_file is set, this is ambiguous.");
        if(mmpol_path != NULL)
            ommp_fatal("xyz_file set but also mmpol_file is set, this is ambiguous.");

        *ommp_sys = ommp_init_xyz(xyz_path, prm_path);
    }
    else if(mmpol_path != NULL){
        if(prm_path != NULL)
            ommp_fatal("prm_file set but it is not needed for mmpol input, this is ambiguous.");
        if(hdf5_path != NULL)
            ommp_fatal("mmpol_file set but also hdf5_file is set, this is ambiguous.");

        *ommp_sys = ommp_init_mmp(mmpol_path);
    }
    else if(hdf5_path != NULL){
        if(prm_path != NULL)
            ommp_fatal("prm_file set but it is not needed for mmpol input, this is ambiguous.");

        ommp_fatal("input from hdf5_file is still not implemented.");
    }
    else{
        ommp_fatal("No input for MM system found in Smart Input file, set one of xyz_file+prm_file, mmpol_file, hdf5_file");
    }
    
    // Set solver in ommp_sys
    ommp_set_default_solver(*ommp_sys, req_solver);
    // Set matv in ommp_sys
    ommp_set_default_matv(*ommp_sys, req_matv);

    // Handle QM part of the system
    if(*ommp_qmh == NULL){
        ommp_message("qm section not present in smart input JSON: QMHelper object is not set.",
                     OMMP_VERBOSE_LOW, "SI");
    }

    // Handle frozen atoms
    if(nfrozen > 0){
        ommp_set_frozen_atoms(*ommp_sys, nfrozen, frozenat);
        free(frozenat);
    }
    // Handle link atoms
    if(nla > 0){
        if(*ommp_qmh == NULL)
            ommp_fatal("Link atoms requested but no qm section is defined!");
        for(int i=0; i<nla; i++){
            if(la_bl[i] < 0) la_bl[i] = OMMP_DEFAULT_LA_DIST;
            if(la_ner[i] == 0) la_ner[i] = OMMP_DEFAULT_LA_N_EEL_REMOVE;

            ommp_create_link_atom(*ommp_qmh,  *ommp_sys, 
                                  la_mm[i], la_qm[i], la_la[i], 
                                  prm_path, 
                                  la_bl[i], la_ner[i]);
        }

        free(la_mm);
        free(la_qm);
        free(la_la);
        free(la_bl);
        free(la_ner);
    }

    return;
}

void *c_json_cherrypick(const char *json_file, char *path, char type){
    if(type != 'd' && type != 's' && type != 'i')
        return NULL;
    
    char msg[OMMP_STR_CHAR_MAX];
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
    
    int32_t *outi = NULL;
    char *outs = NULL;
    double *outd = NULL;

    char *field = strtok(path, "/");
    
    while(field != NULL){
        while(cur != NULL){
            if(strcmp(cur->string, field) == 0)
                break;
            cur = cur->next;
        }
        field = strtok(NULL, "/");
        
        if(cur == NULL){
            ommp_message("Path not found!",
                            OMMP_VERBOSE_LOW, "SI");
            return NULL;
        }
        
        if(cJSON_IsObject(cur)){
            // This is not a leaf
            if(field == NULL){
                ommp_message("Incomplete path provided for SI cherry picking.",
                             OMMP_VERBOSE_LOW, "SI");
                return NULL;
            }
            cur = cur->child;
        }
        else{
            // This is a leaf
            if(field == NULL){
                // Here we are, we found it finally!
                switch(type){
                    case 'i':
                        if(cJSON_IsNumber(cur)){
                            outi = (int32_t *) malloc(sizeof(int32_t));
                            *outi = cur->valueint;
                            return (void *) outi;
                        }
                        break;
                    case 'd':
                        if(cJSON_IsNumber(cur)){
                            outd = (double *) malloc(sizeof(double));
                            *outd = cur->valuedouble;
                            return (void *) outd;
                        }
                        break;
                    case 's':
                        if(cJSON_IsString(cur)){
                            outs = (char *) malloc(sizeof(char) * strlen(cur->valuestring));
                            strcpy(outs, cur->valuestring);
                            return (void *) outs;
                        }
                        break;
                    default:
                        break;
                }
                return NULL;
            }
            else{
                ommp_message("Path not found!",
                             OMMP_VERBOSE_LOW, "SI");
                return NULL;
            }
        }
    }
}

void c_smartinput_cpstr(const char *json_file, char *path, char **s){
    *s = c_json_cherrypick(json_file, path, 's');
    if(*s == NULL){
        ommp_fatal("JSON cherry picking failed.");
    }
}

void ommp_smartinput_cpstr(const char *json_file, char *path, char **s){
    c_smartinput_cpstr(json_file, path, s);
}

void ommp_smartinput(const char *json_file, OMMP_SYSTEM_PRT ommp_sys, OMMP_QM_HELPER_PRT ommp_qmh){
    // Just an interface function to expose same names and functionalities in C and Fortran
    c_smartinput(json_file, ommp_sys, ommp_qmh);
}
