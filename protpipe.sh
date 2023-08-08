#!/usr/bin/env bash

#### CHECK SINGULARITY #############################################################################
if ! command -v singularity &> /dev/null; then
    echo "INFO: singularity command not found, looking for singularity module"
    if ! command -v module &> /dev/null; then
        echo "ERROR: module command not found. Did you mean to run this on an HPC?"
        exit 1
    else
        if $(module avail singularity/3 2>&1 >/dev/null | grep -q 'No module'); then
            echo 'ERROR: singularity cannot be found. Recheck installation?'
            exit 1
        else
            echo 'INFO: module singularity found'
            module load singularity/3
        fi
    fi
else
    echo 'INFO: singularity command found'
fi


#### R ANALYSIS ####################################################################################
r_sif="src/R.sif"
r_sif_sha256_desired='523df8eb1aa5c9eb30ac5cd70c3181431f90af6e12e497ee6d43854e5f578196'
if [ ! -f "${r_sif}" ]; then
    echo "INFO: pulling image from Singularity cloud"
    singularity pull --arch amd64 library://wellerca/protpipe/4.3.1:sha256.523df8eb1aa5c9eb30ac5cd70c3181431f90af6e12e497ee6d43854e5f578196
    mv 4.3.1_sha256.523df8eb1aa5c9eb30ac5cd70c3181431f90af6e12e497ee6d43854e5f578196.sif "${r_sif}"
else
    echo "INFO: ${r_sif} already exists, skipping download"
fi

r_sif_sha256_actual=$(sha256sum "${r_sif}" | awk '{print $1}')

if [ ! "${r_sif_sha256_actual}" == "${r_sif_sha256_desired}" ]; then
    echo "ERROR: ${r_sif} sha256sum does not pass check. Possibly corrupted? Delete or clear singularity cache and try again."
    exit 1
else
    echo "INFO: ${r_sif} sha256sum sum passes check"
fi

#### RUN IN CONTAINER ####################################################################################

singularity exec --cleanenv -H ${PWD} ${r_sif} Rscript $@