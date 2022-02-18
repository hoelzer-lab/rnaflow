process containerGet {
    label 'smallTask'
    tag "$tool"

    //errorStrategy 'retry'
    //maxRetries 10

    if ( workflow.profile.contains('singularity') ) { storeDir "${params.singularityCacheDir}" }
    else if ( workflow.profile.contains('conda') ) { storeDir "${params.condaCacheDir}" }

  input:
    tuple val(tool), val(path), val(conda_env_suffix)

  output:
    path("${img_file_name}.img")

  script:
  img_file_name = path.replace("/", "-").replace(":","-")

  if (workflow.profile.contains('conda') && !tool.contains('rattle')){
    conda_env_file = path.replace('$baseDir','')
    img_file_name = path.replace('$baseDir','').replace('/','_')
  }

  if ( workflow.profile.contains('singularity') || (workflow.profile.contains('conda') && tool.contains('rattle')))
    """
    if [ -e ${params.singularityCacheDir}/${img_file_name}.img ] 
        then
            echo "${tool} singularity image file already exists, skipping."
    else
        singularity pull --name ${img_file_name}.img "docker://${path}"
    fi
    """
  else if ( workflow.profile.contains('docker') || (workflow.profile.contains('conda') && tool.contains('rattle')) )
    """
    if [[ "\$(docker images -q ${path} 2> /dev/null)" == "" ]]; 
    then
      docker pull ${path}
      touch ${img_file_name}.img
	  else
	    echo "${tool} docker image file already exists, skipping."
      touch ${img_file_name}.img
    fi
    """
  else if ( workflow.profile.contains('conda') && !tool.contains('rattle'))
    """
    if [ -e ${params.condaCacheDir}${tool}-${conda_env_suffix} ] 
        then
            echo "${tool} conda environment file already exists, skipping."
    else
        conda env create -f ${projectDir}${conda_env_file} -p ${params.condaCacheDir}/${tool}-${conda_env_suffix}
    fi
    touch ${img_file_name}.img
    """
  else 
    """
    echo 'unknown container engine, exiting'
    exit 1
    """
}