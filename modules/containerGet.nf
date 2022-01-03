process containerGet {
    //label 'basic_tools'
    label 'smallTask'
    tag "$tool"

    errorStrategy 'retry'
    maxRetries 10

    if ( workflow.profile.contains('singularity') ) { storeDir "${params.singularityCacheDir}" }
    else if ( workflow.profile.contains('conda') ) { storeDir "${params.condaCacheDir}" }

  input:
    tuple val(tool), val(path)

  output:
    file("*/*${tool}*")

  script:
  
  if ( workflow.profile.contains('singularity') )
    """
    if [ -e ${params.singularityCacheDir}/${path}*.img ] 
        then
            echo "${tool} singularity image file already exists, skipping."
    else
        singularity pull --name ${path}.img --dir "${params.singularityCacheDir}" "docker://${path}"
    fi
    """
  else if ( workflow.profile.contains('docker') )
    """
    if [[ "\$(docker images -q ${path}* 2> /dev/null)" == "" ]]; 
    then
	    echo "${tool} docker image file already exists, skipping."
	else
      docker pull ${path}
    fi
    """
  else if ( workflow.profile.contains('conda') )
    """
    if [ -e ${params.condaCacheDir}/${tool}* ] 
      then
        echo "${tool} conda environment already exists, skipping."
    else
      yes | conda env create -f ${path}
    fi
    """
  else 
    """
    echo 'unknown container engine, exiting'
    exit 1
    """
}