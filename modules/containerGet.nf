process containerGet {
    if (!params.cloudProcess) { label 'smallTask' }
    tag "$tool"

    errorStrategy 'retry'
    maxRetries 10

    if ( workflow.profile.contains('singularity') ) { storeDir "${params.singularityCacheDir}" }
    else if ( workflow.profile.contains('conda') ) { storeDir "${params.condaCacheDir}" }

  input:
    tuple val(tool), val(path)

  output:
    path("${img_file_name}.img")

  script:
  img_file_name = path.replace("/", "-").replace(":","-")
  if ( workflow.profile.contains('singularity') )
    """
    if [ -e ${params.singularityCacheDir}/${img_file_name}.img ] 
        then
            echo "${tool} singularity image file already exists, skipping."
    else
        singularity pull --name ${img_file_name}.img "docker://${path}"
    fi
    """
  else if ( workflow.profile.contains('docker') )
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
  else if ( workflow.profile.contains('conda') )
    """
    echo "The setup mode currently doesn't support conda environments, skipping."
    touch ${img_file_name}.img
    """
  else 
    """
    echo 'unknown container engine, exiting'
    exit 1
    """
}