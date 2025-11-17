process GET_GEO_TABLE {
    publishDir "data/geo_table", mode: 'copy', overwrite: true

    input:
      val geo_id

      output:
      path "${geo_id}_table.tsv", emit: geo_table

      script:
      """
      # Download the SOFT file from GEO
      geo_parent=\$(echo ${geo_id} | sed -E 's/(GSE[0-9]{3})[0-9]+/\\1nnn/')
      wget -q ftp://ftp.ncbi.nlm.nih.gov/geo/series/\${geo_parent}/${geo_id}/soft/${geo_id}_family.soft.gz -O ${geo_id}.soft.gz
      gunzip -f ${geo_id}.soft.gz
      
      # Extract GSM IDs and conditions from the SOFT file
      awk '
      /^\\^SAMPLE/ {gsm=\$3}
      /^!Sample_description/ {
          desc=\$3
          if (desc ~ /^IP/) cond="IP"
          else if (desc ~ /^ctrl/) cond="control"
          else cond="unknown"
          print gsm"\\t"cond
      }' ${geo_id}.soft > ${geo_id}_table.tsv

      """
}
