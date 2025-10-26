docker {
    enabled = true
    runOptions= "--entrypoint=''"
}

params {
    output_folder = 'results'
}

executor { 
    cache = true
}

report {
    enabled = true
    file = 'reports/report.html'
    overwrite = true
}
trace {
    enabled = true
    file = 'reports/trace.txt'
    overwrite = true
}
timeline {
    enabled = true
    file = 'reports/timeline.html'
    overwrite = true
}
dag {
    enabled = true
    file = 'reports/dag.dot'
    overwrite = true
}

process {
    executor = 'local' 
    //cpus = 2
    //memory = '1 GB'
    maxRetries = 0
    errorStrategy = 'retry'

    //Précision des containers pour les processus
    withName: GET_SRR {
        container="mariemeier/reprohackathon:entrez-direct"
    }
    withName: DOWNLOAD_FASTQ {
        container="mariemeier/reprohackathon:sra-toolkit"
    }
}
