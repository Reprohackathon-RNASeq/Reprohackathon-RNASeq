params {
    output_folder = 'results'
}

process {
    // Exécuteur par défaut: local (votre VM).
    executor = 'local' 
    
    // Configuration des ressources par défaut (optionnel, mais bonne pratique)
    //cpus = 2
    //memory = '1 GB'
    maxRetries = 30
    errorStrategy = 'retry'
}

docker {
    enabled = true
}
