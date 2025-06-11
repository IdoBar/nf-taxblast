#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 
/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
 

def helpMessage() {

     log.info """
      Usage:
      The typical command for running the pipeline is as follows:
      nextflow run nf-taxblast.nf --app blastn --query QUERY.fasta --chunkSize 200 --db "path-of-db/db" -profile conda,blastn_tax
      nextflow run nf-taxblast.nf --app "diamond blastp" --query QUERY.faa --chunkSize 5000 --db "path-of-db/db" -profile docker,diamond_tax

      Mandatory arguments:
       --app <value>                  BLAST/DIAMOND program to use (diamond blastp/x must be quoted!)
                                      Valid options: [blastn, blastp, tblastn, blastx, 'diamond blastp', 'diamond blastx'] 
       --query <file.fatsa>           Query fasta file of sequences you wish to BLAST
       --db <path-of-db/db>           Path of the BLAST or DIAMOND database. 
                                      If BLAST database is provided for DIAMOND and taxonomy information is requested
                                      then a suitable database will be created (see Taxonomy options below). 
                                      Default: [$BLASTDB/nt or $BLASTDB/nr for protein search]   
       -profile <profile1,profile2>   Configuration profile to use. Can use multiple (comma separated)
                                      Available profiles for container systems: [conda/apptainer/singularity/local/docker]
                                      Available profiles with preset database and output formats: [blastn_tax/diamond_tax]
                                      Available profiles with test datasets and databases: [test/test_tax/test_p/test_p_tax/test_d/test_d_tax]

       Optional arguments:
       --out <outfile.outfmt6>        Output filename of final BLAST output. Default: [QUERY.app.db.outfmt6]
       --outDir <path>                Output folder for the results. Default: [results]
       --outCols <'std'>              Output columns (must be quoted!). Default: ['std']
       --headers <false>              Include headers in the output table. Default: false
       --blastOpts <'-evalue 10'>     Additional options for BLAST command (must be quoted!). Default: ['-evalue 1e-10 -max_target_seqs 20']
       --dmndOpts <'-e 10e-10'>       Additional options for BLAST command (must be quoted!). Default: ['-e 1e-10 -k 20'] 
       --chunkSize <num>              Number of fasta records to use in each job when splitting the query fasta file. Default: [250]
                                      This option can also take the size of each subquery (like 200.KB, 5.KB, etc.) 
       --queueSize <num>              Maximum number of jobs to be queued [50]
       --download <false>             Download database before running homology search. Default: false

       Taxonomy options:
       --taxDbDir <path-of-taxdb/db>  Location of taxonomy db files (prot.accession2taxid.FULL.gz, nodes.dmp and names.dmp) to allow DIAMOND 
                                      to return taxonomic information columns. If the required files cannot be found in the path 
                                      they will be automatically downloaded from the NCBI.
                                      Information about the required files and where to download them can be found at 
                                      https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options
                                      Default: [same path as the database]
       --taxListFile <taxid.list>     A file with list of taxonomy IDs to limit the search space.
       --help                         This usage statement.
     """
 }



 // Show help message
 if (params.help) {
     helpMessage()
     exit 0
 }

// def db_prefix = file(param.dbDir).resolve(param.dbName)
def db_name = file(params.db).name
// diamond_db = file(params.db).name
// db_path = (params.app =~ /D$/ && params.outfmtString =~ /staxids|sscinames|scomnames|sskingdoms/ && file("${params.db}.dmnd").exists()) ? "${params.db}.dmnd" : file(params.db).name
def db_basename = file(params.db).baseName
def db_dir = file(params.db).parent ?: "$launchDir/db"
def db_prefix = db_dir.resolve("${db_basename}")
def tax_db_dir = params.taxDbDir ?: db_dir
// out_dir = "${params.outDir}/${params.app}"
// Check if the chunks are provided as Memory units and if not assume KB
// =~ /\d+\.*\w[bB]$/ ? MemoryUnit.of( "${params.chunkSize}.KB" ) : MemoryUnit.of( params.chunkSize )
if (params.chunkSize instanceof String) {
    if (params.chunkSize.isNumber()) {
        int chunk_size = params.chunkSize.toInteger()
    } else if (params.chunkSize ==~ /\d+\.*\w[bB]$/) {
    (mem_value, mem_suffix) = (chunk_size =~ /(\d+)\.*(\w[bB])$/)[0]
    def chunk_size = mem_value + "." + mem_suffix.toUpperCase()
} else {error("`${params.chunkSize}` is not a valid chunk size, please provide an integer (i.e. 200) or a file size (i.e. '200.KB')")}
} else if (params.chunkSize instanceof Number) {
    int chunk_size = params.chunkSize
}



// Format output filename with as "query"
// def output_fmt = params.outfmtString =~ /(\d+) (.+)/ 
def dmnd_opts = params.headers ? params.dmndOpts + " --header simple" : params.dmndOpts
def blast_outfmt = params.headers ? "7 " + params.outCols : "6 " + params.outCols
def tax_filt = file(params.taxListFile)
def app_str = params.app.replaceAll(/\s/, "-")
def std = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
def suffix = params.headers ? 'outfmt7' : 'outfmt6'
def res_name = params.out ?: "${file(params.query).baseName}.${app_str}.${db_basename}.${suffix}"
def db_ready = false
// dag_file = (nextflow.version.matches('22.04+')) ? "${params.outDir}/${params.tracedir}/blast-nf_dag.mmd" : "${params.outDir}/${params.tracedir}/blast-nf_dag.html"
/*
PROCESSES
*/
process blast {
    label 'blast'

    publishDir params.outDir, mode: 'copy', pattern: '*.version'

    input:
        path 'query.fa'
        val dbPrefix 
       // val dbName
        path taxlist
 
    output:
        path 'blast_results'
       // path 'software_versions.txt'

    script:
    def taxid_filt = taxlist.name != 'NO_FILE' ? "-taxidlist ${taxlist}" : ''
    """
    ${params.app} -db ${dbPrefix} -query 'query.fa' -outfmt '${blast_outfmt}' ${taxid_filt} -num_threads ${task.cpus} ${params.blastOpts} > blast_results
    blastn -version | head -n1 > blast.version
    """
}

process taxids_file2list {
    executor 'local'
    input:
        path taxid_file

    output:
        stdout

    shell:
    '''
    tax_str=!( gawk -v ORS="," '1' !{taxid_file} | sed 's/,$//' )
    echo "--taxonlist ${tax_str}"
    '''
}


process extract_db_fasta {
    label 'blast'
    label 'sc_medium'
   // publishDir db_dir, mode: 'copy', overwrite: true
    input:
        val dbPrefix 

    output:
        path "${db_name}.fasta" 
    script:
        def basename = file(dbPrefix).baseName
    shell:
    """
    blastdbcmd -entry all -db ${dbPrefix} -out ${db_name}.fasta
    """

}

process download_db {
    label 'blast'

    publishDir "${db_dir}", mode: 'copy', pattern: "${db_name}.*"
    publishDir params.outDir, mode: 'copy', pattern: '*.version'

    // input:
    //     val dbName
 
    output:
        path "${db_name}.*"
        // path 'software_versions.txt'

    """
    update_blastdb.pl --decompress --verbose --source ncbi --num_threads 0 ${db_name}
    rm ${db_name}*.tar.gz
    update_blastdb.pl --version | head -n1 > update_blastdb.version
    """
}
process download_tax_db {
    label 'download'

    publishDir tax_db_dir, mode: 'copy', overwrite: true 

    output:
        path "*.dmp"
        path "prot.accession2taxid.FULL.gz"

    script:
    """
    aria2c -x5 -c --auto-file-renaming=false https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
    aria2c -x5 -c --auto-file-renaming=false https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    unzip taxdmp.zip
    """
}

process diamond_prep_tax_db {
    label 'diamond'
    publishDir params.outDir, mode: 'copy', pattern: '*.version'
    // publishDir "${tax_db_dir}", mode = 'copy' 

    input:
        val dbPrefix
        path fasta
        path acc2tax
        path names
        path nodes
    
    output:
        val true  
    
    script:
    """
    diamond makedb -p ${task.cpus} --in ${fasta} -d ${dbPrefix} --taxonmap ${acc2tax} --taxonnodes ${nodes} --taxonnames ${names}
    diamond version > diamond.version
    """

}

process diamond_prep_blast_db {
    label 'diamond'
    // publishDir "${dbDir}", mode = 'copy' 
    publishDir params.outDir, mode: 'copy', pattern: '*.version'

    input:
        // val dbDir
        // val dbName
        val dbPrefix // from "${db_name}"
    
    output:
        val true  
    
    script:
    """
    diamond prepdb -p ${task.cpus} -d ${dbPrefix}
    diamond version > diamond.version 
    """

}


process diamond {
    label 'diamond'

    publishDir params.outDir, mode: 'copy', pattern: '*.version'

    input:
        path 'query.fa' // from ch_fasta
        val dbPrefix //from diamond_db_ch
        // val dbName
        val tax_filt_string
    output:
        path 'diamond_results'
       // path 'software_versions.txt'

    script:
    """
    ${params.app} -d ${dbPrefix} -q 'query.fa' ${tax_filt_string} -f 6 ${params.outCols} -p ${task.cpus} ${dmnd_opts} -o diamond_results
    diamond version > diamond.version
    """
} 

/*

workflow diamond_prep_flow { 
    take: 
        fmt_string
        db_prefix
        tax_db_dir
        tax_filt
        // db_ready
    main:
        if (fmt_string =~ /staxids|sscinames|sphylums|skingdoms/) {
            // def diamond_db = db_dir.resolve("${db_basename}.dmnd")
            if (file("${db_prefix}.dmnd").isEmpty()) {
                def acc2taxid = file("${tax_db_dir}/prot.accession2taxid.FULL.gz")
                def names_dmp = file("${tax_db_dir}/names.dmp")
                def nodes_dmp = file("${tax_db_dir}/nodes.dmp")
                db_fasta_ch = file("${db_prefix}.fasta").isEmpty() ? extract_db_fasta(db_prefix) : Channel.fromPath("${db_prefix}.fasta")
                if (acc2taxid.isEmpty() || names_dmp.isEmpty() || nodes_dmp.isEmpty()) 
                    download_tax_db()
                
                diamond_prep_tax_db(db_prefix, db_fasta_ch, Channel.fromPath(acc2taxid), Channel.fromPath(names_dmp), Channel.fromPath(nodes_dmp))
            } 
        } 
        else if (file("${db_prefix}.acc").isEmpty()) 
                diamond_prep_blast_db(db_prefix)
        // process taxonomy filtering string
        def tax_filt_d (tax_filt_file) {
            if (tax_filt_file.name == 'NO_FILE')
                return ''
            tax_file_str = tax_filt_file.text
            tax_str = tax_file_str.replaceAll(/\n/, ",")
            tax_str = "--taxonlist " + tax_str.replaceFirst(/,$/, "")
            return tax_str
            //tax_filt_file.text.replaceAll(/\n/, ",").replaceFirst(/,$/, "")        
        }
        tax_filt_string = tax_filt_d(tax_filt)
        // tax_filt.name != 'NO_FILE' ? "--taxonlist ${taxids_file2list(tax_filt)}" : ''
            
        db_ready = true
        
    emit:
        tax_filt_string
        db_ready

}
*/
workflow {
    
    // Define your channels
    if (chunk_size ==~ /\d+\.\w+[bB]/) {
        Channel.fromPath(params.query)
                .splitFasta(size: MemoryUnit.of( chunk_size ), file: true)
                .set { ch_fasta }
    }
    else {
        Channel.fromPath(params.query)
                .splitFasta(by: chunk_size, file: true)
                .set { ch_fasta }
    }
    
    
    
    // db_ch = Channel.fromPath(params.db)
    // Check if database exists
    if (params.app =~ /diamond|blastp/ && !file("${db_prefix}.pdb").exists() && !params.download) {
        error("Unable to find protein database (${db_prefix}.pdb), please rerun the pipeline with '--download' flag first")

    } else if (params.download) {
        download_tax_db()
    }
    if (params.app =~ /blastn/ && !file("${db_prefix}.ndb").exists() && !params.download) {
        error("Unable to find nucleotide database (${db_prefix}.ndb), please rerun the pipeline with `--download` flag first")

    } else if (params.download) {
        download_tax_db()
    }
    // Run the homology search process 
    if (params.app =~ /diamond/) {
        dmnd_db = "${db_prefix}"
        tax_str_ch = tax_filt.name == 'NO_FILE' ? '' : taxids_file2list(tax_filt) 
        if (params.outCols =~ /staxids|sscinames|sphylums|skingdoms/) {
            // def diamond_db = db_dir.resolve("${db_basename}.dmnd")
            if (!file("${db_prefix}.dmnd").exists()) {
                acc2taxid = file("${tax_db_dir}/prot.accession2taxid.FULL.gz")
                names_dmp = file("${tax_db_dir}/names.dmp")
                nodes_dmp = file("${tax_db_dir}/nodes.dmp")
                // download acc2tax if missing
                if (acc2taxid.isEmpty() || names_dmp.isEmpty() || nodes_dmp.isEmpty()) 
                    download_tax_db()
                db_fasta_ch = file("${db_prefix}.fasta").isEmpty() ? extract_db_fasta(db_prefix) : Channel.fromPath("${db_prefix}.fasta")
                diamond_prep_tax_db(db_prefix, db_fasta_ch, acc2taxid, names_dmp, nodes_dmp)
            }
            dmnd_db = "${db_prefix}.dmnd"
        } 
        else if (!file("${db_prefix}.acc").exists()) {
            diamond_prep_blast_db(db_prefix)
            // dmnd_db = "${db_prefix}.acc"
        } 
     
        ch_hits = diamond(ch_fasta, dmnd_db, tax_str_ch)      
        // tax_filt_ch = Channel.fromPath(params.taxListFile)        
        // run diamond
        // (tax_filt_ch, db_ready_ch) = diamond_prep_flow(Channel.value(db_prefix), Channel.value(tax_db_dir), Channel.value(tax_filt), Channel.value(db_ready))
        // if (diamond_prep_flow.out[1]) {
        
        
        //}
            
        /*
        if (diamond_prep_flow.out[1])
            ch_hits = diamond(ch_fasta, db_prefix, diamond_prep_flow.out[0], diamond_prep_flow.out[1])
            */
    }
    else {
        ch_hits = blast(ch_fasta, db_prefix, tax_filt)
    }
         
    if (params.headers){
        ch_hits
        .collectFile(name: res_name, storeDir: params.outDir, keepHeader: true)
        .subscribe { println "Entries are saved to file: $it" }
    } else {

    }     
    // Collect results
    ch_hits
        .collectFile(name: res_name, storeDir: params.outDir)
        .subscribe { println "Entries are saved to file: $it" }
}

