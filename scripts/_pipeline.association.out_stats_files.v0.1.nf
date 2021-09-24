#!/data/0212964/nextflow/nextflow-0.25.1-all

/* ######################################################################
#                 Association Pipeline v0.1                             #
###################################################################### */

/*
 * Defines pipeline parameters
 */

params.final_chunks = "The final chunking file splited by chromosome"
params.input_pheno = "The phenotype to be used in this analyses, file generated in ESN Windows environment"
params.imputation_dir = "The directory contains the post imputation genotype files (mldose files)"
params.example_mldose = "One of the mldose files in the imputation directory"
params.geno_prefix = "The prefix of the mldose file"
params.working_dirs = "The working directory"
params.out = "The output file name for the stats file"
params.method = "palinear or palogist"

/*
 * probabel_header = "name\tchrom\tposition\tA1\tA2\tFreq1\tMAF\tQuality\tRsq\tn\tMean_predictor_allele\tbeta_SNP_add\tsebeta_SNP_add\tchi2_SNP\tchi\tp\tor_95_percent_ci"
 */

// create chunks channels from chunking file
chunks = Channel
			.from( file(params.final_chunks) )
			.splitCsv(header:['chr', 'chunk', 'start', 'end'], skip: 1, sep: '\t')
chunks.into { chunks_prep_geno; chunks_merge; }

/* *********************************************
 * Step 1: Start Prepare ProbABEL Phenotype File
 */

process prepare_pheno{

	input:
	// input defined by parameters

	output:
	file "probabel_pheno" into probabel_phenotype_file

	memory '8 GB'
	executor 'pbs'
	
    """
	/share/storage/REDS-III/common/software/prepare_probabel_files.pl \
		--in_mldose ${params.example_mldose} \
		--in_pheno ${params.input_pheno} \
		--out_pheno "probabel_pheno"
	"""
}

/* ********************************************
 * Step 2: Start Prepare ProbABEL Genotype File
 */

process prepare_geno{

	input:
	file "probabel_pheno" from probabel_phenotype_file
	set chr, chunk, start, end from chunks_prep_geno

	output:
	file "${params.geno_prefix}${chr}.${chunk}.mach_mldose" into probabel_genotype_files

    memory '20 GB'
    executor 'pbs'

	"""
	/share/storage/REDS-III/common/software/prepare_probabel_files.pl \
		--in_mldose "${params.imputation_dir}/${params.geno_prefix}${chr}.${chunk}.mach.mldose.gz" \
		--in_pheno "${probabel_pheno}" \
		--out_mldose "${params.geno_prefix}${chr}.${chunk}.mach_mldose" 
	"""
}

/* *******************************
 * Step 3: Start ProbABEL Analysis
 */

process do_probabel_analysis{

	input:
	file probabel_pheno from probabel_phenotype_file
	file geno from probabel_genotype_files
	
	output:
	file "${geno.baseName}_add.out.txt" into probabel_out_files

	memory '20 GB'
	executor 'pbs'

	"""
	/share/storage/REDS-III/common/software/${params.method} \
    	--pheno $probabel_pheno \
    	--dose $geno \
    	--info "${params.imputation_dir}/${geno.baseName}.mach.mlinfo" \
    	--map "${params.imputation_dir}/${geno.baseName}.legend" \
    	--out "${geno.baseName}"
	"""
}

/* ***************************
 * Step 4: Start calc chi p OR
 */

process calc_chi_p_or{

	input:
	file result from probabel_out_files
	
	output:
	file "${result.baseName}.stats" into stats_files

	memory '8 GB'
	executor 'pbs'

	"""
	/share/storage/REDS-III/common/software/R/calculate_stats_for_probabel_results_v2.R \
		--remove_missing_p \
		--in_file "${result}" \
		--out_file "${result.baseName}.stats"
	"""
}

/* ***************************
 * Step 5: merge stats files
 */

stats_files
	.collectFile(name: file(params.out), skip: 1)
	.println {"Result saved to file: $it"}
