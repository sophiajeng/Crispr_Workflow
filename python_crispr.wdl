version 1.0

workflow crispr_python_workflow {
	input {
		String sample_name_file
		String library_file
		String library_config_file
		String annot_file
		String script_dir
	}
		#TODO: specify directory where the scripts are

		#library and fastq files
		#this can be used to create a pair that can be transformed into a map once we have 1.1 support
		#Array[String] sample_library_array =  read_lines(sample_library_file)

		Array[Array[String]] sample_name_array = read_tsv(sample_name_file)

		Array[String] library_array = read_lines(library_file)
		Array[Array[String]] library_config_array = read_tsv(library_config_file)

		#Map[String,Array[Array[String]]] sample_library_name_map = collect_by_key(sample_library_name_pair)
		#Map[String, Array[String]] library_config_map = {'yusa':library_config_array}

		#TO DO: specify library
		Map[String, Array[Array[String]]] sample_library_name_map = {"eil_20_21":sample_name_array}

		Array[Pair[String, Array[String]]] library_config_map = zip(library_array, library_config_array)

		scatter (pairs in library_config_map) {
			String library = pairs.left
			Array[String] library_params= pairs.right
			scatter (sample in sample_library_name_map[library]){
				call cutadapt {
					input:
						cur_sample = sample[1],
						sample_file = sample[0],
						five_prime = library_params[1],
						three_prime = library_params[2]
				}
				call count_reads{
					input: 
						cur_sample = sample[1],
						trimmed_reads=cutadapt.trimmed_fastq,
						script_dir = script_dir,
						sgRNA_file = library_params[0]
				}
			}
		}

                      call basicr {
                                input:
                                        counts_files=count_reads.counts_txt,
                                        stats_files=count_reads.stat_txt,
                			all_counts_files=count_reads.all_counts_txt,
                			trimmed_seq_files=count_reads.trimmed_seq_txt,
                                        script_dir = script_dir,
                                        lib_annot = annot_file
                        }
 
                output {        
                	File counts_output = basicr.count_txt
                	File stats_output = basicr.stat_txt
                        File trimmed_output = basicr.trimmed_hist
                        File all_output = basicr.all_length

        	}

}
task cutadapt {
	input {
		String cur_sample
		File sample_file
		String five_prime
		String three_prime
	}
	command {
		set -e
		cutadapt -j 2 -n 2 -g ${five_prime} -a ${three_prime} ${sample_file} > Sample_${cur_sample}_trimmed.fastq
	}
	output {
		File trimmed_fastq="Sample_" + "${cur_sample}" + "_trimmed.fastq"
	}
}
task count_reads {
	input {
		String cur_sample
		File trimmed_reads
		String script_dir
		String sgRNA_file
	}
	
	command {
		set -e
		python ${script_dir}/crispr_counts.py ${trimmed_reads} ${sgRNA_file} Sample_${cur_sample}_statistics.txt Sample_${cur_sample}_counts.txt Sample_${cur_sample}_trimmed_counts.txt Sample_${cur_sample}_sequence_counts.txt
	}
	runtime {
		runtime_minutes: 600
		requested_memory_per_core: "30G"
		cpus: 9
		maxRetries: 1
	}
	output {
                File stat_txt = "Sample_" + "${cur_sample}" + "_statistics.txt"
		File counts_txt = "Sample_" + "${cur_sample}" + "_counts.txt"
                File all_counts_txt = "Sample_" + "${cur_sample}" + "_sequence_counts.txt"
                File trimmed_seq_txt = "Sample_" + "${cur_sample}" + "_trimmed_counts.txt"
	}
}
task basicr {
	input {
		Array[Array[String]] counts_files
		Array[Array[String]] stats_files
                Array[Array[String]] all_counts_files
                Array[Array[String]] trimmed_seq_files
                Array[String] counts = flatten(counts_files)
                Array[String] stats = flatten(stats_files)
                Array[String] all_counts = flatten(all_counts_files)
                Array[String] trimmed_seqs = flatten(trimmed_seq_files)

		String script_dir
		String lib_annot

	}
	command {
		set -e
		R CMD BATCH --vanilla '--args --counts ${sep=' ' counts} --stats ${sep=' ' stats}  --correct_len ${sep=' ' trimmed_seqs} --all_len ${sep=' ' all_counts} --counts_output current_cell_line_crispr_counts.txt --stats_output current_cell_line_crispr_stats.txt --annotation ${lib_annot} --correct_length_output current_cell_line_crispr_trimmed_hist.txt --all_length_output current_cell_line_crispr_all_len_hist.txt' ${script_dir}/basicRPython.R

	}
	runtime {
		runtime_minutes: 600
		requested_memory_per_core: "10G"
		cpus: 1
		maxRetries: 1
	}
	output {
		File count_txt = "current_cell_line_crispr_counts.txt"
		File stat_txt = "current_cell_line_crispr_stats.txt"
		File trimmed_hist = "current_cell_line_crispr_trimmed_hist.txt"
		File all_length = "current_cell_line_crispr_all_len_hist.txt"
	}
}
