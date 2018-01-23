
workflow chia-pet {
    # inputs
    # assuming we have a mandatory configuration file:
    File config_file
    Array[File] fastqs 	# [end_id]

    # splitting two fastq files (one paired end run result)
    # into three types of paired FASTQs:
    # (A) none (no linker detected)
    # (B) linker.single (linker found, but only one side was mappable)
    # (C) linker.paired (linker found and both sides are mappable)
    call fastq_splitting { input :
        fastqs = fastqs
    }

    scatter(fastqs_pair in [fastq_splitting.none_fastqs, fastq_splitting.single_fastqs, fastq_splitting.paired_fastqs]) {
        call mapping { input :
            fastqs = paired_fastqs
        }
    }

    call gather_bams { input: 
        bams = mapping.bam
    }

    call pairs_creation { input: 
        bam = gather_bams.linker_paired_bam
    }

    call loops_clustering { input :
        pairs_file = pairs_creation.linker_paired_bam
    }

    # potentially here will come intermediate task of converting bam into pairs file
    # otherwise it has to be repeated in .hic and .cool creation

    call hic_creation { input :
        pairs_file = pairs_creation.pairs_file
    }

    call cooler_creation { input :
        pairs_file = pairs_creation.pairs_file
    }



    # for peak calling use the merged bam file
    call peak_calling { input :
        pairs_file = pairs_creation.merged_bam

    }
}

# workflow tasks 

task fastq_splitting {
	# Inputs
	Array[File] fastqs 	# [end_id]
    run # the sequencing run ID
    
    # Initialize log file
    command {
        log_file = 1.${run}.filter_linker.log
    }
    
    
    ## Filter linker
    # Reads in the read pairs from separate R1 and R2 files
    # Partitions the read pairs into different categories:
    # 1. no linker, 2. linker - one tag, 3. linker - two tags,
    # (4. conflict), (5. tied)
    # For each category, writes out a FASTQ file containing both
    # R1 and R2 reads
    
    command {
        cpu/cpu stag -W -T 18 -t ${n_thread} -O ${run} \
            ${data_dir}/${r1_fastq} ${data_dir}/${r2_fastq} \
            2>> ${log_file}
    }
    
    # Get linker statistics and write to file
    command {
        cpu/cpu stat -s -p -T 18 -t ${n_thread} ${run}.cpu \
            2>> ${log_file} 1> ${run}.stat
    }
    
    
    # Compress the partitioned FASTQ files
    command {
        # Three standard categories
        pigz -p ${n_thread} ${run}.singlelinker.paired.fastq 2>> ${log_file}
        pigz -p ${n_thread} ${run}.none.fastq 2>> ${log_file}
        pigz -p ${n_thread} ${run}.singlelinker.single.fastq 2>> ${log_file}
        
        # Two accessory categories
        pigz -p ${n_thread} ${run}.conflict.fastq 2>> ${log_file}
        pigz -p ${n_thread} ${run}.tied.fastq 2>> ${log_file}
    }
    
    output {
		# After linker filtering, the output FASTQ file for each
		# category contains both R1 and R2 in the same file 
		
		# Core output FASTQ files
		File none_fastq = ${run}.none.fastq.gz
        File single_fastq = ${run}.singlelinker.single.fastq.gz
        File paired_fastq = ${run}.singlelinker.paired.fastq.gz
        
        # Accessory output FASTQ files
        File conflict_fastq = ${run}.conflict.fastq.gz
        File tied_fastq = ${run}.tied.fastq.gz
        
        # Log and statistics files
		File fastq_splitting_log = 1.${run}.filter_linker.log 
		File fastq_splitting_stat = ${run}.stat
	}

    runtime {
        docker: 'some repo docker for CPU'
    }
}

task mapping {
	File idx_tar 		# reference index tar
	Array[File] fastqs 	# [end_id]
	
	command {

		${bin_dir}/cpu/cpu memaln -T ${map_qual} -t ${n_thread} ${fasta} \
            ${run}.${tag_name}.fastq.gz 1> ${run}.${tag_name}.sam 2>> ${log_file}

	}
	output {
		File bam = glob("*.bam")[0]
		File bai = glob("*.bai")[0]
		File align_log = glob("*.align.log")[0]
		File align_qc = glob("*.align.qc")[0]
	}
	runtime {
        docker: 'some repo docker for CPU'
	}
}

task hic_creation {
    command {

    }
    output {
        File hic_file = glob("*.hic")[0]
    }
}

task cooler_creation {
    command {

    }
    output {
        File cooler_file = glob("*.cool")[0]
    }
}

task tags_clustering {
    command {

    }
    output {
        File clustered_loops = glob("*.bedpe")[0]
    }
}

task peak_calling {
    command {

    }
    output {
        File peaks = glob("*.bed")[0]
    }
}

task coverage {
    command {

    }
    output {
        File coverge = glob("*.bedGraph")[0]
    }
}
# utilities

task gather_bams {
	Array[File] bams

	command <<<
		python <<CODE
		
        # python code that will stdout the bam files in way allowing
        # next steps:
        # we need the linker.paired.bam for loops, hic and cooler
        # we need none+linker.single+linke.paired for peaks and coverage
        # the file names can go into two different txt files that could be read
		with open('single.paired.txt','w') as fp:
		    fp.write("what ever is needed to be filled here")
  		with open('three_bams_merged.txt','w') as fp:
		    fp.write("what ever is needed to be filled here")          

        CODE
	>>>
	output {
		String linker_paired_bam = read_string("single.paired.txt")
        String merged_bam = read_string("three_bams_merged.txt")
	}
}

task pairs_creation {
    File bam_file

    command {
        juicer_shortform2pairs.pl bam_file
    }

    output {
        File pairs_file = glob("*.bedGraph")[0]
    }
}


# Question:
# (1) What are file formats are input for hic and cool creation commands?
# (2) Is FASTQ splitting highly parallelizable step? Why it needs 20 threads?
# (3) Is FASTQ mapping highly parallelizable step? Why it needs 20 threads?
# (4) in 2.map.pbs there is pairs creation, clustering and deduplication - I am not sure where is the paired file created? 
# (5)
#
#
#