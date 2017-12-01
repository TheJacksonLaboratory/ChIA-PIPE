
workflow chia-pet {
    # inputs
    # assuming we have a mandatory configuration file:
    File config_file
    Array[File] fastqs 	# [end_id]


    call split_fastqs { input :
        fastqs = fastqs
    }

    call demultiplex { input :
        split_fastqs = split_fastqs.split_fastqs

    }

    scatter(fastqs_pair in demultiplex.split_fastqs) {
        call mapping { input :
            fastqs = fastq_file
        }
    }

    call gather_bams { input: 
        bams = mapping.bam
    }

    call pairs_creation { input: 
        bam = gather_bams.inker_paired_bam
    }

    call hic_creation { input :
        pairs_file = pairs_creation.pairs_file
    }

}

# workflow tasks 

task split_fastqs {
	Array[File] fastqs 	# [end_id]
    String linker
    
    command {

    }

    output {
        # WDL glob() globs in an alphabetical order
		# so R1 and R2 can be switched, which results in an
		# unexpected behavior of a workflow
		# so we prepend <PREFIX>_fastqs_'end'_ (R1 or R2)
        # the <PREFIX> is demultiplexed by special task (demultiplex)
		# to the basename of original filename
		# this prefix should be later stripped in a later task
		Array[File] split_fastqs = glob("*fastqs_R?_*.fastq.gz")
		File fastq_splitting_log = glob("*.splitting.log")[0]
		File fastq_splitting_qc = glob("*.fastq_splitting.qc")[0]
	}

    runtime {
        docker: 'some repo docker for CPU'
    }
}

task mapping {
	File idx_tar 		# reference index tar
	Array[File] fastqs 	# [end_id]
	
	command {
		python $(which encode_bowtie2.py) \
			${idx_tar} \
			${sep=' ' fastqs} \
			${if paired_end then "--paired-end" else ""} \
			${"--multimapping " + multimapping} \
			${"--score-min " + score_min} \
			${"--nth " + select_first([cpu,4])}
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

task demultiplex {
	Array[File] split_fastqs

	command <<<
		python <<CODE
		
        # python code that will stdout the pairs of splitted fastqs
        # none pair
        # linker.single pair
        # linker.paired pair
		
        CODE
	>>>
	output {
		# pair of true replicates
		Array[Array[File]] pairs = read_tsv(stdout())
	}
}