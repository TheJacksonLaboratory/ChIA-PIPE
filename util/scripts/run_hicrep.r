bin_dir="/projects/capurd/chia_pet/chia_pet_tool_2"
library(hicrep, lib.loc=bin_dir)
options(scipen = 999)


run_hicrep <-  function(
    run_1, run_2, path_1, path_2, chroms, resolut=100000, h=3, maxd=5000000){
    ### 
    #takes 2 matrix with the same index: idx
    #idx has headers: chrom start end
    #resolut = bin resolution in mat1 & mat2
    #chrs should be in idx$CHR

    scc = NULL
    
    for(chrom in chroms){
        print(chrom)
        # Set file names
        file_1 = paste0(path_1, '/', run_1, '.chr', chrom, '.mat')
        file_2 = paste0(path_2, '/', run_2, '.chr', chrom, '.mat')
        
        # Read matrices
        mat_1 = read.table(file_1)
        mat_2 = read.table(file_2)
        
        # Reformat matrices
        if ( nrow(mat_1) > nrow(mat_2)){
            mat_1 = mat_1[1:nrow(mat_2), 1:nrow(mat_2)]
        } else if ( nrow(mat_2) > nrow(mat_1)){
            mat_2 = mat_2[1:nrow(mat_1), 1:nrow(mat_1)]
        }
        
        n_row = nrow(mat_1)
        chrom_col = rep(chrom, n_row)
        start_col = as.integer(seq(0, (resolut*n_row - resolut), resolut))
        end_col = as.integer(seq(resolut, (resolut*n_row), resolut))
        
        mat_1 = cbind(chrom_col, start_col, end_col, mat_1)
        mat_2 = cbind(chrom_col, start_col, end_col, mat_2)
        
        
        # Pre-process data
        mprep =  prep(mat_1, mat_2, resolut, h, maxd)
        
        # Apply HiCRep
        scc_chrom = get.scc(mprep, resolut, maxd)
        
        # Save results
        scc = c(scc, scc_chrom$scc)
    }
    return(scc)
}


## MAIN ######

# Resolution
resolut=100000

# Information for replicate 1
run_1 = "LHH0048H"
run_2 = "LHH0054H"

# Information for replicate 2
path_1 = paste0("/projects/ruan-lab/capurd/processing/results/LHH0048H/",
"LHH0048H_intrachrom_mats_", resolut)

path_2 = paste0("/projects/ruan-lab/capurd/processing/results/LHH0054H/",
"LHH0054H_intrachrom_mats_", resolut)

# Set chromosomes
#chroms = 1:22
#chroms = c(chroms, "X")
chroms = c("21", "22")

# Run HiCRep
scc = run_hicrep(run_1, run_2, path_1, path_2, chroms, resolut)

score = mean(scc)

write(score, file=paste0('hicrep_', run_1, '_', run_2, '_', resolut, ".txt"), sep='\t')

