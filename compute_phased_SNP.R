##########################################
# General settings and parameters passing
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)
PILEUP <- args[1]
SNP <- args[2]
padj <- args[3]
cutoff <- args[4]


##########################################
# fetch number of nuclitides for each location
fetch_number <- function(dat) {
	chrom = as.character(dat[1])
	pos = dat[2]
	ref = as.character(toupper(dat[3]))
	num = as.numeric(dat[4])
	bases = as.character(toupper(dat[5]))
	qual = dat[6]

	type <- c(A = 0, T = 0, C = 0, G = 0)

	i = 1
	while (i <= nchar(bases)) {
		base = substr(bases, i, i)
		# remove '^'(start of a read) and '$'(end of a read)
		if (base == "^" | base == "$") {
			i <- i + 1
		# ignore '-'(deletion from the next base)
		} else if (base == "-") {
			next_base = as.numeric(substr(bases, i+1, i+1))
			if (is.na(next_base)) {
				i = i + 1
			} else {
				i = i + 2 + next_base
			}
			#i = i + 2 + next_base
		# ignore '*'(deletion in the base)
		} else if (base == "*") {
			i <- i + 1
			num <- num - 1
		# ignore '+' (insertion)
		} else if (base == "+") {
			next_base = as.numeric(substr(bases, i+1, i+1))
			if (is.na(next_base)) {
				i = i + 1
			} else {
				i = i + 2 + next_base
			}
			#i = i + 2 + next_base
		# ignore 'N'
		} else if (base == "N") {
			i <- i + 1
		# count number of ref allele
		} else if (base == "." | base == ",") {
			type[ref] <- type[ref] + 1
			i <- i + 1
		# count number of alt allele
		} else if (base %in% names(type)) {
			type[base] <- type[base] + 1
			i <- i + 1
		} else {
			i <- i + 1;
		}
	}
	return(c(chrom, pos, ref, num, type))
}



##########################################
# calculate phased SNP
calculate_bias <- function(count, SNP, qvalue_cutoff = 0.1, count_cutoff = 10) {
	colnames(SNP) <- c("chrom", "pos_minus", "pos", "allele")
	count[[2]] <- as.numeric(count[[2]])
	colnames(count) <- c("chrom", "pos", "ref", "num", "A", "T", "C", "G")
	SNP_count <- merge(count, SNP, by.x = c("chrom", "pos"), all.x = T)

	split <- strsplit(as.character(SNP_count$allele),'|') 
	SNP_count <- data.frame(SNP_count, do.call(rbind, split))

	SNP_count$p_allele <- apply(SNP_count, 1, function(x){p_allele <- as.character(x['X1']); as.numeric(as.character(x[p_allele]))})
	SNP_count$m_allele <- apply(SNP_count, 1, function(x){m_allele <- as.character(x['X3']); as.numeric(as.character(x[m_allele]))})

	SNP_count$pvalue <- apply(SNP_count, 1, function(x){p_allele <- x['X1']; m_allele <- x['X3']; binom.test(x = c(as.numeric(as.character(x[p_allele])), as.numeric(as.character(x[m_allele]))), p = 0.5, alternative = "two.sided")$p.value})
	SNP_count$padj <- p.adjust(SNP_count$pvalue, method = "hochberg")
	# SNP_count$biased <- apply(SNP_count, 1, function(x){padj <- as.numeric(x['padj']); count <- as.numeric(x['num']); if(padj < qvalue_cutoff & count > count_cutoff){"Yes"}else{"No"}})

	return(SNP_count)
}



##########################################
# convert column data types in a data frame
convert_data_type <- function(obj,types){
    out <- lapply(1:length(obj),FUN = function(i){FUN1 <- switch(types[i],character = as.character,numeric = as.numeric,factor = as.factor); FUN1(obj[,i])})
    names(out) <- colnames(obj)
    as.data.frame(out,stringsAsFactors = FALSE)
}



##########################################
# start to fetch number for each type of base and do binomial test
PILEUP <- read.table(PILEUP, header = F, comment.char = "")
PILEUP <- PILEUP[PILEUP$V4 != 0, ]
SNP <- read.table(SNP, header = F)

total_res <- t(apply(PILEUP, 1,  fetch_number))
temp <- data.frame(total_res)
total_res_2 <- convert_data_type(temp, rep("character", 8))
total_res_2 <- convert_data_type(total_res_2, c("character", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric"))

final_count <- calculate_bias(total_res_2, SNP)

if(padj) {
 final_count$biased <- apply(final_count, 1, function(x){padj <- as.numeric(x['padj']); if(padj < cutoff){"Yes"}else{"No"}})
final_count <- final_count[, c("chrom", "pos", "allele", "p_allele", "m_allele", "padj", "biased")]
}else {
final_count$biased <- apply(final_count, 1, function(x){pvalue <- as.numeric(x['pvalue']); if(pvalue < cutoff){"Yes"}else{"No"}})
final_count <- final_count[, c("chrom", "pos", "allele", "p_allele", "m_allele", "pvalue", "biased")]
}
#final_count <- final_count[, c("chrom", "pos", "allele", "p_allele", "m_allele", "padj", "biased")]



##########################################
# output for calculating phased loops
write.table(final_count, "phased_snp_qvalue.txt", quote = F, sep = "\t", col.names = T, row.names = F)
