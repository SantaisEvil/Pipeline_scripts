library(optparse)
library(scales)

##Initialize matrix
mt <- matrix(ncol = 9, nrow = 0)

#Parse file dirs from snakemake
files <- Sys.glob(snakemake@params[[1]])
dup_files <- Sys.glob(snakemake@params[[2]])
rmdup_files <- Sys.glob(snakemake@params[[3]])



for (i in 1:length(files)){

  filename <- basename(files[i])
  
  if (filename != "ucsc") {
    sample_name <- unlist(strsplit(filename,"_"))[[1]]

    
    text_file <- read.delim(paste0(files[i],"/bowtie2/unfiltered_aligned.txt") ,sep ="", header = F)


    seq_depth <- text_file[1,1]
    alignment_rate <- text_file[6,1]
    
    mapped_frags <- as.numeric(text_file[1,1]) - as.numeric(text_file[3,1])
    text_file <- read.delim(dup_files[i]  ,sep = "", header = F)
    read_pair_dups <- as.numeric(text_file[7,8]) + as.numeric(text_file[7,9])
    percent_dup <- text_file[7,10]
    percent_dup <- percent(as.numeric(percent_dup),accuracy = 0.01)
    
    ##Stats after dropping non-primary reads.
    text_file <- read.delim(rmdup_files[i]  ,sep = "", header = F)
    surviving <- text_file[7,4]
    
    final_alignment <- percent((as.numeric(surviving)/as.numeric(seq_depth)) ,accuracy = 0.01)
    per_removed <-  percent(1-(as.numeric(surviving)/as.numeric(seq_depth)) ,accuracy = 0.01)
    data <- c(sample_name,seq_depth,mapped_frags,alignment_rate, read_pair_dups,percent_dup, surviving, final_alignment, per_removed)
    
    mt <- rbind(mt,data)
  }
  
}

df <- data.frame(mt)

rownames(df) <- NULL


colnames(df) <- c("Sample Name","Sequencing Depth","Mapped Fragments", "Alignment Rate", "Read Pair Duplicates", "Percent Duplication","Primary Alignments", "Final Primary Alignment of Total Reads","Dropped Reads %")

df$`Sample Name` <- as.factor(df$`Sample Name`)
df$`Sequencing Depth` <- as.numeric(df$`Sequencing Depth`)
df$`Mapped Fragments` <- as.numeric(df$`Mapped Fragments`)


write.csv(df,snakemake@output[[1]],row.names= F)