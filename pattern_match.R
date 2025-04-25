.libPaths("/home/groups/mcweeney_lab/jengs/rpackages/")
library(Biostrings)
library(argparse)
parser <- ArgumentParser(description='Library Sequence Similiarity')
parser$add_argument('--lib_file',nargs='*',
                   help='library file')
parser$add_argument('--fasta_file',nargs='*',
                   help='fasta file')
parser$add_argument('--num_mismatch',nargs='*',
                   help='maximum number of mismatch allowed')
parser$add_argument('--lib_name',nargs='*',
                   help='library name')

args <- parser$parse_args(commandArgs(TRUE))
print(args$lib_file)
print(args$fasta_file)
print(args$num_mismatch)
print(args$lib_name)

#setwd("/home/exacloud/gscratch/mcweeney_lab/resources/u54/crispr/libraries")
#library_file<-read.csv("yusa_library.csv",stringsAsFactors=F,header=F)
library_file<-read.csv(args$lib_file,stringsAsFactors=F,header=F)

colnames(library_file)<-c("sgRNA","Sequence","Gene")
#library_file<-read.csv("mousa_library.csv",stringsAsFactors=F,header=F)
#library_file<-read.csv("brunello_library1.csv",stringsAsFactors=F,header=F)

sequence<-library_file[,"Sequence"]
#sequence<-sequence[1:10]
genome<-readDNAStringSet(args$fasta_file)
seqPatternMatch<-lapply(sequence,function(x,gen) {
        patternCount<-vcountPattern(x,gen,max.mismatch=as.numeric(as.character(args$num_mismatch)))
        return(length(which(patternCount>0))-1)
	#return(patternCount)
},genome)

names(seqPatternMatch)<-sequence
seqPatternMatch<-do.call(rbind,seqPatternMatch)
seqPatternMatch<-cbind(sequence,seqPatternMatch)
colnames(seqPatternMatch)<-c("Sequence",paste0("num_sequence_le_",args$num_mismatch,"_mismatch"))
seqPatternMatch<-merge(library_file,seqPatternMatch,by.x="Sequence",by.y="Sequence")
setwd("/home/groups/mcweeney_lab/jengs/crispr_library_scores")
write.table(seqPatternMatch,file=paste0(args$lib_name,"_sequence_sim_le_",args$num_mismatch,".txt"),sep="\t",col.names=T,row.names=F)
#save(seqPatternMatch,file="yusaPatternMatchMax1.rda")

