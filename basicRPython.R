#TODO set libpath
.libPaths("/home/groups/mcweeney_lab/jengs/rpackages")

library(data.table)

library(argparse)
parser <- ArgumentParser(description='Process alignment results')
parser$add_argument('--stats',nargs='*',
                   help='python stats files')
parser$add_argument('--counts', nargs='*',help='python counts files')
parser$add_argument('--project_name',nargs='*',help='Project name, prefix for outputfiles')
parser$add_argument('--counts_output', nargs='+',help='counts output file')
parser$add_argument('--stats_output', nargs='+',help='stats output file')
parser$add_argument('--annotation', nargs='+',help='library annotation file')
parser$add_argument('--correct_len',nargs='*',
                   help='python correct length seq files')
parser$add_argument('--all_len',nargs='*',
                   help='python all length count files')

parser$add_argument('--correct_length_output', nargs='+',help='correct length output file')
parser$add_argument('--all_length_output', nargs='+',help='all length output file')

args <- parser$parse_args(commandArgs(TRUE))
print(args$stats)
print(args$counts)
print(args$stats_output)
print(args$counts_output)


stats<- args$stats
correct_length<-args$correct_len
all_length<-args$all_len

align.stats <- do.call(rbind, lapply(stats, function(x){
	tmp <- fread(x)
	tmp[,Sample:=gsub("_statistics.txt","",basename(x))]
	return(tmp)
}))
fwrite(align.stats,file=args$stats_output,sep="\t",col.names=T,row.names=F)


correct_len_seq_ctr <- do.call(rbind, lapply(correct_length, function(x){
        tmp <- fread(x)
        tmp[,Sample:=gsub("_trimmed_counts.txt","",basename(x))]
        return(tmp)
}))
fwrite(correct_len_seq_ctr,file=args$correct_length_output,sep="\t",col.names=T,row.names=F)

length_ctr <- do.call(rbind, lapply(all_length, function(x){
        tmp <- fread(x)
        tmp[,Sample:=gsub("_sequence_counts.txt","",basename(x))]
        return(tmp)
}))
fwrite(length_ctr,file=args$all_length_output,sep="\t",col.names=T,row.names=F)



sample_counts<-lapply(args$counts,function(x) {
	tmp<-fread(x,header=F)
	colnames(tmp)<-c("sgRNA",gsub("_counts.txt","",basename(x)))
	return(tmp)
})
print(colnames(sample_counts[[1]]))
sample_counts<-Reduce(function(...) merge.data.table(...,all=T,by="sgRNA"),sample_counts)
lib_annot<-fread(args$annotation)
colnames(lib_annot)[1:3]<-c("Sequence","sgRNA","Gene")
sample_counts<-merge(lib_annot,sample_counts,by="sgRNA",all=T)
fwrite(sample_counts,file=args$counts_output,sep="\t",col.names=T,row.names=F)

