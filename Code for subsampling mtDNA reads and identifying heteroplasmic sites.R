##subsampling mtDNA reads and identifying heteroplasmic sites
##Reena Debray

#input for this script: a mitocaller output file, uploaded and named "mtcalls"
#output of this script: a data frame named "output" with information on the position, genotypes, and allele frequencies of each heteroplasmic site

#initialize data frame of heteroplasmic sites
output<-data.frame(matrix(nrow=0,ncol=3))
options(stringsAsFactors=F) #rbind command later in the script will not work otherwise

for (line in seq(1,nrow(mtcalls))){
#remove reads with base quality < 20
  reads=unlist(strsplit(substr(mtcalls[line,35],6,nchar(as.character(mtcalls[line,35]))),","))
  base_quals=unlist(strsplit(substr(mtcalls[line,36],7,nchar(as.character(mtcalls[line,36]))),","))
  reads_and_bqs=data.frame(reads,base_quals)
  reads_and_bqs_filtered<-reads_and_bqs[as.numeric(as.character(reads_and_bqs$base_quals))>=20,]

#count number of reads; discard if less than 200
  if (nrow(reads_and_bqs_filtered)>=200){
    
    #downsample sites with more than 200 reads to 200 reads
    if (nrow(reads_and_bqs_filtered)>200) {reads_and_bqs_subsampled<-reads_and_bqs_filtered[sample(seq(1,nrow(reads_and_bqs_filtered)),size=200,replace=F),]}
    
    #keep alleles that are supported by 8 or more reads
    supported_alleles=c()
    for (base in c("A","C","G","T")) {
      if (nrow(reads_and_bqs_subsampled[reads_and_bqs_subsampled$reads==base,])>=8) {supported_alleles<-c(supported_alleles,as.character(reads_and_bqs_subsampled[reads_and_bqs_subsampled$reads==base,"reads"]))}
    }
    
    #if there are multiple alleles (heteroplasmy), recalculate frequencies based on supported alleles
    if (length(table(supported_alleles))>1){
      new_freqs<-data.frame(table(supported_alleles))
      new_freqs$percentreads<-round(new_freqs$Freq/sum(new_freqs$Freq),4)
    
      #re-format results and append to data frame of heteroplasmic sites
      genotypes<-paste(new_freqs$supported_alleles,collapse="/")
      percentreads<-paste(new_freqs$percentreads,collapse="/")
      newline=c(line,genotypes,percentreads)
      output<-rbind(output,newline)
    }
  }
}

colnames(output)=c("Site","Genotype","Frequency")
