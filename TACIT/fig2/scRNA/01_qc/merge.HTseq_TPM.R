#!/bin/env Rscript

args = commandArgs(T)

if (length(args) != 4){
    print('USAGE: Rscript merge.R pattern merge_colname input_dir output_file')
    quit()
}


files = list.files(args[3],pattern = paste(args[1],'$',sep=''))

a <- files[1]
atable <- read.table(a, header=F, check.names=F, sep='\t')
atable <- atable[,c(1,3)]#column 1 is gene name; column 3 is TPM
samplea = strsplit(basename(a), '.', fixed=T)[[1]][1]
colnames(atable) = c(args[2], samplea)

for (i in files[-1]){
        print (i)
        sampleb = strsplit(basename(i), '.', fixed=T)[[1]][1]
    btable <- read.table(i, header=F, check.names=F, sep='\t')
    btable <- btable[,c(1,3)]#column 1 is gene name; column 3 is TPM
    colnames(btable) = c(args[2], sampleb)
        print(colnames(btable))
        #names(btable) <- c('gene_id', fpkm)
        atable <- merge(atable, btable, by = args[2], all =TRUE)
}

row.names(atable) = atable[,1]
atable = atable[,-1]
print(colnames(atable))
write.table(atable, sep='\t', file=args[4], col.names=NA,quote=FALSE, na='NA')