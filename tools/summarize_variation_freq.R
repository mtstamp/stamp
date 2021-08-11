#############################################################################
#
#Updated by: Yiqin Wang
#TODO
#to add similiar functionalities in annot
#
##############################################################################

args <- commandArgs(trailingOnly = TRUE)

var <- read.table(args[1], header = T, sep="\t")
var1 <- read.table(args[2], header = T, sep= "\t")
var2 <- read.table(args[3], header = T, sep= "\t")

#sort sample and position
var = var[order(var$sample, var$pos),]
var1 = var1[order(var1$sample, var1$pos),]
var2 = var2[order(var2$sample, var2$pos),]

#print out the number of variants in each file
print(sprintf("Read %d rows from %s", nrow(var), args[1]))
print(sprintf("Found %d matched records in %s", sum(var$pos == var1$pos & var$sample == var1$sample), args[2]))
print(sprintf("Found %d matched records in %s", sum(var$pos == var2$pos & var$sample == var2$sample), args[3]))

allele.count <- as.matrix(cbind(var1[,c("A","T","C","G")],var2[,c("A","T","C","G")]))

n=0
n1 = 0
#workspace <- 2000000
comp.frac <- function(x){
if (sum(is.na(x))>0 | sum(x[1:4]) == 0 | sum(x[5:8]) == 0){
NA
}else{
#if numbers are too large use chisq test instead
tryCatch({n=n+1;fisher.test(matrix(x,c(2,4),byrow=T))$p}, error = function(err){
print("chisq"); chisq.test(matrix(x,c(2,4),byrow=T))$p.val
})
}
}

fp <- apply(allele.count, 1, comp.frac)
print(n1)
#only record depth, allele fraction, and allele fraction difference information
var[,"f1.depth"] = var1$depth_fwd
var[,"f1.het_allele"] = var1$het_allele
var[,"f1.het_freq"] = var1$het_freq
var[,"f2.depth"] = var2$depth_fwd
var[,"f2.het_allele"] = var2$het_allele
var[,"f2.het_freq"] = var2$het_freq
var[,"fp.freq"] = fp

write.table(var, args[4], col.names= T, row.names = F, sep="\t", quote = F)


