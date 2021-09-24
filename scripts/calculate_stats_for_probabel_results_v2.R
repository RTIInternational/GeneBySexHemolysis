#!/share/apps/R/bin/Rscript

#### old path from ESN cluster (I think? -- EJE)
####/usr/local/bin/R-2.15.0/bin/Rscript

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

removeMissingP = FALSE

args <- commandArgs(TRUE)

loop = TRUE
while (loop) {
	if (args[1] == "--in_file") {
		fileIn = args[2]
	}

	if (args[1] == "--out_file") {
		fileOut = args[2]
	}

	if (args[1] == "--remove_missing_p") {
		removeMissingP = TRUE
	}

	if (length(args)>1) {
		args<-args[2:length(args)]
	} else {
		loop = FALSE
	}
}

cat("Reading ", fileIn, "...\n", sep="")
results = read.table(fileIn, header = TRUE, na.strings = "nan")
#check if file loaded correctly
#x<-nrow(results)
#print(paste0("nrow = ",x))
#x<-colnames(results)
#print(paste0("colnames = ",x))


cat("Calculating stats for SNP...\n")
results$chi_SNP = (results$beta_SNP_add/results$sebeta_SNP_add)^2
results$p_SNP = pchisq(results$chi_SNP,df=1,lower=F)
results$chi_INT = (results$beta_SNP_Gender/results$sebeta_SNP_Gender)^2
results$p_INT = pchisq(results$chi_INT,df=1,lower=F)
if (removeMissingP) {
	results = results[complete.cases(results$p_SNP),]
	results = results[complete.cases(results$p_INT),]
}


cat("Writing ", fileOut, "...\n", sep="")
write.table(results, file = fileOut, row.names = FALSE, quote=FALSE)

cat("Done\n")
