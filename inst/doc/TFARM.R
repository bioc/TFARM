### R code from vignette source 'TFARM.Rnw'

###################################################
### code chunk number 1: style
###################################################

BiocStyle::latex()



###################################################
### code chunk number 2: options
###################################################

options(width=90)


###################################################
### code chunk number 3: preliminaries
###################################################

library(TFARM)
library(GenomicRanges)



###################################################
### code chunk number 4: TFARM.Rnw:76-82
###################################################
# Load and visualize the dataset:

data("MCF7_chr1")
dim(as.data.frame(MCF7_chr1))
head(MCF7_chr1)



###################################################
### code chunk number 5: TFARM.Rnw:136-149
###################################################
# Coming back to the example on the transcription factors of cell line MCF-7,
# in the promotorial regions of chromosome 1.
# Suppose that the user wants to find the most relevant association rules for the
# prediction of the presence of the transcription factor TEAD4 and such that the
# left-hand-side of the rules contains only present transcription factors.
# This means extracting all the association rules with right-hand-side equal to
# {TEAD4=1} setting the parameter type = TRUE; the minimun support and minimum
# confidence thresholds are set, as an example, to 0.005 and 0.62, respectively:

r_TEAD4 <- rulesGen(MCF7_chr1, "TEAD4=1", 0.005, 0.62, TRUE)
dim(r_TEAD4)
head(r_TEAD4)



###################################################
### code chunk number 6: TFARM.Rnw:166-184
###################################################
# Transcription factors present in at least one of the regions in the considered dataset:

c <- names(elementMetadata(MCF7_chr1))
c
lc <- length(c)

names(presAbs(c, r_TEAD4, TRUE))

# Transcription factors present in at least one of the association rules:

p <- presAbs(c, r_TEAD4, TRUE)$pres
p

# Transcription factors absent in all the association rules:

a <- presAbs(c[1:lc], r_TEAD4, TRUE)$abs
a



###################################################
### code chunk number 7: TFARM.Rnw:309-315
###################################################
# To find the subset of rules containing the transcription factor FOSL2:

r_FOSL2 <- rulesTF(TFi  = 'FOSL2=1', rules =  r_TEAD4, verbose = TRUE)
head(r_FOSL2)
dim(r_FOSL2)[1]



###################################################
### code chunk number 8: TFARM.Rnw:319-324
###################################################
# If none of the rules in input to rulesTF contains the given item TFi,
# and verbose = TRUE, a console message warns for an error:

r_CTCF <- rulesTF(TFi = 'CTCF=1', rules = r_TEAD4, verbose = TRUE)



###################################################
### code chunk number 9: TFARM.Rnw:339-343
###################################################
# For example to evaluate FOSL2 importance in the set of rules r_FOSL2:

r_noFOSL2 <- rulesTF0('FOSL2=1', r_FOSL2, r_TEAD4, MCF7_chr1, "TEAD4=1")



###################################################
### code chunk number 10: TFARM.Rnw:346-348
###################################################
head(r_noFOSL2)



###################################################
### code chunk number 11: IComp
###################################################
# Perform the IComp function to compute the Importance Index distribution:

imp_FOSL2 <- IComp('FOSL2=1', r_FOSL2, r_noFOSL2, figures=TRUE)
names(imp_FOSL2)

imp_FOSL2$imp

head(imp_FOSL2$delta)
head(imp_FOSL2$rwi)
head(imp_FOSL2$rwo)



###################################################
### code chunk number 12: TFARM.Rnw:391-413
###################################################
# For the considered example the user could run:

library(plyr)

A <- vector("list", length(p))
B <- vector("list", length(p))
IMP <- matrix(0, length(p), 4)
IMP <- data.frame(IMP)
IMP[,1] <- paste(p)
colnames(IMP) <- c('TF', 'imp', 'sd', 'nrules')
IMP_Z <- vector("list", length(p))

for (i in 1:length(p))  {
	A[[i]] <- rulesTF(p[i], r_TEAD4, FALSE)
	B[[i]] <- rulesTF0(p[i], A[[i]], r_TEAD4, MCF7_chr1, "TEAD4=1")
	IMP_Z[[i]] <- IComp(p[i], A[[i]], B[[i]], figures=FALSE)$imp
}

IMP[, 2] <- sapply(IMP_Z, mean)
IMP[, 3] <- sapply(IMP_Z, sd)
IMP[, 4] <- sapply(IMP_Z, length)



###################################################
### code chunk number 13: TFARM.Rnw:416-422
###################################################

# Sort by imp column of IMP

IMP.ord <- arrange(IMP, desc(imp))
IMP.ord



###################################################
### code chunk number 14: IPCA
###################################################
# Extract the delta variations of support, confidence and lift:

DELTA <- vector("list", length(p))
for (i in 1:length(p)){
DELTA[[i]] <- IComp(p[i], A[[i]], B[[i]], figures=FALSE)$delta
}

# Select the candidate co-regulators and the number of rules associated with them, then
# perform the Principal Component Analysis:

colnames(IMP)
I <- data.frame(IMP$TF, IMP$imp, IMP$nrules)
i.pc <- IPCA(DELTA, I)
names(i.pc)

i.pc$summary

head(i.pc$loadings)

head(i.pc$scores)



###################################################
### code chunk number 15: distribViz
###################################################
# Considering for example the candidate co-regulator transcription factors
# found in the set of rules r_TEAD4:

distribViz(IMP_Z, p)



###################################################
### code chunk number 16: TFARM.Rnw:538-559
###################################################
# Select the index of the list of importances IMP_Z
# containing importance distributions of transcription factor ZNF217
ZNF217_index <- which(p == 'ZNF217=1')

# select outlier rules where ZNF217 has importance greater than 5
o <- IMP_Z[[ZNF217_index]] > 5
rule_o <- A[[ZNF217_index]][o,]
rule_o

# So, ZNF217 is very relevant in the pattern of transcription factors
# {GABPA=1,TCF12=1,ZNF217=1,NR2F2=1}
# for the prediction of the presence of TEAD4.

# To extract support, confidence and lift of the corresponding rule without ZNF217:
B[[ZNF217_index]][o,]

# Since none of the three measures of the rule obtained removing ZNF217 are equal to zero,
# the rule {GABPA=1,TCF12=1,ZNF217=0,NR2F2=1} -> {TEAD4=1},
# obtained removing ZNF217, is found in the relevant rules for the prediction
# of the presence of TEAD4.



###################################################
### code chunk number 17: TFARM.Rnw:591-599
###################################################
# Construct couples as a vector in which all possible combinations of
# transcription factors (present in at least one association rules)
# are included:

couples_0 <- combn(p, 2)
couples <- paste(couples_0[1,], couples_0[2,], sep=',')
head(couples)



###################################################
### code chunk number 18: TFARM.Rnw:602-632
###################################################
# The evaluation of the mean Importance Index of each pair is
# then computed similarly as previously done for single transcription factors:

A_c <- vector("list", length(couples))
B_c <- vector("list", length(couples))
I_c <- matrix(0, length(couples), 2)
I_c <- data.frame(I_c)
I_c[,1] <- paste(couples)
colnames(I_c) <- c('TF', 'imp')
IMP_c <- vector("list", length(couples))

# Compute rulesTF, rulesTF0 and IComp for each pair, avoiding pairs not
# found in the r_TEAD4 set of rules

for (i in 1:length(couples))  {
	A_c[[i]] <- rulesTF(couples[i], r_TEAD4, FALSE)
	if (all(!is.na(A_c[[i]][1]))){
	B_c[[i]] <- rulesTF0(couples[i], A_c[[i]], r_TEAD4, MCF7_chr1, "TEAD4=1")
  	IMP_c[[i]] <- IComp(couples[i], A_c[[i]], B_c[[i]], figures=FALSE)$imp
	}
}

# Delete all NULL elements and compute the mean Importance Index of each pair

null.indexes <- sapply(IMP_c, is.null)
IMP_c <- IMP_c[!null.indexes]
I_c <- I_c[!null.indexes,]

I_c[, 2] <- sapply(IMP_c, mean)



###################################################
### code chunk number 19: TFARM.Rnw:635-641
###################################################
# Select rows with mean Importance Index different from NaN, then order I_c:

I_c <- I_c[!is.na(I_c[,2]),]
I_c_ord <- arrange(I_c, desc(imp))
head(I_c_ord)



###################################################
### code chunk number 20: heatmap
###################################################
# Construction of a vector in which mean Importance Index values of pairs
# of transcription factors are represented.
# These transcription factors are taken from the output of presAbs as
# present in at least one association rules.

# The function rbind is used to combine IMP columns and I_c_ord columns and
# then the function arrange orders the data frame by the imp column.

I_c_2 <- arrange(rbind(IMP[,1:2], I_c_ord), desc(imp))
i.heat <- heatI(p, I_c_2)



