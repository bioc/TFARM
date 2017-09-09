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
### code chunk number 4: TFARM.Rnw:75-81
###################################################
# Load and visualize the dataset:

data("MCF7_chr1")
length(MCF7_chr1)
MCF7_chr1



###################################################
### code chunk number 5: TFARM.Rnw:135-148
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
### code chunk number 6: TFARM.Rnw:165-183
###################################################
# Transcription factors present in at least one of the regions in the considered dataset:

c <- names(mcols(MCF7_chr1))
c
lc <- length(c)

names(presAbs(c, r_TEAD4, TRUE))

# Transcription factors present in at least one of the association rules:

p_TFs <- presAbs(c, r_TEAD4, TRUE)$pres
p_TFs

# Transcription factors absent in all the association rules:

a <- presAbs(c[1:lc], r_TEAD4, TRUE)$abs
a



###################################################
### code chunk number 7: TFARM.Rnw:308-314
###################################################
# To find the subset of rules containing the transcription factor FOSL2:

r_FOSL2 <- rulesTF(TFi  = 'FOSL2=1', rules =  r_TEAD4, verbose = TRUE)
head(r_FOSL2)
dim(r_FOSL2)[1]



###################################################
### code chunk number 8: TFARM.Rnw:318-323
###################################################
# If none of the rules in input to rulesTF contains the given item TFi,
# and verbose = TRUE, a warnig message is reported to the user:

r_CTCF <- rulesTF(TFi = 'CTCF=1', rules = r_TEAD4, verbose = TRUE)



###################################################
### code chunk number 9: TFARM.Rnw:338-342
###################################################
# For example to evaluate FOSL2 importance in the set of rules r_FOSL2:

r_noFOSL2 <- rulesTF0('FOSL2=1', r_FOSL2, r_TEAD4, MCF7_chr1, "TEAD4=1")



###################################################
### code chunk number 10: TFARM.Rnw:345-347
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
### code chunk number 12: TFARM.Rnw:390-412
###################################################
# For the considered example the user could run:

all <- lapply(p_TFs, function(pi) {
  	A <- rulesTF(pi, r_TEAD4, FALSE)
  	B <- rulesTF0(pi, A, r_TEAD4, MCF7_chr1, "TEAD4=1")
    IComp(pi, A, B, figures=FALSE)
  })

for (i in 1:length(p_TFs)) {
    IMP_Z[[i]] <- all[[i]]$imp
# Extract the delta variations of support, confidence and lift:
    DELTA[[i]] <- all[[i]]$delta
}

  IMP <- data.frame(
      TF = p_TFs,
      imp = sapply(IMP_Z, mean),
      sd = sapply(IMP_Z, sd),
      nrules = sapply(IMP_Z, length),
      stringsAsFactors=FALSE
  )



###################################################
### code chunk number 13: TFARM.Rnw:415-422
###################################################

# Sort by imp column of IMP

library(plyr)
IMP.ord <- arrange(IMP, desc(imp))
IMP.ord



###################################################
### code chunk number 14: IPCA
###################################################
# Select the candidate co-regulators and the number of rules associated with them, then
# perform the Principal Component Analysis:

colnames(IMP)
TF_Imp <- data.frame(IMP$TF, IMP$imp, IMP$nrules)
i.pc <- IPCA(DELTA, TF_Imp)
names(i.pc)

i.pc$summary

head(i.pc$loadings)

head(i.pc$scores)



###################################################
### code chunk number 15: distribViz
###################################################
# Considering for example the candidate co-regulator transcription factors
# found in the set of rules r_TEAD4:

distribViz(IMP_Z, p_TFs)



###################################################
### code chunk number 16: TFARM.Rnw:531-552
###################################################
# Select the index of the list of importances IMP_Z
# containing importance distributions of transcription factor ZNF217
ZNF217_index <- which(p_TFs == 'ZNF217=1')

# select outlier rules where ZNF217 has importance greater than 5
o <- IMP_Z[[ZNF217_index]] > 5
rule_o <- all[[ZNF217_index]]$rwi[o,]
rule_o

# So, ZNF217 is very relevant in the pattern of transcription factors
# {FOSL2=1,GABPA=1,MYC=1,ZNF217=1}
# for the prediction of the presence of TEAD4.

# To extract support, confidence and lift of the corresponding rule without ZNF217:
all[[ZNF217_index]]$rwo[o,]

# Since the measure of the rule obtained removing ZNF217 is equal to zero,
# the rule {FOSL2=1,GABPA=1,MYC=1,ZNF217=0} -> {TEAD4=1},
# obtained removing ZNF217, is found in the relevant rules for the prediction
# of the presence of TEAD4.



###################################################
### code chunk number 17: TFARM.Rnw:584-592
###################################################
# Construct couples as a vector in which all possible combinations of
# transcription factors (present in at least one association rules)
# are included:

couples_0 <- combn(p_TFs, 2)
couples <- paste(couples_0[1,], couples_0[2,], sep=',')
head(couples)



###################################################
### code chunk number 18: TFARM.Rnw:595-623
###################################################
# The evaluation of the mean Importance Index of each pair is
# then computed similarly as previously done for single transcription factors:

# Compute rulesTF, rulesTF0 and IComp for each pair, avoiding pairs not
# found in the r_TEAD4 set of rules

IMP_c <- lapply(couples, function(ci) {
  	A_c <- rulesTF(ci, r_TEAD4, FALSE)
  	if (all(!is.na(A_c[[1]][1]))){
	B_c <- rulesTF0(ci, A_c, r_TEAD4, MCF7_chr1, "TEAD4=1")
  	IComp(ci, A_c, B_c, figures=FALSE)$imp
	}
  })


# Delete all NULL elements and compute the mean Importance Index of each pair

I_c <- matrix(0, length(couples), 2)
I_c <- data.frame(I_c)
I_c[,1] <- paste(couples)

null.indexes <- vapply(IMP_c, is.null, numeric(1))
IMP_c <- IMP_c[!null.indexes]
I_c <- I_c[!null.indexes,]

I_c[,2] <- vapply(IMP_c, mean, numeric(1))
colnames(I_c) <- colnames(IMP[,1:2])



###################################################
### code chunk number 19: TFARM.Rnw:626-632
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
i.heat <- heatI(p_TFs, I_c_2)



