# module load gcc/7.1.0  openmpi/3.1.4
# module load R/4.1.1
# use local R in mac
# install.packages("glmnet", repos = "https://cran.us.r-project.org")

library(glmnet)
library(Matrix)

# == step1, determine lambda
dir.create("f2_figs")
cancertypes <- c('BRCA','COAD')
for (cancertype in cancertypes){
	x <- read.csv(paste("f1_data/",cancertype,"_atacseq_sig.csv",sep=""), sep=",",row.names=1) 
	x <- as.matrix(x)
	y <- read.csv(paste("f1_data/",cancertype,"_clinical.csv",sep=""), sep=",",row.names=1) 
	y <- as.matrix(y)

	# cv for the best lambda/s
	for (cvi in 1:10){
		pdf(file=paste("f2_figs/",cancertype,"_cvfit_",cvi,".pdf",sep=""))
		cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C",nfolds=5)
		plot(cvfit)
		dev.off()
		}

	# fit glmnet
	fit <- glmnet(x, y, family = "cox")
	pdf(file=paste("f2_figs/",cancertype,"_fit.pdf",sep=""))
	plot(fit,xvar='lambda')
	dev.off()

	# == step2, after get the cv, use s = e^(-2)
	for (s_value in c(0.01,0.03,0.05,0.07,0.09,0.11,0.13)){
		coef <- coef(fit,s=s_value)
		coef <- coef[coef[,1]!=0,]
		coef <- coef[order(coef,decreasing=TRUE)]
		write.table(as.matrix(coef),paste('f2_figs/',cancertype,'_coef_',s_value,'.csv',sep=''),col.names=FALSE,sep=",")
		}
	}


