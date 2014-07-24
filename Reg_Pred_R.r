#Prep data
#Spp.dat.otos<-read.table('clipboard',header=T)
#save(Spp.dat.otos,file="C://Users//david.lin//documents-export-2014-07-02//Spp_dat_otos.DMP")
load("C://Users//david.lin//documents-export-2014-07-02//Spp_dat_otos.DMP")
#Subset species data
Hake.dat.otos<-subset(Spp.dat.otos,Species=="Hake")
Petrale.dat.otos<-subset(Spp.dat.otos,Species=="Petrale")
Petrale.dat.otos$Sex<-as.factor(toupper(Petrale.dat.otos$Sex))
Sablefish.dat.otos<-subset(Spp.dat.otos,Species=="Sablefish")
Splitnose.dat.otos<-subset(Spp.dat.otos,Species=="Splitnose")
Splitnose.dat.otos$Sex<-as.factor(toupper(Splitnose.dat.otos$Sex))

##########################
### RUN and FIT models ###
##########################

############
### Hake ###
############
#Check summary stats of otolith weight at age
Hake.Bps.dat<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(10,300),steppin=10,breakpts)$Data
par(mfrow=c(2,2))
hake.summ.All<-boxplot(Hake.Bps.dat$A[,1]~Hake.Bps.dat$A[,2])
hake.summ.F<-boxplot(Hake.Bps.dat$F[,1]~Hake.Bps.dat$F[,2])
hake.summ.M<-boxplot(Hake.Bps.dat$M[,1]~Hake.Bps.dat$M[,2])
#Calculate breakpoints
Hake.Bps.A<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(80,90),steppin=0.1,lowbreaks=c(50,90,50),jitter=100)
Hake.Bps.F<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(91,91),steppin=1,lowbreaks=c(50,90,50))
Hake.Bps.M<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(130,140),steppin=0.1,lowbreaks=c(50,90,50))
Hake.Lms<-Oto.Age.Model.fits(Hake.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Hake.Bps.A$Bps$All$par[1],Hake.Bps.F$Bps$Females$par[1],Hake.Bps.M$Bps$Males$par[1]),lowbreaks=c(50,90,50))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Hake.Lms,1)
OtoAge.plot(Hake.Lms,2)
OtoAge.plot(Hake.Lms,3)
#Hake <300mg -- use this one for now
#Check summary stats of otolith weight at age
Hake.dat.otos.300<-subset(Hake.dat.otos,OtoWt<300)
Hake.Bps.10_300.dat<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(10,300),steppin=10)$Data
par(mfrow=c(2,2))
hake300.summ.All<-boxplot(Hake.Bps.10_300.dat$A[,1]~Hake.Bps.10_300.dat$A[,2])
hake300.summ.F<-boxplot(Hake.Bps.10_300.dat$F[,1]~Hake.Bps.10_300.dat$F[,2])
hake300.summ.M<-boxplot(Hake.Bps.10_300.dat$M[,1]~Hake.Bps.10_300.dat$M[,2])
#Calculate breakpoints
Hake.Bps.A.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(135,140),steppin=0.01,lowbreaks=c(50,90,50))
Hake.Bps.F.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(120,130),steppin=0.1,lowbreaks=c(50,90,50))
Hake.Bps.M.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(135,140),steppin=0.01,lowbreaks=c(50,90,50))
Hake.Lms.300<-Oto.Age.Model.fits(Hake.dat.otos.300,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Hake.Bps.A.300$Bps$All$par[1],Hake.Bps.F.300$Bps$Females$par[1],Hake.Bps.M.300$Bps$Males$par[1]),lowbreaks=c(50,90,50))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Hake.Lms.300,1)
OtoAge.plot(Hake.Lms.300,2)
OtoAge.plot(Hake.Lms.300,3)
############

###############
### Petrale ###
###############
#Check summary stats of otolith weight at age
Petrale.Bps.dat<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,100),steppin=10)$Data
par(mfrow=c(2,2))
petrale.summ.All<-boxplot(Petrale.Bps.dat$A[,1]~Petrale.Bps.dat$A[,2])
petrale.summ.F<-boxplot(Petrale.Bps.dat$F[,1]~Petrale.Bps.dat$F[,2])
petrale.summ.M<-boxplot(Petrale.Bps.dat$M[,1]~Petrale.Bps.dat$M[,2])
#Calculate breakpoints
Petrale.Bps.A<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(52,55),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.F<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(20,24),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.M<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(22,24),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Lms<-Oto.Age.Model.fits(Petrale.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Petrale.Bps.A$Bps$All$par[1],Petrale.Bps.F$Bps$Females$par[1],Petrale.Bps.M$Bps$Males$par[1]),lowbreaks=c(0,0,0))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Petrale.Lms,1)
OtoAge.plot(Petrale.Lms,2)
OtoAge.plot(Petrale.Lms,3)
#Petrale <80mg -- use this one for now
#Check summary stats of otolith weight at age
Petrale.dat.otos.80<-subset(Petrale.dat.otos,OtoWt<80)
Petrale.Bps.10_80.dat<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,100),steppin=10)$Data
petrale80.summ.All<-boxplot(Petrale.Bps.10_80.dat$A[,1]~Petrale.Bps.10_80.dat$A[,2])
petrale80.summ.F<-boxplot(Petrale.Bps.10_80.dat$F[,1]~Petrale.Bps.10_80.dat$F[,2])
petrale80.summ.M<-boxplot(Petrale.Bps.10_80.dat$M[,1]~Petrale.Bps.10_80.dat$M[,2])
#Calculate breakpoints
Petrale.Bps.80.A<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(20,26),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.80.F<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(24,30),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Bps.80.M<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(20,26),steppin=0.1,lowbreaks=c(0,0,0))
Petrale.Lms.80<-Oto.Age.Model.fits(Petrale.dat.otos.80,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(0,100),steppin=10,breakpts=c(Petrale.Bps.80.A$Bps$All$par[1],Petrale.Bps.80.F$Bps$Females$par[1],Petrale.Bps.80.M$Bps$Males$par[1]),lowbreaks=c(0,0,0))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Petrale.Lms.80,1)
OtoAge.plot(Petrale.Lms.80,2)
OtoAge.plot(Petrale.Lms.80,3)
###############

#################
### Sablefish ###
#################
#Check summary stats of otolith weight at age
Sablefish.Bps.dat<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,100),steppin=10)$Data
par(mfrow=c(2,2))
sabfish.summ.All<-boxplot(Sablefish.Bps.dat$A[,1]~Sablefish.Bps.dat$A[,2])
sabfish.summ.F<-boxplot(Sablefish.Bps.dat$F[,1]~Sablefish.Bps.dat$F[,2])
sabfish.summ.M<-boxplot(Sablefish.Bps.dat$M[,1]~Sablefish.Bps.dat$M[,2])
#Calculate breakpoints
Sablefish.Bps.A<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(22,24),steppin=0.01,lowbreaks=c(12,12,12))
Sablefish.Bps.F<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(22,24),steppin=0.01,lowbreaks=c(12,12,12))
Sablefish.Bps.M<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(21,23),steppin=0.01,lowbreaks=c(12,12,12))
#Fit linear models
Sablefish.NWFSC.Lms<-Oto.Age.Model.fits(Sablefish.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(15,25),steppin=0.1,breakpts=c(Sablefish.Bps.A$Bps$All$par[1],Sablefish.Bps.F$Bps$Females$par[1],Sablefish.Bps.M$Bps$Males$par[1]),lowbreaks=c(12,12,12))
par(mfrow=c(2,2))
OtoAge.plot(Sablefish.NWFSC.Lms,1)
OtoAge.plot(Sablefish.NWFSC.Lms,2)
OtoAge.plot(Sablefish.NWFSC.Lms,3)
#################

#################
### Splitnose ###
#################
#Check summary stats of otolith weight at age
Splitnose.Bps.dat<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(0,1000),steppin=50,breakpts)$Data
par(mfrow=c(2,2))
splitnose.summ.All<-boxplot(Splitnose.Bps.dat$A[,1]~Splitnose.Bps.dat$A[,2])
splitnose.summ.F<-boxplot(Splitnose.Bps.dat$F[,1]~Splitnose.Bps.dat$F[,2])
splitnose.summ.M<-boxplot(Splitnose.Bps.dat$M[,1]~Splitnose.Bps.dat$M[,2])
#Calculate breakpoints
Splitnose.Bps.A<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(109,111),steppin=0.01,lowbreaks=c(45,45,45))
Splitnose.Bps.F<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(127,129),steppin=0.01,lowbreaks=c(45,45,45))
Splitnose.Bps.M<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="T",rngSplit=c(104,106),steppin=0.01,lowbreaks=c(45,45,45))
Splitnose.Lms<-Oto.Age.Model.fits(Splitnose.dat.otos,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(15,25),steppin=0.1,breakpts=c(Splitnose.Bps.A$Bps$All$par[1],Splitnose.Bps.F$Bps$Females$par[1],Splitnose.Bps.M$Bps$Males$par[1]),lowbreaks=c(45,45,45))
#Fit linear models
par(mfrow=c(2,2))
OtoAge.plot(Splitnose.Lms,1)
OtoAge.plot(Splitnose.Lms,2)
OtoAge.plot(Splitnose.Lms,3)
#################


Spp.Lms<-list(Hake.Lms,Hake.Lms.300,Petrale.Lms,Petrale.Lms.80,Sablefish.NWFSC.Lms,Splitnose.Lms)
names(Spp.Lms)<-c("Hake","Hake300","Petrale","Petrale80","Sablefish","Splitnose")
save(Spp.Lms,file="C://Users//david.lin//documents-export-2014-07-02//Spp_Lms.DMP")
#Get Numbers
Spp.Ns<-Spp.samp.size(Spp.Lms)
N.Hake300.F<-hist(Hake.Lms.300$Data$Females[,2],breaks=c(0,3,6,9,12,20))
N.Hake300.M<-hist(Hake.Lms.300$Data$Males[,2],breaks=c(0,3,6,9,12,20))
N.Petrale80.F<-hist(Petrale.Lms.80$Data$Females[,2],breaks=c(0,3,6,9,12,20))
N.Petrale80.M<-hist(Petrale.Lms.80$Data$Males[,2],breaks=c(0,3,6,9,12,20))
N.Sablefish.F<-hist(Sablefish.NWFSC.Lms$Data$Females[,2],breaks=c(0,10,20,30,50,100))
N.Sablefish.M<-hist(Sablefish.NWFSC.Lms$Data$Males[,2],breaks=c(0,10,20,30,50,100))
N.Splitnose.F<-hist(Splitnose.Lms$Data$Females[,2],breaks=c(0,5,20,40,70,100))
N.Splitnose.M<-hist(Splitnose.Lms$Data$Males[,2],breaks=c(0,5,20,40,70,105))
#Get R^2s
Spp.R2<-Spp.r2(Spp.Lms)
#Get predicted values
spp.low.ages.mat<-list(
matrix(c(50,90,1,2),nrow=2,ncol=2),matrix(c(50,1),nrow=1,ncol=2),
matrix(c(0,1),nrow=1,ncol=2),matrix(c(0,1),nrow=1,ncol=2),
matrix(c(6,12,0,1),nrow=2,ncol=2),matrix(c(6,12,0,1),nrow=2,ncol=2),
matrix(c(25,45,1,2),nrow=2,ncol=2),matrix(c(22,45,1,2),nrow=2,ncol=2))
Spp.age.pred<-Pred.ages.oto(Spp.Lms,spp.low.ages.mat)
hist(Spp.age.pred$Hake$A[,2])

#VBGF parameters
Hake.vbgf<-vbgf.pars(Spp.age.pred$Hake)
Hake300.vbgf<-vbgf.pars(Spp.age.pred$Hake300)
Petrale.vbgf<-vbgf.pars(Spp.age.pred$Petrale)
Petrale80.vbgf<-vbgf.pars(Spp.age.pred$Petrale80)
Sablefish.vbgf<-vbgf.pars(Spp.age.pred$Sablefish)
Splitnose.vbgf<-vbgf.pars(Spp.age.pred$Splitnose)

#Relative errors
Hake.RE.Linf<-Spp.RE(Hake.vbgf,1)
Hake300.RE.Linf<-Spp.RE(Hake300.vbgf,1)
Petrale.RE.Linf<-Spp.RE(Petrale.vbgf,1)
Petrale80.RE.Linf<-Spp.RE(Petrale80.vbgf,1)
Sablefish.RE.Linf<-Spp.RE(Sablefish.vbgf,1)
Splintose.RE.Linf<-Spp.RE(Splitnose.vbgf,1)
Hake.RE.k<-Spp.RE(Hake.vbgf,2)
Hake300.RE.k<-Spp.RE(Hake300.vbgf,2)
Petrale.RE.k<-Spp.RE(Petrale.vbgf,2)
Petrale80.RE.k<-Spp.RE(Petrale80.vbgf,2)
Sablefish.RE.k<-Spp.RE(Sablefish.vbgf,2)
Splintose.RE.k<-Spp.RE(Splitnose.vbgf,2)
Spp.RE.Linf<-cbind(Hake300.RE.Linf,Petrale80.RE.Linf,Sablefish.RE.Linf,Splintose.RE.Linf)
Spp.RE.k<-cbind(Hake300.RE.k,Petrale80.RE.k,Sablefish.RE.k,Splintose.RE.k)

#Plot parameter estimate comps

#VBGF fits

#1:1 plots

#Ageing error
Hake300.F.4AE<-Punt_age_prep(Spp.age.pred$Hake300$F[,c(2,3)])
Hake300.M.4AE<-Punt_age_prep(Spp.age.pred$Hake300$M[,c(2,3)])
Petrale80.F.4AE<-Punt_age_prep(Spp.age.pred$Petrale80$F[,c(2,3)])
Petrale80.M.4AE<-Punt_age_prep(Spp.age.pred$Petrale80$M[,c(2,3)])
Sablefish.F.4AE<-Punt_age_prep(Spp.age.pred$Sablefish$F[,c(2,3)])
Sablefish.M.4AE<-Punt_age_prep(Spp.age.pred$Sablefish$M[,c(2,3)])
Splitnose.F.4AE<-Punt_age_prep(Spp.age.pred$Splitnose$F[,c(2,3)])
Splitnose.M.4AE<-Punt_age_prep(Spp.age.pred$Splitnose$M[,c(2,3)])
Spp.AgeError<-list(Hake300.F.4AE,Hake300.M.4AE,Petrale80.F.4AE,Petrale80.M.4AE,Sablefish.F.4AE,Sablefish.M.4AE,Splitnose.F.4AE,Splitnose.M.4AE)
names(Spp.AgeError)<-c("Hale_F","Hake_M","Petrale_F","Petrale_M","Sabelfish_F","Sablefish_M","Splitnose_F","Splitnose_M")
save(Spp.AgeError,file="C://Users//david.lin//documents-export-2014-07-02//Spp_AgeError.DMP")

