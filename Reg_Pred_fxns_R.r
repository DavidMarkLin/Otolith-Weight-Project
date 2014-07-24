#############################################
#############  FUNCTIONS  ###################
#############################################
library(gplots)
#Fit oto weights and ages
#oto.age.col = columns that contain the otolith weight and ages, respectively
#sextype = "All" means both M and F use any unsexed individuals in the data set.
#BP.find = T- find breakpoints. BP.find = F- fit linear models.
#steppin = intervals to search over
#breakpts = breakpoints used when Bp.find = "F"
#lowbreaks = the otolith weight at which the data set is separated. The low numbers are not used in any of the linear models.
#intcpt = "Y" estimates the intercept; "N" sets it to 0.


Oto.Age.Model.fits<-function(spp.dat.in,oto.age.col=c(5,4),sextype="All",Bp.find="F",rngSplit=c(10,200),steppin=10,breakpts,lowbreaks=c(0,0,0),intcpt="Y",CG="F",jitter=0)
{
  #Data prep
  dat.names<-c("All","Females","Males")
  Spp.dat.AFM<-list()
  Spp.dat.AFM[[1]]<-na.omit(cbind(spp.dat.in[,oto.age.col[1]],spp.dat.in[,oto.age.col[2]]))
  if(sextype=="All")
  {
    Spp.dat.AFM[[2]]<-na.omit(cbind(subset(spp.dat.in,Sex!="M")[,oto.age.col[1]],subset(spp.dat.in,Sex!="M")[,oto.age.col[2]]))
    Spp.dat.AFM[[3]]<-na.omit(cbind(subset(spp.dat.in,Sex!="F")[,oto.age.col[1]],subset(spp.dat.in,Sex!="F")[,oto.age.col[2]]))
    names(Spp.dat.AFM)<-dat.names
    Lts<-list()
    if(length(spp.dat.in[,3])!=length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[as.numeric(-attributes(Spp.dat.AFM[[1]])$na.action),3]}
    if(length(spp.dat.in[,3])==length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[,3]}
    if(length(subset(spp.dat.in,Sex!="M")[,3])!=length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex!="M")[as.numeric(-attributes(Spp.dat.AFM[[2]])$na.action),3]}
    if(length(subset(spp.dat.in,Sex!="M")[,3])==length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex!="M")[,3]}
    if(length(subset(spp.dat.in,Sex!="F")[,3])!=length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex!="F")[as.numeric(-attributes(Spp.dat.AFM[[3]])$na.action),3]}
    if(length(subset(spp.dat.in,Sex!="F")[,3])==length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex!="F")[,3]}
  }
  if(sextype!="All")
  {
    Spp.dat.AFM[[2]]<-na.omit(cbind(subset(spp.dat.in,Sex=="F")[,oto.age.col[1]],subset(spp.dat.in,Sex=="F")[,oto.age.col[2]]))
    Spp.dat.AFM[[3]]<-na.omit(cbind(subset(spp.dat.in,Sex=="M")[,oto.age.col[1]],subset(spp.dat.in,Sex=="M")[,oto.age.col[2]]))
    names(Spp.dat.AFM)<-dat.names
    Lts<-list()
    if(length(spp.dat.in[,3])!=length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[as.numeric(-attributes(Spp.dat.AFM[[1]])$na.action),3]}
    if(length(spp.dat.in[,3])==length(Spp.dat.AFM[[1]][,1])){Lts[[1]]<-spp.dat.in[,3]}
    if(length(subset(spp.dat.in,Sex=="F")[,3])!=length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex=="F")[as.numeric(-attributes(Spp.dat.AFM[[2]])$na.action),3]}
    if(length(subset(spp.dat.in,Sex=="F")[,3])==length(Spp.dat.AFM[[2]][,1])){Lts[[2]]<-subset(spp.dat.in,Sex=="F")[,3]}
    if(length(subset(spp.dat.in,Sex=="M")[,3])!=length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex=="M")[as.numeric(-attributes(Spp.dat.AFM[[3]])$na.action),3]}
    if(length(subset(spp.dat.in,Sex=="M")[,3])==length(Spp.dat.AFM[[3]][,1])){Lts[[3]]<-subset(spp.dat.in,Sex=="M")[,3]}
  }
  names(Lts)<-c("Lt_A","Lt_F","Lt_M")


  #Find breakpoint
  if(Bp.find=="T")
  {
    Spp.bps<-Prof.out<-Low.breaks<-list()
    #Create 2x2 empty panels for graphs
    par(mfrow=c(2,2))
    for(i in 1:length(Spp.dat.AFM))
    {
      Spp.dat.AFM.mod<-Spp.dat.AFM[[i]][Spp.dat.AFM[[i]][,1]>lowbreaks[i],]
      if(CG=="F")
      {
        #minimize the function PReg.obj, given initial guesses for breakpoint, slopes, and y-intercept.
        out.f<- optim(c(sum(rngSplit)/2,0,0.1,0.1),PReg.obj,x.in=Spp.dat.AFM.mod[,1],y.in=Spp.dat.AFM.mod[,2],verbose=F,rngSplit=rngSplit)
        #record results of minimization, record breakpoint, standardize name of original and modified data
      }
      if(CG=="T")
      {
        out.f<- optim(c(sum(rngSplit)/2,0,0.1,0.1),PReg.obj,x.in=Spp.dat.AFM.mod[,1],y.in=Spp.dat.AFM.mod[,2],verbose=F,rngSplit=rngSplit,df,method="CG")
      }
      out.f.jitter<-NA
      if(jitter>0)
      {
        out.f.jitter=list()
        yint<-runif(jitter,0,5)
        m1<-rnorm(jitter,0.1,0.03)
        m2<-rnorm(jitter,0.1,0.03)
        for(x in 1:jitter){
        out.f.jitter[[x]]<- optim(c(sum(rngSplit)/2,yint[x],m1[x],m2[x]),PReg.obj,x.in=Spp.dat.AFM.mod[,1],y.in=Spp.dat.AFM.mod[,2],verbose=F,rngSplit=rngSplit)
        }
        inputy<-inputm1<-inputm2<-c()
        for(i in 1:jitter){
          inputy[i] = out.f.jitter$i$par$yint
          inputm1[i] = out.f.jitter$i$par$m1
          inputm2[i] = out.f.jitter$i$par$m2
        }
        x = c(1:100)
        plot(x,inputy)
        plot(x,inputm1)
        plot(x,inputm2)
      }  
      Spp.bps[[i]]<-out.f
      Low.breaks[[i]]<-lowbreaks[i]
      names(Spp.bps)[[i]]<-names(Low.breaks)[[i]]<-dat.names[i]
      likeProf <- PRegLikeProf.fn(out.f$par,x.in=Spp.dat.AFM.mod[,1],y.in=Spp.dat.AFM.mod[,2],rngSplit=rngSplit,steppin=steppin)
      likeProf$fDiff <- likeProf$f-out.f$value
      likeProf$CI <- likeProf$fDiff-2.5
      Prof.out[[i]]<-tmp <- likeProf
      tmp[tmp$f==min(tmp$f),]
      #plot graphs using the modified data with piecewise curves.
      plot(tmp$split,tmp$fDiff,type="l",xlab="Breakpoint",ylab="-LogLike Difference from minimum",main=dat.names[i],lwd=2)
      abline(h=qchisq(0.95,1)/2,col="red",lwd=1.5,lty=2)
    }
    #Create list of modified data, breakpoints, and minimization results. Standardize names.
    names(Prof.out)<-dat.names
    Dat.Bps.out<-list(Spp.dat.AFM,Spp.bps,Low.breaks,Prof.out)
    names(Dat.Bps.out)<-c("Data","Bps","AgeBreaks","Profile")
    return(Dat.Bps.out)
   }

  if(Bp.find=="F")
  {
  out.f.jitter<-NA
  Lm.out<-list()
  Lm.names<-c("Bp1","Bp2","NoBp")
  for(i in 1:length(Spp.dat.AFM))
  {
    Spp.dat.AFM.mod<-Spp.dat.AFM[[i]][Spp.dat.AFM[[i]][,1]>lowbreaks[i],]
    #Piecewise models
    dat.bp1<-Spp.dat.AFM.mod[Spp.dat.AFM.mod[,1]<=breakpts[i],]
    dat.bp2<-Spp.dat.AFM.mod[Spp.dat.AFM.mod[,1]>breakpts[i],]
    if(intcpt=="Y")
    {
      if(length(dat.bp1)>2) {pw1.lm<-lm(dat.bp1[,2]~dat.bp1[,1])}
      if(length(dat.bp1)<=2) {pw1.lm<-NA}
      if(length(dat.bp2)>2) {pw2.lm<-lm(dat.bp2[,2]~dat.bp2[,1])}
      if(length(dat.bp2)<=2) {pw2.lm<-NA}
      #Linear models
      lm.all<-lm(Spp.dat.AFM.mod[,2]~Spp.dat.AFM.mod[,1])
      Lm.out[[i]]<-list(pw1.lm,pw2.lm,lm.all)
      names(Lm.out[[i]])<-Lm.names
    }
    if(intcpt=="N")
    {
      if(length(dat.bp1)>2) {pw1.lm<-lm(dat.bp1[,2]~dat.bp1[,1]-1)}
      if(length(dat.bp1)<=2) {pw1.lm<-NA}
      if(length(dat.bp2)>2) {pw2.lm<-lm(dat.bp2[,2]~dat.bp2[,1])}
      if(length(dat.bp2)<=2) {pw2.lm<-NA}
      #Linear models
      lm.all<-lm(Spp.dat.AFM.mod[,2]~Spp.dat.AFM.mod[,1]-1)
      Lm.out[[i]]<-list(pw1.lm,pw2.lm,lm.all)
      names(Lm.out[[i]])<-Lm.names
    }
  }
  names(Lm.out)<-dat.names
  colnames(Spp.dat.AFM[[1]])<-colnames(Spp.dat.AFM[[2]])<-colnames(Spp.dat.AFM[[3]])<-c("OtoWt","Age")
  Dat.Bp.Lm.out<-list(Spp.dat.AFM,Lm.out,lowbreaks,breakpts,Lts,out.f.jitter)
  names(Dat.Bp.Lm.out)<-c("Data","LMs","Low_wt_breaks_used","Bps_used","Lt","Jitter")
  return(Dat.Bp.Lm.out)
  }
}

#Fit PW and LM models on same plot
OtoAge.plot<-function(model.in,data.choice,col.bpts=c("orange","blue"))
{
  age.break.low.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]<=model.in$Low_wt_breaks_used[data.choice],]
  age.break.hi.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]>model.in$Low_wt_breaks_used[data.choice],]
  dat.bp1<-age.break.hi.dat[age.break.hi.dat[,1]<=model.in$Bps[data.choice],]
  dat.bp2<-age.break.hi.dat[age.break.hi.dat[,1]>model.in$Bps[data.choice],]
  plot(age.break.hi.dat,pch=21,col="black",bg="gray",xlim=c(0,ceiling(max(model.in$Data[[data.choice]][,1]))),ylim=c(0,ceiling(max(model.in$Data[[data.choice]][,2]))),xlab="Otolith weight",ylab="Age")
  points(age.break.low.dat,pch=21,col="black",bg="lightblue")
  abline(v=model.in$Bps[data.choice],col="red",lty=2)
  abline(v=model.in$Low_wt_breaks_used[data.choice],col="red",lty=3)
  if(dim(dat.bp1)[1]!=0){lines(dat.bp1[,1],model.in$LMs[[data.choice]]$Bp1$fitted,col=col.bpts[1], lwd=2)}
  if(dim(dat.bp2)[1]!=0){lines(dat.bp2[,1],model.in$LMs[[data.choice]]$Bp2$fitted,col=col.bpts[2],lwd=2)}
  if(dim(age.break.hi.dat)[1]!=0){lines(age.break.hi.dat[,1],model.in$LMs[[data.choice]]$NoBp$fitted,col="black",lwd=2)}
}

LM.comp.plots<-function(model.in,data.choice,add.plot=F,col.bpts=c("orange","blue"))
{
  age.break.low.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]<=model.in$Low_wt_breaks_used[data.choice],]
  age.break.hi.dat<-model.in$Data[[data.choice]][model.in$Data[[data.choice]][,1]>model.in$Low_wt_breaks_used[data.choice],]
  dat.bp1<-age.break.hi.dat[age.break.hi.dat[,1]<=model.in$Bps[data.choice],]
  dat.bp2<-age.break.hi.dat[age.break.hi.dat[,1]>model.in$Bps[data.choice],]
  if(add.plot==F)
  {
    if(dim(dat.bp1)[1]!=0){plot(dat.bp1[,1],model.in$LMs[[data.choice]]$Bp1$fitted,col=col.bpts[1], lwd=2,type="l",xlim=c(0,ceiling(max(model.in$Data[[data.choice]][,1]))),ylim=c(0,ceiling(max(model.in$Data[[data.choice]][,2]))),xlab="Otolith weight",ylab="Age")}
    if(dim(dat.bp2)[1]!=0){lines(dat.bp2[,1],model.in$LMs[[data.choice]]$Bp2$fitted,col=col.bpts[2],lwd=2)}
    #if(dim(age.break.hi.dat)[1]!=0){plot(age.break.hi.dat[,1],model.in$LMs[[data.choice]]$NoBp$fitted,col="black",lwd=2,type="l")}
  }
  if(add.plot==T)
  {
    if(dim(dat.bp1)[1]!=0){lines(dat.bp1[,1],model.in$LMs[[data.choice]]$Bp1$fitted,col=col.bpts[1], lwd=2,lty=2)}
    if(dim(dat.bp2)[1]!=0){lines(dat.bp2[,1],model.in$LMs[[data.choice]]$Bp2$fitted,col=col.bpts[2],lwd=2,lty=2)}
    #if(dim(age.break.hi.dat)[1]!=0){lines(age.break.hi.dat[,1],model.in$LMs[[data.choice]]$NoBp$fitted,col="black",lwd=2)}
  }  
}

#Get numbers
Spp.samp.size<-function(dat.list.in)
{
  Spp.N<-matrix(NA,nrow=length(dat.list.in),ncol=3)
  colnames(Spp.N)<-c("A","F","M")
  rownames(Spp.N)<-names(dat.list.in)
  for(i in 1:length(dat.list.in))
    {Spp.N[i,]<-c(dim(dat.list.in[[i]]$Data$All)[1],dim(dat.list.in[[i]]$Data$Females)[1],dim(dat.list.in[[i]]$Data$Males)[1])}
  return(Spp.N)
}

#Get R^2 values
Spp.r2<-function(dat.list.in)
{
  rtwos<-as.data.frame(matrix(NA,nrow=2*length(dat.list.in),ncol=5))
  colnames(rtwos)<-c("Species","Sex","PW1","PW2","LM")
  genders<-c("Female","Male")
  spp.names<-names(dat.list.in)
  x<-1
  for(i in 1:length(dat.list.in))
  {
    for(g in 1:2)
    {
      for(ii in 1:3)
      {
        rtwos[x,1]<-spp.names[i]
        rtwos[x,2]<-genders[g]
        if(is.na(dat.list.in[[i]]$LMs[[g+1]][[ii]])==FALSE){rtwos[x,2+ii]<-summary(dat.list.in[[i]]$LMs[[g+1]][[ii]])$r.squared}
      }
      x<-1+x
    }
  }
  return(rtwos)
}


#Calcualte predicted ages
#low.breaks: matrix of values with breaks for the low ages (col 1= wt; col 2= age)
#intercept: If > 0, use non-zero intercept
Pred.ages.oto<-function(dat.list.in,low.breaks,intercept=999)
{
  Pred.ages<-list()
  for(i in 1:length(dat.list.in))
  {
    AFM.oto.age.temp<-list()
    for(ii in 1:3)
    {
      #Define low breaks
      low.ages.ind<-dat.list.in[[i]]$Data[[ii]][,1]<=low.breaks[[i]][nrow(low.breaks[[i]]),1]
      if(length(dat.list.in[[i]]$Data[[ii]][low.ages.ind,])>2){low.ages.dat<-dat.list.in[[i]]$Data[[ii]][low.ages.ind,]}
      if(length(dat.list.in[[i]]$Data[[ii]][low.ages.ind,])<=2){low.ages.dat<-matrix(dat.list.in[[i]]$Data[[ii]][low.ages.ind,],nrow=1,ncol=2)}
      low.ages<-matrix(NA,nrow=nrow(low.ages.dat),ncol=2)
      bp.ages.ind<-dat.list.in[[i]]$Data[[ii]][,1]>low.breaks[[i]][nrow(low.breaks[[i]]),1]
      if(length(dat.list.in[[i]]$Data[[ii]][bp.ages.ind,])>2){bp.ages.dat<-dat.list.in[[i]]$Data[[ii]][bp.ages.ind,]}
      if(length(dat.list.in[[i]]$Data[[ii]][bp.ages.ind,])<=2){bp.ages.dat<-matrix(dat.list.in[[i]]$Data[[ii]][bp.ages.ind,],nrow=1,ncol=2)}
      bp.ages<-matrix(NA,nrow=nrow(bp.ages.dat),ncol=2)
      #Predict the low end breaks
      for(iii in 1:nrow(low.breaks[[i]]))
      {
       if(iii==1)
       {
        low.ages[low.ages.dat[,1]>0&low.ages.dat[,1]<=low.breaks[[i]][iii,1],]<-low.breaks[[i]][iii,2]

       }
        
       else{low.ages[low.ages.dat[,1]>low.breaks[[i]][iii-1,1]&low.ages.dat[,1]<=low.breaks[[i]][iii,1],]<-low.breaks[[i]][iii,2]}
      }
      #Predict the piecewise breaks
      bp1.temp<-bp.ages.dat[,1]<=dat.list.in[[i]]$Bps[1]
      bp2.temp<-bp.ages.dat[,1]>dat.list.in[[i]]$Bps[1]
      if(intercept==0)
      {
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp1)>1){bp.ages[bp1.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp1$coef*bp.ages.dat[bp1.temp,1]}
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp2)>1){bp.ages[bp2.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp2$coef*bp.ages.dat[bp2.temp,1]}
      }
      if(intercept!=0)
      {
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp1)>1){bp.ages[bp1.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp1$coef[2]*bp.ages.dat[bp1.temp,1]+dat.list.in[[i]]$LMs[[ii]]$Bp1$coef[1]}
        if(length(dat.list.in[[i]]$LMs[[ii]]$Bp2)>1){bp.ages[bp2.temp,1]<-dat.list.in[[i]]$LMs[[ii]]$Bp2$coef[2]*bp.ages.dat[bp2.temp,1]+dat.list.in[[i]]$LMs[[ii]]$Bp2$coef[1] }
      }
      #Predict from linear model
      if(intercept==0){bp.ages[,2]<-dat.list.in[[i]]$LMs[[ii]]$NoBp$coef*bp.ages.dat[,1]}
      if(intercept!=0){bp.ages[,2]<-dat.list.in[[i]]$LMs[[ii]]$NoBp$coef[2]*bp.ages.dat[,1]+dat.list.in[[i]]$LMs[[ii]]$NoBp$coef[1]}
      low.bp.ages<-matrix(NA,nrow=nrow(dat.list.in[[i]]$Data[[ii]]),ncol=2)
      low.bp.ages[low.ages.ind]<-low.ages
      low.bp.ages[bp.ages.ind]<-bp.ages
      names(low.bp.ages)<-c("Bp","noBp")
      AFM.oto.age.temp[[ii]]<-as.data.frame(cbind(dat.list.in[[i]]$Data[[ii]],low.bp.ages,dat.list.in[[i]][[5]][[ii]]))
      colnames(AFM.oto.age.temp[[ii]])<-c("OtoWt","Age","PW_age","Lm_age","Length")
    }
    names(AFM.oto.age.temp)<-c("A","F","M")
    Pred.ages[[i]]<-AFM.oto.age.temp
  }
  names(Pred.ages)<-names(dat.list.in)
  return(Pred.ages)
}

#Piecewise models
PReg.obj <- function(x,x.in,y.in,rngSplit=c(75,180),fixSplt=NULL,verbose=F) {
    #piecewise regression, assuming continuous function at breakpoint.
    #Therefore, same intercept for both equations
    #note: b1 and b2 are slopes and alpha is the y-intercept. splt is the breakpoint.
    splt <- x[1]
    alpha <- x[2]
    b1 <- x[3]
    b2 <- x[4]
    #nobs = number observed
    nobs <- length(x.in)
    #creates an empty vecter with 'nobs' number of elements.
    predY <- rep(NA,nobs)
    #take lower half of data below breakpoint guess
    ind <- x.in <= splt
    
    #plug data into line equation
    predY[ind] <- alpha+b1*(x.in[ind]-splt)
    #difference of squares
    f.left <- sum((predY[ind]-y.in[ind])^2)
    #same as above, but with upper half of data
    ind <- x.in > splt
    predY[ind] <- alpha+b2*(x.in[ind]-splt)
    f.right <- sum((predY[ind]-y.in[ind])^2)
   
    #Sum all the differences of squares (this is to be minimized)
    f <- (f.left+f.right)
    if(verbose){cat(f.left,f.right,f,"\n")}

    #if breakpoint is out of domain, return arbitrarily large number as a flag, and so optim will avoid
    if(splt<rngSplit[1] | splt>rngSplit[2]) {
        return(1e10)
    }
    return((nobs/2)*log(f/nobs))    #the likelihood function (assuming normal errors and constant variance)
}



PRegLikeProf.fn <- function(x,x.in,y.in,rngSplit,steppin) {
  PRegLP.obj <- function(x,breakpt,x.in,y.in)
  {
    #piecewise regression, assuming continuous function at breakpoint.
    #Therefore, same intercept for both equations
    #Likelihood profile over splt
    splt <- breakpt
    alpha <- x[1]
    b1 <- x[2]
    b2 <- x[3]
    nobs <- length(x.in)
    predY <- rep(NA,nobs)
    ind <- x.in <= splt
    predY[ind] <- alpha+b1*(x.in[ind]-splt)
    f.left <- sum((predY[ind]-y.in[ind])^2)
    ind <- x.in > splt
    predY[ind] <- alpha+b2*(x.in[ind]-splt)
    f.right <- sum((predY[ind]-y.in[ind])^2)
    f <- (f.left+f.right)
    return((nobs/2)*log(f/nobs))    #the likelihood function (assuming normal errors and constant variance)
  }
    nobs <- length(x.in)
    splits <- seq(rngSplit[1],rngSplit[2],by=steppin)
    out <- matrix(NA,nrow=length(splits),ncol=5,dimnames=list(NULL,c("split","alpha","b1","b2","f")))
    for(i in 1:length(splits))
    {
        opt <- optim(x[-1],PRegLP.obj,breakpt=splits[i],x.in=x.in,y.in=y.in)
        out[i,] <- c(splits[i],opt$par,opt$value)
    }
    as.data.frame(out)
}

#Plot ages comparisons (1:1)
#gen.type:1=All; 2=Females; 3=Males
#age.comp.type:3=PW; 4=LM
#rd.type:1= none; 2=round; 3=ceiling
plot.age2age<-function(spp.dat,spp.dat.nam,gen.type=1,age.comp.type=3,rd.type=1,col.dots="black" ,lab.nam=c("T","T"),add.plot="F")
{
  name.nums<-c(1:length(names(spp.dat)))
  name.num<-name.nums[names(spp.dat)==spp.dat.nam]
  xlab.nam<-ylab.nam<-""
  if(lab.nam[1]=="T"){xlab.nam<-"Age from otolith annuli"}
  if(lab.nam[2]=="T"){ylab.nam<-"Age from otolith weight"}
  if(add.plot=="F")
  {
    if(rd.type==1)
    {
      plot(Spp.age.pred[[name.num]][[gen.type]][,2],Spp.age.pred[[name.num]][[gen.type]][,age.comp.type],xlim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),ylim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),xlab=xlab.nam, ylab=ylab.nam,pch=21,bg=col.dots)
      abline(lm(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==2)
    {
      plot(Spp.age.pred[[name.num]][[gen.type]][,2],round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),xlim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),ylim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),xlab=xlab.nam, ylab=ylab.nam,"Age from otolith weight",pch=21,bg=col.dots)
      abline(lm(round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==3)
    {
      plot(Spp.age.pred[[name.num]][[gen.type]][,2],ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),xlim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),ylim=c(0,max(Spp.age.pred[[name.num]][[gen.type]][,c(-1,-5)])),xlab=xlab.nam, ylab=ylab.nam,"Age from otolith weight",pch=21,bg=col.dots)
      abline(lm(ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    abline(a=0,b=1,col="black",lwd=2)
  }
  if(add.plot=="T")
  {
    if(rd.type==1)
    {
      points(Spp.age.pred[[name.num]][[gen.type]][,2],Spp.age.pred[[name.num]][[gen.type]][,age.comp.type],pch=21,bg=col.dots)
      abline(lm(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==2)
    {
      points(Spp.age.pred[[name.num]][[gen.type]][,2],round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),pch=21,bg=col.dots)
      abline(lm(round(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    if(rd.type==3)
    {
      points(Spp.age.pred[[name.num]][[gen.type]][,2],ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type]),pch=21,bg=col.dots)
      abline(lm(ceiling(Spp.age.pred[[name.num]][[gen.type]][,age.comp.type])~Spp.age.pred[[name.num]][[gen.type]][,2]),lty=2,lwd=1,col=col.dots)
    }
    abline(a=0,b=1,col="black",lwd=2)
  }
}

### VBGF functions
VBGF<-function(Linf,k,t0,ages)
{
Lengths_exp<-Linf*(1-exp(-k*(ages-t0)))
return(Lengths_exp)
}

VBGF.fit<-function(p,obs,return.type=2)
{
    if(return.type==1)
    {
        exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
        return(exp.lts)
    }
    if(return.type==2)
    {
        exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
        sigma<-sqrt(sum((obs[,2]-exp.lts)^2)/length(obs[,2]))
        neglogsum<-sum(-log(dnorm(exp.lts,obs[,2],sigma)))
        return(neglogsum)
    }
}

VBGF.fit.plot<-function(AG.age.in,AG.lt.in,max.x,max.y,pch.col="red",title.in="",add.plot="F",print.vbgf="T")
{
  fitvbgf<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=AG.age.in,lts=AG.lt.in),start=list(Linf=max(AG.lt.in,na.rm=TRUE),k=0.1,t0=0),control = list(reltol=0.00000000001))
  AG.fit.par<-summary(fitvbgf)$par
  exp.AG<-VBGF(AG.fit.par[1],AG.fit.par[2],AG.fit.par[3],c(0:max(AG.age.in)))
  if(add.plot=="F")
  {
    plot(AG.age.in,AG.lt.in,xlab="Age (years)", ylab="Length (cm)",xlim=c(0,max.x),ylim=c(0,max.y),pch=21,bg=pch.col,col="black",main=title.in)
    lines(c(0:max(AG.age.in)),exp.AG,lwd=2)
  }
  if(add.plot=="T")
  {
    points(AG.age.in,AG.lt.in,pch=21,bg=pch.col,col="black",main=title.in)
    lines(c(0:max(AG.age.in)),exp.AG,lwd=2,lty=2)
  }
  if(print.vbgf=="T")
  {
  text(0.75*max.x,0.25*max.y,labels=bquote(paste(L[infinity],"= ",.(round(AG.fit.par[1],2)))))
  text(0.75*max.x,0.2*max.y,labels=paste("k= ",round(AG.fit.par[2],2)))
  text(0.75*max.x,0.15*max.y,labels=bquote(paste(t[0],"= ",.(round(AG.fit.par[3],2)))))
  }
  return(AG.fit.par)
}
###

vbgf.pars<-function(spp.dat.in, max.x,max.y)
{
  gender<-c("A", "F", "M")
  age.type<-c("OtA","PwA","LmA")
  num.type<-c("","Rd","Cl")
  temp.out<-list(matrix(NA,nrow=21,ncol=3),matrix(NA,nrow=21,ncol=3))
  r.names<-rep(NA,21)
  colnames(temp.out[[1]])<-colnames(temp.out[[2]])<-c("Linf","k","t0")
  x<-1
  for(i in 1:length(gender))
  {
    for(ii in 1:length(age.type))
    {
      if(ii==1)
      {
        temp.out.parms<-VBGF.fit.plot(spp.dat.in[[i]][[ii+1]],spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],sep="")
        x<-x+1
      }
      if(ii>1)
      {
        temp.out.parms<-VBGF.fit.plot(spp.dat.in[[i]][[ii+1]],spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],"_",num.type[1],sep="")
        x<-x+1
        temp.out.parms<-VBGF.fit.plot(round(spp.dat.in[[i]][[ii+1]]),spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],"_",num.type[2],sep="")
        x<-x+1
        temp.out.parms<-VBGF.fit.plot(ceiling(spp.dat.in[[i]][[ii+1]]),spp.dat.in[[i]][[5]],round(max(spp.dat.in[[i]][[ii+1]])/10)*10,round(max(spp.dat.in[[i]][[5]],na.rm=TRUE)/10)*10,"gray","")
        temp.out[[1]][x,]<-temp.out.parms[,1]
        temp.out[[2]][x,]<-temp.out.parms[,2]
        r.names[x]<-paste(gender[i],"_",age.type[ii],"_",num.type[3],sep="")
        x<-x+1
      }
    }
  }
  rownames(temp.out[[1]])<-rownames(temp.out[[2]])<-r.names
  names(temp.out)<-c("Est","Std")
  return(temp.out)
}

vbgf.par.comp.plot<-function(spp.par.in,rows2use=c(1:21),vbgf.num,labelx="Y",dot.col=c("black","orange","orange","orange","blue","blue","blue"))
{
  par.names<-expression(L[inf],"k",t[0])
  plotCI(spp.par.in$Est[rows2use,vbgf.num],uiw=spp.par.in$Std[rows2use,vbgf.num]*1.96,ylim=c(0, max(spp.par.in$Est[rows2use,vbgf.num]+spp.par.in$Est[rows2use,vbgf.num]*0.25)),pch=21,pt.bg=dot.col,xlab="",ylab=par.names[vbgf.num],axes=F)
  box()
  if(labelx=="Y"){axis(1,at=c(1:dim(spp.par.in$Est)[1]),labels=rownames(spp.par.in$Est)[rows2use],las=2,cex=0.75)}
  if(labelx=="N"){axis(1,at=c(1:dim(spp.par.in$Est)[1]),labels=rep("",dim(spp.par.in$Est)[1]))}
  axis(2)
  if(length(rows2use)/7==3){abline(v=c(7.5,14.5),lty=2)}
  if(length(rows2use)/7==2){abline(v=c(7.5),lty=2)}
}

#Relative errors
Spp.RE<-function(spp.par.in,vbgf.num=1)
{
  F_u<-(spp.par.in$Est[c(9,12),vbgf.num]-spp.par.in$Est[8,vbgf.num])/spp.par.in$Est[8,vbgf.num]
  F_uw<-((spp.par.in$Est[c(9,12),vbgf.num]+spp.par.in$Std[c(9,12),vbgf.num]*1.96)-spp.par.in$Est[8,vbgf.num])/spp.par.in$Est[8,vbgf.num]
  M_u<-(spp.par.in$Est[c(16,19),vbgf.num]-spp.par.in$Est[15,vbgf.num])/spp.par.in$Est[15,vbgf.num]
  M_uw<-((spp.par.in$Est[c(16,19),vbgf.num]+spp.par.in$Std[c(16,19),vbgf.num]*1.96)-spp.par.in$Est[15,vbgf.num])/spp.par.in$Est[15,vbgf.num]
  FM.out<-rbind(F_u,F_uw,M_u,M_uw)
  colnames(FM.out)<-c("PW","LM")
  return(FM.out)
}
#############################################
#############################################
#############################################
Punt_age_prep<-function(Spp.Ages.in)
{
  age.in<-table(round(Spp.Ages.in))
  age.out<-matrix(NA,nrow=length(age.in),ncol=5)
  age.out[,c(1,4)]<-0
  age1<-c(1:dim(age.in)[1])
  age2<-c(1:dim(age.in)[2])
  x<-1
  for(i in 1:length(age1))
  {
    for(ii in 1:length(age2))
    {
      age.out[x,c(2,3,5)]<-c(age1[i],age2[ii],age.in[i,ii])
      x<-x+1
    }
  }
  age.return<-subset(age.out,age.out[,5]>0)
  colnames(age.return)<-c("0","Age1","Age2","0","N")
  return(age.return)
}

#plot SD from Punt program
plot.AErr.comp<-function(dat.in,species.in,col.in="orange")
{
dat.plot<-subset(dat.in,Species==species.in)
plot(subset(dat.plot,Reader==1)$SD,subset(dat.plot,Reader==2)$SD,xlim=c(0,round(max(c(subset(dat.plot,Reader==1)$SD,subset(dat.plot,Reader==2)$SD)))),ylim=c(0,round(max(c(subset(dat.plot,Reader==1)$SD,subset(dat.plot,Reader==2)$SD)))),xlab="",ylab="",pch=21,bg=col.in,cex=1.25)
abline(a=0,b=1,lwd=2,col="red")
}

