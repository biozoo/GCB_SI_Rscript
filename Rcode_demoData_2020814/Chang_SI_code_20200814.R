### CHANG_ETAL_SUPPLEMENTARY_INFORMATION_RCODE - Full Code
### Last Updated Aug 14, 2020
rm(list = ls())
library(dplyr) 
library(rEDM) # Empirical dynamical modeling for CCM analysis
library(Kendall) # Kendall's tau test for the convergence of CCM
library(igraph) # Plot network objects

# Set working directory
setwd("D:/data/Meta_analysis/Causality/Manuscript/GCB/Rcode_demoData_2020814")
######################################################################################
###### The reconstruction of causality networks ######################################
######################################################################################
# Loading time series data for Lake Kasumigaura station 9 (Ks9) ####
# Note that the use of Lake Kasumigaura dataset should strictly follow the policy 
# stated in the following website http://db.cger.nies.go.jp/gem/moni-e/inter/GEMS/database/kasumi/
month.dat <- read.csv("Demo_Ks9_tsdata.csv") # Loading dataset with raw variables
colnames(month.dat)=c('Year','Month','Richness','Biomass','Temperature','NO3','PO4')# rename the variables
n=nrow(month.dat) # time series length
seed=25647

####Function for time series standardization (normalization, detrend and deseason)
# x: time series data (vector form)
# normalization: normalization to zero mean & unit variance 
# dseason: logic argument for Deseasonalization by monthly mean
# season_sd: Deseasonalization by both monthly mean and monthly standard deviation
# sea: The period of seasonality (e.g. 12 months for monthly data)
# dtrend: logic argument for detrend
# dTtype: Type of detrend by first difference (first) or linear regression (linear)
nomz=function(x, normalization=T, dseason=T, season_sd=T, sea=12, dtrend=T, dTtype="linear"){
  x=as.numeric(x)
  xt=x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt=diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t=lm(xt~c(1:length(xt)))
    xt=xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt=xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt=xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}


####### Function for generating lag time series 
laf=function(x,y,lagf){
  n <- NROW(x)
  x.t=x;y.t=y
  if(lagf<=0){x.t=x.t[(1-lagf):n];y.t=y.t[1:(n+lagf)]} # if lagf<0, y is leading
  if(lagf>0){x.t=x.t[1:(n-lagf)];y.t=y.t[(1+lagf):n]}  # if lagf>0, x is leading           
  return(cbind(x.t,y.t))
}

# Detrend + deseason time series
sdat=data.frame(month.dat[,c(1:2)],apply(month.dat[,-c(1:2)],2,nomz))
##The index for testing causal links (a total of 12 links)
indmat=matrix(0,12,2);colnames(indmat)=c('Effect','Cause')
indmat[,1]=c('Biomass','NO3','PO4','Richness','NO3','PO4','Richness','Biomass','Richness','Biomass','Richness','Biomass')
indmat[,2]=c('Richness','Richness','Richness','Biomass','Biomass','Biomass','Temperature','Temperature','NO3','NO3','PO4','PO4')
indmat

# Determine the embedding dimensions (En) in CCM by try-and-error with best hindcast (tp=-1) skill
Emax=20
En=NULL
for(i in 1:nrow(indmat)){
  E.test=NULL
  for(E.t in 2:Emax){
    cmxy.t <- ccm(sdat, E = E.t,
                  lib_column = indmat[i,1], target_column = indmat[i,2],
                  lib_sizes = n, tp=-1,random_libs = F) 
    E.test=c(E.test,mean(cmxy.t$rho))
  }
  # Select the embedding dimension that makes the model with the highest hindcast (tp=-1) predictive skill
  En=c(En,which.max(E.test)+1) 
}


################################################################################
### CCM analysis for all causality testlinks
lib_siz=sort(c(5,10,20,30,40,seq(50,n,50),n)) # a sequence of library size
ccmda=NULL
for(i in 1:nrow(indmat)){
  ccmda.t=NULL
  for(j in 0:-3){
    da.t=laf(sdat[,indmat[i,1]],sdat[,indmat[i,2]],lagf=j) # Varying time lags
    colnames(da.t)=indmat[i,]
    # CCM analysis cross-mapping from an effect variable to its cause
    x_xmap_y <- ccm(da.t, E = En[i], # The embedding dimension E for each link were determined in previous step
                    lib_column = indmat[i,'Effect'], target_column = indmat[i,'Cause'],
                    lib_sizes = lib_siz, tp=0,RNGseed = seed,
                    num_samples = 100,replace=F)
    
    # Take average for the predictive skill under each library size
    aveg=cbind(unique(x_xmap_y$lib_size),aggregate(x_xmap_y[,c('rho')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('mae')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('rmse')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'])
    ccm_mean=data.frame(lag=rep(j,nrow(aveg)),x_xmap_y[1:nrow(aveg),]);
    ccm_mean[,c('lib_size','rho','mae','rmse')]=aveg
    ccm_mean[ccm_mean[,'rho']<0,'rho']=0
    
    ###########################
    # Convergence test in CCM
    ###########################
    # Fisher's delta rho Z test
    rho.Lmax=ccm_mean$rho[which.max(ccm_mean$lib_size)]
    rho.Lmin=ccm_mean$rho[1]
    ns=min(sum(!is.na(sdat[,indmat[i,1]])),sum(!is.na(sdat[,indmat[i,2]])))
    delta_rho=rho.Lmax-rho.Lmin
    z=abs(0.5*(log((1+rho.Lmax)/(1-rho.Lmax))-log((1+rho.Lmin)/(1-rho.Lmin)))*(2/(ns-3))^-0.5)
    z.p=(1-pnorm(z))
    # Kendall's tau test
    if(length(ccm_mean$rho)>3){
      kend=MannKendall(ccm_mean$rho)
      kend.tau=kend$tau[1]
      kend.p=kend$sl[[1]]
    }else{
      kend.tau=NA
      kend.p=NA
    }
    # Compile all the testing results
    ccmda.t=rbind(ccmda.t,
                  unlist(c(ccm_mean[1,1:5],rho_Lmax=rho.Lmax,rho_Lmin=rho.Lmin,
                           Z=z,p_Z=z.p,Kendall_tau=kend.tau,Kendall_p=kend.p)))
  }
  # Select the CCM results based on predictive skills 
  ccmda=rbind(ccmda,ccmda.t[which.max(ccmda.t[,'rho_Lmax']),])
}

ccmda=data.frame(indmat,ccmda)
Convergence=ccmda$p_Z<0.05 & ccmda$Kendall_p<=0.05 & ccmda$Kendall_tau>0
(ccmda=data.frame(ccmda,Convergence))

# Standardized linkage strength by dividing the maximal linkage strength within the system
istd=data.frame(system=rep('Ks9'),ccmda[,c('Cause','Effect','rho_Lmax','p_Z','Kendall_tau','Kendall_p','Convergence')],
                Standardized_linkage_strength=ccmda$rho_Lmax/max(ccmda$rho_Lmax[1:12]))

# The matrix of linkage strength (i,j) (i=effect & j=cause) 
linkM=matrix(0,ncol(sdat)-2,ncol(sdat)-2)
colnames(linkM)=rownames(linkM)=colnames(sdat)[-c(1:2)]
for(i in 1:nrow(istd)){
  if(istd$Convergence[i]==TRUE){
    linkM[rownames(linkM)==ccmda[i,1],colnames(linkM)==ccmda[i,2]]=istd$Standardized_linkage_strength[i]
  }
}
linkM

# Loading the library for network plotting 
iweb <- graph.adjacency(t(linkM), mode="directed",weighted=T)
wd=linkM*10
wd=wd[wd>0]
E(iweb)$width=wd
tkplot(iweb,edge.curved=0.2,vertex.label.color="black",vertex.color="white",
       edge.color="black",vertex.size=30,layout=layout_in_circle(iweb))


#################################################################################
####The computation of ecosystem stability (1/CV of the phytoplankton biomass; Chla)
ym=month.dat[,"Biomass"]
n=length(ym)
xm=c(1:n)
xs=matrix(c(ym,rep(NA,12-n%%12)),12,n%/%12+1)
xs=apply(xs,1,mean,na.rm=T)
yrs=ym-rep(xs,n%/%12+1)[1:n]
yrs=yrs-predict(lm(yrs~xm),data.frame(xm=xm))  
mym=mean(ym,na.rm=T)
(Ecosystem_stability.t=1/c(sd(yrs,na.rm=T)/mym))

# The computation of mean species richness
(Mean_richness.t=mean(month.dat[,'Richness'],na.rm=T))

######################################################
##### The computation of warming rate  ###############
library(zyp)# Theil-Sen median based estimator
temp.mon=data.frame(temp=month.dat[,"Temperature"],
                    month=c(1:nrow(month.dat)))
fit<- zyp.sen(temp~month,temp.mon)
(Warming_rate.t=fit$coefficients[2]*12)



############################################
##### Mirrage correlation (Supplemental Materials; Fig. S14)
yr=10 
wdw=yr*12 # size of moving-window
n=nrow(month.dat) # time series length
# Detrend + deseason time series
daa=data.frame(month.dat[,c(1:2)],apply(month.dat[,-c(1:2)],2,nomz))
corr=NULL # correlation coefficients
qur=NULL # quadratic coefficient
slopp=NULL # linear slope
modd=NULL # Linear or quadratic model decided by AIC 
na.t=rep(NA,n)
cormv=c()
quc=lmc=NULL
aici=NULL
for(i in wdw:nrow(daa)){
  x1=daa[(i-wdw+1):i,"Biomass"]
  x2=daa[(i-wdw+1):i,"Richness"]
  sna=sum((!is.na(x1))&(!is.na(x2)))
  if(sna<wdw/3){
    cormv=c(cormv,NA)
    quc=c(quc,NA)
    lmc=c(lmc,NA)
  }else{
    cormv=c(cormv,cor(x1,x2,use="complete.obs"))
    lm1=lm(x2~x1)
    lm2=lm(x2~x1+I(x1^2))
    if(AIC(lm1)<AIC(lm2)){aic.t="L"}else{aic.t="Q"} # select the optimal model based on AIC
    aici=c(aici,aic.t)
    quc=c(quc,lm2$coefficients[3])
    lmc=c(lmc,lm1$coefficients[2])
  }
}

# The time series of regression coefficients
corr=qur=slopp=modd=na.t;
corr[c((wdw:nrow(daa))-wdw/2)]=cormv
qur[c((wdw:nrow(daa))-wdw/2)]=quc
slopp[c((wdw:nrow(daa))-wdw/2)]=lmc
modd[c((wdw:nrow(daa))-wdw/2)]=aici

# plotting moving-window regression coefficients
win.graph(40,25)
y=corr
y1=y;y1[modd=="Q"]=NA
y2=qur; y2[modd=="L"]=0
range.y=range(c(y,y1,y2),na.rm=T)
x=seq(month.dat[1,'Year']+month.dat[1,'Month']/12-1/24,by=1/12,length.out=length(corr)) # Timing (year)
xm=x[which(is.na(y))[which(diff(which(is.na(y)))>10)+1]];xm=xm[length(xm)]
xn=x[which(is.na(y))[which(diff(which(is.na(y)))>10)-1]][1]
plot(y~x,ylab="",main='Mirrage correlation between phytopankton species richness and biomass (Ks9)',
     type="l",ylim=range.y,xlab="",xlim=c(xn,xm),mar=c(1,1,1,1),lty=3)
lines(y1~x,ylim=range.y,xlim=c(xn,xm),lty=1,cex=2,lwd=2)
lines(y2~x,type="l",ylim=range.y,xlim=c(xn,xm),ylab="",xlab="",col=4,lwd=2) 
abline(h=0,lty=1,col="grey80",lwd=2)
legend('bottomleft',c("Linear","Quadratic"),lty=c(1,1),col=c(1,4),bty='n',lwd=2)


####################################################################################
###¡@By repeating all these analyses for all systems, including causality network reconstruction and 
###  the computation of ecosystem stability and warming rate, we can carry out the following cross-system comparison.
###  The dataset from the other systems are also available as explained in Supplementary Table 2 
####################################################################################


##################################################################################
######################### cross system comparison ################################
#    Loading the dataset 
l.strength <- read.csv("Demo_link_strength.csv") # read the CCM data
warm.sb <- read.csv("Demo_warming_stability.csv") # read the data of ecosystem stability and warming rate
attach(warm.sb)

# The abbreviation for variables BD: Richness; EF: Biomass; Temp: Temperature
# The abbreviation for causal links, e.g., "Richness->Biomass" was abbreviated as "BD_EF"
link=c("BD_EF","BD_NO3","BD_PO4","EF_BD","EF_NO3","EF_PO4","Temp_BD","Temp_EF","NO3_BD","NO3_EF","PO4_BD","PO4_EF")
lsda=matrix(0,nrow(warm.sb),12);colnames(lsda)=link; 
lsda[,'BD_EF']=l.strength[l.strength[,'Cause']=='Richness'&l.strength[,'Effect']=='Biomass','Standardized_linkage_strength']
lsda[,'BD_NO3']=l.strength[l.strength[,'Cause']=='Richness'&l.strength[,'Effect']=='NO3','Standardized_linkage_strength']
lsda[,'BD_PO4']=l.strength[l.strength[,'Cause']=='Richness'&l.strength[,'Effect']=='PO4','Standardized_linkage_strength']
lsda[,'EF_BD']=l.strength[l.strength[,'Cause']=='Biomass'&l.strength[,'Effect']=='Richness','Standardized_linkage_strength']
lsda[,'EF_NO3']=l.strength[l.strength[,'Cause']=='Biomass'&l.strength[,'Effect']=='NO3','Standardized_linkage_strength']
lsda[,'EF_PO4']=l.strength[l.strength[,'Cause']=='Biomass'&l.strength[,'Effect']=='PO4','Standardized_linkage_strength']
lsda[,'Temp_BD']=l.strength[l.strength[,'Cause']=='Temperature'&l.strength[,'Effect']=='Richness','Standardized_linkage_strength']
lsda[,'Temp_EF']=l.strength[l.strength[,'Cause']=='Temperature'&l.strength[,'Effect']=='Biomass','Standardized_linkage_strength']
lsda[,'NO3_BD']=l.strength[l.strength[,'Cause']=='NO3'&l.strength[,'Effect']=='Richness','Standardized_linkage_strength']
lsda[,'NO3_EF']=l.strength[l.strength[,'Cause']=='NO3'&l.strength[,'Effect']=='Biomass','Standardized_linkage_strength']
lsda[,'PO4_BD']=l.strength[l.strength[,'Cause']=='PO4'&l.strength[,'Effect']=='Richness','Standardized_linkage_strength']
lsda[,'PO4_EF']=l.strength[l.strength[,'Cause']=='PO4'&l.strength[,'Effect']=='Biomass','Standardized_linkage_strength']
lsda=data.frame(lsda)
rownames(lsda)=System
lsda # Standardized linkage strength for all systems
attach(lsda)

##################################################################################
# Plot Fig. 1 ecosystem stability vs warming rate
dev.new(width=5, height=5); par(mfrow=c(1,1))
y=Ecosystem_stability;x=Warming_rate
plot(y~x,ylab="Ecosystem stability",xlab="Warming rate (oC/yr)",
     pch=16,cex=2,ylim=c(1.2,2.5),xlim=c(-0.02,0.15))
lm2=lm(y~x);summary(lm2)
slv=summary(lm2)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=0}
abline(lm2$coefficients[1],lm2$coefficients[2],lty=lly)
text(y=y*.95,x=x*0.95,System)
text(0.1,2.4,paste("R2=",round(summary(lm2)$r.squared,3)))
text(0.1,2.3,paste("p=",round(summary(lm2)$coefficients[2,4],3)))

# Statistical power of causal links for explaining ecosystem stability 
SP_link=data.frame(Link=colnames(lsda),slope=rep(NA),R2=rep(NA),p_value=rep(NA),sig=rep(NA),rmse=rep(NA),aic=rep(NA))
for(i in 1:ncol(lsda)){
  lm.t=lm(Ecosystem_stability~lsda[,i])
  SP_link[i,'slope']=round((lm.t)$coefficients[2],3)
  SP_link[i,'R2']=round(summary(lm.t)$r.squared,3)
  SP_link[i,'p_value']=round(summary(lm.t)$coefficients[2,4],3)
  SP_link[i,'sig']=sign(lm.t$coefficients[2])
  SP_link[i,'rmse']=sqrt(mean((summary(lm.t)$residuals)^2))
  SP_link[i,'aic']=AIC(lm.t)
}
(SP_link=SP_link[order(SP_link$aic),])

# The link with the lowest AIC is 1st. Richness->Biomass (BD_EF;stabilization), 
#  and followed by 2nd. Richness->NO3 (BD_NO3;stabilization) and 3rd. Temperature->Biomass (Temp_EF; destabilization).


# The optimal pathway in causality network explaining ecosystem stability
# The function for generating all possible pathways
hcom=function(n){
  m.t=matrix(rep(1,n*(2^n-1)),ncol=n)
  cbt=1
  for(i in 1:(n-1)){
    ct=combn(n,i)
    for(j in 1:ncol(ct))m.t[cbt+j,ct[,j]]=0;
    cbt=cbt+ncol(ct)
  }
  return(m.t)
}

# The function for calculating geometric mean
geo.mean <- function(data,na.rm=T) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)],na.rm=na.rm))
  if(any(data==0)){gm=0}
  return(gm)
}

# Ecosystem stability regresses against the geometric mean of linkage strength for all pathways
# All the regression models are equivalent simple linear regression models regardless of the complexity of pathways
phtc=hcom(12); # In total, 12 links generating 2^12-1 = 4095 pathways (link combinations) are examined
r2=aic=rmse=NULL
for(i in 1:nrow(phtc)){ 
  c.t=which(phtc[i,]!=0)
  if(length(c.t)==1){gmt=lsda[,c.t]}else{gmt=apply(lsda[,c.t],1,geo.mean)}
  lm1=lm(Ecosystem_stability~gmt)
  r2=c(r2,summary(lm1)$r.squared)
  aic=c(aic,AIC(lm1))
  rmse=c(rmse,sqrt(mean((summary(lm1)$residuals)^2)))
}


nlk=apply(phtc[,1:12],1,sum)
htcmon=cbind(phtc,AIC=aic,R2=r2,rMSE=rmse)
# Select the optimal pathway with the lowest AIC
# The pathway with minimal AIC is "Rich->NO3->Chla + Rich->Chla"
(optimal_path=link[htcmon[which.min(htcmon[,'AIC']),]==1])
(optimal_path.aic=htcmon[which.min(htcmon[,'AIC']),'AIC'])



##############################################################################
### Plot Fig.3  Mean linkage strength vs ecosystem stability
windows(width=10, height=10); 
# Alternatively, quartz(width=10, height=10) works for Mac; x11(width=10, height=10) works for Linux; 
par(mfrow=c(2,2),mar=c(4,4,2,1))
y=Ecosystem_stability

# The computation of Fig. 3a (rarefacted richness) needs original community data for estimating sampling effort
# For simiplicty, we only demostrate a similar analysis presented in Fig. S5 based on average species richness  
# Fig. S5 y:Ecosystem stability; x:average species richness
x=Mean_richness
plot(y~x,xlab='Mean species Rich', main= '',ylab="Ecosystem stability",cex=2,pch=16)
lm2=lm(y~x);#summary(lm2)
coll=1
slv=summary(lm2)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=3;coll='grey'}
abline(lm2$coefficients[1],lm2$coefficients[2],lty=lly,col=coll)
text(y=y*0.98,x=x*0.98,System)
text(16,2.2,paste("R2=",round(summary(lm2)$r.squared,3)))
text(16,2.1,paste("p=",round(summary(lm2)$coefficients[2,4],3)))

# Fig. 3b y:Ecosystem stability; x:linkage strength (Richness->Biomass)
x=BD_EF
plot(y~x,xlab='Linkage strength', main= '(Rich->Chla)',ylab="Ecosystem stability",cex=2,pch=16)
lm2=lm(y~x);#summary(lm2)
coll=1
slv=summary(lm2)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=3;coll='grey'}
abline(lm2$coefficients[1],lm2$coefficients[2],lty=lly,col=coll)
text(y=y*0.98,x=x*0.98,System)
text(0.4,2.2,paste("R2=",round(summary(lm2)$r.squared,3)))
text(0.4,2.1,paste("p=",round(summary(lm2)$coefficients[2,4],3)))

# Fig. 3c y:Ecosystem stability; x:mean linkage strength (Rich->NO3->Chla+Rich->Chla) 
x=apply(cbind(BD_NO3,NO3_EF,BD_EF),1,geo.mean);
plot(y~x,xlab='Mean linkage strength', main='Rich->NO3->Chla+Rich->Chla',ylab="Ecosystem stability",cex=2,pch=16)
lm2=lm(y~x);
coll=1
slv=summary(lm2)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=3;coll='grey'}
abline(lm2$coefficients[1],lm2$coefficients[2],lty=lly,col=coll)
text(y=y*0.98,x=x*0.98,System)
text(0.4,2.2,paste("R2=",round(summary(lm2)$r.squared,3)))
text(0.4,2.1,paste("p=",round(summary(lm2)$coefficients[2,4],3)))

# Fig. 3d y:Ecosystem stability; x:Linkage strength (NO3->Chla) 
x=NO3_EF
plot(y~x,xlab='Linkage strength', main= '(NO3->Chla)',ylab="Ecosystem stability",cex=2,pch=16)
lm2=lm(y~x);#summary(lm2)
coll=1
slv=summary(lm2)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=3;coll='grey'}
abline(lm2$coefficients[1],lm2$coefficients[2],lty=lly,col=coll)
text(y=y*0.98,x=x*0.98,System)
text(0.4,2.2,paste("R2=",round(summary(lm2)$r.squared,3)))
text(0.4,2.1,paste("p=",round(summary(lm2)$coefficients[2,4],3)))




##############################################################################
### Fig. 4 Warming  effects on linkage strengths #############################
##############################################################################
windows(width=8, height=16); 
# Alternatively, quartz(width=8, height=16) works for Mac; x11(width=8, height=16) works for Linux; 
par(mfcol=c(2,1),mar=c(4,4,2,1))
x=Warming_rate
# Fig. 4a y:Linkage strength (Rich->Chla); x:Warming rate
y=BD_EF;
plot(y~x,ylab="Linkage strength", xlab='',main='Rich->Chla',cex=2,pch=16)
lm1=lm(y~x);summary(lm1)
coll=1;slv=summary(lm1)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=3;coll='grey'}
abline(lm1$coefficients[1],lm1$coefficients[2],lty=lly,col=coll)
text(y=y*1.08,x=x*1.08,System)
text(0.1,0.8,paste("R2=",round(summary(lm1)$r.squared,3)))
text(0.1,0.7,paste("p=",round(summary(lm1)$coefficients[2,4],3)))

# Fig. 4b y:Linkage strength (Rich->NO3->Chla+Rich->Chla); x:Warming rate
y=apply(cbind(BD_NO3,NO3_EF,BD_EF),1,geo.mean)
plot(y~x,ylab="Mean linkage strength", xlab='', 
     main='Rich->NO3->Chla+Rich->Chla',cex=2,pch=16)
lm1=lm(y~x);summary(lm1)
coll=1
slv=summary(lm1)$coefficients[2,4];if(slv<=0.05){lly=1}else if(slv<=0.1){lly=2}else{lly=3;coll='grey'}
abline(lm1$coefficients[1],lm1$coefficients[2],lty=lly,col=coll)

text(y=y*0.95,x=x*0.95,System)
text(0.1,0.8,paste("R2=",round(summary(lm1)$r.squared,3)))
text(0.1,0.7,paste("p=",round(summary(lm1)$coefficients[2,4],3)))
