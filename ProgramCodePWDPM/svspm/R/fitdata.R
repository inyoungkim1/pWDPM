fitdata=function(y,p=1,Q=3,model='pm',Nm=50000,Nb=10000,trace=F){
ptm=proc.time()
if(model!="pm" && model!="dpm" && model!="wdpmp" && model!="wdpmr" && model!="ewdpm" &&
model!="wdpmg" && model!="wdpme"&& model!="wdpmh"){
stop("invalid model name")
}
if(model=="dpm"){
Q=1
p=1
}
if(model=="pm"){
result=cpp3(y,Nm,Nb)
}
if(model=="dpm" || model=="wdpmp"){
result=cpp4(y,p,Q,Nm,Nb)
}
if(model=="wdpmr"){
result=cpp5(y,p,Q,Nm,Nb)
}
if(model=="ewdpm"){
result=cpp6(y,p,Q,3,Nm,Nb)
}

if(model=="wdpmg"){
result=cpp6(y,p,Q,4,Nm,Nb)
}
if(model=="wdpme"){
result=cpp6(y,p,Q,5,Nm,Nb)
}
if(model=="wdpmh"){
result=cpp6(y,p,Q,6,Nm,Nb)
}
len=length(result)
a1=floor(Nb/10)+1
a2=floor(Nm/10)
tsize=a2-a1+1
if(len>8.5){
Macf=as.vector(acf(result$MT[a1:a2],lag.max=1000,plot=F)$acf)
Szacf=as.vector(acf(result$SzT[a1:a2],lag.max=1000,plot=F)$acf)
Snacf=as.vector(acf(result$SnT[a1:a2],lag.max=1000,plot=F)$acf)
Wacf=as.vector(acf(result$WT[a1:a2],lag.max=1000,plot=F)$acf)
fiacf=as.vector(acf(result$fiT[a1:a2],lag.max=1000,plot=F)$acf)
Mu0acf=as.vector(acf(result$Mu0T[a1:a2],lag.max=1000,plot=F)$acf)
if(Q>1){
psiacf=as.vector(acf(result$psiT[a1:a2],lag.max=1000,plot=F)$acf)
}
if(Q==1){
psiacf=c(1)
}
}
if(len<8.5){
Macf=c(1)
Szacf=c(1)
Snacf=as.vector(acf(result$SnT[a1:a2],lag.max=1000,plot=F)$acf)
Wacf=c(1)
psiacf=c(1)
fiacf=as.vector(acf(result$fiT[a1:a2],lag.max=1000,plot=F)$acf)
Mu0acf=as.vector(acf(result$MuT[a1:a2],lag.max=1000,plot=F)$acf)
}
effi_list=cpp7(Macf, Szacf, Snacf, Wacf, psiacf, fiacf, Mu0acf, len,tsize)
fitresult=list(LPS=result$LPS, Sn=result$Sn,effi_Sn=effi_list$Sn, fi=result$fi, effi_fi=effi_list$fi)
if(model=="pm"){
fitresult=c(fitresult,Mu=result$Mu, effi_mu=effi_list$Mu)
}
if(model!="pm"){
fitresult=c(fitresult,mu0=result$mu0,effi_Mu0=effi_list$Mu0,M=result$M,effi_M=effi_list$M,
Sz=result$Sz,effi_Sz=effi_list$Sz,W=result$W,effi_W=effi_list$W)
if(model=="wdpmp" || model=="wdpmr" || model=="ewdpm"){
fitresult=c(fitresult,psi=result$psi,effi_psi=effi_list$psi)
}
}
ptm=proc.time()-ptm
fitresult=c(fitresult,cpu_time=as.vector(ptm)[3])
if(trace==T){
SnT=list(SnT=result$SnT)
fiT=list(fiT=result$fiT)
if(model=="pm"){
MuT=list(MuT=result$MuT)
fitresult=c(fitresult,SnT,fiT,MuT)
}
if(model!="pm"){
Mu0T=list(Mu0T=result$Mu0T)
SzT=list(SzT=result$SzT)
MT=list(MT=result$MT)
WT=list(WT=result$WT)
fitresult=c(fitresult,SnT,fiT,Mu0T,SzT,MT,WT)
}
if(model=="wdpmp" || model=="wdpmr" || model=="ewdpm"){
psiT=list(psiT=result$psiT)
fitresult=c(fitresult,psiT)
}
}
hT=list(hT=result$hT)
fitresult=c(fitresult,hT)
if(model!="pm" && trace==T){
mumat=list(mumat=result$mumat)
fitresult= c(fitresult,mumat)
}
return (fitresult)
}