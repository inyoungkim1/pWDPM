// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::vec  cpp1(double inputn, double inputp, double inputQ, double inputtype){
     int n=floor(inputn);
     int p=floor(inputp);
     int type=floor(inputtype);
     vec zv(1);
     double Q=1;
     mat Xmat(n,p),Candimat(1,1),probmat(1,1),mumat(1,1); 
     double alpha=0.05,M=2,Mu0=-0.2,Sz=3,Sn=0.1,fi=0.96,deflat=1,v1=0,W=0.01,temp=0,probsuc=0;
     v1=(1-alpha)*Sz;
     double hp=R::rnorm(0,1);
     vec h(n);
     for(int t = 0; t<n; ++t){
     if(t==0){
     h[0]=hp*fi+R::rnorm(0,sqrt(Sn));
     }
     if(t>0){
     h[t]=h[t-1]*fi+R::rnorm(0,sqrt(Sn)); 
     }
     }
     if(type>2){
     if(type==3){  
     Q=pow(3,p);
     }
     if(type==4){
     Q=inputQ;
     }
     Candimat.ones(Q,p);
     if(type==3){
     double tempz1=0,tempz2=0;
     int tempz3=0;
     vec tempv6(3);
     tempv6[0]=-3;
     tempv6[1]=0;
     tempv6[2]=3;
     for(int j=0; j<p; j++){
     tempz1=pow(3,j);  
     for(int t=0; t<Q; t++){
     tempz2=t;
     tempz3=floor(tempz2/tempz1);
     Candimat(t,j)=tempv6[tempz3%3];
     }
     }
     }
     if(type==4){
     for(int j = 0; j<p; ++j){
     for(int q = 0; q<Q; ++q){
     Candimat(q,j)=R::runif(-5,5);   
     }
     }
     }
     }
     if(type>1){
     probmat.ones(150,Q);
     mumat.ones(150,Q);
     for(int q = 0; q<Q; ++q){
     deflat=1;
     for(int l = 0; l<150; ++l){
     if (l==0){
     probmat(l,q)=R::rbeta(1,M);
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=1-probmat(l,q);
     }
     if (l>0){
     temp=R::rbeta(1,M);
     probmat(l,q)=deflat*temp;
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=deflat*(1-temp);
     }
     }
     }     
     }  
     vec y(n),distance(Q),a1(1),a2(1),index=linspace<vec>(1,Q,Q),prob(1),candiv(1),zs(n),mu(n),Z(n),
     tempv(1),stempv(1),sortindex(Q),findv(1),qmat(Q),pmat(Q),wmat(Q),qvec(Q),qvec1(Q);
     y.fill(0);
     qvec.fill(1);
     qvec1.fill(1);
     if(type>2){
     for(int t = 0; t<n; ++t){
     if(t==0){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(0,lag)=R::rnorm(0,0.1);
     }
     }
     if(t>0 && t<p){
     for (int lag = 0;  lag< t; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     for (int lag = t;  lag< p; ++lag){
     Xmat(t,lag)=Xmat(t-1,lag-1);
     }
     }
     if(t>=p){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     }
     a1=conv_to< colvec >::from(Xmat.submat(t,0,t,p-1));
     for(int j = 0; j<Q; ++j){
     a2=conv_to< colvec >::from(Candimat.submat(j,0,j,p-1)); 
     distance[j]=sqrt(sum((a1-a2)%(a1-a2)));
     distance[j]=exp(2*log(distance[j]));
     }
     tempv=distance+index/100000000;
     stempv=sort(tempv);
     for(int j = 0; j<Q; ++j){
     findv=index.elem(find(tempv==stempv[j]));
     sortindex[j]=findv[0];
     }
     tempv=1/stempv;
     qvec=qvec1;
     qmat[0]=1;
     for(int j = 0; j<Q; ++j){
     pmat[j]=tempv[j]/sum(tempv.subvec(j,Q-1));
     qvec[j]=1-pmat[j];
     if(j>0){
     qmat[j]=prod(qvec.subvec(0,j-1));
     }
     wmat[j]=(1-0.7*(1-pmat[j]))*qmat[j]*exp(j*log(0.7));
     }
     prob=wmat;
     candiv=sortindex;  
     for (int q = 0;  q< Q; ++q) {
     if (q==0) {
     deflat=1.0;
     }
     if (q>0) {
     deflat=deflat-prob[q-1];
     }
     probsuc=prob[q]/deflat;
     if (R::runif(0,1)<probsuc) {
     Z[t]=candiv[q];
     break;}
     }     
     prob=conv_to< colvec >::from(probmat.submat(0,Z[t]-1,149,Z[t]-1));
     candiv=conv_to< colvec >::from(mumat.submat(0,Z[t]-1,149,Z[t]-1));    
     for (int j = 0;  j< 150; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(alpha*Sz));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<W){
     y[t]=0;
     }
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     } 
     }
     }
     if(type==2){
     for(int t = 0; t<n; ++t){ 
     prob=conv_to< colvec >::from(probmat.submat(0,0,149,0));
     candiv=conv_to< colvec >::from(mumat.submat(0,0,149,0));    
     for (int j = 0;  j< 150; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(alpha*Sz));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<W){
     y[t]=0;
     }
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     }
     }
     }
     if(type==1){
     vec mixmean=NumericVector::create(-10.12999, -3.97281, -8.56686, 
     2.77786, 0.61942, 1.79518, -1.08819);
     mixmean=mixmean-1.2704;
     vec mixprob=NumericVector::create(0.0073, 0.10556, 0.00002, 
     0.04395, 0.34001, 0.24566, 0.25750);
     vec mixsig=NumericVector::create(5.79596, 2.61369,5.17950, 
     0.16735 ,0.64009 ,0.34023, 1.26261);
     prob=mixprob;
     candiv=mixmean;  
     vec varv(n);
     for(int t = 0; t<n; ++t){   
     for (int j = 0;  j< 7; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     varv[t]=mixsig[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(varv[t]));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     } 
     }
     }
     zv=y;
     return zv;
     }

// [[Rcpp::export]]
arma::vec  cpp2(double inputn, double inputp, double inputQ, double inputtype, double method){
     int n=floor(inputn);
     int p=floor(inputp);
     int type=floor(inputtype);
     double Q=1;
     mat Xmat(n,p),Candimat(1,1),probmat(1,1),mumat(1,1); 
     double alpha=0.05,M=2,Mu0=-0.2,Sz=3,Sn=0.1,fi=0.96,deflat=1,v1=0,W=0.01,temp=0,probsuc=0;
     v1=(1-alpha)*Sz;
     double hp=R::rnorm(0,1);
     vec h(n);
     for(int t = 0; t<n; ++t){
     if(t==0){
     h[0]=hp*fi+R::rnorm(0,sqrt(Sn));
     }
     if(t>0){
     h[t]=h[t-1]*fi+R::rnorm(0,sqrt(Sn)); 
     }
     }
     if(type==1){  
     Q=pow(3,p);
     }
     if(type==2){
     Q=inputQ;
     }
     
     Candimat.ones(Q,p);
     if(type==1){
     double tempz1=0,tempz2=0;
     int tempz3=0;
     vec tempv6(3);
     tempv6[0]=-3;
     tempv6[1]=0;
     tempv6[2]=3;
     for(int j=0; j<p; j++){
     tempz1=pow(3,j);  
     for(int t=0; t<Q; t++){
     tempz2=t;
     tempz3=floor(tempz2/tempz1);
     Candimat(t,j)=tempv6[tempz3%3];
     }
     }
     }
     if(type==2){
     for(int j = 0; j<p; ++j){
     for(int q = 0; q<Q; ++q){
     Candimat(q,j)=R::runif(-5,5);   
     }
     }
     }
     probmat.ones(150,Q);
     mumat.ones(150,Q);
     for(int q = 0; q<Q; ++q){
     deflat=1;
     for(int l = 0; l<150; ++l){
     if (l==0){
     probmat(l,q)=R::rbeta(1,M);
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=1-probmat(l,q);
     }
     if (l>0){
     temp=R::rbeta(1,M);
     probmat(l,q)=deflat*temp;
     mumat(l,q)=R::rnorm(Mu0,sqrt(v1));
     deflat=deflat*(1-temp);
     }
     }
     }
     vec y(n),distance(Q),a1(1),a2(1),index=linspace<vec>(1,Q,Q),prob(1),candiv(1),zs(n),mu(n),Z(n),
     tempv(1),stempv(1),sortindex(Q),findv(1),qmat(Q),pmat(Q),wmat(Q),qvec(Q),qvec1(Q);
     y.fill(0);
     qvec.fill(1);
     qvec1.fill(1);
     vec psiv(Q),gamv(Q);
     psiv.fill(3);
     gamv.fill(1);
     for(int j = 0; j<Q; ++j){
     if(method<2){
     gamv[j]=R::rgamma(0.1,0.01);
     }
     if(method==3){
     psiv[j]=R::rgamma(2,1);
     }
     if(method==4){
     psiv[j]=R::rexp(2.42);
     }
     if(method==5){
     psiv[j]=R::rcauchy(0,0.003);
     }
     }
     double tempdist;
     for(int t = 0; t<n; ++t){
     if(t==0){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(0,lag)=R::rnorm(0,0.1);
     }
     }
     if(t>0 && t<p){
     for (int lag = 0;  lag< t; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     for (int lag = t;  lag< p; ++lag){
     Xmat(t,lag)=Xmat(t-1,lag-1);
     }
     }
     if(t>=p){
     for (int lag = 0;  lag< p; ++lag){
     Xmat(t,lag)=y[t-lag-1];
     }
     }
     a1=conv_to< colvec >::from(Xmat.submat(t,0,t,p-1));
     for(int j = 0; j<Q; ++j){ 
     a2=conv_to< colvec >::from(Candimat.submat(j,0,j,p-1)); 
     tempdist=sum((a1-a2)%(a1-a2));
     wmat[j]=log(gamv[j])-psiv[j]*tempdist;
     }
     wmat=wmat-max(wmat);
     wmat=exp(wmat);
     prob=wmat/sum(wmat);
     for (int q = 0;  q< Q; ++q) {
     if (q==0) {
     deflat=1.0;
     }
     if (q>0) {
     deflat=deflat-prob[q-1];
     }
     probsuc=prob[q]/deflat;
     if (R::runif(0,1)<probsuc) {
     Z[t]=q+1;
     break;}
     }  
     prob=conv_to< colvec >::from(probmat.submat(0,Z[t]-1,149,Z[t]-1));
     candiv=conv_to< colvec >::from(mumat.submat(0,Z[t]-1,149,Z[t]-1));    
     for (int j = 0;  j< 150; ++j) {
     if (j==0) {
     deflat=1.0;
     }
     if (j>0) {
     deflat=deflat-prob[j-1];
     }
     probsuc=prob[j]/deflat;
     if (R::runif(0,1)<probsuc) {
     mu[t]=candiv[j];
     break;}
     }
     zs[t]=R::rnorm(mu[t],sqrt(alpha*Sz));
     y[t]=h[t]+zs[t];  
     y[t]=exp(y[t]/2);
     if(R::runif(0,1)<W){
     y[t]=0;
     }
     if(R::runif(0,1)<0.5){
     y[t]=-y[t];
     } 
     }
     return y;
     }

// [[Rcpp::export]]
Rcpp::List cpp3(arma::vec y,double maxiter,double Nb){ 
int n=y.n_elem;
double logc=-10;
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
vec S(n);
vec seven=linspace<vec>(1,7,7);
for(int i = 0; i<n; ++i){
S[i]=seven[i%7];
}
double fi=0.86,sigmaeta=0.2,mu0=0.3;
vec h(n),m(n),C(n);
vec mixprob(7),mixmean(7),mixvar(7),prob(7);
mixprob[0]=0.0073;
mixprob[1]=0.10556;
mixprob[2]=0.00002;
mixprob[3]=0.04395;
mixprob[4]=0.34001;
mixprob[5]=0.24566;
mixprob[6]=0.2575;
mixmean[0]=-10.12999;
mixmean[1]=-3.97281;
mixmean[2]=-8.56686;
mixmean[3]=2.77786;
mixmean[4]=0.61942;
mixmean[5]=1.79518;
mixmean[6]=-1.08819;
mixmean=mixmean-1.2704;
mixvar[0]=5.79596;
mixvar[1]=2.61369;
mixvar[2]=5.17950;
mixvar[3]=0.16735;
mixvar[4]=0.64009;
mixvar[5]=0.34023;
mixvar[6]=1.26261;
prob.fill(0);
double a=1,R=1,f=1,Q=1,A=1,e=1,V=1,C0=sigmaeta/(1-fi*fi),mean1=1,var1=1,state=0;
double deflat=1,probsuc=1,temp=0,lb=0,ub=0,newfi=0;
double logratio=0;
vec tempv3(1),tempv4(1);
vec hstar(1);
int para=0;
vec vart(n),mt(n);
double temp2=0,temp3=0;
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),Mu0T(Nm);
Nb=Nb-1;
double dt=0;
double temps=0;
int draw=0;
vec hT(n);
hT.fill(0);
double reptime=0;
for(int iter = 0; iter < maxiter; ++iter){
dt=mu0*(1-fi);
state=S[0];
V=mixvar[state-1];
R=fi*fi*C0+sigmaeta;
f=mixmean[state-1];
Q=R+V;
A=R/Q;
e=y[0]-f;
m[0]=A*e;
C[0]=R-A*A*Q;
for(int i = 1; i <n; ++i){
state=S[i];
V=mixvar[state-1];
a=fi*m[i-1]+dt;
R=fi*fi*C[i-1]+sigmaeta;
f=a+mixmean[state-1];
Q=R+V;
A=R/Q;
e=y[i]-f;
m[i]=a+A*e;
C[i]=R-A*A*Q;
}
h[n-1]=R::rnorm(m[n-1],sqrt(C[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/C[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+m[n-i-1]/C[n-i-1]-dt*fi/sigmaeta);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int i = 0; i<n; ++i){
for (int j = 0;  j< 7; ++j){
temp=y[i]-h[i];
prob[j]=mixprob[j]*R::dnorm(temp,mixmean[j],sqrt(mixvar[j]),0);
}
prob=prob/sum(prob);
for (int j = 0;  j< 7; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
S[i]=j+1;
break;}
}
}
para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h-mu0;
for (int i = 0; i < n; ++i){
vart[i]=mixvar[S[i]-1];
mt[i]=(y[i]-hstar[i]-mixmean[S[i]-1])/vart[i];
}
var1=1/sum(1/vart);
mean1=sum(mt)*var1;
mu0=R::rnorm(mean1,sqrt(var1));
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv4=hstar.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.99,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.99,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.99,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.99,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025+hstar[0]*hstar[0]*(1-fi*fi)/2;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(hstar[i]-hstar[i-1]*fi)*(hstar[i]-hstar[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
h=hstar+mu0;
}
if(para==0){
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*h[0]+(1-fi)/sigmaeta*sum(h.subvec(1,n-1)-fi*h.subvec(0,n-2)));
mu0=R::rnorm(mean1,sqrt(var1));
logratio=0;
tempv3=h.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=h.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.99,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.99,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.99,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.99,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025+(h[0]-mu0)*(h[0]-mu0)*(1-fi*fi)/2;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((h[i]-mu0)-(h[i-1]-mu0)*fi)*((h[i]-mu0)-(h[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
}
if(iter%10==0){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
Mu0T[draw]=mu0;
draw=draw+1;
if(iter>Nb){
hT=hT+h;
reptime=reptime+1;
}
}
}
Nb=floor(Nb/10.0);
hT=hT/reptime;
double Sn=mean(SnT.subvec(Nb,Nm-1));
fi=mean(fiT.subvec(Nb,Nm-1));
double mu1=mean(Mu0T.subvec(Nb,Nm-1));
double score=0;
vec LPS(n);
LPS.fill(0);
for(int i = 0; i<10000; ++i){
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(mu1,sqrt(Sn/(1-fi*fi)));
}
if(t>0){
h[t]=mu1+fi*(h[t-1]-mu1)+R::rnorm(0,sqrt(Sn));
}
for (int j = 0;  j< 7; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-mixprob[j-1];
}
probsuc=mixprob[j]/deflat;
if (R::runif(0,1)<probsuc) {
temps=j;
break;}
}
LPS[t]=LPS[t]+R::dnorm(y[t]-h[t],mixmean[temps],sqrt(mixvar[temps]),0);
}
}
LPS=LPS/10000;
LPS=LPS.subvec(10,n-1);
score=mean(log(LPS));
return List::create(Named("fi")=fi,Named("Sn")=Sn,Named("Mu")=mu1,
Named("fiT")=fiT,Named("SnT")=SnT,Named("MuT")=Mu0T,Named("LPS")=score,Named("hT")=hT);
}

// [[Rcpp::export]]
Rcpp::List cpp4(arma::vec y,int p,int Q,double maxiter,double Nb){ 
double alpha=0.05,Sz0=0.05,psiub=4;
double logc=-20;
double minp=exp(-50);
double rob1=0.00001;
if(p<3){
rob1=0.001;
}
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
y=y.subvec(p,n1-1);
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
double tempz1=0,tempz2=0,Qmax=pow(3,p);
int tempz3=0;
vec tempv6(3);
tempv6[0]=-3;
tempv6[1]=0;
tempv6[2]=3;
mat candimat(Qmax,p);
for(int j=0; j<p; j++){
tempz1=pow(3,j);  
for(int t=0; t<Qmax; t++){
tempz2=t;
tempz3=floor(tempz2/tempz1);
candimat(t,j)=tempv6[tempz3%3];
}
}
vec tempv4(1),tempv5(1); 
double psi=0.8;
mat distance(n,Qmax),distance1(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)))+rob1;
distance(i,j)=log(distance1(i,j));
}
check[j]=mean(distance1.col(j));
}
if(Q>Qmax){
Q=Qmax;
}
if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
for(int i = 0; i<n; ++i){
distance(i,q)=distance(i,loc);
}
}
}
distance=distance.submat(0,0,n-1,Q-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}
double tempj=1.0;
vec tempv9(1);
double rv=0.7;
mat qmat(n,Q),pmat(n,Q),wmat(n,Q),wmatn(n,Q),odwmat(n,Q),odwmatn(n,Q),marwmat(n,Q),marwmatn(n,Q);
vec qvec(Q),qvec1(Q);
qvec.fill(1);
qvec1.fill(1);
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
qvec=qvec1;
qmat(i,0)=1;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
tempv9=tempv5.subvec(j,Q-1);
tempv9=tempv9-tempv9[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
pmat(i,j)=tempv9[0];
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmat(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmat(i,loc-1)=wmat(i,j);
marwmat(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
if(tempv9[0]==1){
for(int pw = j+1; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmat(i,pw)=0;
odwmat(i,loc-1)=0;
marwmat(i,loc-1)=0;
}
break;
}
}
}
double fi=0.86,sigmaeta=0.08;
vec h(n);
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(0,1)*fi+R::rnorm(0,sqrt(sigmaeta));
}
if(t>0){
h[t]=h[t-1]*fi+R::rnorm(0,sqrt(sigmaeta)); 
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
Sod=unique(S);
int k=Sod.n_elem;
vec mu(n),mustar(k+1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0.4,0.1);
}
mustar[0]=logc;
for (int i = 0; i<n; ++i){
mu[i]=mustar[S[i]];
}
double sigmasq=2.5,W=0.01,M=4,mu0=0;
h=ones<vec>(n);
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
rv=1;
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,V=Sz0,C0=sigmaeta/(1-fi*fi),
mean1=0,var1=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double newpsi=1,logratio=0;
vec tempv(1),gamv(n);
double a0=1,d0=0.2,templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
double newfi=0,lb=0,ub=0;
int para=0;
vec x=linspace<vec>(-30,20,1000);
int nx=1000;
mat predmat(nx,Q);
predmat.fill(0);
double reptime=0;
vec tempv7(Q),tempv8(Q);
int draw=0;
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),SzT(Nm),MT(Nm),WT(Nm),Mu0T(Nm),psiT(Nm);
Nb=Nb-1;
double tempsum=0,tempbr=0;
vec hT(n);
hT.fill(0);
mat mumat(Nm,n);




for(int iter = 0; iter < maxiter; ++iter){
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(logc);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
temp2=a0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=d0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
temp=alpha*sigmasq;
if(S[0]>0){
V=temp;
}
temp2=fi*fi*C0+sigmaeta;
temp3=mu[0];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[0]-temp3;
tempv1[0]=temp5*temp6;
tempv2[0]=temp2-temp5*temp5*temp4;
for(int i = 1; i < n; ++i){
V=Sz0;
if(S[i]>0) {
V=temp;
}
temp1=fi*tempv1[i-1];
temp2=fi*fi*tempv2[i-1]+sigmaeta;
temp3=temp1+mu[i];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[i]-temp3;
tempv1[i]=temp1+temp5*temp6;
tempv2[i]=temp2-temp5*temp5*temp4;
}
h[n-1]=R::rnorm(tempv1[n-1],sqrt(tempv2[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/tempv2[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+tempv1[n-i-1]/tempv2[n-i-1]);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=W*exp(-0.5*(y[position]-h[position]-logc)*(y[position]-h[position]-logc)/Sz0)/sqrt(Sz0);
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=(1-W)*count[j]/sqrt(alpha*sigmasq)/(n2-1+M)*
exp(-0.5*(y[position]-h[position]-tempmustar[j])*(y[position]-h[position]-tempmustar[j])/alpha/sigmasq);
probsum=probsum+prob[j];
}
}
prob[lmu]=(1-W)*M/sqrt(sigmasq)/(n2-1+M)*
exp(-(y[position]-h[position]-mu0)*(y[position]-h[position]-mu0)/sigmasq/2);  
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((y[position]-h[position])/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
if(iter>Nb && iter%10==9){
hT=hT+h;
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<nx; ++i){
for (int q = 0;  q< Q; ++q){
temp6=W*R::dnorm(x[i],logc,sqrt(Sz0),0);
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+(1-W)*basecountcluster(l,q)/(tempv8[q]+M)*
R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+(1-W)*M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
reptime=reptime+1;
}
if(Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
tempv=conv_to< colvec >::from(sortindex.submat(i,0,i,(Q-1)));
locv=index.elem(find(tempv==Z[i]));
loc=locv[0];
temp1=1-pmat(i,loc-1);
if(loc<Q){
rv=R::qbeta(R::runif(0,1)*R::pbeta(temp1,loc,2,1,0),loc,2,1,0)/temp1;
}
if(loc==Q){
rv=exp(log(R::runif(0,1))/Q);
}
gamv[i]=rv;
if(rv==0){
gamv[i]=0.00001;
}
}
psiub=log(4.0/psi);
newpsi=R::qnorm(R::runif(0,R::pnorm(psiub,-0.001,0.05,1,0)),-0.001,0.05,1,0);
newpsi=psi*exp(newpsi);
logratio=0;
if(n<601){
logratio=-R::dnorm(log(psi),0.7,0.05,1)+R::dnorm(log(newpsi),0.7,0.05,1);
}
if(n>600){
logratio=-R::dnorm(log(psi),0.08,0.05,1)+R::dnorm(log(newpsi),0.08,0.05,1);
}
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
rv=gamv[i];
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmatn(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmatn(i,loc-1)=wmatn(i,j);
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmatn(i,pw)=0;
odwmatn(i,loc-1)=0;
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
wmat=wmatn;
odwmat=odwmatn;
marwmat=marwmatn;
} 
}
mustar=S0v.subvec(0,k);
mustar[0]=logc;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(y[i]-h[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
mu[i]=mustar[S[i]];
if(iter%10==0){
mumat(draw,i)=mu[i];
}
}

para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h+mu0;
mustar1=mustar.subvec(1,k)-mu0;
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*hstar[0]+(1-fi)/sigmaeta*sum(hstar.subvec(1,n-1)-fi*hstar.subvec(0,n-2)));
temp6=mu0;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum(mustar1%mustar1)/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=hstar.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi)*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
h=hstar-mu0;
mustar.subvec(1,k)=mustar1+mu0;
for(int i = 0; i<n; ++i){
if(S[i]>0){
mu[i]=mu[i]-temp6+mu0;
}
}
}
if(para==0){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=h.subvec(1,n-1);
tempv4=h.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-h[0]*h[0]*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+h[0]*h[0]*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(h[i]-h[i-1]*fi)*(h[i]-h[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
}
temp2=n-n2;
temp2=temp2+0.1;
temp3=n2;
temp3=temp3+0.9;
W=R::rbeta(temp2,temp3);
if(iter%10==0){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
SzT[draw]=sigmasq;
MT[draw]=M;
WT[draw]=W;
Mu0T[draw]=mu0;
psiT[draw]=psi;
draw=draw+1;
}
}
for(int iter = 0; iter < Nm; ++iter){  
newpsi=R::qnorm(R::runif(0,R::pnorm(log(5.0/psi),-0.5,1,1,0)),-0.5,1,1,0);
newpsi=psi*exp(newpsi);
logratio=5*log(newpsi)-5*log(psi)+2*(psi-newpsi);
for(int i = 0; i<1; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
marwmat=marwmatn;
}
psiT[iter]=psi;
}
predmat=predmat/reptime;
hT=hT/reptime;
hT=exp(hT/2.0);
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,999))+sum(tempv.subvec(0,998));
temp=temp/2.0;
temp=temp*50.0/999;
tempv=tempv/temp;
predmat.col(q)=tempv;
}
Nb=floor(Nb/10.0);
psi=mean(psiT.subvec(Nb,Nm-1));
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
temp=j;
odwmat(i,loc-1)=qmat(i,j)*(1/(temp+1.0)-(1-pmat(i,j))/(temp+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
odwmat(i,loc-1)=0;
}
}
}
temp1=mean(SnT.subvec(Nb,Nm-1));
temp2=mean(fiT.subvec(Nb,Nm-1));
temp3=0;
vec pred0(1),pred1(1),predv(1),odwv(1),LPS(n);
LPS.fill(0);
double err=0;
double step=50.0/999;
for(int i = 0; i<10000; ++i){
if(err>1){
break;
}
for(int t = 0; t<n; ++t){
if(err>1){
break;
}
if(t==0){
h[0]=R::rnorm(0,sqrt(temp1/(1-temp2*temp2)));
}
if(t>0){
h[t]=temp2*h[t-1]+R::rnorm(0,sqrt(temp1));
}
temp4=y[t]-h[t];
if(temp4<x[0]){
err=1000;
}
if(temp4>x[nx-1]){
err=1000;
}
if(err==0){
loc=floor((temp4-x[0])/step);
pred0=conv_to< colvec >::from(predmat.submat(loc,0,loc,Q-1));
pred1=conv_to< colvec >::from(predmat.submat(loc+1,0,loc+1,Q-1));
predv=(temp4-x[loc])/step*(pred1-pred0)+pred0;
odwv=conv_to< colvec >::from(odwmat.submat(t,0,t,Q-1));
LPS[t]=LPS[t]+sum(predv%odwv);
}
}
}
if(err<1){
LPS=LPS/10000;
LPS=LPS.subvec(10-p,n-1);
temp3=mean(log(LPS));
}
if(err>1){
temp3=datum::nan;
}
sigmasq=mean(SzT.subvec(Nb,Nm-1));
W=mean(WT.subvec(Nb,Nm-1));
mu0=mean(Mu0T.subvec(Nb,Nm-1));
M=mean(MT.subvec(Nb,Nm-1));
return List::create(
Named("fi")=temp2,Named("Sn")=temp1,Named("mu0")=mu0,Named("Sz")=sigmasq,Named("M")=M,Named("W")=W,Named("psi")=psi,
Named("fiT")=fiT, Named("SnT")=SnT,Named("Mu0T")=Mu0T,Named("SzT")=SzT,Named("MT")=MT,Named("WT")=WT,Named("psiT")=psiT,
Named("LPS")=temp3,Named("hT")=hT,Named("mumat")=mumat);
} 

// [[Rcpp::export]]
Rcpp::List cpp5(arma::vec y,int p,int Q,double maxiter,double Nb){ 
double alpha=0.05,Sz0=0.05,psiub=4;
double logc=-20;
double minp=exp(-50);
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
y=y.subvec(p,n1-1);
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
double Qmax=Q;
mat candimat(Qmax,p);
for(int j=0; j<p; j++){  
for(int t=0; t<Qmax; t++){
candimat(t,j)=R::runif(-5,5);
}
}
vec tempv4(1),tempv5(1); 
double psi=0.8;
mat distance(n,Qmax),distance1(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)));
distance(i,j)=log(distance1(i,j));
}
check[j]=mean(distance1.col(j));
}
if(Q>Qmax){
Q=Qmax;
}
if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
for(int i = 0; i<n; ++i){
distance(i,q)=distance(i,loc);
}
}
}
distance=distance.submat(0,0,n-1,Q-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}
double tempj=1.0;
vec tempv9(1);
double rv=0.7;
mat qmat(n,Q),pmat(n,Q),wmat(n,Q),wmatn(n,Q),odwmat(n,Q),odwmatn(n,Q),marwmat(n,Q),marwmatn(n,Q);
vec qvec(Q),qvec1(Q);
qvec.fill(1);
qvec1.fill(1);
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
qvec=qvec1;
qmat(i,0)=1;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
tempv9=tempv5.subvec(j,Q-1);
tempv9=tempv9-tempv9[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
pmat(i,j)=tempv9[0];
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmat(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmat(i,loc-1)=wmat(i,j);
marwmat(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
if(tempv9[0]==1){
for(int pw = j+1; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmat(i,pw)=0;
odwmat(i,loc-1)=0;
marwmat(i,loc-1)=0;
}
break;
}
}
}
double fi=0.86,sigmaeta=0.08;
vec h(n);
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(0,1)*fi+R::rnorm(0,sqrt(sigmaeta));
}
if(t>0){
h[t]=h[t-1]*fi+R::rnorm(0,sqrt(sigmaeta)); 
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
Sod=unique(S);
int k=Sod.n_elem;
vec mu(n),mustar(k+1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0.4,0.1);
}
mustar[0]=logc;
for (int i = 0; i<n; ++i){
mu[i]=mustar[S[i]];
}
double sigmasq=2.5,W=0.01,M=4,mu0=0;
h=ones<vec>(n);
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
rv=1;
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,V=Sz0,C0=sigmaeta/(1-fi*fi),
mean1=0,var1=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double newpsi=1,logratio=0;
vec tempv(1),gamv(n);
double a0=1,d0=0.2,templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
double newfi=0,lb=0,ub=0;
int para=0;
vec x=linspace<vec>(-30,20,1000);
int nx=1000;
mat predmat(nx,Q);
predmat.fill(0);
double reptime=0;
vec tempv7(Q),tempv8(Q);
int draw=0;
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),SzT(Nm),MT(Nm),WT(Nm),Mu0T(Nm),psiT(Nm);
Nb=Nb-1;
double tempsum=0,tempbr=0;
vec hT(n);
hT.fill(0);
double stp=0.2;
mat mumat(Nm,n);


for(int iter = 0; iter < maxiter; ++iter){
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(logc);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
temp2=a0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=d0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
temp=alpha*sigmasq;
if(S[0]>0){
V=temp;
}
temp2=fi*fi*C0+sigmaeta;
temp3=mu[0];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[0]-temp3;
tempv1[0]=temp5*temp6;
tempv2[0]=temp2-temp5*temp5*temp4;
for(int i = 1; i < n; ++i){
V=Sz0;
if(S[i]>0) {
V=temp;
}
temp1=fi*tempv1[i-1];
temp2=fi*fi*tempv2[i-1]+sigmaeta;
temp3=temp1+mu[i];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[i]-temp3;
tempv1[i]=temp1+temp5*temp6;
tempv2[i]=temp2-temp5*temp5*temp4;
}
h[n-1]=R::rnorm(tempv1[n-1],sqrt(tempv2[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/tempv2[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+tempv1[n-i-1]/tempv2[n-i-1]);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=W*exp(-0.5*(y[position]-h[position]-logc)*(y[position]-h[position]-logc)/Sz0)/sqrt(Sz0);
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=(1-W)*count[j]/sqrt(alpha*sigmasq)/(n2-1+M)*
exp(-0.5*(y[position]-h[position]-tempmustar[j])*(y[position]-h[position]-tempmustar[j])/alpha/sigmasq);
probsum=probsum+prob[j];
}
}
prob[lmu]=(1-W)*M/sqrt(sigmasq)/(n2-1+M)*
exp(-(y[position]-h[position]-mu0)*(y[position]-h[position]-mu0)/sigmasq/2);  
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((y[position]-h[position])/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
if(iter>Nb && iter%10==9){
hT=hT+h;
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<nx; ++i){
for (int q = 0;  q< Q; ++q){
temp6=W*R::dnorm(x[i],logc,sqrt(Sz0),0);
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+(1-W)*basecountcluster(l,q)/(tempv8[q]+M)*
R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+(1-W)*M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
reptime=reptime+1;
}
if(Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(odwmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
tempv=conv_to< colvec >::from(sortindex.submat(i,0,i,(Q-1)));
locv=index.elem(find(tempv==Z[i]));
loc=locv[0];
temp1=1-pmat(i,loc-1);
if(loc<Q){
rv=R::qbeta(R::runif(0,1)*R::pbeta(temp1,loc,2,1,0),loc,2,1,0)/temp1;
}
if(loc==Q){
rv=exp(log(R::runif(0,1))/Q);
}
gamv[i]=rv;
if(rv==0){
gamv[i]=0.00001;
}
}
newpsi=R::qnorm(R::runif(R::pnorm(0,psi,stp,1,0),R::pnorm(psiub,psi,stp,1,0)),psi,stp,1,0);
logratio=0;
logratio=logratio+3*log(newpsi)-3*log(psi)+2*(psi-newpsi);;
lb=R::pnorm(0,psi,stp,1,0);
ub=R::pnorm(psiub,psi,stp,1,0);
logratio=logratio+log(ub-lb);
lb=R::pnorm(0,newpsi,stp,1,0);
ub=R::pnorm(psiub,newpsi,stp,1,0);
logratio=logratio-log(ub-lb);
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
rv=gamv[i];
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
wmatn(i,j)=(1-rv*(1-pmat(i,j)))*qmat(i,j)*exp(j*log(rv));
odwmatn(i,loc-1)=wmatn(i,j);
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
wmatn(i,pw)=0;
odwmatn(i,loc-1)=0;
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
wmat=wmatn;
odwmat=odwmatn;
marwmat=marwmatn;
} 
}
mustar=S0v.subvec(0,k);
mustar[0]=logc;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(y[i]-h[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
mu[i]=mustar[S[i]];
if(iter%10==0){
mumat(draw,i)=mu[i];
}
}
para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h+mu0;
mustar1=mustar.subvec(1,k)-mu0;
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*hstar[0]+(1-fi)/sigmaeta*sum(hstar.subvec(1,n-1)-fi*hstar.subvec(0,n-2)));
temp6=mu0;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum(mustar1%mustar1)/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=hstar.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi)*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
h=hstar-mu0;
mustar.subvec(1,k)=mustar1+mu0;
for(int i = 0; i<n; ++i){
if(S[i]>0){
mu[i]=mu[i]-temp6+mu0;
}
}
}
if(para==0){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=h.subvec(1,n-1);
tempv4=h.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-h[0]*h[0]*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+h[0]*h[0]*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(h[i]-h[i-1]*fi)*(h[i]-h[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
}
temp2=n-n2;
temp2=temp2+0.1;
temp3=n2;
temp3=temp3+0.9;
W=R::rbeta(temp2,temp3);
if(iter%10==0){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
SzT[draw]=sigmasq;
MT[draw]=M;
WT[draw]=W;
Mu0T[draw]=mu0;
psiT[draw]=psi;
draw=draw+1;
}
}
for(int iter = 0; iter < Nm; ++iter){  
newpsi=R::qnorm(R::runif(0,R::pnorm(log(5.0/psi),-0.5,1,1,0)),-0.5,1,1,0);
newpsi=psi*exp(newpsi);
logratio=5*log(newpsi)-5*log(psi)+2*(psi-newpsi);
for(int i = 0; i<1; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-newpsi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
tempj=j;
marwmatn(i,loc-1)=qmat(i,j)*(1/(tempj+1.0)-(1-pmat(i,j))/(tempj+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
marwmatn(i,loc-1)=0;
}
}
logratio=logratio+log(marwmatn(i,Z[i]-1))-log(marwmat(i,Z[i]-1));
}
if(log(R::runif(0,1))<logratio){
psi=newpsi;
marwmat=marwmatn;
}
psiT[iter]=psi;
}
predmat=predmat/reptime;
hT=hT/reptime;
hT=exp(hT/2.0);
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,999))+sum(tempv.subvec(0,998));
temp=temp/2.0;
temp=temp*50.0/999;
tempv=tempv/temp;
predmat.col(q)=tempv;
}
Nb=floor(Nb/10.0);
psi=mean(psiT.subvec(Nb,Nm-1));
for(int i = 0; i<n; ++i){
tempv5=conv_to< colvec >::from(sortdistance.submat(i,0,i,Q-1));
tempv5=-psi*tempv5;
tempv9=tempv5-tempv5[0];
tempv9=exp(tempv9);
tempv9=tempv9/sum(tempv9);
qvec=qvec1;
qmat(i,0)=1;
tempsum=0;
tempbr=0;
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
pmat(i,j)=tempv9[j]/(1-tempsum);
qvec[j]=1-pmat(i,j);
if(j>0){
qmat(i,j)=prod(qvec.subvec(0,j-1));
}
temp=j;
odwmat(i,loc-1)=qmat(i,j)*(1/(temp+1.0)-(1-pmat(i,j))/(temp+2.0));
tempsum=tempsum+tempv9[j];
if(tempsum==1 || qmat(i,j)<minp){
tempsum=j+1;
tempbr=1;
break;
}
}
if(tempbr>0){
for(int pw = tempsum; pw<Q; ++pw){  
loc=sortindex(i,pw);
odwmat(i,loc-1)=0;
}
}
}
temp1=mean(SnT.subvec(Nb,Nm-1));
temp2=mean(fiT.subvec(Nb,Nm-1));
temp3=0;
vec pred0(1),pred1(1),predv(1),odwv(1),LPS(n);
LPS.fill(0);
double err=0;
double step=50.0/999;
for(int i = 0; i<10000; ++i){
if(err>1){
break;
}
for(int t = 0; t<n; ++t){
if(err>1){
break;
}
if(t==0){
h[0]=R::rnorm(0,sqrt(temp1/(1-temp2*temp2)));
}
if(t>0){
h[t]=temp2*h[t-1]+R::rnorm(0,sqrt(temp1));
}
temp4=y[t]-h[t];
if(temp4<x[0]){
err=1000;
}
if(temp4>x[nx-1]){
err=1000;
}
if(err==0){
loc=floor((temp4-x[0])/step);
pred0=conv_to< colvec >::from(predmat.submat(loc,0,loc,Q-1));
pred1=conv_to< colvec >::from(predmat.submat(loc+1,0,loc+1,Q-1));
predv=(temp4-x[loc])/step*(pred1-pred0)+pred0;
odwv=conv_to< colvec >::from(odwmat.submat(t,0,t,Q-1));
LPS[t]=LPS[t]+sum(predv%odwv);
}
}
}
if(err<1){
LPS=LPS/10000;
LPS=LPS.subvec(10-p,n-1);
temp3=mean(log(LPS));
}
if(err>1){
temp3=datum::nan;
}
sigmasq=mean(SzT.subvec(Nb,Nm-1));
W=mean(WT.subvec(Nb,Nm-1));
mu0=mean(Mu0T.subvec(Nb,Nm-1));
M=mean(MT.subvec(Nb,Nm-1));
return List::create(
Named("fi")=temp2,Named("Sn")=temp1,Named("mu0")=mu0,Named("Sz")=sigmasq,Named("M")=M,Named("W")=W,Named("psi")=psi,
Named("fiT")=fiT, Named("SnT")=SnT,Named("Mu0T")=Mu0T,Named("SzT")=SzT,Named("MT")=MT,Named("WT")=WT,Named("psiT")=psiT,
Named("LPS")=temp3,Named("hT")=hT,Named("mumat")=mumat);
} 

// [[Rcpp::export]]
Rcpp::List cpp6(arma::vec y,int p,int Q, int wdp,double maxiter,double Nb){   
double alpha=0.05,Sz0=0.05;
double logc=-20;
double caup=0.003*0.003;
int n1=y.n_elem;
int n=n1-p;
mat Xmat(n,p);
for(int t = 0; t<n; ++t){
for (int lag = 0;  lag< p; ++lag){
Xmat(t,lag)=y[t+p-lag-1];
}
}
y=y.subvec(p,n1-1);
for(int t = 0; t<n; ++t){
y[t]=log(y[t]*y[t]+exp(logc));
} 
double tempz1=0,tempz2=0,Qmax=pow(3,p);
int tempz3=0;
vec tempv6(3);
tempv6[0]=-3;
tempv6[1]=0;
tempv6[2]=3;
mat candimat(Qmax,p);
for(int j=0; j<p; j++){
tempz1=pow(3,j);  
for(int t=0; t<Qmax; t++){
tempz2=t;
tempz3=floor(tempz2/tempz1);
candimat(t,j)=tempv6[tempz3%3];
}
}
vec tempv4(1),tempv5(1); 
double psi=3;
mat distance(n,Qmax),distance1(n,Qmax);
vec index=linspace<vec>(1,Qmax,Qmax),check(Qmax),findv1(1);
int loc=0;
for(int j = 0; j<Qmax; ++j){
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(Xmat.submat(i,0,i,p-1));
tempv5=conv_to< colvec >::from(candimat.submat(j,0,j,p-1)); 
distance1(i,j)=sqrt(sum((tempv4-tempv5)%(tempv4-tempv5)));
}
check[j]=mean(distance1.col(j));
}
if(Q>Qmax){
Q=Qmax;
}
if(Q<Qmax){
tempv4=sort(check);
tempv4=tempv4.subvec(0,Q-1);
for(int q = 0; q<Q; ++q){
findv1=index.elem(find(check==tempv4[q]));
loc=findv1[0]-1;
for(int i = 0; i<n; ++i){
distance(i,q)=distance1(i,loc);
}
}
}
if(Q==Qmax){
distance=distance1;
}
distance=distance.submat(0,0,n-1,Q-1);
double temp=0;
mat sortdistance(n,Q),sortindex(n,Q),sortdistance1(n,Q);
vec Z(n),countz(Q),S(n),stempv(1);
countz.fill(0);
index=linspace<vec>(1,Q,Q);
for(int i = 0; i<n; ++i){
tempv4=conv_to< colvec >::from(distance.submat(i,0,i,Q-1));
tempv4=tempv4+index/100000000;
stempv=sort(tempv4);
for(int j = 0; j<Q; ++j){
sortdistance(i,j)=stempv[j];
sortdistance1(i,j)=-log(sortdistance(i,j))/psi;
temp=sortdistance(i,j);
findv1=index.elem(find(tempv4==temp));
sortindex(i,j)=findv1[0];
}
Z[i]=sortindex(i,0);
countz[Z[i]-1]=countz[Z[i]-1]+1;
S[i]=R::rpois(6);
}

mat wmat(n,Q),wmatn(n,Q);
double fi=0.86,sigmaeta=0.08;
vec h(n);
for(int t = 0; t<n; ++t){
if(t==0){
h[0]=R::rnorm(0,1)*fi+R::rnorm(0,sqrt(sigmaeta));
}
if(t>0){
h[t]=h[t-1]*fi+R::rnorm(0,sqrt(sigmaeta)); 
}
}
S=10000*Z+S;  
int od=1,len=0;
vec tempa(1),Sod(n);
uvec findv2(1);
for (int i = 0; i < n; ++i) {
if(i==0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(S[i]==0){
Sod[i]=0;
}
if(S[i]>0){
tempa=S.subvec(0,i-1);
findv2=find(tempa==S[i]);
len=findv2.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv4=Sod.subvec(0,i-1);
findv1=tempv4.elem(find(tempa==S[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
mat odwmat(n,Q);
for(int i = 0; i<n; ++i){
for(int j = 0; j<Q; ++j){
loc=sortindex(i,j);
odwmat(i,loc-1)=wmat(i,j);
}
}
Sod=unique(S);
int k=Sod.n_elem;
vec mu(n),mustar(k+1);
for(int i = 1; i<k+1; ++i){
mustar[i]=R::rnorm(0.4,0.1);
}
mustar[0]=logc;
for (int i = 0; i<n; ++i){
mu[i]=mustar[S[i]];
}
double sigmasq=2.5,W=0.01,M=4,mu0=0;
sortdistance=sortdistance1;
h=ones<vec>(n);
temp=0;
countz=ones<vec>(1);
findv1=ones<vec>(1);
loc=0;
len=0;
k=0;
od=1;
tempv4=ones<vec>(1);
tempv5=ones<vec>(1);
tempa=ones<vec>(1);
Sod=ones<vec>(n);
index=linspace<vec>(1,Q,Q);
vec tempv1(n),tempv2(n); 
double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,V=Sz0,C0=sigmaeta/(1-fi*fi),
mean1=0,var1=0;
vec countzp(1),countzus(1),locv(1),sp0(1),sp1(1),sp(1),pos(1),ord(1),maxS(1),Zuniq(1);
vec q0v(Q);
q0v.fill(0);
int zt=0,st=0,locs=0,kz=0;
uvec findv(1);
vec index1=linspace<vec>(0,n-1,n);
mat basecountcluster(n,Q),basecluster(n,Q),baseclusterS0(n,Q),
baseclusterS(n,Q),baseconfigS(n,Q),basemu(n,Q);
int position=0,n2=0,lmu=0,addnew=0,nz=0;
double probsum=0,deflat=1,probsuc=0,newmu=1;
vec newS(n),tempS(1),location(1),count(1),tempmustar(1),prob(1),
trueS(1),lamz(1),countz1(1),counts(1);
vec S0v(n+1);
S0v.fill(0);
vec tempmean(1),tempvar(1);
vec zprob(1),sz(1),scount(1),zcluster(1);
mat Q0k(Q,n),zprobmat(1,1);
Q0k.fill(0);
double sm=1e-320;
double newpsi=3.1,logratio=0;
vec tempv(1);
double a0=1,d0=0.2,templogb=0,tempinvb=1;
vec mix(1),coef(1),mixgam(1);
vec hstar(1),mustar1(1);
vec tempv3(1);
double newfi=0,lb=0,ub=0;
int para=0;
vec x=linspace<vec>(-30,20,1000);
int nx=1000;
mat predmat(nx,Q);
predmat.fill(0);
double reptime=0;
vec tempv7(Q),tempv8(Q);
int draw=0;
double logL1=0,logL2=0;
mat KS(n,Q),SS(n,Q);
vec sumv(n);
vec psiv(Q);
psiv.fill(1.2);
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
KS(t,q)=exp(-psiv[q]*distance(t,q));
}
sumv[t]=sum(KS.row(t));
}
vec sumpv(1);
int change=0;
vec psivt(Q);
psivt.fill(0);
double Nm=maxiter/10.0;
vec fiT(Nm),SnT(Nm),SzT(Nm),MT(Nm),WT(Nm),Mu0T(Nm),psiT(Nm);
Nb=Nb-1;
wmat.fill(1.0/n);
vec probv1(Q),probv2(Q);
double repsi=0;
mat mumat(Nm,n);
vec hT(n);
hT.fill(0);

for(int iter = 0; iter < maxiter; ++iter){ 
temp=0;
countz=q0v;
Zuniq=countz;
for(int i = 0; i<n; ++i){
zt=Z[i];
countz[zt-1]=countz[zt-1]+1;
}
for(int j = 0; j<Q; ++j){
if(countz[j]>0){
temp=temp+1;
Zuniq[j]=temp;
}
}
kz=temp;
countzp=q0v.subvec(0,kz-1);
countzus=countzp;
pos=countzp;
ord=countzp;
maxS=countzp;
basecountcluster.fill(0);
basecluster.fill(-1);
baseclusterS.fill(-1);
baseclusterS0.fill(-1);
baseconfigS.fill(0);
basemu.fill(logc);
for(int i = 0; i<n; ++i){
loc=Zuniq[Z[i]-1]-1;
st=S[i];
baseclusterS0(ord[loc],loc)=st;
baseconfigS(ord[loc],loc)=i+1;
if(st==0){
baseclusterS(ord[loc],loc)=0;
basecountcluster(0,loc)=basecountcluster(0,loc)+1;
basecluster(0,loc)=0;
}
if(st>0){
countzp[loc]=countzp[loc]+1;
if(ord[loc]==0){
baseclusterS(0,loc)=1;
maxS[loc]=maxS[loc]+1;
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1; 
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
}
if(ord[loc]>0){
sp0=baseclusterS0.col(loc);
sp0=sp0.subvec(0,ord[loc]-1);
findv=find(sp0==st);
len=findv.n_elem;
if(len==0){
countzus[loc]=countzus[loc]+1;
pos[loc]=pos[loc]+1;
basecluster(pos[loc],loc)=st;
basemu(pos[loc],loc)=mustar[st];
basecountcluster(pos[loc],loc)=basecountcluster(pos[loc],loc)+1;
maxS[loc]=maxS[loc]+1;
baseclusterS(ord[loc],loc)=maxS[loc];
}
if(len>0){
sp1=baseclusterS.col(loc);
sp1=sp1.subvec(0,ord[loc]-1);
findv1=sp1.elem(find(sp0==st));
baseclusterS(ord[loc],loc)=findv1[0];
sp=basecluster.col(loc);
locv=index1.elem(find(sp==st));
locs=locv[0];
basecountcluster(locs,loc)=basecountcluster(locs,loc)+1;
}
}
}
ord[loc]=ord[loc]+1;
}
countz1=ones<vec>(kz);
od=-1;
for(int i = 0; i < Q; ++i) {
if(countz[i]>0){
od=od+1;
countz1[od]=countz[i];
}
}
temp2=a0+k-kz;
mix=ones<vec>(kz+1);
coef=mix-1;
mix[0]=lgamma(temp2);
temp3=M+1;
temp4=d0;
temp4=temp4-log(R::rbeta(temp3,countz1[0]));
coef[0]=countz1[0];
coef[1]=1;
if(kz>1){
for(int i = 1; i<kz; ++i){ 
temp4=temp4-log(R::rbeta(temp3,countz1[i]));
mix[i]=lgamma(temp2+i);
tempv=coef.subvec(0,i);
coef.subvec(1,i+1)=tempv;
coef[0]=0;
coef.subvec(0,i)=coef.subvec(0,i)+tempv*countz1[i];
coef.subvec(0,i+1)=coef.subvec(0,i+1)/max(coef.subvec(0,i+1));
}
}
mix[kz]=lgamma(temp2+kz);
templogb=log(temp4);
mix[0]=mix[0]+log(coef[0])-temp2*templogb;
temp5=mix[0];
mixgam=mix; 
tempinvb=1/temp4;   
mixgam[0]=R::rgamma(temp2,tempinvb);
for(int i = 1; i<kz+1; ++i){ 
temp1=temp2+i;
mix[i]=mix[i]+log(coef[i])-temp1*templogb;
if(mix[i]>temp5){
temp5=mix[i];
}
mixgam[i]=R::rgamma(temp1,tempinvb);
}
mix=mix-temp5;
mix=exp(mix);
mix=mix/sum(mix);
for (int q = 0;  q< kz+1; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-mix[q-1];
}
probsuc=mix[q]/deflat;
if (R::runif(0,1)<probsuc) {
M=mixgam[q];
break;
}
}
temp=alpha*sigmasq;
if(S[0]>0){
V=temp;
}
temp2=fi*fi*C0+sigmaeta;
temp3=mu[0];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[0]-temp3;
tempv1[0]=temp5*temp6;
tempv2[0]=temp2-temp5*temp5*temp4;
for(int i = 1; i < n; ++i){
V=Sz0;
if(S[i]>0) {
V=temp;
}
temp1=fi*tempv1[i-1];
temp2=fi*fi*tempv2[i-1]+sigmaeta;
temp3=temp1+mu[i];
temp4=temp2+V;
temp5=temp2/temp4;
temp6=y[i]-temp3;
tempv1[i]=temp1+temp5*temp6;
tempv2[i]=temp2-temp5*temp5*temp4;
}
h[n-1]=R::rnorm(tempv1[n-1],sqrt(tempv2[n-1]));
for(int i = 1; i < n; ++i){
var1=1/(fi*fi/sigmaeta+1/tempv2[n-i-1]);
mean1=var1*(fi*h[n-i]/sigmaeta+tempv1[n-i-1]/tempv2[n-i-1]);
h[n-i-1]=R::rnorm(mean1,sqrt(var1));
}
for(int l = 0; l<kz; ++l){
nz=countz1[l];
tempS=conv_to< colvec >::from(baseclusterS.submat(0,l,nz-1,l));
location=conv_to< colvec >::from(baseconfigS.submat(0,l,nz-1,l));
count=conv_to< colvec >::from(basecountcluster.submat(0,l,countzus[l],l));
tempmustar=conv_to< colvec >::from(basemu.submat(0,l,countzus[l],l));  
lmu=countzus[l]+1;
n2=countzp[l]; 
lamz=ones<vec>(nz);
for(int i = 0; i<nz; ++i){
addnew=0;
st=tempS[i];
position=location[i]-1;
count[st]=count[st]-1;
if(st>0){
n2=n2-1;
}
if (count[st]==0 && st>0){
count.shed_row(st);
tempmustar=tempmustar.elem(find(tempmustar!=tempmustar[st]));
tempS.elem(find(tempS>st))=tempS.elem(find(tempS>st))-1;
lmu=lmu-1;
}
prob=ones<vec>(lmu+1);
prob[0]=W*exp(-0.5*(y[position]-h[position]-logc)*(y[position]-h[position]-logc)/Sz0)/sqrt(Sz0);
probsum=prob[0];
if(lmu>1){
for (int j = 1;  j< lmu; ++j) {
prob[j]=(1-W)*count[j]/sqrt(alpha*sigmasq)/(n2-1+M)*
exp(-0.5*(y[position]-h[position]-tempmustar[j])*(y[position]-h[position]-tempmustar[j])/alpha/sigmasq);
probsum=probsum+prob[j];
}
}
prob[lmu]=(1-W)*M/sqrt(sigmasq)/(n2-1+M)*
exp(-(y[position]-h[position]-mu0)*(y[position]-h[position]-mu0)/sigmasq/2);  
probsum=probsum+prob[lmu];
prob=prob/probsum;
for (int j = 0;  j< lmu+1; ++j) {
if (j==0) {
deflat=1.0;
}
if (j!=0) {
deflat=deflat-prob[j-1];
}
probsuc=prob[j]/deflat;
if (R::runif(0,1)<probsuc) {
tempS[i]=j;
break;
}
}
if(tempS[i]>0){
n2=n2+1;
if(tempS[i]>lmu-1){
addnew=1;
}
}
if(addnew==1){
var1=1/(1 / (alpha * sigmasq) + 1 / ((1 - alpha) * sigmasq));
mean1=var1*((y[position]-h[position])/alpha/sigmasq+mu0/(1-alpha)/sigmasq);
newmu=R::rnorm(mean1,sqrt(var1));
tempv4=tempmustar.subvec(0,lmu-1);
tempmustar=ones<vec>(lmu+1);
tempmustar.subvec(0,lmu-1)=tempv4;
tempmustar[lmu]=newmu;
tempv4=count;
count=ones<vec>(lmu+1);
count.subvec(0,lmu-1)=tempv4;
count[lmu]=1;
lamz[i]=newmu;
lmu=lmu+1;
}
if(addnew==0){
st=tempS[i];
count[st]=count[st]+1;
lamz[i]=tempmustar[st];
}
}
od=0;
trueS=ones<vec>(nz);
trueS.fill(1);
position=location[0]-1;
if(tempS[0]==0){
newS[position]=0; 
}  
if(tempS[0]>0){
newS[position]=1+10000*l;
od=od+1;
trueS[0]=od;
}
if(nz>1){
for(int i = 1; i<nz; ++i){
position=location[i]-1;  
if(tempS[i]==0){
newS[position]=0;
}
if(tempS[i]>0){
tempv5=lamz.subvec(0,i);
findv1=trueS.elem(find(tempv5==lamz[i]));
len=findv1.n_elem;
if(len==1){
od=od+1;
trueS[i]=od;
}
if(len>1){
trueS[i]=findv1[0];
}
newS[position]=trueS[i]+10000*l;
}
}       
}
}
tempS=newS;
od=1;
len=0;
for (int i = 0; i < n; ++i) {
if(i==0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
Sod[i]=1;
od=od+1;
}
}
if(i>0){
if(tempS[i]==0){
Sod[i]=0;
}
if(tempS[i]>0){
tempa=tempS.subvec(0,i-1);
findv=find(tempa==tempS[i]);
len=findv.n_elem;
if(len==0){
Sod[i]=od;
od=od+1;
}
if(len>0){
tempv5=Sod.subvec(0,i-1);
findv1=tempv5.elem(find(tempa==tempS[i]));
Sod[i]=findv1[0];
}
}
}
}
S=Sod;
counts=S0v.subvec(0,od-1);  
for (int i = 0; i < n; ++i) {
counts[S[i]]=counts[S[i]]+1;
}
n2=n-counts[0];
k=od-1;
if(iter>Nb && iter%10==9){
hT=hT+h;
tempv7.fill(0);
tempv8.fill(0);
od=-1;
for (int q = 0;  q< Q; ++q){
if(countz[q]>0){
od=od+1;
tempv7[q]=countzus[od];
tempv8[q]=countzp[od];
}
}
for(int i = 0; i<nx; ++i){
for (int q = 0;  q< Q; ++q){
temp6=W*R::dnorm(x[i],logc,sqrt(Sz0),0);
if(tempv7[q]>0){ 
for (int l = 1;  l< tempv7[q]+1; ++l){
temp6=temp6+(1-W)*basecountcluster(l,q)/(tempv8[q]+M)*
R::dnorm(x[i],basemu(l,q),sqrt(alpha*sigmasq),0);
}
}
temp6=temp6+(1-W)*M/(tempv7[q]+M)*R::dnorm(x[i],mu0,sqrt(sigmasq),0);
predmat(i,q)=predmat(i,q)+temp6;      
}          
}
reptime=reptime+1;
}

if(Q>1){
if(k>0.5){ 
scount=counts.subvec(1,k); 
sz=S0v.subvec(0,k-1);
zcluster=sz;
zprobmat=Q0k.submat(0,0,Q-1,k-1);
tempv4=q0v;
for(int i = 0; i<n; ++i){
if(S[i]>0){
sz[S[i]-1]=Z[i];
zprob=conv_to< colvec >::from(wmat.submat(i,0,i,Q-1));
zprob=zprob+sm;
zprobmat.col(S[i]-1)=zprobmat.col(S[i]-1)+log(zprob);
tempv4[Z[i]-1]=tempv4[Z[i]-1]+1;
}
}
for(int j = 0; j<k; ++j){
tempv4[sz[j]-1]=tempv4[sz[j]-1]-scount[j];
for(int q = 0; q<Q; ++q){
zprobmat(q,j)=zprobmat(q,j)+lgamma(M+tempv4[q])-lgamma(M+tempv4[q]+scount[j]);
}
zprob=conv_to< colvec >::from(zprobmat.submat(0,j,Q-1,j));
zprob=zprob-max(zprob);
zprob=exp(zprob);
zprob=zprob/sum(zprob);
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
zcluster[j]=q+1;
break;}
}
tempv4[zcluster[j]-1]=tempv4[zcluster[j]-1]+scount[j];
}
} 
for(int i = 0; i<n; ++i){
if(S[i]==0){
zprob=conv_to< colvec >::from(wmat.submat(i,0,i,Q-1));
for (int q = 0;  q< Q; ++q) {
if (q==0) {
deflat=1.0;
}
if (q!=0) {
deflat=deflat-zprob[q-1];
}
probsuc=zprob[q]/deflat;
if (R::runif(0,1)<probsuc) {
Z[i]=q+1;
break;
}
}
}
if(S[i]>0){
Z[i]=zcluster[S[i]-1];
}
}



if(wdp<4){
logL1=0;
logL2=0;
logratio=0;
newpsi=psi*exp(R::rnorm(0,0.1));
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
probv1[q]=-psi*distance(t,q);
probv2[q]=-newpsi*distance(t,q);
}
probv1=probv1-max(probv1);
probv1=exp(probv1);
probv1=probv1/sum(probv1);
wmat.row(t)=conv_to< rowvec >::from(probv1);
probv2=probv2-max(probv2);
probv2=exp(probv2);
probv2=probv2/sum(probv2);
wmatn.row(t)=conv_to< rowvec >::from(probv2);
logL1=logL1+log(wmat(t,Z[t]-1));
logL2=logL2+log(wmatn(t,Z[t]-1));
}
logL1=logL1+R::dnorm(log(psi),1.03,0.02,1);
logL2=logL2+R::dnorm(log(newpsi),1.03,0.02,1); 
logratio=exp(logL2-logL1);
if(logratio>R::runif(0,1)){
psi=newpsi;
wmat=wmatn;
}
}
if(wdp>3){
sumpv=sumv;
logL1=0;
logL2=0;
logratio=0;
change=0;
for(int j = 0; j<Q; ++j){
if(wdp<6){
newpsi=R::qnorm(R::runif(R::pnorm(0,psiv[j],0.01,1,0),R::pnorm(4,psiv[j],0.01,1,0)),psiv[j],0.01,1,0);
}
if(wdp==6){
newpsi=R::qnorm(R::runif(R::pnorm(0,psiv[j],0.01,1,0),R::pnorm(2,psiv[j],0.01,1,0)),psiv[j],0.01,1,0);
}
logL1=0;
logL2=0;
logratio=0;
for(int t = 0; t<n; ++t){  
SS(t,j)=exp(-newpsi*distance(t,j));
sumv[t]=sumpv[t]+SS(t,j)-KS(t,j);
logL1=logL1+log(KS(t,Z[t]-1)/sumpv[t]);
if(j==Z[t]-1){
logL2=logL2+log(SS(t,Z[t]-1)/sumv[t]);
}
if(j!=Z[t]-1){
logL2=logL2+log(KS(t,Z[t]-1)/sumv[t]);
}
}
if(wdp<6){
logL1=logL1-log(R::pnorm(4,psiv[j],0.01,1,0)-R::pnorm(0,psiv[j],0.01,1,0));
logL2=logL2-log(R::pnorm(4,newpsi,0.01,1,0)-R::pnorm(0,newpsi,0.01,1,0));
}
if(wdp==6){
logL1=logL1-log(R::pnorm(2,psiv[j],0.01,1,0)-R::pnorm(0,psiv[j],0.01,1,0));
logL2=logL2-log(R::pnorm(2,newpsi,0.01,1,0)-R::pnorm(0,newpsi,0.01,1,0));
}
if(wdp==4){
logL1=logL1-psiv[j]+log(psiv[j]);
logL2=logL2-newpsi+log(newpsi);
}
if(wdp==5){
logL1=logL1-0.413*psiv[j];
logL2=logL2-0.413*newpsi;
}
if(wdp==6){
logL1=logL1-log(caup+psiv[j]*psiv[j]);
logL2=logL2-log(caup+newpsi*newpsi);
}
logratio=exp(logL2-logL1);
if(logratio>R::runif(0,1)){
psiv[j]=newpsi;
KS.col(j)=SS.col(j);
sumpv=sumv;
change=1;
}
if (change==0){
SS.col(j)=KS.col(j);
sumv=sumpv;
}
}
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
wmat(t,q)=exp(-psiv[q]*distance(t,q));
}
wmat.row(t)=wmat.row(t)/sum(wmat.row(t));
}
}
countz=q0v;
for(int i = 0; i<n; ++i){
countz[Z[i]-1]=countz[Z[i]-1]+1;
}
}


mustar=S0v.subvec(0,k);
mustar[0]=logc;
if(k>0.5){
tempmean=S0v.subvec(0,k-1);
tempvar=tempmean;
for(int i = 0; i < n; ++i) {
st=S[i];
if (st>0){
tempmean[st-1]=(y[i]-h[i])*(1-alpha)/(alpha+counts[st]*(1-alpha))+tempmean[st-1];
}
}
for(int j = 0; j < k; ++j) {
tempmean[j]=tempmean[j]+mu0*alpha/(alpha+counts[j+1]*(1-alpha));
tempvar[j]=sigmasq/(counts[j+1]/alpha+1/(1-alpha)); 
mustar[j+1]=R::rnorm(tempmean[j],sqrt(tempvar[j]));
}
}
for(int i = 0; i < n; ++i) {
mu[i]=mustar[S[i]];
if(iter%10==0){
mumat(draw,i)=mu[i];
}
}
para=0;
if(R::runif(0,1)<0.5){
para=1;
}
if(para==1){
hstar=h+mu0;
mustar1=mustar.subvec(1,k)-mu0;
var1=sigmaeta/((n-1)*(1-fi)*(1-fi)+(1-fi*fi));
mean1=var1*((1-fi*fi)/sigmaeta*hstar[0]+(1-fi)/sigmaeta*sum(hstar.subvec(1,n-1)-fi*hstar.subvec(0,n-2)));
temp6=mu0;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum(mustar1%mustar1)/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=hstar.subvec(1,n-1);
tempv3=tempv3-mu0;
tempv4=hstar.subvec(0,n-2);
tempv4=tempv4-mu0;
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi)*((hstar[i]-mu0)-(hstar[i-1]-mu0)*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
h=hstar-mu0;
mustar.subvec(1,k)=mustar1+mu0;
for(int i = 0; i<n; ++i){
if(S[i]>0){
mu[i]=mu[i]-temp6+mu0;
}
}
}
if(para==0){
mean1=(sum(mustar)-mustar[0])/k;
var1=(1-alpha)*sigmasq/k;
mu0=R::rnorm(mean1,sqrt(var1));
temp2=n2+k;
temp2=temp2/2+1;
temp3=3;
for(int i = 0; i<n; ++i){
if(S[i]>0){
temp3=temp3+(y[i]-h[i]-mu[i])*(y[i]-h[i]-mu[i])/2/alpha;
}
}
temp3=temp3+sum((mustar.subvec(1,k)-mu0)%(mustar.subvec(1,k)-mu0))/2/(1-alpha);
sigmasq=R::rgamma(temp2,1/temp3);
sigmasq=1/sigmasq;
logratio=0;
tempv3=h.subvec(1,n-1);
tempv4=h.subvec(0,n-2);
mean1=sum(tempv3%tempv4)/sum(tempv4%tempv4);
var1=sigmaeta/sum(tempv4%tempv4);
if(var1<0.0001){
var1=0.0001;
}
lb=R::pnorm(-0.995,mean1,sqrt(var1),1,0);
ub=R::pnorm(0.995,mean1,sqrt(var1),1,0);
newfi=R::qnorm(R::runif(lb,ub),mean1,sqrt(var1),1,0);
logratio=log(R::dnorm(fi,mean1,sqrt(var1),0)/(ub-lb))-log(R::dnorm(newfi,mean1,sqrt(var1),0)/(ub-lb));
lb=R::pnorm(-0.995,0.95,sqrt(0.5),1,0);
ub=R::pnorm(0.995,0.95,sqrt(0.5),1,0);
logratio=logratio+log(R::dnorm(newfi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio-(hstar[0]-mu0)*(hstar[0]-mu0)*(1-newfi*newfi)/2/sigmaeta+0.5*log(1-newfi*newfi);
logratio=logratio-log(R::dnorm(fi,0.95,sqrt(0.5),0)/(ub-lb));
logratio=logratio+(hstar[0]-mu0)*(hstar[0]-mu0)*(1-fi*fi)/2/sigmaeta-0.5*log(1-fi*fi);
if(log(R::runif(0,1))<logratio){
fi=newfi;
}
temp2=n;
temp2=temp2/2+2.5;
temp3=0.025;
for(int i = 1; i<n; ++i){
temp3=temp3+0.5*(h[i]-h[i-1]*fi)*(h[i]-h[i-1]*fi);
}
sigmaeta=R::rgamma(temp2,1/temp3);
sigmaeta=1/sigmaeta;
}
temp2=n-n2;
temp2=temp2+0.1;
temp3=n2;
temp3=temp3+0.9;
W=R::rbeta(temp2,temp3);
if(iter%10==0){
fiT[draw]=fi;
SnT[draw]=sigmaeta;
SzT[draw]=sigmasq;
MT[draw]=M;
WT[draw]=W;
Mu0T[draw]=mu0;
psiT[draw]=psi;
if(wdp>3 && iter>Nb){
psivt=psivt+psiv;
repsi=repsi+1;
}
draw=draw+1;
}
}

predmat=predmat/reptime;
hT=hT/reptime;
hT=exp(hT/2.0);
for(int q = 0; q<Q; ++q){
tempv=predmat.col(q);
temp=sum(tempv.subvec(1,999))+sum(tempv.subvec(0,998));
temp=temp/2.0;
temp=temp*50.0/999;
tempv=tempv/temp;
predmat.col(q)=tempv;
}
Nb=floor(Nb/10.0);
if(wdp<4){
psi=mean(psiT.subvec(Nb,Nm-1));
psivt.fill(psi);
}
if(wdp>3){
psivt=psivt/repsi;
}
for(int t = 0; t<n; ++t){
for(int q = 0; q<Q; ++q){
odwmat(t,q)=exp(-psiv[q]*distance(t,q));
}
odwmat.row(t)=odwmat.row(t)/sum(odwmat.row(t));
}


temp1=mean(SnT.subvec(Nb,Nm-1));
temp2=mean(fiT.subvec(Nb,Nm-1));
temp3=0;
vec pred0(1),pred1(1),predv(1),odwv(1),LPS(n);
LPS.fill(0);
double err=0;
double step=50.0/999;
for(int i = 0; i<10000; ++i){
if(err>1){
break;
}
for(int t = 0; t<n; ++t){
if(err>1){
break;
}
if(t==0){
h[0]=R::rnorm(0,sqrt(temp1/(1-temp2*temp2)));
}
if(t>0){
h[t]=temp2*h[t-1]+R::rnorm(0,sqrt(temp1));
}
temp4=y[t]-h[t];
if(temp4<x[0]){
err=1000;
}
if(temp4>x[nx-1]){
err=1000;
}
if(err==0){
loc=floor((temp4-x[0])/step);
pred0=conv_to< colvec >::from(predmat.submat(loc,0,loc,Q-1));
pred1=conv_to< colvec >::from(predmat.submat(loc+1,0,loc+1,Q-1));
predv=(temp4-x[loc])/step*(pred1-pred0)+pred0;
odwv=conv_to< colvec >::from(odwmat.submat(t,0,t,Q-1));
LPS[t]=LPS[t]+sum(predv%odwv);
}
}
}
if(err<1){
LPS=LPS/10000;
LPS=LPS.subvec(10-p,n-1);
temp3=mean(log(LPS));
}
if(err>1){
temp3=datum::nan;
}
sigmasq=mean(SzT.subvec(Nb,Nm-1));
W=mean(WT.subvec(Nb,Nm-1));
mu0=mean(Mu0T.subvec(Nb,Nm-1));
M=mean(MT.subvec(Nb,Nm-1));
return List::create(
Named("fi")=temp2,Named("Sn")=temp1,Named("mu0")=mu0,Named("Sz")=sigmasq,Named("M")=M,Named("W")=W,Named("psi")=psi,
Named("fiT")=fiT, Named("SnT")=SnT,Named("Mu0T")=Mu0T,Named("SzT")=SzT,Named("MT")=MT,Named("WT")=WT,Named("psiT")=psiT,
Named("LPS")=temp3,Named("hT")=hT,Named("mumat")=mumat);
}

// [[Rcpp::export]]
Rcpp::List cpp7(arma::vec Macf, arma::vec Szacf, arma::vec Snacf, arma::vec Wacf, 
arma::vec psiacf, arma::vec fiacf, arma::vec Mu0acf, double len, double tsize){
double temp=2/sqrt(tsize);
double a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0;
double b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0;
for(int i = 0; i<1001; ++i){
if(b1<0.5){
if(Snacf[i]<temp){
b1=1;
}
if(b1<0.5){
a1=a1+Snacf[i];
}
}
if(b2<0.5){
if(fiacf[i]<temp){
b2=1;
}
if(b2<0.5){
a2=a2+fiacf[i];
}
}
if(b3<0.5){
if(Mu0acf[i]<temp){
b3=1;
}
if(b3<0.5){
a3=a3+Mu0acf[i];
}
}
}
if(len>8.5){
for(int i = 0; i<1001; ++i){
if(b4<0.5){
if(Szacf[i]<temp){
b4=1;
}
if(b4<0.5){
a4=a4+Szacf[i];
}
}
if(b5<0.5){
if(Wacf[i]<temp){
b5=1;
}
if(b5<0.5){
a5=a5+Wacf[i];
}
}
if(b6<0.5){
if(Macf[i]<temp){
b6=1;
}
if(b6<0.5){
a6=a6+Macf[i];
}
}
if(b7<0.5){
if(psiacf[i]<temp){
b7=1;
}
if(b7<0.5){
a7=a7+psiacf[i];
}
}
}
}
a1=2*a1-1;
a2=2*a2-1;
a3=2*a3-1;
a4=2*a4-1;
a5=2*a5-1;
a6=2*a6-1;
a7=2*a7-1;
if(len>8.5){
a8=datum::nan;
}
if(len<8.5){
a6=datum::nan;
a4=datum::nan;
a7=datum::nan;
a5=datum::nan;
a8=a3;
a3=datum::nan;
}
return List::create(Named("M")=a6,Named("Sz")=a4,Named("fi")=a2,Named("Sn")=a1,
Named("W")=a5,Named("Mu0")=a3,Named("psi")=a7,Named("Mu")=a8);
}