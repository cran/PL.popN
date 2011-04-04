PL.popN0 <-
function(model,tau)
  {

## Inverse Logit Transformation
  logit<-function(xb)    
    {
    p<-1/(1+exp(-xb));
    return(p)
    }

## Population size estimation (with standard error) for a GLM.
 
  VarNhat.glm<-function(m,tau)   
    {
    P<-m$fitted.values;
    Pi<-(1-(1-P)^tau);
    Nhat<-sum(1/Pi);    # H-T estimator.
    var.beta<-vcov(m);
    X<-model.matrix(m);
    gdash.beta<-t(X)%*%(Pi^(-2)*(1-P)^tau*tau*P);
    varA<-sum((1-Pi)/Pi^2);
    varB<-(t(gdash.beta)%*%var.beta)%*%gdash.beta;
    varN<-as.vector(varA+varB);
    Se.Nhat<-sqrt(varN);
    list(Nhat=Nhat,Se.Nhat=Se.Nhat)
    }

## Population size estimation (with standard error) for a GAM.

  VarNhat.gam<-function(m,tau)   
    {
    P<-m$fitted.values ;
    Pi<-(1-(1-P)^tau);
    Nhat<-sum(1/Pi);    # H-T estimator.
    var.beta<-m$Ve;
    X<-model.matrix(m);
    gdash.beta<-t(X)%*%(Pi^(-2)*(1-P)^tau*tau*P);
    varA<-sum((1-Pi)/Pi^2);
    varB<-(t(gdash.beta)%*%var.beta)%*%gdash.beta;
    varN<-as.vector(varA+varB);
    Se.Nhat<-sqrt(varN);
    list(Nhat=Nhat,Se.Nhat=Se.Nhat)
    }

## Population size estimation (with standard error) for a GLMM.

  VarNhat.glmmPQL<-function(m,tau)      # Using glmmPQL()
    {
    P<-logit(fitted(m));
    Pi<-(1-(1-P)^tau);
    Nhat<-sum(1/Pi);         # H-T estimator.
    var.beta<-vcov(m);
    X<-model.matrix(m);
    gdash.beta<-t(X)%*%(Pi^(-2)*(1-P)^tau*tau*P);
    varA<-sum((1-Pi)/Pi^2);
    varB<-(t(gdash.beta)%*%var.beta)%*%gdash.beta;
    varN<-as.vector(varA+varB);
    Se.Nhat<-sqrt(varN);
    list(Nhat=Nhat,Se.Nhat=Se.Nhat)
    }

## Population size estimation (with standard error) for a SIMEX model.
  
  VarNhat.simex<-function(m,m.glm,tau)   
    {
    P<-m$fitted.values;
    Pi<-(1-(1-P)^tau);
    Nhat<-sum(1/Pi);      # H-T estimator.
    var.beta<-m$variance.asymptotic;
    X<-model.matrix(m.glm);
    gdash.beta<-t(X)%*%(Pi^(-2)*(1-P)^tau*tau*P);
    varA<-sum((1-Pi)/Pi^2);
    varB<-(t(gdash.beta)%*%var.beta)%*%gdash.beta;
    varN<-as.vector(varA+varB);
    Se.Nhat<-sqrt(varN);
    list(Nhat=Nhat,Se.Nhat=Se.Nhat)
    }

## Check for model type and evaluate population size estimate and standard error. 
  
  if(model$call[1]=="glm()")
    {
    D<-length(model$y);
    est<-VarNhat.glm(model,tau);
    N.hat<-est$Nhat;
    N.hat.se<-est$Se.Nhat;
    model.type<-"GLM";
    }
  
  if(model$call[1]=="gam()")
    {
    D<-length(model$y);
    est<-VarNhat.gam(model,tau);
    N.hat<-est$Nhat;
    N.hat.se<-est$Se.Nhat;
    model.type<-"GAM";
    }
    
  if(model$call[1]=="glmmPQL()")
    {
    D<-nrow(model$fitted);
    est<-VarNhat.glmmPQL(model,tau);
    N.hat<-est$Nhat;
    N.hat.se<-est$Se.Nhat;
    model.type<-"GLMM";
    }
  
  if(model$call[1]=="simex()")
    {
    D<-length(model$fitted.values);
    est<-VarNhat.simex(model,model$model,tau);
    N.hat<-est$Nhat
    N.hat.se<-est$Se.Nhat;
    model.type<-"SIMEX";
    }
      
  list(model.type=model.type,D=D,tau=tau,N.hat=N.hat,N.hat.se=N.hat.se)
  }

