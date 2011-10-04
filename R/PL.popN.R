PL.popN <-
function(model,tau,D)
  {
   est1<-PL.popN0(model,tau,D);
   class(est1)<-"PL.popN"
   est1;
  }

