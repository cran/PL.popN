PL.popN <-
function(model,tau)
  {
   est1<-PL.popN0(model,tau);
   class(est1)<-"PL.popN"
   est1;
  }

