ocpt.mean.initialise=function(data,penalty="Manual",pen.value=length(data),Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,shape=1,minseglen=1,alpha=1,verbose=FALSE){
    checkData(data)
    if(test.stat=="ECP"){
        ecpans = ocpt.np.initialise(Z=data, K=Q, delta=minseglen+1, alpha=alpha,verbose=verbose)
        return(ecpans)
    }
  if(minseglen<1){minseglen=1;warning('Minimum segment length for a change in mean is 1, automatically changed to be 1.')}
  if(length(data)<=2*minseglen){stop('Data length must be larger than 2*minseglen.')}
  #if(!((test.stat=="Normal")||(test.stat=="CUSUM"))){ stop("Invalid test statistic, must be Normal or CUSUM") #}
  if(penalty == "CROPS"){
    stop('Use cpt.mean from changepoint as CROPS is not available with changepoint.online')
  }
  if(penalty =="MBIC"){
      if(test.stat=="Normal"){
      cost_func="mean.norm.mbic"
  } else if(test.stat=="Exponential"){
      cost_func="meanvar.exp.mbic"
  } else if(test.stat=="Gamma"){
      cost_func="meanvar.gamma.mbic"
  }else if(test.stat=="Poisson"){
    cost_func="meanvar.poisson.mbic"
  }
  }else{
      if(test.stat=="Normal"){
          cost_func="mean.norm"
      } else if(test.stat=="Exponential"){
          cost_func="meanvar.exp"
      } else if(test.stat=="Gamma"){
          cost_func="meanvar.gamma"
      }else if(test.stat=="Poisson"){
          cost_func="meanvar.poisson"
      }else{
          stop("Not a valid test statistic. Must be Normal, Exponential, Gamma or Poisson")
      }
  }


          mu=mean(data)
          sumstat=cbind(c(0,cumsum(coredata(data))),c(0,cumsum(coredata(data)^2)),cumsum(c(0,(coredata(data)-mu)^2)))

          pen.value = online.decision(penalty, pen.value, length(data), diffparam=1, asymcheck=cost_func, method="AMOC")
          method = "PELT"
          ans=PELT.online.initialise(sumstat,pen=pen.value,cost_func = cost_func, shape = shape, minseglen = minseglen)

          return(online.class_input(sumstat=sumstat, cpttype="mean", method=method, test.stat=test.stat, penalty=penalty, pen.value=ans$penalty, minseglen=minseglen, param.estimates=param.estimates, out=sort(ans$cptsout[ans$cptsout>0]),shape=ans$shape,Q=Q,lastchangelike=ans$lastchangelike,lastchangecpts=ans$lastchangecpts,checklist=ans$checklist,nchecklist=ans$nchecklist,ndone=ans$ndone,nupdate=ans$nupdate,cost_func=ans$cost_func))
}

ocpt.mean.initialize=function(data,penalty="Manual",pen.value=length(data),Q=5,test.stat="Normal",class=TRUE,param.estimates=TRUE,shape=1,minseglen=1,alpha=1,verbose=FALSE){
return(ocpt.mean.initialise(data,penalty,pen.value,Q,test.stat,class,param.estimates,shape,minseglen,alpha,verbose))
}

ocpt.mean.update=function(previousanswer,newdata){

    checkData(newdata)
    if(class(previousanswer)=="ecp.ocpt"){
        ecpans = ocpt.np.update(previousanswer,newdata,K=2*(previousanswer@number))
        return(ecpans)
    }
    if(class(previousanswer@param.est) == "list"){
        param.estimates = TRUE
    }else{
        param.estimates = FALSE
    }

    nextans = PELT.online.update(previousanswer=previousanswer,newdata=newdata)

    return(online.class_input(sumstat=nextans$sumstat, cpttype="mean", method=previousanswer@method, test.stat=previousanswer@test.stat, penalty=previousanswer@pen.type, pen.value=nextans$penalty, minseglen=nextans$minseglen, param.estimates=param.estimates, out=sort(nextans$cptsout[nextans$cptsout>0]),shape = previousanswer@shape, lastchangelike=nextans$lastchangelike,lastchangecpts=nextans$lastchangecpts,checklist=nextans$checklist,nchecklist=nextans$nchecklist,ndone=nextans$ndone,nupdate=nextans$nupdate,cost_func=previousanswer@cost_func))

}
