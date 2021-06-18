
calcularbnotimizada<-function(md,itmax=200,cent=-1){
  n<-ncol(md)
  ass<-vector()
  ass_ant<-rep(2,n)
  Ass<-matrix(0,ncol=itmax,nrow=n)
  GuardaBn<-vector()

  if(cent==-1){
    cent<-sample(n,3)
  }

  for(i in 1:n){
    if(min(md[i,cent[1]],md[i,cent[2]],md[i,cent[3]])==md[i,cent[1]]){
      ass[i]<-1
    } else{
      if(md[i,cent[2]]<md[i,cent[3]]){
        ass[i]<-2
      }else{
        ass[i]<-3
      }
    }
  }
  it<-1
  Ass[,it]<-ass

  while(it<itmax && !prod(ass==ass_ant)){
    ass_ant<-ass
    ord<-sample(n,n)
    for(i in ord){
      ass[i]<-1
      bn1<--bn3(ass,md)
      ass[i]<-2
      bn2<--bn3(ass,md)
      ass[i]<-3
      bn33<--bn3(ass,md)
      if(min(bn1,bn2,bn33)==bn1){
        ass[i]<-1
      }else{
        if(bn2<bn33){
          ass[i]<-2
        }
      }
    }
    GuardaBn[it]<-bn3(ass,md)
    it<-it+1
    Ass[,it]<-ass
    it
  }

  ans<-list(which(ass==1),which(ass==2),which(ass==3),GuardaBn,it-1)
  names(ans)=c("grupo 1","grupo 2","grupo 3","Bn", "numIt")

  return(ans)
}

repeteBnmax<-function(mdist,repBn=50){
  Bn<-vector()
  iter<-vector()
  grupos<-matrix(0,ncol=dim(mdist)[1],nrow<-repBn)
  for(k in 1:repBn){
    resultado<-calcularbnotimizada(mdist)
    Bn[k]<-resultado[[4]][resultado[[5]]]
    grupos[k,resultado[[1]]]<-1;grupos[k,resultado[[2]]]<-2;grupos[k,resultado[[3]]]<-3
    iter[k]<-resultado[[5]]
  }
  Bn.max<-max(Bn)
  ind.max<-which(Bn==Bn.max)
  if(is.vector(ind.max)){ind.max<-ind.max[1]}
  grupos.max1<-which(grupos[ind.max,]==1);grupos.max2<-which(grupos[ind.max,]==2);grupos.max3<-which(grupos[ind.max,]==3)
  iter.max<-iter[ind.max]
  resposta<-list(grupos.max1,grupos.max2,grupos.max3,Bn.max,iter.max)
  names(resposta)=c("grupo 1","grupo 2","grupo 3","Bn", "numIt")
  return(resposta)
}

maxBn3size1<-function(md,itmax=200){
  n<-ncol(md)
  ass_ant<-rep(2,n)
  Ass<-matrix(0,ncol=itmax,nrow=n)
  GuardaBn<-vector()
  GuardaBn.max<-vector()
  ass<-rep(3,n)
  Ass.unit<-matrix(0,ncol=n,nrow=n)

  ass[1]<-1
  ass[sample(2:n,n%/%2)]<-2
  it.unit<-1

  while(it.unit<=n){

    it<-1
    while(it<itmax && !prod(ass==ass_ant)){
      ass_ant<-ass
      ord<-sample((1:n)[-it.unit],(n-1))
      for(i in ord){
        ass[i]<-2
        bn2<--bn3(ass,md)
        ass[i]<-3
        bn33<--bn3(ass,md)
        if(bn2<bn33){
          ass[i]<-2
        }else{
          ass[i]<-3
        }
      }

      GuardaBn[it]<--bn3(ass,md)
      it<-it+1
      Ass[,it]<-ass
    }

    GuardaBn.max[it.unit]<--bn3(ass,md)
    Ass.unit[,it.unit]<-ass
    it.unit<-it.unit+1
    ass<-rep(3,n)
    ass[it.unit]<-1
    ass[sample((1:n)[-it.unit],(n-1)%/%2)]<-2
    GuardaBn<-vector()
    Ass<-matrix(0,ncol=itmax,nrow=n)
  }

  maxBn<-min(GuardaBn.max)
  ind.max<-which(GuardaBn.max==maxBn)
  ass<-Ass.unit[,ind.max]

  ans<-list(-maxBn,which(ass==1),which(ass==2),which(ass==3))
  names(ans)=c("Bn","grupo 1", "grupo 2","grupo 3")

  return(ans)
}

repeteBnmaxsize1<-function(md,rep=15){
  Bn<-vector()
  iter<-vector()
  grupos<-matrix(0,ncol=dim(md)[1],nrow=rep)
  for(k in 1:rep){
    resultado<-maxBn3size1(md)
    Bn[k]<-resultado[[1]]
    grupos[k,resultado[[2]]]<-1;grupos[k,resultado[[3]]]<-2;grupos[k,resultado[[4]]]<-3

  }
  Bn.max<-max(Bn)
  ind.max<-which(Bn==Bn.max)
  if(is.vector(ind.max)){ind.max<-ind.max[1]}
  grupos.max1<-which(grupos[ind.max,]==1);grupos.max2<-which(grupos[ind.max,]==2);grupos.max3<-which(grupos[ind.max,]==3)
  iter.max<-iter[ind.max]
  resposta<-list(Bn.max,grupos.max1,grupos.max2,grupos.max3)
  names(resposta)=c("Maximized Bn","grupo 1","grupo 2","grupo 3")
  return(resposta)
}

calculaBnmaxrestrito<-function(md,n1_max,n1_min,n2,it_max=1000){
  if(n1_min>n1_max){
    print("ERROR: n1_max  must be larger than n1_min")
  }

  n<-dim(md)[1]
  it=1
  ass<-rep(3,n)
  ass_ant<-rep(2,n)
  Ass<-matrix(0,ncol=it_max,nrow=n)

  GuardaBnmax<-vector()
  n2_temp<-sample(n,n2)
  n1_temp<-sample((1:n)[-n2_temp],sample(n1_min:min(n-(n2+2),n1_max),1))
  ass[n1_temp]<-1;ass[n2_temp]<-2

  if(n1_min==n1_max){
    count=0
    while(it<it_max && count<max(n/2,10)){
      ass_ant<-ass

      ord<-sample(n,n)
      for(i in ord){
        v=(1:n)[-i]
        j=sample(v,1)
        f0= -bn3(ass,md)
        temp_ass=ass
        ass[j]=ass[i]
        ass[i]=temp_ass[j]

        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==n2) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==n2)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==n2) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==n2)){
          f1 = -bn3(ass,md)
        }else{
          f1=Inf
        }

        if(f0 < f1){
          ass=temp_ass
        }
      }

      GuardaBnmax[it]<--bn3(ass,md)
      it<-it+1
      Ass[,it]<-ass

      if(prod(ass==ass_ant)){
        count=count+1
      }
      else{
        count=0
      }

    }

  }else{

    while(it<it_max){
      ass_ant<-ass
      ord<-sample(n,n)
      for(i in ord){
        ass[i]<-1
        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==n2) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==n2)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==n2) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==n2)){
          f1<--bn3(ass,md)
        }else{
          f1<-Inf
        }

        ass[i]<-2
        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==n2) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==n2)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==n2) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==n2)){
          f2<--bn3(ass,md)
        }else{
          f2<-Inf
        }

        ass[i]<-3
        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==n2) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==n2) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==n2)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==n2) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==n2)){
          f3<--bn3(ass,md)
        }else{
          f3<-Inf
        }

        if(sum(is.infinite(c(f1,f2,f3)))==2){
          ind.aux<-which(is.finite(c(f1,f2,f3)))
          f0<-c(f1,f2,f3)[ind.aux]

          grupos<-(1:3)[-ind.aux]

          grupo.interesse<-which(ass==grupos[1])
          Bn.temp<-vector()
          cont.aux<-1
          ass[i]<-grupos[1]
          for(j in grupo.interesse){
            ass[j]<-ind.aux
            Bn.temp[cont.aux]<--bn3(ass,md)
            cont.aux<-cont.aux+1
            ass[j]<-grupos[1]
          }
          f4<-min(Bn.temp)
          ind.min.1<-grupo.interesse[which(Bn.temp==f4)]

          grupo.interesse2<-which(ass==grupos[2])
          Bn.temp.2<-vector()
          cont.aux<-1
          ass[i]<-grupos[2]
          for(j in grupo.interesse2){
            ass[j]<-ind.aux
            Bn.temp.2[cont.aux]<--bn3(ass,md)
            cont.aux<-cont.aux+1
            ass[j]<-grupos[2]
          }
          f5<-min(Bn.temp.2)
          ind.min.2<-grupo.interesse2[which(Bn.temp.2==f5)]

          if(min(f0,f4,f5)==f0){
            ass[i]<-ind.aux
          }else{
            if(f4<f5){
              ass[i]<-grupos[1]
              ass[ind.min.1]<-ind.aux
            }else{
              ass[i]<-grupos[2]
              ass[ind.min.2]<-ind.aux
            }
          }


        }else{
          if(f1==min(f1,f2,f3)){
            ass[i]<-1
          }else{
            if(f2==min(f2,f3)){
              ass[i]<-2
            }else{
              ass[i]<-3
            }
          }
        }
      }

      GuardaBnmax<--bn3(ass,md)
      Ass[,it]<-ass
      it<-it+1
    }

  }

  ans<-list(which(ass==1),which(ass==2),which(ass==3),-GuardaBnmax,it-1)
  names(ans)=c("grupo 1","grupo 2","grupo 3","Bn", "numIt")

  return(ans)
}

repeteBnmaxrestrito<-function(md,n1_max,n1_min,n2,rep=15){
  Bn<-vector()
  iter<-vector()
  grupos<-matrix(0,ncol=dim(md)[1],nrow=rep)
  for(k in 1:rep){
    resultado<-calculaBnmaxrestrito(md,n1_max,n1_min,n2)
    Bn[k]<-resultado[[4]]
    grupos[k,resultado[[1]]]<-1;grupos[k,resultado[[2]]]<-2;grupos[k,resultado[[3]]]<-3
    iter[k]<-resultado[[5]]
  }
  Bn.max<-max(Bn)
  ind.max<-which(Bn==Bn.max)
  if(is.vector(ind.max)){ind.max<-ind.max[1]}
  grupos.max1<-which(grupos[ind.max,]==1);grupos.max2<-which(grupos[ind.max,]==2);grupos.max3<-which(grupos[ind.max,]==3)
  iter.max<-iter[ind.max]
  resposta<-list(grupos.max1,grupos.max2,grupos.max3,Bn.max,iter.max)
  names(resposta)=c("grupo 1","grupo 2","grupo 3","Bn Padronizado", "numIt")
  return(resposta)
}

calculaBnmaxrestritogrupo1<-function(md,n1_max,n1_min,it_max=1000){

  if(n1_min>n1_max){
    print("ERROR: n1_max  must be larger than n1_min")
  }

  n<-dim(md)[1]
  it=1
  ass<-rep(3,n)
  ass_ant<-rep(2,n)
  Ass<-matrix(0,ncol=it_max,nrow=n)

  GuardaBnmax<-vector()

  n2_temp<-sample(n,1)
  n1_temp<-sample((1:n)[-n2_temp],sample(n1_min:min(n-3,n1_max),1))
  ass[n1_temp]<-1;ass[n2_temp]<-2

  if(n1_min==n1_max){
    count=0
    while(it<it_max && count<max(n/2,10)){
      ass_ant<-ass

      ord<-sample(n,n)
      for(i in ord){
        v=(1:n)[-i]
        j=sample(v,1)
        f0= -bn3(ass,md)
        temp_ass=ass
        ass[j]=ass[i]
        ass[i]=temp_ass[j]

        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==1) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==1)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==1) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==1)){
          f1 = -bn3(ass,md)
        }else{
          f1=Inf
        }

        if(f0 < f1){
          ass=temp_ass
        }
      }

      GuardaBnmax[it]<--bn3(ass,md)
      it<-it+1
      Ass[,it]<-ass

      if(prod(ass==ass_ant)){
        count=count+1
      }
      else{
        count=0
      }

    }

  }else{

    while(it<it_max){
      ass_ant<-ass
      ord<-sample(n,n)
      for(i in ord){
        ass[i]<-1

        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==1) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==1)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==1) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==1)){
          f1<--bn3(ass,md)
        }else{
          f1<-Inf
        }

        ass[i]<-2

        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==1) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==1)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==1) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==1)){
          f2<--bn3(ass,md)
        }else{
          f2<-Inf
        }

        ass[i]<-3

        if ((length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==2))==1) ||(length(which(ass==1))<=n1_max && length(which(ass==1))>=n1_min && length(which(ass==3))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==1))==1) ||(length(which(ass==2))<=n1_max && length(which(ass==2))>=n1_min && length(which(ass==3))==1)||(length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==1))==1) || (length(which(ass==3))<=n1_max && length(which(ass==3))>=n1_min && length(which(ass==2))==1)){
          f3<--bn3(ass,md)
        }else{
          f3<-Inf
        }

        if(sum(is.infinite(c(f1,f2,f3)))==2){

          ind.aux<-which(is.finite(c(f1,f2,f3)))
          f0<-c(f1,f2,f3)[ind.aux]

          grupos<-(1:3)[-ind.aux]

          grupo.interesse<-which(ass==grupos[1])
          Bn.temp<-vector()
          cont.aux<-1
          ass[i]<-grupos[1]
          for(j in grupo.interesse){
            ass[j]<-ind.aux
            Bn.temp[cont.aux]<--bn3(ass,md)
            cont.aux<-cont.aux+1
            ass[j]<-grupos[1]
          }
          f4<-min(Bn.temp)
          ind.min.1<-grupo.interesse[which(Bn.temp==f4)]

          grupo.interesse2<-which(ass==grupos[2])
          Bn.temp.2<-vector()
          cont.aux<-1
          ass[i]<-grupos[2]
          for(j in grupo.interesse2){
            ass[j]<-ind.aux
            Bn.temp.2[cont.aux]<--bn3(ass,md)
            cont.aux<-cont.aux+1
            ass[j]<-grupos[2]
          }
          f5<-min(Bn.temp.2)
          ind.min.2<-grupo.interesse2[which(Bn.temp.2==f5)]

          if(min(f0,f4,f5)==f0){
            ass[i]<-ind.aux
          }else{
            if(f4<f5){

              ass[i]<-grupos[1]
              ass[ind.min.1]<-ind.aux
            }else{

              ass[i]<-grupos[2]
              ass[ind.min.2]<-ind.aux
            }
          }


        }else{
          if(f1==min(f1,f2,f3)){
            ass[i]<-1
          }else{
            if(f2==min(f2,f3)){
              ass[i]<-2
            }else{
              ass[i]<-3
            }
          }
        }
      }

      GuardaBnmax<--bn3(ass,md)
      Ass[,it]<-ass
      it<-it+1
    }

  }

  ans<-list(which(ass==1),which(ass==2),which(ass==3),-GuardaBnmax,it-1)
  names(ans)=c("grupo 1","grupo 2","grupo 3","Bn", "numIt")

  return(ans)
}

repeteBnmaxrestritogrupo1<-function(md,n1_max,n1_min,rep=15){
  Bn<-vector()
  iter<-vector()
  grupos<-matrix(0,ncol=dim(md)[1],nrow=rep)
  for(k in 1:rep){
    resultado<-calculaBnmaxrestritogrupo1(md,n1_max,n1_min)
    Bn[k]<-resultado[[4]]
    grupos[k,resultado[[1]]]<-1;grupos[k,resultado[[2]]]<-2;grupos[k,resultado[[3]]]<-3
    iter[k]<-resultado[[5]]
  }
  Bn.max<-max(Bn)
  ind.max<-which(Bn==Bn.max)
  if(is.vector(ind.max)){ind.max<-ind.max[1]}
  grupos.max1<-which(grupos[ind.max,]==1);grupos.max2<-which(grupos[ind.max,]==2);grupos.max3<-which(grupos[ind.max,]==3)
  iter.max<-iter[ind.max]
  resposta<-list(grupos.max1,grupos.max2,grupos.max3,Bn.max,iter.max)
  names(resposta)=c("grupo 1","grupo 2","grupo 3","Bn Padronizado", "numIt")
  return(resposta)
}

funcaovar<-function(posicao,BootB,BootB1){
  n1<-length(which(posicao==1))
  n2<-length(which(posicao==2))
  n3<-length(which(posicao==3))

  n<-n1+n2+n3

  naux<-n%/%3

  somaeta<-function(ng){
    nng<-0
    for(i in 1:3){
      nng<-nng+((n-ng[i])/(2*(ng[i]-1)))
    }
    soma<-2*choose(n,2)*(1+(1/n)*(nng))
    return(soma)
  }

  eta2<-function(n){
    n.igual.1<-which(n==1)
    n2<-n[-n.igual.1][1];n3<-n[-n.igual.1][2]
    N<-1+n2+n3
    soma<-(2/((N^2)*(N-1)))*(2+((2*n2*n3)/(N-1))+((n2*(n3+2)^2)/((N-1)*(n2-1)))+((n3*(n2+2)^2)/((N-1)*(n3-1))))
    return(soma)
  }

  somaetacentral<-somaeta(c(naux,naux,n-2*naux))
  eta2central<-eta2(c(1,naux,n-(naux+1)))
  tauestimado<-BootB/(somaetacentral*(2/n*(n-1))^2)

  if(n1 >1 && n2>1 && n3>1){
    varBn<-(somaeta(c(n1,n2,n3))/somaetacentral)*BootB
  }else{

    if((n1 == 1 && n2>1 && n3>1)| (n2==1 && n1>1 && n3>1) | (n3==1 && n1>1 && n2>1)){
      varBn<- BootB1-(eta2central-eta2(c(n1,n2,n3)))*tauestimado
    } else{
      varBn<-1
    }
  }
  ans<-list(varBn,BootB,BootB1)
  names(ans)<-c("varbn","Bootb","Bootb1")
  return(ans)
}

calcularbnotimizadapad<-function(md,itmax=200,cent=-1){
  n<-ncol(md)

  ass<-vector()
  ass_ant<-rep(2,n)
  Ass<-matrix(0,ncol=itmax,nrow=n)
  GuardaBnPad<-vector()

  naux<-n%/%3

  posicaocentral<-c(rep(1,naux),rep(2,naux),rep(3,(n-2*naux)))
  Bncentral<-vector()
  for(fi in 1:1000){
    posicaofisher<-sample(posicaocentral)
    Bncentral[fi]<--bn3(posicaofisher,md)
  }
  BootB<-var(Bncentral)

  posicaocentralgrupo1<-c(1,rep(2,naux),rep(3,(n-naux-1)))
  Bncentralgrupo1<-vector()
  for(fi in 1:1000){
    posicaofisher<-sample(posicaocentralgrupo1)
    Bncentralgrupo1[fi]<--bn3(posicaofisher,md)
  }
  BootB1<-var(Bncentralgrupo1)

  if(cent==-1){
    cent<-sample(n,3)
  }

  for(i in 1:n){
    if(min(md[i,cent[1]],md[i,cent[2]],md[i,cent[3]])==md[i,cent[1]]){
      ass[i]<-1
    } else{
      if(md[i,cent[2]]<md[i,cent[3]]){
        ass[i]<-2
      }else{
        ass[i]<-3
      }
    }
  }

  it<-1
  Ass[,it]<-ass

  while(it<itmax && !prod(ass==ass_ant)){
    ass_ant<-ass

    ord<-sample(n,n)
    for(i in ord){
      ass[i]<-1
      bnpad1<-(-bn3(ass,md)/sqrt(funcaovar(ass,BootB,BootB1)$varbn))
      ass[i]<-2
      bnpad2<--bn3(ass,md)/sqrt(funcaovar(ass,BootB,BootB1)$varbn)
      ass[i]<-3
      bnpad3<--bn3(ass,md)/sqrt(funcaovar(ass,BootB,BootB1)$varbn)
      if(min(bnpad1,bnpad2,bnpad3)==bnpad1){
        ass[i]<-1
      }else{
        if(bnpad2<bnpad3){
          ass[i]<-2
        }
      }

    }
    GuardaBnPad[it]<--bn3(ass,md)/sqrt(funcaovar(ass,BootB,BootB1)$varbn)
    it<-it+1
    Ass[,it]<-ass

  }

  varBn<-funcaovar(ass,BootB,BootB1)$varbn

  ans<-list(which(ass==1),which(ass==2),which(ass==3),-GuardaBnPad,it-1,varBn,BootB,BootB1)
  names(ans)=c("grupo 1","grupo 2","grupo 3","Bn Padronizado", "numIt","VarBn","BootB","BootB1")

  return(ans)
}

repeteBnpadmax<-function(mdist,repBn=15){
  Bn<-vector()
  iter<-vector()
  varBn<-vector()
  BootB<-vector()
  BootB1<-vector()
  grupos<-matrix(0,ncol=dim(mdist)[1],nrow<-repBn)
  for(k in 1:repBn){
    resultado<-calcularbnotimizadapad(mdist)
    Bn[k]<-resultado[[4]][resultado[[5]]]
    grupos[k,resultado[[1]]]<-1;grupos[k,resultado[[2]]]<-2;grupos[k,resultado[[3]]]<-3
    iter[k]<-resultado[[5]]
    varBn[k]<-resultado[[6]]
    BootB[k]<-resultado[[7]]
    BootB1[k]<-resultado[[8]]
  }
  Bn.max<-max(Bn)
  ind.max<-which(Bn==Bn.max)
  if(is.vector(ind.max)){ind.max<-ind.max[1]}
  grupos.max1<-which(grupos[ind.max,]==1);grupos.max2<-which(grupos[ind.max,]==2);grupos.max3<-which(grupos[ind.max,]==3)
  iter.max<-iter[ind.max]
  varBn.max<-varBn[iter.max]
  BootB.max<-BootB[iter.max]
  BootB1.max<-BootB1[iter.max]
  resposta<-list(grupos.max1,grupos.max2,grupos.max3,Bn.max,iter.max,varBn.max,BootB.max,BootB1.max)
  names(resposta)=c("grupo 1","grupo 2","grupo 3","Bn Padronizado", "numIt","varBn","BootB","BootB1")
  return(resposta)
}

gama3<-function(n){
  se3<-(243*(3^(n-6))+(1+n+(n^2))-(2+n)*(2^(n-1)))/2
  delta3<-(2^(n-2)-n)*n
  return(se3+delta3)
}

