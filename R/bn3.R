#' Computes Bn Statistic for 3 Groups.
#'
#'Returns the value for the Bn statistic that measures the degree of separation between three groups.
#'The statistic is computed as a combination of differences of average within group and between group distances.
#'Large values of Bn indicate large group separation. Under overall sample homogeneity we have E(Bn)=0.
#'
#' Either \code{data} OR \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean distance.
#'
#' For more detail see
#' Bello, Debora Zava,  Marcio Valk and Gabriela Bettella Cybis.
#'  "Clustering inference in multiple groups." arXiv preprint arXiv:2106.09115 (2021).
#'
#' @param group_id A vector of 1s, 2s and 3s indicating to which group the samples belong. Must be in the same order as data or md.
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @return Value of the Bn3 statistic.
#'
#' @examples
#' n=7
#' set.seed(1234)
#' x=matrix(rnorm(n*10),ncol=10)
#' bn3(c(1,2,2,2,3,3,3),data=x)     # option (a) entering the data matrix directly
#' md=as.matrix(dist(x))^2
#' bn3(c(1,2,2,2,3,3,3),md)         # option (b) entering the distance matrix
#' @export


bn3<-function(group_id,md=NULL,data=NULL){
  if (is.null(md)) {
    if (is.null(data)) {
      stop("No data provided")
    }
    md <- as.matrix(dist(data)^2)
  }
  if (sum("matrix" %in% class(md)) == 0) {
    stop("md is not of class matrix")
  }

  grupo<-list()
  grupo[[1]]<-which(group_id==1) ;grupo[[2]]<-which(group_id==2) ; grupo[[3]]<-which(group_id==3)
  ng<-NULL
  ng[1]<-length(grupo[[1]]);ng[2]<-length(grupo[[2]]);ng[3]<-length(grupo[[3]])
  n<-sum(ng)
  dim_d <- dim(md)[1]
  if(dim_d != n){
    stop("Incorrect dimension or group_id")
  }

  if(ng[1]==0 | ng[2]==0 |ng[3]==0){
    return(-Inf)
  }

  if(ng[1]>1 && ng[2]>1 && ng[3]>1){
    Ug<-matrix(0,ncol=3,nrow=3)
    for(k in 1:3){
      i<-1
      while(i<ng[k]){
        j<-i+1
        while(j<=ng[k]){
          Ug[k,k]<-Ug[k,k]+md[grupo[[k]][i],grupo[[k]][j]]
          j<-j+1
        }
        i<-i+1
      }
    }
    for(k in 1:3){
      Ug[k,k]<-Ug[k,k]*(choose(ng[k],2)^(-1))
    }
    for(k in 1:2){
      for(l in (k+1):3){
        i<-1
        while(i<=ng[k]){
          j<-1
          while(j<=ng[l]){
            Ug[k,l]<-Ug[k,l]+md[grupo[[k]][i],grupo[[l]][j]]
            j<-j+1
          }
          i<-i+1
        }
      }
    }

    for(k in 1:2){
      for(l in (k+1):3){
        Ug[k,l]<-Ug[k,l]*(ng[k]*ng[l])^(-1)
      }
    }

    Bn<-(ng[1]*ng[2])/(n*(n-1))*(2*Ug[1,2]-Ug[1,1]-Ug[2,2])+(ng[1]*ng[3])/(n*(n-1))*(2*Ug[1,3]-Ug[1,1]-Ug[3,3])+(ng[2]*ng[3]/(n*(n-1)))*(2*Ug[2,3]-Ug[2,2]-Ug[3,3])

    return(Bn)
  }

  if(length(which(ng==1))==1){
    indicealone<-which(ng==1)

    Ug<-matrix(0,ncol=3,nrow=3)
    for(k in 1:3){
      if(k==indicealone){
        Ug[k,k]<-0
      }else{
        i<-1
        while(i<ng[k]){
          j<-i+1
          while(j<=ng[k]){
            Ug[k,k]<-Ug[k,k]+md[grupo[[k]][i],grupo[[k]][j]]
            j<-j+1
          }
          i<-i+1
        }
      }
    }
    for(k in 1:3){
      if(k != indicealone){
        Ug[k,k]<-Ug[k,k]*(choose(ng[k],2))^(-1)
      }
    }
    for(k in 1:2){
      for(l in (k+1):3){
        i<-1
        while(i<=ng[k]){
          j<-1
          while(j<=ng[l]){
            Ug[k,l]=Ug[l,k]<-Ug[k,l]+md[grupo[[k]][i],grupo[[l]][j]]
            j<-j+1
          }
          i<-i+1
        }
      }
    }
    for(k in 1:2){
      for(l in (k+1):3){
        Ug[k,l]=Ug[l,k]<-Ug[k,l]*(ng[k]*ng[l])^(-1)
      }
    }
    Bn<-0
    k<-1
    while(k<3){
      l<-k+1
      while(l<=3){
        if(k==indicealone){
          Bn<-Bn+(2*ng[l]/(n*(n-1)))*(Ug[k,l]-Ug[l,l])
        }else{
          if(l==indicealone){
            Bn<-Bn+(2*ng[k]/(n*(n-1)))*(Ug[k,l]-Ug[k,k])
          }else{
            Bn<-Bn+((ng[k]*ng[l])/(n*(n-1)))*(2*Ug[k,l]-Ug[k,k]-Ug[l,l])
          }
        }
        l<-l+1
      }
      k<-k+1
    }
    return(Bn)

  }
  if(length(which(ng==1))>1){
    return(-Inf)

  }

}
