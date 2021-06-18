#########################################################
###   Main Uclust function
#########################################################
#' U-statistic based significance clustering for three way partitions
#'
#' Partitions data into three groups only when these partitions are statistically significant.
#'  If no significant partition exists, the test will return "homogeneous".
#'
#' This is the significance clustering procedure of Bello et al. (2021).
#' The method first performs a homogeneity test to verify whether the data can be significantly
#' partitioned. If the hypothesis of homogeneity is rejected, then the method will search, among all
#' the significant partitions, for the partition that better separates the data, as measured by larger
#' \code{bn} statistic. This function should be used in high dimension small sample size settings.
#'
#'
#' Either \code{data} or \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean distance.
#'
#' Variance of \code{bn} is estimated through resampling, and thus, p-values may vary a bit in different runs.
#'
#' For more detail see
#' Bello, Debora Zava,  Marcio Valk and Gabriela Bettella Cybis.
#'  "Clustering inference in multiple groups." arXiv preprint arXiv:2106.09115 (2021).
#'See also \code{is_homo3}, \code{uclust}.
#'
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param alpha Significance level.
#' @param rep Number of times to repeat optimization procedures. Important for problems with
#' multiple optima.
#' @return  Returns a list with the following elements:\describe{
#'   \item{groups}{List with elements of final three groups}
#'   \item{p.value}{P-value for the test that renders the final partition, if heterogeneous.
#'   Homogeneity test p-value, if homogeneous.}
#'   \item{alpha_corrected}{Bonferroni corrected significance level for the test that renders the final
#'   partition, if heterogeneous. Homogeneity test significance level, if homogeneous.}
#'   \item{ishomo}{Logical, returns \code{TRUE} when the sample is homogeneous.}
#'   \item{Bn}{Value of Bn statistic for the final partition, if heterogeneous.
#'   Value of Bn statistic for the maximal homogeneity test partition, if homogeneous.}
#'   \item{varBn}{Variance estimate for final partition, if heterogeneous.
#'   Variance estimate for the maximal homogeneity test partition, if homogeneous.}
#' }
#'
#' @examples
#' set.seed(123)
#' x = matrix(rnorm(70000),nrow=7)  #creating homogeneous Gaussian dataset
#' res = uclust3(data=x)
#' res
#'
#' # uncomment to run
#' # x = matrix(rnorm(15000),nrow=15)
#' # x[1:6,] = x[1:6,]+1.5 #Heterogeneous dataset (first 5 samples have different mean)
#' # x[7:12,] = x[7:12,]+3
#' # res = uclust3(data=x)
#' # res$groups
#'
#'@export


uclust3<-function(md=NULL,data=NULL,alpha=0.05,rep=15){

  if (is.null(md)) {
    if (is.null(data)) {
      stop("No data provided")
    }
    md <- as.matrix(dist(data)^2)
  }
  if (sum("matrix" %in% class(md)) == 0) {
    stop("md is not of class matrix")
  }


  n<-dim(md)[1]

  is.h<-is_homo3(rep=rep,md) #  max std Bn
  ResultadoTesteIsHomo<-is.h   #remove?

  BootB<-is.h$BootB
  BootB1<-is.h$BootB1

  if(is.h$p.value< alpha){ # is_homo return homogeneous?
    oBn<-repeteBnmax(md,rep) # rep to find max Bn
    maxBn<-oBn[[4]]
    posicaoBnmax<-vector()
    posicaoBnmax[oBn[[1]]]<-1;posicaoBnmax[oBn[[2]]]<-2;posicaoBnmax[oBn[[3]]]<-3

    # for size one group
    n1<-length(oBn[[1]])

    # Bonferroni test for max Bn
    varbn<-funcaovar(posicaoBnmax,BootB,BootB1)$varbn
    pvalor.Bonf<-pnorm(maxBn/sqrt(varbn),lower.tail = FALSE)
    alpha.Bonf<-alpha/gama3(n)

    if(pvalor.Bonf<alpha.Bonf){
      Bn1.temp<-repeteBnmaxsize1(md)
      maxBn1.temp<-Bn1.temp[[1]]
      if(maxBn>maxBn1.temp){
        clust1<-oBn[[1]];clust2<-oBn[[2]];clust3<-oBn[[3]]
        p.valor<-pvalor.Bonf
      }else{
        clust1<-Bn1.temp[[2]];clust2<-Bn1.temp[[3]];clust3<-Bn1.temp[[4]]
        # Bonferroni test for max  Bn with group size 1
        posicaoBn1max.temp<-vector()
        posicaoBn1max.temp[clust1]<-1;posicaoBn1max.temp[clust2]<-2;posicaoBn1max.temp[clust3]<-3
        varbn1.temp<-funcaovar(posicaoBn1max.temp,BootB,BootB1)$varbn
        pvalor.Bonf<-pnorm(maxBn1.temp/sqrt(varbn1.temp),lower.tail = FALSE)
        p.valor<-pvalor.Bonf
      }
    }else{
      n2<-length(oBn[[2]]);n3<-length(oBn[[3]])
      minsize<-((n1==floor(n/3)&&n2==floor(n/3))||(n1==floor(n/3)&&n3==floor(n/3))||(n2==floor(n/3)&&n3==floor(n/3)))# os tamanhos com menor variancia para os quais nao queremos fazer otimizacao restrita

      # Bn for group size 1
      oBn1<-repeteBnmaxsize1(md) #
      maxBn1<-oBn1[[1]]

      # for centralized groups
      if(minsize==TRUE){
        clust1<-oBn1[[2]];clust2<-oBn1[[3]];clust3<-oBn1[[4]]
        # Bonferroni test for max Bn with size one group
        posicaoBn1max<-vector()
        posicaoBn1max[clust1]<-1;posicaoBn1max[clust2]<-2;posicaoBn1max[clust3]<-3
        varbn1<-funcaovar(posicaoBn1max,BootB,BootB1)$varbn
        pvalor.Bonf.<-pnorm(maxBn1/sqrt(varbn1),lower.tail = FALSE)
        alpha.Bonf<-alpha/gama3(n)
        n1<-1
      }

      # while  Bonferroni test significance is not TRUE and the groups are not central

      n1.1<-min(n1,n2,n3)
      while(minsize==FALSE &&pvalor.Bonf<alpha.Bonf){
       # repeat procedure for two groups in uhclust
        while(minsize==FALSE && pvalor.Bonf<alpha.Bonf){
          n2m<-min(n1,n2,n3)
          n2_min<-n2m+1
          n2_max<-(n-n1.1-n2_min)

          oBn<-repeteBnmaxrestrito(md,n2_max,n2_min,n1,rep=rep)

          maxBn<-oBn[[4]] # compare with Bn for size one group
          maxBnmax<-max(maxBn1,maxBn)

          # atualizing n2 value
          if (maxBnmax==maxBn1){ #
            n1<-1 # if group one's Bn is bigger
            pvalor.Bonf<-pnorm(maxBnmax/varbn1,lower.tail = TRUE)
            if(pvalor.Bonf>=alpha/gama3(n)){ # If Bonferroni's test is not significant, keep the interactions
              n2=n2_min
            }
          }else{ # If the larger  Bn was the restrict
            # Define n2 based on max Bn
            n1.temp<-length(oBn$grupo1);n2.temp<-length(oBn$grupo2);n3.temp<-length(oBn$grupo3)
            #
            if(n1.temp==n1.1){
              n2<-n2.temp
            }else{
              if(n2.temp==n1.1){
                n2<-n1.temp
              }
            }
            n3<-n-(n1.1+n2)
            pvalor.Bonf<-pnorm(maxBnmax/funcaovar(c(rep(1,n1),rep(2,n2),rep(3,n3)[[1]]),BootB,BootB1)$varbn,lower.tail = TRUE)
          }

          # update the 3 groups
          n3<-n-(n1.1+n2)

          # With new group sizes, we need to update the central configuration
          minsize<-((n1.1==floor(n/3)&&n2==floor(n/3))||(n1.1==floor(n/3)&&n3==floor(n/3))||(n2==floor(n/3)&&n3==floor(n/3)))# os tamanhos com menor variancia para os quais nao queremos fazer otimizacao restrita


          if (pvalor.Bonf<alpha.Bonf){
            clust1<-oBn[[2]];clust2<-oBn[[3]];clust3<-oBn[[4]]}
          if (minsize==TRUE && pvalor.Bonf<alpha.Bonf){
            clust1<-oBn1[[2]];clust2<-oBn1[[3]];clust3<-oBn1[[4]]
            posicaoBn1max[clust1]<-1;posicaoBn1max[clust2]<-2;posicaoBn1max[clust3]<-3
            varbn1<-funcaovar(posicaoBn1max,BootB,BootB1)$varbn
            pvalor.Bonf.<-pnorm(maxBn1/sqrt(varbn1),lower.tail = FALSE)
            n1=1
          }

        }

        n1.1<-n1.1+1

      }

    }

    p.valor<-pvalor.Bonf
    alpha_correct<-alpha.Bonf

  }else{ # If classic test returned homogeneous
    clust1<-1:n;clust2<-0;clust3<-0
    p.valor<-is.h$p.value
    n1<-n
    alpha_correct<-alpha
  }



  # answer

  if(n1!=n){ # If we have groups, then non-homogeneous elements
    if(pvalor.Bonf>alpha.Bonf){ # for the rare case of not being homogeneous but not finding significant groups
      n1<-length(is.h$grupo1)
      n2<-length(is.h$grupo2)
      n3<-length(is.h$grupo3)
      ishomo=FALSE
      p<-is.h$p.value
      alpha_correct<-alpha
      posicaoBn<-vector()
      posicaoBn[is.h$grupo1]<-1;posicaoBn[is.h$grupo2]<-2;posicaoBn[is.h$grupo3]<-3
      Bn<-bn3(posicaoBn,md)
      varbn<-funcaovar(posicaoBn,BootB,BootB1)$varbn
      clust1<-is.h$grupo1;clust2<-is.h$grupo2;clust3<-is.h$grupo3
      if(maxBn1>Bn){
        Bn=maxBn1
        clust1<-oBn1[[2]];clust2<-oBn1[[3]];clust3<-oBn1[[4]]
        posicaoBn1max.temp<-vector()
        posicaoBn1max.temp[clust1]<-1;posicaoBn1max.temp[clust2]<-2;posicaoBn1max.temp[clust3]<-3
        varbn<-funcaovar(posicaoBn1max.temp,BootB,BootB1)$varbn
        #fazer o teste
        p<-pnorm(Bn/sqrt(varbn),lower.tail = FALSE)
        if(p>alpha.Bonf){
          ishomo=TRUE
        }
      }
    }else{

      ishomo=FALSE
      #Calculate Bn with Group Configuration
      posicaoBn<-vector()
      posicaoBn[clust1]<-1;posicaoBn[clust2]<-2;posicaoBn[clust3]<-3
      Bn<-bn3(posicaoBn,md)
      varbn<-funcaovar(posicaoBn,BootB,BootB1)$varbn
    }
  }else{ # If we only have one group, then homogeneous elements
    ishomo=TRUE
    posicaoBn<-vector()
    clust1<-is.h[[2]];clust2<-is.h[[3]];clust3<-is.h[[4]]
    posicaoBn[clust1]<-1;posicaoBn[clust2]<-2;posicaoBn[clust3]<-3
    Bn<--bn3(posicaoBn,md)
    varbn<-funcaovar(posicaoBn,BootB,BootB1)$varbn
  }

  ans=list(list(cluster1=clust1,cluster2=clust2,cluster3=clust3),p.valor,alpha_correct,ishomo,Bn,varbn)#, ResultadoTesteIsHomo)

  names(ans)=c("groups","p.value","alpha_corrected","ishomo","Bn","varBn")#, "ResultadoTesteIsHomo")
  return(ans)
}
