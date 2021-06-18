#############################################################
#### Homogeneity test main function
#############################################################
#' U-statistic based homogeneity test for 3 groups
#'
#' Homogeneity test based on the statistic bn3. The test assesses whether there exists a data
#' partition for which three group separation is statistically significant according to utest3.
#' The null hypothesis is overall sample homogeneity, and a sample is considered homogeneous if
#'  it cannot be divided into three groups with at least one significantly different from the others.
#'
#' This is the homogeneity test of Bello et al. (2021).
#' The test is performed through two steps: an optimization procedure that finds the data partition that
#' maximizes the standardized Bn and a test for the resulting maximal partition.
#' Should be used in high dimension small sample size settings.
#'
#'
#'
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
#'
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param rep Number of times to repeat optimization procedure. Important for problems with
#' multiple optima.
#' @param test_max Logical indicating whether to employ the max test
#' @param alpha Significance level
#' @return Returns a list with the following elements:\describe{
#'   \item{stdBn}{Test statistic. Maximum standardized Bn.}
#'   \item{group1}{Elements in group 1 in the maximal partition. (obs: this is not the best
#'   partition for the data, see \code{uclust3})}
#'   \item{group2}{Elements in group 2 in the maximal partition.}
#'   \item{group3}{Elements in group 3 in the maximal partition.}
#'   \item{pvalue.Bonferroni}{P-value for the homogeneity test.}
#'   \item{alpha_Bonferroni}{Alpha after Bonferroni correction}
#'   \item{bootB}{Resampling variance estimate for partitions with central group sizes.}
#'   \item{bootB1}{Resampling variance estimate for partitions with one group of size 1.}
#'   \item{varBn}{Estimated variance of Bn for maximal standardized Bn configuration.}
#'
#'
#' }
#'
#'
#' @examples
#' set.seed(123)
#' x = matrix(rnorm(70000),nrow=7)  #creating homogeneous Gaussian dataset
#' res = is_homo3(data=x)
#' res
#'
#' #uncomment to run
#' # x = matrix(rnorm(18000),nrow=18)
#' # x[1:5,] = x[1:5,]+0.5 #Heterogeneous dataset (first 5 samples have different mean)
#' # x[6:9,] = x[6:9,]+1.5
#' # res = is_homo3(data=x)
#' #  res
#' # md = as.matrix(dist(x)^2) #squared Euclidean distances for the same data
#' # res = is_homo3(md)       # uncomment to run

#'
#' # Multidimensional sacling plot of distance matrix
#' #fit <- cmdscale(md, eig = TRUE, k = 2)
#' #x <- fit$points[, 1]
#' #y <- fit$points[, 2]
#' #plot(x,y, main=paste("Homogeneity test: p-value =",res$p.MaxTest))
#'
#'@export


is_homo3<-function(md = NULL, data = NULL,rep=20,test_max=TRUE,alpha=0.05){
  if (is.null(md)) {
    if (is.null(data)) {
      stop("No data provided")
    }
    md <- as.matrix(dist(data)^2)
  }
  if (sum("matrix" %in% class(md)) == 0) {
    stop("md is not of class matrix")
  }
  n <- dim(md)[1]
  if (n <= 4) {
    stop("samples size n is too small for homogeneity test")
  }


  nt<-gama3(n)

  resultado<-repeteBnpadmax(md,rep)
  Bn<-resultado[[4]]
  BootB<-resultado[[7]]
  BootB1<-resultado[[8]]

  if(test_max==TRUE){

    if(nt<=2^(28)){
      pvalor<-1-exp(nt*pnorm(Bn,log.p=TRUE))
    }else{
      bn=sqrt(2*log(nt))-(log(log(nt))+log(4*pi*log(2)^2))/(2*sqrt(2*log(nt)))
      an=(log((4*log(2)^2)/log(4/3)^2))/(2*sqrt(2*log(nt)))
      bgumbel=(1/an)*(Bn-bn)
      pvalor<-1-exp(-exp(-bgumbel))
    }
    g1<-resultado[[1]];g2<-resultado[[2]];g3<-resultado[[3]]
    varBn<-resultado[[6]]
    BootB<-resultado[[7]]
    BbootB1<-resultado[[8]]

    ans<-list(Bn,g1,g2,g3,pvalor,BootB,BootB1,varBn)
    names(ans)=list("stdBn","grupo1","grupo2","grupo3","p.value","BootB","BootB1","varBn")

  }else{

    pvalor.Bonf<-pnorm(Bn,lower.tail=FALSE)
    alpha.Bonf<-alpha/nt
    ans=list(Bn,g1,g2,g3,pvalor.Bonf,alpha.Bonf,BootB,BootB1,varBn)
    names(ans)=list("stdBn","grupo1","grupo2","grupo3","pvalue.Bonferroni","alpha_Bonferroni","BootB","BootB1","varBn")

  }
  return(ans)
}




