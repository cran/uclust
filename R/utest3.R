
#' U-test for three groups
#'
#'Test for the separation of three groups.
#' The null hypothesis states that the groups are homogeneous and the alternative
#' hypothesis states that at least one is separated from the others.
#'
#'
#'
#' Either \code{data} or \code{md} should be provided.
#' If data are entered directly, Bn will be computed considering the squared Euclidean
#'  distance.
#'
#' For more detail see
#' Bello, Debora Zava,  Marcio Valk and Gabriela Bettella Cybis.
#'  "Clustering inference in multiple groups." arXiv preprint arXiv:2106.09115 (2021).
#'
#' @param group_id A vector of 1s, 2s and 3s indicating to which group the samples belong.
#' Must be in the same order as data or md.
#' @param md Matrix of distances between all data points.
#' @param data Data matrix. Each row represents an observation.
#' @param alpha Significance level
#' @param numB Number of resampling iterations.
#'
#' @return Returns a list with the following elements:\describe{
#'   \item{is.homo}{Logical of whether test indicates that data is homogeneous}
#'   \item{Pvalue}{Replication based p-value}
#'   \item{Bn}{Test Statistic}
#'   \item{sdBn}{Standard error for Bn statistic computed through resampling}
#'
#' }
#'
#' @seealso \code{\link{bn3}},\code{\link{utest}},\code{\link{is_homo3}}
#'
#' @examples
#'
#'# Simulate a dataset with two separate groups,
#'# the first row has mean -4, the next 5 rows have mean 0 and the last 5 rows have mean 4.
#' data <- matrix(c(rnorm(15, -4),rnorm(75, 0), rnorm(75, 4)), nrow = 11, byrow=TRUE)
#' # U test for mixed up groups
#' utest3(group_id=c(1,2,3,1,2,3,1,2,3,1,2), data=data, numB=3000)
#' # U test for correct group definitions
#' utest3(group_id=c(1,2,2,2,2,2,3,3,3,3,3), data=data, numB=3000)
#'
#'
#'@export
## Public version of U test for 3 groups






utest3<-function(group_id,md=NULL,data=NULL,alpha=0.05,numB=1000){
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

  Bn<-bn3(group_id,md)
  Bnboot<-vector()
  for(b in 1:numB){
    group_idboot<-sample(group_id)
    Bnboot[b]<-bn3(group_idboot,md)
  }
  sdBn<-sd(Bnboot)

  pvalor<-2*pnorm(abs(Bn/sdBn),lower.tail=FALSE)

  if(pvalor<alpha){
    ans<-list(FALSE,pvalor,Bn,sdBn)
    names(ans)<-c("is.homo","pvalue","Bn","sdBn")
  }else{
    ans<-list(TRUE,pvalor,Bn,sdBn)
    names(ans)<-c("is.homo","pvalue","Bn","sdBn")
  }
   # message(paste("\t U-test for group separation  \n\nTest Statistic Bn =",
   # round(ans$Bn, digits = 4), "\t p-value = ", round(ans$Pvalue,
   # digits = 4), "\nAlternative hypothesis: The groups are not homogeneous,
   # \nthere exists some separation between groups. "))
 return(ans)
}
