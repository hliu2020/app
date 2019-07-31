

dfHclust = function(df) {
  # validate input
  stopifnot(inherits(df, "data.frame"))
  stopifnot(ncol(df)>1)
  # obtain software
  require(shiny)
  require(cluster)
}
# global variables ... 
cmeths = c("ward.D", "ward.D2",
           "single", "complete", "average", "mcquitty",
           "median", "centroid")
dmeths = c("euclidean", "maximum", "manhattan", "canberra",
           "binary")