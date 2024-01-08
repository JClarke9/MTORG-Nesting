
# Creating a script for age nest analysis ---------------------------------

# This is the original Code that I modified slightly - it is from PhiDot
create.stage.var=function(data,agevar.name,stagevar.name,time.intervals,cutoff1, cutoff2)
{
  nocc=length(time.intervals)
  age.mat=matrix(data[,agevar.name],nrow=dim(data)[1],ncol=nocc-1)
  age.mat=t(t(age.mat)+cumsum(c(0,time.intervals[1:(nocc-2)])))
  stage.mat=t(apply(age.mat,1,function(x) as.numeric(x<=cutoff1)))
  stage.mat=data.frame(stage.mat)
  names(stage.mat)=paste(stagevar.name,1:(nocc-1),sep="")
  return(stage.mat)
}

# This is the new code - potentially use this for multiple stages
# create.stage.var=function(data,agevar.name,stagevar.name,time.intervals, cutoff1, cutoff2)
# {
#   nocc=length(time.intervals)
#   age.mat=matrix(data[,agevar.name],nrow=dim(data)[1],ncol=nocc-1)
#   age.mat=t(t(age.mat)+cumsum(c(0,time.intervals[1:(nocc-2)])))
#   stage.mat=apply(age.mat, c(1,2), function(row) {
#   ifelse(row < cutoff1, 1, ifelse(row >= cutoff2, 2, NA))
# })
#   stage.mat=data.frame(stage.mat)
#   names(stage.mat)=paste(stagevar.name,1:(nocc-1),sep="")
#   return(stage.mat)
# }
