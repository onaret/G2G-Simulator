# Function: get.g2stat
# Computes the G2 test statistic given two data matrices
# This is used internally by the G2 function
# Inputs: Z, W: both square, symmetric matrices with an equal number of rows
# Output: test statistic (single value)
# This is the most efficient way I found of computing this, making use of a trace property

get.g2stat.had <- function(W,Z)
{
 g2tstat <- sum( W * Z )
 g2tstat
}

#
# Function: G2
### Input 
### dep data and indep data with samples on the rows and genes on the columns
### grouping: Either a logical value = F or A matrix with a single column and same number of rows as samples. 
###         Column name should be defined.
###         Contains clinical information of the samples. 
###         Should have two groups only. 
### stand : scaling the columns of the data, logical value
### nperm : number of permutations 

### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$Sg


G2 <- function(dep.data,indep.data,stand = TRUE, nperm=1000,grouping=F){
      
    nperm = nperm
      ## check for the number of samples in dep and indep data

      if (nrow(dep.data)!=nrow(indep.data)){
         cat("number of samples not same in dep and indep data","\n")
      }

     

      ### check for centering and standardizing the data
       
      dep.data = scale(dep.data,center=T,scale=stand)
      indep.data = scale(indep.data,center=T,scale=stand)

        

     
#################################################################################

   ### If grouping is T divide the dataset into two groups and get G2 pval for each group.
   
   if(length(grouping) > 1 ){

    ## check the number of groups
  
      if (length(table(grouping))>2){
       cat("number of groups cannot exceed 2, please check grouping variable","\n")
      }
   ## Group 1
   dep = dep.data[grouping == names(table(grouping))[1],]
   indep = indep.data[grouping == names(table(grouping))[1],]
   
   ## Calculate Z = XX' and W = YY'
   Z = tcrossprod(indep)
   W = tcrossprod(dep)
   ### G2 for Group1
   samp_names = rownames(indep)
   Sg1 = get.g2stat.had(W,Z)
  

   ### Permutations

   perm_samp = matrix(0,nrow=nrow(indep),ncol=nperm)   ## generate the permutation matrix
       for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
        }

   ## permutation starts
   for (perm in 1:nperm){
        permX = Z[perm_samp[,perm],]
        permX = permX[,perm_samp[,perm]]

        Sg1 = c( Sg1, get.g2stat.had(W,permX) )
      }
 
      perm_samp.g1 = perm_samp  ## permutation matrix for output
##############################################################################

   ## Group 2
   dep = dep.data[grouping == names(table(grouping))[2],]
   indep = indep.data[grouping == names(table(grouping))[2],]

   ## Calculate Z = XX' and W = YY'
   Z = tcrossprod(indep)
   W = tcrossprod(dep)
   
   ### G2 for Group2
   samp_names = rownames(indep)
   Sg2 = get.g2stat.had(W,Z)
  

   ### Permutations

   perm_samp = matrix(0,nrow(indep),nperm)   ## generate the permutation matrix
       for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
        }

   ## permutation starts
   for (perm in 1:nperm){
        permX = Z[perm_samp[,perm],]
        permX = permX[,perm_samp[,perm]]

        Sg2 = c( Sg2, get.g2stat.had(W,permX) )
      }

       perm_samp.g2 = perm_samp   #### permutation matrix for output

########################################################################

    #### G2 test statistic
    Sg1 = t(as.matrix(Sg1))
    Sg2 = t(as.matrix(Sg2))
    colnames(Sg1) = rep("g1",length(Sg1))
    colnames(Sg2) = rep("g2",length(Sg2))

    ### calculate the p-values for g1 and g2

    G2p1 = mean(Sg1[1]<= c(Inf , Sg1[2:(nperm+1)]))
    G2p2 = mean(Sg2[1]<= c(Inf , Sg2[2:(nperm+1)]))
    
    ### prepare for output
    G2p = c(G2p1,G2p2)
    Sg = cbind(Sg1,Sg2)
    perm_samp = list(perm_samp.g1, perm_samp.g2)
    names(perm_samp) = c("g1.perm.mat","g2.perm.mat")
    names(G2p) = c(paste(colnames(grouping)," status",":",names(table(grouping))[1],sep=""),paste(colnames(grouping)," status",":",names(table(grouping))[2],sep=""))
  
   } else {
         #### No  grouping of the samples.
         
         samp_names = rownames(indep.data)
         
         ## Calculate Z = XX' and W = YY'
         Z = tcrossprod(indep.data)
         W = tcrossprod(dep.data)
         Sdi = (sum(diag(W)))*(sum(diag(Z)))  ### for standardization

   
         ### G2 
         samp_names = rownames(indep.data)
         Sg = get.g2stat.had(W,Z)
         Sd = Sg/Sdi
         ### Permutations
         perm_samp = matrix(0,nrow(indep.data),nperm)   ## generate the permutation matrix
         for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
         }
         ## permutation starts
         for (perm in 1:nperm){
            permX = Z[perm_samp[,perm],]
            permX = permX[,perm_samp[,perm]]
            Sg = c(Sg, get.g2stat.had(W,permX) )
            Sd = c(Sd,(Sg[perm+1]/Sdi))

         }


########################################################################

         #### G2 test statistic
         Sg = t(as.matrix(Sg))
         Sd = t(as.matrix(Sd))

         ### Calculte G2 pval
         G2p = mean(Sg[1]<= c(Inf , Sg[2:(nperm+1)]))
         names(G2p) = "G2 pval"
         
         G2pD = mean(Sd[1]<= c(Inf , Sd[2:(nperm+1)]))
         names(G2pD) = "G2 pval stand"
      }
     
return (list(perm = perm_samp,G2p = G2p,test_stat = Sg,std_test = Sd, std_pval = G2pD))
}
   
