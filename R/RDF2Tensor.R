## 31/1/2018
## represent RDF data as sparse matrices representing frontal slices 
## to be a class

#  Copyright (C) 2018 Abdelmoneim Amer Desouki, 
#   Data Science Group, Paderborn University, Germany.
#  All right reserved.
#  Email: desouki@mail.upb.de
#
#  This file is part of RDFTensor.
#
#  RDFTensor is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RDFTensor is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RDFTensor.  If not, see <http://www.gnu.org/licenses/>.
	
getTensor<-function(trp,SO=NULL,P=NULL){
# trp: dataframe of strings with 3 columns S,P,O
	if(is.null(SO)) SO=unique(c(trp[,1],trp[,3]))
	if(is.null(P))  P=unique(trp[,2])
	modes=c(length(SO),length(SO),length(P))
	# print(sprintf("Expected memory size:%f Gb",4*prod(modes)/1024/1024/1024))#GB
	print("The modes:")
	print(modes)
	#dogfood 1.8 TB
	# s,o,p
	i=match(trp[,1],SO)
	j=match(trp[,3],SO)
	k=match(trp[,2],P)
	## unique(i,j,k)
	G =cbind(i,j,k)
	# RDFTensor <-list(n,m,X)#T_nxnxm
	n=length(SO)
	m=length(P)
	X <- list()
	for(p in 1:m){
		Xp=Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(n, n))
		Xp[G[G[,3]==p,1:2,drop=FALSE]]=1
		X[[p]]=methods::as(Xp, "TsparseMatrix")#i,j,x
	}
	return(list(X=X,n=n,m=m,SO=SO,P=P))
}

getTensor3m <- function(trp, S = NULL, P = NULL, O = NULL){
#subject and object values may have different indices
	if(is.null(S)) S=unique(trp[,1])
	if(is.null(O)) O=unique(trp[,3])
	if(is.null(P))  P=unique(trp[,2])
	modes=c(length(S),length(O),length(P))
	# print(sprintf("Expected memory size:%f Gb",4*prod(modes)/1024/1024/1024))#GB
	print("The modes:")
	print(modes)
	#dogfood 1.8 TB
	# s,o,p
	i=match(trp[,1],S)
	j=match(trp[,3],O)
	k=match(trp[,2],P)
	## unique(i,j,k)
	G =cbind(i,j,k)
	# RDFTensor <-list(n,m,X)#T_nxnxm
	n=length(S)
	m=length(O)
	l=length(P)
	X <- list()
	for(p in 1:l){
		Xp=Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(n, m))
		Xp[G[G[,3]==p,1:2,drop=FALSE]]=1
		X[[p]]=methods::as(Xp, "TsparseMatrix")#i,j,x
	}
	return(list(X=X,n=n,m=m,l=l,S=S,O=O,P=P))
}

# 29/4/2018
getTnsrijk <- function(X){
    ijk=NULL
	for(k in 1:length(X)){#NB:
		if(any(X[[k]]@x==1)){#numeric tensor
		   tt<- methods::as(X[[k]], "TsparseMatrix")
		   ijk=rbind(ijk,cbind(tt@i[X[[k]]@x==1]+1,k,tt@j[X[[k]]@x==1]+1))
		}
	}
	return(ijk)
}
###
getmode1slices <- function(X){
	## returns X as mode1 slices(Horizontal slices)
	n=nrow(X[[1]])
	m=ncol(X[[1]])
	l=length(X)
	X1=list(n)
	G=getTnsrijk(X)
	for(s in 1:n){
		Xp=Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(l, m))
		Xp[G[G[,1]==s,2:3,drop=FALSE]]=1
		X1[[s]]=methods::as(Xp, "TsparseMatrix")#i,j,x
	}
	return(X1)
}
###
tnsr2trp<-function(tnsr){
  ijk=getTnsrijk(tnsr$X)
  trp = cbind(unlist(tnsr$SO)[ijk[,1]],unlist(tnsr$P)[ijk[,2]],unlist(tnsr$SO)[ijk[,3]])
  return(trp)
}
# Xh=getmode1slices(tnsr$X)
# 14/3/2018
# from RDFTensor to xyz format to be read by Python module.

#7/5/2018
truncTnsr <- function(X,H=0,L=0,F=0){
	#Horiz slices
	#Lateral slices
	#Frontal slices
	ijk=NULL
	for(k in 1:length(X)){#NB:
		if(any(X[[k]]@x==1)){#numeric tensor
		   tt<- methods::as(X[[k]], "TsparseMatrix")
		   ijk=rbind(ijk,cbind(tt@i[X[[k]]@x==1]+1,k,tt@j[X[[k]]@x==1]+1))
		}
	}
	
	print('Calculating fibers nonzeros...')
	trp = cbind(ijk[,1],ijk[,2],ijk[,3])
	pf1=aggregate(trp[,2],by=list(trp[,1]),FUN="length")
	sf1=aggregate(trp[,1],by=list(trp[,2]),FUN="length")
	of1=aggregate(trp[,3],by=list(trp[,1],trp[,2]),FUN="length")
    
	if(F>0){
	
	}
}
# remove slices with nonzeros less than tnnz

tnsr2xyz<-function(tnsr,fname){
	# count from 1
	con <- file(fname,"w")
	for(k in 1:tnsr$l){
		print(sprintf("predicate:%d",k))
		# ij=which(tnsr$X[[i]]==1,arr.ind=TRUE)
		s1=methods::as(tnsr$X[[k]],"TsparseMatrix")
		i=s1@i+1
		j=s1@j+1
		if(tnsr$X[[k]]@x[1]==0){# caused by bad initialization
			i=i[-1]
			j=j[-1]
		}
		writeLines(as.character(paste(i,k,j)),con)
	}
	close(con)
}	

#12/1/2020
tensor2mat<-function(X,binary=FALSE,symmetrize=FALSE){
##X list of matrices (lateral slices)
    # ijk=getTnsrijk(X)
    alli=NULL
    allj=NULL
    for(k in 1:length(X)){#NB:
		if(any(X[[k]]@x==1)){#numeric tensor
		   # tt<- methods::as(X[[k]], "TsparseMatrix")
		   # ijk=rbind(ijk,cbind(tt@i[X[[k]]@x==1]+1,k,tt@j[X[[k]]@x==1]+1))
           alli=append(alli,X[[k]]@i[X[[k]]@x==1]+1)
           allj=append(allj,X[[k]]@j[X[[k]]@x==1]+1)
		}
	}
    n=nrow(X[[1]])
    
    if(binary && !symmetrize){
        Xp=Matrix::sparseMatrix(i=alli,j=allj,x=1,dims=c(n, n))
        Xp@x[Xp@x > 0]=1 
    }else{ 
        if(symmetrize){
           ii=append(alli,allj); jj=append(allj,alli)
           if(binary){
              # tmpij=cbind(ii,jj)[!duplicated(paste(ii,jj)),]
              Xp=Matrix::sparseMatrix(i=ii,j=jj,x=1,dims=c(n, n))
               Xp@x[Xp@x > 0]=1               
            }else{
               # tmpc=table(paste(ii,jj))
               # tmpij=matrix(as.integer(unlist(strsplit(names(tmpc),split=' '))),ncol=2,byrow=TRUE)
               Xp=Matrix::sparseMatrix(i=ii,j=jj,x=1,dims=c(n, n)) 
           }
        }else{#not symmetrize & not binary
               Xp=Matrix::sparseMatrix(i=alli,j=allj,x=1,dims=c(n, n))            
        }
    }
return(Xp)
}

getPrdcnts<-function(X){
  ijk=getTnsrijk(X)
  tt=table(ijk[,2])
  return(as.vector(tt[as.character(sort(as.integer(names(tt))))]))
}

# Usage
 # fname='sider_xyz.txt'
 # print(load('sider_tensor_noltr.RData'))
# tnsr2xyz(tnsr,fname)

 # fname='dogfood_xyz.txt'
 # print(load('D:\\RDF\\mats\\t2i_dogfood.RData'))
 # tnsr=getTensor(T2I)
 # tnsr2xyz(tnsr,fname)

# ------------------------------------------
#OhneLtr
 # print(load('D:\\RDF\\mats\\t2i_dogfood.RData'))
 # zz=table(unlist(T2I[substring(T2I[,3],1,1)!='\"',2]))
 # tnsr=getTensor(T2I[substring(T2I[,3],1,1)!='\"',])
 # tnsr2xyz(tnsr,'dogfoodOhneLtr_xyz.txt')
	# name='dogfood'
	# org_nt=cbind(T2I[substring(T2I[,3],1,1)!='\"',],dot='.')
	# write.table(file=sprintf('org_%s_OhneLtr.nt',name),org_nt,sep=' ',row.names=FALSE)
# ------------------------------------------
#	'Tensor Times Matrix (m-Mode Product)
#'Contracted (m-Mode) product between a Tensor of arbitrary number of modes and a matrix. The result is folded back into Tensor.
