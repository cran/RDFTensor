#15/4/2020

#Get scores from RESCAL decomposition (R & A) for pairs

#  Copyright (C) 2020 Abdelmoneim Amer Desouki, 
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


rescal_SO_Val<-function(R,A,Subj,P,Obj){
#Get scores from RESCAL decomposition (R & A) for pairs
#Subj: integer vector containing Subjects 
#Obj: integer vector containing Objects
# P: predicate index ( slice in tensor)
 
    if(length(Subj)!=length(Obj)){
        stop("rescal_SO_Val: List of subjects and objects must be of the same length")
    }
    
    val=rowSums((A[Subj,,drop=FALSE]%*%R[[P]])*(A[Obj,,drop=FALSE])) #Note * NOT %*%
    
    return(data.frame(Subj,P,Obj,val))
    
}

rescal_Trp_Val<-function(R,A,ntnsr,Sl_ix=NULL,verbose=1){
# Calculate score values corresponding to triples in tensor
# Sl_ix:index of slices (optional) default all

    ijkv=NULL
    # tmp=foreach(s =1:length(R),.packages='Matrix', .combine="rbind") %dopar%{
    if(is.null(Sl_ix)) Sl_ix=1:length(R)
    for(s in Sl_ix){
      if(verbose>1)print(sprintf('-------------=======s=%d==========---------',s))
      t1=proc.time()
      # t1n=A%*%R[[s]]%*%Matrix::t(A)
      ix=cbind(ntnsr$X[[s]]@i[ntnsr$X[[s]]@x==1]+1,ntnsr$X[[s]]@j[ntnsr$X[[s]]@x==1]+1) 
      #print(dim(ix))
      ResS=rescal_SO_Val(R, A, Subj=ix[,1], P=s , Obj=ix[,2] )
      
        ijkv=rbind(ijkv,ResS)
      
      t2=proc.time()
      if(verbose>0)print(sprintf('==== Slice %d #triples: %d == in %.2f secs',s,nrow(ix),(t2-t1)[3]))
    }
    
    return(ijkv)
}
###-------------------------------------------------------

