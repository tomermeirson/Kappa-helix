#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# This is the code to assign k-helix 
# Author: Tomer Meirson, 2022
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


library(bio3d)
library(dplyr)
library(magrittr)
library(stringr)

# Setting the path to STRIDE is required to run the function
# STRIDE (stride_WIN32.exe) can be downloaded from the following URL:
# http://webclu.bio.wzw.tum.de/stride/

assign.ss.khelix.stride = function(pdb,chain.id,K.config=NULL,path.stride,path.pdb=NULL){
  ## input:
  # pdb - Protein Data Bank ID
  # chain.id - chain to be assigned
  # K.config - configurations of thresholds for k-helix assignment. Default can be accessed through default.K.config()
  # path.stride - path to STRIDE software. Required
  # path.pdb - optional: path for download PDB files. Default is working directory
  
  ## output: data.frame of residues with assigned secondary structure based on STRIDE along with updated k-helix assignment
  
  if(is.null(K.config)) K.config = default.K.config() # if empty use the default settings
  if(is.null(path.pdb)) path.pdb=getwd() # if no path for PDB use the current wordking directory
  
  # if PDB does not exist download 
  if(!file.exists(paste0(path.pdb,pdb,'.pdb'))) get.pdb(pdb,path = path.pdb)
  
  # load PDB file
  raw = system(paste0('"',path.stride,'/stride_WIN32.exe" ','-r', chain.id,paste0(' "',path.pdb,'/',pdb,'.pdb','"')),intern = T) # use only the relevant chain.id
  
  if(length(raw)>1){ # if stride fails - probably status 5 
    str = grep('^ASG',raw,value = T) %>% as.matrix()
    
    str.list = lapply(str, function(x) str_split(x,'\\s+',simplify = T))
    
    
    data = do.call(rbind,str.list) %>% as.data.frame()
    
    colnames(data) = c('row','res.id','chain','no','seq','ss','old_type','phi','psi','area','pdb')
    
    data %<>% mutate(phi= as.double(as.character(phi)))
    data %<>% mutate(psi= as.double(as.character(psi)))
    data %<>% mutate(chain= as.character(chain))
    data %<>% mutate(no= as.character(no))

    # at terminal residues or  breaks, only diherdral angle is unambiguous (not 360), thus use only the relevant angle
    data %<>% mutate(dif = ifelse(phi==360.0,sqrt((psi-K.config$psi)^2),
                                  ifelse(psi==360.0,sqrt((phi-K.config$phi)^2),
                                         (sqrt((phi-K.config$phi)^2+(psi-K.config$psi)^2))/2)))
    
    data %<>% mutate(kappa = 0) # initialize 
    data %<>% mutate(old_type = as.character(old_type))
    data %<>% mutate(new_type = as.character(old_type))
    
    # filter-in only the relevant chain
    data %<>% filter(chain == chain.id)
    
    K.res = K.config$minK.res - 1 # reduce by 1 because of indexing
    K.iterations = with(K.config,maxK.res-minK.res) # how many averaging windows of increasing size to test  
    
    # iterate over increasing averaging windows
    for(iter in 0:(K.iterations)){
      K.res = K.res+iter
      # Assign kappa helix
      for(i in 1:(nrow(data)-K.res)){
        # print(i)
        if(mean(data$dif[i:(i+K.res)]) < K.config$K.thresh){
          data$kappa[i:(i+K.res)] = 1
        }
      }
    }
    
    data$no %<>% as.integer()

    data$new_type[data$kappa==1] = 'Kappa'
    
    return(data)
    
  }
}


default.K.config = function(K.thresh = 17,K.phi = -78,K.psi = 146,minK.res = 2,maxK.res=4){
  
  ## input:
  # Based on parameters used in doi.org/10.1093/bioinformatics/btz527
  # K.phi - reference phi dihedral angle. Default is to -78. In DSSP-PPII the default is -75
  # K.psi - reference psi dihedral angle. Default is to 146. In DSSP-PPII the default is 145
  # K.thresh - kappa assignment threshold. Default is 17
  # minK.res - minimum successive residues required to assign k-helix. Default is 2
  # maK.res - maximum size of averaging window of residues  to assign k-helix. Default is 4
  
  ## output: default K-configuration

  
  # in Bioinformatics I used phi/psi (-78/+146)
  K.phi = -78 # kappa phi  reference value
  K.psi = 146 # kappa psi reference value

  K.config = data.frame(phi=K.phi,psi=K.psi,K.thresh=K.thresh,
                        minK.res=minK.res,maxK.res)
  return(K.config)
}

