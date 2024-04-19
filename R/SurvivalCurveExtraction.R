# Copyright (c) 2023 Mathyn Vervaart
# Licensed under the MIT License

#####################################################################################
# Load packages 
#####################################################################################
library(compiler)
enableJIT(3)

#####################################################################################
# Function for extracting curves from postscript file, taken from: 
# Liu Z, Rich B, Hanley JA. Recovering the raw data behind a non-parametric survival curve. Systematic Reviews. 2014
#####################################################################################
ExtractCurvesFromLargeFileFunction <- function(filename,n.lines.max,n.lines.at.once) {
  
  parse =function(Line){
    part = unlist(strsplit( Line , " " )) ; 
    x.val = as.numeric( part[1] ) ; 
    y.val = as.numeric( part[2] )
    #mode = 1*(part[3]=="mo") +0*(part[3]=="cv")
    mode = 1*(part[3]=="m") +0*(part[3]=="l")
    return(c(x.val,y.val,mode))
  }
  
  TFile = file(filename, "r")
  
  N.lines=0
  X.v=rep(NA,n.lines.max) ; Y.v=rep(NA,n.lines.max) ; M.v=rep(NA,n.lines.max) ;
  eof = 0;
  msg="";
  block=1
  n.start=1
  last.line.so.far = n.lines.at.once
  n.blocks = ceiling(n.lines.max / 10000)  
  
  while (eof == 0 & last.line.so.far <= n.lines.max) { 
    
    lines = try( readLines(TFile, n=n.lines.at.once) )  ;
    lines[substr(lines,1,11)=="% Copyright"] = "blank" # avoid trouble
    
    print(lines[1:2])
    
    if( is.na( lines[1] ) ) { eof = 1; msg = "lines do not exist" } 
    
    if(block<10) print(noquote( paste("block",
                                      toString(block),
                                      toString( (block-1)*n.lines.at.once+1 ),
                                      toString( block*n.lines.at.once ),
                                      msg)
    ) )
    if(!eof) {
      n.lines= length(lines)
      N.lines=N.lines+n.lines
      L= grep('[0-9]+ l$', lines) ; n.L=length(L) ; 
      MO=grep('[0-9]+ m$', lines) ; n.MO=length(MO)  

      EL = function(i) ( (i+1) %in% L ) | (i == MO[n.MO] )
      WH = sapply(MO, EL )
      elig.MO = MO[ WH ]
      
      all.lines = sort( c( elig.MO, L) )
      #print("post all.lines")
      kept.lines = lines[all.lines]
      vector=unlist( lapply(kept.lines,parse) ) ;
      x.vals = vector[ seq(1,length(vector),3) ];
      y.vals = vector[ seq(2,length(vector),3) ];
      m.vals = vector[ seq(3,length(vector),3) ];
      n.stop = n.start + length(x.vals) -1
      X.v[n.start:n.stop] = x.vals
      Y.v[n.start:n.stop] = y.vals
      M.v[n.start:n.stop] = m.vals
      n.start=n.stop+1
      block=block+1
      last.line.so.far = last.line.so.far + n.lines.at.once
    }  
  } # end of while ...
  
  close(TFile) 
  
  X.v=X.v[!is.na(X.v)] ; Y.v=Y.v[!is.na(Y.v)] ; M.v=M.v[!is.na(M.v)]
  
  Y.max=max(Y.v); 
  
  # Y.v = Y.max-Y.v # flip vertically
  
  Y.min=min(Y.v);  Y.max=max(Y.v) # revised
  X.min=min(X.v);  X.max=max(X.v)
  
  l.n = cumsum(M.v) # so can split X.v and Y.v on this vector
  
  X = split(X.v,l.n) ; Y = split(Y.v,l.n) ;
  
  LENGTH = as.vector(table(l.n))
  
  X = X[LENGTH>1] ; Y = Y[LENGTH>1]
  LENGTH = LENGTH[LENGTH>1]
  
  return(list(X,Y,LENGTH,X.min,Y.min,X.max,Y.max))
  
}


#####################################################################################
# Function for extracting lines from postscript file, taken from: 
# Liu Z, Rich B, Hanley JA. Recovering the raw data behind a non-parametric survival curve. Systematic Reviews. 2014
#####################################################################################
ExtractLinesFromLargeFileFunction = function(filename,n.lines.max,n.lines.at.once) {
  
  parse =function(Line){
    part = unlist(strsplit( Line , " " )) ; 
    x.val = as.numeric( part[1] ) ; 
    y.val = as.numeric( part[2] )
    mode = 1*(part[3]=="mo") +0*(part[3]=="li")
    return(c(x.val,y.val,mode))
  }
  
  TFile = file(filename, "r")
  
  N.lines=0
  X.v=rep(NA,n.lines.max) ; Y.v=rep(NA,n.lines.max) ; M.v=rep(NA,n.lines.max) ;
  eof = 0;
  msg="";
  block=1
  n.start=1
  last.line.so.far = n.lines.at.once
  n.blocks = ceiling(n.lines.max / 10000)  
  
  while (eof == 0 & last.line.so.far <= n.lines.max) { 
    
    lines = try( readLines(TFile, n=n.lines.at.once) )  ;
    #print(lines[1:2])
    
    if( is.na( lines[1] ) ) { eof = 1; msg = "lines do not exist" } 
    
    if(block<10) print(noquote( paste("block",
                                      toString(block),
                                      toString( (block-1)*n.lines.at.once+1 ),
                                      toString( block*n.lines.at.once ),
                                      msg)
    ) )
    if(!eof) {
      n.lines= length(lines)
      N.lines=N.lines+n.lines
      L= grep('[0-9]+ li$', lines) ; n.L=length(L) ; 
      MO=grep('[0-9]+ mo$', lines) ; n.MO=length(MO)  
      elig.MO = MO[sapply(MO,function(i) ( (i+1) %in% L ) | (i == MO[n.MO] ) )]
      all.lines = sort( c( elig.MO, L) )
      kept.lines = lines[all.lines]
      vector=unlist( lapply(kept.lines,parse) ) ;
      x.vals = vector[ seq(1,length(vector),3) ];
      y.vals = vector[ seq(2,length(vector),3) ];
      m.vals = vector[ seq(3,length(vector),3) ];
      n.stop = n.start + length(x.vals) -1
      X.v[n.start:n.stop] = x.vals
      Y.v[n.start:n.stop] = y.vals
      M.v[n.start:n.stop] = m.vals
      n.start=n.stop+1
      block=block+1
      last.line.so.far = last.line.so.far + n.lines.at.once
    }  
  } # end of while ...
  
  close(TFile) 
  
  X.v=X.v[!is.na(X.v)] ; Y.v=Y.v[!is.na(Y.v)] ; M.v=M.v[!is.na(M.v)]
  
  Y.max=max(Y.v); 
  
  Y.v = Y.max-Y.v # flip vertically
  
  Y.min=min(Y.v);  Y.max=max(Y.v) # revised
  X.min=min(X.v);  X.max=max(X.v)
  
  l.n = cumsum(M.v) # so can split X.v and Y.v on this vector
  
  X = split(X.v,l.n) ; Y = split(Y.v,l.n) ;
  
  LENGTH = as.vector(table(l.n))
  
  X = X[LENGTH>1] ; Y = Y[LENGTH>1]
  LENGTH = LENGTH[LENGTH>1]
  
  return(list(X,Y,LENGTH,X.min,Y.min,X.max,Y.max))
  
}