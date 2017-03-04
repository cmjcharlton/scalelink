# To run type source('Scalingmodelnew.R')
# Program for record linkage scaling algorithm. Version 1.0.  October 2016
#read files from current working directory
#The files are specified as follows:
#FOI	Name of File of interest (NFOI x PFOI). NFOI is no. of records
#LDF	Linking data file (NLDF x PLDF). NLDF is no. of records
#IDFOI	1 if first column of FOI is ID. 
#IDLDF	1 if first column of LDF is ID. 
#IDENTIFIERS_FOI	A row vector length PFOI with a 1 for an identifier, 0 else
#IDENTIFIERS_LDF	A row vector length PLDF with a 1 for an identifier, 0 else
#NCATS	Each identifier has k(j) categories specified in this row vector, coded 0 to k(j)-1. 
#MATCHALL	Set to 1 if records with complete disagreement or complete agreement are to be included in estimation: disabled this version
#NREPORT set the interval for printing the FOI record being accessed while computin A matrix. Default det to 100. 
#NOTE Current version 1.0 October 2016 only allows two (0,1) categories for each identifier
# Identifiers should be in same order in each of the two files.
# If any identifier is missing use the missing value code -9.999e+029    A random replacement for the whole record
#where this occurs will be used.

options(digits=10)

FOI<-read.table("FOI")
LDF<-read.table("LDF")
#MATCHALL<-read.table("MATCHALL")
#print(MATCHALL)
#<-read.table("FOI")
#LDF_ID<-read.table("LDF")
IDENTIFIERS_FOI<-read.table('IDENTIFIERS_FOI')
IDENTIFIERS_LDF<-read.table('IDENTIFIERS_LDF')

NEXISTS<-file.exists("NREPORT")
if(NEXISTS==TRUE)
{NREPORT<-read.table("NREPORT")}
if(NEXISTS==FALSE)
{NREPORT<-100}


MISSINGVALUE<-(-9.999e+029)
SIZEFOI<-dim(FOI)
SIZELDF<-dim (LDF)
NFOI<-SIZEFOI[1]
PFOI<-SIZEFOI[2]
NLDF<-SIZELDF[1]
PLDF<-SIZELDF[2]
IDENTPOSFOI<-matrix(0,0,0)
IDENTPOSLDF<-matrix(0,0,0)


for (I in 1:SIZEFOI[2])
{
    if(IDENTIFIERS_FOI[I]==1)        IDENTPOSFOI<-c(IDENTPOSFOI, I)
    if(IDENTIFIERS_LDF[I]==1)        IDENTPOSLDF<-c(IDENTPOSLDF, I)
}

NIDENT<-sum(IDENTIFIERS_FOI)
NEXISTS<-file.exists("NCATS")

if(NEXISTS==TRUE)
{
  NCATS<-read.table("NCATS")
  
}
#Assume binary identifiers if no NCATS file
if(NEXISTS==FALSE)
{
  
  
  NCATS<-matrix(1,1,NIDENT)*2
}



print(("Numbers of categories for identifiers"))
print(NCATS)

NCATSTOTAL<-sum(NCATS)

#End input

# randomly insert for missing values - first LDF then FOI
#create LDFNEW FOINEW if any missing values
#calc sum of each row - if <0 there is a missvalue

RECSUMLDF<-(rowSums(LDF))
LDFNEW<-matrix(0,0,SIZELDF[2])
for (NRECLDF in 1:SIZELDF[1])
{
  if( RECSUMLDF[NRECLDF]>-10000)
  {
    
    LDFNEW<-rbind(LDFNEW,LDF[NRECLDF,])
  }
}
SIZELDFNEW<-dim (LDFNEW)
NLDFNEW<-SIZELDFNEW[1]
for (NRECLDF in 1:SIZELDF[1])
{
  if( RECSUMLDF[NRECLDF]<(-10000))
  {
    
    
    N<-sample(NLDFNEW,size= 1, replace = TRUE, prob = NULL) 
    LDF[NRECLDF,]<-LDFNEW[N,]
    
  }
}


RECSUMFOI<-(rowSums(FOI))
FOINEW<-matrix(0,0,SIZEFOI[2])
for (NRECFOI in 1:SIZEFOI[1])
{
 if( RECSUMFOI[NRECFOI]>(-10000))
  {
    
    FOINEW<-rbind(FOINEW,FOI[NRECFOI,])
}
}
SIZEFOINEW<-dim (FOINEW)
NFOINEW<-SIZEFOINEW[1]
for (NRECFOI in 1:SIZEFOI[1])
{
  if( RECSUMFOI[NRECFOI]<(-10000))
  {
    
   
    N<-sample(NFOINEW,size= 1, replace = TRUE, prob = NULL) 
    FOI[NRECFOI,]<-FOINEW[N,]
    
  }
}


 #now re-use FOINEW LDFNEW just containg indicators
 LDFNEW<-matrix(0,NLDF,NIDENT)
 FOINEW<-matrix(0,NFOI,NIDENT)
 for (I in 1:NLDF){
   for(J in 1:NIDENT)
   {LDFNEW[I,J]<-LDF[I,IDENTPOSLDF[J]]}
 }
 for (I in 1:NFOI){
   for(J in 1:NIDENT)
 {FOINEW[I,J]<-FOI[I,IDENTPOSFOI[J]]}
 }
SIZEFOI<-SIZEFOI[1]
SIZELDF<-SIZELDF[1]

A<-matrix(0,NCATSTOTAL,NCATSTOTAL)
TOTMATCH<-0
TOTNOTMATCH<-0

#  Accumulates matrix values. Note: assumes all identifiers have just 2 categories.
FOIEXP<-matrix(0,1,NIDENT*2)
AGREE<-matrix(0,1,NIDENT)
for (NRECFOI in 1:NFOI){ 
for (NRECLDF in 1:NLDF){
  FOIEXP<-FOIEXP*0
 
   for (I in 1:NIDENT)
  {J<-{{FOINEW[NRECFOI,I]-LDFNEW[NRECLDF,I]}!=0}
  J<-I*2-J
  FOIEXP[1,J]<-1}
 
  A<-A+crossprod(FOIEXP,FOIEXP)
  }
    out<-NRECFOI%%NREPORT
  if (out==0)
  {
    print(paste( "FOI record number",NRECFOI))
  }
  }

for (I in 1:NCATSTOTAL)
{
A[I,I]<-A[I,I]*(NIDENT-1)/(NIDENT*NIDENT)
for (J in 1:NCATSTOTAL)
{
if(J!=I)
{
  A[J,I]<- (-A[J,I]/(NIDENT*NIDENT))
}

}
}



TOTAL<-NFOI*NLDF


#now A*
  q<-matrix(0,NCATSTOTAL,1)
  r<-matrix(0,NCATSTOTAL,1)
  # Note assumes 2 cats only for each indicator
  for (I in 1:NIDENT)
  {  q[(I-1)*2+1,1]<-1
    r[2*I,1]<-1
}

  q<-q/NIDENT
  r<-r/NIDENT
  B<-array(0,c(1,NCATSTOTAL))
  B<-cbind(B,0,1/NIDENT)

  B<-t(B)
q<-q*TOTAL
r<-r*TOTAL
B<-B*TOTAL
ASTAR<-cbind(2*A,-q,-r)
qa<-cbind(t(q),0,0)
ra<-cbind(t(r),0,0)
qr<-rbind(qa,ra)
ASTAR<-rbind(ASTAR,qr)


SCORES<-solve(ASTAR,B)
SCORES<-t(SCORES)

    L<-1
    for (I in 1:NIDENT)
    {
    J<-NCATS[1,I] 
    
    ADD<-SCORES[1,L]
    for (K in 1:J)
    {
    SCORES[1,L]<-SCORES[1,L]-ADD
    L=L+1
    }
    }
   # SCORES<-ENDPOINT_SCORES[1:NCATSTOTAL]
    SCORES<-100*SCORES
    options(digits=3)
    print("Endpoint scores")
    print(SCORES)
   
  
    #The folowing to be translated to save scoresdec 2016
    
  
    
  
    # parameterfile=fullfile(subfolder,'ASTAR');
     #csvwrite('parameterfile',ASTAR);
    # parameterfile=fullfile(subfolder,'SCORES');
    # csvwrite('parameterfile',SCORES);
    
    
    
    
    
    
    
   
    
    print(paste("Number_of_primary_file_records",NFOI))
   
    