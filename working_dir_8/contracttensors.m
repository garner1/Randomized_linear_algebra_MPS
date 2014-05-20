function [X,numindX]=contracttensors(X,numindX,indX,Y,numindY,indY)

Xsize=ones(1,numindX); Xsize(1:length(size(X)))=size(X); 
Ysize=ones(1,numindY); Ysize(1:length(size(Y)))=size(Y);

indXl=1:numindX; indXl(indX)=[]; %the meaning is clear from the use in line 34
indYr=1:numindY; indYr(indY)=[];

sizeXl=Xsize(indXl); 
sizeX=Xsize(indX); 
sizeYr=Ysize(indYr); 
sizeY=Ysize(indY);

if prod(sizeX)~=prod(sizeY)
    error('indX and indY are not of same dimension.');
end

%special cases
if isempty(indYr)  %if one is contracting all the indexes of Y
    if isempty(indXl)  %if one is contracting all the indexes of X
        X=permute(X,[indX]); 
        X=reshape(X,[1,prod(sizeX)]);
        
       
        Y=permute(Y,[indY]); 
        Y=reshape(Y,[prod(sizeY),1]);
        
        X=X*Y; 
        Xsize=1;
        
        return;
    
    else   %if one is NOT contracting all the indexes of X
        X=permute(X,[indXl,indX]); 
        X=reshape(X,[prod(sizeXl),prod(sizeX)]);
        
        Y=permute(Y,[indY]); 
        Y=reshape(Y,[prod(sizeY),1]);
        
        X=X*Y; 
        Xsize=Xsize(indXl);
        
        X=reshape(X,[Xsize,1]);
        
        return 
    end
end
%%%%%%%%%

X=permute(X,[indXl,indX]); 
X=reshape(X,[prod(sizeXl),prod(sizeX)]);

Y=permute(Y,[indY,indYr]); 
Y=reshape(Y,[prod(sizeY),prod(sizeYr)]);

X=X*Y; 
Xsize=[Xsize(indXl),Ysize(indYr)]; 
numindX=length(Xsize); 
X=reshape(X,[Xsize,1]);