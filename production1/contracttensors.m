function [Xout,numindX]=contracttensors(Xin,numindX,indX,Yin,numindY,indY)

Xsize=ones(1,numindX); Xsize(1:length(size(Xin)))=size(Xin); 
Ysize=ones(1,numindY); Ysize(1:length(size(Yin)))=size(Yin);

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
        Xin=permute(Xin,indX); 
        Xin=reshape(Xin,[1,prod(sizeX)]);
        
        Yin=permute(Yin,indY); 
        Yin=reshape(Yin,[prod(sizeY),1]);
        
        Xout=Xin*Yin; 
%         Xsize=1;
        
        return;
    
    else   %if one is NOT contracting all the indexes of X
        Xin=permute(Xin,[indXl,indX]); 
        Xin=reshape(Xin,[prod(sizeXl),prod(sizeX)]);
        
        Yin=permute(Yin,indY); 
        Yin=reshape(Yin,[prod(sizeY),1]);
        
        Xout=Xin*Yin; 
        Xsize=Xsize(indXl);
        
        Xout=reshape(Xout,[Xsize,1]);
        
        return 
    end
end
%%%%%%%%%

Xin1=permute(Xin,[indXl,indX]); 
Xin2=reshape(Xin1,[prod(sizeXl),prod(sizeX)]);

Yin1=permute(Yin,[indY,indYr]); 
Yin2=reshape(Yin1,[prod(sizeY),prod(sizeYr)]);

Xout=Xin2*Yin2; 
Xsize=[Xsize(indXl),Ysize(indYr)]; 
numindX=length(Xsize); 
Xout=reshape(Xout,[Xsize,1]);