% ************************ one-site optimization **************************

function [A,E,Heff]=minimizeE_onesite(hsetj,Hleft,Hright,P)

DAl=size(Hleft{1},1); 
DAr=size(Hright{1},1); 
d=size(hsetj{1},1);

% calculation of Heff
M=size(hsetj,1);

Heff=0; 
for m=1:M
    Heffm=contracttensors(Hleft{m},3,2,Hright{m},3,2); 
    Heffm=contracttensors(Heffm,5,5,hsetj{m},3,3); 
    Heffm=permute(Heffm,[1,3,5,2,4,6]); 
    Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]); 
    Heff=Heff+Heffm;
end

% projection on orthogonal subspace 
if ~isempty(P), Heff=P'*Heff*P; end

%
% optimization
%
options.disp=0; 
[A,E]=eigs(Heff,1,'sr',options); 
if ~isempty(P), A=P*A; end 


A=reshape(A,[DAl,DAr,d]);
