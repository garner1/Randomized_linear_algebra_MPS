function [B,U,DB]=prepare_onesite(A,direction,factorization)

[D1,D2,d]=size(A); 
if strcmp(factorization,'svd'),
    switch direction
        case 'lr'
            A=permute(A,[3,1,2]); A=reshape(A,[d*D1,D2]); 
            [B,S,U]=svd2(A); DB=size(S,1);
            B=reshape(B,[d,D1,DB]); B=permute(B,[2,3,1]);
            U=S*U; 
        case 'rl'
            A=permute(A,[1,3,2]); A=reshape(A,[D1,d*D2]); 
            [U,S,B]=svd2(A); DB=size(S,1); 
            B=reshape(B,[DB,d,D2]); B=permute(B,[1,3,2]); 
            U=U*S;
    end
end

if strcmp(factorization,'qr'),
    switch direction
        case 'lr'
            A=permute(A,[3,1,2]); A=reshape(A,[d*D1,D2]); 
            [B,U]=qr2(A); DB=size(B,2);
            B=reshape(B,[d,D1,DB]); B=permute(B,[2,3,1]); 
        case 'rl'
            A=permute(A,[1,3,2]); A=reshape(A,[D1,d*D2]); 
            [B,U]=qr2(A'); DB=size(B',1);
            B=reshape(B',[DB,d,D2]); B=permute(B,[1,3,2]);
            U = U';
    end
end
