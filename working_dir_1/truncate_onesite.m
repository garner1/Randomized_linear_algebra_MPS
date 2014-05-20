function [B,U,DB]=truncate_onesite(A,direction,maxD,factorization)
%used in truncate
[D1,D2,d]=size(A); 
if strcmp(factorization,'svd'),
    switch direction
        case 'lr'
            A=permute(A,[3,1,2]); A=reshape(A,[d*D1,D2]); 
            [B,S,U]=svd2(A); DB=size(S,1); D2=size(B,2);
            if DB>maxD,
                S=diag(S);   %transform the diag matrix S into a vector
                S=S(1:maxD); %truncate the singular values

                % THIS NORMALIZATION PRODUCES WRONG RESULTS; INSTEAD KEEP TRACK
                % OF U AS NORMALIZATION OF THE STATE;
    %             S=S./norm(S);  %rinormalize such that \sum_i s_i^2=1

                S=diag(S);   %transform the vector S into a diag matrix
                B=B(:,1:maxD); %truncate B
                U=U(1:maxD,:);  %truncate U
                B=reshape(B,[d,D1,maxD]);
            else B=reshape(B,[d,D1,D2]);
            end
            B=permute(B,[2,3,1]);
            U=S*U; 
        case 'rl'
            A=permute(A,[1,3,2]); A=reshape(A,[D1,d*D2]); 
            [U,S,B]=svd2(A); DB=size(S,1); D1=size(B,1);
            if DB>maxD,
                S=diag(S);
                S=S(1:maxD);
    %             S=S./norm(S);
                S=diag(S);
                B=B(1:maxD,:);
                U=U(:,1:maxD);
                B=reshape(B,[maxD,d,D2]);
            else B=reshape(B,[D1,d,D2]);
            end
            B=permute(B,[1,3,2]); 
            U=U*S;
    end
end

if strcmp(factorization,'fastsvd'),
    its = 0;
    switch direction
        case 'lr'
            A = permute(A,[3,1,2]); 
            A = reshape(A,[d*D1,D2]); 
            k = min([maxD,size(A)]);
            l = k;
            [B,S,U] = svd_fast(A,k,its,l); 
            B = reshape(B,[d,D1,k]);
            B = permute(B,[2,3,1]);
            U = S*U; 
        case 'rl'
            A = permute(A,[1,3,2]); 
            A = reshape(A,[D1,d*D2]); 
            k = min([maxD,size(A)]);
            l = k;
            [U,S,B] = svd_fast(A,k,its,l); 
            B = reshape(B,[k,d,D2]);
            B = permute(B,[1,3,2]); 
            U = U*S;
    end
end

if strcmp(factorization,'fastqr'),
    its = 0;
    switch direction
        case 'lr'
            A = permute(A,[3,1,2]); 
            A = reshape(A,[d*D1,D2]); 
            k = min([maxD,size(A)]);disp([maxD,size(A)])
            l = k;
%
%SUBSAMPLED RFFT
%
            [B,U] = qr_fast_srft(A,k,l); 
%
%RAND GAUSSIAN
%
%             [B,U] = qr_fast(A,k,its,l); 
            B = reshape(B,[d,D1,k]);
            B = permute(B,[2,3,1]);
        case 'rl'
            A = permute(A,[1,3,2]); 
            A = reshape(A,[D1,d*D2]); 
            k = min([maxD,size(A)]);disp([maxD,size(A)])
            l = k;
%
%SUBSAMPLED RFFT
%
            [B,U] = qr_fast_srft(A',k,l);
%
%RAND GAUSSIAN
%
%             [B,U] = qr_fast(A',k,its,l);
            B = B'; U = U';
            B = reshape(B,[k,d,D2]);
            B = permute(B,[1,3,2]);
    end
end
