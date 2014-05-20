function [Q,listn,j] = adaptive_rand(mps,mpo,site,r,epsilon,iter)
%
% algo 4.2 in REF
%
[D,~,~] = size(mps{site});
x = size(mpo{site},1);
omega = cell(1,r);
for ind = 1:r
   omega{ind} = randn(D,x,D)+1i*randn(D,x,D); 
end
y = cell(1,r);
listn = [];
for ind = 1:r
    y{ind} = mpoxrand(mpo,mps,site,'right',omega{ind},'normal');
    for ii=1:iter
        y{ind} = mpoxrand(mpo,mps,site,'right',y{ind},'dagger');
        y{ind} = mpoxrand(mpo,mps,site,'right',y{ind},'normal');
    end    
    y{ind} = reshape(y{ind},[D*x*D,1]);
    listn = [listn,norm(y{ind})];
end
for j = 1:2*D*x*D,
% disp(j)
    if j>1, y{j} = y{j}-fast_projection(Q,y{j},D,x);end
    q = y{j}/norm(y{j});
    if j==1, Q = q;else Q = [Q q];end
    omega{j+r} = randn(D,x,D)+1i*randn(D,x,D);
    Aomega = mpoxrand(mpo,mps,site,'right',omega{j+r},'normal');
    for ii=1:iter
        Aomega = mpoxrand(mpo,mps,site,'right',Aomega,'dagger');
        Aomega = mpoxrand(mpo,mps,site,'right',Aomega,'normal');
    end    
    Aomega = reshape(Aomega,[D*x*D,1]);
    y{j+r} = Aomega-fast_projection(Q,Aomega,D,x);
    listn(j+r) = norm(y{j+r});
    for ind = j+1:j+r-1
        qaux = reshape(q',[D,x,D]);
        yaux = reshape(y{ind},[D,x,D]);
        y{ind} = y{ind}-q*contracttensors(qaux,3,[1,2,3],yaux,3,[1,2,3]);
        listn(ind) = norm(y{ind});
    end
    if j>=r,
        if max(listn(j+1:j+r)) <= epsilon/(10*sqrt(2/pi)), 
            break;
        end
    end
end
