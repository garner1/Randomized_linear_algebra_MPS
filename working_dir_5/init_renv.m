function [U] = init_renv(mps,mpo,U_env,tag)
N = size(mpo,2);
U = cell(1,N);
if strcmp(tag,'two'),
    U{N/2} = U_env;
    for ind = N/2:-1:2
        U{ind-1} = contracttensors(conj(mps{ind}),3,2,U{ind},4,1);
        U{ind-1} = contracttensors(U{ind-1},5,[2,3],mpo{ind},4,[3,2]);
        U{ind-1} = contracttensors(U{ind-1},5,[5,2],mps{ind},3,[3,2]);
        U{ind-1} = permute(U{ind-1},[1 3 4 2]);
    end
else
    U{N} = U_env;
    for ind = N:-1:2+N/2
        U{ind-1} = contracttensors(conj(mps{ind}),3,2,U{ind},4,1);
        U{ind-1} = contracttensors(U{ind-1},5,[2,3],mpo{ind},4,[3,2]);
        U{ind-1} = contracttensors(U{ind-1},5,[5,2],mps{ind},3,[3,2]);
        U{ind-1} = permute(U{ind-1},[1 3 4 2]);
    end
end
    