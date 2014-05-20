function vec = mpoxrand(mpo,mps,dir,vec,dagger)

N = size(mpo,2);
[x,~,~,~] = size(mpo{1});
[m,~,~] = size(mps{1});
if strcmp(dagger,'normal'),
    if strcmp(dir,'right'),
        for ind = N:-1:1   
            T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
            T = contracttensors(T,5,[2,5],vec,3,[2,3]);
            T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
            vec = permute(T,[3,1,2]);
        end
    end
    if strcmp(dir,'left'),
        vec = reshape(vec,[m,x,m]);
        for ind = 1:N    
            T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
            T = contracttensors(T,5,[1,4],vec,3,[2,3]);
            T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
            vec = permute(T,[3,1,2]);
        end
    end
end
if strcmp(dagger,'dagger'),
    if strcmp(dir,'right'),
        for ind = 1:N
            mpo{ind} = permute(conj(mpo{ind}),[2,1,3,4]);
            mps{ind} = permute(conj(mps{ind}),[2,1,3]);
            T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
            T = contracttensors(T,5,[2,5],vec,3,[2,3]);
            T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
            vec = permute(T,[3,1,2]);
        end
    end
    if strcmp(dir,'left'),
        vec = reshape(vec,[m,x,m]);
        for ind = N:-1:1  
            mpo{ind} = permute(conj(mpo{ind}),[2,1,3,4]);
            mps{ind} = permute(conj(mps{ind}),[2,1,3]);
            T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
            T = contracttensors(T,5,[1,4],vec,3,[2,3]);
            T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
            vec = permute(T,[3,1,2]);
        end
    end
end

    

