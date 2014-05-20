function p = mpoxrand(mpo,mps,site,dir,p,dagger)

N = size(mpo,2);
[x,~,~,~] = size(mpo{site});
[m,~,~] = size(mps{site});
if strcmp(dagger,'normal'),
    if site==1, 
        if strcmp(dir,'right'),
            for ind = N:-1:2   
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[2,5],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
                p = permute(T,[3,1,2]);
            end
        end
        if strcmp(dir,'left'),
            p = reshape(p,[m,x,m]);
            for ind = 2:N    
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[1,4],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
                p = permute(T,[3,1,2]);
            end
        end
    else
        if strcmp(dir,'right'),
            for ind = site-1:-1:1  
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[2,5],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
                p = permute(T,[3,1,2]);
            end
            for ind = N:-1:site+1
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[2,5],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
                p = permute(T,[3,1,2]);
            end
        end
        if strcmp(dir,'left'),
            p = reshape(p,[m,x,m]);
            for ind = site+1:N    
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[1,4],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
                p = permute(T,[3,1,2]);
            end
            for ind = 1:site-1    
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[1,4],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
                p = permute(T,[3,1,2]);
            end
        end
    end
end
if strcmp(dagger,'dagger'),
    if site==1, 
        if strcmp(dir,'right'),
            for ind = 2:1:N
                mpo{ind} = permute(conj(mpo{ind}),[2,1,3,4]);
                mps{ind} = permute(conj(mps{ind}),[2,1,3]);
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[2,5],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
                p = permute(T,[3,1,2]);
            end
        end
        if strcmp(dir,'left'),
            p = reshape(p,[m,x,m]);
            for ind = 2:N    
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[1,4],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
                p = permute(T,[3,1,2]);
            end
        end
    else
        if strcmp(dir,'right'),
            for ind = site+1:1:N
                mpo{ind} = permute(conj(mpo{ind}),[2,1,3,4]);
                mps{ind} = permute(conj(mps{ind}),[2,1,3]);
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[2,5],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
                p = permute(T,[3,1,2]);
            end
            for ind = 1:1:site-1
                mpo{ind} = permute(conj(mpo{ind}),[2,1,3,4]);
                mps{ind} = permute(conj(mps{ind}),[2,1,3]);
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[2,5],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,2]);
                p = permute(T,[3,1,2]);
            end
        end
        if strcmp(dir,'left'),
            p = reshape(p,[m,x,m]);
            for ind = site+1:N    
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[1,4],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
                p = permute(T,[3,1,2]);
            end
            for ind = 1:site-1    
                T = contracttensors(mpo{ind},4,4,mps{ind},3,3);
                T = contracttensors(T,5,[1,4],p,3,[2,3]);
                T = contracttensors(T,4,[2,4],conj(mps{ind}),3,[3,1]);
                p = permute(T,[3,1,2]);
            end
        end
    end
end

    

