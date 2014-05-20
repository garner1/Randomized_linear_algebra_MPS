clear all;
x = [];
y = [];
% for N = linspace(10,100,10)
    N = 20;
    DB = 2^7;
    d= 2;
    precision = 1e-5;

    mpoX = mpo_hh(N,1,1,1,0,0);
    mpsA = createrandommps(N,DB,d); 
%     tic
    mps = mpomps(mpoX,mpsA);
%     toc

    tic
    mpsB1 = truncate(mps,'rl',DB,'fastsvd'); 
    t(1) = toc
    tic
    mpsB2 = truncate(mps,'rl',DB,'svd');
    t(2) = toc
    x = [x,N];
    y = [y,t(2)-t(1)]
%     plot(x,y);drawnow;
% end


