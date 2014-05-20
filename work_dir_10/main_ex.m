%cellfun
%Heisenberg
%
clear all;
N = 5*50;
D = 2^5;
size_bulk = 50;
site_in = 1;
iter = 0;
sweeps = round(5*N/size_bulk);
m_in = 10;
        
phi = 0.0;
J = 0.25;
jx = 1.0*J;
jy = 1.0*J;
jz = 1.0*J;
hz = 0; h = hz;
% tic
[elist1,mps1] = heisenberg(N,D,size_bulk,site_in,iter,m_in,sweeps,jx,jy,jz,h,phi,'fast');
% toc
% %
% % Ising
% % example: 
% % N = 2^2*2^6;
% % D = 2^5;
% % size_bulk= 2^6;
% % sweeps = round(5*N/size_bulk);
% % precision = 1e-5; 
% % e = -1.273233707327088
% N = 2^2*2^6;
% D = 2^5;
% size_bulk= 2^6;
% sweeps = round(5*N/size_bulk);
% precision = 1*1e-5; %used in the outer loop of the solver
% 
% [elist2,mps2] = ising_critical(N,D,size_bulk,sweeps,precision);
