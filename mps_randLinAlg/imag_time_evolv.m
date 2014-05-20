clear;
N=10; D=2;
nreal=5;
beta=10;steps=100;
precision=1e-5; dt=beta/steps; jflipped=5;

% magnetization in z-direction
oset=cell(1,N);
sx=0*[0,1;1,0]; sy=0*[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2); 
for j=1:N, oset{1,j}=id; end;
oset{1,jflipped}=sz;

% time evolution operator
h=kron(sx,sx)+kron(sy,sy)+kron(sz,sz); %[d^2,d^2]~[(s1,s2);(s1',s2')]
w=expm(-0.5*dt*h);  %[d^2,d^2]~[(s1,s2);(s1',s2')]
w=reshape(w,[2,2,2,2]); %[s1,s2,s1',s2']
w=permute(w,[1,3,2,4]); %[s1,s1',s2,s2']
w=reshape(w,[4,4]); %[(s1,s1'); (s2,s2')]
[U,S,V]=svd2(w);
eta=size(S,1);
U=U*sqrt(S);
V=sqrt(S)*V;
U=reshape(U,[2,2,eta]); %[s1,s1',k]
U=permute(U,[4,3,2,1]); %add a dummy index to make a proper mpo [:,k,s1',s1]
V=reshape(V,[eta,2,2]); %[k,s2,s2']
V=permute(V,[1,4,3,2]); %[:,k,s2',s2]
I=reshape(id,[1,1,2,2]);
mpo_even=cell(1,N);
mpo_odd=cell(1,N);
for j=1:N, mpo_even{j}=I; mpo_odd{j}=I; end
for j=1:2:(N-1), mpo_odd{j}=U; mpo_odd{j+1}=V; end 
for j=2:2:(N-1), mpo_even{j}=U; mpo_even{j+1}=V; end

% starting state (one spin flipped) 
% mps0=cell(1,N);
% for j=1:N
%     if j==jflipped, state=[0; 1]; else state=[1; 0]; end
%     mps0{j}=reshape(state,[1,1,2]); 
% end
mzvalues=zeros(nreal,steps);
for indr=1:nreal
    disp(indr)
mps0=randmps_norm(N,D);
% time evolution 
mps=mps0; 
for step=1:steps
    [mps]=reduceD(mps,mpo_even,D,precision); 
    [mps]=reduceD(mps,mpo_odd,D,precision); 
    mzvalues(indr,step)=real(expectationvalue(mps,oset));
end
end
plot(mean(mzvalues,1))