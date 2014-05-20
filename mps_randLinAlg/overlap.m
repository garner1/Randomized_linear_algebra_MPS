function overlap=overlap(mps1,mps2)

%evaluates <MPS1|MPS2> [hence the norm is sqrt(overlap)] for open boundary conditions
%the complexity is O(N*D^3*d)
%!!!!REMEMBER TO PROVIDE INPUT STATE THAT ARE NORMALIZED!!!!

N1=length(mps1);N2=length(mps2);
d=size(mps1{1},3);
if N1~=N2, error('the two MPSs do not have the same number of qubits'); end

overlap=1;
X=eye(d); X=reshape(X,[1,1,d,d]); 
for j=N1:-1:1
    overlap=updateCright(overlap,mps1{j},X,mps2{j}); 
end


