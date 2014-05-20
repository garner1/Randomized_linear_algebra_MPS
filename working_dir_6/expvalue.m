function value = expvalue(mpo,mps,Vh,Uh,Vn,Un)

num = contracttensors(conj(mps),3,[1,2,3],matvectn(Vh,Uh,mpo,mps),3,[1,2,3]);
den = contracttensors(conj(mps),3,[1,2,3],matvectn(Vn,Un,reshape(eye(2),[1 1 2 2]),mps),3,[1,2,3]);
value = real(num/den);
