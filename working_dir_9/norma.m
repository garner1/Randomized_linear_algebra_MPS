function out = norma(mps,u,v)

% [D,x,D,pn] = size(u);

out = contracttensors(conj(mps),3,2,u,4,1);
out = contracttensors(out,5,[2,4],mps,3,[3,2]);
out = contracttensors(out,4,[1,2,3,4],v,4,[2,3,1,4]);


