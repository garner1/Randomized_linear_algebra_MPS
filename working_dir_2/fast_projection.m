function out = fast_projection(Q,y,D,x)

[m,n] = size(Q);
out = zeros(m,1);
for col = 1:n
    aux = reshape(Q(:,col)',[D,x,D]);
    y = reshape(y,[D,x,D]);
    out = out + contracttensors(aux,3,[1,2,3],y,3,[1,2,3])*Q(:,col);
end
