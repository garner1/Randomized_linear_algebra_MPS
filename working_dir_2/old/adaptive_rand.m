%
% algo 4.2 in REF
%
clear all;
r = 10;epsilon = 0.1;
D = 100;
A = randn(D);
omega = randn(D,r); 
y = zeros(D,r);
listn = [];
y = A*omega;
for col = 1:r
    listn = [listn,norm(y(:,col))];
end
j = 0;
Q = zeros(D,1);
while max(listn(j+1:j+r)) > epsilon/(10*sqrt(2/pi))
    j = j+1;
    if j>1, y(:,j) = (eye(D)-Q*Q')*y(:,j);end
    q(:,j) = y(:,j)/norm(y(:,j));
    Q = [Q q(:,j)];
    omega(:,j+r) = randn(D,1);
    y(:,j+r) = (eye(D)-Q*Q')*A*omega(:,j+r);
    listn(j+r) = norm(y(:,j+r));
    for ind = j+1:j+r-1
        y(:,ind) = y(:,ind)-q(:,j)*q(:,j)'*y(:,ind);
        listn(ind) = norm(y(:,ind));
    end
disp(listn(j+1:j+r))
end
