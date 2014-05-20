%
% algo 4.2 in REF
%
clear all;
r = 100;
D = 100;
A = randn(D);
Q = [];
nlist = [];
for ind = 1:r
    omega = randn(D,1)/sqrt(D);
    y(:,ind) = A*omega;
    if ind ==1, 
        qq(:,ind) = y(:,ind);
        Q = [Q qq(:,ind)/norm(qq(:,ind))];
    else
        qq(:,ind) = (eye(D)-Q*Q')*y(:,ind);
        nlist = [nlist,norm(qq(:,ind))];
        disp(norm(qq(:,ind)))
        Q = [Q qq(:,ind)/norm(qq(:,ind))];
        plot(nlist,'x');drawnow;
    end
end

