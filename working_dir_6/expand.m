function [mps] = expand(mps0,D)

N = size(mps0,2);
[D0,~,d] = size(mps0{1});

mps = cell(1,N);
for ii = 1:N
    mps{ii} = zeros(D,D,d);
    mps{ii}(1:D0,1:D0,1:d) = mps0{ii};
end
