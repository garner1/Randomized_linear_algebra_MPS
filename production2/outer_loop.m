function [mpsma,mpsmb,ee] = outer_loop(V1,mpo1a,mpo1b,U1,V2,mpo2,U2,mpsma,mpsmb,m_in,m_out,precision,N)

[D,~,d] = size(mpsma);

mpsm = permute(contracttensors(mpsma,3,2,mpsmb,3,1),[1,3,2,4]);
mpsm = reshape(mpsm,[D,D,d*d]);

mpo1aux = contracttensors(mpo1a,4,2,mpo1b,4,1);
mpo1aux = permute(mpo1aux,[1,4,2,5,3,6]);
mpo1aux = reshape(mpo1aux,[5,5,d*d,d*d]);

mpo2aux = contracttensors(mpo2,4,2,mpo2,4,1);
mpo2aux = permute(mpo2aux,[1,4,2,5,3,6]);
mpo2aux = reshape(mpo2aux,[1,1,d*d,d*d]);

num = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V1,U1,mpo1aux,mpsm),3,[1,2,3]);
den = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V2,U2,mpo2aux,mpsm),3,[1,2,3]);
ee = real(num/den);
list = ee;
ii = 0;
while 1,
    ii = ii+1;
    [mpsma,mpsmb,ee] = inner_loop(V1,mpo1a,mpo1b,U1,V2,mpo2,U2,mpsma,mpsmb,ee,m_in);
    list = [list,ee];
%     disp([std(list)/abs(mean(list)),ee/N])
    plot(list/N);drawnow
    if ii>m_out && (std(list(end-m_out:end))/abs(mean(list(end-m_out:end)))<precision), break;end
    if ii>50, break;end
end



