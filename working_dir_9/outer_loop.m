function [mpsm,ee,list] = outer_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,m_in,m_out,precision,N)
% outer_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,m_in,m_out,precision,ee)

num = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V1,U1,mpo1,mpsm),3,[1,2,3]);
den = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V2,U2,mpo2,mpsm),3,[1,2,3]);
ee = real(num/den);
list = ee;
ii = 0;
while 1,
    ii = ii+1;
    [mpsm,ee] = inner_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,ee,m_in);
    list = [list,ee];
%     disp([std(list)/abs(mean(list)),ee/N])
%     plot(list/N);drawnow
    if ii>m_out && (std(list(end-m_out:end))/abs(mean(list(end-m_out:end)))<precision), break;end
%     if ii>50 && list(end)>list(1), break;end
    if ii>100, break;end
end

