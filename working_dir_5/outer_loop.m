function [mpsm,ee] = outer_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,m_in,m_out)

num = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V1,U1,mpo1,mpsm),3,[1,2,3]);
den = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V2,U2,mpo2,mpsm),3,[1,2,3]);
ee = num/den;
elist = [];
for ii = 1:m_out
    [mpsm,ee] = inner_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,ee,m_in);
%     elist = [elist,real(ee)];
%     plot(elist);drawnow;
%     disp(ee)
end


