function [mpsma,mpsmb,ee] = outer_loop(V1,mpo1a,mpo1b,U1,V2,mpo2,U2,mpsma,mpsmb,ee,m_in,m_out)

for ii = 1:m_out
    [mpsma,mpsmb,ee] = inner_loop(V1,mpo1a,mpo1b,U1,V2,mpo2,U2,mpsma,mpsmb,ee,m_in);
end


