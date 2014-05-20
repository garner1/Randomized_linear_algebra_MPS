function [mps_out,ee_out] = outer_loop(V1,mpo1,U1,V2,mpo2,U2,mpsm,m_in,m_out,precision,N)

num = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V1,U1,mpo1,mpsm),3,[1,2,3]);
den = contracttensors(conj(mpsm),3,[1,2,3],matvectn(V2,U2,mpo2,mpsm),3,[1,2,3]);
ee = real(num/den);

list = zeros(50,1);
mps = cell(50,1);
list(1) = ee;
mps{1} = mpsm;

ii = 1;
while 1,
    ii = ii+1;
    [mps{ii},list(ii)] = inner_loop(V1,mpo1,U1,V2,mpo2,U2,mps{ii-1},list(ii-1),m_in);
%   plot(list(1:ii)/N);drawnow
    if ii>m_out && (std(list(ii-m_out:ii))/abs(mean(list(ii-m_out:ii)))<precision), 
    	mps_out=mps{ii};ee_out=list(ii);
    	if list(ii)>list(1), mps_out = mps{1}; ee_out = list(1); end	
	break;
    end
    if ii>=50, mps_out = mps{1}; ee_out = list(1);break; end
end

%     disp([std(list)/abs(mean(list)),ee/N])

