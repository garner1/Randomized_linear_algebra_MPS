function V_env = update_lenv(mpsA,mpox,V_env)


V_env = contracttensors(conj(mpsA),3,1,V_env,4,2);
V_env = contracttensors(V_env,5,[2,4],mpox,4,[3,1]);
V_env = contracttensors(V_env,5,[5,3],mpsA,3,[3,1]);
V_env = permute(V_env,[2 1 3 4]);