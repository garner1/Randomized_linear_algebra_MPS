function v = matvectn(L,R,mpox,v0)

v = contracttensors(L,4,4,v0,3,1);

v = contracttensors(v,5,[3,5],mpox,4,[1,4]);

v = contracttensors(v,5,[1,4,3],R,4,[4,2,3]);

v = permute(v,[1,3,2]);


