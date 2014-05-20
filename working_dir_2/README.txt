ALGORITHM
Variational optimization of energy for MPS with PBC. 


INPUT variables
Hamiltonian as MPO {X_1,...,X_N}
X_i = [l,r,b,t] of dim {5,5,d,d}
initial state {A_1,...,A_N}
A_i = [l,r,t] of dim {m,m,d}

STAGE A [cost O(N*p*m^3)]
repeat for #p of random vectors V
start from site #1
generate a random normal vector V of size [5*m^2,1] and reshape to [m,5,m]
contract from left X_{N} with P, costing O(m^3) 
repeat contraction till X_{2}
reshape final state to [5*m^2,1] column vector

collect #p of final vectors into a matrix of size Y=[5*m^2,p]
perform QR decomposition of Y at cost O(p^2*m^2): [Q,~]=qr(Y); Q[5*m^2,p]

STAGE B [cost O(N*p*m^3)]
construct B[p,5*m^2]=Q^dag * A, multiplying row by column
SVD(B,0)=[U1,S,V] at cost O(p^2*m^2)
U[5*m^2,k]=Q*U1
S[k,k]
V^dag[k,5*m^2]
constract S with V^dag; cost O(5*k^2*m^2)

Construct effective Hamiltonian[1]
input data: {U[5*m^2,k],SV[k,5*m^2],X_1[5,5,d,d]}
contract U[m,5,m,k] with X_1[5,5,d,d] into UX_1[m,m,k,5,d,d]
contract UX_1[m,m,k,5,d,d] with SV[m,5,m,k] into UXSV[m,m,d,d,m,m] 

Repeat the same construction for N_eff[1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
The most intensive computational part is the solution of the 
generalized eigenvalue problem. 
One can simplify it first by observing that both Heff and Neff are 
tensor product of a network and the local mpo, which is Id in the 
case of Neff. Hence the diagonalization in the local mpo sector 
is trivial, involving simply the diagonalization of the local mpo
operator. The network sector could be simplified if Neff is positive 
definite. Given the iterative nature in the solution of the problem one 
could adopt randomized linear algebra schemes in the solution of the 
generalized eigenvalue problem: like Nystrom transform, Lanczos or similar.





