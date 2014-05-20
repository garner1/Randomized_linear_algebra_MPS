% mps = randmps_norm(N,D) PREPARE A LEFT NORMALIZED INHOMOGENEOUS RANDOM MPS STATE. 

% THE INPUT ARGUMENTS ARE: 
% N = THE NUMBER OF QUBITS AND 
% D = THE MAX AUX DIMENSION

% THE OUTPUT IS:
% mps = IS A CELL WITH N ENTRIES, EACH ONE A TENSOR OF 3 INDICES, WITH THE
% LAST INDEX BEING THE PHYSICAL INDEX; <mps|mps> = 1;

% TO HAVE THE STATE IN LEFT CANONICAL FORM THE A MATRICES ARE TAKEN FROM THE FIRST \CHI COLUMN OF RANDOM UNITARIES

function [mps] = randmps_norm(N,D)

% initialize the mps cell
mps = cell(1,N); 

% prepare the fhe other A matrices except for the last one
U = randU(2*D);
U = U(:,1:D);
U = reshape(U,[2,D,D]);
for ind = 1:N
    mps{ind} = permute(U,[2,3,1]);
end
