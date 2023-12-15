function   [v , maxEig] = eigenvector_centrality_und_fab(CIJ)
%   EIGENVECTOR_CENTRALITY_UND      Spectral measure of centrality
%
%   v = eigenvector_centrality_und(CIJ)
%
%   EigenVector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector
%   corresponding to the largest eigenvalue of the adjacency matrix.
%
%   Inputs:     CIJ,        binary/weighted undirected adjacency matrix.
%
%   Outputs:      v,        eigenvector associated with the largest
%                           eigenvalue of the adjacency matrix CIJ.
%
%   Reference: Newman, MEJ (2002). The mathematics of networks.
%
%   Contributors:
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012
%   Mika Rubinov, University of Cambridge, 2015

%   MODIFICATION HISTORY
%   2010/2012: original (XNZ, RB)
%   2015: ensure the use of leading eigenvector (MR)



[V,D] = eigs(CIJ , 1 , 'largestreal');

% [M,idx] = max(diag(D));
ec = abs(V);
v = reshape(ec, length(ec), 1);

maxEig = abs(D);
end