function [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell)
% find the index to read out delta matrix
% N: number of cells for calculation
% nOrbitals: clear!
% iCell: cell for first index i
% jCell: cell for second index j
if numel(N)==1
   N=[N,N];
end
iRange = (((iCell(1)-1)*N(2) + iCell(2) - 1)*nOrbitals + 1):((iCell(1)-1)*N(2) + iCell(2))*nOrbitals;
jRange = (((jCell(1)-1)*N(2) + jCell(2) - 1)*nOrbitals + 1):((jCell(1)-1)*N(2) + jCell(2))*nOrbitals;
