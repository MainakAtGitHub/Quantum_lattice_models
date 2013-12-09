function [iRange, jRange] = find_lattice_translation_index(N, nOrbitals, iCell, jCell)

iRange = (((iCell(1)-1)*N + iCell(2) - 1)*nOrbitals + 1):((iCell(1)-1)*N + iCell(2))*nOrbitals;
jRange = (((jCell(1)-1)*N + jCell(2) - 1)*nOrbitals + 1):((jCell(1)-1)*N + jCell(2))*nOrbitals;