function cellvector=periodic_latticevectors(cellvector,N)
% so far only square systems; to do make it work for non-square
% (rectangular) systems
cellvector(1)=mod(cellvector(1)-1,N)+1;
cellvector(2)=mod(cellvector(2)-1,N)+1;