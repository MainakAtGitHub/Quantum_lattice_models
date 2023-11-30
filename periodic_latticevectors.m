function cellvector=periodic_latticevectors(cellvector,N)
% so far only square systems; to do make it work for non-square (ready for rectangular systems, 23 Dec 2021)
% (rectangular) systems
if numel(N)==1
    N=[N,N];
end
cellvector(1)=mod(cellvector(1)-1,N(1))+1;
cellvector(2)=mod(cellvector(2)-1,N(2))+1;
