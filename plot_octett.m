% lowres
energy=[0.001:0.002:0.03];
% highres
energy=[0.001:0.003:0.036];
[kx,ky]=banana_1band('~/docs/real_space/BdG/bscco/Z3/input_SC_U015_N35_Z3.txt',energy);
symm=true;
[o,figure]=plot_octett_1band(kx,ky,energy,[],symm);
print_pdf('/tmp/bscco_octett2.pdf')
% make individual plots for each energy;
for n=1:length(energy);
    [kx,ky]=banana_1band('~/itp/docs/real_space/BdG/bscco/Z3/input_SC_U015_N35_Z3.txt',energy(n));
    [o,figureq]=plot_octett_1band(kx,ky,energy(n),[],symm);
    hoffman_map(figureq)
    %[o,figure]=octett_1band(kx,ky,energy(n),symm);
    print_pdf(['/tmp/bscco_octett_n_',sprintf('%03d',n)]);
    close all;
end;