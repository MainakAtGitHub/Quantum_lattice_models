% lowres
energy=[0.001:0.002:0.03];
% highres
energy=[0.001:0.0001:0.036];
[kx,ky]=banana_1band('bscco/Z3/input_SC_U015_N35_Z3.txt',energy);
octett_1band(kx,ky,energy);
print_pdf('/tmp/bscco_octett2.pdf')