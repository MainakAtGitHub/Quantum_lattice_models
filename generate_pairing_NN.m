function [ Gamma,latticeVectorsSC ] = generate_pairing_NN( pairing )
% script to set up pairing interactions with NN pairing
% argument: pairing gives the value of the NN pairing interaction
% output: pairing interaction and lattice vectors for the last entry
format long;
latticeVectorsSC = [-1,0; 0,-1; 0,0; 0,1; 1,0];

% allow the argument "pairing" to be a vector with multiple data, i.e. do a
% loop here

for ind=1:numel(pairing)
    Gamma(:,:,1:5) = pairing(ind);
    Gamma(:,:,3) = 0;
    pairingString = num2str(pairing(ind));
    filename = ['Gamma_',pairingString,'_NN.mat'];
    % not a good idea to hardcode paths that might not exist on other computers
    % or for other users, better to just use the local directory, or some
    % default directory which should exist on all systems, for example /tmp/
    % filename = sprintf('/home/UFAD/mainak.pal/Desktop/summer_2019_office/andreas/BdG/Tutorial/matlabtest/Gamma_%s_NN.mat',pairingString);
    % this is probably not what is intended, because also the "irrelevant"
    % variables 'filename', 'pairing', 'pairingString' would end up saved in
    % the file
    % save(filename);
    save(filename,'Gamma','latticeVectorsSC');
    % clearing any variables at the end of a function has no meaning because
    % with the function to end, all variables are cleared anyhow!
    % clear pairing pairingString filename
end

end

