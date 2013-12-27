function f=make_plot_impurity_dos(inputfile,smooth)
if nargin <1
    inputfile='LDOS_FeSe_Milan_Gamma_Vimp_4_N_9_M_40_ita_0.003'
end;
if nargin < 2
    if isempty(strfind(inputfile, 'tetra'))
        smooth=0;
    else
        %default: smoothing if "tetra" in filename
        smooth=1;
    end
end
range=[-0.2 0.2];
orb={'orbital1','orbital2','orbital3','orbital4','orbital5','total'};

% load input
load(inputfile,'-mat')
% some smoothing if necessary
orbitalLDOSFarAway=sg_smooth(orbitalLDOSFarAway,smooth);
orbitalLDOSImp=sg_smooth(orbitalLDOSImp,smooth);
orbitalLDOSImpNN=sg_smooth(orbitalLDOSImpNN,smooth);
orbitalLDOSImpNNN=sg_smooth(orbitalLDOSImpNNN,smooth);
% fix negative values (numerical error, smoothing artefacts)
orbitalLDOSFarAway(orbitalLDOSFarAway<0)=0;
orbitalLDOSImp(orbitalLDOSImp<0)=0;
orbitalLDOSImpNN(orbitalLDOSImpNN<0)=0;
orbitalLDOSImpNNN(orbitalLDOSImpNNN<0)=0;
% DOS far away (without impurity)
figure1=figure;
plot1=plot(energy,[orbitalLDOSFarAway;sum(orbitalLDOSFarAway)]);
setlabels(plot1,orb,range);
print_pdf('/tmp/LDOS_far_away.pdf');

% impurity DOS
figure2=figure;
plot2=plot(energy,[orbitalLDOSImp;sum(orbitalLDOSImp)]);
setlabels(plot2,orb,range);
print_pdf('/tmp/LDOS_Imp.pdf');

% NN dos
figure3=figure;
plot3=plot(energy,[orbitalLDOSImpNN;sum(orbitalLDOSImpNN)]);
setlabels(plot3,orb,range);
print_pdf('/tmp/LDOS_Imp_NN.pdf');

% NNN dos
figure4=figure;
plot4=plot(energy,[orbitalLDOSImpNNN;sum(orbitalLDOSImpNNN)]);
setlabels(plot4,orb,range);
print_pdf('/tmp/LDOS_Imp_NNN.pdf');

% compare total dos
figure5=figure;
plot5=plot(energy,[sum(orbitalLDOSFarAway);sum(orbitalLDOSImp);sum(orbitalLDOSImpNN);sum(orbitalLDOSImpNNN)]);
set(plot5(1),'DisplayName','tot far away');
set(plot5(2),'DisplayName','tot impurity');
set(plot5(3),'DisplayName','tot NN');
set(plot5(4),'DisplayName','tot NNN');
xlim(range);
xlabel({'\omega'});

% Create ylabel
ylabel({'DOS [1/eV]'});

% Create legend
legend show
print_pdf('/tmp/LDOS_tot.pdf');

end

function setlabels(plot,orb,range)
for i=1:6
set(plot(i),'DisplayName',orb{i});
end;
xlim(range);
% Create xlabel
xlabel({'\omega'});

% Create ylabel
ylabel({'DOS [1/eV]'});

% Create legend
legend show
end
