function ubs_bmp_VASP
% Prepare a binary file to plot the undolded band structure as a bitmap figure
%
% Execution:
%    octave ubs_bmp_VASP
%
% Generated file:
%    WAVECAR_spin1.f2b.bin
%
% Band structure plotting:
%    gnuplot f2b-band-structure.plt
%
% (c) Oleg Rubel modified 05 Dec 2018

%% Init. parameters
KPATH = [0 0 0;...
         1.0/2 0 0; ...
         1.0/3 1.0/3 0;...
         0.0 0.0 0.0]; % k-point path
FOLDS = [4 4 1 ]; % multiplicity in the corresponding directions used when constructing the super-cell
KLABEL = {'G'; 'M'; 'K'; 'G'};
finpt = './WAVECAR_spin1.f2b'; % input file name
Ef = -1.30469679; % Fermi energy (Ry)
ERANGE = [Ef-12.0 Ef+10.0]; % energy range for plot (eV)
pwr = 1/1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.
         % 0 - folded bands (needs wth = 0)
nK = 280; % pixels along k
nE =  400; % pixels along Energy axis
% sK = 0.0025; % smearing factor in k-space
% sE = 0.025; % smearing factor in energy
sK = 0.02; % smearing factor in k-space
sE = 0.2; % smearing factor in energy
G = [0.101273670  0.058470381  0.000000000
0.000000000  0.116940762  0.000000000
 0.000000000  0.000000000  0.050000000];    % Reciprocal latt. vect. from *.outputkgen

%% INITIALIZATION
[KEIG, EIG, W] = readinput(finpt); % read input data from file
% EIG - energy eigenvalues
% KEIG - k-list for eigenvalues
% W - list of characters

%% MAIN
L = [];
ENE = [];
WGHT = [];
%G = G'; % transpose G matrix (need for Wien2k)
for i=1 : 3
    G(i,:)=G(i,:)*FOLDS(i); % rescale reciprocal lattice vectors 
end                         % from supercell to primitive cell
dl = 0; % cumulative length of the path
KPATH = coordTransform(KPATH,G);
KEIG = coordTransform(KEIG,G);
XTICKS = [0];
for ikp = 1 : size(KPATH,1)-1
    B = KPATH(ikp,:) - KPATH(ikp+1,:);
    dk = sqrt(dot(B,B));
    XTICKS = [XTICKS; XTICKS(ikp)+dk];
    for j = 1 : length(EIG)
        if EIG(j) > ERANGE(1) && EIG(j) < ERANGE(2) && W(j) > 0
            dist = dp2l( KEIG(j,:) , KPATH(ikp,:) , KPATH(ikp+1,:) );
            if dist < eps % k-point is on the path
                A = KPATH(ikp,:) - KEIG(j,:);
                x = dot(A,B)/dk;
                if x >= 0  &&  x <= dk+eps % k-point is within the path range
                    L = [L; x+dl]; % append k-point coordinate along the path
                    ENE = [ENE; EIG(j)]; % append energy list
                    WGHT = [WGHT; W(j)];
                end
            end
        end
    end
    dl = dl + dk;
end
if isempty(L)
    msg = ['No eigenvalues are selected for the plot. ', ...
        'The likely reason is that the energy range is ', ...
        'too restrictive (check ERANGE), or no k-points are located ', ...
        'on the path selected (check KPATH)'];
    error(msg);
end


WGHTRS = rescale(WGHT,pwr);
Li = linspace(XTICKS(1),XTICKS(end),240);
Ei = linspace(ERANGE(1)-Ef,ERANGE(2)-Ef,200);
nLi = length(Li);
nEi = length(Ei);
Wi = zeros(nLi,nEi);
ENE = ENE-Ef;
[XTICKS(1) XTICKS(end)]
for l = 1:length(L) % loop over all unfolded points
    disp(100*l/length(L));
    for i = 1:nLi
    if abs(Li(i)-L(l)) > 3*sK
        continue; % skip if > 3*sigma
    end
        for j = 1:nEi
            if abs(Ei(j)-ENE(l)) > 3*sE
                continue; % skip if > 3*sigma
            end
            wK = gauss(Li(i),L(l),sK);
            wE = gauss(Ei(j),ENE(l),sE);
            Wi(i,j) = Wi(i,j) + wK*wE*WGHTRS(l);
        end
    end
end
Wi = Wi/(gauss(0,0,sK)*gauss(0,0,sE)); % normalize weights
save_binary_matrix([finpt,'.bin'],Li,Ei,Wi');
XTICKS
KLABEL

%set(gca,'XTick',XTICKS);
%set(gca,'XTickLabel',KLABEL);

