%% Data import and wrangling 

Path = [path 'Collaborators/EmiliLab/EMT.TGFB_scProt_1/001-SCPipeline_Output/'];

[dat, txt] = xlsread( [Path 'EpiToMesen.TGFB.nPoP_trial1_1PercFDR.xlsx']);

[~, IDs] = xlsread( [Path 'cellIDToTimepoint.xlsx']);
ind.d0 = IDs( ismember( IDs(:,2), 'd0' ), 1 );
ind.d3 = IDs( ismember( IDs(:,2), 'd3' ), 1 );
ind.d9 = IDs( ismember( IDs(:,2), 'd9' ), 1 );

Cell_Ids = txt(1,2:end);
ins = zeros( size(Cell_Ids) );
ins ( ismember( Cell_Ids, ind.d0 ) ) = 1;
ins ( ismember( Cell_Ids, ind.d3 ) ) = 2;
ins ( ismember( Cell_Ids, ind.d9 ) ) = 3;

%% Comute pairwise protein correlations for days 0, 3, and 9
r1 = pear ( dat(:, ins == 1 )'  );
r2 = pear ( dat(:, ins == 2 )'  );
r3 = pear ( dat(:, ins == 3 )'  );

%% Compute correlartions between correlation vectors corresponding to the same protein
rr1 = pear2i( r1, r3);

in = rr1 < -0.1;

rr2 = pear2i( r1(in,in), r3(in,in) );

subplot(1,2,1), hist( rr1, 20 );
subplot(1,2,2), hist( rr2, 20 );

%% Visually inspect the agreement between pairwise protein correlations 
%  for the proteins in set in

plot( re(r1(in,in)), re(r3(in,in)), '.' )

%% Select proteins from both tails of the distribution rr1
% in set Tails

Tails = rr1 < -0.15 | rr1 > 0.6;

% Compute permutaion Clu based in the avarage correlations
Clu = clu( r1(Tails,Tails)+r3(Tails,Tails), 3 ); rb(1)

%% Visualize the clustered correlations 
% from set Tails

% Make correlation matrices having just Tails proteins 
r1_in = r1(Tails,Tails);
r2_in = r2(Tails,Tails);
r3_in = r3(Tails,Tails);

% Plot the correlations for all proteins ordred by permutaion Clu 
imagesc( [r1_in(Clu,Clu) zeros(sum(Tails),5) r2_in(Clu,Clu) zeros(sum(in),5) r3_in(Clu,Clu)]  ); 
set(gca, 'FontSize', 16 );
rb(1);

ht(1) = text(20,-2, 'Day 0' );
ht(2) = text(90,-2, 'Day 3' );
ht(3) = text(160,-2, 'Day 9' );
set(ht, 'FontSize', 22 );

%%
pdf( [path 'NS/PDFs/EMT_correlations'], [12 7], 1 );



%% Plot change in protein abundance



lowTail = rr1 < -0.15; 
Clu = clu ( [dat(lowTail, ins == 1 )  zeros(sum(lowTail),5) ...
             dat(lowTail, ins == 2 )  zeros(sum(lowTail),5) ...
             dat(lowTail, ins == 3 ) ], 3  ); 
set(gca, 'FontSize', 16 );
rb(1);


%%



































