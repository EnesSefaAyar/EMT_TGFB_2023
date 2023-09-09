% Scripts modified from EMT_2
% Forked on Sept 9th
% Goal: Load DART-ID updated data and compute correlations without imputaion

%% Section 1: Data import and wrangling 
% =========================================

Path = [path 'Collaborators/EmiliLab/EMT.TGFB_scProt_1/001-SCPipeline_Output/'];

[dat_L, txt] = xlsread( [Path 'EpiToMesen.TGFB.nPoP_trial1_1PercDartFDRTMTBulkDIA.WallE_unimputed.xlsx']);
Cell_Ids = txt(1,2:end);
ProtIDs_L = txt(2:end,1);

[~, IDs] = xlsread( [Path 'cellIDToTimepoint.xlsx']);
ind.d0 = IDs( ismember( IDs(:,2), 'd0' ), 1 );
ind.d3 = IDs( ismember( IDs(:,2), 'd3' ), 1 );
ind.d9 = IDs( ismember( IDs(:,2), 'd9' ), 1 );


ins = zeros( size(Cell_Ids) );
ins ( ismember( Cell_Ids, ind.d0 ) ) = 1;
ins ( ismember( Cell_Ids, ind.d3 ) ) = 2;
ins ( ismember( Cell_Ids, ind.d9 ) ) = 3;

%% Filter proteins to a subset with less missing data

in1 = sum(~isnan(dat_L),2) > 30;
in2 = sum( ~isnan(dat_L(:, ins==1 )), 2) > 4 & ...
      sum( ~isnan(dat_L(:, ins==2 )), 2) > 4 & ...
      sum( ~isnan(dat_L(:, ins==3 )), 2) > 4;
in3 = (in1 + in2)>1 ;

fprintf( 'Selected Proteins: %d\n\n', sum( in3 ) );
ProtIDs = ProtIDs_L(in3,:);
dat_nan = dat_L(in3,:);


%% Without IMPUTAION comute pairwise protein correlations for days 0, 3, and 9
dat_Lin = dat_L(in3,:);
r1 = pearNaN ( dat_Lin(:, ins == 1 )'  );
r2 = pearNaN ( dat_Lin(:, ins == 2 )'  );
r3 = pearNaN ( dat_Lin(:, ins == 3 )'  );


%% oOad SwissProt Annotation 
load ( [path 'bin/Human_ONLYsp_Swiss_2_GeneName'] );
[Lia,Locb] = ismember( ProtIDs,  HumanONLYspSwiss2GeneName(:, 1 ) );
%%



%% Section 2: Identify and plot correlation CLUSTERS
% ====================================================


rd = r3 - r1;

% Select proteins with complete pairwise correlations
Mean = mean ( isnan( rd ) );
in_NAN =  find( Mean == 0 ); 
Rd = rd( in_NAN, in_NAN );


% Select correlations without missing data and large difference 
SumR = sum( Rd.^2 );
in = SumR > prctile( SumR, 50 ); 
fprintf( 'Number of selected proteins: %d, \n\n', sum(in) );
ProtIDs2 = ProtIDs( in_NAN(in) );

% Cluster the matrix of correlation differences 
Clu = clu( Rd(in,in), 3); rb(1);
%%

In = in_NAN(in);
In = In(Idx); %(Clu);


Lim = 1;
subplot(2,3,1)
imagesc( r1(In,In)); 
set( gca, 'CLim', [-Lim Lim] );
set( gca, 'Xtick', [], 'Ytick', [] );
h(1) = title( 'Day 0' );


subplot(2,3,2)
imagesc( r2(In,In)); 
set( gca, 'CLim', [-Lim Lim] );
set( gca, 'Xtick', [], 'Ytick', [] );
h(2) = title( 'Day 3' );

subplot(2,3,3)
imagesc( r3(In,In)); hb(1) = rb(Lim); 
set( gca, 'Xtick', [], 'Ytick', [] );
h(3) = title( 'Day 9' );


% Differences 
Lim = 0.7;
subplot(2,3,4)
imagesc( r2(In,In) - r1(In,In) );
set( gca, 'Xtick', [], 'Ytick', [] );
h(4) = title( 'Day 3 - Day 0' );

subplot(2,3,5)
imagesc( r2(In,In) - r1(In,In) ); 
set( gca, 'Xtick', [], 'Ytick', [] );
h(5) = title( 'Day 9 - Day 3' );

subplot(2,3,6)
imagesc( r3(In,In) - r1(In,In)  ); hb(2) = rb(Lim); 
set( gca, 'Xtick', [], 'Ytick', [] );
h(6) = title( 'Day 9 - Day 0' );


set(h,  'FontSize', 18 );
set(hb, 'FontSize', 12 );

set( hb(1), 'Ytick', [-1 0 1] );
set( hb(2), 'Ytick', [-0.7 0 0.7] );

%%
pdf( [path 'NS/PDFs/EMT_Global_correlation-Patterns_Kmeans'], [12 7], 1 );




%% Analyze protein clusters

%[Lia,Locb] = ismember( ProtIDs2,  HumanONLYspSwiss2GeneName(:, 1 ) );

fprf( ProtIDs2(Clu), 3, [path 'NS/txt/Prot_CorrClusters.txt']  );
fprf( ProtIDs2(flip(Clu,2)), 3, [path 'NS/txt/Prot_CorrClusters_flip.txt']  );




%% Cluster the matrix of correlation differences using K-means

Rd2 = Rd(in,in);
[idx,C] = kmeans( Rd2, 3); 

Idx = [find(idx==1); find(idx==2);  find(idx==3) ];

imagesc( Rd2(Idx,Idx) );  rb(1);

%%
fprf( ProtIDs2, 3, [path 'NS/txt/Prot_CorrClusters_Background.txt']  );
fprf( ProtIDs2(idx==1), 3, [path 'NS/txt/Prot_CorrCluster-1.txt']  );
fprf( ProtIDs2(idx==2), 3, [path 'NS/txt/Prot_CorrCluster-2.txt']  );
fprf( ProtIDs2(idx==3), 3, [path 'NS/txt/Prot_CorrCluster-3.txt']  );




















