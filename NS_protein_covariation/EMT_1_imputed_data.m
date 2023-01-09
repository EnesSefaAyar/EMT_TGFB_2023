%% Data import and wrangling
% Goal: Expore the EMT data using the complete data matrix exported from the SCoPE2 pipline
% Started Dec 26 2022 

Path = [path 'Collaborators/EmiliLab/EMT.TGFB_scProt_1/001-SCPipeline_Output/'];

[dat, txt] = xlsread( [Path 'EpiToMesen.TGFB.nPoP_trial1_1PercFDR.xlsx']);
Cell_Ids = txt(1,2:end);
ProtIDs = txt(2:end,1);

[~, IDs] = xlsread( [Path 'cellIDToTimepoint.xlsx']);
clear ind
ind.d0 = IDs( ismember( IDs(:,2), 'd0' ), 1 );
ind.d3 = IDs( ismember( IDs(:,2), 'd3' ), 1 );
ind.d9 = IDs( ismember( IDs(:,2), 'd9' ), 1 );


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

in = rr1 < -0.15;

rr2 = pear2i( r1(in,in), r3(in,in) );

subplot(1,2,1), hist( rr1, 20 );
subplot(1,2,2), hist( rr2, 20 );

%% Visually inspect the agreement between pairwise protein correlations
%  for the proteins in set in

subplot(1,2,1)
scatter( re(r1(in,in)), re(r3(in,in)), 'MarkerFaceColor', 'b' )

subplot(1,2,2)
clu( r1(in,in) ); rb(1);

%% Select example proteins with changing correlations

r1_norm = sum( r2.^2 ); %hist(r1_norm);

left_Tail = find( rr1 < -0.15 );
right_Tail = find( rr1 > 0.6 );

    [val_l, ind_l] = sort( sum( r1(left_Tail,left_Tail).^2 ), 'descend' );
    [val_r, ind_r] = sort( sum( r1(right_Tail,right_Tail).^2 ), 'descend' );

subplot(1,2,1)
imagesc( r1( right_Tail(ind_r), left_Tail(ind_l) ) ); rb(1);

subplot(1,2,2)
imagesc( r3( right_Tail(ind_r), left_Tail(ind_l) ) ); rb(1);

Rdiff = abs( r1( right_Tail(ind_r), left_Tail(ind_l) ) - r3( right_Tail(ind_r), left_Tail(ind_l) ));

[val, ind] = sort( Rdiff(:), 'descend' );
[x, y] = ind2sub(size(Rdiff), ind);
%%
[Lia,Locb] = ismember( ProtIDs,  HumanONLYspSwiss2GeneName(:, 1 ) );
%%
clear hh ht hl
for i = 24 % 14 9
    close all
    ind_x = right_Tail( ind_r(x(i) ));
    ind_y = left_Tail( ind_l(y(i) ));
xx = dat(ind_x, : );
yy = dat(ind_y, : );
color = 'k';

subplot(1,3,1)
    scatter( xx(ins == 1 ), yy(:, ins == 1 ), 'MarkerFaceColor', color ); set(gca, 'FontSize', 14)
    [corr1, p1] = corr( xx(ins == 1 )', yy(:, ins == 1 )' );
    ht(1) = title( ['Day 0, $\rho = ' sprintf( '%1.2f', corr1) '$'] );
    hl(1) = ylabel(regexprep( HumanONLYspSwiss2GeneName{Locb(ind_y),3}, '_HUMAN', '' ) );


subplot(1,3,2)
    scatter( xx(ins == 2 ), yy(:, ins == 2 ), 'MarkerFaceColor', color ); set(gca, 'FontSize', 14)
    corr2 = corr( xx(ins == 2 )', yy(:, ins == 2 )' );
    ht(2) = title( ['Day 3, $\rho = ' sprintf( '%1.2f', corr2) '$'] );
        hl(2) = xlabel( regexprep( HumanONLYspSwiss2GeneName{ Locb(ind_x) ,3}, '_HUMAN', '' ) );


subplot(1,3,3)
    scatter( xx(ins == 3 ), yy(:, ins == 3 ), 'MarkerFaceColor', color );  set(gca, 'FontSize', 14)
    corr3 = corr( xx(ins == 3 )', yy(:, ins == 3 )' );
    ht(3) = title( ['Day 9, $\rho = ' sprintf( '%1.2f', corr3) '$'] );

set( ht, 'FontSize', 22, 'Interpreter', 'latex' );
set( hl, 'FontSize', 20 );
set( gcf, 'position', [139  197  1125   490] );
pause(1)
end

%%
pdf( [path 'NS/PDFs/EMT_Exmp_CorrPair_' num2str(i)], [15 6], 1 );



%% Select proteins from both tails of the distribution rr1
% in set Tails

Tails = rr1 < -0.1 | rr1 > 0.56;

% Compute permutaion Clu based in the avarage correlations
Clu = clu( r1(Tails,Tails)+r3(Tails,Tails), 3 ); rb(1)

%% Visualize the clustered correlations
% from set Tails

% Make correlation matrices having just Tails proteins
r1_in = r1(Tails,Tails);
r2_in = r2(Tails,Tails);
r3_in = r3(Tails,Tails);

% Plot the correlations for all proteins ordred by permutaion Clu
imagesc( [r1_in(Clu,Clu) zeros(sum(Tails),5) ...
          r2_in(Clu,Clu) zeros(sum(Tails),5)...
          r3_in(Clu,Clu)]  );
set(gca, 'position', [0.07  0.12  0.9  0.75] );
set(gca, 'FontSize', 16); %, 'Xtick', [], 'Ytick', [] );
h = rb(1);
set(h, 'location', 'NorthOutside', 'FontSize', 20, 'Xtick', -1:0.5:1 );


ht(1) = text(20,-3, 'Day 0' );
ht(2) = text(90,-3, 'Day 3' );
ht(3) = text(160,-3, 'Day 9' );
set(ht, 'FontSize', 26 );

hh(1) = xlabel( 'Proteins' );
hh(2) = ylabel( 'Proteins' );
set(hh, 'FontSize', 28 );


%%
pdf( [path 'NS/PDFs/EMT_correlations_2'], [12 8], 1 );


%% Plot change in correlations & protein abundances
close all
Nm = sum(Tails);
tails = find(Tails);
tails = tails(Clu);
Ones = ones(1,4);

% Plot the correlations for all proteins ordred by permutaion Clu
imagesc( [r1_in(Clu,Clu)  zeros(Nm,4)  mean(dat(tails, ins == 1),2)*Ones zeros(Nm,8) ...
          r2_in(Clu,Clu)  zeros(Nm,4)  mean(dat(tails, ins == 2),2)*Ones zeros(Nm,8) ...
          r3_in(Clu,Clu)  zeros(Nm,4)  mean(dat(tails, ins == 3),2)*Ones       ]     );
set(gca, 'FontSize', 12 );
rb(1);

%%
pdf( [path 'NS/PDFs/EMT_correlations_Abundance'], [12 7], 1 );







%% Plot change in protein abundance
clear hh

close all
Nm = sum(Tails);
tails = find(Tails);
tails = tails(Clu);
%lowTail = rr1 < -0.15; tails
d1 = dat(tails, ins == 1 );  c1 = clu(d1', 3, 1);
d2 = dat(tails, ins == 2 );  c2 = clu(d2', 3, 1);
d3 = dat(tails, ins == 3 );  c3 = clu(d3', 3, 1);

  imagesc ( [d1(Clu,c1)  zeros(Nm,10) ...
             d2(Clu,c2)  zeros(Nm,10) ...
             d3(Clu,c3) ]  );
set(gca, 'FontSize', 16, 'Xtick', [], 'Ytick', [] );

h = rb(1); colormap( yellowblue(50) );
set(h, 'location', 'NorthOutside', 'FontSize', 20, 'Xtick', -1:0.5:1 );
set(gca, 'position', [0.07  0.12  0.9  0.75] );
hh(1) = xlabel( 'Single Cells' );
set(hh, 'FontSize', 28 );

ht(1) = text(40,-3, 'Day 0' );
ht(2) = text(200,-3, 'Day 3' );
ht(3) = text(350,-3, 'Day 9' );
set(ht, 'FontSize', 26 );
%%
pdf( [path 'NS/PDFs/EMT_Abundance'], [12 7], 1 );



%%
[vals,inds] = sort( rr1 );

fprf( ProtIDs(inds), 3, 'EMT_prots.txt' )

%%
%Annot = getFunH( ProtIDs(inds(1:20)) );
% load ( [path 'bin/Human_ONLYsp_Swiss_2_GeneName'] );
% tails(50:75)
[c,ai,bi] = intersect( ProtIDs(in),  HumanONLYspSwiss2GeneName(:, 1 ) );

Annot = HumanONLYspSwiss2GeneName( bi, : );

%%
fprf( Annot(:,3), 1 )
