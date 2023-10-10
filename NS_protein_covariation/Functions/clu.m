 function [inds, r, Clust_Tree_n, T] = clu ( x, TYPE, NoPlot )


if size( x, 1 ) > 3e4 && TYPE > 1
   warning( 'More than 30,000 variables. I give up :(' ); %#ok<*WNTAG>
   return
end

if nargin >= 2 && TYPE<1
     Conditions = (1:size(x,2))';
     r = corr( x', Conditions, 'type', 'Pearson' );%Spearman
     [val, ind] = sort( r );
     inds = ind( abs(val)>TYPE );
     im( nor_mean( x(inds,:), 2 ) ); fgf,  rg(1);
     %colormap( redblue );
     set( gca, 'FontSize', 16, 'FontWeight', 'Bold' );
     set( gca, 'Xtick', Conditions ) 
   return
end

if nargin < 2, TYPE =  'corr'; end; 

switch TYPE
    case { 'corr', 1 }, Type = 'Correlation';              TYPE = 'corr';
    case { 'euc',  2 }, Type = 'Euclidean Distance';       TYPE = 'euc';
    case { 'cos',  3 }, Type = 'Unnormalized Correlation'; TYPE = 'cosine';
    case { 'sort', 4 }, Type = 'Sorting';                  TYPE = 'SORT';
end;

fprintf( 'I will use %s as a similarity measure ... \n' ,  Type );
%%

if ~strcmp( TYPE, 'SORT' )
    
    dist_n = pdist( x, TYPE );  
    Clust_Tree_n = linkage( dist_n, 'complete' ); %average
    
    
    if nargout >= 4
         T = cluster( Clust_Tree_n ,'maxclust', 3 ); 
    end

    %fig
    [~, ~, inds ] = dendrogram (Clust_Tree_n, 0); 
    clear H T; % Clust_Tree_n;
    %close all;

else 
    %[~, inds] = sort( sum( x, 2) );
    [~, inds] = sort( x(:,1) );
end

if nargin>2 && NoPlot>0, return, end
    
%[row clm] = size( x );

if isequal( x, x' ), 
    r = x( inds, inds );
else
    r = x( inds,  : );
end;
%%

imagesc( r );               


% colormap( redgreencmap ); colorbar 
% 
% Lim = input( 'Please, Pass Me a Limit (the default is [-1 1]):  ' );
% 
% if isempty(Lim), Lim = 1; end 
% 
% 
% set( gca, 'CLim', [-Lim Lim] );













