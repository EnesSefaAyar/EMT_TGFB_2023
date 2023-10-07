function r = pear2i( m, n, Type  )



% Inds_keep = sum(isnan(m),2)==0  &  sum(isnan(n),2)==0 &...
%             sum(abs(m) == Inf,2)==0  &  sum(abs(n) == Inf,2)==0; 
% if sum(Inds_keep==0)>0
%     m = m(Inds_keep,:);
%     n = n(Inds_keep,:);    
%     fprintf( 'Size: %d x %d\n', size( n ) );
% end


if nargin <= 2
    Type = 'Pearson'; 
else
    Type = 'Spearman';
    [~, m] = sort( m );   
    [~, n] = sort( n );   
end


[row, clm] = size( m );  

%l   =   1 / (row-1);  

r = zeros( clm, 1 ); 
            
for i = 1:clm
    
       Inds = ~isnan( m( :, i ) ) & ~isnan( n( :, i ) ) & abs( m( :, i ) ) < 1e100 & abs( n( :, i) ) < 1e100; 
         
       r(i) = norv( m( Inds, i) )' *  norv( n( Inds, i) ) / ( sum( Inds ) - 1 );         
end

end
    

function n = norv ( n )
      
   n = n - mean( n );                   
   n = n * ( 1 / std ( n ) ) ; 

end
 