function R = pearNaN ( m )





M = size(m,2);
R = zeros( M, M ); 

            
for i=1:M                   
    for j=i:M 
     
        R(i,j) = pearv( m(:,i), m(:,j)  );
        R(j,i) = R(i,j);
    
    end
end

 
 


function r = pearv( m, n  )

Inds_keep = ~isnan(m) & ~isnan(n); 

if sum(Inds_keep) > 2
    m = m(Inds_keep);
    n = n(Inds_keep); 

    r = ( norv( m(:) )' * norv( n(:) )  ) / (numel(n)-1);
else
    r = NaN;
end
    

% if numel(m) ~= numel(n), 
%     r = []; 
%     fprintf( 'Different Size Vectosrs !\n'  );
% end



 
 
    

function x = norv ( x )
      
   x = x -  mean( x );                   
   x = x * (1/std( x ) ) ; 
   
   
   
   