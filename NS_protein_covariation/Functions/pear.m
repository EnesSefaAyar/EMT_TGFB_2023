function [r, m] = pear ( m, Filter_Tresh )



Inds_2_remove = sum( isnan(m), 2 ) > 0  |  sum(abs(m) == Inf,2) > 0; 
if sum(Inds_2_remove)>0
    m = m( ~Inds_2_remove ,:);   
    fprintf( 'Rows removed: %d \n', sum(Inds_2_remove)  );
end




            [row, clm] = size( m );  
            
for i=1:clm,                     m(:,i) = m(:,i) -  mean( m(:,i) );               
                                 m(:,i) = m(:,i) /   std( m(:,i) );              
end


r   =   m' * m / (row-1);   
 
    

if nargin >= 2  &&  Filter_Tresh ~=0   
    
    r   = ( abs(r) > Filter_Tresh ) .* r;
end

 
 