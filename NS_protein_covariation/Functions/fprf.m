function fprf( Names,  TYPE, File_Name  )


if ~exist( 'TYPE', 'var' ), TYPE=0; end
if nargin  <3, File_Name = 'dump'; end



if nargin ==2 && TYPE == 1
    fprintf( '%s\n',  Names{:} );    
end

if TYPE ==2
    for i=1:numel(Names)
       
        annot = getFunH( Names{i}  );
        if isempty( annot{1} )
            annot = getFun( Names{i}, 2 );
        end
        fprintf( '%s--%s\n',  annot{2}, annot{3}  );  
    end
end


if nargin ==3
    
    fid = fopen(  File_Name, 'w' );

        fprintf( fid, '%s\n',  Names{:} );

    fclose( fid ); 
    
end