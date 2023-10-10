function pdf( File_Name, Paper_Size, OPEN, Flag )



%set( gca, 'Position', [0.13 0.10 0.78 0.84] );
if nargin > 1
   
    if isempty( Paper_Size )    
        load( 'C:\Program Files\MATLAB\PaperSize.mat' );
        Paper_Size = PaperSize;
    end
    
   set( gcf, 'PaperSize', Paper_Size,...
             'PaperPositionMode', 'auto' );
end



if nargin >= 4 && Flag >= 1
    set( gcf, 'Position', [440  359  692   419] )
    set( gca, 'Position', [0.14 0.14 0.80  0.78] );
end 
    
File_Name = [ File_Name  '.pdf' ];
%print( '-dpdf',  File_Name ); %'-cmyk' '-bestfit', '-fillpage',
print( '-dpdf', '-fillpage',   File_Name ); %'-cmyk' '-bestfit', '-fillpage',


if nargin >= 3 && OPEN >= 1
    %system( sprintf( 'start  "%s"',  File_Name ) );
    PWD = pwd;
    if strmatch( 'C', PWD(1)  )
        system( ['start  C:'  File_Name ] );
    else
        system( ['open  '  File_Name ] ); 
    end
    
    %open(  File_Name  );
    %system( ['open  -a /app/adobe  ' File_Name ] );  
    %system( File_Name );    "C:\Program Files (x86)\Adobe\Reader 9.0\Reader AcroRd32.exe"
end