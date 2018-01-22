% ConvertWdfToHDF5.m
%
% Matlab script because of available Matlab libraries for Yokogawa WDF file
% format.
%
% Libraries downloaded from Yokogawa website. Trial version but I do not
% know what the limitation is.
%
%
% Based on:
%  MATLAB M-File : WdfSample4.m
%  WDF file access sample script- 4
%  Copyright (C) Yokogawa Meters & Instruments Corporation
%  Software Japan. All rights reserved.
%
% 
%
% 2017-10-05 QuTech/Hansonlab, Wouter Westerveld
%

% Input dataDir
dataDir = 'D:\measuring\data\20171012\scope\' 
 
% 
files = dir( strcat( dataDir, '*.wdf' ) );

for i = 1:length(files)
    filename = strcat( dataDir, files(i).name );
    filenameOut = strcat( filename, '.h5'  );    
    if exist( filenameOut, 'file' )
        disp( strcat( filename, '... skipping.' ) );
        continue
    else
        disp( strcat( filename, '...' ) );
    end
        
    %
    %   Input the Filename
    %
    [ ret, chNum ] = mexWdfGetChNum( filename );
    block = 0;

    for ch = 0 : double( chNum ) - 1       
        %   Get file parameters                
        [ ret, traceName ] = mexWdfItemRead ( filename, 'TraceName', ch, 0 );     
        [ ret, vResolution ] = mexWdfItemRead ( filename, 'VResolution', ch, block );
        [ ret, vOffset ] = mexWdfItemRead ( filename, 'VOffset', ch, block );  
        [ ret, vScaleUpper ] = mexWdfItemRead ( filename, 'VScaleUpper', ch, block );
        [ ret, vScaleLower ] = mexWdfItemRead ( filename, 'VScaleLower', ch, block );                                        
        [ ret, hResolution ] = mexWdfItemRead ( filename, 'HResolution', ch, block );
        [ ret, hOffset ] = mexWdfItemRead ( filename, 'HOffset', ch, block );  
        [ ret, vUnit ] = mexWdfItemRead ( filename, 'VUnit', ch, block );
        [ ret, hUnit ] = mexWdfItemRead ( filename, 'HUnit', ch, block );                  
                        
        %   Get waveform data        
        clear data;
        [ ret, param, data ] = mexWdfDataRead( filename, ch, block );             
        h5create( filenameOut, sprintf('/CH%.0f', ch), size( data ), 'Datatype', 'int16' );
        % Store data
        h5write( filenameOut, sprintf('/CH%.0f', ch), data );             
        % Store attributes
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'TraceName', traceName );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'VResolution', vResolution );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'VOffset', vOffset );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'VScaleUpper', vScaleUpper );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'VScaleLower', vScaleLower );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'HResolution', hResolution );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'HOffset', hOffset );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'VUnit', vUnit );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'HUnit', hUnit );
        h5writeatt( filenameOut, sprintf('/CH%.0f', ch), 'Note', 'Data in V is: data * VResolution + VOffset' );
                                
    end   
    
end

disp( 'finished :-) ' )
