function write_perseus_cubtop( M, output )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function write_perseus_sparse( matfile, output )
%
%    Author: Jesse Berwald
%    Created: 22 January 2013
%
%    Convert a 2D matrix (digital image) stored in a text file to
%    format readable by Perseus. See
%    http://www.math.rutgers.edu/~vidit/perseus.html.
%     
%      M -- n x m 2D matrix
%      
%      output -- name of file to output for perseus input
% 
%      p is a return code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% someday we should change this to an input argument
num_dim = 2;
dim = num2str( num_dim );

mat_size = size( M );

fh = fopen( output, 'w' );

% write dimension to top line
fprintf( fh, dim );
fprintf( fh, '\n' );

% non-sparse cubical format requires matrix dimensions for
% allocating memory 
fprintf( fh, num2str( mat_size( 1 ) ) );
fprintf( fh, '\n' );
fprintf( fh, num2str( mat_size( 2 ) ) );
fprintf( fh, '\n' );

% now we unravel the matrix into one column, row 1, followed by row
% 2, etc.
% fix the row
for i = 1:mat_size( 1 )
    % traverse the column
    for j = 1:mat_size( 2 )
        fprintf( fh, '%d\n', M( i, j ) );
    end
end       
        
fclose( fh );
