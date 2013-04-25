function write_perseus_sparse( fname, output )
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
%      fname -- name of text file storing N x M matrix in row/column
%      format. File can use any common delimiter, dlmread will
%      determine the appropriate one from the context.
%      
%      output -- name of file to output for perseus input
% 
%      p is a return code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% someday we should change this to an input argument
num_dim = 2;
dim = num2str( num_dim );

M = dlmread( fname );
mat_size = size( M );

fh = fopen( output, 'w' );

% write dimension to top line
fprintf( fh, dim );
fprintf( fh, '\n' );

for i = 1:mat_size( 1 )
    for j = 1:mat_size( 2 )
        fprintf( fh, '%d %d %d\n', i, j, M( i, j ) );
    end
end       
        
fclose( fh );
