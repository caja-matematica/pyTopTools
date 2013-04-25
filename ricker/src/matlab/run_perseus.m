function run_perseus( fname, output, pers_type )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function run_perseus
%
%    Author: Jesse Berwald
%    Created: 22 January 2013
%
%    Call perseus (must be in path) on fname. Assumes that fname is in
%    'sparse cubical toplex' format (ie. 'scubtop'). See
%    http://www.math.rutgers.edu/~vidit/perseus.html.
%
%      fname -- full path to file
%
%      output -- prefix for perseus output. Will be appended with
%      output_*.txt, where * denotes homological dimension.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( pers_type, 'sparse' )
    command = sprintf( 'perseus scubtop %s %s > p.log', fname, output );
else
    command = sprintf( 'perseus cubtop %s %s > p.log', fname, output );
end
   
unix( command );
