function process_perseus_sub_super( N, tmpPerseus, scale, max_dim, ...
                                    pers_type )

 if nargin < 3
     error( 'process_perseus takes at least 3 arguments.' )
 end
 if nargin < 4
     max_dim = 1;
 else
     max_dim = max_dim;
 end
 if nargin < 5
     p = 'cubtop';
 else
     p = pers_type;
 end
        
 % convert population matrix to integers
 M = scale * N;
 M = ceil( M )
 
 % convert M to super-level format
 max_val = max( max( M ) );
 Msup = -M + max_val + 1;  % conversion function
 max_msup = max( max( Msup ) );
 min_sup = min( min( Msup ) );
% $$$  
% $$$  figure;
% $$$  bar3( M );
% $$$  figure;
% $$$  bar3( Msup );
 

 % Convert mat_text to Perseus file format. 
 perseus_in = strcat( tmpPerseus, '_sub.txt' );
 perseus_in_sup = strcat( tmpPerseus, '_super.txt' );

 
 % assume if not sparse then full cubtop complex
 if strcmp( p, 'sparse' )
     write_perseus_sparse_mat( M, perseus_in );
     write_perseus_sparse_mat( Msup, perseus_in_sup );
 else
     disp( 'Writing Perseus file!' )
     write_perseus_cubtop( M, perseus_in );
     write_perseus_cubtop( Msup, perseus_in_sup );
 end
 
 % call assumes that 'perseus' is in your path. Perseus will
 % append '_*.txt' to perseus_out
 disp( 'Running Perseus on:' );
 disp( perseus_in );
 disp( perseus_in_sup );
 %
 perseus_out = strcat( tmpPerseus, '_sub' );
 perseus_out_sup = strcat( tmpPerseus, '_super' );

 % run perseus on sub- and super-level versions of pop matrix
 run_perseus( perseus_in, perseus_out, p );
 run_perseus( perseus_in_sup, perseus_out_sup, p );
 
 % Now combine the output files into one file for diagram plotting
 for d = 0:max_dim
     
     fprintf( 'DIMENSION %d\n', d );
     
     tmp1 = strcat( perseus_out, '_', num2str( d ), '.txt' );
     tmp2 = strcat( perseus_out_sup, '_', num2str( d ), '.txt' );
     try
         P = load( tmp1 );
         % for dim > 0, sometimes we get [] matrices
         if size( P ) == [0,0]
             got_sub = 0;
         else
             got_sub = 1;
         end
     catch
         got_sub = 0;
         disp( strcat( 'No perseus file found for',{' '}, tmp1 ) );
     end
     try
         Psup = load( tmp2 );
         if size( Psup ) == [0,0]
             got_super = 0;
         else
             got_super = 1;          
         end
     catch
         got_super = 0;
         disp( strcat( 'No perseus file found for', {' '}, tmp2 ) );
     end
   
     % convert super-level Perseus data back using max_val and account for
     % -1's
     if ( got_super == 1 )
         % Shift births by max_val since we are considering
         % super-level sets.
         inf_idx = find( Psup(:,2) == -1 );
         non_inf = find( Psup(:,2) ~= -1 );
         Psup = -Psup + max_val + 1;
         
         % designate superlevel inf's differently and undo scaling
         Psup( inf_idx, 2 ) = -2;
         Psup( non_inf, : ) = Psup( non_inf, : ) / scale;
         Psup( :, 1 ) = Psup( :, 1 ) / scale; % this undoes scale
                                              % twice in frist col?!
         if ( got_sub == 1 ) 
             non_inf = find( P(:,2) ~= -1 );
             P( non_inf, : ) = P( non_inf, : ) / scale;
             P( :, 1 ) = P( :, 1 ) / scale;
             full_pers = [ P; Psup ];    
         else
             % already rescaled above, so no need here
             full_pers = Psup;
         end
     % We just have sublevels 
     else
         non_inf = find( P(:,2) ~= -1 );
         P( non_inf, : ) = P( non_inf, : ) / scale;
         P( :, 1 ) = P( :, 1 ) / scale;
         full_pers = P;
     end
     
     full_pers_fname = strcat( tmpPerseus, '_full_', num2str( d ), ...
                               '.pdia' );
     dlmwrite( full_pers_fname, full_pers, 'delimiter', ...
               ' ');
 end
 % end for
 
 
