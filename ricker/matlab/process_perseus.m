function process_perseus( N, tmpPerseus, scale, pers_type )

 if nargin < 3
     error( 'process_perseus takes at least 3 arguments.' )
 end
 if nargin < 4
     p = 'cubtop';
 else
     p = pers_type;
 end
 
        
 % convert population matrix to integers
 M = scale * N;
 M = floor( M );
 
 % write scaled population matrix to disk
 % mat_text = strcat( tmpPerseus, '.txt' );
 % dlmwrite( mat_text, M, 'delimiter', ' ', 'precision', '%d');
 
 % Convert mat_text to Perseus file format. 
 perseus_in = strcat( tmpPerseus, '_pers.txt' );
 
 % assume if not sparse then full cubtop complex
 if strcmp( p, 'sparse' )
     write_perseus_sparse_mat( M, perseus_in );
 else
     write_perseus_cubtop( M, perseus_in );
 end
 % call assumes that 'perseus' is in your path. Perseus will
 % append '_*.txt' to perseus_out
 perseus_out = strcat( tmpPerseus, '_pers' );
 
 run_perseus( perseus_in, perseus_out, p );
 
