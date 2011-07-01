
function M = flip_strand_a( M, strand );

if( ischar(strand) )
  switch( strand )
   case '+'
    % nothing to do
   case '-'
    M = -M( end:-1:1, end:-1:1 );
   otherwise
    error( 'illegal strand' );
  end;
else
  assert( isnumeric(strand) );
  switch( strand )
   case '+'
    % nothing to do
   case '-'
    M = -M( end:-1:1, end:-1:1 );
   otherwise
    error( 'illegal strand' );
  end;
end;


