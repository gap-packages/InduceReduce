gap> START_TEST("indredtest.tst");

# Neither the order of the classes of a character table nor the order of
# its irreducibles is specified.
# Compare the values at given group elements, as sets.
gap> ValuesAt:= function( irr, elms )
>      local ccl, pos;
>      ccl:= ConjugacyClasses( UnderlyingCharacterTable( irr[1] ) );
>      pos:= List( elms, g -> PositionProperty( ccl, C -> g in C ) );
>      return Set( List( irr, chi -> ValuesOfClassFunction( chi ){ pos } ) );
>    end;;

#############################################################################
# testing the main function
gap> G:= AlternatingGroup( 6 );;
gap> reps:= List( ConjugacyClasses( G ), Representative );;
gap> ct:= CharacterTableUnger( G );
CharacterTable( Alt( [ 1 .. 6 ] ) )
gap> ValuesAt( Irr( ct ), reps ) = ValuesAt( IrrDixonSchneider( G ), reps );
true

#
gap> G:= GeneralLinearGroup( 2, 3 );;
gap> reps:= List( ConjugacyClasses( G ), Representative );;
gap> ct:= CharacterTableUnger( G );
CharacterTable( GL(2,3) )
gap> ValuesAt( Irr( ct ), reps ) = ValuesAt( IrrDixonSchneider( G ), reps );
true
gap> OrdersClassRepresentatives( ct )
>    = List( ConjugacyClasses( ct ), C -> Order( Representative( C ) ) );
true
gap> SizesCentralizers( ct )
>    = List( ConjugacyClasses( ct ), C -> Size( G ) / Size( C ) );
true

# The classes of the table may be ordered differently from those of the
# group: store the identity class last in the group, the table has it first.
gap> H:= GroupWithGenerators( GeneratorsOfGroup( G ) );;
gap> SetConjugacyClasses( H,
>        List( Concatenation( reps{ [ 2 .. Length( reps ) ] }, reps{ [ 1 ] } ),
>              g -> ConjugacyClass( H, g ) ) );
gap> ConjugacyClasses( CharacterTable( H ) ) = ConjugacyClasses( H );
false
gap> ValuesAt( IrrUnger( H ), reps ) = ValuesAt( IrrDixonSchneider( G ), reps );
true
gap> ct:= CharacterTableUnger( H );;
gap> ValuesAt( Irr( ct ), reps ) = ValuesAt( IrrDixonSchneider( G ), reps );
true
gap> OrdersClassRepresentatives( ct )
>    = List( ConjugacyClasses( ct ), C -> Order( Representative( C ) ) );
true

# testing the `Irr` method from the package;
# call it directly, `Irr( G )` may select a method from elsewhere
gap> G:= DihedralGroup( 24 );;
gap> ct:= CharacterTable( G );;
gap> ForAny( ApplicableMethod( Irr, [ G, 0 ], 0, "all" ),
>            f -> IsIdenticalObj( f, IndRed.Irr ) );
true
gap> irr:= IndRed.Irr( G, 0 );;
gap> irr = Irr( ct );
true
gap> HasInfoText( ct );
false
gap> HasIsSupersolvableGroup( G );
true
gap> ForAll( PrimeDivisors( Size( ct ) ),
>            p -> IsBound( ComputedPowerMaps( ct )[p] ) );
true
gap> G:= MathieuGroup( 12 );;
gap> ct:= CharacterTable( G );;
gap> irr:= IndRed.Irr( G, 0 );;
gap> irr = Irr( ct );
true
gap> InfoText( ct );
"origin: Unger's algorithm"
gap> ForAll( PrimeDivisors( Size( ct ) ),
>            p -> IsBound( ComputedPowerMaps( ct )[p] ) );
true
gap> reps:= List( ConjugacyClasses( G ), Representative );;
gap> ValuesAt( irr, reps ) = ValuesAt( IrrDixonSchneider( G ), reps );
true

#
gap> STOP_TEST("indredtest.tst");
