
DeclareAttribute("NaturalCharacters", IsGroup);

DeclareGlobalFunction ("pPrimeDecompositionPower");

DeclareAttribute( "pPrimeDecompositions", IsCharacterTable, "mutable" );
DeclareOperation( "pPrimeDecomposition", [ IsCharacterTable, IsPosInt ]  );
DeclareOperation( "pPrimeDecomposition", [ IsGroup, IsPosInt ]  );


DeclareGlobalFunction ("pPrimeCharacter");
DeclareGlobalFunction ("pPrimeRestriction");
