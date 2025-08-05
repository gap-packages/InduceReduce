##  (C) 2025 Frank Lübeck, Lehrstuhl für Algebra und Zahlentheorie, RWTH Aachen
##  
##  This file provides function for computing certain (generalized)
##  characters of a finite group G.
##  
##  NaturalCharacters(G) depending on the representation of the group:
##      permutation groups: permutation characters on orbits of points minus
##              trivial
##      matrix groups in charactistic zero: traces of elements
##      matrix groups over finite fields: permutation characters 
##              on non-zero vectors over field and quadratic extension and
##              on one dimensional subspaces
##      otherwise: regular character
##      
##  Using power maps we can construct from known generalized characters chi:
##      pPrimeCharacter(chi, p)  with values chi(x_{p'}) for primes p
##      pPrimeRestriction(chi, p) with values 0 in p-singular classes
##                               |G|_p chi(x) for x p-regular
##  

# fallback
InstallMethod(NaturalCharacters, ["IsGroup"], function(G)
  local t, ch;
  # regular character
  t := CharacterTable(G);
  ch := 0*[1..NrConjugacyClasses(t)];
  ch[1] := Size(G);
  return [Character(t, ch)];
end);

InstallMethod(NaturalCharacters, ["IsPermGroup"], function(G)
  local t, cls, reps, orb, res, ch, o;
  t := CharacterTable(G);
  cls := ConjugacyClasses(G){IdentificationOfConjugacyClasses(t)};
  reps := List(cls, Representative);
  orb := Orbits(G, MovedPoints(G));
  res := [];
  for o in orb do
    # permutation character on o minus 1
    ch := List(reps, x-> Number(o, i-> i^x = i)) - 1;
    Add(res, Character(t, ch));
  od;
  return res;
end);
InstallMethod(NaturalCharacters, ["IsMatrixGroup"], function(G)
  local t, cls, reps, res, F, q, n, rk, ch, x;
  t := CharacterTable(G);
  cls := ConjugacyClasses(G){IdentificationOfConjugacyClasses(t)};
  reps := List(cls, Representative);
  res := [];
  if Characteristic(G) = 0 then
    Add(res, Character(t, List(reps, TraceMat)));
  elif IsFinite(FieldOfMatrixGroup(G)) then
    # permutation character on non-zero vectors
    # over F_q and over F_(q^2)
    F := FieldOfMatrixGroup(G);
    q := Size(F);
    n := DimensionOfMatrixGroup(G);
    rk := List(reps, x-> RankMat(x - One(G)));
    Add(res, Character(t, List(rk, k-> q^(n-k))-1));
    Add(res, Character(t, List(rk, k-> (q^2)^(n-k))-1));
    # permutation character of one dimensional subspaces of F_q^n minus 1
    ch := [];
    for x in reps do
      Add(ch, Sum(Eigenvalues(F, x), e-> q^(n-RankMat(x - e*One(G)))-1));
    od;
    Add(res, Character(t, ch/(q-1)-1));
  fi;
  return res;
end);

# Let p be a prime.
# An element x of a finite group has a unique decomposition
# x = y z = z y with y being of p'-order and z being a p-element.
# Both, y and z are powers of x, this function returns k with
# y = x^k and z = x/y. This only depends on the order o of x.

InstallGlobalFunction(pPrimeDecompositionPower, function(o, p)
  local q, m, g, res;
  if o mod p <> 0 then
    return 1;
  fi;
  q := p;
  m := o/p;
  while m mod p = 0 do
    q := q*p;
    m := m/p;
  od;
  if m = 1 then
    return 0;
  fi;
  g := GcdRepresentation(q, m);
  res := q*g[1];
  if 2*res > o then
    res := res - o;
  fi;
  return res;
end);

# missing method for groups
InstallOtherMethod(OrdersClassRepresentatives, ["IsGroup"], 
  G-> List(ConjugacyClasses(G), c-> Order(Representative(c))));

# empty on first call
InstallMethod( pPrimeDecompositions, [ IsCharacterTable ],
function(tbl)
    return rec ();
end);
# delegate to character table
InstallMethod(pPrimeDecomposition, [IsGroup, IsPosInt],
function(G, p)
    return pPrimeDecomposition( CharacterTable(G), p );
end);
# 
InstallMethod(pPrimeDecomposition, 
              [IsCharacterTable and HasUnderlyingGroup, IsPosInt],
function(tbl, p)
  local G, map, pms, res, ords, k, i;
  G := UnderlyingGroup(tbl);
  map := pPrimeDecompositions(tbl);
  if not IsPrimeInt(p) then
    Error("pPrimeDecomposition: second argument must be prime.");
  fi;
  if not IsBound(map.(p)) then
      pms := PowerMapsOfAllClasses(G);
      res := [];
      ords := OrdersClassRepresentatives(tbl);
      for i in [1..Length(ords)] do
          k := pPrimeDecompositionPower(ords[i], p);
          Add(res, pms[i][(k mod ords[i])+1]);
      od;
      map.(p) := res;
  fi;
  return map.(p);
end);

# for a generalized chi and prime p the class function f(g) with values 
# chi(h) where h is the p-prime part of g is again a generalized character
# [Isaacs, problem (8.3)]

InstallGlobalFunction(pPrimeCharacter,
function(chi, p)
    local tbl;
    tbl := UnderlyingCharacterTable(chi);
    return Character(tbl, chi{pPrimeDecomposition(tbl, p)});
end);

# for a generalized character chi and prime p the restriction to p'-elements
# (otherwise 0) multiplied by the p-part of the group order
# is again a generalized character [Isaacs, problem (15.3)]

InstallGlobalFunction(pPrimeRestriction,
function(chi, p)
  local tbl, s, q, ch, ords, i;
  tbl := UnderlyingCharacterTable(chi);
  s := Size(tbl);
  q := 1;
  while s mod p = 0 do
    s := s/p;
    q := q*p;
  od;
  ch := [];
  ords := OrdersClassRepresentatives(tbl);
  for i in [1..Length(ords)] do
    if ords[i] mod p = 0 then
      Add(ch, 0);
    else
      Add(ch, q*chi[i]);
    fi;
  od;
  # return original character if is was 0 on p-regular classes
  if q*chi = ch then
    ch := chi;
  fi;
  return Character(tbl, ch);
end);

