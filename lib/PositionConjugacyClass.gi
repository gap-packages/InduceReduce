##  (C) 2025 Frank Lübeck, Lehrstuhl für Algebra und Zahlentheorie, RWTH Aachen
##  
##  This file contains functions to find the number of the conjugacy class
##  of an element in a group G. 
##  
##  ConjugacyClassInvariants(G) produces a record r with a recursive structure 
##      r.tree
##  This is a list with 4 entries:
##      poss: a list of class positions
##      f: a function with arguments r and a group element
##         that returns a class invariant that can distinguish classes poss
##      res: the sorted list of results of f for representatives of classes poss
##      ts: a list of the same length as res, if res[i] is an invariant of
##          a unique class j in poss then ts[i] = j; otherwise ts[i] is again
##          a tree for the classes with invariant res[i].
##  
##  
##  This structure is used by PositionConjugacyClass(G, x) to identify the
##  number of the conjugacy class of x by computing class invariants until the
##  number is found.
##
##  An optional second argument of ConjugacyClassInvariants can be a list 
##  of functions of form g(r, tree) ##  that extend a tree that so far 
##  only contains the list poss.  There are default functions for permutation
##  and matrix groups.
##  
##  Note that a similar scheme could also be used for identifying other
##  equivalence classes.
##  

# refine the tree with function fu
CCInvFuncs.refine := function(r, tree, fu)
  local invs, s, ts, pos, a;
  invs := List(r.reps{tree[1]}, x-> fu(r, x));
  s := Set(invs);
  if Length(s) > 1 then
    Add(tree, fu);
    Add(tree, s);
    ts := [];
    for a in s do
      pos := Positions(invs, a);
      if Length(pos) = 1 then
        Add(ts, tree[1][pos[1]]);
      else
        Add(ts, [tree[1]{pos}]);
      fi;
    od;
    Add(tree, ts);
  fi;
end;

# In many practical cases it is likely that conjugacy tests are
# made with class representatives. Here is a cheap test to recognize such
# cases.
CCInvFuncs.NoticeReps := function(r, tree)
  local s, len, l, fu;
  if tree[1] = [1..Length(r.classes)] then
    s := ShallowCopy(r.reps);
    len := Length(s);
    l := [1..len];
    SortParallel(s, l);
    fu := function(r, x)
      local pos;
      pos := PositionSorted(s, x);
      if not IsBound(s[pos]) or s[pos] <> x then
        return fail;
      else
        return l[pos];
      fi;
    end;
    tree[2] := fu;
    tree[3] := [1..len];
    tree[4] := [1..len];
    Add(tree[3], fail);
    Add(tree[4], [[1..len]]);
  fi;
end;

# fallback: try IsConjugate with all classes  but the last
CCInvFuncs.ConjTest := function(r, tree)
  local fu;
  fu := function(r, x)
    local l, i;
    l := tree[1];
    for i in [1..Length(l)-1] do
      if x = r.reps[l[i]] or x in r.classes[l[i]] then
        return l[i];
      fi;
    od;
    return l[Length(l)];
  end;
  Add(tree, fu);
  Add(tree, Set(tree[1]));
  Add(tree, tree[3]);
end;

# use characteristic polynomial in matrix groups
CCInvFuncs.CharPol := function(r, tree)
  local fu;
  fu := function(r, x)
    return CharacteristicPolynomial(x);
  end;
  CCInvFuncs.refine(r, tree, fu);
end;

# use minimal polynomial in matrix groups
CCInvFuncs.MinPol := function(r, tree)
  local fu;
  fu := function(r, x)
    return MinimalPolynomial(x);
  end;
  CCInvFuncs.refine(r, tree, fu);
end;


# for permutation groups: cycle structures on orbits of points
CCInvFuncs.CycStruct := function(r, tree)
  local o, fu;
  o := Orbits(r.G, MovedPoints(r.G));
  if Length(o) = 1 then
    fu := function(r, x)
      return CycleStructurePerm(x);
    end;
  else
    fu := function(r, x)
      return List(o, a-> CycleStructurePerm(RestrictedPerm(x,a)));
    end;
  fi;
  CCInvFuncs.refine(r, tree, fu);
end;

BindGlobal("ConjugacyClassInvariants", function(G, args...)
  local r, funcs, tree, find;
  if IsBound(G!.ConjugacyClassInvariants) then
    return G!.ConjugacyClassInvariants;
  fi;
  r := rec(G := G);
  r.classes := ConjugacyClasses(G);
  r.reps := List(r.classes, Representative);
  r.tree := [[1..Length(r.reps)]];
  if Length(args) > 0 then
    funcs := ShallowCopy(args);
    Add(funcs, CCInvFuncs.ConjTest);
  elif IsPermGroup(G) then
    funcs := [CCInvFuncs.NoticeReps, CCInvFuncs.CycStruct, CCInvFuncs.ConjTest];
  elif IsMatrixGroup(G) then
    funcs := [CCInvFuncs.NoticeReps, CCInvFuncs.CharPol, CCInvFuncs.MinPol, CCInvFuncs.ConjTest];
  else
    funcs := [CCInvFuncs.NoticeReps, CCInvFuncs.ConjTest];
  fi;
  tree := r.tree;
  find := function(tree, funcs)
    local f2, a;
    funcs[1](r, tree);
    f2 := funcs{[2..Length(funcs)]};
    if Length(tree) = 1 then
      find(tree, f2);
    else
      for a in tree[4] do
        if IsList(a) then
          find(a, f2);
        fi;
      od;
    fi;
  end;
  find(tree, funcs);
  G!.ConjugacyClassInvariants := r;
  return r;
end);

# applying this is easy:
BindGlobal("PositionConjugacyClass", function(G, x)
  local r, tree, inv, pos;
  r := ConjugacyClassInvariants(G);
  tree := r.tree;
  while IsList(tree) do
    inv := tree[2](r, x);
    pos := PositionSorted(tree[3], inv);
    tree := tree[4][pos];
  od;
  return tree;
end);

