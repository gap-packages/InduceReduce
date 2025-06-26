#############################################################################
##  
##  PackageInfo.g for the package `InduceReduce'              Jonathan Gruber
##                                                            
##  This file contains meta-information on the package. It is used by
##  the package loading mechanism and the upgrade mechanism for the
##  redistribution of the package via the GAP website.
##

SetPackageInfo( rec(

PackageName := "InduceReduce",
Subtitle := "Unger's algorithm to compute character tables of finite groups",
Version := "1.0",
Date := "04/09/2018", # dd/mm/yyyy format

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Jonathan",
    LastName := "Gruber",
    Email := "gruber@mathemaik.uni-kl.de",
    PostalAddress := Concatenation( [
                       "Fachbereich Mathematik\n",
                       "Technische Universität Kaiserslautern\n",
                       "Erwin Schrödinger Straße\n",
						"67663 Kaiserslautern\n",
                       "Germany" ] ),
    Place := "Kaiserslautern",
    Institution := "TU Kaiserslautern",
  ),
],

#SourceRepository := rec( Type := "TODO", URL := "URL" ),
#IssueTrackerURL := "TODO",
#SupportEmail := "TODO",

PackageWWWHome :="https://jonathan-144.github.io/InduceReduce",

PackageInfoURL := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "/README.md" ),
ArchiveURL     := Concatenation( ~.PackageWWWHome,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "This package provides an implementation of Unger's algorithm\
 to compute the character table of a finite group.",

PackageDoc := rec(
  BookName  := "InduceReduce",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Computing Character Tables using Unger's algorithm",
),

Dependencies := rec(
  GAP := ">= 4.0",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.5" ] ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

Keywords := [ "character table", "elementary subgroups", "induced characters" ],

AutoDoc := rec(
    entities := rec(
        VERSION := ~.Version,
        RELEASEYEAR := ~.Date{[7..10]},
        RELEASEDATE := function(date)
          local day, month, year, allMonths;
          day := Int(date{[1,2]});
          month := Int(date{[4,5]});
          year := Int(date{[7..10]});
          allMonths := [ "January", "February", "March", "April", "May", "June", "July",
                         "August", "September", "October", "November", "December"];
          return Concatenation(String(day)," ", allMonths[month], " ", String(year));
        end(~.Date),
    ),
),

));


