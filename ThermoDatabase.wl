(* ::Package:: *)

(* ::Title:: *)
(*  THERMODATABASE*)
(*  Thermochemical Properties Module for Gas-Turbine Combustor Modeling*)


(* ::Text:: *)
(*Purpose*)
(**)
(*  Import the thermochemical data required by the DRM19 kinetic mechanism from  THERM.DAT and evaluate standard-state ideal-gas thermodynamic properties  using NASA 7-coefficient polynomials. *)
(*  *)
(*  Properties provided*)
(**)
(*      ThermoCp[species, T]   Cp^0(T)   [J mol^-1 K^-1]*)
(*      ThermoH [species, T]   H^0(T)    [J mol^-1]*)
(*      ThermoS [species, T]   S^0(T)    [J mol^-1 K^-1]*)
(*      ThermoG [species, T]   G^0(T)    [J mol^-1]*)
(**)
(*  Project*)
(*  -------*)
(*  Gas-Turbine Combustor Model*)
(*  One-Dimensional Reactive-Flow and Chemical-Kinetics Framework*)
(**)
(*  - Author*)
(*  Hiago David Zogbi Silva Oliveira*)
(*  Engenharia Qu\[IAcute]mica/Mec\[AHat]nica - UFF*)
(*  *)
(*  - Orientador*)
(*  Dr . Isabela Florindo *)
(**)
(*  - Required data*)
(*  THERM.DAT*)
(**)
(*  - Scope*)
(**)
(*  This module intentionally reads ONLY the thermochemical records required by the combustor model (DRM19 species plus the N and S species of the NOx/SOx subsets, listed in $CombustorThermoSpecies). THERM.DAT contains many additional gas, liquid and solid records   some of which use formats that are not relevant to the present gas-phase  combustor model. Restricting the import to the required records avoids  ambiguity and prevents unrelated condensed-phase entries from being parsed as NASA-7 gas records.*)
(**)
(*  NASA 7-coefficient formulation*)
(**)
(*      Cp^0/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4*)
(**)
(*      H^0/(R T) = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T*)
(**)
(*      S^0/R = a1 Log[T] + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7*)
(**)
(*      G^0 = H^0 - T S^0*)
(**)
(*  Internal unit convention*)
(**)
(*      Temperature        K*)
(*      Cp                 J mol^-1 K^-1*)
(*      H                  J mol^-1*)
(*      S                  J mol^-1 K^-1*)
(*      G                  J mol^-1*)
(**)
(*  Main references*)
(**)
(*  [1] McBride, B. J.; Gordon, S.; Reno, M. A. Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species.  NASA Technical Memorandum 4513, 1993.*)
(**)
(*  [2] Gordon, S.; McBride, B. J. Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. NASA Reference Publication 1311, 1994.*)
(**)
(*  [3] Kee, R. J.; Rupley, F. M.; Miller, J. A. CHEMKIN-II: A FORTRAN Chemical  Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics.  Sandia National Laboratories, SAND89-8009B.*)
(**)
(*  Notes*)
(**)
(*  - The names used by DRM19 do not always match the descriptive labels used in  THERM.DAT. Exact 18-character header-field selectors are therefore defined  below. Examples: *)
(*     HCO -> "CHO", CH2 -> "CH2 TRIPLET",*)
(*    CH2(S) -> "CH2 SINGLET", CH2O -> "HCHO Formaldehyde".*)
(*  - No kinetic calculations are performed in this module.*)
(*  - Extrapolation outside each NASA polynomial validity interval is rejected.*)
(**)
(*  Version 2.0*)
(*  September 2026*)


(* ::Section:: *)
(*GLOBAL CONSTANTS AND DATABASE STATE *)


ClearAll[
  $ThermoR,
  $ThermoDatabase,
  $ThermoSourcePath,
  $ThermoSelectors
];

$ThermoR = 8.31446261815324;  (* J mol^-1 K^-1 *)

$ThermoDatabase = <||>;
$ThermoSourcePath = None;

$CombustorThermoSpecies = {
  "H2", "H", "O", "O2", "OH", "H2O", "HO2",
  "CH2", "CH2(S)", "CH3", "CH4", "CO", "CO2",
  "HCO", "CH2O", "CH3O", "C2H4", "C2H5", "C2H6",
  "N2", "AR","H2O2",

  (* Nitrogen species *)
  "N", "NO", "NO2", "N2O",

  (* Sulfur species *)
  "SO", "SO2", "SO3", "S", "SH", "H2S", "COS"
};


(* ::Section:: *)
(*EXACT THERM.DAT RECORD SELECTORS *)


(*
   THERM.DAT uses a fixed-width 18-character species/description field.
   These selectors were chosen explicitly for the species used by DRM19.
   Exact matching prevents CH3 from being confused with CH3O, CH3OH, etc.
*)

$ThermoSelectors = <|
  "H2"     -> "H2 REF ELEMENT    ",
  "H"      -> "H                 ",
  "O"      -> "O                 ",
  "O2"     -> "O2 REF ELEMENT    ",
  "OH"     -> "OH HYDROXYL RADI  ",
  "H2O"    -> "H2O               ",
  "HO2"    -> "HO2               ",
  "CH2"    -> "CH2 TRIPLET       ",
  "CH2(S)" -> "CH2 SINGLET       ",
  "CH3"    -> "CH3               ",
  "CH4"    -> "CH4    RRHO       ",
  "CO"     -> "CO                ",
  "CO2"    -> "CO2               ",
  "HCO"    -> "CHO               ",
  "CH2O"   -> "HCHO Formaldehyde ",
  "CH3O"   -> "CH3O  METHOXY RAD ",
  "C2H4"   -> "C2H4              ",
  "C2H5"   -> "C2H5 ethyl radic  ",
  "C2H6"   -> "C2H6              ",
  "N2"     -> "N2 REF ELEMENT    ",
  "AR"     -> "AR REF ELEMENT    ",
  "H2O2"   -> "H2O2              ",
  (* Nitrogen species *)
  "N"      -> "N                 ",
  "NO"     -> "NO                ",
  "NO2"    -> "NO2               ",
  "N2O"    -> "N2O               ",

  (* Sulfur species *)
  "SO"     -> "SO                ",
  "SO2"    -> "SO2               ",
  "SO3"    -> "SO3               ",
  "S"      -> "S                 ",
  "SH"     -> "SH                ",
  "H2S"    -> "H2S               ",
  "COS"    -> "COS               "
|>;


(* ::Section:: *)
(*MESSAGES*)

ClearAll[LoadThermoDatabase, ThermoRecord, ThermoCp, ThermoH, ThermoS, ThermoG];

LoadThermoDatabase::nofile =
  "THERM.DAT was not found at `1`.";

LoadThermoDatabase::read =
  "THERM.DAT could not be imported from `1`.";

LoadThermoDatabase::nomatch =
  "No exact THERM.DAT record was found for DRM19 species `1` using selector `2`.";

LoadThermoDatabase::multiple =
  "More than one exact THERM.DAT record was found for DRM19 species `1` using selector `2`: lines `3`.";

LoadThermoDatabase::phase =
  "The selected THERM.DAT record for species `1` at line `2` is not a gas-phase (G) record.";

LoadThermoDatabase::marker =
  "Invalid NASA record line markers for species `1` beginning at THERM.DAT line `2`.";

LoadThermoDatabase::number =
  "A numeric field could not be parsed for species `1` beginning at THERM.DAT line `2`.";

LoadThermoDatabase::failed =
  "THERM.DAT loading failed. The thermodynamic database was not modified.";

ThermoRecord::noload =
  "No thermo database is loaded. Call LoadThermoDatabase[path] first.";

ThermoRecord::species =
  "Species `1` is not present in the loaded thermodynamic database.";

ThermoRecord::temp =
  "Temperature `1` K is outside the valid NASA range [`2`, `3`] K for species `4`.";

ThermoRecord::tpos =
  "Temperature must be a positive numeric value; received `1`.";


(* ::Section:: *)
(*LOW-LEVEL FIXED-WIDTH PARSER*)


ClearAll[
  normalizeThermoLines,
  fixedField,
  thermoNumber,
  nasaCoefficientFields,
  nasaLineMarker,
  exactSelectorHits,
  parseThermoCard
];

normalizeThermoLines[text_String] :=
  StringSplit[
    StringReplace[text, {"\r\n" -> "\n", "\r" -> "\n"}],
    "\n"
  ];

(* Return a fixed-width field, padding short lines with spaces if necessary. *)
fixedField[line_String, i_Integer, j_Integer] :=
  StringTake[StringPadRight[line, Max[j, StringLength[line]]], {i, j}];

(* Convert FORTRAN-style E/D scientific notation to a Wolfram real number. *)
thermoNumber[s_String] := Module[{x, cleaned},
  cleaned = StringTrim[s];

  If[cleaned === "", Return[$Failed]];

  cleaned = StringReplace[
    cleaned,
    RegularExpression["[EeDd]([+-]?[0-9]+)"] :> "*^$1"
  ];

  x = Quiet @ Check[ToExpression[cleaned], $Failed];

  If[NumericQ[x], N[x], $Failed]
];

(* NASA coefficient lines contain 15-character numeric fields. *)
nasaCoefficientFields[line_String, n_Integer] :=
  Table[
    thermoNumber @ fixedField[line, 15 (k - 1) + 1, 15 k],
    {k, 1, n}
  ];

nasaLineMarker[line_String] :=
  StringTrim @ fixedField[line, 80, 80];

(* Match ONLY the exact 18-character first field of a THERM.DAT header. *)
exactSelectorHits[lines_List, selector_String] :=
  Select[
    Range[1, Length[lines]],
    fixedField[lines[[#]], 1, 18] === selector &
  ];

parseThermoCard[
  lines_List,
  species_String,
  selector_String
] := Module[
  {
    hits, i, header, phase, tmin, tmax, tmid,
    c2, c3, c4, high, low, markers
  },

  hits = exactSelectorHits[lines, selector];

  Which[
    Length[hits] == 0,
      Message[LoadThermoDatabase::nomatch, species, selector];
      Return[$Failed],

    Length[hits] > 1,
      Message[LoadThermoDatabase::multiple, species, selector, hits];
      Return[$Failed]
  ];

  i = First[hits];

  If[i + 3 > Length[lines],
    Message[LoadThermoDatabase::marker, species, i];
    Return[$Failed]
  ];

  header = lines[[i]];

  (* NASA/CHEMKIN fixed-width header positions. *)
  
  phase = StringTrim @ fixedField[header, 45, 45];
  tmin  = thermoNumber @ fixedField[header, 46, 55];
  tmax  = thermoNumber @ fixedField[header, 56, 65];
  tmid  = thermoNumber @ fixedField[header, 66, 75];

  If[phase =!= "G",
    Message[LoadThermoDatabase::phase, species, i];
    Return[$Failed]
  ];

  markers = nasaLineMarker /@ lines[[i ;; i + 3]];

  If[markers =!= {"1", "2", "3", "4"},
    Message[LoadThermoDatabase::marker, species, i];
    Return[$Failed]
  ];

  c2 = nasaCoefficientFields[lines[[i + 1]], 5];
  c3 = nasaCoefficientFields[lines[[i + 2]], 5];
  c4 = nasaCoefficientFields[lines[[i + 3]], 4];

  If[
    ! And @@ (NumericQ /@ Join[{tmin, tmax, tmid}, c2, c3, c4]),
    Message[LoadThermoDatabase::number, species, i];
    Return[$Failed]
  ];

  (*
     Standard NASA 7 layout:
       line 2 : high a1 ... a5
       line 3 : high a6 a7, low a1 a2 a3
       line 4 : low  a4 ... a7
  *)
  high = Join[c2, Take[c3, 2]];
  low  = Join[Drop[c3, 2], c4];

  If[Length[high] =!= 7 || Length[low] =!= 7,
    Message[LoadThermoDatabase::number, species, i];
    Return[$Failed]
  ];

  <|
    "Species"        -> species,
    "Selector"       -> selector,
    "Header"         -> header,
    "SourceLine"     -> i,
    "Phase"          -> phase,
    "Tmin"           -> tmin,
    "Tmax"           -> tmax,
    "Tmid"           -> tmid,
    "Low"            -> low,
    "High"           -> high
  |>
];


(* ::Section:: *)
(*DATABASE LOADING*)


ClearAll[LoadThermoDatabase];

LoadThermoDatabase[path_String] := Module[
  {text, lines, candidate},

  If[! FileExistsQ[path],
    Message[LoadThermoDatabase::nofile, path];
    Return[$Failed]
  ];

  text = Quiet @ Check[Import[path, "Text"], $Failed];

  If[text === $Failed || ! StringQ[text],
    Message[LoadThermoDatabase::read, path];
    Return[$Failed]
  ];

  lines = normalizeThermoLines[text];

  candidate = Association @ Table[
    sp -> parseThermoCard[lines, sp, $ThermoSelectors[sp]],
    {sp, $CombustorThermoSpecies}
  ];

  If[AnyTrue[Values[candidate], # === $Failed &],
    Message[LoadThermoDatabase::failed];
    Return[$Failed]
  ];

  (* Database is replaced only after ALL records have passed validation. *)
  $ThermoDatabase = candidate;
  $ThermoSourcePath = ExpandFileName[path];

  <|
    "Status"     -> "Loaded",
    "Source"     -> $ThermoSourcePath,
    "Count"      -> Length[$ThermoDatabase],
    "Species"    -> Keys[$ThermoDatabase]
  |>
];


(* ::Section:: *)
(*DATABASE ACCESS*)


ClearAll[ThermoLoadedQ, ThermoSpecies, ThermoRecord, ThermoDatabaseSummary];

ThermoLoadedQ[] := AssociationQ[$ThermoDatabase] && Length[$ThermoDatabase] > 0;

ThermoSpecies[] := Module[{},
  If[! ThermoLoadedQ[],
    Message[ThermoRecord::noload];
    Return[$Failed]
  ];
  Keys[$ThermoDatabase]
];

ThermoRecord[species_String] := Module[{},
  If[! ThermoLoadedQ[],
    Message[ThermoRecord::noload];
    Return[$Failed]
  ];

  If[! KeyExistsQ[$ThermoDatabase, species],
    Message[ThermoRecord::species, species];
    Return[$Failed]
  ];

  $ThermoDatabase[species]
];

ThermoDatabaseSummary[] := Module[{rows},
  If[! ThermoLoadedQ[],
    Message[ThermoRecord::noload];
    Return[$Failed]
  ];

  rows = Table[
    With[{r = $ThermoDatabase[sp]},
      <|
        "Species"    -> sp,
        "SourceLine" -> r["SourceLine"],
        "Tmin [K]"   -> r["Tmin"],
        "Tmid [K]"   -> r["Tmid"],
        "Tmax [K]"   -> r["Tmax"]
      |>
    ],
    {sp, Keys[$ThermoDatabase]}
  ];

  Dataset[rows]
];


(* ::Section:: *)
(*NASA-7 PROPERTY KERNEL*)


ClearAll[
  nasaCp,
  nasaH,
  nasaS,
  nasaG,
  thermoCoefficients
];

nasaCp[a_List, T_?NumericQ] :=
  $ThermoR (a[[1]] + a[[2]] T + a[[3]] T^2 + a[[4]] T^3 + a[[5]] T^4);

nasaH[a_List, T_?NumericQ] :=
  $ThermoR T (
    a[[1]] +
    a[[2]] T/2.0 +
    a[[3]] T^2/3.0 +
    a[[4]] T^3/4.0 +
    a[[5]] T^4/5.0 +
    a[[6]]/T
  );

nasaS[a_List, T_?NumericQ] :=
  $ThermoR (
    a[[1]] Log[T] +
    a[[2]] T +
    a[[3]] T^2/2.0 +
    a[[4]] T^3/3.0 +
    a[[5]] T^4/4.0 +
    a[[7]]
  );

nasaG[a_List, T_?NumericQ] := nasaH[a, T] - T nasaS[a, T];

thermoCoefficients[species_String, T_?NumericQ] := Module[{r},

  If[T <= 0,
    Message[ThermoRecord::tpos, T];
    Return[$Failed]
  ];

  r = ThermoRecord[species];
  If[r === $Failed, Return[$Failed]];

  If[T < r["Tmin"] || T > r["Tmax"],
    Message[
      ThermoRecord::temp,
      T,
      r["Tmin"],
      r["Tmax"],
      species
    ];
    Return[$Failed]
  ];

  If[T <= r["Tmid"], r["Low"], r["High"]]
];


(* ::Section:: *)
(*PUBLIC THERMODYNAMIC FUNCTIONS *)


ClearAll[ThermoCp, ThermoH, ThermoS, ThermoG];

ThermoCp[species_String, T_?NumericQ] := Module[{a},
  a = thermoCoefficients[species, T];
  If[a === $Failed, $Failed, nasaCp[a, T]]
];

ThermoH[species_String, T_?NumericQ] := Module[{a},
  a = thermoCoefficients[species, T];
  If[a === $Failed, $Failed, nasaH[a, T]]
];

ThermoS[species_String, T_?NumericQ] := Module[{a},
  a = thermoCoefficients[species, T];
  If[a === $Failed, $Failed, nasaS[a, T]]
];

ThermoG[species_String, T_?NumericQ] := Module[{a},
  a = thermoCoefficients[species, T];
  If[a === $Failed, $Failed, nasaG[a, T]]
];


(* ::Section:: *)
(*VALIDATION UTILITIES  *)


ClearAll[
  relativeJump,
  ValidateThermoContinuity,
  ThermoPropertyTable,
  ThermoBasicChecks
];

relativeJump[x_?NumericQ, y_?NumericQ] :=
  Abs[x - y]/Max[1.0, Abs[x], Abs[y]];

ValidateThermoContinuity[] := Module[{rows},
  If[! ThermoLoadedQ[],
    Message[ThermoRecord::noload];
    Return[$Failed]
  ];

  rows = Table[
    With[
      {
        r = $ThermoDatabase[sp],
        tm = $ThermoDatabase[sp]["Tmid"]
      },
      <|
        "Species" -> sp,
        "Tmid [K]" -> tm,
        "RelJumpCp" -> relativeJump[
          nasaCp[r["Low"], tm],
          nasaCp[r["High"], tm]
        ],
        "RelJumpH" -> relativeJump[
          nasaH[r["Low"], tm],
          nasaH[r["High"], tm]
        ],
        "RelJumpS" -> relativeJump[
          nasaS[r["Low"], tm],
          nasaS[r["High"], tm]
        ]
      |>
    ],
    {sp, Keys[$ThermoDatabase]}
  ];

  Dataset[rows]
];

ThermoPropertyTable[
  species_List : {"H2", "O2", "CH4", "CO", "CO2", "H2O"},
  T_?NumericQ : 298.15
] := Dataset @ Table[
  <|
    "Species" -> sp,
    "T [K]" -> T,
    "Cp [J/mol/K]" -> ThermoCp[sp, T],
    "H [kJ/mol]" -> ThermoH[sp, T]/1000.0,
    "S [J/mol/K]" -> ThermoS[sp, T],
    "G [kJ/mol]" -> ThermoG[sp, T]/1000.0
  |>,
  {sp, species}
];

(*
  Basic sanity checks for this THERM.DAT implementation.
  These are not a substitute for external validation; they are intended to
  reveal gross parsing or coefficient-order errors immediately.
*)

ThermoBasicChecks[] := Module[{h2h, o2h, ch4h, h2ocp},
  If[! ThermoLoadedQ[],
    Message[ThermoRecord::noload];
    Return[$Failed]
  ];

  h2h   = ThermoH["H2", 298.15]/1000.0;
  o2h   = ThermoH["O2", 298.15]/1000.0;
  ch4h  = ThermoH["CH4", 298.15]/1000.0;
  h2ocp = ThermoCp["H2O", 298.15];

  <|
    "LoadedSpeciesCount" -> Length[$ThermoDatabase],
    "ExpectedSpeciesCount" -> Length[$CombustorThermoSpecies],
    "H2_H_298_kJmol" -> h2h,
    "O2_H_298_kJmol" -> o2h,
    "CH4_H_298_kJmol" -> ch4h,
    "H2O_Cp_298_JmolK" -> h2ocp,
    "CountCheck" ->(Length[$ThermoDatabase] == Length[$CombustorThermoSpecies]),
    "ReferenceStateCheck" -> (Abs[h2h] < 0.1 && Abs[o2h] < 0.1),
    "CH4FormationEnthalpyCheck" -> (Abs[ch4h + 74.6] < 0.5)
  |>
];


(* ::Section::Closed:: *)
(*RECOMMENDED INITIALIZATION *)


(* ::Text:: *)
(*  If THERM.DAT is in the same directory as the notebook, evaluate:*)
(**)
(*      thermoPath = FileNameJoin[{NotebookDirectory[], "THERM.DAT"}];*)
(*      LoadThermoDatabase[thermoPath]*)
(**)
(*  Expected result:*)
(**)
(*      <|"Status" -> "Loaded", ..., "Count" -> 33, ...|>*)
(**)
(*  Recommended first tests:*)
(**)
(*      ThermoRecord["CH4"]*)
(*      ThermoDatabaseSummary[]*)
(*      ThermoBasicChecks[]*)
(*      ThermoPropertyTable[]*)
(*      ValidateThermoContinuity[]*)
