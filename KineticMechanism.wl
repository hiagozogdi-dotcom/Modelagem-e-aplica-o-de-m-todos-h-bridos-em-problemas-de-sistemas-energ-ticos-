(* ::Package:: *)

(* ::Title:: *)
(*KINETICMECHANISM*)
(*Chemical-Mechanism Data Library for Gas-Turbine Combustor Modeling*)


(* ::Text:: *)
(*Description*)
(**)
(*  Wolfram Language library for importing and organizing the reduced*)
(*  DRM19 chemical mechanism and appending simplified NOx and SOx subsets.*)
(**)
(*  IMPORTANT:*)
(*  *)
(*  This module is a MECHANISM DATA LIBRARY. It does NOT evaluate:*)
(*      - forward/reverse rate constants k(T)*)
(*      - reaction progress rates q*)
(*      - species production rates omega_k*)
(*      - reactor balances*)
(**)
(*  These calculations are intentionally deferred to the combustor model.*)
(**)
(*  The library provides:*)
(*      - species list*)
(*      - reaction equations*)
(*      - reactant/product stoichiometric coefficients*)
(*      - net stoichiometric coefficients*)
(*      - Arrhenius parameters A, beta, Ea*)
(*      - reaction reversibility*)
(*      - third-body flags and efficiencies*)
(*      - LOW/falloff parameters*)
(*      - TROE parameters*)
(*      - source/submechanism metadata*)
(*      - stoichiometric matrices*)
(**)
(*  Project*)
(*  Gas-Turbine Combustor Model*)
(*  One-Dimensional Reactive Flow and Chemical-Kinetics Framework*)
(**)
(*  Author*)
(*  Hiago David Zogbi Silva Oliveira*)
(**)
(*  Universidade Federal Fluminense (UFF)*)
(*  Department of Mechanical Engineering*)
(*  Laboratory of Thermosciences (LATERMO)*)
(**)
(*  Required mechanism file*)
(*  drm19.dat*)
(**)
(*  The supplied DRM19 file identifies itself as: "Reduced version of GRI-MECH 1.2. 19 species (+ N2, AR); 84 reactions. PennState Dec, 1994."*)
(**)
(*  Extended mechanism used by this library*)
(*  Core:*)
(*      DRM19                    21 species / 84 reactions*)
(**)
(*  Simplified NOx subset: N, NO, NO2 Extended-Zeldovich thermal-NO reactions plus reduced NO/NO2 interconversion reactions.*)
(**)
(*  Simplified SOx subset: SO, SO2, SO3 Reduced oxidation scheme representing SO -> SO2 and pressure-dependent SO2 -> SO3 conversion.*)
(**)
(*  Total after extension:*)
(*      27 species*)
(*      92 reactions*)
(**)
(*  NOTE:*)
(*  N and SO3 must also be available in the thermodynamic database.*)
(**)
(*  CHEMKIN parameter convention*)
(*  *)
(*  The supplied DRM19 REACTIONS block does not explicitly override the traditional CHEMKIN unit convention. Raw parameters are therefore stored using the mechanism convention:*)
(**)
(*      A     : CHEMKIN concentration units (order dependent)*)
(*      beta  : dimensionless*)
(*      Ea    : cal/mol*)
(*      T     : K*)
(*      time  : s*)
(*      concentration basis : mol/cm^3*)
(**)
(*  No SI conversion of A is performed in this library. This is deliberate: conversion of A depends on reaction order and belongs to the kinetics evaluation layer.*)
(**)
(*  Main references*)
(*  *)
(*  DRM19 / CHEMKIN:*)
(*  [1] Kee, R. J.; Rupley, F. M.; Miller, J. A. CHEMKIN-II: A FORTRAN Chemical Kinetics Package for the Analysis of Gas-Phase Chemical Kinetics. Sandia National Laboratories, SAND89-8009B.*)
(*  [2] DRM19 mechanism file supplied with this project: Reduced GRI-MECH 1.2 mechanism, Penn State, Dec. 1994.*)
(*  NOx subset:*)
(*  [3] Smith, G. P. et al. GRI-Mech 3.0.*)
(*  *)
(*      http://www.me.berkeley.edu/gri_mech/*)
(**)
(*      The thermal-NO subset retains:*)
(*         N + NO  <=> N2 + O*)
(*         N + O2  <=> NO + O*)
(*         N + OH  <=> NO + H*)
(**)
(*      Reduced NO2 chemistry retains:*)
(*         HO2 + NO <=> NO2 + OH*)
(*         NO + O + M <=> NO2 + M*)
(*         NO2 + O <=> NO + O2*)
(**)
(*  SOx subset:*)
(*  *)
(*  [4] Garland, N. L. Temperature dependence of the reaction SO + O2. Chemical Physics Letters, 290 (1998), 385-390. doi:10.1016/S0009-2614(98)00553-3.*)
(**)
(*  [5] Mueller, M. A.; Yetter, R. A.; Dryer, F. L.  Kinetic modeling of the CO/H2O/O2/NO/SO2 system: Implications for high-pressure fall-off in the SO2 + O(+M) = SO3(+M) reaction. International Journal of Chemical Kinetics, 32 (2000), 317-339.*)
(**)
(*  [6] Rasmussen/Glarborg sulfur-chemistry parameterization as compiled in later detailed sulfur-combustion mechanisms. The reduced  SO2 + O(+M) <=> SO3(+M) entry retains high-pressure, low-pressure and Troe parameters.*)
(**)
(*  Version: 1.0*)
(*  Date: September 2026*)


(* ::Section::Closed:: *)
(*PUBLIC SYMBOLS AND GLOBAL STORAGE  *)


ClearAll[
  LoadKineticMechanism,
  MechanismLoadedQ,
  MechanismSummary,
  MechanismSpecies,
  MechanismReactions,
  MechanismReactionCount,
  ReactionRecord,
  ReactionReactants,
  ReactionProducts,
  ReactionArrhenius,
  ReactionThirdBodyQ,
  ReactionEfficiencies,
  ReactionFalloffQ,
  ReactionLowPressure,
  ReactionTroe,
  ReactantStoichiometricMatrix,
  ProductStoichiometricMatrix,
  StoichiometricMatrix,
  MechanismUnits,
  MechanismReferences,
  MissingMechanismSpecies,
  ReactionsBySubset,
  $KineticMechanism
];

$KineticMechanism = <||>;


(* ::Section::Closed:: *)
(*MESSAGES  *)


LoadKineticMechanism::nofile =
  "Mechanism file `1` was not found.";

LoadKineticMechanism::section =
  "Could not locate a valid `1` section in the mechanism file.";

LoadKineticMechanism::rxn =
  "Could not parse reaction line: `1`.";

LoadKineticMechanism::count =
  "Expected 84 DRM19 reactions, but `1` reactions were parsed.";

ReactionRecord::noload =
  "No kinetic mechanism is loaded. Call LoadKineticMechanism[path] first.";

ReactionRecord::id =
  "Reaction index `1` is outside the valid range.";

MechanismSpecies::noload = ReactionRecord::noload;
MechanismReactions::noload = ReactionRecord::noload;
MechanismSummary::noload = ReactionRecord::noload;
MechanismUnits::noload = ReactionRecord::noload;


(* ::Section::Closed:: *)
(*BASIC UTILITIES    *)


MechanismLoadedQ[] :=
  AssociationQ[$KineticMechanism] &&
  KeyExistsQ[$KineticMechanism, "Reactions"];

normalizeLine[s_String] :=
  StringTrim @ First @ StringSplit[s, "!", 2];

numberStringQ[s_String] :=
  StringMatchQ[
    StringTrim[s],
    RegularExpression[
      "[+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[EeDd][+-]?\\d+)?"
    ]
  ];

toNumber[s_String] := Module[
  {x},

  x = StringTrim[s];

  x = StringReplace[
    x,
    RegularExpression["[EeDd]([+-]?\\d+)"] :> "*^$1"
  ];

  ToExpression[x]
];


(* ::Section::Closed:: *)
(*STOICHIOMETRIC-SIDE PARSER *)


ClearAll[parseStoichToken, parseStoichSide];

parseStoichToken[token_String] := Module[
  {t, m, coeff, sp},
  t = StringTrim[token];
  If[
  t === "" ||
  t === "M" ||
  t === "m",
  Return[Nothing]
  ];

  m = StringCases[
    t,
    RegularExpression[
      "^([0-9]+(?:\\.[0-9]+)?)?(.+)$"
    ] -> {"$1", "$2"}
  ];

  If[m === {},
    Return[$Failed]
  ];

  coeff =
    If[
      m[[1, 1]] === "",
      1.,
      ToExpression[m[[1, 1]]]
    ];

  sp = StringTrim[m[[1, 2]]];

  sp -> N[coeff]
];

parseStoichSide[side_String] := Module[
  {clean, tokens, rules},

  (* Remove pressure-dependent third-body notation
     BEFORE splitting species by "+" *)
     
  clean = StringReplace[
    side,
    {
      "(+M)" -> "",
      "(+m)" -> ""
    }
  ];

  tokens = StringSplit[clean, "+"];

  rules = DeleteCases[
    parseStoichToken /@ tokens,
    Nothing | $Failed
  ];

  Merge[rules, Total]
];


(* ::Section::Closed:: *)
(*REACTION-LINE PARSER *)


ClearAll[parseReactionHeader];

parseReactionHeader[line_String] := Module[
  {clean, parts, eqn, a, beta, ea, reversible, lr, lhs, rhs},

  clean = normalizeLine[line];

  parts = StringSplit[clean];

  If[Length[parts] < 4,
    Message[LoadKineticMechanism::rxn, line];
    Return[$Failed]
  ];

  (* Last three whitespace-separated fields are A, beta, Ea. *)
  If[!And @@ (numberStringQ /@ Take[parts, -3]),
    Message[LoadKineticMechanism::rxn, line];
    Return[$Failed]
  ];

  {a, beta, ea} = toNumber /@ Take[parts, -3];
  eqn = StringJoin @ Drop[parts, -3];

  reversible = StringContainsQ[eqn, "<=>"];

  lr = Which[
    reversible, StringSplit[eqn, "<=>"],
    StringContainsQ[eqn, "=>"], StringSplit[eqn, "=>"],
    True,
      Message[LoadKineticMechanism::rxn, line];
      Return[$Failed]
  ];

  If[Length[lr] != 2,
    Message[LoadKineticMechanism::rxn, line];
    Return[$Failed]
  ];

  {lhs, rhs} = lr;

  <|
    "Equation" -> eqn,
    "Reactants" -> parseStoichSide[lhs],
    "Products" -> parseStoichSide[rhs],
    "Reversible" -> reversible,

    "Arrhenius" -> <|
      "A" -> N[a],
      "Beta" -> N[beta],
      "Ea" -> N[ea]
    |>,

    "ThirdBody" ->
      (
        StringContainsQ[eqn, "+M"] ||
        StringContainsQ[eqn, "(+M)"]
      ),

    "Falloff" -> StringContainsQ[eqn, "(+M)"],

    "Efficiencies" -> <||>,
    "LowPressure" -> None,
    "Troe" -> None
  |>
];


(* ::Section::Closed:: *)
(*CONTINUATION-LINE PARSERS*)


ClearAll[
  parseTripleSlashLine,
  parseTroeLine,
  parseEfficiencyLine
];

parseTripleSlashLine[line_String, keyword_String] := Module[
  {inside, vals},

  inside = StringCases[
    line,
    RegularExpression[
      "(?i)" <> keyword <> "\\s*/\\s*([^/]*)/"
    ] -> "$1"
  ];

  If[inside === {}, Return[None]];

  vals = Select[StringSplit[StringTrim @ First[inside]], numberStringQ];

  If[Length[vals] < 3, Return[None]];

  vals = toNumber /@ Take[vals, 3];

  <|
    "A" -> N[vals[[1]]],
    "Beta" -> N[vals[[2]]],
    "Ea" -> N[vals[[3]]]
  |>
];

parseTroeLine[line_String] := Module[
  {inside, vals},

  inside = StringCases[
    line,
    RegularExpression["(?i)TROE\\s*/\\s*([^/]*)/"] -> "$1"
  ];

  If[inside === {}, Return[None]];

  vals = Select[StringSplit[StringTrim @ First[inside]], numberStringQ];
  vals = N[toNumber /@ vals];

  Which[
    Length[vals] == 3,
      <|
        "Alpha" -> vals[[1]],
        "T3" -> vals[[2]],
        "T1" -> vals[[3]],
        "T2" -> None
      |>,

    Length[vals] >= 4,
      <|
        "Alpha" -> vals[[1]],
        "T3" -> vals[[2]],
        "T1" -> vals[[3]],
        "T2" -> vals[[4]]
      |>,

    True,
      None
  ]
];

parseEfficiencyLine[line_String] := Module[
  {pairs},

  pairs = StringCases[
    line,
    RegularExpression[
      "([A-Za-z0-9()]+)\\s*/\\s*([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[EeDd][+-]?\\d+)?)\\s*/"
    ] :> ("$1" -> "$2")
  ];

  Association[
    (#1 -> N[toNumber[#2]]) & @@@ pairs
  ]
];


(* ::Section::Closed:: *)
(*DRM19 FILE PARSER *)


ClearAll[parseDRM19];

parseDRM19[text_String] := Module[
  {
    lines, norm, speciesStart, speciesEnd, reactionStart, reactionEnd,
    speciesLines, species, rxnLines, reactions, current, i, line,
    low, troe, eff
  },

  lines = StringSplit[
    StringReplace[text, {"\r\n" -> "\n", "\r" -> "\n"}],
    "\n"
  ];

  norm = normalizeLine /@ lines;

  speciesStart = FirstPosition[norm, "SPECIES", Missing["NotFound"]];
  If[MissingQ[speciesStart],
    Message[LoadKineticMechanism::section, "SPECIES"];
    Return[$Failed]
  ];
  speciesStart = First[speciesStart];

  speciesEnd = FirstPosition[
    norm[[speciesStart + 1 ;;]],
    "END",
    Missing["NotFound"]
  ];
  If[MissingQ[speciesEnd],
    Message[LoadKineticMechanism::section, "SPECIES/END"];
    Return[$Failed]
  ];
  speciesEnd = speciesStart + First[speciesEnd];

  speciesLines = Select[
    norm[[speciesStart + 1 ;; speciesEnd - 1]],
    # =!= "" &
  ];
  species = Flatten[StringSplit /@ speciesLines];

  reactionStart = FirstPosition[norm, "REACTIONS", Missing["NotFound"]];
  If[MissingQ[reactionStart],
    Message[LoadKineticMechanism::section, "REACTIONS"];
    Return[$Failed]
  ];
  reactionStart = First[reactionStart];

  reactionEnd = FirstPosition[
    norm[[reactionStart + 1 ;;]],
    "END",
    Missing["NotFound"]
  ];
  If[MissingQ[reactionEnd],
    Message[LoadKineticMechanism::section, "REACTIONS/END"];
    Return[$Failed]
  ];
  reactionEnd = reactionStart + First[reactionEnd];

  rxnLines = norm[[reactionStart + 1 ;; reactionEnd - 1]];
  reactions = {};
  current = None;

  Do[
    line = StringTrim[rxnLines[[i]]];

    If[line === "", Continue[]];

    Which[
      StringContainsQ[line, "<=>"] || StringContainsQ[line, "=>"],

        If[current =!= None,
          AppendTo[reactions, current]
        ];

        current = parseReactionHeader[line];

        If[current === $Failed,
          Return[$Failed]
        ],

      StringMatchQ[line, RegularExpression["(?i)^LOW\\s*/.*"]],

        If[current =!= None,
          low = parseTripleSlashLine[line, "LOW"];
          current["LowPressure"] = low;
          current["Falloff"] = True;
          current["ThirdBody"] = True;
        ],

      StringMatchQ[line, RegularExpression["(?i)^TROE\\s*/.*"]],

        If[current =!= None,
          troe = parseTroeLine[line];
          current["Troe"] = troe;
          current["Falloff"] = True;
          current["ThirdBody"] = True;
        ],

      StringContainsQ[line, "/"],

        If[current =!= None,
          eff = parseEfficiencyLine[line];
          If[Length[eff] > 0,
            current["Efficiencies"] =
              Join[current["Efficiencies"], eff]
          ];
        ],

      True,
        Null
    ],
    {i, Length[rxnLines]}
  ];

  If[current =!= None,
    AppendTo[reactions, current]
  ];

  reactions = MapIndexed[
    Join[
      #1,
      <|
        "GlobalID" -> First[#2],
        "LocalID" -> First[#2],
        "Subset" -> "DRM19",
        "Source" -> "DRM19 (reduced GRI-MECH 1.2, Penn State, Dec. 1994)"
      |>
    ] &,
    reactions
  ];

  <|
    "Species" -> species,
    "Reactions" -> reactions
  |>
];


(* ::Section::Closed:: *)
(*SIMPLIFIED NOx SUBMECHANISM *)


ClearAll[makeReactionRecord];

makeReactionRecord[
  eqn_String,
  a_?NumericQ,
  beta_?NumericQ,
  ea_?NumericQ,
  subset_String,
  source_String,
  opts___Rule
] := Module[
  {base, rules},

  base = parseReactionHeader[
    eqn <> " " <>
    ToString[FortranForm @ N[a]] <> " " <>
    ToString[FortranForm @ N[beta]] <> " " <>
    ToString[FortranForm @ N[ea]]
  ];

  rules = Association[{opts}];

  Join[
    base,
    <|
      "Subset" -> subset,
      "Source" -> source
    |>,
    rules
  ]
];


$NOxReactions = {

  makeReactionRecord[
    "N+NO<=>N2+O",
    2.700*^13, 0.0, 355.0,
    "NOx",
    "GRI-Mech 3.0 / extended Zeldovich"
  ],

  makeReactionRecord[
    "N+O2<=>NO+O",
    9.000*^9, 1.0, 6500.0,
    "NOx",
    "GRI-Mech 3.0 / extended Zeldovich"
  ],

  makeReactionRecord[
    "N+OH<=>NO+H",
    3.360*^13, 0.0, 385.0,
    "NOx",
    "GRI-Mech 3.0 / extended Zeldovich"
  ],

  makeReactionRecord[
    "HO2+NO<=>NO2+OH",
    2.110*^12, 0.0, -480.0,
    "NOx",
    "GRI-Mech 3.0 reduced NO/NO2 chemistry"
  ],

  makeReactionRecord[
    "NO+O+M<=>NO2+M",
    1.060*^20, -1.410, 0.0,
    "NOx",
    "GRI-Mech 3.0 reduced NO/NO2 chemistry",
    "ThirdBody" -> True,
    "Efficiencies" -> <|
      "H2" -> 2.0,
      "H2O" -> 6.0,
      "CH4" -> 2.0,
      "CO" -> 1.5,
      "CO2" -> 2.0
    |>
  ],

  makeReactionRecord[
    "NO2+O<=>NO+O2",
    3.900*^12, 0.0, -240.0,
    "NOx",
    "GRI-Mech 3.0 reduced NO/NO2 chemistry"
  ]
};


(* ::Section::Closed:: *)
(*SIMPLIFIED SOx SUBMECHANISM*)


$SOxReactions = {

  (* SO + O2 -> SO2 + O is the dominant reduced SO oxidation step.
     Parameters are retained in CHEMKIN-style mol-cm-s units. *)
     
  makeReactionRecord[
    "SO+O2<=>SO2+O",
    8.91*^6, 1.4, 3712.0,
    "SOx",
    "Garland (1998); sulfur mechanism compilation"
  ],

  (* Pressure-dependent SO2 + O association.
     High-pressure Arrhenius:
         A = 3.70E11, beta = 0, Ea = 1689 cal/mol
     Low-pressure limit:
         A = 2.40E27, beta = -3.60, Ea = 5186 cal/mol
     Troe:
         alpha = 0.442, T3 = 316 K, T1 = 7442 K
     Reduced third-body efficiencies retained from sulfur mechanisms. *)
     
  makeReactionRecord[
    "SO2+O(+M)<=>SO3(+M)",
    3.70*^11, 0.0, 1689.0,
    "SOx",
    "Mueller et al. (2000); Rasmussen/Glarborg sulfur chemistry",
    "ThirdBody" -> True,
    "Falloff" -> True,
    "LowPressure" -> <|
      "A" -> 2.40*^27,
      "Beta" -> -3.60,
      "Ea" -> 5186.0
    |>,
    "Troe" -> <|
      "Alpha" -> 0.442,
      "T3" -> 316.0,
      "T1" -> 7442.0,
      "T2" -> None
    |>,
    "Efficiencies" -> <|
      "N2" -> 0.0,
      "CO2" -> 2.5,
      "SO2" -> 10.0,
      "H2O" -> 10.0
    |>
  ]
};


(* ::Section::Closed:: *)
(*STOICHIOMETRIC MATRICES *)


ClearAll[buildStoichMatrices];

buildStoichMatrices[species_List, reactions_List] := Module[
  {nr, ns, nuR, nuP, nuN},

  ns = Length[species];
  nr = Length[reactions];

  nuR = Table[
    Lookup[reactions[[r, "Reactants"]], species[[k]], 0.0],
    {k, ns}, {r, nr}
  ];

  nuP = Table[
    Lookup[reactions[[r, "Products"]], species[[k]], 0.0],
    {k, ns}, {r, nr}
  ];

  nuN = nuP - nuR;

  <|
    "Reactants" -> nuR,
    "Products" -> nuP,
    "Net" -> nuN
  |>
];


(* ::Section::Closed:: *)
(*MECHANISM LOADER*)


Options[LoadKineticMechanism] = {
  "IncludeNOx" -> True,
  "IncludeSOx" -> True
};

LoadKineticMechanism[
  path_String,
  OptionsPattern[]
] := Module[
  {
    text, drm, coreReactions, extraReactions, reactions,
    species, matrices, includeNOx, includeSOx, allReactions
  },

  If[!FileExistsQ[path],
    Message[LoadKineticMechanism::nofile, path];
    Return[$Failed]
  ];

  text = Import[path, "Text"];

  drm = parseDRM19[text];

  If[drm === $Failed,
    Return[$Failed]
  ];

  coreReactions = drm["Reactions"];

  If[Length[coreReactions] != 84,
    Message[LoadKineticMechanism::count, Length[coreReactions]];
    Return[$Failed]
  ];

  includeNOx = TrueQ[OptionValue["IncludeNOx"]];
  includeSOx = TrueQ[OptionValue["IncludeSOx"]];

  extraReactions = Join[
    If[includeNOx, $NOxReactions, {}],
    If[includeSOx, $SOxReactions, {}]
  ];

  allReactions = Join[coreReactions, extraReactions];

  (* Assign global IDs and sequential local IDs inside each subset. *)
  Module[{localCounters = <||>},
    allReactions = MapIndexed[
      Function[{rxn, idx},
        Module[{sub, local},
          sub = Lookup[rxn, "Subset", "Unspecified"];
          local = Lookup[localCounters, sub, 0] + 1;
          AssociateTo[localCounters, sub -> local];
          Join[
            rxn,
            <|
              "GlobalID" -> First[idx],
              "LocalID" -> local
            |>
          ]
        ]
      ],
      allReactions
    ];
  ];

  species = DeleteDuplicates @ Join[
    drm["Species"],
    If[includeNOx, {"N", "NO", "NO2"}, {}],
    If[includeSOx, {"SO", "SO2", "SO3"}, {}]
  ];

  matrices = buildStoichMatrices[species, allReactions];

  $KineticMechanism = <|

    "Name" -> "DRM19 + reduced NOx/SOx",

    "CoreMechanism" -> "DRM19",

    "Species" -> species,

    "Reactions" -> allReactions,

    "Matrices" -> matrices,

    "Units" -> <|
      "Temperature" -> "K",
      "Time" -> "s",
      "Concentration" -> "mol/cm^3",
      "ActivationEnergy" -> "cal/mol",
      "ArrheniusA" -> "CHEMKIN order-dependent units"
    |>,

    "Options" -> <|
      "IncludeNOx" -> includeNOx,
      "IncludeSOx" -> includeSOx
    |>,

    "References" -> <|
      "DRM19" ->
        "Reduced GRI-MECH 1.2 mechanism, Penn State, Dec. 1994.",
      "NOx" ->
        "GRI-Mech 3.0 thermal-NO and reduced NO/NO2 subset.",
      "SOx_SO_O2" ->
        "Garland, Chemical Physics Letters 290 (1998), 385-390.",
      "SOx_SO3" ->
        "Mueller, Yetter & Dryer, Int. J. Chem. Kinet. 32 (2000), 317-339; Rasmussen/Glarborg sulfur chemistry."
    |>
  |>;

  MechanismSummary[]
];


(* ::Section::Closed:: *)
(*PUBLIC QUERY FUNCTIONS*)


MechanismSpecies[] := Module[{},
  If[!MechanismLoadedQ[],
    Message[MechanismSpecies::noload];
    Return[$Failed]
  ];
  $KineticMechanism["Species"]
];

MechanismReactions[] := Module[{},
  If[!MechanismLoadedQ[],
    Message[MechanismReactions::noload];
    Return[$Failed]
  ];
  $KineticMechanism["Reactions"]
];

MechanismReactionCount[] :=
  If[MechanismLoadedQ[],
    Length[$KineticMechanism["Reactions"]],
    Message[MechanismReactions::noload]; $Failed
  ];

ReactionRecord[id_Integer] := Module[{nr},
  If[!MechanismLoadedQ[],
    Message[ReactionRecord::noload];
    Return[$Failed]
  ];

  nr = Length[$KineticMechanism["Reactions"]];

  If[id < 1 || id > nr,
    Message[ReactionRecord::id, id];
    Return[$Failed]
  ];

  $KineticMechanism["Reactions"][[id]]
];

ReactionReactants[id_Integer] :=
  Lookup[ReactionRecord[id], "Reactants", $Failed];

ReactionProducts[id_Integer] :=
  Lookup[ReactionRecord[id], "Products", $Failed];

ReactionArrhenius[id_Integer] :=
  Lookup[ReactionRecord[id], "Arrhenius", $Failed];

ReactionThirdBodyQ[id_Integer] :=
  TrueQ @ Lookup[ReactionRecord[id], "ThirdBody", False];

ReactionEfficiencies[id_Integer] :=
  Lookup[ReactionRecord[id], "Efficiencies", <||>];

ReactionFalloffQ[id_Integer] :=
  TrueQ @ Lookup[ReactionRecord[id], "Falloff", False];

ReactionLowPressure[id_Integer] :=
  Lookup[ReactionRecord[id], "LowPressure", None];

ReactionTroe[id_Integer] :=
  Lookup[ReactionRecord[id], "Troe", None];

ReactantStoichiometricMatrix[] :=
  If[MechanismLoadedQ[],
    $KineticMechanism["Matrices"]["Reactants"],
    Message[ReactionRecord::noload]; $Failed
  ];

ProductStoichiometricMatrix[] :=
  If[MechanismLoadedQ[],
    $KineticMechanism["Matrices"]["Products"],
    Message[ReactionRecord::noload]; $Failed
  ];

StoichiometricMatrix[] :=
  If[MechanismLoadedQ[],
    $KineticMechanism["Matrices"]["Net"],
    Message[ReactionRecord::noload]; $Failed
  ];

MechanismUnits[] :=
  If[MechanismLoadedQ[],
    $KineticMechanism["Units"],
    Message[MechanismUnits::noload]; $Failed
  ];

MechanismReferences[] :=
  If[MechanismLoadedQ[],
    $KineticMechanism["References"],
    Message[ReactionRecord::noload]; $Failed
  ];

ReactionsBySubset[subset_String] :=
  If[MechanismLoadedQ[],
    Select[
      $KineticMechanism["Reactions"],
      Lookup[#, "Subset", ""] === subset &
    ],
    Message[ReactionRecord::noload]; $Failed
  ];

MissingMechanismSpecies[availableSpecies_List] :=
  If[MechanismLoadedQ[],
    Complement[MechanismSpecies[], availableSpecies],
    Message[ReactionRecord::noload]; $Failed
  ];


(* ::Section::Closed:: *)
(*SUMMARY *)


MechanismSummary[] := Module[
  {rxns, subsets},

  If[!MechanismLoadedQ[],
    Message[MechanismSummary::noload];
    Return[$Failed]
  ];

  rxns = $KineticMechanism["Reactions"];
  subsets = Counts @ Lookup[rxns, "Subset"];

  <|
    "Status" -> "Loaded",
    "Name" -> $KineticMechanism["Name"],
    "SpeciesCount" -> Length[$KineticMechanism["Species"]],
    "ReactionCount" -> Length[rxns],
    "DRM19Reactions" -> Lookup[subsets, "DRM19", 0],
    "NOxReactions" -> Lookup[subsets, "NOx", 0],
    "SOxReactions" -> Lookup[subsets, "SOx", 0],
    "ThirdBodyReactions" ->
      Count[Lookup[rxns, "ThirdBody"], True],
    "FalloffReactions" ->
      Count[Lookup[rxns, "Falloff"], True],
    "TroeReactions" ->
      Count[
        Lookup[rxns, "Troe"],
        x_ /; AssociationQ[x]
      ],
    "Units" -> $KineticMechanism["Units"]
  |>
];
