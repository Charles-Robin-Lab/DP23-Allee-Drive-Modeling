#!/usr/bin/env wolframscript



lociEffectPMFR[0, n_] := {{1,PDF[BinomialDistribution[2*n, p], 0]}}

lociEffectPMFR[k_, n_] := Append[lociEffectPMFR[k-1,n],{1-(k/2/n)^2,PDF[BinomialDistribution[2*n, p], k]}]
lociEffectPMF[individuals_] := lociEffectPMFR[2*individuals, individuals]



w :=lociEffectPMF[4]

collapse[c1_,c2_] := {First[c1],Last[c1]+Last[c2]}

crossDists[x_,y_] := 
    Module[{values,probabilitiesuncollapsedDist},
    values := Flatten[Outer[Times,Map[First,x],Map[First,y]]];
    probabilities := Flatten[Outer[Times,Map[Last,x],Map[Last,y]]];
    uncollapsedDist := MapThread[List, {values, probabilities}];
    Map[Fold[collapse],Values[GroupBy[uncollapsedDist, First]]]]

cross[x_] := crossDists[x,w]
Print[w]
Print[Nest[cross,w,1]]
