# Bacterial diversity during salinity-conveyed thermotolerance in the coral model Aiptasia

This repository contains the scripts used to analyze data and create figures for the manuscript “Salinity-conveyed thermotolerance in the coral model Aiptasia is accompanied by distinct changes of the bacterial microbiome”.

## Summary

Coral bleaching, i.e. the loss of photosynthetic algal endosymbionts, caused by ocean warming is now among the main factors driving global reef decline, making the elucidation of factors that contribute to thermotolerance important. Recent studies implicate high salinity as a contributing factor, potentially explaining the high thermotolerance of corals from the Arabian Seas. Here we characterized bacterial community composition under heat stress at different salinities using the coral model Aiptasia. Exposure of two Aiptasia host-algal symbiont pairings (H2-SSB01 and CC7-SSA01) to ambient (25°C) and heat stress (34°C) temperatures at low (36 PSU), intermediate (39 PSU), and high (42 PSU) salinities showed that bacterial community composition at high salinity was significantly different, concomitant with reduced bleaching susceptibility in H2-SSB01, not observed in CC7-SSA01. Elucidation of bacteria that showed increased relative abundance at high salinity, irrespective of heat stress, revealed candidate taxa that could potentially contribute to the observed increased thermotolerance. We identified 4 (H2-SSB01) and 3 (CC7-SSA01) bacterial taxa belonging to the orders Alteromondales, Rhizobiales, and Rhodobacterales, suggesting that only few bacterial taxa are potential contributors to an increase in thermal tolerance at high salinities. These taxa have previously been implicated in nitrogen and DMSP cycling, processes that are considered to affect thermotolerance. Our study is the first to demonstrate microbiome restructuring in symbiotic cnidarians under heat stress at different salinities. As such, it underlines the putative role of bacterial communities adapting to prevailing environmental conditions with consequences for the environmental stress tolerance of the associated animal host.

## Workflow

1. Identification and removal of contaminant OTUs, and overal taxonomic profiles were done using the script `Aipasia_salinity_barplots.R`

2. Ordination plots and PERMANOVAs were done using the script `Aiptasia_salinity_ordination_and_Permanovas.R`

3. Alpha diversity indexes were calculated using the script `Aiptasia_salinity_AlphaDiversity.R`

4. The Similarity Percentage analysis (SIMPER) was used to identify OTUs that best discriminated between salinities and temperatures as described in the script `Aiptasia_salinity_SIMPER.R`

