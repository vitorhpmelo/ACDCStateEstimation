# ACDCStateEstimation

[![Build Status](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## WARNINGS

- Currently dc loads and generators are not supported
- We don't allow the DC linearization for the AC power flow equations - we are not vandals - all `DC's in this repository refer to actual direct current lines/features

## TODOs

- PRIORITY: add ACR formulations for AC side -> first update to new PowerModelsMCDC!!
- PRIORITY: add missing measurements and measurement conversions (see `synthetic_measurements.jl`) - below is a list of measurements to be added
- PRIORITY: add tests
- DISCUSS: which network cases are relevant? keep shortcuts to all?
- add docs
- add github workflows

## Notes

Based on the PowerModelsMCDC data format:
1. all dc buses have three terminals, even when some of them are not actually connected to anything. the order is (positive, negative, neutral)
2. for dc branches, you need to check first `line_confi`. If this is 2, all terminals are connected. Only if it is 1, go look at `connect_at`
3. the above holds for converters, too

NOTE: in the opf results, dc powers are always zero because they are dummy placeholders. multiply current and voltage in post-processing for the actual values!!

## Open questions

- what to do if we have multiple generators/slackbuses? do we just fix angle of one of them and let the estimator find the others, or do we fix them all?

## Measurements that still need to be added

### All, prioritary: enable ``to" flow measurements

### For ACPolar, AC side
- current injections (generators and loads)
- current flows (branches)

### For ACRectangular, AC side
- ????? DISCUSS

### DC measurements
- power and current injections to dc buses (from converter)
- other converter stuff????

## Changes to PMMCDC that might be nice
- do not model non-utilized terminals of DC buses
- make the mapping of multi-conductor components a bit more streamlines (a bit like PMD), instead of the lookup settings like right now
- remove the dc_power_flow variable and power = 0 placeholder
- inconsistency in naming, e.g., fbusdc vs f_bus
- lots of commented code in `opf_dcp`
- matteo mentioning something about grounding (although not a problem here apparently)