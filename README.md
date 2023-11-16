# ACDCStateEstimation

[![Build Status](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml?query=branch%3Amaster)


## WARNINGS

- Currently dc loads are not supported
- We don't allow the DC formulation for the AC side, we are not vandals

## TODOs

- import export PMMCDC parsers so that they are not needed in the examples. add preferential functions for 5-bus ntw and 67-bus ntw parsing
- find a pretty way to choose to "drop" neutral/ground dc voltage, current and power measurements ????
- add missing measurements and measurement conversions (see `synthetic_measurements.jl`) - below is a list of measurements to be added

## Notes

Based on the PowerModelsMCDC data format:
1. all dc buses have three terminals, even when some of them are not actually connected to anything. the order is (positive, negative, neutral)
2. for dc branches, you need to check first `line_confi`. If this is 2, all terminals are connected. Only if it is 1, go look at `connect_at`
3. the above holds for converters, too

NOTE: in the opf results, dc powers are always zero because they are dummy placeholders. multiply current and voltage in post-processing for the actual values!!

## Open questions

- what to do if we have multiple generators/slackbuses? do we just fix angle of one of them and let the estimator find the others, or do we fix them all?
- no ACRPowerModel for the AC side in PMMCDC? why?

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