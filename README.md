# ACDCStateEstimation

[![Build Status](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MartaVanin/ACDCStateEstimation.jl/actions/workflows/CI.yml?query=branch%3Amaster)


## WARNINGS

- Currently dc loads are not supported
- We don't allow the DC formulation for the AC side, we are not vandals

## TODOs

- import export PMMCDC parsers so that they are not needed in the examples. add preferential functions for 5-bus ntw and 67-bus ntw parsing
- remove `constraint_voltage_dc` altogether from problem description?

## Notes

Based on the PowerModelsMCDC data format:
1. all DC buses have three terminals, even unused ones
2. for dc branches, you need to check first `line_confi`. If this is 2, all terminals are connected. Only if it is 1, go look at `connect_at`
3. above is the same for converters, too

NOTE: in the opf results, dc powers are always zero because they are dummy placeholders. multiply current and voltage in post-processing for the actual values!!

## Open questions

- what to do if we have multiple generators/slackbuses? do we just fix angle of one of them and let the estimator find the others, or do we fix them all?