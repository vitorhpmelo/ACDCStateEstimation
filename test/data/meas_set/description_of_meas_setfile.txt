Each collumn of the meas_set file correspond to a type of measurement  

the column name is orgonized as (var)-(source)

-var = the monitored variables
-source = where the measurement comes from



AC side measurements (unsyc)
- vm-scada - Buses with voltage maginitude measurements (AC) 
- pinj-scada - Buses with power injections measurements (AC)
- p-scad - Branches with power flow measurements (AC), notation (bus from, bus to) 

AC side measurements (sync)
- v-pmu - Buses with voltage phasor measurements (AC)
- c-pmu - Branches with current measurements (AC), notation (bus from, bus to)

Converter measurements (unsyc)
- pconv_ac-conv — AC converter active power measurements
- pconv_pr_fr-conv — Converter power on the primary (from-side)
- pconv_tf_fr-conv — Transformer power on the from-side of the converter
- pconv_tf_to-conv — Transformer power on the to-side of the converter
- vmc-conv — Converter AC bus voltage magnitude
- vmf-conv — Converter filter bus voltage magnitude
- mconv-conv — Converter modulation index
- pconvdc-conv — Converter DC power measurement
- pconv_dc-conv — Converter DC power measurement (alternative notation)
- pconv_dcg-conv — DC grid-side converter power measurement
- pconv_dcg_shunt-conv — Converter shunt DC power measurement

DC grid measurement (unsyc) 
- i_dcgrid-dc — DC grid current measurements (bus from, bus to)
- vdcm-dc — DC grid midpoint voltage measurements
- p_dcgrid-dc — DC grid power measurements (bus from, bus to)

DC syncronized terminal unit (DTU) measurements (sync)
- i_dcgrid-dtu — DC current measurements (bus from, bus to)
- vdcm-dtu — DC voltage measurements
- p_dcgrid-dtu — DC power measurements (bus from, bus to)