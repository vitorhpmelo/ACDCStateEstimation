%
%   AC/DC grid OPF test case based on:
%   Cigre B4 DC grid benchmark model
%   
%  
%
%

function mpc = cigreb4
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type		Pd		Qd		  Gs	 Bs	   area		Vm	   Va		   baseKV	zone	  Vmax			 Vmin
mpc.bus = [
    1   3   0	0	0   0   1       1.06	0	400     1       1.1     0.9;
    2   3   0	0	0   0   1       1.06   	0	400     1       1.1     0.9;
    3   2   10	0	0   0   1       1.06   	0	220     1       1.1     0.9;
    4   2   10	0	0   0   1       1.06   	0	220     1       1.1     0.9;
];

%% generator data
%	bus	   Pg	  Qg		Qmax		 Qmin	  Vg	   mBase	status		Pmax	  Pmin	 Pc1	 Pc2   Qc1min	Qc1max    Qc2min	Qc2max	ramp_agc	    ramp_10	      ramp_30	  ramp_q	   apf         alpha
mpc.gen = [
    1	0       0	500      -500    1.06	100       1       2000     0 0 0 0 0 0 0 0 0 0 0 0;
    2	40      0	300      -300    1.06	100       1       2000     0 0 0 0 0 0 0 0 0 0 0 0;
    3   0      0	300      -300    1.06	100       1       2000     0 0 0 0 0 0 0 0 0 0 0 0;
    4   0      0	300      -300    1.06	100       1       2000     0 0 0 0 0 0 0 0 0 0 0 0;
];

%% branch data
%	fbus	tbus	r			 x			 b				  rateA	   rateB	   rateC	 ratio	  angle	    status	angmin	    angmax
mpc.branch = [
];

%% dc grid topology
%colunm_names% dcpoles
mpc.dcpol=2;
% numbers of poles (1=monopolar grid, 2=bipolar grid)
%% bus data
%column_names%  busdc_i    grid    Pdc     Vdc   basekVdc    Vdcmax   Vdcmin   Cdc
mpc.busdc = [
    1   1       0       1       525         1.1     0.9     0;
    2   1       0       1       525         1.1     0.9     0;
    3   1       0       1       525         1.1     0.9     0;
    4   1       0       1       525         1.1     0.9     0;
    5   1       0       1       525         1.1     0.9     0;
    6   1       0       1       525         1.1     0.9     0;
];



%% converters
%column_names% busdc_i busac_i type_dc type_ac P_g Q_g islcc Vtar  rtf xtf transformer tm   bf filter   rc   xc reactor basekVac Vmmax Vmmin Imax status LossA LossB LossCrec LossCinv  droop   Pdcset Vdcset dVdcset Pacmax Pacmin Qacmax Qacmin conv_confi connect_at ground_type ground_z status_p status_n
mpc.convdc = [
     1      1       2       1  -100 -40     0    1 0.00 0.01          1  1 0.01      1 0.00 0.01       1       138   1.1   0.9  1.1      1 1.103 0.887    2.885    1.885 0.0050 -58.6274 1.0079       0    100   -100     50    -50          2          0           1      0.5        1        1;
     2      2       3       1  -100   0     0    1 0.00 0.01          1  1 0.01      1 0.00 0.01       1      138   1.1   0.9  1.1      1 1.103 0.887    2.885    2.885 0.0070  21.9013 1.0000       0    100   -100     50    -50          2          0           0      0.5        1        1;
     5      3       3       1  -100   0     0    1 0.00 0.01          1  1 0.01      1 0.00 0.01       1      138   1.1   0.9  1.1      1 1.103 0.887    2.885    2.885 0.0070  21.9013 1.0000       0    100   -100     50    -50          2          0           0      0.5        1        1;
     6      4       3       1   50   0     0    1 0.00 0.01          1  1 0.01      1 0.00 0.01       1      138   1.1   0.9  1.1      1 1.103 0.887    2.885    2.885 0.0070  21.9013 1.0000       0    100   -100     50    -50          2          0           0      0.5        1        1;
];

%% dc branches
%column_names%  fbusdc  tbusdc  r  l   c   rateA   rateB   rateC   status   line_confi return_type return_z connect_at status_p status_n status_r
mpc.branchdc = [
    1   3   0.0689   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
    2   4   0.0689   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
    3   5   2.0667   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
    5   6   1.3778   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
    4   6   1.3778   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
    3   4   2.0667   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
    3   6   2.7556   0    0   2500   2500   2500   1  2           2    0.052          0        1        1        1; % bipolar
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
2	0	0	3	0  1	0;
2	0	0	3   0  2	0;
2	0	0	3   0  0.5	0;
2	0	0	3   0  0.5	0;
];
