slurmstepd: error: task/cgroup: unable to add task[pid=30699] to memory cg '(null)'
MANDATORY INPUT VARIABLES:
  NStates =
		 2
  NSteps =
		 2
  dtI =
		 0.25 (fs)
  ESTEPS =
		 50
  ISTATE =
		 1
  FUNCTIONAL =
		 WB97XD
  BASIS_SET =
		 6-311++G**
  CHARGE =
		 0
  MULTIPLICITY =
		 1
  MEMORY =
		 5 (GB)
  PARALLEL_FORCES =
		 True
  NAMD_METHOD =
		 EH
  EL_PROP =
		 VV
  MD_ENSEMBLE =
		 NVE
  VELOC =
		 MB
Input looks good.
Submitting electronic structure for step 0.
	Gaussian terminated normally in 3.27 s. (GS_NEW)
	Gaussian terminated normally in 9.57 s. (TD_NEW_S1)
Submitting electronic structure for step 1.
	Gaussian terminated normally in 2.96 s. (GS_NEW)
	Gaussian terminated normally in 11.12 s. (TD_NEW_S1)
	Gaussian terminated normally in 1.63 s. (DIMER)
	TD_NEW_S1 DIM_INFO TIME: 0.05 s
	TD_NEW_S1 GET_MO TIME: 0.01 s
	TD_NEW_S1 GET_CI TIME: 0.01 s
	TD_OLD_S1 GET_MO TIME: 0.01 s
	TD_OLD_S1 GET_CI TIME: 0.01 s
	DIMER TIME: 2.01 s
	G16_NAC.py main: 2.08 s
	gau_nac() (G16_NAC.py): 0.0 s
gau_nac prepare 0.01 s
WFN OVERLAP FORTRAN TOOK 0.13 s
gau_nac dump 0.86 s
	gau_nac.worker() (G16_NAC.py): 1.0 s
	GET_OVERLAP OVERALL TIME (G16_TD.py): 3.11 s
	GET_NACT OVERALL TIME (G16_TD.py): 0.0 s
	GET_APPROX NACR OVERALL TIME (G16_TD.py): 0.0 s
Total QM took 19.11 s.
Electronic propagation took 0.00 s.
Electronic Norm.: 1.0 1.0
Total MD Step took 19.11 s.
Submitting electronic structure for step 2.
	Gaussian terminated normally in 3.00 s. (GS_NEW)
	Gaussian terminated normally in 10.95 s. (TD_NEW_S1)
	Gaussian terminated normally in 1.80 s. (DIMER)
	TD_NEW_S1 DIM_INFO TIME: 0.05 s
	TD_NEW_S1 GET_MO TIME: 0.01 s
	TD_NEW_S1 GET_CI TIME: 0.01 s
	TD_OLD_S1 GET_MO TIME: 0.01 s
	TD_OLD_S1 GET_CI TIME: 0.01 s
	DIMER TIME: 1.96 s
	G16_NAC.py main: 2.05 s
	gau_nac() (G16_NAC.py): 0.0 s
gau_nac prepare 0.01 s
WFN OVERLAP FORTRAN TOOK 0.09 s
gau_nac dump 0.86 s
	gau_nac.worker() (G16_NAC.py): 0.96 s
	GET_OVERLAP OVERALL TIME (G16_TD.py): 3.04 s
	GET_NACT OVERALL TIME (G16_TD.py): 0.0 s
	GET_APPROX NACR OVERALL TIME (G16_TD.py): 0.0 s
Total QM took 19.13 s.
Electronic propagation took 0.00 s.
Electronic Norm.: 1.0 1.0
Total MD Step took 19.14 s.
