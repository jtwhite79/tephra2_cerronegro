ptf ~
VENT_EASTING 532596.00
VENT_NORTHING 1381862.00
VENT_ELEVATION 1
#
# Note: UTM coordinates are used (add 10,000,000 m in 
#      northern hemisphere
#
PLUME_HEIGHT ~PLUME_HEA           ~ 15000
ALPHA ~ALPHAA              ~ 100.0
BETA ~BETAA               ~ 100.0
ERUPTION_MASS ~ERUPTIONA           ~ 5.0e11
MAX_GRAINSIZE -7  -7
MIN_GRAINSIZE 7    7
MEDIAN_GRAINSIZE ~MEDIAN_GA           ~ 3.0
STD_GRAINSIZE ~STD_GRAIA           ~ 2.0

/*eddy diff for small particles in m2/s (400 cm2/s) */
EDDY_CONST ~EDDY_CONA           ~ 0.009

# diffusion coeff for large particles (m2/s)
DIFFUSION_COEFFICIENT ~DIFFUSIOA           ~ 10000

# threshold for change in diffusion (seconds fall time)
FALL_TIME_THRESHOLD ~FALL_TIMA           ~ 5000

# density model for the pyroclasts
LITHIC_DENSITY 	~  LITHDENA    ~  2600.0  
PUMICE_DENSITY 	~  PUMDENA    ~  1000.0

#define column integration steps
COL_STEPS 400  75
PART_STEPS 140  75

# Note: 
# 0 = uniform distribution using threshold at PLUME_RATIO (no longer used)
# 1 = log-normal distribution using beta (no longer used)
# 2 = beta distribution using parameters alpha and beta (set below)
PLUME_MODEL 2
