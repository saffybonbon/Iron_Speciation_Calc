# Iron_Speciation_Calc

Python tool to estimate melt FeO and Fe2O3 from FeOt, with error propagation via monte carlo simulation.

#### Input Spreadsheet:

  - Temperature is given in °K
  - Pressure is given in Pa
  - Oxygen fugacity is given as log10fO2
  
 #### Available models:
 
  - Irvine & Barager (1971)
  - Le Maitre (1976)
  - Sack et al. (1980)
  - Kress & Carmichael (1991): Eqn 6 & 7
  - Jayasuriya et al. (2004)
  - Putirka (2016b)
  - O'Neill et al. (2018)
  - Borisov et al. (2018)
  
### Requirements

  - pandas
  - numpy
  - periodictable
  
  For plotting:
  - matplotlib
  - seaborn
