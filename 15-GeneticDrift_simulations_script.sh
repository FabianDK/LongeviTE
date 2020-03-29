#Genetic drift simulations of TE frequency change using parameters used by experimental evolution studies
#TE frequencies from Kofler et al. 2015, PLoS Gen.

#Required scripts
simulation_Rscript_not_round.R
simulation_Fabian_Rscript_not_round.R

#IMPORTANT: Open the two requried scripts and adjust folder names


#Rscript for running simulations
popsizeC = as.integer(args[1]) #population size of controls, e.g. 640 in Remolina2012
popsizeS = as.integer(args[2]) #population size of selected, e.g. 320 in Remolina2012
genC = as.integer(args[3]) #generations of controls, e.g. 80 in Remolina2012
genS = as.integer(args[4]) #generations of selected populations, e.g. 50 in Remolina2012
reps = as.integer(args[5]) #number of replicate populations, e.g. 3 in Remolina2012
direction.cutoff = as.integer(args[6]) #cut-off to define S>C and C>S TEs, e.g. 0.5: define a TE family that is of larger frequency in S than C populations in >50% of the simulations to be S>C.
run.numb = as.integer(args[7]) #number of simulations: e.g. 100
study = as.character(args[8]) #Only select TEs that were also detected in each study. Options: Remolina, Hoedjes, Carnes or choose. For Fabian, a separate script was used 

#If study is set to choose, two other options are available:
setfreq = as.character(args[9]) #Set frequency between 0 and 1.
ins = as.integer(args[10]) #Set the number of TE insertions
#This was done as a control step to check if the number of TE insertions makes a difference on whether a TE is classified as S>C or C>S.


############################
#Remolina: upper boundary of population size
Rscript simulation_Rscript_not_round.R 640 320 80 50 3 "0.5" 5000 remolina
[1] "Average Prop. S>C TEs: 0.474 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.477 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.048 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 4.281 || sd = 2.04"
[1] "Average TEs Consistent C>S: 93.049 || sd = 2.54"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0.804"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0188"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "46.58% of sims with more S>C than C>S"
[1] "48.96% of sims with more C>S than S>C"
[1] "4.46% of sims with Equal numbers of S>C and C>S"

#Remolina: lower boundary of population size
Rscript simulation_Rscript_not_round.R 440 220 80 50 3 "0.5" 10000 remolina
[1] "Average Prop. S>C TEs: 0.467 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.473 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.06 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 3.99 || sd = 1.95"
[1] "Average TEs Consistent C>S: 91.033 || sd = 2.68"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0.7687"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0123"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "45.73% of sims with more S>C than C>S"
[1] "50.29% of sims with more C>S than S>C"
[1] "3.98% of sims with Equal numbers of S>C and C>S"

#Remolina: upper boundary of population size, with 50% bottleneck in S lines
Rscript simulation_Rscript_not_round.R 640 160 80 50 3 "0.5" 5000 remolina
[1] "Average Prop. S>C TEs: 0.459 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.486 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.055 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 3.9 || sd = 1.93"
[1] "Average TEs Consistent C>S: 91.47 || sd = 2.44"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0.757"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0072"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "36.94% of sims with more S>C than C>S"
[1] "58.86% of sims with more C>S than S>C"
[1] "4.2% of sims with Equal numbers of S>C and C>S"

#Remolina: upper boundary of population size, with 75% bottleneck in S lines 
Rscript simulation_Rscript_not_round.R 640 80 80 50 3 "0.5" 5000 remolina
[1] "Average Prop. S>C TEs: 0.438 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.502 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.06 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 4.417 || sd = 1.98"
[1] "Average TEs Consistent C>S: 89.299 || sd = 2.41"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0.8372"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0024"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "22.2% of sims with more S>C than C>S"
[1] "74.66% of sims with more C>S than S>C"
[1] "3.14% of sims with Equal numbers of S>C and C>S"

#Remolina: lower boundary of population size, with 50% bottleneck in S lines
Rscript simulation_Rscript_not_round.R 440 110 80 50 3 "0.5" 5000 remolina
[1] "Average Prop. S>C TEs: 0.445 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.485 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.071 || sd = 0.02"
[1] ""
[1] "Average TEs Consistent S>C: 3.433 || sd = 1.84"
[1] "Average TEs Consistent C>S: 89.347 || sd = 2.62"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0.6724"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0034"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "31.9% of sims with more S>C than C>S"
[1] "64.68% of sims with more C>S than S>C"
[1] "3.42% of sims with Equal numbers of S>C and C>S"

#Remolina: lower boundary of population size, with 75% bottleneck in S lines
Rscript simulation_Rscript_not_round.R 440 55 80 50 3 "0.5" 5000 remolina
[1] "Average Prop. S>C TEs: 0.414 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.508 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.078 || sd = 0.02"
[1] ""
[1] "Average TEs Consistent S>C: 4.004 || sd = 1.92"
[1] "Average TEs Consistent C>S: 86.991 || sd = 2.62"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0.7732"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 4e-04"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "14.36% of sims with more S>C than C>S"
[1] "83.32% of sims with more C>S than S>C"
[1] "2.32% of sims with Equal numbers of S>C and C>S"

############################
#Carnes: upper boundary of population size
Rscript simulation_Rscript_not_round.R 2000 2000 850 170 5 "0.5" 5000 carnes
[1] "Average Prop. S>C TEs: 0.493 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.473 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.034 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 3.432 || sd = 1.81"
[1] "Average TEs Consistent C>S: 96.732 || sd = 2.46"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0.5086"
[1] "56.66% of sims with more S>C than C>S"
[1] "39.66% of sims with more C>S than S>C"
[1] "3.68% of sims with Equal numbers of S>C and C>S" 

#Carnes: lower boundary of population size
Rscript simulation_Rscript_not_round.R 1000 1000 850 170 5 "0.5" 5000 carnes
[1] "Average Prop. S>C TEs: 0.497 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.452 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.052 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 3.348 || sd = 1.77"
[1] "Average TEs Consistent C>S: 86.635 || sd = 3.39"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0.9298"
[1] "66.54% of sims with more S>C than C>S"
[1] "30.28% of sims with more C>S than S>C"
[1] "3.18% of sims with Equal numbers of S>C and C>S"

#Carnes: lower boundary of population size, with 50% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 1000 500 850 170 5 "0.5" 5000 carnes
[1] "Average Prop. S>C TEs: 0.485 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.442 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.073 || sd = 0.02"
[1] ""
[1] "Average TEs Consistent S>C: 0.777 || sd = 0.87"
[1] "Average TEs Consistent C>S: 84.921 || sd = 3.45"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0.9962"
[1] "65.3% of sims with more S>C than C>S"
[1] "30.72% of sims with more C>S than S>C"
[1] "3.98% of sims with Equal numbers of S>C and C>S"

#Carnes: upper boundary of population size, with 50% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 2000 1000 850 170 5 "0.5" 5000 carnes
[1] "Average Prop. S>C TEs: 0.487 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.466 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.047 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 1.196 || sd = 1.07"
[1] "Average TEs Consistent C>S: 95.573 || sd = 2.42"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0.8696"
[1] "56.94% of sims with more S>C than C>S"
[1] "39.14% of sims with more C>S than S>C"
[1] "3.92% of sims with Equal numbers of S>C and C>S"

#Carnes: lower boundary of population size, with 75% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 1000 250 850 170 5 "0.5" 5000 carnes
[1] "Average Prop. S>C TEs: 0.455 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.441 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.104 || sd = 0.02"
[1] ""
[1] "Average TEs Consistent S>C: 0.213 || sd = 0.46"
[1] "Average TEs Consistent C>S: 82.737 || sd = 3.46"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 1"
[1] "52.98% of sims with more S>C than C>S"
[1] "42.96% of sims with more C>S than S>C"
[1] "4.06% of sims with Equal numbers of S>C and C>S"

#Carnes: upper boundary of population size, with 75% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 2000 500 850 170 5 "0.5" 5000 carnes
[1] "Average Prop. S>C TEs: 0.472 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.466 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.061 || sd = 0.01"
[1] ""
[1] "Average TEs Consistent S>C: 0.337 || sd = 0.58"
[1] "Average TEs Consistent C>S: 93.794 || sd = 2.39"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0.9852"
[1] "50.82% of sims with more S>C than C>S"
[1] "45.34% of sims with more C>S than S>C"
[1] "3.84% of sims with Equal numbers of S>C and C>S"


############################
#HOEDJES
#4 populations per 3 diet regimes - ignored diet regime!

#Hoedjes: lower boundary of population size
Rscript simulation_Rscript_not_round.R 2000 2000 115 58 12 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.492 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.49 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 0 || sd = 0"
[1] "Average TEs Consistent C>S: 106.984 || sd = 0.12"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0258"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "50.54% of sims with more S>C than C>S"
[1] "49.22% of sims with more C>S than S>C"
[1] "0.24% of sims with Equal numbers of S>C and C>S"

#Hoedjes: upper boundary of population size
Rscript simulation_Rscript_not_round.R 4000 4000 115 58 12 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.492 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.489 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 0 || sd = 0"
[1] "Average TEs Consistent C>S: 107 || sd = 0.02"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.027"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "50.68% of sims with more S>C than C>S"
[1] "49.22% of sims with more C>S than S>C"
[1] "0.1% of sims with Equal numbers of S>C and C>S"

#Hoedjes: lower boundary of population size, with 50% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 2000 1000 115 58 12 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.492 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.489 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 0 || sd = 0"
[1] "Average TEs Consistent C>S: 107 || sd = 0.02"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.027"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "50.68% of sims with more S>C than C>S"
[1] "49.22% of sims with more C>S than S>C"
[1] "0.1% of sims with Equal numbers of S>C and C>S"

#Hoedjes: upper boundary of population size, with 50% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 4000 2000 115 58 12 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.491 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.49 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 0 || sd = 0.01"
[1] "Average TEs Consistent C>S: 106.985 || sd = 0.12"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.025"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "50.36% of sims with more S>C than C>S"
[1] "49.42% of sims with more C>S than S>C"
[1] "0.22% of sims with Equal numbers of S>C and C>S"


#Hoedjes: lower boundary of population size, with 75% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 2000 500 115 58 12 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.49 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.491 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 0 || sd = 0.01"
[1] "Average TEs Consistent C>S: 106.104 || sd = 0.83"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0232"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "49% of sims with more S>C than C>S"
[1] "50.18% of sims with more C>S than S>C"
[1] "0.82% of sims with Equal numbers of S>C and C>S"

#Hoedjes: upper boundary of population size, with 75% bottleneck in S lines:
Rscript simulation_Rscript_not_round.R 4000 1000 115 58 12 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.489 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.492 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 0.001 || sd = 0.02"
[1] "Average TEs Consistent C>S: 106.795 || sd = 0.43"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.026"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "48.64% of sims with more S>C than C>S"
[1] "51.24% of sims with more C>S than S>C"
[1] "0.12% of sims with Equal numbers of S>C and C>S"

#Hoedjes: lower boundary of population size, with 4 replicates for consistent differences:
Rscript simulation_Rscript_not_round.R 2000 2000 115 58 4 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.491 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.488 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.021 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 1.784 || sd = 1.31"
[1] "Average TEs Consistent C>S: 104.668 || sd = 1.47"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0256"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "50.34% of sims with more S>C than C>S"
[1] "48.02% of sims with more C>S than S>C"
[1] "1.64% of sims with Equal numbers of S>C and C>S"

#Hoedjes: upper boundary of population size, with 4 replicates for consistent differences:
Rscript simulation_Rscript_not_round.R 4000 4000 115 58 4 "0.5" 5000 hoedjes
[1] "Average Prop. S>C TEs: 0.492 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.489 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.019 || sd = 0"
[1] ""
[1] "Average TEs Consistent S>C: 1.703 || sd = 1.29"
[1] "Average TEs Consistent C>S: 105.19 || sd = 1.33"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.0288"
[1] "Prop. of sims with larger/equal C>S than observed: 1"
[1] "Prop. of sims with larger/equal Equal than observed: 0"
[1] "51.26% of sims with more S>C than C>S"
[1] "48.38% of sims with more C>S than S>C"
[1] "0.36% of sims with Equal numbers of S>C and C>S"

############################
#For simulations of Fabian2018 data, a different script was used (this is because the setup is slightly different with 2 controls and 4 selected populations. Also the number of generations vary between all populations).

#Fabian: Generations and reps fixed as different within/between Regimes
popsizeC = as.integer(args[1]) #population size of C, i.e. 300
popsizeS = as.integer(args[2]) #population size of S, i.e. 300
run.numb = as.integer(args[3]) #Number of simulations, e.g. 5000

#Fabian2018 with N = 300:
Rscript simulation_Fabian_Rscript_not_round.R 300 300 "0.5" 5000
[1] "Average Prop. S>C TEs: 0.494 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.37 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.135 || sd = 0.03"
[1] ""
[1] "Average TEs Consistent S>C: 4.591 || sd = 1.99"
[1] "Average TEs Consistent C>S: 59.248 || sd = 4.27"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.8232"
[1] "Prop. of sims with larger/equal C>S than observed: 0.874"
[1] "Prop. of sims with larger/equal Equal than observed: 0.0018"
[1] "91.6% of sims with more S>C than C>S"
[1] "6.38% of sims with more C>S than S>C"
[1] "2.02% of sims with Equal numbers of S>C and C>S"


#Fabian2018, with 50% bottleneck in S populations:
Rscript simulation_Fabian_Rscript_not_round.R 300 150 "0.5" 5000
[1] "STARTING SIMS FOR PROPORTION OF S>C AND C>S"
[1] "Average Prop. S>C TEs: 0.429 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.356 || sd = 0.05"
[1] "Average Prop. Equal TEs: 0.215 || sd = 0.03"
[1] ""
[1] "Average TEs Consistent S>C: 2.246 || sd = 1.47"
[1] "Average TEs Consistent C>S: 57.456 || sd = 4.21"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.3076"
[1] "Prop. of sims with larger/equal C>S than observed: 0.7912"
[1] "Prop. of sims with larger/equal Equal than observed: 0.5"
[1] "78.28% of sims with more S>C than C>S"
[1] "18.84% of sims with more C>S than S>C"
[1] "2.88% of sims with Equal numbers of S>C and C>S"

#Fabian2018, with 75% bottleneck in S populations:
Rscript simulation_Fabian_Rscript_not_round.R 300 75 "0.5" 5000
[1] "Average Prop. S>C TEs: 0.361 || sd = 0.05"
[1] "Average Prop. C>S TEs: 0.346 || sd = 0.04"
[1] "Average Prop. Equal TEs: 0.293 || sd = 0.04"
[1] ""
[1] "Average TEs Consistent S>C: 1.57 || sd = 1.23"
[1] "Average TEs Consistent C>S: 55.065 || sd = 4.23"
[1] ""
[1] "Prop of Sims with more or equal Consistent S>C: 0"
[1] "Prop of Sims with more or equal Consistent C>S: 1"
[1] "Prop of Sims with Consistent S>C more than consistent C>S: 0"
[1] ""
[1] "Prop. of sims with larger/equal S>C than observed: 0.021"
[1] "Prop. of sims with larger/equal C>S than observed: 0.7228"
[1] "Prop. of sims with larger/equal Equal than observed: 0.9898"
[1] "55.74% of sims with more S>C than C>S"
[1] "39.9% of sims with more C>S than S>C"
[1] "4.36% of sims with Equal numbers of S>C and C>S"


