-R package needed:
	# truncnorm, doParallel,mgcv
	#install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable", dep=TRUE)
	#install.packages() install these package first
#install r package with conda on unix

-input files:
	"initial_par_input.csv"
	"initial_parrange_input.csv"
	"data.csv"


-output files;
	all_parm_cl.csv #parameter and confidence level
	range_fit.png #the fitted data
	outfile="optim_parameter"---including some intermediate results files.



-Main.R is the main function 
	using sbatch run Main.R on HPC (see batch file example, example.sh)
	$ sbatch example.sh #to run model

--setting in Main.R 
First part:
files related with re-sample the observation data using a gam model implemented in INLA for posterior predictive analyses, initialized with warm starts.
	--("simulate.R")
	--("lme_inla_formula.R")
	--("proc.formula.R")
	--("getVars.R")
	--("s2inla.R")
Second part:
find the optimal parameters(initial parameters) for dammm.
	--optimizeP.R
		find the optimal parameters(initial parameters) for dammm with R optim function.
	{
	## start the iterations with optimizing O2 while fixing CH4,N2O
	## continue iteration with CH4,N2O  while fixing O2
	## repeat iteration with O2 while fixing CH4.N2O
	
	optim() try to find the best parameter that minimize residual sum of squares.
	the minimize residual sum of squares caculate from ---lulik.rh.optim, lulik.ch4.n2o.optim
	}

	--lulik_utils.R 
		utility function to setting parameters run lulik() and the conventional Nelder-Mead algorithm using optim

	--lulik.R
 		log likelihood function given observation data.

	--damm.R
 		damm {
			 First creat numbers of microsite at each obervation 
			 then calling damm_work at each microsite
			 last combine the average over one obervation with all microsites.
		 }
	--damm_work
		calculate each microsites with given parameters with Arrhenius and Michaelis Menten equstions

	--soilm_pdf_3.R
		build microsite bivariate distribution
	--damm_utils.R 
		log tranfer parameters and back the logged values to normal


