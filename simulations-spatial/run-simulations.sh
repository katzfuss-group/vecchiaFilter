Rscript generate_scenarios.R gauss &> gauss.log
Rscript generate_scenarios.R poisson &> poisson.log
Rscript generate_scenarios.R logistic &> logistic.log
Rscript generate_scenarios.R gamma &> gamma.log
rm results
cat gauss.results > results
tail -n +2 poisson.results >> results
tail -n +2 logistic.results >> results
tail -n +2 gamma.results >> results

python generate-plot.py
