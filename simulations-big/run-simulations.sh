#!/bin/bash

Rscript sim-adv-diff-big.r gauss &> gauss.log &
Rscript sim-adv-diff-big.r poisson &> poisson.log &
Rscript sim-adv-diff-big.r logistic &> logistic.log &
Rscript sim-adv-diff-big.r gamma &> gamma.log &

python plot-scores.py
