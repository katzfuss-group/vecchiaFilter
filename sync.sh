#!/bin/bash

rsync -av --rsync-path="rsync" --log-file="/home/marcin/HVLF/rsync.log" --exclude=simulations-lorenz/Lorenz* ~/HVLF/ marcin.jurek@z10.stat.tamu.edu:/home/grad/marcin.jurek/HVLF
