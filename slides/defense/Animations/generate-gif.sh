python /home/marcin/AgProject/plots-code/generateAgPlots.py
cp /home/marcin/AgProject/plots/truth_* ~/MRF/slides/JSM-Summer18/Animations/
cp /home/marcin/AgProject/plots/obs_* ~/MRF/slides/JSM-Summer18/Animations/

rm truth.gif obs.gif
convert -delay 20 -loop 0 truth_* truth.gif
convert -delay 20 -loop 0 obs_* obs.gif
