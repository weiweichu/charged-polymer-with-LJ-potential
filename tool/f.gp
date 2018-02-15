set term x11 enhanced font "Hevitica, 20"

set xlabel "psi" font "Hevitica, 20"
set ylabel "free energy (128 x 32)" font "Hevitica, 20"

w = 1
plot [0:][-0.1:10] \
   '0/doshistory/weight.19' w l lw w ti "chiN = 11",\
   '1/doshistory/weight.18' w l lw w ti "       12",\
   '2/doshistory/weight.18' w l lw w ti "       13",\
   '3/doshistory/weight.18' w l lw w ti "       14",\
   '4/doshistory/weight.18' w l lw w ti "       15",\
   '5/doshistory/weight.18' w l lw w ti "       16",\
   '6/doshistory/weight.18' w l lw w ti "       17",\
   '7/doshistory/weight.18' w l lw w ti "       18"
