#!/bin/sh
#call the programm: ./LB2D_h10_atmos_entropic[nx, Re, v_max, h(obstacle), H/h, width(obstacle), dist(obstacle), pertubation,read_state, entropic, input_file]
#export OMP_NUM_THREADS=16

export OMP_NUM_THREADS=48
#== run kbc jobs ==
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 5 40 5 12 0 1 0 states/h10_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 5 40 5 25 0 1 0 states/h10_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 5 40 5 38 0 1 0 states/h10_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 5 40 5 50 0 1 0 states/h10_kbc_atmos

bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 10 20 10 12 0 1 0 states/h10_kbc_atmos
bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 10 20 10 25 0 1 0 states/h10_kbc_atmos
bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 10 20 10 38 0 1 0 states/h10_kbc_atmos
bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 10 20 10 50 0 1 0 states/h10_kbc_atmos

#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 20 20 20 12 0 1 0 states/h20_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 20 20 20 25 0 1 0 states/h20_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 20 20 20 38 0 1 0 states/h20_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 20 20 20 50 0 1 0 states/h20_kbc_atmos

#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 30 20 30 12 0 1 0 states/h30_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 30 20 30 25 0 1 0 states/h30_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 30 20 30 38 0 1 0 states/h30_kbc_atmos
#bsub -n 48 -W 08:00 ./LB2D_kbc_obstac 100 80000 0.08 30 20 30 50 0 1 0 states/h30_kbc_atmos

#== run entropic jobs ==
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 5 40 5 12 0 1 1 states/h10_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 5 40 5 25 0 1 1 states/h10_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 5 40 5 38 0 1 1 states/h10_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 5 40 5 50 0 1 1 states/h10_entropic_atmos

#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 10 20 10 12 0 1 1 states/h10_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 10 20 10 25 0 1 1 states/h10_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 10 20 10 38 0 1 1 states/h10_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 10 20 10 50 0 1 1 states/h10_entropic_atmos

#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 20 20 20 12 0 1 1 states/h20_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 20 20 20 25 0 1 1 states/h20_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 20 20 20 38 0 1 1 states/h20_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 20 20 20 50 0 1 1 states/h20_entropic_atmos

#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 30 20 30 12 0 1 1 states/h30_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 30 20 30 25 0 1 1 states/h30_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 30 20 30 38 0 1 1 states/h30_entropic_atmos
#bsub -n 48 -W 08:00 ./LB2D_entropic_obstac 100 80000 0.08 30 20 30 50 0 1 1 states/h30_entropic_atmos
