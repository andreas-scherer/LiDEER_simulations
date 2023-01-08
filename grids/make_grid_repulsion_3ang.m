clear all;
close all;
clc;

npoints=50;

repulsion(npoints,2,1000) ,


[alphas,betas,gammas,weights]=repulsion(npoints,2,1000);



save(['rep_1ang_', num2str(npoints), 'pts_hem.mat'],'alphas','betas','gammas','weights');