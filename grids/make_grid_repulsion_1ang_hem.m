clear all;
close all;
clc;

npoints=400;


betas=linspace(0,pi/2,npoints).';
alphas=0*betas;
gammas=0*betas;
weights=sin(betas)/sum(sin(betas));


save(['rep_1ang_', num2str(npoints), 'pts_hem.mat'],'alphas','betas','gammas','weights');