stat_rand_0_05 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/d554058970f57b9c04ac3f4683df1d68.mat');
stat_rand_0_1 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/a1abdba9ceba03b7669685f628719b13.mat');
stat_rand_0_5 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/db07946d36a15c1f47ce2f98f9d48e87.mat');
stat_rand_1 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/7fd6e73f324588699b2c1cde1f550ffe.mat');
stat_rand_5 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/671f61f0e9a451912284ffc3ab3636a7.mat');
stat_rand_50 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/80c9cdfda75c4ecd0cb00d2c54856194.mat');
stat_rand_10 = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/99b65f2926d90d8522769902ed793a03.mat');

stat_all_excl = matfile ('/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS/pre_calculated_stats/e83af7da3fb7b4d899cd401a386b648f.mat');

s_rand_0_05__ = stat_rand_0_05.stats(1,1);
s_rand_0_1__ = stat_rand_0_1.stats(1,1);
s_rand_0_5__ = stat_rand_0_5.stats(1,1);
s_rand_1__ = stat_rand_1.stats(1,1);
s_rand_5__ = stat_rand_5.stats(1,1);
s_rand_10__ = stat_rand_10.stats(1,1);
s_rand_50__ = stat_rand_50.stats(1,1);

s_all__ = stat_all_excl.stats(1,1);

n_elements = numel (s_all__.Cov_Psi_C_train);
mse = [norm(s_all__.Cov_Psi_C_train - s_rand_0_05__.Cov_Psi_C_train,'fro')^2, ...
       norm(s_all__.Cov_Psi_C_train - s_rand_0_1__.Cov_Psi_C_train,'fro')^2, ...
       norm(s_all__.Cov_Psi_C_train - s_rand_0_5__.Cov_Psi_C_train,'fro')^2, ...
       norm(s_all__.Cov_Psi_C_train - s_rand_1__.Cov_Psi_C_train,'fro')^2, ...
       norm(s_all__.Cov_Psi_C_train - s_rand_5__.Cov_Psi_C_train,'fro')^2, ...
       norm(s_all__.Cov_Psi_C_train - s_rand_10__.Cov_Psi_C_train,'fro')^2, ...
       norm(s_all__.Cov_Psi_C_train - s_rand_50__.Cov_Psi_C_train,'fro')^2] / n_elements;
   
close all
subplot (1, 2, 1)
plot ([0.05, 0.1, 0.5, 1, 5, 10, 50], mse, '*-');
grid()
title ('MSE: Approximation of the input-feature covariance matrix.')
xlabel ('% of random candidates');
ylabel ('MSE');

subplot (1, 2, 2)
loglog ([0.05, 0.1, 0.5, 1, 5, 10, 50], mse, '*-');
grid()
title ('MSE: Approximation of the input-feature covariance matrix.')
xlabel ('% of random candidates');
ylabel ('MSE');