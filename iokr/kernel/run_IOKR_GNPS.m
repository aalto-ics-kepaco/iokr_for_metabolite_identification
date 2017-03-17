
addpath('general_functions');
addpath('preprocessing')

%--------------------------------------------------------------
% Load Data
%--------------------------------------------------------------

inchi = readtext('../data/inchi.txt'); % Inchi
mf_corres = load('../data/matching_mf_train.txt');  % correspondance with the set of unique mf
load ../data/fp.mat Y; % fingerprints
Y = full(Y);
[~,n] = size(Y);

% Candidates description
load ../data/GNPS_cand.mat cand;

% Input kernels
load ../data/input_kernels/kernels_computed_by_myself/KX_list.mat KX_list;

% Indices of test examples used in the evaluation
eval = load('../data/ind_eval.txt'); 

%--------------------------------------------------------------
% Cross-validation
%--------------------------------------------------------------

iokr_param = struct('center',1,'mkl','unimkl');
select_param = struct('cv_type','loocv','lambda',[1e-5 1e-4 1e-3 1e-2 5e-2 1e-1 5e-1 1 10 100]);
ky_param = struct('type','gaussian','base_kernel','tanimoto','param_selection','entropy');
output_param = struct('representation','kernel','kernel_param',ky_param);

rank = zeros(n,1);
cand_num = zeros(n,1); % vector containing the number of candidates for each test example

n_folds = 10; % number of folds
ind_fold = load('../data/cv_ind.txt'); % indices of the different folds

for i = 1:n_folds
    
    % Create training and test sets
    test_set = find(ind_fold == i);
    train_set = setdiff(1:n,test_set);
    test_set = intersect(test_set,eval);

    % training
    KX_list_train = cellfun(@(x) x(train_set,train_set), KX_list, 'UniformOutput', false);
    Y_train = Y(:,train_set);
    
    train_model = Train_IOKR(KX_list_train, Y_train, output_param, select_param, iokr_param);
    
    % Prediction
        
    KX_list_train_test = cellfun(@(x) x(train_set,test_set), KX_list, 'UniformOutput', false);
    KX_list_test = cellfun(@(x) x(test_set,test_set), KX_list, 'UniformOutput', false);
    %Y_C_test = cellfun(@(x) full(x.fp), cand(mf_corres(test_set)),'UniformOutput',false);
    Y_C_test = cellfun(@(x) full(x.fp)', cand(mf_corres(test_set(1:10))),'UniformOutput',false);
    
    score = Test_IOKR(KX_list_train_test, KX_list_test, train_model, Y_train, Y_C_test, iokr_param.center);
        
end
