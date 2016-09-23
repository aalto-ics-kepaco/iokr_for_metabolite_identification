%% Make 'cand' from cell-array to struct-array

% tic; 
% inputDir = '/scratch/cs/kepaco/bache1/data/metabolite-identification/GNPS';
% cand = load (strcat (inputDir, '/candidates/GNPS_cand.mat'));
% cand = cand.cand;
% disp ('cand loaded')
% toc; 

% cand_struct = struct ([]);
% 
% tic; 
% for i = 1:size (cand, 1)
%     disp (i)
%     cand_struct(i, 1).data = cand{i}.fp;
%     cand_struct(i, 1).id = cand{i}.inchi;
%     cand_struct(i, 1).num = cand{i}.num;
% end % for
% disp ('put into structure')
% toc;
% 
% save (strcat (inputDir, '/candidates/GNPS_cand_as_struct.mat'), 'cand_struct', '-v7.3')

% for i = 1:size (cand, 1)
%     disp (i)
%     assert (isequal (cand_struct(i, 1).data, cand{i}.fp));
%     assert (isequal (cand_struct(i, 1).id, cand{i}.inchi));
%     assert (isequal (cand_struct(i, 1).num, cand{i}.num));
% end % for

for i = 1:size (cand, 1)
    disp(i)
    cand(i, 1).data = cand(i, 1).data';
    cand(i, 1).id = cand(i, 1).id';
end % for

cand_struct = cand;
save (strcat (inputDir, '/candidates/GNPS_cand_as_struct_transp.mat'), 'cand_struct', '-v7.3')