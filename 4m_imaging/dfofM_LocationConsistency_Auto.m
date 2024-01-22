%%

clear all; close all; clc

base = 'D:\Mice';
cd(base)

folds = findSubF('pcaica',3,[],0);


%%

for ff = 1:length(folds)
    cd(folds{ff})

    try
        copyfile([base '\dfofM_LocationConsistency.m'],'dfofM_LocationConsistency.m');
        dfofM_LocationConsistency
    catch err
        disp(err);
    end
end

cd(base)

