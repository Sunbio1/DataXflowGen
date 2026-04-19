clear; clc;
clear functions

root = pwd;
pattern = 'r_*';

nStarts = 50;
L1penalty = 50;

outCsv = fullfile(root, 'batch_compare_noL1_vs_L1refit_allLHS.csv');

D = dir(fullfile(root, pattern));
D = D([D.isdir]);
names = sort({D.name});

T = table();

protocols = { ...
    struct('name','NO_L1',    'useL1',false,'refitWithoutL1',false,'freezeZeroAfterL1',false,'zeroTol',1e-6), ...
    struct('name','L1_REFIT', 'useL1',true, 'refitWithoutL1',true, 'freezeZeroAfterL1',true, 'zeroTol',1e-6) ...
};

for i = 1:numel(names)
    runDir = fullfile(root, names{i});

    for pIdx = 1:numel(protocols)
        P = protocols{pIdx};
        seed = 1000 + 10*i + pIdx;

        R = fit_one_runfolder_finalfit_struct(runDir, ...
            'mode','random', ...
            'nStarts', nStarts, ...
            'L1penalty', L1penalty, ...
            'seed', seed, ...
            'useL1', P.useL1, ...
            'refitWithoutL1', P.refitWithoutL1, ...
            'freezeZeroAfterL1', P.freezeZeroAfterL1, ...
            'zeroTol', P.zeroTol, ...
            'forceStartValue', true, ...
            'startValue', 0.2);

        R.protocol = string(P.name);
        T = [T; struct2table(R)]; 
        writetable(T, outCsv);
    end
end

disp(['DONE: ', outCsv]);
