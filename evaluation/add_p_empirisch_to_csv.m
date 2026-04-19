clear; clc;
clear functions

csvIn  = 'batch_compare_noL1_vs_L1refit_allLHS.csv';
csvOut = 'batch_compare_noL1_vs_L1refit_with_pEmp.csv';

Nboot = 10;                 
L1penalty = 50;
zeroTol = 1e-6;

forceStartValue = true;
startValue = 0.2;

% For random mode: multi-start like in your benchmark
nStarts_random = 50;

% For reproducibility
baseSeedFit  = 1234;
baseSeedBoot = 9000;

T = readtable(csvIn);

% OPTIONAL: compute only r_000 (for testing)
% T = T(contains(string(T.runDir), "/r_000"), :);

% If p_empirisch column is missing / comes in as text, initialize/clean it
if ~ismember('p_empirisch', T.Properties.VariableNames)
    T.p_empirisch = NaN(height(T),1);
end

for r = 1:height(T)
    if T.ok(r) ~= 1
        continue;
    end

    if ~isnan(T.p_empirisch(r))
        continue;
    end

    runDir   = string(T.runDir(r));
    protocol = string(T.protocol(r));

    fprintf('\n=== %d/%d ===\n', r, height(T));
    fprintf('runDir: %s\n', runDir);
    fprintf('protocol: %s\n', protocol);

    % Decide fit settings from protocol
    if protocol == "NO_L1"
        useL1 = false;
        refitWithoutL1 = false;
        freezeZeroAfterL1 = false;
    else
        useL1 = true;
        refitWithoutL1 = true;
        freezeZeroAfterL1 = true;
    end

    % Decide mode/nStarts: r_000 is "AI structure" in your setup, but this is fitting.
    % If you do NOT want to LHS-fit r_000: set nStarts=1 for r_000.
    if contains(runDir, "/r_000")
        mode = "ai";    % no LHS
        nStarts = 1;
    else
        mode = "random";
        nStarts = nStarts_random;
    end

    fitSeed  = baseSeedFit  + r;
    bootSeed = baseSeedBoot + r;

    try
        % 1) load + fit to get FINAL parameters in global ar
        fit_one_row(runDir, mode, nStarts, ...
            useL1, L1penalty, refitWithoutL1, freezeZeroAfterL1, zeroTol, ...
            forceStartValue, startValue, fitSeed);

        % 2) empirical p via parametric bootstrap (refit each replicate with arFit only)
        [p_emp, chi2_boot] = empirical_p_bootstrap(Nboot, bootSeed);

        T.p_empirisch(r) = p_emp;

        fprintf('p_empirisch(N=%d) = %.4g\n', Nboot, p_emp);
        fprintf('boot chi2 (min/med/max) = %.4g / %.4g / %.4g\n', ...
            min(chi2_boot), median(chi2_boot,'omitnan'), max(chi2_boot));

        writetable(T, csvOut);

    catch ME
        fprintf('FAIL row %d: %s\n', r, ME.message);
        % conservative: keep NaN, but still write intermediate results
        writetable(T, csvOut);
    end
end

disp(['DONE: ', csvOut]);

%% ======================= LOCAL FUNCTIONS =======================

function fit_one_row(runDir, mode, nStarts, useL1, L1penalty, refitWithoutL1, freezeZeroAfterL1, zeroTol, forceStartValue, startValue, seed)
    old = pwd;
    cleanup = onCleanup(@() cd(old));
    cd(runDir);

    ensure_d2d_layout(pwd);
    [modelname, dataname] = autodetect_model_data(pwd);

    clear global ar
    global ar
    ar = [];

    rng(seed);

    arInit;
    arLoadModel(modelname);
    arLoadData(dataname, 1);
    arCompileAll;

    % Robustness
    if isfield(ar,'config')
        if isfield(ar.config,'cvodes_maxsteps'), ar.config.cvodes_maxsteps = 200000; end
        if isfield(ar.config,'cvodes_mxstep'),   ar.config.cvodes_mxstep   = 200000; end
        if isfield(ar.config,'mxstep'),          ar.config.mxstep          = 200000; end
        if isfield(ar.config,'mxsteps'),         ar.config.mxsteps         = 200000; end
        if isfield(ar.config,'rtol'),            ar.config.rtol            = 1e-6; end
        if isfield(ar.config,'atol'),            ar.config.atol            = 1e-8; end
    end

    arInitValues;

    % L1: do NOT touch qFit, only set qL1reg
    if ~isfield(ar,'qL1reg') || isempty(ar.qL1reg)
        ar.qL1reg = zeros(size(ar.p));
    end
    ar.qL1reg(:) = 0;

    if useL1
        ar.config.useL1 = 1;
        ar.config.l1penalty = L1penalty;

        target_prefixes = {'a_','b_','h_','delta_'};
        for ii = 1:numel(ar.pLabel)
            pName = ar.pLabel{ii};
            if any(cellfun(@(prefix) startsWith(pName, prefix), target_prefixes))
                ar.qL1reg(ii) = 1;
            end
        end
    else
        ar.config.useL1 = 0;
        ar.config.l1penalty = 0;
    end

    % Enforce start values (all fitted parameters start the same)
    if forceStartValue
        idx = find(ar.qFit(:)==1);
        ar.p(idx) = startValue;
        ar.p = max(min(ar.p, ar.ub), ar.lb);
    end

    % Optimizer
    optimizer_name = 'lsqnonlin';
    idxOpt = find(strcmp(ar.config.optimizers, optimizer_name), 1);
    if isempty(idxOpt), error('Optimizer lsqnonlin not found in ar.config.optimizers'); end
    ar.config.optimizer = idxOpt;

    ar.config.maxsteps = 20000;
    ar.config.optimoptions = optimoptions(optimizer_name, ...
        'MaxIter', 20000, ...
        'MaxFunctionEvaluations', 500000, ...
        'Display', 'off');

    % FIT
    if mode == "ai"
        ar.config.useInitialGuess = 1;
        arFit;                         % exactly 1 start, no LHS
    else
        % Random: arFit from start + additional LHS (keep the best)
        ar.config.useInitialGuess = 1;
        arFit;
        arCalcMerit(true);
        chi2_best = get_data_chi2(ar);
        p_best = ar.p;

        if nStarts > 1
            ar.config.useInitialGuess = 0;
            arFitLHS(nStarts);
            arCalcMerit(true);
            chi2_lhs = get_data_chi2(ar);
            if isfinite(chi2_lhs) && chi2_lhs < chi2_best
                chi2_best = chi2_lhs;
                p_best = ar.p;
            end
            ar.p = p_best;
            arCalcMerit(true);
        end
    end

    % optional refit without L1 (fair "pure" fit)
    if useL1 && refitWithoutL1
        ar.config.useL1 = 0;
        ar.config.l1penalty = 0;

        if freezeZeroAfterL1 && isfield(ar,'qL1reg') && ~isempty(ar.qL1reg)
            idxReg  = (ar.qL1reg(:)==1);
            idxZero = idxReg & (abs(ar.p(:)) < zeroTol);

            ar.p(idxZero)    = 0;
            ar.qFit(idxZero) = 0;
            ar.lb(idxZero)   = 0;
            ar.ub(idxZero)   = 0;
        end

        ar.config.useInitialGuess = 1;
        arFit;
    end

    arCalcMerit(true);
end

function [p_emp, chi2_boot] = empirical_p_bootstrap(Nboot, seed)
    rng(seed);
    global ar

    arCalcMerit(true);
    chi2_obs = get_data_chi2(ar);

    backup = backup_yExp(ar);
    p_backup = ar.p;

    chi2_boot = NaN(Nboot,1);

    for b = 1:Nboot
        ok = false;
        for attempt = 1:3
            try
                ar.p = p_backup;       % always start from the same final fit
                arCalcMerit(true);

                [yPredCells, sigmaCells] = get_pred_and_sigma(ar);
                ar = set_yExp_from_cells(ar, yPredCells, sigmaCells);

                % refit (fast): arFit only, no LHS
                ar.config.useInitialGuess = 1;
                arFit;

                arCalcMerit(true);
                chi2_boot(b) = get_data_chi2(ar);
                ok = true;
                break;
            catch
                ok = false;
            end
        end
        if ~ok
            % conservative: if bootstrap fit fails, count as "worse"
            chi2_boot(b) = Inf;
        end
    end

    ar = restore_yExp(ar, backup);
    ar.p = p_backup;
    arCalcMerit(true);

    % plus-one correction
    p_emp = (1 + sum(chi2_boot >= chi2_obs)) / (Nboot + 1);
end

function chi2 = get_data_chi2(ar)
    if isfield(ar,'chi2fit') && ~isempty(ar.chi2fit) && isfinite(ar.chi2fit)
        chi2 = ar.chi2fit;
    else
        chi2 = ar.chi2;
    end
end

function backup = backup_yExp(ar)
    backup = struct();
    backup.yExp = cell(numel(ar.model),1);
    for m = 1:numel(ar.model)
        backup.yExp{m} = cell(numel(ar.model(m).data),1);
        for d = 1:numel(ar.model(m).data)
            backup.yExp{m}{d} = ar.model(m).data(d).yExp;
        end
    end
end

function ar = restore_yExp(ar, backup)
    for m = 1:numel(ar.model)
        for d = 1:numel(ar.model(m).data)
            ar.model(m).data(d).yExp = backup.yExp{m}{d};
        end
    end
end

function [yPredCells, sigmaCells] = get_pred_and_sigma(ar)
    yPredCells = cell(numel(ar.model),1);
    sigmaCells = cell(numel(ar.model),1);

    for m = 1:numel(ar.model)
        yPredCells{m} = cell(numel(ar.model(m).data),1);
        sigmaCells{m} = cell(numel(ar.model(m).data),1);
        for d = 1:numel(ar.model(m).data)
            dd = ar.model(m).data(d);

            if ~isfield(dd,'y')
                error('Missing prediction field ar.model(%d).data(%d).y', m, d);
            end
            ypred = dd.y;

            % sigma heuristics (depends on D2D version)
            if isfield(dd,'yExpStd') && ~isempty(dd.yExpStd)
                sig = dd.yExpStd;
            elseif isfield(dd,'ystd') && ~isempty(dd.ystd)
                sig = dd.ystd;
            elseif isfield(dd,'yStd') && ~isempty(dd.yStd)
                sig = dd.yStd;
            elseif isfield(dd,'yExp')
                sig = ones(size(dd.yExp));   % smoke-test fallback
            else
                sig = ones(size(ypred));
            end

            sig(~isfinite(sig)) = 1;
            sig(sig<=0) = 1;

            yPredCells{m}{d} = ypred;
            sigmaCells{m}{d} = sig;
        end
    end
end

function ar = set_yExp_from_cells(ar, yPredCells, sigmaCells)
    for m = 1:numel(ar.model)
        for d = 1:numel(ar.model(m).data)
            ypred = yPredCells{m}{d};
            sig   = sigmaCells{m}{d};

            yboot = ypred + randn(size(ypred)).*sig;
            ar.model(m).data(d).yExp = yboot;
        end
    end
end

function ensure_d2d_layout(runDir)
    cd(runDir);

    if ~exist('Models','dir') && exist('model','dir')
        if isunix
            system('rm -f Models; ln -s model Models');
        else
            [st,~] = system('cmd /c mklink /J Models model');
            if st ~= 0
                mkdir('Models'); copyfile(fullfile('model','*'), 'Models');
            end
        end
    end

    if ~exist('Data','dir') && exist('data','dir')
        if isunix
            system('rm -f Data; ln -s data Data');
        else
            [st,~] = system('cmd /c mklink /J Data data');
            if st ~= 0
                mkdir('Data'); copyfile(fullfile('data','*'), 'Data');
            end
        end
    end

    assert(exist('Models','dir')==7, 'Neither Models nor model found.');
    assert(exist('Data','dir')==7, 'Neither Data nor data found.');
end

function [modelname, dataname] = autodetect_model_data(runDir)
    modelsDir = fullfile(runDir,'Models');
    dataDir   = fullfile(runDir,'Data');

    mdef = dir(fullfile(modelsDir,'*.def'));
    assert(~isempty(mdef), 'No .def found in Models folder.');
    [~, modelname] = fileparts(mdef(1).name);

    ddef = dir(fullfile(dataDir,'*.def'));
    assert(~isempty(ddef), 'No .def found in Data folder.');
    [~, dataname] = fileparts(ddef(1).name);
end
