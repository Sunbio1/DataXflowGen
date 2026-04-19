function R = fit_one_runfolder_finalfit_struct(runDir, varargin)
% Fit one run folder and report the final fit plus basic structure metrics.
% - Start values: all fitted (qFit) parameters start at startValue (default 0.2).
% - mode="ai": run a single fit (arFit), no multi-start (no LHS).
% - mode="random": fit once from the start value, then run LHS multi-start (arFitLHS(nStarts));
%   keep the best result.
% - useL1: use L1 regularization during the structure fit; optionally refit without L1 for scoring.
%
% Returns (main outputs):
% chi2 at start, after first fit, and final; number of fitted params (k); ndata; degrees of freedom (dof);
% p-value (theory); AIC/AICc/BIC; model hash, etc.


p = inputParser;
p.addParameter('nStarts', 20, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('L1penalty', 50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('seed', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('useL1', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('mode', 'random', @(s)ischar(s)||isstring(s)); % "random" or "ai"
p.addParameter('refitWithoutL1', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('freezeZeroAfterL1', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('zeroTol', 1e-6, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Enforce start values
p.addParameter('forceStartValue', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('startValue', 0.2, @(x)isnumeric(x)&&isscalar(x));

p.parse(varargin{:});
opt = p.Results;

t0 = tic;

R = struct();
R.runDir = string(runDir);
R.modelname = "";
R.dataname  = "";
R.modelDef  = "";
R.dataDef   = "";
R.modelHash = "";
R.modelHashShort = "";

R.ok = 0;
R.errmsg = "";
R.seconds = NaN;

R.mode = string(opt.mode);

R.chi2_start      = NaN;
R.chi2_afterFirst = NaN;   % after first fit
R.chi2_final      = NaN;   % after optional refit without L1

R.k_start = NaN;
R.k       = NaN;
R.ndata   = NaN;
R.dof     = NaN;

R.p_theory    = NaN;
R.p_empirisch = NaN;

R.AIC  = NaN;
R.AICc = NaN;
R.BIC  = NaN;

R.nRegParams_total = NaN;
R.nRegParams_fit   = NaN;
R.nSdParams_total  = NaN;
R.nSdParams_fit    = NaN;
R.nScaleParams_total = NaN;
R.nScaleParams_fit   = NaN;

R.nFrozenZero = 0;

old = pwd;
cleanup = onCleanup(@() cd(old));
cd(runDir);

try
    ensure_d2d_layout(pwd);
    [modelname, dataname] = autodetect_model_data(pwd);
    R.modelname = string(modelname);
    R.dataname  = string(dataname);

    modelDef = fullfile(pwd, 'Models', [modelname '.def']);
    dataDef  = fullfile(pwd, 'Data',   [dataname  '.def']);
    R.modelDef = string(modelDef);
    R.dataDef  = string(dataDef);

    R.modelHash = string(md5_of_files({modelDef, dataDef}));
    if strlength(R.modelHash) >= 8
        R.modelHashShort = extractBetween(R.modelHash, 1, 8);
    else
        R.modelHashShort = R.modelHash;
    end

    clear global ar
    global ar
    ar = [];

    rng(opt.seed);

    arInit;
    arLoadModel(modelname);
    arLoadData(dataname, 1);
    arCompileAll;

    % ---- CVODES robustness ----
    if isfield(ar,'config')
        if isfield(ar.config,'cvodes_maxsteps'), ar.config.cvodes_maxsteps = 200000; end
        if isfield(ar.config,'cvodes_mxstep'),   ar.config.cvodes_mxstep   = 200000; end
        if isfield(ar.config,'mxstep'),          ar.config.mxstep          = 200000; end
        if isfield(ar.config,'mxsteps'),         ar.config.mxsteps         = 200000; end
        if isfield(ar.config,'rtol'),            ar.config.rtol            = 1e-6; end
        if isfield(ar.config,'atol'),            ar.config.atol            = 1e-8; end
    end

    % ---- Bounds sanity ----
    for i = 1:numel(ar.pLabel)
        pnam = ar.pLabel{i};
        lb = ar.lb(i); if ~isfinite(lb), lb = -Inf; end
        ub = ar.ub(i); if ~isfinite(ub), ub =  Inf; end

        if startsWith(pnam, 'h_')
            lb = max(lb, 1); ub = min(ub, 10);
        end
        if strcmp(pnam, 'b_6_2')
            lb = max(lb, 0); ub = min(ub, 5);
        end
        if startsWith(pnam, 'a_')
            ub = min(ub, 50);
        end
        if ~(lb < ub)
            epsBW = 1e-9;
            if isfinite(ub), lb = ub - epsBW; else, ub = lb + epsBW; end
        end

        ar.lb(i) = lb;
        ar.ub(i) = ub;
    end

    arInitValues;

    % ---- L1 setup: do NOT touch qFit, only set qL1reg ----
    if ~isfield(ar,'qL1reg') || isempty(ar.qL1reg)
        ar.qL1reg = zeros(size(ar.p));
    end
    ar.qL1reg(:) = 0;

    if opt.useL1
        ar.config.useL1 = 1;
        ar.config.l1penalty = opt.L1penalty;

        target_prefixes = {'a_','b_','h_','delta_'};
        for ii = 1:length(ar.pLabel)
            pName = ar.pLabel{ii};
            if any(cellfun(@(prefix) startsWith(pName, prefix), target_prefixes))
                ar.qL1reg(ii) = 1;
            end
        end
    else
        ar.config.useL1 = 0;
        ar.config.l1penalty = 0;
    end

    % ---- Structure counting (based on D2D qFit) ----
    [R.nRegParams_total, R.nRegParams_fit] = count_by_prefix(ar.pLabel, ar.qFit, {'a_','b_','h_','delta_'});
    [R.nSdParams_total,  R.nSdParams_fit]  = count_by_prefix(ar.pLabel, ar.qFit, {'sd_','sigma_'});
    [R.nScaleParams_total, R.nScaleParams_fit] = count_by_prefix(ar.pLabel, ar.qFit, {'scale_','sc_','offset_'});

    % ---- Enforce start values (really all the same) ----
    if opt.forceStartValue
        idxFit = find(ar.qFit(:)==1);
        ar.p(idxFit) = opt.startValue;
        ar.p = max(min(ar.p, ar.ub), ar.lb); % clamp into bounds
    end

    % ---- START chi2 ----
    arCalcMerit(true);
    R.chi2_start = get_data_chi2(ar);
    R.k_start    = sum(ar.qFit);

    % ---- Optimizer ----
    optimizer_name = 'lsqnonlin';
    idx = find(strcmp(ar.config.optimizers, optimizer_name), 1);
    assert(~isempty(idx), ['Optimizer "' optimizer_name '" not found!']);
    ar.config.optimizer = idx;

    ar.config.maxsteps = 20000;
    ar.config.optimoptions = optimoptions(optimizer_name, ...
        'MaxIter', 20000, ...
        'MaxFunctionEvaluations', 500000, ...
        'Display', 'off');

    % ---- FIT: AI vs RANDOM ----
    mode = string(opt.mode);

    if mode == "ai"
        % exactly one start, no LHS
        ar.config.useInitialGuess = 1;
        arFit;

    else
        % Random: first fit from the same start (more stable baseline)
        ar.config.useInitialGuess = 1;
        arFit;
        arCalcMerit(true);
        chi2_best = get_data_chi2(ar);
        p_best = ar.p;

        % then add multi-start
        if opt.nStarts > 1
            ar.config.useInitialGuess = 0;
            arFitLHS(opt.nStarts);

            arCalcMerit(true);
            chi2_lhs = get_data_chi2(ar);

            if isfinite(chi2_lhs) && chi2_lhs < chi2_best
                chi2_best = chi2_lhs;
                p_best = ar.p;
            end

            % restore best result
            ar.p = p_best;
            arCalcMerit(true);
        end
    end

    % after first fit (possibly with L1)
    arCalcMerit(true);
    R.chi2_afterFirst = get_data_chi2(ar);

    % ---- optional: refit WITHOUT L1 (pure fit for scoring) ----
    if opt.useL1 && opt.refitWithoutL1
        ar.config.useL1 = 0;
        ar.config.l1penalty = 0;

        if opt.freezeZeroAfterL1 && isfield(ar,'qL1reg') && ~isempty(ar.qL1reg)
            tol = opt.zeroTol;
            idxReg  = (ar.qL1reg(:)==1);
            idxZero = idxReg & (abs(ar.p(:)) < tol);

            R.nFrozenZero = sum(idxZero);

            ar.p(idxZero)    = 0;
            ar.qFit(idxZero) = 0;
            ar.lb(idxZero)   = 0;
            ar.ub(idxZero)   = 0;
        end

        ar.config.useInitialGuess = 1;
        arFit;
    end

    % ---- FINAL ----
    arCalcMerit(true);
    R.chi2_final = get_data_chi2(ar);

    R.k     = sum(ar.qFit);
    R.ndata = ar.ndata;
    R.dof   = R.ndata - R.k;

    if isfinite(R.chi2_final) && isfinite(R.dof) && (R.dof > 0)
        R.p_theory = 1 - chi2cdf(R.chi2_final, R.dof);
    end

    chi2 = R.chi2_final;
    k = R.k;
    n = R.ndata;

    R.AIC = chi2 + 2*k;
    R.BIC = chi2 + k*log(n);

    den = (n - k - 1);
    if den > 0
        R.AICc = R.AIC + (2*k*(k+1))/den;
    else
        R.AICc = NaN;
    end

    R.ok = 1;
    R.seconds = toc(t0);

catch ME
    R.ok = 0;
    R.errmsg = string(ME.message);
    R.seconds = toc(t0);
end

end

% -------------------------------------------------------------------------
function chi2 = get_data_chi2(ar)
if isfield(ar,'chi2fit') && ~isempty(ar.chi2fit) && isfinite(ar.chi2fit)
    chi2 = ar.chi2fit;
else
    chi2 = ar.chi2;
end
end

% -------------------------------------------------------------------------
function [nTotal, nFit] = count_by_prefix(pLabel, qFit, prefixes)
mask = false(size(pLabel));
for j = 1:numel(prefixes)
    mask = mask | startsWith(string(pLabel), string(prefixes{j}));
end
nTotal = sum(mask);
nFit   = sum(mask(:) & (qFit(:)==1));
end

% -------------------------------------------------------------------------
function h = md5_of_files(files)
import java.security.*;
md = MessageDigest.getInstance('MD5');

for i=1:numel(files)
    f = files{i};
    if exist(f,'file')~=2, continue; end
    txt = fileread(f);
    md.update(uint8(txt));
end

d = typecast(md.digest(),'uint8');
h = lower(reshape(dec2hex(d)',1,[]));
end

% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
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
