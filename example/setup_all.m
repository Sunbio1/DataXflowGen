clear;

% -----------------------------------------------------------
% 1) Define model name & initialize D2D
modelname = 'neu_model';
arInit;

% -----------------------------------------------------------
% 2) Load model & data and compile
arLoadModel(modelname);
arLoadData('neu_data', 1);
arCompileAll;

% -----------------------------------------------------------
% 3) Sanity bounds for critical parameters (before first fit)
% 3) Sanity bounds (before arInitValues)
for i = 1:numel(ar.pLabel)
    pnam = ar.pLabel{i};

    % Read current bounds safely
    lb = ar.lb(i); if ~isfinite(lb), lb = -Inf; end
    ub = ar.ub(i); if ~isfinite(ub), ub =  Inf; end

    % h_* in [1..10]
    if startsWith(pnam, 'h_')
        lb = max(lb, 1);
        ub = min(ub, 10);
    end

    % b_6_2 >= 0 (optional UB=5)
    if strcmp(pnam, 'b_6_2')
        lb = max(lb, 0);
        ub = min(ub, 5);   % comment out if too restrictive
    end

    % set upper bound for a_* (to avoid outliers like 1000)
    if startsWith(pnam, 'a_')
        ub = min(ub, 50);  % if scientifically reasonable, also 10/100
        % no LB -> negative a_* still allowed if the model supports it
    end

    % Safety: enforce lb < ub
    if ~(lb < ub)
        epsBW = 1e-9;
        if isfinite(ub)
            lb = ub - epsBW;
        else
            ub = lb + epsBW;
        end
    end

    ar.lb(i) = lb;
    ar.ub(i) = ub;
end

% Only then:
arInitValues;

arFit;

% -----------------------------------------------------------
% 4) Export before structure reduction
arExportPEtab('export_full');

% Chi² table before reduction
output_filename = 'chi2_full.csv';
writeChi2Table(output_filename);

% -----------------------------------------------------------
% 5) Enable L1 regularization
if ~isfield(ar, 'qL1reg') || isempty(ar.qL1reg)
    ar.qL1reg = zeros(size(ar.p));
end

ar.config.useL1 = 1;
ar.config.l1penalty = 50;

target_prefixes = {'a_', 'b_', 'h_', 'delta_'};
count = 0;
for i = 1:length(ar.pLabel)
    pName = ar.pLabel{i};
    if any(cellfun(@(prefix) startsWith(pName, prefix), target_prefixes))
        ar.qFit(i) = 1;
        ar.qL1reg(i) = 1;
        count = count + 1;
    end
end
save('qFit_reference.mat', 'ar');
fprintf('L1 regularization enabled for %d parameters (a_, b_, h_, delta_)\n', count);

% -----------------------------------------------------------
% 6) Configure the optimizer
if ~isfield(ar, 'config') || ~isstruct(ar.config)
    fprintf('ar.config was corrupted. Re-initialized.\n');
    ar.config = struct();
end

optimizer_name = 'lsqnonlin';
idx = find(strcmp(ar.config.optimizers, optimizer_name), 1);
assert(~isempty(idx), ['❌ Optimizer "' optimizer_name '" not found!']);
ar.config.optimizer = idx;
ar.config.maxsteps = 20000;
ar.config.useInitialGuess = 0; % set to 1 if the fit should reuse previous start values
ar.config.optimoptions = optimoptions(optimizer_name, ...
    'MaxIter', 20000, ...
    'Display', 'off');

% -----------------------------------------------------------
% 7) Fit with L1
disp('DEBUG: ar.config before arFit:');
disp(ar.config);

arFit;
arChi2Test;
arCalcMerit(true);
dof = ar.ndata - sum(ar.qFit);
fprintf('Chi² = %.2f | nData = %d | nParam = %d | Degrees of freedom = %d\n', ...
        ar.chi2fit, ar.ndata, sum(ar.qFit), dof);

% -----------------------------------------------------------
% 8) Structure reduction
threshold = 10;
count_removed = 0;
for i = 1:length(ar.p)
    if ar.qL1reg(i) && abs(ar.p(i)) < threshold
        ar.qFit(i) = 0;
        count_removed = count_removed + 1;
    end
end
fprintf('%d parameters below threshold %.1e were fixed (value kept, not set to 0).\n', count_removed, threshold);

% -----------------------------------------------------------
% 9) Refit without L1
ar.config.useL1 = 0;
ar.config.l1penalty = 0;
ar.config.useInitialGuess = 1;   % start from the L1 result

fprintf('Refit without L1 regularization using the reduced model structure...\n');
arFit;
arChi2Test;
arCalcMerit(true);

dof = ar.ndata - sum(ar.qFit);
fprintf('Refit result: Chi² = %.2f | nData = %d | nParam = %d | Degrees of freedom = %d\n', ...
        ar.chi2fit, ar.ndata, sum(ar.qFit), dof);

%-------------------------------------------------------------------
export_kept_interactions('kept_interactions.csv', 1e-6);

% -----------------------------------------------------------
% 10) AIC/BIC
chi2   = arGetMerit();
nParam = sum(ar.qFit);
nData  = ar.ndata;

AIC = 2 * nParam + chi2;
BIC = nParam * log(nData) + chi2;

fprintf('Chi² = %.2f | nParam = %d | nData = %d | AIC = %.2f | BIC = %.2f\n', ...
        chi2, nParam, nData, AIC, BIC);

% -----------------------------------------------------------
% 11) Export after reduction
arExportPEtab('export_reduced');

% Chi² table after reduction
output_filename = 'chi2_reduced.csv';
writeChi2Table(output_filename);

% -----------------------------------------------------------
% 12) Auto-lever diagnostics (robust, always write a summary)
diagOut = fullfile(pwd,'diag_out',datestr(now,'yyyymmdd_HHMMSS'));
if ~exist(diagOut,'dir'), mkdir(diagOut); end
summtxt = fullfile(diagOut,'summary.txt');
fid = fopen(summtxt,'w');
if fid<0, error('Cannot open summary.txt: %s', summtxt); end
fprintf(fid,'Diagnostics folder: %s\n', diagOut);

try
    addpath(genpath(fullfile(pwd,'Util')));  % add auto_* to path
    arSimu(false,true,true);

    has_from = exist('auto_from_csv','file')==2;
    has_lev  = exist('auto_find_levers','file')==2;
    fprintf(fid,'auto_from_csv=%d | auto_find_levers=%d\n', has_from, has_lev);

    % Prefer external CSV, otherwise create on the fly
    externalCsv = fullfile(pwd,'l_kras6','l6','gene_with_regulation.csv');
    if exist(externalCsv,'file')==2
        geneCsv = externalCsv; usedCsv = 'external';
    else
        geneCsv = fullfile(diagOut,'gene_with_regulation.csv');
        % IMPORTANT: onlyStates=false so measurement names != state names are included
        build_gene_with_regulation(ar, geneCsv, 'minChi2',0, 'onlyStates',false, 'decision','median');
        usedCsv = 'generated';
    end
    fprintf(fid,'gene_with_regulation.csv: %s (%s)\n', geneCsv, usedCsv);

    % Options
    opt_measure     = 'max';
    opt_graphScope  = 'upstream';
    opt_factorUp    = 2.0;
    opt_factorDown  = 0.2;
    opt_tlim        = [];
    opt_onlyFree    = true;
    opt_includeInit = false;
    opt_applyTop    = false;

    % Diagnostics (prefer batch, otherwise single EGFR, or Gene of interest!!!)
    out = []; usedMode = '(none)';
    if has_from && exist(geneCsv,'file')==2
        usedMode = 'batch';
        out = auto_from_csv(ar, geneCsv, ...
            'measure',opt_measure,'graphScope',opt_graphScope, ...
            'factorUp',opt_factorUp,'factorDown',opt_factorDown, ...
            'tlim',opt_tlim,'onlyFree',opt_onlyFree,'includeInit',opt_includeInit, ...
            'applyTop',opt_applyTop,'outdir',diagOut,'writeCSV',true);
    elseif has_lev
        usedMode = 'single';
        out = auto_find_levers(ar, ...
            'target','EGFR','desired','up', ...
            'measure',opt_measure,'graphScope',opt_graphScope, ...
            'factorUp',opt_factorUp,'factorDown',opt_factorDown, ...
            'tlim',opt_tlim,'onlyFree',opt_onlyFree,'includeInit',opt_includeInit, ...
            'outdir',diagOut,'writeCSV',true,'applyTop',opt_applyTop);
    else
        error('Diagnostics scripts are not on the path.');
    end
    fprintf(fid,'Mode: %s\n', usedMode);

    % If no CSVs were produced: force a fallback
    csvList = listAllCSVs(diagOut);
    if isempty(csvList)
        fprintf(fid,'WARN: No CSVs produced. Fallback (global scope, onlyFree=false)…\n');
        if has_from && exist(geneCsv,'file')==2
            out = auto_from_csv(ar, geneCsv, ...
                'measure',opt_measure,'graphScope','global', ...
                'factorUp',opt_factorUp,'factorDown',opt_factorDown, ...
                'tlim',opt_tlim,'onlyFree',false,'includeInit',opt_includeInit, ...
                'applyTop',false,'outdir',diagOut,'writeCSV',true);
            usedMode = 'batch-fallback';
        elseif has_lev
            out = auto_find_levers(ar, ...
                'target','EGFR','desired','up', ...
                'measure',opt_measure,'graphScope','global', ...
                'factorUp',opt_factorUp,'factorDown',opt_factorDown, ...
                'tlim',opt_tlim,'onlyFree',false,'includeInit',opt_includeInit, ...
                'outdir',diagOut,'writeCSV',true,'applyTop',false);
            usedMode = 'single-fallback';
        end
        fprintf(fid,'Mode (after fallback): %s\n', usedMode);
        csvList = listAllCSVs(diagOut);
    end

    % Short report
    if isstruct(out) && isfield(out,'rows')
        for i=1:numel(out.rows)
            r = out.rows(i);
            if ~isempty(r.rep) && isfield(r.rep,'top') && ~isempty(r.rep.top)
                fprintf(fid,'[%d] %s (%s, Chi2=%.3f): TOP %s %s x%.3g | Δmax=%.3g ΔAUC=%.3g score=%.3g\n', ...
                    i, r.gene, r.regulation, r.chi2, r.rep.top.param, r.rep.top.direction, ...
                    r.rep.top.factor, r.rep.top.delta_max, r.rep.top.delta_auc, r.rep.top.score);
            else
                fprintf(fid,'[%d] %s: no stable recommendation.\n', i, r.gene);
            end
        end
    elseif isstruct(out) && isfield(out,'top')
        fprintf(fid,'Single: TOP %s %s x%.3g | Δmax=%.3g ΔAUC=%.3g score=%.3g\n', ...
            out.top.param, out.top.direction, out.top.factor, ...
            out.top.delta_max, out.top.delta_auc, out.top.score);
    else
        fprintf(fid,'WARN: No output struct available.\n');
    end

    % Log CSV list
    csvList = listAllCSVs(diagOut);
    fprintf(fid,'\nCSV files (%d):\n', numel(csvList));
    for k=1:numel(csvList), fprintf(fid,'- %s\n', csvList{k}); end

catch ME
    fprintf(fid,'ERROR: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf(fid,'STACK:\n');
        for s=1:numel(ME.stack)
            fprintf(fid,'  %s:%d\n', ME.stack(s).file, ME.stack(s).line);
        end
    end
end
fclose(fid);
fprintf('[AUTO] Summary: %s\n', summtxt);

% ===========================================================
% Local helper functions (must be at the end)
function b = getBound(v, i, def)
    if i<=numel(v) && ~isempty(v(i)), b = v(i); else, b = def; end
end

% -----------------------------------------------------------
% Function to export Chi² (with protection against yStd==0)
function writeChi2Table(filename)
    global ar
    results = [];
    for m = 1:length(ar.model)
        model = ar.model(m);
        for d = 1:length(model.data)
            data = model.data(d);
            if ~isfield(data,'tExp') || isempty(data.tExp), continue; end
            tExp = data.tExp;
            yExp = data.yExp;
            ySim = data.yExpSimu;
            yStd = data.yExpStd;
            yNames = data.yNames;

            if isempty(yStd), yStd = ones(size(yExp)); end
            yStd = max(yStd, eps);

            chi2vec = ((yExp - ySim) ./ yStd).^2;
            for iy = 1:length(yNames)
                T = table();
                T.Time           = tExp;
                T.Variable       = repmat(yNames(iy), length(tExp), 1);
                T.ExpData        = yExp(:, iy);
                T.SimulatedData  = ySim(:, iy);
                T.Std            = yStd(:, iy);
                T.Chi2           = chi2vec(:, iy);
                T.ModelIndex     = repmat(m, length(tExp), 1);
                T.DataIndex      = repmat(d, length(tExp), 1);
                results = [results; T];
            end
        end
    end
    writetable(results, filename);
    fprintf('Chi² table saved as %s\n', filename);
end

% -----------------------------------------------------------
% On-the-fly: build gene_with_regulation.csv (Variable,Chi2,Regulation)
function build_gene_with_regulation(ar, outCsv, varargin)
    ip = inputParser;
    ip.addParameter('minChi2',0,@(x)isnumeric(x)&&isscalar(x));
    ip.addParameter('onlyStates',true,@islogical);
    ip.addParameter('decision','median',@(x)ischar(x)||isstring(x));
    ip.parse(varargin{:});
    opt = ip.Results;

    % Collect states
    stateSet = string([]);
    for m=1:numel(ar.model)
        stateSet = unique([stateSet; string(ar.model(m).x(:))]);
    end

    S = struct(); % aggregate per variable
    for m = 1:numel(ar.model)
        mdl = ar.model(m);
        for d = 1:numel(mdl.data)
            dat = mdl.data(d);
            if ~isfield(dat,'yExp') || isempty(dat.yExp) || ~isfield(dat,'yExpSimu') || isempty(dat.yExpSimu)
                continue;
            end
            yNames = string(dat.yNames);
            Yexp = dat.yExp; Ysim = dat.yExpSimu;
            Ystd = dat.yExpStd; if isempty(Ystd), Ystd = ones(size(Yexp)); end
            Ystd = max(Ystd, eps);

            for iy = 1:numel(yNames)
                vname = char(yNames(iy));
                res   = Yexp(:,iy) - Ysim(:,iy);

                wres  = res ./ Ystd(:,iy);
                chi2v = nansum(wres.^2);

                if ~isfield(S, vname), S.(vname) = struct('chi2',0,'res',[]); end
                S.(vname).chi2 = S.(vname).chi2 + chi2v;
                S.(vname).res  = [S.(vname).res; res(~isnan(res))];
            end
        end
    end

    vars = fieldnames(S);
    outRows = {};
    for i=1:numel(vars)
        v = vars{i};
        if opt.onlyStates && ~any(stateSet == string(v)), continue; end
        chi2v = S.(v).chi2; if chi2v < opt.minChi2, continue; end
        r = S.(v).res; if isempty(r), continue; end

        if strcmpi(opt.decision,'median'), score = median(r,'omitnan'); else, score = mean(r,'omitnan'); end
        regulation = "Activation"; if score < 0, regulation = "Inhibition"; end

        outRows(end+1,:) = {v, chi2v, char(regulation)}; %#ok<AGROW>
    end

    T = cell2table(outRows, 'VariableNames', {'Variable','Chi2','Regulation'});
    writetable(T, outCsv);
    fprintf('gene_with_regulation.csv written: %s  (n=%d)\n', outCsv, height(T));
end

% -----------------------------------------------------------
% List CSVs in a folder (recursive)
function L = listAllCSVs(root)
    pp = regexp(genpath(root), pathsep, 'split'); L = {};
    for ii=1:numel(pp)
        if isempty(pp{ii}), continue; end
        d = dir(fullfile(pp{ii},'*.csv'));
        for jj=1:numel(d), L{end+1} = fullfile(d(jj).folder,d(jj).name); end %#ok<AGROW>
    end
end
