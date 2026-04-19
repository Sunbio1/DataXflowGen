function kept = export_kept_interactions(csvfile, threshold_lin)
    % Lists all remaining interactions (a_*) with linear gain > threshold_lin
    % and writes them to a CSV (Parameter, GainLinear).
    if nargin < 1 || isempty(csvfile),        csvfile = 'kept_interactions.csv'; end
    if nargin < 2 || isempty(threshold_lin),  threshold_lin = 1e-6;              end

    global ar
    % Candidates: a_* and currently free (qFit==1)
    idx = find(ar.qFit==1 & startsWith(ar.pLabel,'a_'));

    labels  = ar.pLabel(idx);
    gains   = ar.p(idx);               % may be log10-scaled
    islog10 = ar.qLog10(idx)==1;

    % Convert to linear scale
    gains_lin = gains;
    gains_lin(islog10) = 10.^gains_lin(islog10);

    % Filter: truly "active" (linear > threshold_lin)
    keepmask = gains_lin > threshold_lin;

    kept = labels(keepmask);           % return for console use

    % Write CSV (keep only active ones)
    T = table(labels(keepmask)', gains_lin(keepmask)', ...
              'VariableNames', {'Parameter','GainLinear'});
    writetable(T, csvfile);

    fprintf('✅ %d of %d interactions were kept. File: %s\n', ...
            sum(keepmask), numel(idx), csvfile);

    % Optional: print to console for verification
    disp(kept);
end
