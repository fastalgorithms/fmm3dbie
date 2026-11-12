%% Sweep all supported two-sphere AMR options
% Run this script directly. Edit the combinations below.
% Calls amr_lapl_density_new.m for each case; saves completed rows to CSV.
% results and reference remain in the workspace when the sweep finishes.

smoke = true; % true: one strong/tails_vals case; false: full sweep

% Supported error/estimator pairs (weak and functional select their own estimator).
pairs = [ ...
    "strong",     "tails_vals";
    "strong",     "tails_coefs";
    "strong",     "upsample";
    "strong",     "nystrom";
    "weak",       "weak";
    "functional", "functional"];

lptypes = [0 1 2];             % Linf, L1, L2
etas = [0 1 2];                % patch-length scaling exponents
tols = [1e-5 1e-6 1e-7 1e-8];      % full experiment

% Smoke mode runs one smaller two-sphere problem with the full-experiment gap.
% The example still builds/checks a fresh reference and runs AMR.
if smoke
    pairs = ["strong", "tails_vals"];
    lptypes = 2;
    etas = 1;
    tols = 1e-3;
    fprintf('Smoke mode: one smaller case, same sphere gap as full mode.\n');
end

script_dir = fileparts(mfilename('fullpath'));
library_dir = fileparts(script_dir);
addpath(library_dir, fullfile(library_dir,'src'), script_dir);
assert(isfile(fullfile(script_dir,'amr_lapl_density_new.m')), ...
    'Place this script beside amr_lapl_density_new.m.');

valid_pairs = ["strong","tails_vals"; "strong","tails_coefs"; ...
    "strong","upsample"; "strong","nystrom"; "weak","weak"; ...
    "functional","functional"];
pairs = string(pairs);
assert(~isempty(pairs) && ismatrix(pairs) && size(pairs,2)==2, ...
    'Supply one [errtype,esttype] pair per row.');
for k = 1:size(pairs,1)
    assert(any(all(valid_pairs==pairs(k,:),2)), ...
        'Unsupported pair: %s / %s.',pairs(k,1),pairs(k,2));
end
validateattributes(lptypes,{'numeric'},{'vector','nonempty','real','finite'});
assert(all(ismember(lptypes,[0 1 2])),'Use 0 (Linf), 1, or 2 for lptype.');
validateattributes(etas,{'numeric'},{'vector','nonempty','real','finite'});
validateattributes(tols,{'numeric'},{'vector','nonempty','real','finite','positive'});

% Keep completed rows if a later solve fails or the sweep is interrupted.
csv_file = [tempname(pwd) '_amr_sweep.csv'];
ncase = size(pairs,1)*numel(lptypes)*numel(etas)*numel(tols);
fprintf('Running %d cases. Results file: %s\n',ncase,csv_file);
results = table();
reference = []; % verified once, then shared by matching cases
for k = 1:size(pairs,1)
    for lptype = lptypes(:).'
        for eta_dens = etas(:).'
            for tol_amr = tols(:).'
                settings = struct('errtype',char(pairs(k,1)), ...
                    'esttype',char(pairs(k,2)),'lptype',lptype, ...
                    'eta_dens',eta_dens,'tol_amr',tol_amr, ...
                    'reference_tol',min(tols)/10,'reference',reference, ...
                    'smoke',smoke);
                fprintf('\nCase %d/%d: %s, %s, p=%g, eta=%g, tol=%g\n', ...
                    height(results)+1,ncase,settings.errtype,settings.esttype, ...
                    lptype,eta_dens,tol_amr);
                row = table(pairs(k,1),pairs(k,2),lptype,eta_dens,tol_amr, ...
                    'VariableNames',{'errtype','esttype','lptype','eta_dens','tol_amr'});
                metrics = array2table(nan(1,7),'VariableNames', ...
                    {'irefine','npatches','npts','actual_abs','actual_rel', ...
                     'estimated_abs','estimated_rel'});
                error_message = "";
                try
                    [values,reference] = run_case(settings);
                    metrics{1,:} = values;
                catch ME
                    error_message = string(getReport(ME,'extended','hyperlinks','off'));
                    fprintf(2,'%s\n',error_message);
                    if startsWith(ME.identifier,'two_sphere_reference:') || ...
                            strcmp(ME.identifier,'refine_close_patches:PointLimit')
                        rethrow(ME)
                    end
                end
                row = [row,metrics,table(error_message)];
                if isempty(results)
                    results = row;
                    writetable(row,csv_file);
                else
                    results = [results;row]; %#ok<AGROW>
                    writetable(row,csv_file,'WriteMode','append','WriteVariableNames',false);
                end
            end
        end
    end
end

disp(results);

function [values,reference] = run_case(amr_sweep_case) %#ok<INUSD>
% Isolate each example run so its variables cannot overwrite the sweep loop.
reference = amr_sweep_case.reference; %#ok<NASGU>
run(fullfile(fileparts(mfilename('fullpath')),'amr_lapl_density_new.m'));
if irefine == 0
    values = [irefine,S.npatches,S.npts,fun_monitor_lp_abs, ...
        fun_monitor_lp_abs/sig_p_rhs_norm_lp,NaN,NaN];
else
    values = [irefine,S.npatches,S.npts,err_monitor_lp_abs_comp, ...
        err_monitor_lp_abs_comp/fun_monitor_norm_S,err_monitor_est_abs,err_monitor_est_rel];
end
end
