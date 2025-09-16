%% Clear Workspace

clear; clc; close all;

%% P2T1

% Store case files based on file structure (change this if you change
% location or even type)
case_paths = {'cases/Case_1.csv','cases/Case_2.csv','cases/Case_3.csv'};

% Iterate through case files, read tables with *preserved names*
for case_num = 1:numel(case_paths)
    cases(case_num).tbl = readtable(case_paths{case_num});
end

% Names of dependent variables and the loads
deps = {'F0_lbf_','F1_lbf_','F2_lbf_','F3D_lbf_', 'LVDT_in_'};
loads = 'LoadingCase_lbs_';

% Loop for evaluating polyfit/polyval and plotting them against raw data
for dep_num = 1:numel(deps)

    % Set current dependent variable
    curr_dep = [deps{dep_num}];

    % Enumerate through all cases, fit/val them
    for case_num = 1:numel(cases)
        x = cases(case_num).tbl.(loads);
        y = cases(case_num).tbl.(curr_dep);
        cases(case_num).fit.(deps{dep_num}) = polyfit(x, y, 1);
        cases(case_num).fit_coef.(deps{dep_num}) = polyval(cases(case_num).fit.(deps{dep_num}), x);
    end

    % Set up individual plots for each dependent variable, not case
    figure();
    hold on;
    grid on;

    % Loop through cases, plot relevant dependent variable altogether!
    for case_num = 1:numel(cases)
        x = cases(case_num).tbl.(loads);
        y = cases(case_num).tbl.([deps{dep_num}]);
        plot(x, y, 'o', 'DisplayName', sprintf('Case %d, raw', case_num));
        plot(x, cases(case_num).fit_coef.(deps{dep_num}), 'DisplayName', sprintf('Case %d, fit', case_num));
    end

    % Make it pretty
    title(sprintf('All cases of %s w/ linear fits', deps{dep_num}), 'Interpreter', 'none');
    xlabel('Loading case (lb)');
    ylabel(sprintf('%s (lbf)', deps{dep_num}), 'Interpreter', 'none');
    legend('Location','best');

    % Save figure into figures folder (hidden on git)
    savefig(sprintf('figs/Dep%s.fig', deps{dep_num}));
    saveas(gcf, sprintf('figs/Dep%s.png', deps{dep_num}));

end