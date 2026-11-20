function [ResultsTable, OuterLog, IterationHistory] = shaft_design_from_csv(csvFile)
% ========================================================================
% SHAFT DIAMETER DESIGN (FATIGUE) - AUTOMATIC FILLET UPDATE
% ------------------------------------------------------------------------
% Input CSV columns:
%   Section, Shaft, r_initial_mm, d_guess_mm, Ma_Nm, Mm_Nm, Ta_Nm, Tm_Nm
%
% Rows must be ordered from LEFT to RIGHT on each shaft.
% For every non-terminal section:
%   r_geom_max = abs(D_left - D_right)/2
% The rightmost section of each shaft has r = 0.
% Outer loop uses continuous radius relaxation. At the end, fillet radii
% are selected from a preferred standard-radius series, always rounding
% DOWN so the selected fillet does not exceed its geometric maximum.
% OUTPUT:
%   Clean, fixed-width report in the MATLAB Command Window.
%   Tables are also returned to the MATLAB workspace if requested.
% ========================================================================

    if nargin < 1 || strlength(string(csvFile)) == 0
        csvFile = input('Enter input CSV filename: ', 's');
    end

    if exist(csvFile, 'file') ~= 2
        error('Input CSV file was not found: %s', csvFile);
    end

    %% =========================
    % READ AND VALIDATE INPUT CSV
    %% =========================
    InputTable = readtable(csvFile, 'TextType', 'string');

    requiredColumns = {'Section','Shaft','r_initial_mm','d_guess_mm', ...
                       'Ma_Nm','Mm_Nm','Ta_Nm','Tm_Nm'};
    for k = 1:numel(requiredColumns)
        if ~ismember(requiredColumns{k}, InputTable.Properties.VariableNames)
            error('Missing required CSV column: %s', requiredColumns{k});
        end
    end

    secName = string(InputTable.Section);
    shaftID = string(InputTable.Shaft);
    rInitial = InputTable.r_initial_mm;
    r = rInitial;
    dGuess = InputTable.d_guess_mm;
    Ma_Nm = InputTable.Ma_Nm;
    Mm_Nm = InputTable.Mm_Nm;
    Ta_Nm = InputTable.Ta_Nm;
    Tm_Nm = InputTable.Tm_Nm;

    nSec = height(InputTable);
    if nSec < 4
        error('The input must contain the sections of both shafts.');
    end
    if any(ismissing(secName)) || any(strlength(secName) == 0)
        error('Every row must have a Section label.');
    end
    if any(ismissing(shaftID)) || any(strlength(shaftID) == 0)
        error('Every row must have a Shaft label.');
    end

    numericInputs = [rInitial; dGuess; Ma_Nm; Mm_Nm; Ta_Nm; Tm_Nm];
    if any(~isfinite(numericInputs))
        error('Numeric CSV inputs must all be finite numbers.');
    end
    if any(rInitial < 0) || any(dGuess <= 0)
        error('Fillet radii must be >= 0 and diameter guesses must be > 0.');
    end

    shaftNames = unique(shaftID, 'stable');
    if numel(shaftNames) ~= 2
        error('The input CSV must contain exactly two shaft groups.');
    end

    shaftRows = cell(numel(shaftNames),1);
    for s = 1:numel(shaftNames)
        shaftRows{s} = find(shaftID == shaftNames(s));
        if numel(shaftRows{s}) < 2
            error('Each shaft must contain at least two ordered sections.');
        end
    end

    % Convert N.m to N.mm, as in the original code.
    Ma = Ma_Nm * 1000;
    Mm = Mm_Nm * 1000;
    Ta = Ta_Nm * 1000;
    Tm = Tm_Nm * 1000;

    %% =========================
    % MATERIAL PROPERTIES (UNCHANGED)
    %% =========================
    Sut = 630;     % MPa
    Sy  = 370;     % MPa (kept from the original code)
    n   = 2;       % factor of safety

    ka = 4.51 * Sut^(-0.265);   % machined surface factor
    Se0 = 0.5 * Sut;            % unmodified endurance limit

    %% =========================
    % ITERATION SETTINGS
    %% =========================
    innerTol = 0.1;             % mm; original Bisection tolerance
    outerTol = 0.1;             % mm; outer diameter convergence tolerance
    maxOuterIter = 200;
    radiusRelaxation = 0.30;    % continuous update factor for outer loop

    % Preferred fillet-radius values [mm].
    % A radius is selected DOWN to the greatest listed value <= geometric max.
    standardRadii = [0 0.5 0.8 1.0 1.2 1.5 2.0 2.5 3.0 4.0 5.0 ...
                     6.0 8.0 10.0 12.0 16.0 20.0 25.0 32.0];

    %% =========================
    % PREALLOCATE RESULTS
    %% =========================
    dRaw = zeros(nSec,1);
    dReport = zeros(nSec,1);
    dStd = zeros(nSec,1);
    iterNo = zeros(nSec,1);

    KtOut = zeros(nSec,1);
    KtsOut = zeros(nSec,1);
    KfOut = zeros(nSec,1);
    KfsOut = zeros(nSec,1);
    kbOut = zeros(nSec,1);
    SeOut = zeros(nSec,1);
    qOut = zeros(nSec,1);
    qsOut = zeros(nSec,1);

    OuterLog = table('Size',[0 4], ...
        'VariableTypes', {'double','double','double','logical'}, ...
        'VariableNames', {'OuterIteration','MaxDiameterChange_mm', ...
                          'MaxRadiusResidual_mm','Converged'});

    IterationHistory = table('Size',[0 10], ...
        'VariableTypes', {'string','double','string','string','double', ...
                          'double','double','double','double','double'}, ...
        'VariableNames', {'Stage','OuterIteration','Shaft','Section', ...
                          'FilletRadiusUsed_mm','FilletRadiusTarget_mm', ...
                          'D_raw_mm','D_reported_mm','D_standard_mm', ...
                          'InnerBisectionIterations'});

    prevDRaw = [];
    converged = false;
    maxDiameterChange = NaN;
    targetR = zeros(nSec,1);

    %% ====================================================================
    % OUTER LOOP: CONTINUOUS DIAMETER <-> FILLET-RADIUS UPDATE
    %% ====================================================================
    for outer = 1:maxOuterIter

        % Inner Bisection calculations using the current radius values.
        for i = 1:nSec
            [dRaw(i), KtOut(i), KtsOut(i), KfOut(i), KfsOut(i), ...
                kbOut(i), SeOut(i), iterNo(i), qOut(i), qsOut(i)] = ...
                solveSection(r(i), dGuess(i), Ma(i), Mm(i), Ta(i), Tm(i), ...
                             Sut, n, ka, Se0, innerTol);

            dReport(i) = ceil(dRaw(i) * 10) / 10;  % original rule
            dStd(i) = nextStandardDiameter(dReport(i));
        end

        % Compute geometric maximum radii from left to right on each shaft.
        targetR = calculateGeometricRadii(dRaw, shaftRows);

        if isempty(prevDRaw)
            maxDiameterChange = NaN;
        else
            maxDiameterChange = max(abs(dRaw - prevDRaw));
        end
        maxRadiusResidual = max(abs(targetR - r));
        convergedThisPass = ~isempty(prevDRaw) && ...
                            (maxDiameterChange <= outerTol + 1e-9);

        OuterLog = [OuterLog; table(outer, maxDiameterChange, ...
            maxRadiusResidual, convergedThisPass, ...
            'VariableNames', {'OuterIteration','MaxDiameterChange_mm', ...
                              'MaxRadiusResidual_mm','Converged'})]; %#ok<AGROW>

        historyRows = table(repmat("Continuous",nSec,1), repmat(outer,nSec,1), ...
            shaftID, secName, r, targetR, dRaw, dReport, dStd, iterNo, ...
            'VariableNames', {'Stage','OuterIteration','Shaft','Section', ...
            'FilletRadiusUsed_mm','FilletRadiusTarget_mm','D_raw_mm', ...
            'D_reported_mm','D_standard_mm','InnerBisectionIterations'});
        IterationHistory = [IterationHistory; historyRows]; %#ok<AGROW>

        if convergedThisPass
            converged = true;
            break;
        end

        % Under-relaxed continuous update; end section of each shaft = 0.
        r = r + radiusRelaxation * (targetR - r);
        for s = 1:numel(shaftRows)
            idx = shaftRows{s};
            r(idx(end)) = 0;
        end
        prevDRaw = dRaw;
    end

    if ~converged
        warning('Maximum outer iterations reached before diameter convergence.');
    end

    %% ====================================================================
    % STANDARDIZE THE FILLET RADII AND RE-RUN THE FINAL DESIGN
    %% ====================================================================
    % First choose the greatest standard radius not exceeding the geometric
    % target from the last continuous outer-loop pass.
    rStandard = zeros(nSec,1);
    for i = 1:nSec
        rStandard(i) = nextStandardRadius(targetR(i), standardRadii);
    end
    for s = 1:numel(shaftRows)
        idx = shaftRows{s};
        rStandard(idx(end)) = 0;
    end

    % Recalculate with standard radii. If any selected radius is larger than
    % the geometric maximum generated by the recalculated diameters, step
    % that radius DOWN to the next permissible standard value. Radii are
    % never increased in this correction loop, so it cannot oscillate.
    maxRadiusCorrectionPasses = numel(standardRadii) + 2;
    for correctionPass = 1:maxRadiusCorrectionPasses

        for i = 1:nSec
            [dRaw(i), KtOut(i), KtsOut(i), KfOut(i), KfsOut(i), ...
                kbOut(i), SeOut(i), iterNo(i), qOut(i), qsOut(i)] = ...
                solveSection(rStandard(i), dGuess(i), Ma(i), Mm(i), Ta(i), Tm(i), ...
                             Sut, n, ka, Se0, innerTol);

            dReport(i) = ceil(dRaw(i) * 10) / 10;
            dStd(i) = nextStandardDiameter(dReport(i));
        end

        finalGeometricMaxR = calculateGeometricRadii(dRaw, shaftRows);
        rCorrected = rStandard;

        for i = 1:nSec
            if rStandard(i) > finalGeometricMaxR(i) + 1e-9
                rCorrected(i) = nextStandardRadius(finalGeometricMaxR(i), standardRadii);
            end
        end
        for s = 1:numel(shaftRows)
            idx = shaftRows{s};
            rCorrected(idx(end)) = 0;
        end

        if all(abs(rCorrected - rStandard) < 1e-9)
            break;
        end
        rStandard = rCorrected;
    end

    % One final table/history row set for the actual standardized-radius run.
    finalGeometricMaxR = calculateGeometricRadii(dRaw, shaftRows);
    finalHistoryRows = table(repmat("Final standardized",nSec,1), ...
        repmat(outer + correctionPass,nSec,1), shaftID, secName, ...
        rStandard, finalGeometricMaxR, dRaw, dReport, dStd, iterNo, ...
        'VariableNames', {'Stage','OuterIteration','Shaft','Section', ...
        'FilletRadiusUsed_mm','FilletRadiusTarget_mm','D_raw_mm', ...
        'D_reported_mm','D_standard_mm','InnerBisectionIterations'});
    IterationHistory = [IterationHistory; finalHistoryRows]; %#ok<AGROW>

    %% =========================
    % FINAL RESULTS TABLE (RETURNED TO WORKSPACE, NOT WRITTEN TO CSV)
    %% =========================
    ResultsTable = table(shaftID, secName, dRaw, dReport, dStd, ...
        finalGeometricMaxR, rStandard, Ma_Nm, Mm_Nm, Ta_Nm, Tm_Nm, ...
        KtOut, KtsOut, KfOut, KfsOut, qOut, qsOut, kbOut, SeOut, ...
        iterNo, repmat(outer,nSec,1), repmat(converged,nSec,1), ...
        'VariableNames', {'Shaft','Section','D_raw_mm','D_reported_mm', ...
        'D_standard_mm','FilletGeometricMax_mm','FilletStandardUsed_mm', ...
        'Ma_Nm','Mm_Nm','Ta_Nm','Tm_Nm','Kt','Kts','Kf','Kfs','q','qs', ...
        'kb','Se_MPa','InnerBisectionIterations','OuterIterations','Converged'});

    %% =========================
    % CLEAN COMMAND WINDOW REPORT
    %% =========================
    printFinalReport(shaftID, secName, dRaw, dStd, ...
        finalGeometricMaxR, rStandard, Ma_Nm, Mm_Nm, Ta_Nm, Tm_Nm, ...
        KtOut, KtsOut, KfOut, KfsOut, qOut, qsOut, kbOut, SeOut, ...
        iterNo, outer, converged, maxDiameterChange, outerTol, ...
        innerTol, radiusRelaxation, OuterLog);

    % Keep Sy declared as in the supplied original source.
    %#ok<NASGU>
    Sy = Sy;
end

%% ========================================================================
% GEOMETRIC FILLET-RADIUS TARGETS
%% ========================================================================
function targetR = calculateGeometricRadii(d, shaftRows)
    targetR = zeros(size(d));
    for s = 1:numel(shaftRows)
        idx = shaftRows{s};
        for j = 1:(numel(idx)-1)
            iLeft = idx(j);
            iRight = idx(j+1);
            targetR(iLeft) = abs(d(iLeft) - d(iRight)) / 2;
        end
        targetR(idx(end)) = 0;
    end
end

function rStd = nextStandardRadius(rGeom, standardRadii)
    % Select the greatest standard radius that does not exceed rGeom.
    idx = find(standardRadii <= rGeom + 1e-10, 1, 'last');
    if isempty(idx)
        rStd = 0;
    else
        rStd = standardRadii(idx);
    end
end

%% ========================================================================
% ORIGINAL INNER SOLVER (BISECTION)
%% ========================================================================
function [dFinal, KtF, KtsF, KfF, KfsF, kbF, SeF, iter, q, qs] = ...
    solveSection(r, d0, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, tol)

    % Original special case: no shoulder / no fillet.
    if r <= 0
        q = 0;
        qs = 0;
    else
        % Original Neuber/Peterson notch sensitivity.
        Sut_ksi = Sut / 6.894;
        sqrt_a_in = 0.246 - 3.08e-3*Sut_ksi + 1.51e-5*Sut_ksi^2 - 2.67e-8*Sut_ksi^3;
        sqrt_a_in = max(sqrt_a_in, 1e-6);
        sqrt_a_mm = sqrt_a_in * sqrt(25.4);  % sqrt(in) -> sqrt(mm)

        q = 1 / (1 + sqrt_a_mm / sqrt(r));
        qs = q;
    end

    % Original bracketing around the initial diameter guess.
    dL = max(0.5 * d0, 0.1);
    dU = max(2.0 * d0, dL + 0.5);

    fL = gRes(dL, r, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, q, qs);
    fU = gRes(dU, r, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, q, qs);

    while fL > 0 && fU > 0
        dU = 1.5 * dU;
        fU = gRes(dU, r, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, q, qs);
        if dU > 5000
            error('Could not bracket the root on the upper side.');
        end
    end

    while fL < 0 && fU < 0
        dL = 0.5 * dL;
        if dL < 0.05
            dL = 0.05;
        end
        fL = gRes(dL, r, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, q, qs);
        if dL <= 0.05
            error('Could not bracket the root on the lower side.');
        end
    end

    if sign(fL) == sign(fU)
        error('Root not bracketed. Check inputs.');
    end

    % Original Bisection method and iteration limit.
    iter = 0;
    maxIter = 250;
    while (dU - dL) > tol
        iter = iter + 1;
        dM = 0.5 * (dL + dU);
        fM = gRes(dM, r, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, q, qs);

        if fM > 0
            dL = dM;
        else
            dU = dM;
        end

        if iter >= maxIter
            warning('Maximum inner Bisection iterations reached.');
            break;
        end
    end

    dFinal = dU;

    % Original final stress-concentration and Marin-factor values.
    kbF = sizeFactorKb(dFinal);
    SeF = ka * kbF * Se0;

    if r <= 0
        KtF = 1;
        KtsF = 1;
        KfF = 1;
        KfsF = 1;
    else
        % Peterson shoulder-fillet correlations, simplified for t/r = 1.
        alpha = 2*r / (dFinal + 2*r);
        KtF  = 2.022 - 2.468*alpha + 2.091*alpha^2 - 0.648*alpha^3;
        KtsF = 1.613 - 1.853*alpha + 2.052*alpha^2 - 0.804*alpha^3;
        KfF = 1 + q * (KtF - 1);
        KfsF = 1 + qs * (KtsF - 1);
    end
end

function res = gRes(d, r, Ma, Mm, Ta, Tm, Sut, n, ka, Se0, q, qs)
    % Original Goodman residual with kb updated using current d.
    kb = sizeFactorKb(d);
    Se = ka * kb * Se0;

    if r <= 0
        Kf = 1;
        Kfs = 1;
    else
        % Peterson shoulder-fillet correlations, simplified for t/r = 1.
        alpha = 2*r / (d + 2*r);
        Kt  = 2.022 - 2.468*alpha + 2.091*alpha^2 - 0.648*alpha^3;
        Kts = 1.613 - 1.853*alpha + 2.052*alpha^2 - 0.804*alpha^3;
        Kf = 1 + q * (Kt - 1);
        Kfs = 1 + qs * (Kts - 1);
    end

    termA = sqrt(4*(Kf*Ma)^2 + 3*(Kfs*Ta)^2);
    termM = sqrt(4*(Kf*Mm)^2 + 3*(Kfs*Tm)^2);

    lhs = (16 / (pi * d^3)) * (termA / Se + termM / Sut);
    res = lhs - 1/n;
end

function kb = sizeFactorKb(d)
    if d < 2.79
        d = 2.79;
    end

    if d <= 51
        kb = 1.24 * d^(-0.107);
    else
        kb = 1.51 * d^(-0.157);
    end

    kb = min(kb, 1.0);
end

function dStd = nextStandardDiameter(d)
    % Original common shaft-diameter series in mm.
    stdList = [5 6 8 10 12 14 16 18 20 22 25 28 30 32 35 40 45 50 ...
               55 60 65 70 75 80 85 90 95 100 110 120 130 140 150 160 180 200];

    idx = find(stdList >= d, 1, 'first');
    if isempty(idx)
        dStd = NaN;
    else
        dStd = stdList(idx);
    end
end

%% ========================================================================
% FIXED-WIDTH COMMAND WINDOW REPORT (NO DISP(TABLE), NO CSV EXPORT)
%% ========================================================================
function printFinalReport(shaftID, secName, dRaw, dStd, ...
    rGeom, rUsed, Ma_Nm, Mm_Nm, Ta_Nm, Tm_Nm, ...
    Kt, Kts, Kf, Kfs, q, qs, kb, Se, innerIter, ...
    outerIter, converged, maxDChange, outerTol, innerTol, ...
    radiusRelaxation, OuterLog)

    line = repmat('=',1,108);
    thin = repmat('-',1,108);

    fprintf('\n%s\n', line);
    fprintf('                    TWO-SHAFT FATIGUE DESIGN - FINAL REPORT\n');
    fprintf('%s\n', line);

    %% TABLE 1: FOUR MOST IMPORTANT FINAL RESULTS
    fprintf('\n%s\n', thin);
    fprintf('TABLE 1 - FINAL DESIGN RESULTS\n');
    fprintf('%-8s %-10s %12s %12s %12s %12s\n', ...
        'Shaft','Section','D_calc [mm]','D_std* [mm]','r_calc [mm]','r_std* [mm]');
    fprintf('%s\n', thin);

    for i = 1:numel(secName)
        % Format cells BEFORE adding ANSI style codes so columns stay aligned.
        dStdText = emphasizeStandard(sprintf('%12.0f', dStd(i)));
        rStdText = emphasizeStandard(sprintf('%12.2f', rUsed(i)));

        fprintf('%-8s %-10s %12.2f %s %12.2f %s\n', ...
            char(shaftID(i)), char(secName(i)), dRaw(i), dStdText, ...
            rGeom(i), rStdText);
    end
    fprintf('%s\n', thin);
    fprintf('* Standardized dimensions are emphasized (bold where supported).\n');

    %% TABLE 2: NUMERICAL COMPUTATION / OUTER-LOOP REPORT
    fprintf('\n%s\n', thin);
    fprintf('TABLE 2 - NUMERICAL COMPUTATION REPORT (OUTER LOOP)\n');
    fprintf('Outer-loop iterations completed : %d\n', outerIter);
    fprintf('Convergence status              : %s\n', yesNo(converged));
    fprintf('Diameter convergence tolerance  : %.3f mm\n', outerTol);
    fprintf('Final maximum diameter change   : %.3f mm\n', maxDChange);
    fprintf('Radius relaxation factor        : %.2f\n', radiusRelaxation);
    fprintf('\nOUTER-LOOP CONVERGENCE HISTORY\n');
    fprintf('%8s %22s %24s %14s\n', ...
        'Pass','Max |Delta D| [mm]','Max |r_target-r| [mm]','Converged');
    fprintf('%s\n', thin);

    for k = 1:height(OuterLog)
        if isnan(OuterLog.MaxDiameterChange_mm(k))
            dText = '---';
        else
            dText = sprintf('%.3f', OuterLog.MaxDiameterChange_mm(k));
        end
        fprintf('%8d %22s %24.3f %14s\n', ...
            OuterLog.OuterIteration(k), dText, ...
            OuterLog.MaxRadiusResidual_mm(k), ...
            yesNo(OuterLog.Converged(k)));
    end

    %% TABLE 3: INNER-LOOP / BISECTION REPORT
    fprintf('\n%s\n', thin);
    fprintf('TABLE 3 - INNER-LOOP REPORT (BISECTION METHOD)\n');
    fprintf('Bisection stopping tolerance : %.3f mm\n', innerTol);
    fprintf('Maximum allowed iterations  : 250\n');
    fprintf('%-8s %-10s %16s %24s\n', ...
        'Shaft','Section','r_used [mm]','Bisection iterations');
    fprintf('%s\n', thin);
    for i = 1:numel(secName)
        fprintf('%-8s %-10s %16.2f %24d\n', ...
            char(shaftID(i)), char(secName(i)), rUsed(i), innerIter(i));
    end
    fprintf('%s\n', thin);
    fprintf('Total inner Bisection iterations: %d\n', sum(innerIter));

    %% TABLE 4: APPLIED LOADS (ORIGINAL LOAD TABLE)
    fprintf('\n%s\n', thin);
    fprintf('TABLE 4 - APPLIED LOADS\n');
    fprintf('%-8s %-10s %13s %13s %13s %13s\n', ...
        'Shaft','Section','Ma [N.m]','Mm [N.m]','Ta [N.m]','Tm [N.m]');
    fprintf('%s\n', thin);
    for i = 1:numel(secName)
        fprintf('%-8s %-10s %13.2f %13.2f %13.2f %13.2f\n', ...
            char(shaftID(i)), char(secName(i)), ...
            Ma_Nm(i), Mm_Nm(i), Ta_Nm(i), Tm_Nm(i));
    end

    %% TABLE 5: STRESS-CONCENTRATION / FATIGUE FACTORS (ORIGINAL FACTOR TABLE)
    fprintf('\n%s\n', thin);
    fprintf('TABLE 5 - STRESS-CONCENTRATION AND FATIGUE FACTORS\n');
    fprintf('%-8s %-10s %8s %8s %8s %8s %7s %7s %9s %11s\n', ...
        'Shaft','Section','Kt','Kts','Kf','Kfs','q','qs','kb','Se [MPa]');
    fprintf('%s\n', thin);
    for i = 1:numel(secName)
        fprintf('%-8s %-10s %8.3f %8.3f %8.3f %8.3f %7.3f %7.3f %9.4f %11.1f\n', ...
            char(shaftID(i)), char(secName(i)), Kt(i), Kts(i), ...
            Kf(i), Kfs(i), q(i), qs(i), kb(i), Se(i));
    end

    %% NOTES (RETAINED)
    fprintf('\n%s\n', line);
    fprintf('NOTES\n');
    fprintf('1) D_calc is the converged Bisection diameter from the final standardized-fillet run.\n');
    fprintf('2) D_std is the next standard shaft diameter at or above the reported calculated diameter.\n');
    fprintf('3) r_calc is the geometric shoulder-radius limit: abs(D_left-D_right)/2.\n');
    fprintf('4) The rightmost section of each shaft is assigned r = 0.\n');
    fprintf('5) r_std is the selected preferred fillet radius and does not exceed r_calc.\n');
    fprintf('6) The outer loop stops when the maximum change in calculated diameters is within tolerance.\n');
    fprintf('7) Table 3 reports inner Bisection iterations for the final standardized-fillet run.\n');
    fprintf('8) No output CSV files are generated; this report is printed in the Command Window.\n');
    fprintf('%s\n\n', line);
end

function out = emphasizeStandard(txt)
    % MATLAB R2026b and newer render ANSI bold in the Command Window.
    % Older releases receive asterisks as a visible fallback.
    persistent supportsAnsiBold
    if isempty(supportsAnsiBold)
        try
            supportsAnsiBold = isMATLABReleaseAtLeast('R2026b');
        catch
            supportsAnsiBold = false;
        end
    end

    if supportsAnsiBold
        esc = char(27);
        out = [esc '[1m' txt esc '[0m'];
    else
        % Keep the same 12-character cell width while marking the value.
        out = ['*' sprintf('%10s*', strtrim(txt))];
    end
end

function out = yesNo(tf)
    if tf
        out = 'YES';
    else
        out = 'NO';
    end
end
