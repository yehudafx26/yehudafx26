function [fn, vb, xf, ff] = generatePhrases(SR, phr_num)
%PLAYNOTES Generate bowing and fingering control signals.
%
% Outputs:
%   fn : signed bow normal force signal [N]
%        Positive and negative signs are used here to select different strings.
%
%   vb : signed bow velocity signal [m/s]
%        The sign is determined by bow_str.
%
%   xf : finger position signal, normalized by string length [-]
%
%   ff : finger force signal [N]
%        In this version, the finger force is kept constant.
%
% Main modifications:
%   1. Bow velocity is shaped by sinusoidal attack and release.
%   2. Bow normal force is shaped by sinusoidal / raised-cosine transitions.
%   3. Finger force is constant over the whole phrase.
%   4. Finger position transitions use a minimum-jerk curve, which gives
%      smoother and more natural motion than linear interpolation.

    %% Phrase definition
    switch phr_num
        case 1
            notes   = '6655006655332200225500443322110022216622776655';
            bow_str = '1100001100110000110000110011000011101100110011';
            durationSec = 8;

        case 2
            notes   = '12161212565356586865323523212323565356532435232178676763232178765653565653561213231232352123535767656768565356532435232178676763232178765653565678687876565356567868787653561217676545435357656856535652';
            bow_str = '10101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010101010';
            durationSec = 30;

        case 3
            notes   = '333323443300666600112211223311110033232177665555';
            bow_str = '110010110000110000110011001100001100101011001111';
            durationSec = 10;

        case 4
            notes   = '11112222333344445555666677778888';
            bow_str = '11110000111100001111000011110000';
            durationSec = 3;

        case 5
            notes   = '323532316161312312353213616131265321612655115614';
            bow_str = '101010101010101010101010101010101010101010101010';
            durationSec = 8;

        case 6
            notes   = '58653353568';
            bow_str = '10101101010';
            durationSec = 3;

        otherwise
            error('Unsupported phrase number.');
    end

    %% Global signal length
    Nf = round(SR * durationSec);

    %% Musical-position mapping
    % The values are normalized finger positions.
    % Here, 1/5 are treated as open-string positions.
    do_pos  = 0;
    re_pos  = 0.035 / 0.38;
    mi_pos  = 0.063 / 0.38;
    fa_pos  = 0.078 / 0.38;
    sol_pos = 0.093 / 0.38; %#ok<NASGU>
    % sol_pos is kept here for possible later extension.

    %% Control amplitudes
    bowForceAmp    = 0.4;    % Absolute bow normal force amplitude [N]
    bowVelocityAmp = 0.2;    % Absolute bow velocity amplitude [m/s]
    fingerForceAmp = 0.4;    % Constant finger force [N]

    %% Transition ratios
    bowRampRatio     = 0.9; % Ratio of each note used for bow attack/release
    forceRampRatio   = 0.9; % Ratio of each note used for bow-force transition
    fingerMoveRatio  = 0.30; % Ratio of each note used for finger-position movement

    %% Convert character strings to numeric arrays
    melody = notes - '0';
    bow_dir_raw = bow_str - '0';

    num_notes = length(melody);

    % Make the bow-direction string robust to length mismatch.
    % If bow_str is shorter than notes, extend it using the last direction.
    % If bow_str is longer than notes, truncate it.
    if length(bow_dir_raw) < num_notes
        bow_dir_raw(end+1:num_notes) = bow_dir_raw(end);
    elseif length(bow_dir_raw) > num_notes
        bow_dir_raw = bow_dir_raw(1:num_notes);
    end

    %% Pre-compute target finger positions and string directions
    note_pos = zeros(1, num_notes);
    note_string_dir = zeros(1, num_notes);

    last_pos = do_pos;

    for i = 1:num_notes
        m = melody(i);

        switch m
            case 0
                % Rest: keep the previous finger position.
                note_pos(i) = last_pos;
                note_string_dir(i) = 0;

            case 1
                % Open note on string 1.
                note_pos(i) = do_pos;
                note_string_dir(i) = 1;

            case 5
                % Open note on string 2.
                note_pos(i) = do_pos;
                note_string_dir(i) = -1;

            case 2
                note_pos(i) = re_pos;
                note_string_dir(i) = 1;

            case 6
                note_pos(i) = re_pos;
                note_string_dir(i) = -1;

            case 3
                note_pos(i) = mi_pos;
                note_string_dir(i) = 1;

            case 7
                note_pos(i) = mi_pos;
                note_string_dir(i) = -1;

            case 4
                note_pos(i) = fa_pos;
                note_string_dir(i) = 1;

            case 8
                note_pos(i) = fa_pos;
                note_string_dir(i) = -1;

            otherwise
                % Undefined pitch symbol: keep previous finger position.
                note_pos(i) = last_pos;
                note_string_dir(i) = 0;
        end

        if m ~= 0
            last_pos = note_pos(i);
        end
    end

    %% Allocate output vectors
    fn = zeros(Nf, 1);
    vb = zeros(Nf, 1);
    xf = zeros(Nf, 1);

    % Finger force is constant in this version.
    ff = fingerForceAmp * ones(Nf, 1);

    %% Use rounded note boundaries to avoid leaving unused samples at the end
    note_edges = round(linspace(1, Nf + 1, num_notes + 1));

    % Initial finger position.
    % This avoids an artificial slide before the first note.
    currentFingerPos = note_pos(1);

    %% Main synthesis loop
    for i = 1:num_notes
        idx = note_edges(i) : note_edges(i+1) - 1;
        segLen = length(idx);

        if segLen <= 0
            continue;
        end

        m = melody(i);
        isRest = (m == 0);

        % Robust ramp lengths.
        bowRampLen    = getRampLength(segLen, bowRampRatio);
        forceRampLen  = getRampLength(segLen, forceRampRatio);
        fingerMoveLen = max(1, min(segLen, round(fingerMoveRatio * segLen)));

        %% Previous and next note information
        if i > 1
            prev_m = melody(i-1);
            prev_bow_dir = bowDirectionSign(bow_dir_raw(i-1));
            prev_fn_target = note_string_dir(i-1) * bowForceAmp;
        else
            prev_m = 0;
            prev_bow_dir = 0;
            prev_fn_target = 0;
        end

        if i < num_notes
            next_m = melody(i+1);
            next_bow_dir = bowDirectionSign(bow_dir_raw(i+1));
        else
            next_m = 0;
            next_bow_dir = 0;
        end

        %% Current targets
        bow_dir = bowDirectionSign(bow_dir_raw(i));

        if isRest
            targetBowVelocity = 0;
            targetBowForce = 0;
        else
            targetBowVelocity = bow_dir * bowVelocityAmp;
            targetBowForce = note_string_dir(i) * bowForceAmp;
        end

        targetFingerPos = note_pos(i);

        %% A. Bow velocity vb
        % Bow velocity uses sinusoidal attack and release.
        % This replaces linear ramps and avoids sharp slope changes.
        v_seg = zeros(segLen, 1);

        if ~isRest
            v_seg(:) = targetBowVelocity;

            % Start ramp:
            % Use attack if this is the first note, the previous note is a rest,
            % or the bow direction changes.
            needAttack = (i == 1) || (prev_m == 0) || (bow_dir ~= prev_bow_dir);

            % End ramp:
            % Use release if this is the last note, the next note is a rest,
            % or the bow direction changes at the next note.
            needRelease = (i == num_notes) || (next_m == 0) || (bow_dir ~= next_bow_dir);

            if needAttack
                r = sinRampUp(bowRampLen);
                v_seg(1:bowRampLen) = targetBowVelocity * r;
            end

            if needRelease
                r = cosRampDown(bowRampLen);
                v_seg(end-bowRampLen+1:end) = targetBowVelocity * r;
            end
        end

        vb(idx) = v_seg;

        %% B. Bow normal force fn
        % Bow force uses smooth trigonometric transitions.
        % The sign of fn still follows the selected string.
        f_seg = zeros(segLen, 1);

        if ~isRest
            f_seg(:) = targetBowForce;

            % If entering from rest, use a sinusoidal attack from zero.
            if (i == 1) || (prev_m == 0)
                r = sinRampUp(forceRampLen);
                f_seg(1:forceRampLen) = targetBowForce * r;

            % If switching string, use raised-cosine interpolation between
            % the previous signed force and the current signed force.
            elseif prev_fn_target ~= targetBowForce
                f_seg(1:forceRampLen) = raisedCosineInterp( ...
                    prev_fn_target, targetBowForce, forceRampLen);
            end

            % If leaving into rest or reaching the end, use a sinusoidal release.
            if (i == num_notes) || (next_m == 0)
                r = cosRampDown(forceRampLen);
                f_seg(end-forceRampLen+1:end) = targetBowForce * r;
            end
        end

        fn(idx) = f_seg;

        %% C. Finger position xf
        % Finger motion uses a minimum-jerk trajectory.
        % This gives zero velocity and zero acceleration at both endpoints,
        % which is more natural than linear interpolation.
        x_seg = targetFingerPos * ones(segLen, 1);

        if targetFingerPos ~= currentFingerPos
            x_seg(1:fingerMoveLen) = minimumJerkInterp( ...
                currentFingerPos, targetFingerPos, fingerMoveLen);
        end

        xf(idx) = x_seg;

        % Store the target position for the next transition.
        currentFingerPos = targetFingerPos;
    end
    fn = fn(:);
    vb = vb(:);
    xf = xf(:);
    ff = ff(:);
end

%% ------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------

function s = bowDirectionSign(rawValue)
%BOWDIRECTIONSIGN Convert bow direction symbol to signed direction.
%
% rawValue = 1 gives +1.
% rawValue = 0 gives -1.

    if rawValue == 1
        s = 1;
    else
        s = -1;
    end
end

function n = getRampLength(segLen, ratio)
%GETRAMPLENGTH Compute a safe ramp length for a note segment.
%
% The ramp length is limited to half of the segment length so that attack
% and release do not strongly overlap.

    if segLen <= 2
        n = 1;
    else
        n = round(ratio * segLen);
        n = max(1, n);
        n = min(n, floor(segLen / 2));
    end
end

function y = sinRampUp(n)
%SINRAMPUP Smooth attack curve from 0 to 1.
%
% y(1)   = 0
% y(end) = 1

    if n <= 1
        y = 1;
        return;
    end

    theta = linspace(0, pi/2, n).';
    y = sin(theta);
end

function y = cosRampDown(n)
%COSRAMPDOWN Smooth release curve from 1 to 0.
%
% y(1)   = 1
% y(end) = 0

    if n <= 1
        y = 0;
        return;
    end

    theta = linspace(0, pi/2, n).';
    y = cos(theta);
end

function y = raisedCosineInterp(a, b, n)
%RAISEDCOSINEINTERP Smooth trigonometric interpolation from a to b.
%
% This curve has zero slope at both endpoints, which is useful for avoiding
% abrupt force changes when switching strings.

    if n <= 1
        y = b;
        return;
    end

    t = linspace(0, 1, n).';
    w = 0.5 - 0.5 * cos(pi * t);
    y = a + (b - a) * w;
end

function y = minimumJerkInterp(a, b, n)
%MINIMUMJERKINTERP Minimum-jerk interpolation from a to b.
%
% The polynomial
%
%     w(t) = 10 t^3 - 15 t^4 + 6 t^5
%
% has zero velocity and zero acceleration at t = 0 and t = 1.
% It is often used as a simple model of natural human motion.

    if n <= 1
        y = b;
        return;
    end

    t = linspace(0, 1, n).';
    w = 10*t.^3 - 15*t.^4 + 6*t.^5;
    y = a + (b - a) * w;
end
