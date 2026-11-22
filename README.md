# Fixed-point-FFT
%% 1024-point radix-2 FFT with data-dependent fixed-point word lengths
% Requires Fixed-Point Designer (fi). Each stage measures the actual
% butterfly results and chooses the smallest signed word length that holds
% them at the original fractional precision. This is block-floating style:
% the entire stage shares one word length, selected from that stage's peak.

clear; close all; clc;

N = 1024;
inputWordLength = 16;
inputFractionLength = 15;
twiddleWordLength = 16;
twiddleFractionLength = 15;
numStages = log2(N);

assert(numStages == fix(numStages), 'N must be a power of two.');

% Demonstration input: complex sinusoid plus a second, smaller tone.
% Replace this block with your own real or complex input vector as needed.
n = 0:N-1;
x = 0.70*cos(2*pi*37*n/N) + 1i*0.55*sin(2*pi*91*n/N);
x = x(:);

% Signed Q1.15 input. The real and imaginary parts each use 16 bits.
Tinput = numerictype(1, inputWordLength, inputFractionLength);
Ttwiddle = numerictype(1, twiddleWordLength, twiddleFractionLength);
F = fimath('RoundingMethod', 'Nearest', ...
           'OverflowAction', 'Saturate', ...
           'ProductMode', 'FullPrecision', ...
           'SumMode', 'FullPrecision');
x = fi(x, Tinput, F);

% Bit-reverse the input to use an iterative decimation-in-time FFT.
numBits = numStages;
indices = uint32((0:N-1).');
reversed = zeros(N, 1, 'uint32');
for b = 1:numBits
    reversed = bitor(bitshift(reversed, 1), bitand(bitshift(indices, -(b-1)), 1));
end
y = x(double(reversed) + 1);

stageWordLengths = zeros(1, numStages + 1);
stageWordLengths(1) = inputWordLength;
stageFractionLengths = inputFractionLength * ones(1, numStages + 1);

for stage = 1:numStages
    span = 2^stage;
    halfSpan = span / 2;
    % Calculate unquantized-to-stage-width butterfly results first.
    upperRAll = zeros(N, 1);
    upperIAll = zeros(N, 1);
    lowerRAll = zeros(N, 1);
    lowerIAll = zeros(N, 1);

    for blockStart = 1:span:N
        for j = 0:halfSpan-1
            upperIndex = blockStart + j;
            lowerIndex = upperIndex + halfSpan;

            % W_span^j = exp(-2*pi*i*j/span), quantized to signed Q1.15.
            angle = -2*pi*j/span;
            wr = fi(cos(angle), Ttwiddle, F);
            wi = fi(sin(angle), Ttwiddle, F);

            upperValue = y(upperIndex);
            ar = real(upperValue);
            ai = imag(upperValue);

            lowerValue = y(lowerIndex);
            br = real(lowerValue);
            bi = imag(lowerValue);

            % Complex twiddle multiply. Full-precision intermediate products
            % and sums avoid accidental narrowing before stage assignment.
            tr = br * wr - bi * wi;
            ti = br * wi + bi * wr;

            upperR = ar + tr;
            upperI = ai + ti;
            lowerR = ar - tr;
            lowerI = ai - ti;

            upperRAll(upperIndex) = double(upperR);
            upperIAll(upperIndex) = double(upperI);
            lowerRAll(lowerIndex) = double(lowerR);
            lowerIAll(lowerIndex) = double(lowerI);
        end
    end

    % Find minimum signed width for the measured peak, with fixed FL.
    peakMagnitude = max(abs([upperRAll; upperIAll; lowerRAll; lowerIAll]));
    if peakMagnitude == 0
        requiredIntegerBits = 1; % signed zero-only stage still needs a sign bit
    else
        % Signed range with I integer bits is [-2^(I-1), 2^(I-1)).
        requiredIntegerBits = max(1, ceil(log2(peakMagnitude + 2^(-inputFractionLength))) + 1);
        while peakMagnitude >= 2^(requiredIntegerBits - 1)
            requiredIntegerBits = requiredIntegerBits + 1;
        end
    end
    outWordLength = requiredIntegerBits + inputFractionLength;
    stageWordLengths(stage + 1) = outWordLength;
    Tout = numerictype(1, outWordLength, inputFractionLength);

    % Quantize once to the selected stage format and store as fi values.
    % Butterfly outputs occupy the upper and lower halves of each block.
    stageR = zeros(N, 1);
    stageI = zeros(N, 1);
    for blockStart = 1:span:N
        for j = 0:halfSpan-1
            upperIndex = blockStart + j;
            lowerIndex = upperIndex + halfSpan;
            stageR(upperIndex) = upperRAll(upperIndex);
            stageI(upperIndex) = upperIAll(upperIndex);
            stageR(lowerIndex) = lowerRAll(lowerIndex);
            stageI(lowerIndex) = lowerIAll(lowerIndex);
        end
    end
    y = complex(fi(stageR, Tout, F), fi(stageI, Tout, F));
end

% Plot actual selected width of each stage's stored data type.
figure('Color', 'w');
stairs(0:numStages, stageWordLengths, 'o-', 'LineWidth', 1.8, ...
       'MarkerFaceColor', [0.1 0.45 0.8]);
grid on;
xlim([0 numStages]);
ylim([inputWordLength-1 stageWordLengths(end)+1]);
xticks(0:numStages);
yticks(inputWordLength:stageWordLengths(end));
xlabel('FFT stage (0 = input)');
ylabel('Stored word length (bits)');
title('1024-point radix-2 FFT: data-dependent bit growth');
for k = 1:numel(stageWordLengths)
    text(k-1, stageWordLengths(k)+0.12, sprintf('%d bits', stageWordLengths(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

fprintf('FFT size: %d points (%d radix-2 stages)\n', N, numStages);
fprintf('Input: %d-bit signed fixed point with %d fractional bits\n', ...
        inputWordLength, inputFractionLength);
fprintf('Data-dependent stage widths (input through output): %s bits\n', ...
        mat2str(stageWordLengths));
fprintf('Output stage word length: %d bits (%d fractional bits)\n', ...
        stageWordLengths(end), inputFractionLength);
fprintf('Output numeric range is approximately [-2^%d, 2^%d).\n', ...
        stageWordLengths(end)-inputFractionLength-1, ...
        stageWordLengths(end)-inputFractionLength-1);

% y contains the unnormalized, bit-true fixed-point FFT result.
