%% YEHU GESTURE EXAMPLES
% Champ C. Darabundit and Zhen Zhang, CAML 2026

clear;
close all;
addpath('../matlab')
addpath('../cpp/release/mex')

%% EXAMPLE CONTROL
runCpp = false; % run C++ code, build me first!
audioWrite = false; % Audio write flag, if false plays audio
aniPlot = false;    % Animation flag, only with Matlab engine

% Gesture selection
%   1 : Open string and press string
%   2 : A scale
%   3 : Finger up and down
%   4 : A Teochew scale
%   5 : Slide + vibrato
gestureNum = 5;

%% SIMULATION PARAMETERS
SR = 44100;
t = 2;
Nf = ceil (SR * t);

%% GENERATE CONTROL SIGNALS

[fn, vb, xf, ff] = generateGestures(gestureNum, Nf, SR);

%% NUMERICAL SIMULATION
[fStringParams, cStringParams, bodyParams, fingerParams, ...
  bowParams, frictionParams] = yehuParams();

if runCpp
  [u, Fbri] = yehuMex(vb, fn, xf, ff, fStringParams, cStringParams, ...
    bodyParams, fingerParams, bowParams, frictionParams, SR);
else
  [u, Fbri] = yehu(vb, fn, xf, ff, fStringParams, cStringParams, ...
    bodyParams, fingerParams, bowParams, frictionParams, SR, gestureNum, aniPlot);
end

%% RADIATION PROCESSING
record = load("IIR_record.mat");
Bmr = record.Bm;
Amr = record.Am;
FIRr = record.FIR;
outr = parfilt(Bmr ,Amr ,FIRr ,-Fbri);
outr = outr / max(abs(outr)) * 0.8;

player = load("IIR_player.mat");
Bmp = player.Bm;
Amp = player.Am;
FIRp = player.FIR;
outp = parfilt(Bmp ,Amp ,FIRp ,-Fbri);
outp = outp / max(abs(outp)) * 0.5;

%% AUDIO WRITE
if audioWrite
  if ~exist('../samples', 'dir')
    mkdir('../samples')
  end
  filename_u = sprintf('../samples/gesture_%d_u.wav', gestureNum);
  filename_f = sprintf('../samples/gesture_%d_F.wav', gestureNum);
  filename_outr = sprintf('../samples/gesture_%d_pos1_out.wav', gestureNum);
  filename_outp = sprintf('../samples/gesture_%d_pos2_out.wav', gestureNum);

  audiowrite(filename_u, u, SR);
  audiowrite(filename_f, Fbri, SR);
  audiowrite(filename_outr, outr, SR);
  audiowrite(filename_outp, outp, SR);
else
  soundsc( outr, SR )
  pause(length(fn)/SR)
end