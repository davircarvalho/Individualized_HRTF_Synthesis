% test AKtools distance shift 


clear all; clc; close all;
Obj = SOFAload([pwd, '\..\Datasets\CIPIC\subject_003.sofa']);
addpath(genpath([pwd, '\..\Functions']));

sg = Obj.SourcePosition;
ear=[90 0];
offCenter = false;
a = Obj.ReceiverPosition(3);
r_0 = sg(1,3);
Nsh = 100; 
Nsamples = size(Obj.Data.IR, 3);
fs = Obj.Data.SamplingRate;
c = 343;

sg = sg(1,:);
sg(2,:) = sg;
sg(2,3) = .2; % objective

tic
h = AKsphericalHead(sg, ear, offCenter, a, r_0, Nsh, Nsamples, fs, c);
toc

%Get distance variation functions for left and right ear
H_eq = fft(h1)./fft(h0);

startingHRTF = fft(shiftdim(Obj.Data.IR, 2));
%Apply distance variation functions to HRTFs
shifted_H= shiftdim(ifft(H_eq .*startingHRTF), 1);

Obj.Data.IR = shifted_H;

% SOFAsave('20cm.sofa', Obj);

