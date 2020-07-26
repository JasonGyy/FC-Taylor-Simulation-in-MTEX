clc;
clear;
cs = crystalSymmetry('432'); % Specify the crystal Symmetry
ss = specimenSymmetry('1'); % Specify sample Symmetry
sS = slipSystem.bcc(cs); % Define slip system (change from bcc to fcc if required)
q = 0.0;epsilon = strainTensor(diag([1 -q -(1-q)])); % Plain strain condition 
% Uniaxial condition 
%q = 0.5;epsilon = strainTensor(diag([1 -q -(1-q)])); % Choose any q to mimic different strain paths
% Cross rolling about ND
% q = 0.0;epsilon = strainTensor(diag([-q 1 -(1-q)]));
% For Multiaxial Plain strain compression
% epsilon_1 = strainTensor(diag([1 0 -1])); % compression along z
% epsilon_2= strainTensor(diag([ 0 -1 1])); % compression along y
% epsilon_3= strainTensor(diag([-1 1 0]));  % compression along x

% Select appropriate phi2 Sections
sP = phi2Sections(cs,specimenSymmetry('222'),'phi2',45*degree,'resolution',1*degree);
% Keep low resolution for better match
oriGrid = sP.makeGrid('resolution',1.25*degree);
oriGrid.SS = specimenSymmetry;
% Select appropriate strain tensor
[~,~,W] = calcTaylor(inv(oriGrid)*epsilon,sS.symmetrise)
sP.plot(W.angle./degree)
mtexColorbar;