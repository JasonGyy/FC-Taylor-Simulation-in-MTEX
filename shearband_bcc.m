clc;
clear;

cs = crystalSymmetry('432');
ss = specimenSymmetry('2');

sS = slipSystem.bcc(cs);

% % gamma fiber
phi1 = linspace(0*degree,90*degree);
Phi =  55*degree;
phi2 = 45*degree;
ori = orientation('euler',phi1,Phi,phi2,cs);

odf_gammafiber = unimodalODF(ori,'halfwidth',15*degree);
ori1 = calcOrientations(odf_gammafiber,500);
plot(odf_gammafiber,'phi2',[0 45]*degree,'contour','antipodal','linewidth',2);

 epsilon= 2* strainTensor([0 0.0 1.0; 0.0 0 0;0 0 0]);
 numIter = 20;
 progress(0,numIter);
 
 for sas=1:numIter

  % compute the Taylor factors and the orientation gradients
  [M,~,W] = calcTaylor(inv(ori1) * epsilon ./ numIter, sS.symmetrise,'silent');

  % rotate the individual orientations
  ori1 = ori1 .* orientation(-W);
  progress(sas,numIter);
 end

ODF=calcODF(ori1);
plot(ODF,'phi2',[0 45]*degree,'contour','antipodal','linewidth',2);
figure ();
%plotPDF(ori1,Miller({1,1,1},cs),'contourf','complete','upper');

%  figure()

 rot = rotation('axis', yvector, 'angle',40*degree);
 ODF1=rotate(ODF,rot);
 plot(ODF1,'phi2',[0 45]*degree,'contour','antipodal','linewidth',2);
%  figure()
%  plotPDF(ODF1,Miller({1,1,1},cs),'contourf','complete','upper');
