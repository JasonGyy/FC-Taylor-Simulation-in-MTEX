clc;
clear;

cs = crystalSymmetry('432');
ss = specimenSymmetry('1');

sS = slipSystem.fcc(cs);


   copper = orientation('Euler',[10,0,0]*degree,cs,ss);
   odf_copper = unimodalODF(copper,'halfwidth',15*degree);
   ori1 = calcOrientations(odf_copper,500);
   %ori1 = orientation.rand(1000,cs);
  
%    plot(odf_copper,'phi2',[0 45 65]*degree,'contour','antipodal','linewidth',2);
 plot(odf_copper,'phi2',[0 45 65]*degree)
 
   epsilon= 2* strainTensor([0 0.0 1.0; 0.0 0 0;0 0 0]);
   numIter = 10;
   progress(0,numIter);
 
 for sas=1:numIter

  % compute the Taylor factors and the orientation gradients
  [M,b,W] = calcTaylor(inv(ori1)*epsilon./numIter, sS.symmetrise);

  % rotate the individual orientations
  ori1 = ori1 .* orientation(-W);
  progress(sas,numIter);
 end

ODF=calcODF(ori1);
% plot(ODF,'phi2',[0 45]*degree,'contour','antipodal','linewidth',2);
% figure ();
% plotPDF(ori1,Miller({1,1,1},cs),'contourf','complete','upper');

 figure()

 rot = rotation('axis', yvector, 'angle',35*degree);
 ODF1=rotate(ODF,rot);
 plot(ODF1,'phi2',[0 45 65]*degree,'contour','antipodal','linewidth',2);
 
  plot(ODF1,'phi2',[0]*degree)
  plot(ODF1,'phi2',[45]*degree)
  plot(ODF1,'phi2',[65]*degree)
 
%  figure()
%  plotPDF(ODF1,Miller({1,1,1},cs),'contourf','complete','upper');