
%This script for simualating pseudo cross rolling or two step cross rolling textures in f.c.c and b.c.c
%materials. Pseudo cross rolling refers here to change in rolling direction by 90
%degrees after certain number of user defined pass. 

clc;
clear all;

cs = crystalSymmetry('432');
ss = specimenSymmetry('triclinic');

% define a family of slip systems
sS = slipSystem.fcc(cs); %change from slipSystem.fcc to Slipsystems.bcc if b.c.c cross rolling is of interest

% some plane strain
q = 0.0;

%Define any texture by using appropriate Euler angles

cube = orientation('Euler',[0,0,0]*degree,cs,ss);
odf_Cube = unimodalODF(cube,'halfwidth',15*degree);
ori1 = calcOrientations(odf_Cube,1000);
ODF=calcODF(ori1);
plot(ODF,'phi2',[0 45 65]*degree,'contour','antipodal','linewidth',3);

%In this case we have selected a random texture of 1000 single orientations

% ori1 = orientation.rand(1000,cs);


%define the strain matrix for rolling along RD1
epsilon = strainTensor(diag([1 -q -(1-q)]));

%define the strain matrix for rolling along RD2
epsilon_cross= strainTensor(diag([-q 1 -(1-q)]));

%define the strain increment per pass
inc=0.1;

%define the numebr of iterations. e.g. 20 iterations with 0.1 increment
%corresponds to true strain of 2.0
numIter = 20;

%Perform the actual calculations

for sas=1:numIter

    if sas <= 14 %first 14 steps viz. till strain of 1.4 rolling is carried out along RD1
    
    [M,~,W] = calcTaylor(inv(ori1) * epsilon .*inc, sS.symmetrise,'silent');
    ori1 = ori1 .* orientation(-W);
    
    end
    
    if sas > 14 && sas <=20 %last 6 steps viz. a strain of 0.6 (total cumulative strain 2.0) rolling is carried out along RD2
    [M,~,W] = calcTaylor(inv(ori1) * epsilon_cross .*inc, sS.symmetrise,'silent');
    ori1 = ori1 .* orientation(-W);
    end
    
end

%plot the texture results uisng ODF and pole figures after simualations
ODF1=calcODF(ori1);
rot = rotation('axis', zvector, 'angle',0*degree);
ODF1=rotate(ODF1,rot);
plot(ODF1,'phi2',[0 45 65]*degree,'contour','antipodal','linewidth',3);
setColorRange('equal');
mtexColorbar ('FontSize',25,'Fontweight','bold');
figure ()
plotPDF(ori1,Miller({1,1,1},cs),'contour','antipodal','complete','linewidth',3,'complete','upper');
mtexColorbar ('FontSize',25,'Fontweight','bold');
setColorRange('equal') % set equal color range for all subplots
annotate([xvector, yvector], 'label', {'RD','TD'}, 'BackgroundColor', 'w',...
    'FitBoxToText','on','FontSize',15,'LineStyle','none','Fontname','Times New Roman','Fontweight','bold');
