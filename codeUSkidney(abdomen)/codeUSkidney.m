clc;
clear all;
close all;
warning off;
addpath ('subcodes');
fprintf('Running Code.m...\n'); 

%% INITIAL ASSUMPTIONS
id = 5;        
algorithm = 0;  
nrun2 = 500;    % max iteration 
nconv = 50;     % convergence condition
pstep = 0.01;
nPop=100;
SearchAgents=30;  
x0=10;
y0=10;
width=700;
height=600;

%% LOADING INPUT DATASET
[I,path]=uigetfile('*.png','select a input image');
str=strcat(path,I);
s=imread(str);
figure;
set(gcf,'position',[x0,y0,width,height]);
subplot(321)
imshow(s);                                                                                                                                                                                                                                                                                                     
title('Input Image');
    
%% PREPROCESSING(Alpha Trimmed Mean Filter)
ma= uint8(255 * mat2gray(s));
data=rgb2gray(ma);
data=im2double(data);
masksize=2;
d=4;
[ro col]=size(data);
temp1=[];
graber=0;
akkumulator=[];

for i=1:ro;
    for j=1:col;
        for m=-masksize:masksize;
            for n=-masksize:masksize;
                if (i+m>0 && i+m<ro && j+n>0 && j+n<col && ...      % To keep indices in limit
                        masksize+m>0 && masksize+m<ro && ...
                        masksize+n>0 && masksize+n<col) 
                    temp1=[temp1 data(i+m,j+n)];               
                end
            end
        end
        temp1=sort(temp1);
        lenth=length(temp1);
        for k=((d/2)-1):(lenth-(d/2))
            akkumulator=[akkumulator temp1(k)];
        end
        akkumulator=sum(akkumulator);
        reformedimage(i,j)=(akkumulator) / (25-d);
        akkumulator=[];
        temp1=[];
    end
end

subplot(322);
imshow(reformedimage);
title('Filtered Image');

%% SUPERPIXEL SEGMENTATION

[L,g] = superpixels(reformedimage,100);
subplot(323);
BW = boundarymask(L);
imshow(imoverlay(reformedimage,BW,'cyan'),'InitialMagnification',67);
title('Pixel Segmentation');
outputImage = zeros(size(ma),'like',ma);
idx = label2idx(L);
numRows = size(ma,1);
numCols = size(ma,2);
for labelVal = 1:g
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = mean(ma(redIdx));
    outputImage(greenIdx) = mean(ma(greenIdx));
    outputImage(blueIdx) = mean(ma(blueIdx));
end    
subplot(324);
imshow(outputImage,'InitialMagnification',67);
title({'Super Pixel Segmented Image' ,'used for Feature Extraction'});

%% affinity propagation clustering 
truelabels=Main_AP_cluster(id,algorithm,nrun2,nconv,pstep);

%% THRIPS OPTIMIZATION
Feature='F1';
Feature1='F10';
[lowerbound,upperbound,dimension,fitness]=fun_info(Feature);
[Best_score,Best_pos,TSO_curve]=tso(SearchAgents,nrun2,lowerbound,upperbound,dimension,fitness);
disp(['The best solution obtained by TSO is : ', num2str(Best_pos)]);
disp(['The best optimal value of the objective funciton found by TSO is : ', num2str(Best_score)]);
[AGPSO1_cg_curve,PSO_cg_curve,IPSO_cg_curve]=Optcomparison(nrun2,SearchAgents);
%% CLASSIFICATION by CNN & Zebra Optimization

[Best_score,Best_pos,ZOA_curve]=ZOA(SearchAgents,nrun2,lowerbound,upperbound,dimension,fitness);
[abnormal,abnormalOutlineInserted,max_area]=ZOAclassUS(s);
disp(['The best solution obtained by ZOA is : ', num2str(Best_pos)]);
disp(['The best optimal value of the objective funciton found by ZOA is : ', num2str(Best_score)]);
[CNVG]=HHO_opt(Feature,nrun2,SearchAgents);
[Convergence_curve]=HPO(nrun2,SearchAgents);
if max_area>50
subplot(325);
imshow(abnormalOutlineInserted);
title('Affected Area','FontSize',12);

subplot(326);
imshow(abnormal);
title('Detected Abnormality','FontSize',12);
   
    else
subplot(325);
imshow(outputImage);
title('No Affected Area detected','FontSize',12);

subplot(326);
imshow(outputImage);
title('No Abnormalities','FontSize',12);
end

figure;
semilogx(TSO_curve,'-r','linewidth',2);
hold on
semilogx(ZOA_curve,'-g','linewidth',2);
hold on
semilogx(AGPSO1_cg_curve,'Color','b','linewidth',2);
hold on
semilogx(PSO_cg_curve,'Color','k','linewidth',2);
hold on
semilogx(IPSO_cg_curve,'Color','m','linewidth',2);
hold on
semilogx(CNVG,'Color','y','linewidth',2);
hold on
semilogx(Convergence_curve,'Color','c','linewidth',2);
title('Objective space')
xlabel('Iterations');
ylabel('Best solution');
legend('TSO','ZOA','AGPSO','PSO','IPSO','HHO','HPO');
axis tight
grid on
box on
InpRes50=[ZOA_curve];
TarRes50=ones(size(InpRes50,1),1);
[Acc,sensitivity,specificity,f1score]= Res50ClassificationUS(InpRes50,TarRes50);
disp('%%% ===== CNN(Resnet 50)Classifier PERFORMANCE=====%%%')
disp([' Accuracy = ' num2str(Acc)])
disp([' Precision = ' num2str(sensitivity)])
disp([' Specificity = ' num2str(specificity)])
disp([' F1 score = ' num2str(f1score)])

[Acc,sensitivity,specificity,f1score]= CNNclassificationUS(InpRes50,TarRes50);
disp('%%% ===== CNN Classifier PERFORMANCE=====%%%')
disp([' Accuracy = ' num2str(Acc)])
disp([' Precision = ' num2str(sensitivity)])
disp([' Specificity = ' num2str(specificity)])
disp([' F1 score = ' num2str(f1score)])
[precision,specificity,Acc,f1score] = DeepLearningUS(InpRes50,TarRes50);

disp('========= Performance Analysis of DeepLearning ========')
disp(['Precision = ' num2str(precision)])
disp(['Specificity = ' num2str(specificity)])
disp(['Accuracy = ' num2str(Acc)])
disp(['F1score = ' num2str(f1score)])

[precision,specificity,Acc,f1score]= ANNUS(InpRes50,TarRes50);
disp('========= Performance Analysis of ANN ========')
disp(['Precision = ' num2str(precision)])
disp(['Specificity = ' num2str(specificity)])
disp(['Accuracy = ' num2str(Acc)])
disp(['F1score = ' num2str(f1score)])
