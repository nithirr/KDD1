function truelabels=Main_AP_cluster(id,algorithm,nrun2,nconv,pstep)

lam = 0.5;      % damping factor
cut = 3;        % after clustering, drop an cluster with number of samples < cut
splot = 'noplot'; % observing a clustering process when it is on

switch id
  case 1
     sw='3k2lap.txt';               % true number of clusters is 3
  case 2
     sw='5k8close.txt';           % true NC=5
  case 3
     sw='14k10close.txt';         % true NC=14
  case 4
     sw='22k10far.txt';            % true NC=22
  case 5
     sw='ionosphere.txt';        % true NC=2
  case 6
     sw='wine.txt';               % true NC=3
  case 7
     sw='yourdata.txt';
 
  case 11   
     sw='FaceClusteringSimilarities.txt'; nrow = 900; % number of samples
  case 12
     sw='DocumentSummarization.mat'; nrow = 125;
  case 13
     sw='TravelRouting.mat'; nrow = 456; lam = 0.8;
  case 14
     sw='GeneFindingProblem.mat';
     nrow = 75067; cut = 1; nsubset = 3500; % a subset of 75067
  case 15
     sw='yourdata.txt';      
    

  case 21
     sw='yeast.txt';                  % true NC=4
  case 22
     sw='nci60.txt';                 % true NC=8
  case 23
     sw='yourdata.txt'; 

  case 31
     sw='yourdata.txt'; 
end

type = 1;       % 1: Euclidean distances
if id > 20
   type = 2;    % 2: Pearson correlation coefficients
end
simatrix = 0;   % 0: data as input; 1: similarity matrix as input
if id > 10 && id <15
    simatrix = 1;
end
data_load      % loading a data file or similarity matrix


disp(' '); disp([' Clustering is in process  ' ,' please wait ...']);
if algorithm
   tic;
   if simatrix
      [labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(M,type,...
        p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);
   else
      [labels,NCs,labelid,iend,Sp,Slam,NCfixs] = adapt_apcluster(data,type,...
        p,pstep,simatrix,'convits',nconv,'maxits',nrun,'dampfact',lam,splot);
   end
  [NC,Sil,Silmin] = solution_evaluation(data,M,labels,NCs,...
      NCfixs,simatrix,nrow,type,cut);
  trun = toc;
  if id == 12 || id == 13
      NCs = unique(labelid);
  end
  
  % finding an optimal clustering solution
  solution_findK
  
else
    tic;
    if ~simatrix
       M = simatrix_make(data,type,nrow);
    end
    if ~length(p)
        dn = find(M(:,3)>-realmax);
        p = median(M(dn,3));        
    end
    [labels,netsim,iend,unconverged] = apcluster(M,p,'convits',...
        nconv,'maxits',nrun2,'dampfact',lam,splot);
    trun = toc;
    
    solution_findK
end


truek = unique(truelabels);
truek = length(truek);
if truek > 1
    C = valid_external(labels(:,Sid), truelabels);
end
if NCopt == truek
   valid_errorate(labels(:,Sid), truelabels);
end
if id == 12 || id == 13
    for j = 1:length(NCs)
       disp(name{NCs(j)});
    end
end
