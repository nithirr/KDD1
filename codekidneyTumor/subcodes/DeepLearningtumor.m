function  [precision,specificity,Acc,f1score] = DeepLearningtumor(InpAnn,TarAnn)

outputnum = 1; 
hiddennum = 8; 
inputnum = size(InpAnn,2); 
inputdata = InpAnn;
outputdata = TarAnn;
dbn = randDBN([inputnum, hiddennum, outputnum]); 
dbn = pretrainDBN( dbn, inputdata );  
dbn = trainDBN( dbn, inputdata, outputdata );
label = v2h( dbn,inputdata);
Sortdata = label;
 
 tp = 1;
 tn = 1;
 fp = 1;
 fn = 1;
  for i = 1:10
     for j = 1:length(Sortdata)
         if (TarAnn(j) == i) && (Sortdata(j) == i)
             tp = tp+1;        
         elseif (TarAnn(j) == i) && (Sortdata(j) ~= i)
             fp = fp+1;
        elseif (TarAnn(j) ~= i) && (Sortdata(j) == i)
             fn = fn+1;
         elseif (TarAnn(j) ~= i) && (Sortdata(j) ~= i)
             tn = tn+1;
         end
     end
 end
  
 Acc = (tp+tn)/(tp+tn+fp+fn)*100;
precision=tp/(tp+fp)*100;
recall=(tp+1)/(tp+fn)*100;
f1score=2*((precision*recall)/(precision+recall));
specificity=tn/(tn+fp)*100;