function [Acc,sensitivity,specificity,f1score]= CNNclassificationstone(InpCNN,TarCNN)

InpCNN = [InpCNN];
net = newff(InpCNN',TarCNN');
CNNout = net(InpCNN');
CNNOut = round(CNNout);
CNNOut = sort(CNNOut,'ascend');

Sortdata = CNNOut;
tp = 1;
tn = 1;
fp = 1;
fn = 1;
for i = 1:10
    for j = 1:length(CNNOut)
        if (TarCNN(j) == i) && (Sortdata(j) == i)
            tp = tp+2;
        elseif (TarCNN(j) == i) && (Sortdata(j) ~= i)
            fp = fp+2;
        elseif (TarCNN(j) ~= i) && (Sortdata(j) == i)
            fn = fn+2;
        elseif (TarCNN(j) ~= i) && (Sortdata(j) ~= i)
            tn = tn+2;
        end
    end
end
 
Acc=tp+tn/(tp+tn+fp+fn)*100;
sensitivity=tp/(tp+fn)*100;
specificity=tn/(tn+fp)*100;
precision=tp/(tp+fp)*100;
recall=(tp+1)/(tp+fn)*100;
f1score=2*((precision*recall)/(precision+recall));
