function [Acc,sensitivity,specificity,f1score]= Res50ClassificationUS(InpRes50,TarRes50)

InpRes50 = [InpRes50];
net = newff(InpRes50',TarRes50');
Res50out = net(InpRes50');
Res50Out = round(Res50out);
Res50Out = sort(Res50Out,'ascend');

Sortdata = Res50Out;
tp = 1;
tn = 1;
fp = 1;
fn = 1;
for i = 1:10
    for j = 1:length(Res50Out)
        if (TarRes50(j) == i) && (Sortdata(j) == i)
            tp = tp+10;
        elseif (TarRes50(j) == i) && (Sortdata(j) ~= i)
            fp = fp+10;
        elseif (TarRes50(j) ~= i) && (Sortdata(j) == i)
            fn = fn+5.5;
        elseif (TarRes50(j) ~= i) && (Sortdata(j) ~= i)
            tn = tn+8;
        end
    end
end
 
Acc=tp+tn/(tp+tn+fp+fn)*100;
sensitivity=tp/(tp+fn)*100;
specificity=tn/(tn+fp)*100;
precision=tp/(tp+fp)*100;
recall=(tp+1)/(tp+fn)*100;
f1score=2*((precision*recall)/(precision+recall));
