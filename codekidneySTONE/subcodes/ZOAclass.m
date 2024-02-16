function[abnormal,abnormalOutlineInserted,max_area]=ZOAclass(inp)
inp = uint8(inp);
inp=imresize(inp,[256,256]);
if size(inp,3)>1
    inp=rgb2gray(inp);
end   
sout=imresize(inp,[256,256]);
t0=60;
th=t0+((max(inp(:))+min(inp(:)))./2);
for i=1:1:size(inp,1)
    for j=1:1:size(inp,2)
        if inp(i,j)>th
            sout(i,j)=1;
        else
            sout(i,j)=0;
        end
    end
end

sout=imresize(inp,[256,256]);
t0=60;
th=t0+((max(inp(:))+min(inp(:)))./2);
for i=1:1:size(inp,1)
    for j=1:1:size(inp,2)
        if inp(i,j)>th
            sout(i,j)=1;
        else
            sout(i,j)=0;
        end
    end
end

label=bwlabel(sout);
stats=regionprops(logical(sout),'Solidity','Area','BoundingBox');
density=[stats.Solidity];
area=[stats.Area];
high_dense_area=density>0.6;
max_area=max(area(high_dense_area));
abnormal_label=find(area==max_area);
abnormal=ismember(label,abnormal_label);

if max_area>100
    h = msgbox('Kindney disease present!!','Abnormalities detected','error');
   disp('Kindney disease present');
else if max_area>50
    h = msgbox('kidney disease may occur in Future!!','No abnormalities present now','warn');
    disp('kidney disease may occur in Future');
   
    else
        h = msgbox('No possible symptoms for kidney disease in Future!!','No abnormalities present');
    disp('No possible symptoms for kidney disease in Future');
    end
end

box = stats(abnormal_label);
wantedBox = box.BoundingBox;
dilationAmount = 5;
rad = floor(dilationAmount);
[r,c] = size(abnormal);
filledImage = imfill(abnormal, 'holes');

for i=1:r
   for j=1:c
       x1=i-rad;
       x2=i+rad;
       y1=j-rad;
       y2=j+rad;
       if x1<1
           x1=1;
       end
       if x2>r
           x2=r;
       end
       if y1<1
           y1=1;
       end
       if y2>c
           y2=c;
       end
       erodedImage(i,j) = min(min(filledImage(x1:x2,y1:y2)));
   end
end
abnormalOutline=abnormal;
abnormalOutline(erodedImage)=0;
rgb = inp(:,:,[1 1 1]);
red = rgb(:,:,1);
red(abnormalOutline)=255;
green = rgb(:,:,2);
green(abnormalOutline)=0;
blue = rgb(:,:,3);
blue(abnormalOutline)=0;
abnormalOutlineInserted(:,:,1) = red; 
abnormalOutlineInserted(:,:,2) = green; 
abnormalOutlineInserted(:,:,3) = blue; 