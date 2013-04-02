% function descriptorArray = getDescriptor(image,corners,BINNUM,PATCHDIM)
% 
% [rowNum, colNum] = size(corners);
% [imageHeight, imageWidth] = size(image);
% descriptorArray = cell(rowNum,1);
% 
% for n=1:rowNum
%     coord = corners(n,:);
%     x=coord(1);
%     y=coord(2);
%     descriptor = zeros(BINNUM+1,1);
%     gap = 1/BINNUM;
% 
% 
%     if (x-PATCHDIM>0 && x+PATCHDIM <= imageWidth && y-PATCHDIM>0 && y+PATCHDIM <= imageHeight)
%         for i=x-PATCHDIM:x+PATCHDIM
%             for j=y-PATCHDIM:y+PATCHDIM
%                 bin = uint32(image(j,i)/gap)+1;
%                 descriptor(bin)=descriptor(bin)+1;
%             end
%         end
% 
%         descriptorArray{n}=descriptor;
%     else
%         descriptorArray{n}=zeros(BINNUM+1,1);
% 
%     end
% 
% end
% 
% end


function descriptorArray = getDescriptor(image,corners,BINNUM,PATCHDIM)

rowNum = size(corners,1);
[imageHeight, imageWidth] = size(image);
descriptorArray = cell(rowNum,1);

for n=1:rowNum
    coord = corners(n,:);
    x=coord(1);
    y=coord(2);
    descriptor = zeros((2*PATCHDIM+1)^2,1);
    
    index=1;
    if (x-PATCHDIM>0 && x+PATCHDIM <= imageWidth && y-PATCHDIM>0 && y+PATCHDIM <= imageHeight)
        for i=x-PATCHDIM:x+PATCHDIM
            for j=y-PATCHDIM:y+PATCHDIM
                descriptor(index)=image(j,i);
                index=index+1;
            end
        end
    end
 
    descriptorArray{n}=descriptor;
    
end

end