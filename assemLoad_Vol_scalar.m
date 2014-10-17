function L = assemLoad_Vol_scalar(Mesh,h,QuadRule,FHandle)

n = 1/h;
L = zeros((n^2)/4,1);

for i = 1:1/h/2
    for j = 1:1/h/2
        idx = i+(j-1)/h/2;
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);      % bottom left corner vertice of a 2 by 2 block
        
        x = ones(size(QuadRule.x,1),1)*Mesh.Coordinates(ref_no,:)+QuadRule.x*2*h;
        FVal = FHandle(x);
        L(idx) = sum(QuadRule.w.*FVal)*(2*h)^2;
    end
end