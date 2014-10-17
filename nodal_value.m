function U = nodal_value(Mesh,U_aux,mark,h,UD1,UD2)

nElements = size(Mesh.Elements,1);
U = zeros(2*nElements,1);
for i = 1:size(mark,2)
    Elements = Mesh.Edge2Elem(mark(i),:);
    if sum(Elements==0)==0
        if abs(Elements(1)-Elements(2))~=1
            U(2*Elements(1)-1) = U_aux(i);
            U(2*Elements(2)-1) = U_aux(i);
        else
            U(2*Elements(1)) = U_aux(i);
            U(2*Elements(2)) = U_aux(i);
        end
    else
        Element = Elements(Elements~=0);
        if abs(Mesh.Edges(mark(i),1)-Mesh.Edges(mark(i),2))~=1
            U(2*Element-1) = U_aux(i);
        else
            U(2*Element) = U_aux(i);
        end
    end
end

for i = 1:1/h/2
    for j = 1:1/h/2
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);      % bottom left corner vertice of a 2 by 2 block
        Dof = Mesh.Vert2Edge(ref_no,ref_no+1/h+1);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no,ref_no+1/h+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element) = UD2(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+1/h+1);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2,ref_no+1/h+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element) = UD2(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+2/h+2+1);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2,ref_no+2/h+2+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element-1) = UD1(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+2/h+2+1);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2+2,ref_no+2/h+2+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element-1) = UD1(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+1/h+1+2);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2+2,ref_no+1/h+1+2],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element) = UD2(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no+2,ref_no+1/h+1+2);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no+2,ref_no+1/h+1+2],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element) = UD2(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no+2,ref_no+1);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no+2,ref_no+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element-1) = UD1(mid_pt);
        end
        Dof = Mesh.Vert2Edge(ref_no,ref_no+1);
        Element = sum(Mesh.Edge2Elem(Dof,:));
        mid_pt = sum(Mesh.Coordinates([ref_no,ref_no+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                U(2*Element-1) = UD1(mid_pt);
        end
    end
end