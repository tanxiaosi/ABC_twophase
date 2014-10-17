function g = assemLoad_bnd2(Mesh,h,UD1,UD2)

% nCoordinates = size(Mesh.Coordinates,1);
n=1/h;
g = zeros(n^2/4,1);
for i = 1:1/h/2
    for j = 1:1/h/2
        idx = i+(j-1)/h/2;
        g(idx) = 0;
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);      % bottom left corner vertice of a 2 by 2 block
        Dof = Mesh.Vert2Edge(ref_no,ref_no+1/h+1);
        mid_pt = sum(Mesh.Coordinates([ref_no,ref_no+1/h+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)-UD2(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+1/h+1);
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2,ref_no+1/h+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)-UD2(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+2/h+2+1);
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2,ref_no+2/h+2+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)+UD1(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+2/h+2+1);
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2+2,ref_no+2/h+2+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)+UD1(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+1/h+1+2);
        mid_pt = sum(Mesh.Coordinates([ref_no+2/h+2+2,ref_no+1/h+1+2],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)+UD2(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no+2,ref_no+1/h+1+2);
        mid_pt = sum(Mesh.Coordinates([ref_no+2,ref_no+1/h+1+2],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)+UD2(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no+2,ref_no+1);
        mid_pt = sum(Mesh.Coordinates([ref_no+2,ref_no+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)-UD1(mid_pt)*h;
        end
        Dof = Mesh.Vert2Edge(ref_no,ref_no+1);
        mid_pt = sum(Mesh.Coordinates([ref_no,ref_no+1],:))/2;
        if Mesh.BdFlags(Dof) <0
                g(idx)=g(idx)-UD1(mid_pt)*h;
        end
    end
end