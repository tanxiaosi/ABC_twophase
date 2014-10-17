function [B,mark] = assemMat_Inn2(Mesh,h)

% nCoordinates = size(Mesh.Coordinates,1);
n=1/h;
B = sparse(n^2/4,n^2-2*n);
mark = zeros(1,n^2-2*n);
cnt = 1;
for i = 1:1/h/2
    for j = 1:1/h/2
        idx = i+(j-1)/h/2;
        ref_no = 1+2*(j-1)+(2/h+2)*(i-1);      % bottom left corner vertice of a 2 by 2 block
        Dof = Mesh.Vert2Edge(ref_no,ref_no+1/h+1);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = -h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=-h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+1/h+1);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = -h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=-h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2,ref_no+2/h+2+1);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+2/h+2+1);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no+2/h+2+2,ref_no+1/h+1+2);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no+2,ref_no+1/h+1+2);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no+2,ref_no+1);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = -h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=-h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
        Dof = Mesh.Vert2Edge(ref_no,ref_no+1);
        if Mesh.BdFlags(Dof) >=0
            if sum(Dof==mark)==0
                B(idx,cnt) = -h;
                mark(1,cnt) = Dof;
                cnt = cnt+1;
            else
                B(idx,Dof==mark)=-h;
            end
%         else
%             mark(1,cnt) = Dof;
%             cnt = cnt+1;
        end
    end
end