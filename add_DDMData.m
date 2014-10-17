function New_Mesh = add_DDMData(Old_Mesh,Hx,Hy,Nx_subdom,Ny_subdom,gp)

New_Mesh = Old_Mesh;

% Initialize constants

nEdges = size(Old_Mesh.Edges,1);
% nCoordinates_Old = nCoordinates-size(Old_Mesh.Elements,1)/3;
Hx_subdom = Hx/Nx_subdom;
Hy_subdom = Hy/Ny_subdom;

% Allocate memory

New_Mesh.DDMData = sparse(1,nEdges);

for i = 1:nEdges
    
    if Old_Mesh.BdFlags(i) >= 0
        
        % Extract edge data
        
        EdgeNo = Old_Mesh.Edges(i,:);
        Edge = (Old_Mesh.Coordinates(EdgeNo(1),:)+Old_Mesh.Coordinates(EdgeNo(2),:))/2;
        
        % Mark the Dofs of vector funciton in non-overlapping subdomain
        
        for j = 1:Nx_subdom*Ny_subdom
            xy = [ceil(j/Ny_subdom),mod(j,Ny_subdom)];
            if xy(2)==0
                xy(2)=Ny_subdom;
            end
            if abs(Edge(2)-(xy(2)-1)*Hy_subdom)<= 10*eps && Edge(1) < xy(1)*Hx_subdom && Edge(1) > (xy(1)-1)*Hx_subdom
                New_Mesh.DDMData(1,i)=gp(j,1);
                break;
            elseif abs(Edge(2)-xy(2)*Hy_subdom)<= 10*eps && Edge(1) < xy(1)*Hx_subdom && Edge(1) > (xy(1)-1)*Hx_subdom
                New_Mesh.DDMData(1,i)=gp(j,2);
                break;
            elseif abs(Edge(1)-(xy(1)-1)*Hx_subdom)<= 10*eps && Edge(2) < xy(2)*Hy_subdom && Edge(2) > (xy(2)-1)*Hy_subdom
                New_Mesh.DDMData(1,i)=gp(j,3);
                break;
            elseif abs(Edge(1)-xy(1)*Hy_subdom)<= 10*eps && Edge(2) < xy(2)*Hy_subdom && Edge(2) > (xy(2)-1)*Hy_subdom
                
                New_Mesh.DDMData(1,i)=gp(j,4);
                break;
            elseif Edge(1) < xy(1)*Hx_subdom && Edge(1) > (xy(1)-1)*Hx_subdom && Edge(2) < xy(2)*Hy_subdom && Edge(2) > (xy(2)-1)*Hy_subdom
                New_Mesh.DDMData(1,i)=Ny_subdom*(Nx_subdom-1)+Nx_subdom*(Ny_subdom-1)+j;
                break;
            end
        end
    end
end

% for i = 1:Nx_subdom*Ny_subdom
%     
%     
%     
%     xy = [ceil(i/Ny_subdom),mod(i,Ny_subdom)];
%     if xy(2)==0
%         xy(2)=Ny_subdom;
%     end
%     
%     Mark the Dofs of vector funciton in overlapping subdomain
%     
%     for j = 1:nEdges
%         
%         if Old_Mesh.BdFlags(j) >= 0
%             
%             Extract edge data
%             
%             EdgeNo = Old_Mesh.Edges(j,:);
%             Edge = (Old_Mesh.Coordinates(EdgeNo(1),:)+Old_Mesh.Coordinates(EdgeNo(2),:))/2;
%             
%             if Edge(2)==(xy(2)-1)*Hy_subdom-extend && Edge(1) < xy(1)*Hx_subdom+extend && Edge(1) > (xy(1)-1)*Hx_subdom-extend
%                 New_Mesh.DDMData(i,j)=2;
%             elseif Edge(2)==xy(2)*Hy_subdom+extend && Edge(1) < xy(1)*Hx_subdom+extend && Edge(1) > (xy(1)-1)*Hx_subdom-extend
%                 New_Mesh.DDMData(i,j)=3;
%             elseif Edge(1)==(xy(1)-1)*Hx_subdom-extend && Edge(2) < xy(2)*Hy_subdom+extend && Edge(2) > (xy(2)-1)*Hy_subdom-extend
%                 New_Mesh.DDMData(i,j)=4;
%             elseif Edge(1)==xy(1)*Hy_subdom+extend && Edge(2) < xy(2)*Hy_subdom+extend && Edge(2) > (xy(2)-1)*Hy_subdom-extend
%                 New_Mesh.DDMData(i,j)=5;
%             elseif Edge(1) < xy(1)*Hx_subdom+extend && Edge(1) > (xy(1)-1)*Hx_subdom-extend && Edge(2) < xy(2)*Hy_subdom+extend && Edge(2) > (xy(2)-1)*Hy_subdom-extend
%                 New_Mesh.DDMData(i,j)=1;
%             end
%         end
%     end
% end
