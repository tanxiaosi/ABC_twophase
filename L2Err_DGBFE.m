function err = L2Err_DGBFE(Mesh,u,QuadRule,FHandle,varargin)
% L2ERR_BFE Discretization error in L2 norm for bilinear finite elements.
 
% Chak Shing Lee
% Texas A&M University 2013

  % Intialize constants
  
  nElements = size(Mesh.Elements,1);
  
  % Precompute shape function values at thhe quedrature points
  
  N = shap_BFE(QuadRule.x);
  grad_N = grad_shap_BFE(QuadRule.x);  

  % Compute discretization error

  err = 0;
  for i = 1:nElements
       
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element mapping  
      
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    P4 = Mesh.Coordinates(vidx(4),:);
    
    x = N(:,1)*P1+N(:,2)*P2+N(:,3)*P3+N(:,4)*P4;
    z1 = P1(1)*grad_N(:,1:2)+P2(1)*grad_N(:,3:4) + ...
         P3(1)*grad_N(:,5:6)+P4(1)*grad_N(:,7:8);
    z2 = P1(2)*grad_N(:,1:2)+P2(2)*grad_N(:,3:4) + ...
         P3(2)*grad_N(:,5:6)+P4(2)*grad_N(:,7:8);
    det_DPhi_K = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));
  
    % Evaluate solutions
      
    u_EX = FHandle(x,varargin{:});
    u_FE = u(i)*N(:,1)+u(i)*N(:,2) + ...
           u(i)*N(:,3)+u(i)*N(:,4);
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*abs(u_EX-u_FE).^2.*det_DPhi_K);
      
  end
  err = sqrt(err);
  
return