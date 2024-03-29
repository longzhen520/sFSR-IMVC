function [X,S,Convergence_curve] = FSR_IMVC(X,Omega,T, opts)
%-----------------------------------------------------------------------------------------------------
% sample missing
% Model: min sum ||Hv||_{F}+lambda ||E||_{2,1}+beta* Mrank(Y)
%      
%        s.t.  Xv(Omega)=Tv(Omega),Xv=Pv*Hv+Ex,Pv'*Pv=I; Hv=Hv*Zv+Eh
%            Y=Z,Z=(Z1,...,ZN)
%             E=[E1;...;EN]
%
% 
%--------------------------------INPUTs----------------------------------------------------------



rho1=1e-3;
rho2=1e-3;
rho3=1e-3;


  rX= opts.rX;   
% init. variables 
sizeD=size(X{1});
% N    = prod(sizeD);
% R_mera=opts.R;
K=opts.K; % subspce dimension
miu   =opts.miu;
%%  check and initilization
if ~exist('opts','var'); opts=[]; end
if isfield(opts,'maxIter'); maxIter = opts.maxIter; else maxIter =50; end
if isfield(opts,'tol');  tol = opts.tol; else tol = 1e-6; end
if isfield(opts,'lambda'); lambda = opts.lambda; else lambda =0.01; end
if isfield(opts,'beta');  beta = opts.beta; else beta = 0.01; end
if isfield(opts,'gamma');  gamma = opts.gamma; else gamma = 0.01; end
if isfield(opts,'I');  I = opts.I;  end

% if isfield(opts,'beta'); beta = opts.beta; else beta(1)=0.01; beta(2)=0.001; end
if isfield(opts,'isCont')
    isCont = opts.isCont; contpar = opts.contpar; contfactor = opts.contfactor;
 else
    isCont = true; contpar = 0.95; contfactor = 2;
end

V = length(X); 
N = size(X{1},2); %sample number


for v=1:V
    sX =size(X{v});
    X{v}=ones(sX).*mean(X{v}(:));
    ind_1 = find(Omega(:,v) == 1);
    X{v}(:,ind_1)=T{v}(:,ind_1);
    Ex{v}=zeros(sX);
    Lambda1{v}=zeros(sX);
    H{v} = zeros(K,N);  % Low-rank coeffs
    P{v} = zeros(size(X{v},1),K);     % Low-rank dictionaries
    
    Eh{v}=zeros(size(H{v}));
    Lambda2{v} = zeros(size(H{v}));
    
    Z{v} = zeros(N,N); 
    y{v}=zeros(N,N);
    lambda3{v}=zeros(N,N);

end
Y=zeros(N,N,V);
Lambda3=Y;
% Y=zeros(rX);
% Lambda3=zeros(rX);
    


 




% gamma    = 2;
display  = 1;
tic;
for iter = 1 : maxIter
    
%     fprintf('\n*****************************iter: %d ******************************\n', iter');
     Z_pre=Z; X_pre=X; H_pre=H;
     
        for v=1:V
               %% update Ex
          tempE=X{v}-P{v}*H{v}+Lambda1{v}/rho1;
          Ex{v} = max(abs(tempE) - lambda/rho1, 0).*sign(tempE);  
%           Ex{v}= max(0,tempE-lambda/rho1)+min(0,tempE+lambda/rho1);
%           solve_l1l2(X{v}-P{v}*H{v}+Lambda1{v}/rho1,lambda/rho1); 
      clear tempE
          %% update Pv
            tempD=H{v}*(Lambda1{v}'+rho1*(X{v}'-Ex{v}'));
            [U1,~,V1]=svd(tempD,'econ');
           P{v}=V1*U1';   
             clear tempD
         %% update Hv
        tmpH=(rho2*Eh{v}-Lambda2{v})*(eye(N,N)-Z{v}')+P{v}'*(rho1*X{v}-rho1*Ex{v}+Lambda1{v});
        H{v}=tmpH*inv((gamma+rho1)*eye(N,N)+rho2*(eye(N,N)-Z{v})*(eye(N,N)-Z{v}'));
         clear tmpH
           %% update Xv
         X{v}=P{v}*H{v}+Ex{v}-Lambda1{v}/rho1;
         tmpx{v}=X{v};
         ind_1 = find(Omega(:,v) == 1);
         X{v}(:,ind_1)=T{v}(:,ind_1);

          %% update Zv
         tmpz = H{v}'*(Lambda2{v} + rho2*H{v} - rho2*Eh{v})  +  rho3*y{v}+ lambda3{v};
         Z{v}=max(inv(rho3*eye(N,N)+ rho2*H{v}'*H{v})*tmpz,0);     
         %% update Eh
%         Eh{v}=solve_l1l2(H{v}-H{v}*Z{v}+Lambda2{v}/rho2,lambda/rho2); 
          tempE=H{v}-H{v}*Z{v}+Lambda2{v}/rho2;
          Eh{v} = max(abs(tempE) - lambda/rho2, 0).*sign(tempE);  

        %% update Lambda
        Lambda1{v}=Lambda1{v}+rho1*(X{v}-P{v}*H{v}-Ex{v});
        Lambda2{v}=Lambda2{v}+rho2*(H{v}-H{v}*Z{v}-Eh{v});
        end
         

    %% update M
    Z_tensor = cat(3, Z{:,:});
    L_tensor = cat(3, lambda3{:,:});

    Y=reshape(full_tr(tensor_ring_als(reshape(Z_tensor - 1/rho3*L_tensor,rX),beta*ones(length(rX),1))),[N,N,V]);   

    Lambda3=Lambda3+rho3*(Y-Z_tensor);
    T_Y=reshape(Y,[N,N,V]);
    T_Lambda3=reshape(Lambda3,[N,N,V]);
   
       for v=1:V
         y{v}=T_Y(:,:,v);
         lambda3{v}=T_Lambda3(:,:,v);
       end
     
     rho1=min(rho1*miu,1e12);
     rho2=min(rho2*miu,1e12);
     rho3=min(rho3*miu,1e12);


  
    diff_Z = 0;
    diff_X = 0;
    diff_H = 0;
    res    = 0;
  
    %% check convergence
    for v = 1:V
         ind_1 = find(Omega(:,v) == 1);
        Rec_error = norm( tmpx{v}(:,ind_1)-T{v}(:,ind_1),'fro')/norm(T{v}(:,ind_1),'fro');
        res = max(res,max(abs(Rec_error(:))));
        diff_Z = max(diff_Z,max(abs(Z{v}(:)-Z_pre{v}(:))));
        diff_X = max(diff_X,max(abs(X{v}(:)-X_pre{v}(:))));
        diff_H = max(diff_H,max(abs(H{v}(:)-H_pre{v}(:))));        
    end

%     err = max([leqm1,leqm2,diff_Z,diff_E,diff_B,diff_P]);
%     err = max([leqm1,leqm2]);
%     fprintf('iter = %d,\t difZ = %.2d, \t difX = %.2d,\t difH = %.2d,\t err = %.6f \n'...
%             ,iter, diff_Z, diff_X, diff_H,res);
%     obj(iter) = err;  
    
%     Rec_sample = reshape(E{2}(:,1),50,40);
%     imshow(Rec_sample,[]);title('E')
    if ((res < 1e-6)|| (diff_Z < 1e-6)|| (diff_X < 1e-6) || (diff_H < 1e-6))&& iter>10     
        break
    end
    
    Res(iter)=res;
    Diff_Z(iter)=diff_Z;
    Diff_H(iter)=diff_H;
    Diff_X(iter)=diff_X;
end
    
   Convergence_curve.res=Res;
  Convergence_curve.diff_Z=Diff_Z;
  Convergence_curve.diff_H=Diff_H;
  Convergence_curve.diff_X=Diff_X;   
  
S = 0;
for v=1:V
    S = S + abs(Z{v})+abs(Z{v}');
end


%  C = SpectralClustering(S,cls_num);


% %% ouput
% out.time     = toc;
% out.iter     = iter;
% out.nois     = nois;
% out.L       = L; 
% out.relChgZPath = relChgZPath(1:iter);
% out.relChgXPath = relChgXPath(1:iter);
% return;

end


