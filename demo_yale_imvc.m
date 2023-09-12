 clc;
 clear;
%% cho0se methods:
% set to 0 for turning off;
EN_FSR_IMVC    = 1;
EN_s_FSR_IMVC  = 1; 
methodname ={'FSR_IMVC','sFSR_IMVC'};

%% load datasets:
 load('yale.mat')
 %% setting rX
   %Yale 
   rX=[11,15,11,15,3];
   % ALOI rX=[108 100 108 100 4];
   % scene15  rX=[65,69,65,69,3];
   % CCV rX=[67 100 67 100 3];
   % Caletch rX=[114,80,114,80,5];
   % Reuters   rX=[150 125 150 125 5];
   % Coil20  rX=[36,40,36,40,3];
   % BDGP   rX=[50,50,50,50,4]
   
  gt=double(gt);
  numClust = length(unique(gt));
%% set random observed ratio (n_p is the observes sample)
MR=[0.1,0.3,0.5,0.7];
V=length(X);
for v=1:V
    X{v}=X{v}';
end
N=size(X{1},1);

 
%% Note: needs to re-randomly generate missing set by 10 times
for m=1:length(MR)
Omega=zeros(N,1);
enList =[];
% n_p=0.45;
for v = 1:V-1
    ind_folds(:,v)=ones(N,1);
    rng('default');
    rng('shuffle');
    ind = randsample(N,floor(N*MR(m))); 
%     ind=find(rand(N,1)< MR(m));
    ind_folds(ind,v)=0;
    Omega=Omega|ind_folds(:,v);
end
indv=find(Omega);
if length(indv)> floor(N*MR(m))
% ind=find(rand(length(indv),1)< MR(m));
  ind = randsample(length(indv),floor(N*MR(m)));
else
  ind = randsample(length(indv),length(indv));
end
  
ind_folds(:,V)=ones(N,1);
ind_folds(indv(ind),V)=0;

i=0;


%% FSR_IMVC
i=i+1;
if EN_FSR_IMVC
     disp(['performing ',methodname{i}, ' ... ']);
     t0=tic;
    [Result{i}]=FSR_IMVC_clustering(X,ind_folds,gt,rX);  
     Time(i)=toc(t0);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' ACC ' num2str(Result{i}.ACC), ' .'])
    disp([methodname{i}, ' NMI ' num2str(Result{i}.NMI), ' .'])
    disp([methodname{i}, ' Purity ' num2str(Result{i}.Purity), ' .'])
    disp('...')
    enList = [enList,i]; 
end

%% sFSR_IMVC
i=i+1;
if EN_s_FSR_IMVC
     disp(['performing ',methodname{i}, ' ... ']);
     
%% parameter setting
 %% yale
   j=8;k=1e-4;a=0.1;
   %% coil j=4;k=1e-3;a=0.05;
   %% BDGP j=4;k=1e-5;a=0.002; 
   %% Scence15 j=4;k=1e-3;a=0.002
   %% CCV      j=4;k=10;a=0.002
   %% Calteh-all j=10;k=1e-4;a=0.002
   %% Animal   j=6;k=1e-2;a=0.002
   %% AOLI    j=6;k=1e-3;a=0.001
   %% Reuters j=2;k=1e-3;a=0.0001

 paras.lambda=1;
   paras.beta=j;
   paras.gamma=k;
   K = length(X);  %sample number
      paras.M=floor(a*N);
       paras.rX   = [floor(a*N),rX(3),rX(4),rX(5)];
paras.K=numClust+1;% subspce dimension
paras.miu=1.8;
t0=tic;
Result{i}=sFSR_IMVC_clustering(X,ind_folds,gt,paras); 
Time(i)=toc(t0);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' ACC ' num2str(Result{i}.ACC), ' .'])
    disp([methodname{i}, ' NMI ' num2str(Result{i}.NMI), ' .'])
    disp([methodname{i}, ' Purity ' num2str(Result{i}.Purity), ' .'])
    disp('...')
    enList = [enList,i]; 
end
%% Show result
fprintf('\n');
Missing_Ratio=MR(m)
fprintf('================== Result =====================\n');
fprintf(' %8.8s \t   %5.4s \t   %5.4s \t  %5.6s  \t   %5.6s   \n','method','Time','ACC', 'NMI','Purity' );
for i = 1:length(enList)
    fprintf(' %8.8s \t  %5.3f \t %5.3f \t   %5.3f \t    %5.3f    \n',...
    methodname{enList(i)},Time(i),Result{enList(i)}.ACC, Result{enList(i)}.NMI,Result{enList(i)}.Purity);
end

end


function Clu_result=sFSR_IMVC_clustering(X,ind_folds,gt,paras)
%% parameter setting
% lambda=[1e-2,5e-1,1e-1,1,5,10];
    cls_num = length(unique(gt));
for iv = 1:length(X)
%     X{iv}(end-3:end,:)=[];
    X1 = X{iv}';   
    X1 = NormalizeFea(X1,0);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = 0;    % 缺失视角补0
    Y{iv} = X1;         % 一列一个样本
end
clear X X1 ind_0
X = Y;
clear Y  

  %% completion and clustering
    [X_mera,S,Convergence_curve] = FR_IMVC(X,ind_folds,X, paras); 
    
     
      [U,~,~]=svd(S','econ');
      for c=1:10
%      C = SpectralClustering(S,cls_num);% C = kmeans(U,numClust,'EmptyAction','drop'); 
      C=kmeans(U, cls_num);
     result_CLU = ClusteringMeasure(double(gt), C)*100;
     acc(c)=result_CLU(1);
     nmi(c)=result_CLU(2);
     purity(c)=result_CLU(3);
      end
      
  ACC=mean(acc);
  NMI=mean(nmi);
  Purity=mean(purity);
                 
%                 
      for c=1:10
%      C = SpectralClustering(S,cls_num);% C = kmeans(U,numClust,'EmptyAction','drop'); 
     C=kmeans(S', cls_num);
     result_CLU = ClusteringMeasure(double(gt), C)*100;
     acc(c)=result_CLU(1);
     nmi(c)=result_CLU(2);
     purity(c)=result_CLU(3);
      end
      
  ACC1=mean(acc);
  NMI1=mean(nmi);
  Purity1=mean(purity);

 [max_ACC,id] = max(ACC(:));
Clu_result.ACC = max_ACC;
Clu_result.NMI = NMI(id);
Clu_result.Purity = Purity(id);

end


function [Clu_result,parameters]=FSR_IMVC_clustering(X,ind_folds,gt,rX)
%% parameter setting
% lambda=[1e-2,5e-1,1e-1,1,5,10];
    cls_num = length(unique(gt));
for iv = 1:length(X)
%     X{iv}(end-3:end,:)=[];
    X1 = X{iv}';
    
    X1 = NormalizeFea(X1,0);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = 0;    % 缺失视角补0
    Y{iv} = X1;         % 一列一个样本
end
clear X X1 ind_0
X = Y;
clear Y




%% yale:
i=100;j=8; k=1e-5;
%% BDGP: i=1; j=6;k=1e-5;
%% coil20:i=10;j=6;k=1e-5
%% scene15: i=1;j=4;k=1e-5;
 paras.lambda=i;
 paras.beta=j;
 paras.gamma=k;
% rX=[35,42,70,126];

N = size(X{1},2); %sample number
paras.rX   = rX;
paras.K=cls_num+1; % subspce dimension
paras.miu=1.5;
       
% -------------------0。5------------------- clustering 

  %% completion and clustering
     [X_mera,S,Convergence_curve] = FSR_IMVC(X,ind_folds,X, paras);   
     for c=1:10
     C = SpectralClustering(S,cls_num);% C = kmeans(U,numClust,'EmptyAction','drop');  
     result_CLU = ClusteringMeasure(gt, C)*100;
     acc(c)=result_CLU(1);
     nmi(c)=result_CLU(2);
     purity(c)=result_CLU(3);
     end
Clu_result.ACC = mean(acc);
Clu_result.NMI = mean(nmi);
Clu_result.Purity = mean(purity);
end