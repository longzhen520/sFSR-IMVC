
function [t] = tensor_ring_als(c,Rank,mu)


% Tensor Ring toolbox
% written by Qibin Zhao, BSI, RIKEN

if (nargin == 0)
    t.d    = 0;
    t.r    = 0;
    t.n    = 0;
    t.node = 0;                    % empty tensor
    t = class(t, 'tensor_ring');
  
    return;
end
% ip = inputParser;
% ip.addParamValue('Tol', 1e-6, @isscalar);
% ip.addParamValue('Rank', [], @ismatrix);
% ip.addParamValue('MaxIter', 1000, @isscalar);
% ip.addParamValue('weight', 0.01, @isscalar);
% ip.addParamValue('lambda', 0.01, @isscalar);
%  ip.parse(varargin{2:end});
% 
% Tol = ip.Results.Tol;
% Rank = ip.Results.Rank;
% MaxIter = ip.Results.MaxIter;
% lambda=ip.Results.lambda;


 maxit=8;
 Tol=1e-5;
 
%  Rank=[8,8,8,8,10];
        n = size(c);
        n = n(:);
        d = numel(n);
        node=cell(1,d);
        r=Rank(:);
        for i=1:d-1
            node{i}=randn(r(i),n(i),r(i+1));
        end
        node{d}=randn(r(d),n(d),r(1));
        od=[1:d]';
        err=1;
        for it=1:maxit
            err0=err;
            if it>1
                c=shiftdim(c,1);
                od=circshift(od,-1);
            end
            c=reshape(c,n(od(1)),numel(c)/n(od(1)));
            b=node{od(2)};
            for k=3:d
                j=od(k);
                br=node{j};
                br=reshape(br,[r(j),numel(br)/r(j)]);
                b=reshape(b,[numel(b)/r(j),r(j)]);
                b=b*br;
            end
            b=reshape(b,[r(od(2)),prod(n(od(2:end))),r(od(1))]);
            b=permute(b,[1,3,2]);
            b=reshape(b,[r(od(2))*r(od(1)), prod(n(od(2:end)))]);
%              a=c*b'*inv(b*b'+mu*eye(r(od(2))*r(od(1))));
               a=c/b;
              % Un(i,:)=Un(i,:).*(Xn(i,idx)*Uneq(idx,:))./(Un(i,:)*(Uneq(idx,:)'*Uneq(idx,:))); 

%             inv(Uneq(idx,:)'*Uneq(idx,:)+para.mu*eye(R(n)*R(n+1)))
            
            err=norm(c-a*b,'fro')/norm(c(:));
            a=reshape(a,[n(od(1)),r(od(2)),r(od(1))]);
            node{od(1)}=permute(a,[3,1,2]);
            s=norm(node{od(1)}(:));
            node{od(1)}=node{od(1)}./s;
            
%             fprintf('it:%d, err=%f\n',it,err);
%             error(it)=err;
            if abs(err0-err)<=1e-5 && it>=2*d && err<=Tol
                break;
            end
            c=reshape(c,n(od)');
          
            
       
            
        end

          node{od(1)}=node{od(1)}.*s;
        t.node=node;
        t.d=d;
        t.n=n;
        t.r=r;    
        return;
    end