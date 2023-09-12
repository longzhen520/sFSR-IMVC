function [a] = full_tr(tr,flag)
%Converts TR-tensor to multi array
%
%
%
%---------------------------

d=tr.d; n=tr.n;  node=tr.node; r=tr.r;



if 1        
    [~,id]=min(r);
    
    
    a=node{id};
    for i=[id+1:d,1:id-1]
        cr=node{i};
        if i==d
            cr=reshape(cr,[r(i),n(i)*r(1)]);
        else
            cr=reshape(cr,[r(i),n(i)*r(i+1)]);
        end
        a=reshape(a,[numel(a)/r(i),r(i)]);
        a=a*cr;
    end
    
    a=reshape(a,[r(id),prod(n),r(id)]);
    a=permute(a,[2,3,1]);
    a=reshape(a,[prod(n),r(id)*r(id)]);
    temp=eye(r(id),r(id));
    a=a*temp(:);
    
    a=reshape(a,n([id:d,1:id-1])');
    a=shiftdim(a,d-id+1);
    a=a(:);
    
    
end


if 0    
    a=node{1};
    for i=2:d
        cr=node{i};
        if i==d
            cr=reshape(cr,[r(i),n(i)*r(1)]);
        else
            cr=reshape(cr,[r(i),n(i)*r(i+1)]);
        end
        a=reshape(a,[numel(a)/r(i),r(i)]);
        a=a*cr;
    end
    
    a=reshape(a,[r(1),prod(n),r(1)]);
    a=permute(a,[2,3,1]);
    a=reshape(a,[prod(n),r(1)*r(1)]);
    temp=eye(r(1),r(1));
    a=a*temp(:);
    
end

if 0
    a=node{2};
    for i=3:d
        cr=node{i};
        if i==d
            cr=reshape(cr,[r(i),n(i)*r(1)]);
        else
            cr=reshape(cr,[r(i),n(i)*r(i+1)]);
        end
        a=reshape(a,[numel(a)/r(i),r(i)]);
        a=a*cr;
    end
    
    
    a=reshape(a,[r(2),prod(n(2:end)),r(1)]);
    a=permute(a,[1,3,2]);
    a=reshape(a,[r(2)*r(1), prod(n(2:end))]);
    
    b=node{1};
    b=permute(b,[2,3,1]);
    b=reshape(b,[n(1),r(2)*r(1)]);
    
    a=b*a;
    a=a(:);
    
    
    
end






if (nargin>1)&& (flag==1)
    a = reshape(a, n');
end;



return;
