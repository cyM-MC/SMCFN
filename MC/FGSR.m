%baseline: FGSR 
%regula_B = L2 norm
function [AB, Time] = FGSR(Smiss,M, d, alpha, u)
%Smiss: input incomplete matrix, M: mask [0,1], d: rank, 
%hyperpara alpha: = u: 0.001, maxIter:iteration 
[m,n]=size(Smiss);
X=Smiss; 
tolFGSR = [];

%disp(['alpha=' num2str(alpha)]) 
if min(m,n)<5000
    [U,S,V]=svd(Smiss*m*n/sum(M(:)),'econ');
%     switch regul_B
%         case 'L2'
            A=alpha^(1/3)*U(:,1:d)*(S(1:d,1:d).^(2/3));
            B=alpha^(-1/3)*(S(1:d,1:d).^(1/3))*V(:,1:d)';
%         case 'L21'
%             A=U(:,1:d)*(S(1:d,1:d).^(1/2));
%             B=(S(1:d,1:d).^(1/2))*V(:,1:d)';
%     end
else
    A=randn(m,d);
    B=randn(d,m);
end    
% if isfield(options,'tol')==0
     tol=1e-4;
% else
%     tol=options.tol;
% end
%
iter=0;
Q=zeros(m,n);
%

% Initialize time variables 
Time =0;
startTime = cputime;

%while iter<maxIter
while iter<10000 % Pg7 default setting
    iter=iter+1;
    % A_new
    tau=1*u*normest(B*B');
    temp=A-u*((X-A*B+Q/u)*(-B'))/tau;
    A_new=solve_l1l2(temp,1/tau);
    % B_new
%     switch regul_B
%         case 'L2'
            B_new=inv(u*A_new'*A_new+alpha*eye(d))*(u*A_new'*(X+Q/u));
%         case 'L21'
%             tau=3*u*normest(A_new)^2;
%             temp=B-u*(-A_new')*((X-A_new*B+Q/u))/tau;
%             B_new=solve_l1l2(temp',alpha/tau)';
%         case 'L1'
%             tau=1.5*u*normest(A_new)^2;
%             temp=B-u*(-A_new')*((X-A_new*B+Q/u))/tau;
%             B_new=max(0,temp-alpha/tau)+min(0,temp+alpha/tau);
%     end
    % X_new
    AB=A_new*B_new;
    X_new=AB-Q/u;
    X_new=X_new.*(1-M)+Smiss.*M;
    % 
    Q=Q+u*(X_new-AB);    %
    et=[norm(A_new-A,'fro')/norm(A,'fro') norm(B_new-B,'fro')/norm(B,'fro')  norm(X_new-X,'fro')/norm(X,'fro')];
    %Convergence test
    tolFGSRt = norm(X_new - X, 'fro')/norm(X_new,'fro');
    tolFGSR = [tolFGSR,tolFGSRt];
    stopC=max(et);
    %
    isstopC=stopC<tol||et(3)<tol/10;
    if mod(iter,50)==0||isstopC||iter<=5
        %disp(['iteration=' num2str(iter) '/' num2str(maxIter) '  stopC=' num2str(stopC) ...
        %disp(['iteration=' num2str(iter) '/' num2str(1000) '  stopC=' num2str(stopC) ...
        %    '  e_A=' num2str(et(1)) '  e_B=' num2str(et(2)) '  e_X=' num2str(et(3)) ' d=' num2str(d)])
    end
    if isstopC
        %disp('converged')
        break;
    end
    B=B_new;
    A=A_new;
    X=X_new; 
    

    % extract nnz
    id1=find(sum(A.^2)>1e-5);
    id2=find(sum((B').^2)>1e-5);
    id=intersect(id1,id2);
    A=A(:,id);
    B=B(id,:);
    d=length(id);
    %     
end
save('tolFGSR.mat','tolFGSR');
 
Time = cputime - startTime;
end
%%
function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end
 