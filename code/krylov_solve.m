function [x,its,normr,data] = krylov_solve(A,b,n,m,tol,maxits,krylov_method,premeth,conv_data)
% Pick a krylov method....

%!! add this option later....
const = sqrt(2)*conv_data.z_max + conv_data.max_sigma* ...
        conv_data.x_max;
newtol = tol*((conv_data.mu)^0.5)/const;
%!!


switch krylov_method
  case 1 % Saunders' MINRES implementation
    damp = 0;
    check = 0;
    show = 0;
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
    [x, istop_minres, its, normr, normAr, ...
     normA, condA, normx] = saunders_minres( ...
         -A, -b, indef_pre, damp, show, check,...
         maxits, tol);
    if istop_minres==-1 || istop_minres==5 || istop_minres==6   % conlim or itnlim
        fprintf('\n    MINRES stopped early:  istop = %3d \n', istop_minres)
    end
  case 2 % My MINRES 
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
    [x, its, resvec] = minres_t(-A,-b,indef_pre,zeros(size(b)),...
                                               tol, maxits, 0);
    data.resvec = resvec;
    normr = resvec(its);
  case 3 % Matlab's minres 
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
    here = pwd; cd  ~/code/ifiss3.3/matlab704/ % avoid calling pdco
                                               % function with same name...
    [x,flag,normr,its,data.resvec] = minres(-A,-b,tol,maxits,indef_pre);
    cd(here)
    if flag > 0
        fprintf('\n   MINRES failed, flag = %3d \n ',flag);
        keyboard
    end
  case 4 % my GMRES on the indefinite system
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
        
    [x,normr,its,resvec]=...
        mpgmres(-A,-b,{indef_pre},'trunc',tol,maxits);
    data.resvec = resvec;
    data.denom = const;
    data.tol = newtol;
  case 5 % matlab's GMRES on the indefinite system
    predata = generate_pre(A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
    fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    fprintf('\n Warning: needs maxits issue needs fixing!')
    fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    [x,flag,normr,itvec] = gmres(A,b,[],tol,[], ...
                               indef_pre);
    its = itvec(2);
    if flag ~= 0
        fprintf('\n\n**  Error, flag = %d **\n', flag)
        keyboard    
    end
  case 6 % (my) Conjugate gradients
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
          
    [x, its, resvec] = ...
        preconjgrad(-A,-b,maxits,zeros(size(b)),newtol,indef_pre);
        
    normr = norm(resvec(its));
    data.resvec = resvec;
    data.tol = newtol;
  case 7 % tfqmr
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
          
    [x,flag,normr,its] = ...
        tfqmr(-A,-b,tol,maxits,indef_pre);
    if flag ~= 0
        fprintf('\n\n**  Error, flag = %d **\n', flag)
        keyboard    
    end
  case 8 % symmlq
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
          
    [x,flag,normr,its] = ...
        symmlq(-A,-b,tol,maxits,indef_pre);
    if flag ~= 0
        fprintf('\n\n**  Error, flag = %d **\n', flag)
        keyboard    
    end
  case 9 
    %BiCGStab
    predata = generate_pre(-A,n,m,premeth);
    indef_pre = @(x) minres_preconditioner(x,predata);
          
    [x,flag,normr,its] = ...
        bicgstab(-A,-b,tol,maxits,indef_pre);
    if flag ~= 0
        fprintf('\n\n**  Error, flag = %d **\n', flag)
        keyboard    
    end
  case 10 
    % a full-blown constraint preconditioner...
    predata = generate_pre(-A,n,m,4);
    indef_pre = @(x) minres_preconditioner(x,predata);
    xhat = indef_pre(-b);
    bhat = zeros(n+m,1);
    bhat(1:n) =  -b(1:n) + A(1:n,1:n)*xhat(1:n);
    [x, its, resvec] = ...
        preconjgrad(-A,bhat,maxits,zeros(size(b)),1e-10,indef_pre);
    data.resvec = resvec;
    x(1:n) = x(1:n) + xhat(1:n);
    %    keyboard
    normr = norm(resvec(its));
  case 11
    % CG on the normal equations
    S = -A(n+1:n+m,1:n)*(A(1:n,1:n)\A(1:n,n+1:n+m));
    rhs = -(A(n+1:n+m,1:n)*(A(1:n,1:n)\b(1:n)) - b(n+1:n+m));
    x = zeros(n+m,1);
    
    
    try
        pc= hsl_mi28_precond(S);
        pre = @(z) pc.apply(z);
    catch
        keyboard
    end
    exact = 1;
    if exact
        x_ex = S\rhs;
        [sol, its, resvec] = ...
            preconjgrad_exact(S,rhs,maxits,zeros(m,1),newtol,pre,x_ex);
    else
        [sol, its, resvec] = ...
            preconjgrad(S,rhs,maxits,zeros(m,1),newtol,pre);
    end

    clear pc;
    
    x(n+1:n+m) = sol;  %l;
    
    %    if ans > 1e-10
    %        keyboard 
    %    end
    %    keyboard
    %    x(n+1:n+m) = S\rhs;
    
    
    x(1:n) = A(1:n,1:n)\(b(1:n)- A(1:n,n+1:n+m)*x(n+1:n+m));
    
    its = 0;
    normr = 0;
    data.bleugh = 0;
    
end

