function pre = generate_pre(K,n,m,number)

% $$$ if nargin<5
% $$$     gamma = norm(K(one,one),'fro'); 
% $$$     % else we have gamma = mu;
% $$$ end

pre.n = n;
pre.m = m;
pre.number = number;
one = 1:pre.n;
two = pre.n+1:pre.n+pre.m;

switch pre.number
  case 1
    % Aug lag
    %    pre.gamma = max(max(abs(K(one,one))));
    %    keyboard
    % will calculate ||H + rhoI + X^{-1}Z||_F (Morini, Simoncini, Tani)
    %    normA = norm(K(one,one),'fro');  
   
    pre.gamma = norm(K(one,one),'fro'); 
    %    pre.gamma = gamma;%(normA);
    pre.W = (pre.gamma*speye(m) + K(two,two));
    pre.Aaug = K(one,one) + K(one,two)*((1/pre.gamma)*K(two,one));

  case 2
    % MI30
    control.iscale = 3; % mc64
    pre.pc = hsl_mi30_precond(K,m,control); 
    
  case 3 
    % ilu
    %    setup.type = 'nofill';
    setup.type = 'crout';
    [pre.L,pre.U,pre.P] = ilu(K,setup);
  case 4
    % constraint preconditioner
    pre.P = K;
    pre.P(one,one) = diag(diag(K(one,one)));
  case 5
    % a simple diagonal preconditioner...
    d = diag(K(one,one));
    d(d==0) = 1e-8;
    D = diag(d);

    pre.A = K(one,one) + D - diag(diag(K(one,one)));
    
    pre.SD = K(two,one)*(D\K(one,two));
    
end
