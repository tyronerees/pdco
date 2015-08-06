function y = minres_preconditioner(x,pre);

% set up the Preconditioner
    
    y(pre.n+pre.m,1) = 0.0;
    one = 1:pre.n;
    two = pre.n+1:pre.n+pre.m;
    
    switch pre.number
      case (1)
        % Augmented Lagrangian
        y(two) = pre.W\x(two);
        y(one) = pre.Aaug\x(one);
      case (2)
        % MI30
        y = pre.pc.apply(x);
      case (3) 
        %        keyboard
        y = pre.U\(pre.L\(pre.P*x));
      case (4) 
        y = pre.P\x;
      case (5)
        y(one) = pre.A \x(one);
        y(two) = pre.SD\x(two);
    end