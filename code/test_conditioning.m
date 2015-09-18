% test_conditioning.m

one = 1:n;
two = n+1 : n+m;

d = -diag(A(one,one));
di = 1.0./d;
di(abs(d) < eps ) = 1.0;

AD = A(two,one) * diag(di);
S = AD * AD';
%S = A(two,one) * (-A(one,one)\A(one,two));
PS = zeros(m,m);

for i = 1:m
   vec_in = zeros(n+m,1);
   vec_in(n+1:n+m) = S(:,i);
   vec_out = indef_pre(vec_in);
   PS(i,:) = vec_out(n+1:n+m);
end

min(eig(full(PS)))