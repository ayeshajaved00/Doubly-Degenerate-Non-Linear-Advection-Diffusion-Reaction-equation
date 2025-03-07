function [y] = tridiagQQ( AA, B, C, F )
%This Function Solves Tridiagonal Matrices
%  Solve the  n x n  tridiagonal system for y:
%
%  [ a(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
%  [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
%  [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
%  [                                 b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
%
%  f must be a vector (row or column) of length n
%  a, b, c must be vectors of length n (note that b(1) and c(n) are not used)

% some additional information is at the end of the file

n = length(F);
v = zeros(n,1);   
y = v;
w = AA(1);
y(1) = F(1)/w;
for i=2:n
    v(i-1) = C(i-1)/w;
    w = AA(i) - B(i).*v(i-1);
    y(i) = ( F(i) - B(i).*y(i-1) )/w;
end

for j=n-1:-1:1    
   y(j) = y(j) - v(j).*y(j+1);
end



end