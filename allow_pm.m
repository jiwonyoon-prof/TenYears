function B = allow_pm(A, q)
B = A;
Jind = find(A>q/2);
B = A;
B(Jind) = B(Jind)-q;
