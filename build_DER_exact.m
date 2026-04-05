function [Aj,Bj,Cj] = build_DER_exact()

% parámetros del paper (~15 ms)
tau = [0.012 0.013 0.015 0.014 0.015 0.014];

kj1 = 1; 
kj2 = 1;

Aj = -diag(1./tau);

Bj = zeros(6,2);

% según ecuación del paper
Bj(1,1) = 1/tau(1);
Bj(2,2) = 1/tau(2);

Bj(3,1) = kj1/(tau(3)*(tau(3)-tau(4)));
Bj(4,1) = kj1/(tau(4)*(tau(4)-tau(3)));

Bj(5,2) = -kj2/(tau(5)*(tau(5)-tau(6)));
Bj(6,2) = -kj2/(tau(6)*(tau(6)-tau(5)));

Cj = [1 0 1 1 0 0;
      0 1 0 0 1 1];

end