function [A,B,Ctilde] = build_global_model(Aj,Bj,Cj,l)

A = kron(eye(l), Aj);
B = kron(eye(l), Bj);
Ctilde = kron(eye(l), Cj);

end