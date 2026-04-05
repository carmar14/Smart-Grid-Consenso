function [Ad,Bd] = discretize_model(A,B,Ts)

Ad = eye(size(A)) + Ts*A;
Bd = Ts*B;

end