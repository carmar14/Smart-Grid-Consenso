function xhat = local_observer(xhat,u,y,Ad,Bd,C,K)

xhat = Ad*xhat + Bd*u + K*(y - C*xhat);

end