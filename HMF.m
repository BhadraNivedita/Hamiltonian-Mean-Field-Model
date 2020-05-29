% Simulation of the HMF model

function dy=HMF(t,y,N,K)

% Setting the parameters
dy=zeros(2*N,1);

% equations for the many body Kapitza pendula
for i=1:N
    dy(i) = y(N+i);
end

for i=1:N
    for j=1:N
        dy(N+i) = dy(N+i) + (K/N)*sin(y(j)-y(i));
    end
end

end