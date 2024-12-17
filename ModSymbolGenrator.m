function modulated_symbols = ModSymbolGenrator(npsk,n)
%modulated_symbols: generate n modulated signals
As=[1 sqrt(2) 0 sqrt(10) 0 sqrt(42)];
M=2^npsk;
A=As(npsk);
data=randi([0,M-1],n,1);
modulated_symbols=qammod(data,M)./A;
end