u_0 = initialize_rpl_n_prior(Y,X,W,mu_rpl_tv,par_rpl);

%% OLD 

f = band(:);
p = 3; q = 2; r = 2;

band = linspace(10,20,r);

X = kron( speye(p*q), ones(size(f)) );
Z = kron( speye(p*q), log(f) );
W = kron( speye(p*q), -4*f );

figure, 
subplot(131)
imagesc(X), title('X');

subplot(132),
imagesc(Z), title('Z');

subplot(133),
imagesc(W), title('W');
%% Attempt EMZ (coila style


X = kron( ones(size(f)) , speye(p*q) );
Z = kron( log(f) , speye(p*q) );
W = kron( f , speye(p*q) );


figure, 
subplot(131)
imagesc(X), title('X');

subplot(132),
imagesc(Z), title('Z');

subplot(133),
imagesc(W), title('W');
