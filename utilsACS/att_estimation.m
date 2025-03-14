function [att_RED, att_RSLD, SNR, x, z, ind_x, ind_z] = att_estimation(att_ref_dB, BW, fs, c0, ROI, transducer_footprint, blocksize, ref, sam, mu_rsld)

%% ------------------------------------------------------------------------
% DESCRIPTION
% Provides:
% * att_RED     attenuation coefficient map calculate by RED
% * att_RSLD    attenuation coefficient map calculate by RSLD
% * SNR         SNR map

%% Parameters
% Attenuation reference:    att_ref_dB   % [dB/cm/MHz]      [alpha_2 alpha_1 alpha_0]
% Bandwidth:                BW           % [MHz]            [freq_L freq_H]
% Sampling Frequency:       fs           % [Hz]             integer
% Sound velocity:           c0           % [m/s]            integer 
% Region of interest:       ROI          % [mm]             [x1, x2, z1, z2]
% Bmode Width:              transducer_footprint % [mm]     integer
% Box size:                 blocksize    % [lambda]
% Reference Phantom         ref          3D matrix
% Sample Tissue             rf           2D matrix

%%% Reference phantom
% RFdata_ref_liver = load('LiverPhantom_RFdata2D_2mm.mat');
% ref = RFdata_ref_liver.RFdata;

%%% Sample Tissue
% RFdata_liver = load("mpUS_2021_2\2week\591\RFdata2D_1.mat");
% sam = RFdata_liver.RF(:,:,frame);
env = abs(hilbert(sam));
%% Settings
winsize     = 0.5;
window_type = 5;
saran_layer = 0;
overlap_pc = 0.5; %0.8;

%% Hyperparameters
% Coila
% mu_rsld = 10^0.6; %10^0.6;%1e1;
filter_win = [9 9];

% Anthony
mu_red = 10^7.7; %1e6; %1e8; %1e4;
median_filter_win = 5; %5
% sigma  = 1e3; %1e3;
% alpha  = 1.5;
% beta   = 2.4;
max_iter = 1000; %2000; %50

%% Initialization
freq_L = BW(1)*1e6;   % [Hz]
freq_H = BW(2)*1e6;   % [Hz]
wl = c0/mean([freq_L freq_H]);   % Wavelength (m)

% x and z at the starting time
x = 0:size(ref,2)-1;
x = x*(transducer_footprint*(10^(-3))/size(ref,2));
z = (0:size(sam,1)-1)*(c0/2)/fs;
last = min([length(z), length((0:size(ref,1)-1)*(c0/2)/fs), length(sam), length(ref)]);
z = z(:,1:last);

% Original axes
% x_ori = x;
% z_ori = z;

% Lateral (x) and axial (z) deltas                
dx=(x(end)-x(1))/(length(x)-1);
dz=(z(end)-z(1))/(length(z)-1);

% x and z in [cm]
x = x*1e2;  
z = z*1e2;

%%% Region of interest selection
x_inf=ROI(1)*1e-1; % [cm]
x_sup=ROI(2)*1e-1; % [cm]
z_inf=ROI(3)*1e-1; % [cm]
z_sup=ROI(4)*1e-1; % [cm]

ind_x = x_inf <= x & x <= x_sup;
% x = x(ind_x);
ind_z = z_inf <= z & z <= z_sup;
% z = z(ind_z);

sam = sam(ind_z,ind_x,:);
ref = ref(ind_z,ind_x,:);
env = env(ind_z,ind_x,:);

% Number of repetitions of datablocks within a single datablock
rpt = 1/(1-overlap_pc);
rpt = round(rpt);
% overlap = 1 - (1/rpt);

% Part without lateral overlap
wx = round((blocksize*wl)/(dx*rpt));

% Number of lines laterally = repetitions * block without repetition
nx = rpt*wx;
% new_blocksize = round(nx*dx/(wl));

% RF data columns
L2 = size(sam,2);

% Blocksize colums
ncol = floor((L2-(rpt-1)*wx)/wx);
sam = sam(:,1:wx*(ncol+rpt-1),:);
env = env(:,1:wx*(ncol+rpt-1),:);

% Actual rf data columns
L2 = size(sam,2);   
% x  = x(1:L2);

xi = 1;
xf = L2;
x0 = (xi:wx:xf+1-nx);
n  = length(x0);
wz = floor(nx*dx/(dz*rpt));
nz = rpt*wz;

% winsize: Percentage of window (max 0.5)
% nw: Samples of each window axially
nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;  
L = (nz - nw)*dz*100;   % [cm]

NFFT = 2^(nextpow2(nw)+2);
band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
% p = length(rang);

L1   = size(sam,1);   % RF data: rows
nrow = floor((L1-(rpt-1)*wz)/wz);        % Blocksize rows
sam = sam(1:wz*(nrow+rpt-1),:,:);
env = env(1:wz*(nrow+rpt-1),:,:);

L1   = size(sam,1);   % RF data: rows
% z    = z(1:L1);
zi = 1;
zf = L1;
z0 = (zi:wz:zf+1-nz);
m  = length(z0);
z0p = z0 + (nw-1)/2;
z0d = z0 + (nz-1) - (nw-1)/2;

ref = ref(1:L1,1:L2,:);

% disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
% num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
% disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ', num2str(overlap*1e2,'%4.0f'),'%']);
% disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
% disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
% disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);

switch window_type
    case 5
        windowing = tukeywin(nw,0.25);  % Tukey Window. Parameter 0.25
    case 6
        windowing = hamming(nw);        % Hamming
    case 7
        windowing = rectwin(nw);        % Boxcar
end

% Windowing neccesary before Fourier transform
windowing = windowing*ones(1,nx);

%% Spectra Processing
Sp = zeros(m,n,length(rang));
Sd = zeros(m,n,length(rang));
SNR = zeros(m,n);

for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);
        
        sub_block_p = sam(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1,:);
        sub_block_d = sam(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1,:);
        env_block =   env(zp-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1,:);
        
        % SNR
        SNR(ii,jj) = mean(env_block(:))./std(env_block(:));
        
        % Sample
        [tempSp,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

att_ref = (att_ref_dB(1)*(f.^2) + att_ref_dB(2)*f + att_ref_dB(3))/8.686;

%% Diffraction compensation
% Attenuation reference
att_ref_map = zeros(m,n,length(rang)); %att_ref_map = zeros(n,m,length(rang));
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Reference phantom Spectrum
Sp_ref = zeros(m,n,length(rang));
Sd_ref = zeros(m,n,length(rang));
for jj=1:n
    for ii=1:m                                
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);
        
        % Reference
        sub_block_p = ref(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1,:);
        sub_block_d = ref(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1,:);

        [tempSp,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        
        Sp_ref(ii,jj,:) = (tempSp(rang));
        Sd_ref(ii,jj,:) = (tempSd(rang));                           
    end
end

diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;

%% Au = b
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

%% Coila
% Without Regularization
% rng('default')
[u,~] = cgs2(A'*A,A'*b(:),1e-6,20);

% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
att_RSLD = u(1:end/2); %CS = u(end/2+1:end);
att_RSLD = 8.686*att_RSLD;   % [dB.cm^{-1}.MHz^{-1}]
att_RSLD = reshape(att_RSLD,m,n);

% disp("Conjugate Gradients Squared")
% disp("beta median: "+median(BS,'all'))
             
% att_RSLD = movmedian(att_RSLD,filter_win);

% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
att_RSLD = IRLS_TV(att_RSLD(:),speye(m*n),mu_rsld,m,n,1e-3,ones(m,n),ones(m*n,1));
att_RSLD = reshape(att_RSLD,m,n);

% disp("IRLS TV - mu: "+mu_tvd)
% disp("beta median: "+median(BS,'all'))

%% Anthony
rng('default')
% [~,~,u_med]  =
% admm_red4(A'*A,A'*b(:),mu_red,sigma,0.01,size(A'*b(:),1),max_iter,4,5,m,n,alpha,beta); v0
% [~,~,u_med]  = admm_red(A'*A,A'*b(:),mu_red,0.001,size(A'*b(:),1),max_iter,10,1,median_filter_win,m,n,mu_red/1.5); v1
[~,~,u_med] = admm_red(A'*A,A'*b(:),mu_red,0.001,size(A'*b(:),1),max_iter,4,1,median_filter_win,m,n,mu_red/1.5);

att_RED = u_med(1:end/2); %CS = u(end/2+1:end);
att_RED = 8.686*att_RED;   % [dB.cm^{-1}.MHz^{-1}]
att_RED = reshape(att_RED,m,n);
% disp("ADMM RED - mu: "+mu_aux)
% disp("beta median: "+median(BS,'all'))

end
%% Local Functions
% A*x = b
function [x,ite_cgs] = cgs2(A,b,tol,maxit,varargin)

    if length(varargin) == 1
        x = varargin{1};
    else
        x= zeros(size(A,2),1);
    end
    r = A*x - b;
    p = -r;
    ite_cgs = 0;
    
    while norm(r,2) > tol && ite_cgs < maxit
        
        alpha = (r'*r)/(p'*(A*p));
        x = x + alpha*p;
        rn = r + alpha*(A*p);
        beta = (rn'*rn)/(r'*r);
        p = -rn + beta*p;
        r = rn;
        ite_cgs = ite_cgs + 1;
    end

end

% TV Andres Leonel Coila
function u = IRLS_TV(b,A,mu,M,N,tol,~,minimask)
% function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)

    [u,~] = cgs2(A'*A,A'*b,1e-6,20);
    
    % Create two variables instead of having a changing size G -> that
    % would be helpful to see the optimization error, but as it is not an
    % output of interest now, we would change that
    
%     G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
    G_0 = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
    
    D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
    D(:,end) = [];
    D(M,M) = 0;
    Dx = kron(speye(N),D);
    
    D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
    D(:,end) = [];
    D(N,N) = 0;
    Dy = kron(D,speye(M));
    
    D = [Dx' Dy']';
    
%     ite_irls = 0;
    error = 1;
    
    while error > tol
        
        X = reshape(u,M,N);
%         ite_irls = ite_irls + 1;
%         Dh = diff(X,[],1);
%         Dh = [Dh;zeros(1,N)];
        Dh = [diff(X,[],1); zeros(1,N)];
        
%         Dv = diff(X,[],2);
%         Dv = [Dv zeros(M,1)];
        Dv = [diff(X,[],2) zeros(M,1)];

        
        %Dx*X(:) - Dh(:);
        %Dy*X(:) - Dv(:);
        
        P = Dh.^2 + Dv.^2;
        eps = 0.01;
        P = 2*sqrt(P.^2 + eps^2);
        P = P.^(-0.5);
        P = P(:).*minimask;
        omega = speye(M*N);
        omega = spdiags(P,0,omega);
        W = kron(speye(2),omega);
        
        AtA = A'*A;
        %mu=5000;
        %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
        [u] = cgs2( AtA + mu*D'*W*D, A'*b, 1e-6 , 20, u );
        
%         G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
        G_1 = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
        
%         error = abs(G(ite_irls+1) - G(ite_irls));
        error = abs(G_1 - G_0);
        G_0 = G_1;
    end
end

% TV Andres Leonel Coila
function [TV] = TVcalc(B,M,N,mask)

    mask(isnan(mask)) = 0;
    mask = mask(:);
    
    X = reshape(B,M,N);
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];
    
    P = Dh.^2 + Dv.^2;
    P = sqrt(P);
    TV = norm(P(:).*mask,1);
end

function [res_admm,n_out , out ] = admm_red(A,y,lambda,error,N,max_iter,inner_iters,inner_iters2,median,m,ni,beta)
   
    x_est= ones([N 1]);
    v_est=ones([N 1]);
    u_est=ones([N 1]);
   
%     x_est= 0.5*abs(x_est)+0.75;
%     v_est= 0.5*abs(v_est)+0.75;
%     u_est= 0.5*abs(u_est)+0.75;
   
    n_out = [x_est v_est u_est];
    alpha=1.5;
   
   
    for k = 1:1:max_iter

        % Part1 of the ADMM, approximates the solution of:
        % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2

            for j = 1:1:inner_iters
%                e = A'*(A*x_est-y) + beta*(x_est-(v_est-u_est));
%                 r = A'*A*e + beta*e;
%                 mu_opt = (e'*e)/(e'*r);
%                 x_est = x_est + mu_opt*e;
%                 %x_est = max( min(x_est, 10), 0);


                 b = (A*y) + beta*(v_est - u_est);
                A_x_est = (A'*A*x_est) + beta*x_est;
                res = b - A_x_est;
                %res = (1 / (input_sigma^2))*(A'*(A*x_est-y))+beta*(x_est -(v_est-u_est));
                a_res = (A'*A*res) + beta*res;
                mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
                x_est = x_est + mu_opt*res;
                x_est = max( min(x_est, 0.5), 0);

            end


        % relaxation
         x_hat = alpha*x_est + (1-alpha)*v_est;

        %x_hat = x_est;
        % Part2 of the ADMM, approximates the solution of
        % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
        % using gradient descent
        for j = 1:1:inner_iters2
           
%              v_est = reshape(v_est,30,100);
%           f_v_est = imnlmfilt(v_est,'DegreeOfSmoothing',10,'SearchWindowSize',15);
%           f_v_est = f_v_est(:);
           
             v_est = reshape(x_est,m,2*ni);
            v_est1=v_est((1:m),(1:ni));
            v_est2=v_est((1:m),(ni+1:2*ni));
           f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
            f_v_est2 = medfilt2(v_est2, [median median],'symmetric');
%            


%             v_est1aux = v_est1;
%             v_est2aux = v_est2;
%             filtroadaptativo(v_est1aux,7);
%             filtroadaptativo(v_est2aux,7);
%            
%             f_v_est1=v_est1aux;
%              f_v_est2=v_est2aux;
             
             
%             f_v_est1 = BM3D(v_est1, median);
%             f_v_est2 = BM3D(v_est2, median);


%             f_v_est1 = imnlmfilt(v_est1,'DegreeOfSmoothing',median);
%             f_v_est2 = imnlmfilt(v_est2,'DegreeOfSmoothing',median);


            f_v_est = [f_v_est1 f_v_est2];
            f_v_est = f_v_est(:);    


           
            v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
        end

        % Part3 of the ADMM, update the dual variable
        u_est = u_est + x_hat - v_est;

         r = norm(A*x_hat-y)^2/norm(y)^2;
       
         if r < error
%              out = x_hat;
             %n_out=k;
%              res_admm = r;
             break
         end
       
       
    end
    out = x_hat;
    % n_out=k;
    res_admm = r;
end

% % RED Anthony v2
% function [res_admm, n_out ,out] = admm_red(A,y,lambda,error,N,max_iter,inner_iters,inner_iters2,median,m,ni,beta)
%    
%     x_est= ones([N 1]);
%     v_est=ones([N 1]);
%     u_est=ones([N 1]);
%     
% %     x_est= 0.5*abs(x_est)+0.75;
% %     v_est= 0.5*abs(v_est)+0.75;
% %     u_est= 0.5*abs(u_est)+0.75;
%     
%     n_out = [x_est v_est u_est];
%     alpha=1.5;
%     
%     
%     for k = 1:1:max_iter
% 
%         % Part1 of the ADMM, approximates the solution of:
%         % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2
% 
%             for j = 1:1:inner_iters
% %                e = A'*(A*x_est-y) + beta*(x_est-(v_est-u_est));
% %                 r = A'*A*e + beta*e;
% %                 mu_opt = (e'*e)/(e'*r);
% %                 x_est = x_est + mu_opt*e;
% %                 %x_est = max( min(x_est, 10), 0);
% 
% 
%                  b = (A*y) + beta*(v_est - u_est);
%                 A_x_est = (A'*A*x_est) + beta*x_est;
%                 res = b - A_x_est;
%                 %res = (1 / (input_sigma^2))*(A'*(A*x_est-y))+beta*(x_est -(v_est-u_est));
%                 a_res = (A'*A*res) + beta*res;
%                 mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
%                 x_est = x_est + mu_opt*res;
%                 %x_est = max( min(x_est, 5), 0);
% 
%             end
% 
% 
%         % relaxation
%          x_hat = alpha*x_est + (1-alpha)*v_est;
% 
%         %x_hat = x_est;
%         % Part2 of the ADMM, approximates the solution of
%         % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
%         % using gradient descent
%         for j = 1:1:inner_iters2
%             
% %              v_est = reshape(v_est,30,100);
% %           f_v_est = imnlmfilt(v_est,'DegreeOfSmoothing',10,'SearchWindowSize',15);
% %           f_v_est = f_v_est(:);
%             
%             v_est = reshape(x_est,m,2*ni);
%             v_est1=v_est((1:m),(1:ni));
%             v_est2=v_est((1:m),(ni+1:2*ni));
%             
% %             f_v_est1 = medfilt2(v_est1, [median median],'symmetric'); % original
% %             f_v_est2 = medfilt2(v_est2, [median median],'symmetric'); % original
% 
% % Kernel variable
% % Filtro de mediana adaptivo -> u otro que tenga adaptabilidad que pueda
% % evitar el efecto de blur
% % Difusión anisotrópica
% 
%             addpath(genpath('bm3d\'))
%             f_v_est1 = BM3D(v_est1, median);
%             f_v_est2 = BM3D(v_est2, median);
% 
% 
% %             f_v_est1 = imnlmfilt(v_est1,'DegreeOfSmoothing',median);
% %             f_v_est2 = imnlmfilt(v_est2,'DegreeOfSmoothing',median);
% 
%             f_v_est = [f_v_est1 f_v_est2];
%             f_v_est = f_v_est(:);    
% 
% 
%             
%             v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
%         end
% 
%         % Part3 of the ADMM, update the dual variable
%         u_est = u_est + x_hat - v_est;
% 
%          r = norm(A*x_hat-y)^2/norm(y)^2;
%        
%          if r < error
% %              out = x_hat;
%              %n_out=k;
% %              res_admm = r;
%              break
%          end
%     end
%     out = x_hat;
%     % n_out=k;
%     res_admm = r;
% end

% % RED Anthony v1
% function [res_admm, n_out, out ] = admm_red4(A,y,lambda,input_sigma,error,N,max_iter,inner_iters,median,m,ni,alpha,beta)
% 
%     x_est= rand([N 1]);
%     v_est=rand([N 1]);
%     u_est=rand([N 1]);
%     
%     x_est= abs(x_est);
%     v_est= abs(v_est);
%     u_est= abs(u_est);
%     
%     n_out = [x_est v_est u_est];
% %     effective_sigma=0.5;
% 
%     for k = 1:1:max_iter
% 
%         % Part1 of the ADMM, approximates the solution of:
%         % x = argmin_z 1/(2sigma^2)||Hz-y||_2^2 + 0.5*beta||z - v + u||_2^2
% 
%             for j = 1:1:inner_iters
%                 b = ((A*y)/(input_sigma^2)) + beta*(v_est - u_est);
%                 A_x_est = ((A'*A*x_est)/(input_sigma^2)) + beta*x_est;
%                 res = b - A_x_est;
%                 %res = (1 / (input_sigma^2))*(A'*(A*x_est-y))+beta*(x_est -(v_est-u_est));
%                 a_res = (A'*A*res)/(input_sigma^2) + beta*res;
%                 mu_opt = mean(res(:).*res(:))/mean(res(:).*a_res(:));
%                 x_est = x_est + mu_opt*res;
%                 %x_est = max( min(x_est, 10), 0);
%             end
% 
% 
%         % relaxation
%         x_hat = alpha*x_est + (1-alpha)*v_est;
% 
% 
%         % Part2 of the ADMM, approximates the solution of
%         % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
%         % using gradient descent
%         for j = 1:1:inner_iters
%             
% %              v_est = reshape(v_est,30,100);
% %           f_v_est = imnlmfilt(v_est,'DegreeOfSmoothing',10,'SearchWindowSize',15);
% %           f_v_est = f_v_est(:);
%             
%             v_est = reshape(x_est,m,2*ni);
%             v_est1=v_est((1:m),(1:ni));
%             v_est2=v_est((1:m),(ni+1:2*ni));
%             f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
%             f_v_est2 = medfilt2(v_est2, [median median],'symmetric');
% %             
% %             f_v_est1 = BM3D(v_est1, median);
% %             f_v_est2 = BM3D(v_est2, median);
% 
% 
%             f_v_est = [f_v_est1 f_v_est2];
%             f_v_est = f_v_est(:);
%             
%             v_est = (beta*(x_hat + u_est) + lambda*f_v_est)/(lambda + beta);
%         end
% 
%         % Part3 of the ADMM, update the dual variable
%         u_est = u_est + x_hat - v_est;
% 
% %          r(k) = norm(A*x_hat-y)^2/norm(y)^2;
%          r = norm(A*x_hat-y)^2/norm(y)^2;
%        
% %          if r(k) < error
%          if r < error
% %              out = x_hat;
%              %n_out=k;
% %              res_admm = r;
%              break
%          end
%     end
%     out = x_hat;
%     % n_out=k;
%     res_admm = r;
% end